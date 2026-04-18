#!/usr/bin/env python
"""
Convert GEOM 3-fragment dummy-labeled CSV (+ optional SDF)
into DeLinker stage-1 input (5-column format):

    <molecule> <pair_linker> <pair_fragments> <abs_dist> <angle>

This keeps two fragments as `pair_fragments` and merges the third fragment
into the linker so that linker/fragments each have exactly two attachment dummies.
"""

import argparse
import csv
import json
import random
import re
from rdkit import Chem
from rdkit import RDLogger
from data.frag_utils import compute_distance_and_angle

RDLogger.DisableLog("rdApp.error")
ENABLE_SANITIZE = True
ENABLE_KEKULIZE = True


def load_sdf(path):
    if not path:
        return {}
    sup = Chem.SDMolSupplier(path, sanitize=ENABLE_SANITIZE, removeHs=False)
    out = []
    for idx, mol in enumerate(sup):
        if mol is None:
            out.append(None)
        else:
            out.append(mol)
    return out


def build_mol_lookup_by_smiles(mols):
    lookup = {}
    for mol in mols:
        if mol is None:
            continue
        try:
            key = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)), isomericSmiles=True)
            lookup[key] = mol
        except Exception:
            continue
    return lookup


def find_dummy_idx_by_mapnum(mol, map_num):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == map_num:
            return atom.GetIdx()
    return None


def mol_from_any(text):
    mol = Chem.MolFromSmiles(text, sanitize=ENABLE_SANITIZE)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(text, sanitize=False)
    if mol is not None:
        return mol
    q = Chem.MolFromSmarts(text)
    if q is None:
        return None
    try:
        smi = Chem.MolToSmiles(q, isomericSmiles=True)
        return Chem.MolFromSmiles(smi, sanitize=ENABLE_SANITIZE)
    except Exception:
        return None


def merge_on_mapnum(base_smi, frag_smi, map_num):
    base = mol_from_any(base_smi)
    frag = mol_from_any(frag_smi)
    if base is None or frag is None:
        return None

    idx_b = find_dummy_idx_by_mapnum(base, map_num)
    idx_f = find_dummy_idx_by_mapnum(frag, map_num)
    if idx_b is None or idx_f is None:
        return None

    nei_b = base.GetAtomWithIdx(idx_b).GetNeighbors()[0].GetIdx()
    nei_f = frag.GetAtomWithIdx(idx_f).GetNeighbors()[0].GetIdx()
    bond_b = base.GetBondBetweenAtoms(idx_b, nei_b).GetBondType()
    bond_f = frag.GetBondBetweenAtoms(idx_f, nei_f).GetBondType()
    bond_type = bond_b if bond_b == bond_f else Chem.BondType.SINGLE

    combo = Chem.CombineMols(base, frag)
    rw = Chem.RWMol(combo)
    offset = base.GetNumAtoms()

    rw.AddBond(nei_b, nei_f + offset, bond_type)

    to_remove = sorted([idx_b, idx_f + offset], reverse=True)
    for idx in to_remove:
        rw.RemoveAtom(idx)

    mol = rw.GetMol()
    if ENABLE_SANITIZE:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
    return Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=False)


def split_frags_by_mapnum(frags_smi):
    frags = frags_smi.split(".")
    mapping = {}
    for frag in frags:
        labels = re.findall(r"\[\*:([0-9]+)\]", frag)
        if len(labels) != 1:
            continue
        mapping[int(labels[0])] = frag
    return mapping


def heavy_atom_count(smiles):
    mol = mol_from_any(smiles)
    if mol is None:
        return 0
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def normalize_for_delinker(smiles):
    mol = mol_from_any(smiles)
    if mol is None:
        return None
    if ENABLE_KEKULIZE:
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            pass
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=ENABLE_KEKULIZE)
    except Exception:
        return Chem.MolToSmiles(mol, isomericSmiles=True)


def parse_anchor_indices(text):
    if not text:
        return []
    return [int(x) for x in text.split("-") if x != ""]


def add_dummy_to_fragment(fragment_smiles, full_mol, anchor_idx, label):
    frag = mol_from_any(fragment_smiles)
    if frag is None:
        return None
    matches = list(full_mol.GetSubstructMatches(frag, uniquify=False))
    chosen = None
    for m in matches:
        if anchor_idx in m:
            chosen = m
            break
    if chosen is None:
        return None
    local_idx = list(chosen).index(anchor_idx)
    rw = Chem.RWMol(frag)
    dummy = Chem.Atom(0)
    dummy.SetAtomMapNum(label)
    d_idx = rw.AddAtom(dummy)
    rw.AddBond(local_idx, d_idx, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    if ENABLE_SANITIZE:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
    return Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=ENABLE_KEKULIZE)


def add_dummy_to_linker(linker_smiles, full_mol, anchors, label_map):
    linker = mol_from_any(linker_smiles)
    if linker is None:
        return None
    matches = list(full_mol.GetSubstructMatches(linker, uniquify=False))
    if not matches:
        return None
    best = max(matches, key=lambda m: sum(1 for a in anchors if a in m))
    rw = Chem.RWMol(linker)
    # add in descending local idx order to keep indices stable
    add_ops = []
    for a in anchors:
        if a not in best:
            continue
        local_idx = list(best).index(a)
        add_ops.append((local_idx, label_map[a]))
    for local_idx, label in sorted(add_ops, key=lambda x: x[0], reverse=True):
        dummy = Chem.Atom(0)
        dummy.SetAtomMapNum(label)
        d_idx = rw.AddAtom(dummy)
        rw.AddBond(local_idx, d_idx, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    if ENABLE_SANITIZE:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
    return Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=ENABLE_KEKULIZE)


def build_dummy_strings_from_nondummy(row, molecule_col, fragments_nondummy_col, linker_nondummy_col, anchors_col):
    full_mol = mol_from_any(row[molecule_col])
    if full_mol is None:
        return None, None
    anchors = parse_anchor_indices(row.get(anchors_col, ""))
    if len(anchors) < 2:
        return None, None
    # label by anchor order: first anchor->1, second->2, ...
    label_map = {a: idx + 1 for idx, a in enumerate(anchors)}

    frag_list = row[fragments_nondummy_col].split(".")
    remaining_anchors = set(anchors)
    frag_with_dummy = []
    for frag in frag_list:
        chosen_anchor = None
        for a in list(remaining_anchors):
            trial = add_dummy_to_fragment(frag, full_mol, a, label_map[a])
            if trial is not None:
                chosen_anchor = a
                frag_with_dummy.append(trial)
                break
        if chosen_anchor is None:
            return None, None
        remaining_anchors.remove(chosen_anchor)

    linker_with_dummy = add_dummy_to_linker(row[linker_nondummy_col], full_mol, anchors, label_map)
    if linker_with_dummy is None:
        return None, None
    return ".".join(frag_with_dummy), linker_with_dummy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_txt", required=True, help="DeLinker 5-column input text")
    parser.add_argument("--output_plan", required=True, help="stage bookkeeping json")
    parser.add_argument("--mol_sdf", default="", help="optional mol SDF for distance/angle; row order should match csv")
    parser.add_argument("--linker_sdf", default="", help="optional linker SDF (order check only)")
    parser.add_argument("--frag_sdf", default="", help="optional fragment SDF (order check only)")
    parser.add_argument("--uuid_col", default="uuid", help="uuid column name in csv")
    parser.add_argument("--molecule_col", default="molecule", help="full molecule smiles column")
    parser.add_argument("--fragments_nondummy_col", default="fragments", help="fragment smiles column (without dummy)")
    parser.add_argument("--linker_nondummy_col", default="linker", help="linker smiles column (without dummy)")
    parser.add_argument("--anchors_col", default="anchors", help="anchor indices column from full molecule")
    parser.add_argument("--fragments_col", default="fragments_with_dummy", help="fragment smiles column (dummy-labeled)")
    parser.add_argument("--linker_col", default="linker_with_dummy", help="linker smiles column (dummy-labeled)")
    parser.add_argument("--build_dummy_from_nondummy", action="store_true",
                        help="rebuild dummy-labeled fragments/linker from non-dummy columns + anchors")
    parser.add_argument("--disable_sanitize", action="store_true",
                        help="disable RDKit sanitize in parsing/SDMolSupplier")
    parser.add_argument("--disable_kekulize", action="store_true",
                        help="disable kekulize during output normalization")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--default_abs_dist", default="0.0")
    parser.add_argument("--default_angle", default="0.0")
    args = parser.parse_args()
    global ENABLE_SANITIZE, ENABLE_KEKULIZE
    ENABLE_SANITIZE = not args.disable_sanitize
    ENABLE_KEKULIZE = not args.disable_kekulize

    rng = random.Random(args.seed)
    mol_sdf = load_sdf(args.mol_sdf)
    mol_lookup = build_mol_lookup_by_smiles(mol_sdf)
    linker_sdf = load_sdf(args.linker_sdf)
    frag_sdf = load_sdf(args.frag_sdf)
    plan = []
    rows_out = []
    skip_bad_dummy = 0
    skip_merge_fail = 0
    skip_normalize_fail = 0
    dist_angle_fail = 0
    dummy_rebuild_fail = 0

    with open(args.input_csv, "r") as f:
        reader = csv.DictReader(f)
        for row_idx, row in enumerate(reader):
            case_id = row[args.uuid_col]
            mol = row[args.molecule_col]
            if args.build_dummy_from_nondummy:
                frags_dummy, linker_dummy = build_dummy_strings_from_nondummy(
                    row, args.molecule_col, args.fragments_nondummy_col, args.linker_nondummy_col, args.anchors_col
                )
                if frags_dummy is None or linker_dummy is None:
                    dummy_rebuild_fail += 1
                    continue
                linker = linker_dummy
                frags_map = split_frags_by_mapnum(frags_dummy)
            else:
                linker = row[args.linker_col]
                frags_map = split_frags_by_mapnum(row[args.fragments_col])
            if set(frags_map.keys()) != {1, 2, 3}:
                skip_bad_dummy += 1
                continue

            pair = sorted(rng.sample([1, 2, 3], 2))
            third = [k for k in [1, 2, 3] if k not in pair][0]

            pair_frags = frags_map[pair[0]] + "." + frags_map[pair[1]]
            pair_linker = merge_on_mapnum(linker, frags_map[third], third)
            if pair_linker is None:
                skip_merge_fail += 1
                continue

            mol_norm = normalize_for_delinker(mol)
            linker_norm = normalize_for_delinker(pair_linker)
            frags_norm = normalize_for_delinker(pair_frags)
            if None in (mol_norm, linker_norm, frags_norm):
                skip_normalize_fail += 1
                continue

            abs_dist = args.default_abs_dist
            angle = args.default_angle
            conf = mol_sdf[row_idx] if row_idx < len(mol_sdf) else None
            if conf is None and mol_norm is not None:
                mol_key = Chem.MolToSmiles(Chem.MolFromSmiles(mol_norm), isomericSmiles=True)
                conf = mol_lookup.get(mol_key)
            if conf is not None:
                d, a = compute_distance_and_angle(conf, linker_norm, frags_norm)
                if d is not None and a is not None:
                    abs_dist = str(float(d))
                    angle = str(float(a))
                else:
                    dist_angle_fail += 1

            rows_out.append((mol_norm, linker_norm, frags_norm, abs_dist, angle))
            original_linker_heavy = heavy_atom_count(linker)
            stage1_target = max(1, original_linker_heavy // 2)
            stage2_target = max(1, original_linker_heavy - stage1_target)
            plan.append({
                "uuid": int(case_id),
                "pair": pair,
                "third": third,
                "pair_linker": pair_linker,
                "pair_frags": pair_frags,
                "third_frag_with_dummy": frags_map[third],
                "original_linker_heavy_atoms": original_linker_heavy,
                "stage1_target_linker_heavy_atoms": stage1_target,
                "stage2_target_linker_heavy_atoms": stage2_target,
                "row_idx": row_idx,
            })

    with open(args.output_txt, "w") as f:
        for r in rows_out:
            f.write("%s %s %s %s %s\n" % r)

    with open(args.output_plan, "w") as f:
        json.dump(plan, f, indent=2)

    if args.linker_sdf:
        print("Loaded linker_sdf entries: %d" % len(linker_sdf))
    if args.frag_sdf:
        print("Loaded frag_sdf entries: %d" % len(frag_sdf))
    print("Skipped rows (dummy parse): %d" % skip_bad_dummy)
    print("Skipped rows (merge fail): %d" % skip_merge_fail)
    print("Skipped rows (normalize fail): %d" % skip_normalize_fail)
    print("Rows with dist/angle fallback defaults: %d" % dist_angle_fail)
    if args.build_dummy_from_nondummy:
        print("Skipped rows (dummy rebuild fail): %d" % dummy_rebuild_fail)
    print("Wrote %d rows -> %s" % (len(rows_out), args.output_txt))
    print("Wrote plan -> %s" % args.output_plan)


if __name__ == "__main__":
    main()
