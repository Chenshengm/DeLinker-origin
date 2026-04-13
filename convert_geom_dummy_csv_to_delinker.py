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
from rdkit import Chem
from data.frag_utils import compute_distance_and_angle


def load_sdf(path):
    if not path:
        return {}
    sup = Chem.SDMolSupplier(path)
    out = []
    for idx, mol in enumerate(sup):
        if mol is None:
            out.append(None)
        else:
            out.append(mol)
    return out


def find_dummy_idx_by_mapnum(mol, map_num):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == map_num:
            return atom.GetIdx()
    return None


def merge_on_mapnum(base_smi, frag_smi, map_num):
    base = Chem.MolFromSmiles(base_smi)
    frag = Chem.MolFromSmiles(frag_smi)
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
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def split_frags_by_mapnum(frags_smi):
    frags = frags_smi.split(".")
    mapping = {}
    for frag in frags:
        m = Chem.MolFromSmiles(frag)
        if m is None:
            continue
        labels = [a.GetAtomMapNum() for a in m.GetAtoms() if a.GetAtomicNum() == 0]
        if len(labels) != 1:
            continue
        mapping[labels[0]] = frag
    return mapping


def heavy_atom_count(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return mol.GetNumHeavyAtoms()


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
    parser.add_argument("--fragments_col", default="fragments_with_dummy", help="fragment smiles column (dummy-labeled)")
    parser.add_argument("--linker_col", default="linker_with_dummy", help="linker smiles column (dummy-labeled)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--default_abs_dist", default="0.0")
    parser.add_argument("--default_angle", default="0.0")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    mol_sdf = load_sdf(args.mol_sdf)
    linker_sdf = load_sdf(args.linker_sdf)
    frag_sdf = load_sdf(args.frag_sdf)
    plan = []
    rows_out = []

    with open(args.input_csv, "r") as f:
        reader = csv.DictReader(f)
        for row_idx, row in enumerate(reader):
            case_id = row[args.uuid_col]
            mol = row[args.molecule_col]
            linker = row[args.linker_col]
            frags_map = split_frags_by_mapnum(row[args.fragments_col])
            if set(frags_map.keys()) != {1, 2, 3}:
                continue

            pair = sorted(rng.sample([1, 2, 3], 2))
            third = [k for k in [1, 2, 3] if k not in pair][0]

            pair_frags = frags_map[pair[0]] + "." + frags_map[pair[1]]
            pair_linker = merge_on_mapnum(linker, frags_map[third], third)
            if pair_linker is None:
                continue

            abs_dist = args.default_abs_dist
            angle = args.default_angle
            conf = mol_sdf[row_idx] if row_idx < len(mol_sdf) else None
            if conf is not None:
                d, a = compute_distance_and_angle(conf, pair_linker, pair_frags)
                if d is not None and a is not None:
                    abs_dist = str(float(d))
                    angle = str(float(a))

            rows_out.append((mol, pair_linker, pair_frags, abs_dist, angle))
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
    print("Wrote %d rows -> %s" % (len(rows_out), args.output_txt))
    print("Wrote plan -> %s" % args.output_plan)


if __name__ == "__main__":
    main()
