#!/usr/bin/env python
"""
Convert GEOM 3-fragment dummy-labeled CSV (+ optional SDF with _Name=id)
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


def load_sdf_by_id(path):
    if not path:
        return {}
    sup = Chem.SDMolSupplier(path)
    out = {}
    for mol in sup:
        if mol is None:
            continue
        if not mol.HasProp("_Name"):
            continue
        out[mol.GetProp("_Name")] = mol
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_txt", required=True, help="DeLinker 5-column input text")
    parser.add_argument("--output_plan", required=True, help="stage bookkeeping json")
    parser.add_argument("--sdf", default="", help="optional SDF with _Name=id for distance/angle")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--default_abs_dist", default="0.0")
    parser.add_argument("--default_angle", default="0.0")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    sdf_by_id = load_sdf_by_id(args.sdf)
    plan = []
    rows_out = []

    with open(args.input_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            case_id = row["id"]
            mol = row["molecule"]
            linker = row["linker"]
            frags_map = split_frags_by_mapnum(row["fragments"])
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
            conf = sdf_by_id.get(str(case_id))
            if conf is not None:
                d, a = compute_distance_and_angle(conf, pair_linker, pair_frags)
                if d is not None and a is not None:
                    abs_dist = str(float(d))
                    angle = str(float(a))

            rows_out.append((mol, pair_linker, pair_frags, abs_dist, angle))
            plan.append({
                "id": int(case_id),
                "pair": pair,
                "third": third,
                "pair_linker": pair_linker,
                "pair_frags": pair_frags,
            })

    with open(args.output_txt, "w") as f:
        for r in rows_out:
            f.write("%s %s %s %s %s\n" % r)

    with open(args.output_plan, "w") as f:
        json.dump(plan, f, indent=2)

    print("Wrote %d rows -> %s" % (len(rows_out), args.output_txt))
    print("Wrote plan -> %s" % args.output_plan)


if __name__ == "__main__":
    main()
