#!/usr/bin/env python
"""
Build a two-stage sampling plan for GEOM-style 3-fragment inputs.

Input line format (whitespace separated):
    <full_smi> <frag1.frag2.frag3>

Example:
    COC1... O=CNC1CCCC1.C1=CSC=C1.COC1=C(OC)C=CC=C1

Output:
    A JSON list where each item stores:
      - full_smi
      - frags (3 fragments)
      - randomly selected stage-1 pair
      - remaining third fragment for stage-2

This script only builds the *plan* and does not run model inference.
"""

import argparse
import json
import random


def parse_line(line):
    toks = line.strip().split()
    if len(toks) < 2:
        return None
    full_smi = toks[0]
    frag_smi = toks[1]
    frags = frag_smi.split('.')
    if len(frags) != 3:
        return None
    return full_smi, frags


def build_record(case_id, full_smi, frags, rng):
    pair = sorted(rng.sample([0, 1, 2], 2))
    remaining = [idx for idx in [0, 1, 2] if idx not in pair][0]
    return {
        "case_id": case_id,
        "full_smi": full_smi,
        "frags": frags,
        "stage1_pair_indices": pair,
        "stage1_frag_smi": ".".join([frags[pair[0]], frags[pair[1]]]),
        "stage2_remaining_index": remaining,
        "stage2_remaining_frag": frags[remaining],
        "stage2_input_template": "<STAGE1_GENERATED_MOL>." + frags[remaining],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="path to txt/smi input")
    parser.add_argument("--output", required=True, help="path to json output")
    parser.add_argument("--seed", type=int, default=0, help="random seed")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    records = []
    with open(args.input, "r") as f:
        for idx, line in enumerate(f):
            parsed = parse_line(line)
            if parsed is None:
                continue
            full_smi, frags = parsed
            records.append(build_record(idx, full_smi, frags, rng))

    with open(args.output, "w") as f:
        json.dump(records, f, indent=2)

    print("Saved %d GEOM two-stage records to %s" % (len(records), args.output))


if __name__ == "__main__":
    main()
