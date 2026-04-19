#!/usr/bin/env python
"""
Build DeLinker-compatible test inputs from a GEOM two-stage plan.

Stage-1 output format (for `data/prepare_data.py --test_mode`):
    <frag_a.frag_b> <abs_dist> <angle>

Stage-2 output format (for `data/prepare_data.py --test_mode`):
    <stage1_generated_mol.frag_c> <abs_dist> <angle>
"""

import argparse
import json


def load_plan(path):
    with open(path, "r") as f:
        return json.load(f)


def write_stage1(plan, out_path, abs_dist, angle):
    with open(out_path, "w") as f:
        for item in plan:
            f.write("%s %s %s\n" % (item["stage1_frag_smi"], abs_dist, angle))
    print("Wrote stage-1 input: %s (%d lines)" % (out_path, len(plan)))


def parse_stage1_gen_line(line):
    toks = line.strip().split()
    # default DeLinker generated line format:
    #   <smiles_in> <smiles_out> <generated_smiles>
    if len(toks) < 3:
        return None
    return toks[2]


def load_stage1_generated(path):
    generated = []
    with open(path, "r") as f:
        for line in f:
            gen = parse_stage1_gen_line(line)
            if gen is not None:
                generated.append(gen)
    return generated


def write_stage2(plan, stage1_generated, out_path, abs_dist, angle, stage1_per_input=1):
    required = len(plan) * stage1_per_input
    if len(stage1_generated) < required:
        raise ValueError("stage1 generated count (%d) < required (%d)" % (len(stage1_generated), required))

    with open(out_path, "w") as f:
        line_count = 0
        for idx, item in enumerate(plan):
            remaining = item.get("stage2_remaining_frag", item.get("third_frag_with_dummy"))
            if remaining is None:
                raise ValueError("plan item missing remaining fragment field")
            start = idx * stage1_per_input
            end = start + stage1_per_input
            for gen in stage1_generated[start:end]:
                stage2_frag = "%s.%s" % (gen, remaining)
                f.write("%s %s %s\n" % (stage2_frag, abs_dist, angle))
                line_count += 1
    print("Wrote stage-2 input: %s (%d lines)" % (out_path, line_count))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("stage1")
    p1.add_argument("--plan", required=True)
    p1.add_argument("--output", required=True)
    p1.add_argument("--abs_dist", default="0.0")
    p1.add_argument("--angle", default="0.0")

    p2 = sub.add_parser("stage2")
    p2.add_argument("--plan", required=True)
    p2.add_argument("--stage1_generated_smi", required=True)
    p2.add_argument("--output", required=True)
    p2.add_argument("--abs_dist", default="0.0")
    p2.add_argument("--angle", default="0.0")
    p2.add_argument("--stage1_per_input", type=int, default=1,
                    help="how many stage-1 generated molecules correspond to each plan item")

    args = parser.parse_args()
    plan = load_plan(args.plan)

    if args.cmd == "stage1":
        write_stage1(plan, args.output, args.abs_dist, args.angle)
    elif args.cmd == "stage2":
        stage1_generated = load_stage1_generated(args.stage1_generated_smi)
        write_stage2(plan, stage1_generated, args.output, args.abs_dist, args.angle, args.stage1_per_input)


if __name__ == "__main__":
    main()
