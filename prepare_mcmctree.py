#!/usr/bin/env python3
"""Turn a filled-in PhyloSlide calibration file into MCMCtree inputs.

PhyloSlide's --dating_template writes a directory containing:

    node_key.txt        every internal node, in plain English
    tree.template.nwk   rooted tree with @N1@.. placeholders
    calibrations.txt    for you to fill in
    mcmctree.ctl        control file with ndata pre-filled

Fill in calibrations.txt, then run this script on that directory. It substitutes
your calibrations into the tree, strips the unused placeholders, and writes the
two control files for MCMCtree's two-pass approximate-likelihood workflow:

    tree.nwk                    calibrated, rooted tree
    mcmctree_step1_hessian.ctl  usedata=3, runs baseml -> out.BV
    mcmctree_step2_mcmc.ctl     usedata=2, reads in.BV

Then, in that directory:

    cp mcmctree_step1_hessian.ctl mcmctree.ctl && mcmctree mcmctree.ctl
    mv out.BV in.BV
    cp mcmctree_step2_mcmc.ctl mcmctree.ctl && mcmctree mcmctree.ctl

baseml must be on PATH for step 1: mcmctree calls it through the shell, and if it
is missing every locus fails with "file rst2 not found!" while mcmctree still
exits 0 and writes an out.BV that is silently useless.

Usage:
    prepare_mcmctree.py <template_dir> [--burnin N] [--sampfreq N] [--nsample N]
                                       [--clock {1,2,3}] [--ndata N]
"""
import argparse
import re
import sys
from pathlib import Path

CALIB_RE = re.compile(r"^\s*(?:B|L|U|G|SN|ST)\s*\(", re.I)


def parse_calibrations(path: Path):
    calib, seen = [], {}
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            sys.exit(f"{path}:{lineno}: expected '<NODE>  <calibration>', got: {raw.strip()}")
        node, cal = parts[0].strip(), parts[1].strip()
        if not re.fullmatch(r"N\d+", node):
            sys.exit(f"{path}:{lineno}: node must look like N1, N2, ...; got '{node}'")
        if not CALIB_RE.match(cal):
            sys.exit(
                f"{path}:{lineno}: '{cal}' does not look like an MCMCtree calibration.\n"
                "Expected one of B(...), L(...), U(...), G(...), SN(...), ST(...).\n"
                "Remember times are in units of 100 Myr: 17.2 Ma -> 0.172"
            )
        if cal.count("(") != cal.count(")"):
            sys.exit(f"{path}:{lineno}: unbalanced parentheses in '{cal}'")
        if node in seen:
            sys.exit(f"{path}:{lineno}: {node} already calibrated on line {seen[node]}")
        seen[node] = lineno
        calib.append((node, cal))
    return calib


def clade_tips(tree_str):
    """Return {node_label: [tips]} for every @Nk@ label in a newick string."""
    out, stack, buf = {}, [], ""
    i = 0
    while i < len(tree_str):
        c = tree_str[i]
        if c == "(":
            stack.append([]); buf = ""
        elif c == ",":
            if buf.strip():
                stack[-1].append(buf.strip().split(":")[0])
            buf = ""
        elif c == ")":
            if buf.strip():
                stack[-1].append(buf.strip().split(":")[0])
            buf = ""
            grp = stack.pop()
            j = i + 1
            lab = ""
            while j < len(tree_str) and tree_str[j] not in ",);":
                lab += tree_str[j]; j += 1
            lab = lab.split(":")[0].strip()
            if lab.startswith("@") and lab.endswith("@"):
                out[lab.strip("@")] = sorted(grp)
            if stack:
                stack[-1].extend(grp)
            i = j - 1
        else:
            buf += c
        i += 1
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Build MCMCtree inputs from a filled-in PhyloSlide calibration file.")
    ap.add_argument("template_dir", type=Path)
    ap.add_argument("--ndata", type=int, default=None,
                    help="Override the ndata from the template control file.")
    ap.add_argument("--clock", choices=["1", "2", "3"], default=None,
                    help="1 strict, 2 uncorrelated, 3 autocorrelated.")
    ap.add_argument("--burnin", type=int, default=None)
    ap.add_argument("--sampfreq", type=int, default=None)
    ap.add_argument("--nsample", type=int, default=None)
    a = ap.parse_args()

    d = a.template_dir
    tmpl, calf, ctlf = d / "tree.template.nwk", d / "calibrations.txt", d / "mcmctree.ctl"
    for f in (tmpl, calf, ctlf):
        if not f.exists():
            sys.exit(f"ERROR: {f} not found. Run PhyloSlide with --dating_template first.")

    calib = parse_calibrations(calf)
    if not calib:
        sys.exit(
            f"ERROR: no calibrations in {calf}.\n"
            "Add at least one line like:  N1   B(0.172, 0.195, 1e-300, 0.025)\n"
            "See node_key.txt for which node is which."
        )

    header, tree = tmpl.read_text().strip().split("\n", 1)
    tree = tree.strip()
    available = set(re.findall(r"@(N\d+)@", tree))

    # Read back which clade each calibrated node actually is, so the user sees what
    # they constrained rather than trusting a node number. Node numbering is
    # specific to this tree; the same number means something different elsewhere.
    clades = clade_tips(tree)

    for node, cal in calib:
        if node not in available:
            sys.exit(
                f"ERROR: {node} is not a node in this tree. "
                f"Valid nodes: N1..N{len(available)} (see node_key.txt)"
            )
        tree = tree.replace(f"@{node}@", f"'{cal}'")

    tree = re.sub(r"@N\d+@", "", tree)              # drop uncalibrated placeholders
    tree = re.sub(r":0\.?0*(?=[,)])", "", tree)     # tidy dummy branch lengths
    tree = re.sub(r":0\.?0*;", ";", tree)
    (d / "tree.nwk").write_text(f"{header}\n{tree}\n")

    ctl = ctlf.read_text()
    subs = {"ndata": a.ndata, "clock": a.clock,
            "burnin": a.burnin, "sampfreq": a.sampfreq, "nsample": a.nsample}
    for k, v in subs.items():
        if v is not None:
            ctl = re.sub(rf"^(\s*{k}\s*=\s*)\S+", rf"\g<1>{v}", ctl, flags=re.M)

    step1 = re.sub(r"^(\s*usedata\s*=\s*)\S+", r"\g<1>3", ctl, flags=re.M)
    step1 = re.sub(r"^(\s*outfile\s*=\s*)\S+", r"\g<1>out_hessian.txt", step1, flags=re.M)
    step2 = re.sub(r"^(\s*usedata\s*=\s*)\S+", r"\g<1>2", ctl, flags=re.M)
    step2 = re.sub(r"^(\s*outfile\s*=\s*)\S+", r"\g<1>out_dates.txt", step2, flags=re.M)
    (d / "mcmctree_step1_hessian.ctl").write_text(step1)
    (d / "mcmctree_step2_mcmc.ctl").write_text(step2)

    m = re.search(r"^\s*ndata\s*=\s*(\d+)", ctl, flags=re.M)
    print(f"calibrations applied : {len(calib)}\n")
    for node, cal in calib:
        tips = clades.get(node, [])
        if len(tips) <= 5:
            desc = ", ".join(tips)
        else:
            desc = f"{len(tips)} taxa: {', '.join(tips[:3])}, ... {tips[-1]}"
        print(f"    {node:5} {cal}")
        print(f"          -> {desc}")
    print("\n    CHECK these are the clades you meant. Node numbers are specific to")
    print("    this tree and are not comparable with any other analysis.")
    print(f"ndata                : {m.group(1) if m else '?'}")
    print(f"\nwritten to {d}:")
    print("    tree.nwk")
    print("    mcmctree_step1_hessian.ctl   (usedata=3, baseml -> out.BV)")
    print("    mcmctree_step2_mcmc.ctl      (usedata=2, approx.-lik. MCMC)")
    print(f"\nNext, in {d}:")
    print("    cp mcmctree_step1_hessian.ctl mcmctree.ctl && mcmctree mcmctree.ctl")
    print("    grep -c 'rst2 not found' mcmctree.log   # must be 0, else baseml is not on PATH")
    print("    mv out.BV in.BV")
    print("    cp mcmctree_step2_mcmc.ctl mcmctree.ctl && mcmctree mcmctree.ctl")


if __name__ == "__main__":
    main()
