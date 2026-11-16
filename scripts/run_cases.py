"""Build and run all OpenFOAM cases (needs OpenFOAM v2412 in the shell environment).

Usage:  python scripts/run_cases.py            # all cases, 2 at a time (skips finished ones)
        python scripts/run_cases.py sudden     # only cases whose name contains 'sudden'
"""
import sys, pathlib, subprocess, time, json
from concurrent.futures import ThreadPoolExecutor
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from case_builder import DuctCase, write_case

CASES = {
    # mesh-independence study on the sudden expansion (baseline mesh = scale 1.0)
    "sudden_coarse": DuctCase(diffuser_angle_deg=180, mesh_scale=0.6, end_iterations=5000),
    "sudden_medium": DuctCase(diffuser_angle_deg=180, mesh_scale=0.8, end_iterations=5000),
    "sudden":        DuctCase(diffuser_angle_deg=180, mesh_scale=1.0, end_iterations=5000),
    # design variants: planar diffusers (total included angle), baseline mesh
    "diffuser_60deg": DuctCase(diffuser_angle_deg=60, end_iterations=3000),
    "diffuser_30deg": DuctCase(diffuser_angle_deg=30, end_iterations=3000),
    "diffuser_20deg": DuctCase(diffuser_angle_deg=20, end_iterations=3000),
    "diffuser_10deg": DuctCase(diffuser_angle_deg=10, end_iterations=3000),
}


def run(name, c):
    d = ROOT / "cases" / name
    if (d / "done.json").exists():
        return name, json.loads((d / "done.json").read_text())
    write_case(c, d)
    t0 = time.time()
    for cmd in (["blockMesh"], ["simpleFoam"], ["postProcess", "-func", "writeCellCentres", "-latestTime"],
                ["postProcess", "-func", "writeCellVolumes", "-latestTime"]):
        with open(d / f"log.{cmd[0]}", "w") as log:
            subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
    log = (d / "log.simpleFoam").read_text()
    iters = log.count("\nTime = ")
    info = {"runtime_s": round(time.time() - t0, 1), "iterations": iters,
            "converged": "SIMPLE solution converged" in log, "cells": int(
                [l for l in (d / "log.blockMesh").read_text().splitlines() if "nCells:" in l][-1].split()[-1])}
    (d / "done.json").write_text(json.dumps(info))
    return name, info


if __name__ == "__main__":
    sel = {k: v for k, v in CASES.items() if len(sys.argv) < 2 or sys.argv[1] in k}
    with ThreadPoolExecutor(max_workers=2) as ex:
        for name, info in ex.map(lambda kv: run(*kv), sel.items()):
            print(name, info, flush=True)
