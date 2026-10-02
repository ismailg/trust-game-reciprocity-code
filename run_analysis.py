#!/usr/bin/env python3
"""Run the code in a separate private directory; never populate the release with data."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
from verify_package import verify

ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mode", choices=("fresh", "reuse-fits"), default="fresh")
    parser.add_argument("--private-archive", type=Path,
                        help="Authorised original analysis directory, only for reuse-fits mode")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument("--no-render", action="store_true")
    args = parser.parse_args()
    work = args.work_dir.resolve()
    data = args.data_dir.resolve()
    if work == ROOT or ROOT in work.parents or work in ROOT.parents:
        parser.error("The private work directory must be outside the release directory")
    if work.exists():
        parser.error("Choose a new work directory; existing results are never overwritten")
    if args.workers < 1 or args.bootstrap < 8 or args.bootstrap % 8:
        parser.error("workers must be positive; bootstrap must be a positive multiple of eight")
    if args.mode != "fresh" and args.private_archive is None:
        parser.error("This mode requires --private-archive")
    if args.mode == "fresh" and args.private_archive is not None:
        parser.error("Fresh mode takes only the two data files, not archived results")
    if not args.no_render and args.bootstrap != 10000:
        parser.error("Small bootstrap checks cannot render a manuscript containing publication counts")
    for name in ("full_RTG_data.csv", "demographics.csv"):
        if not (data / name).is_file():
            parser.error(f"Required input missing: {name}")
    verify(ROOT)
    # Only package code/documentation are copied here. Private inputs are added
    # to this disposable workspace, never to ROOT or a distribution archive.
    shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git"))
    for name in ("full_RTG_data.csv", "demographics.csv"):
        shutil.copy2(data / name, work / "Data" / name)
    (work / ".gitignore").write_text("*\n")
    (work / "PRIVATE_WORKSPACE.txt").write_text(
        "Contains restricted participant data and analysis objects. Do not distribute.\n")
    archive = args.private_archive.resolve() if args.private_archive else None
    if args.mode == "reuse-fits":
        (work / "private_fits").mkdir()
        fits = {
            "investor.rds": "results/HMM/corrected_submission/standardized_investor_zero_return_2026-10-02_v1/standardized_search_result.rds",
            "trustee.rds": "results/HMM/corrected_submission/standardized_trustee_search_2026-07-21_v2/standardized_search_result.rds",
        }
        for name, source in fits.items():
            shutil.copy2(archive / source, work / "private_fits" / name)
    env = os.environ.copy()
    # Do not inherit unrelated run overrides from an interactive session.
    env = {k: v for k, v in env.items() if not k.startswith(("VTC_", "RTG_", "TRUSTEE_"))}
    env.update(RTG_RUN_MODE=args.mode, RTG_WORKERS=str(args.workers),
               RTG_BOOTSTRAP_DRAWS=str(args.bootstrap), OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    commands = [["Rscript", "scripts/recompute.R"]]
    if not args.no_render:
        commands.append(["Rscript", "render.R", "both"])
    status = {"mode": args.mode, "bootstrap_draws": args.bootstrap, "commands": []}
    try:
        for index, command in enumerate(commands, 1):
            print("Running", " ".join(command), flush=True)
            with (work / f"stage_{index}.stdout.log").open("w") as out, (work / f"stage_{index}.stderr.log").open("w") as err:
                result = subprocess.run(command, cwd=work, env=env, stdout=out, stderr=err)
            status["commands"].append({"command": command, "exit_code": result.returncode})
            if result.returncode:
                raise RuntimeError(f"Stage {index} failed; inspect its log in the private work directory")
        status["completed"] = True
    except Exception:
        status["completed"] = False
        raise
    finally:
        (work / "RUN_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")
    print("Completed. Generated files are private and must not be added to the code release.")


if __name__ == "__main__":
    main()
