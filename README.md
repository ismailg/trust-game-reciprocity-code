# Trust Game analysis code

Code accompanying **Beyond Behavioural Averages: Clinical Differences in Trust and Reciprocity Dynamics**. Submission version: 2 October 2026.

The package contains analysis code, manuscript source, run instructions, and software versions. Participant data and fitted models are excluded. Reproducing the empirical results requires authorised access to `full_RTG_data.csv` and `demographics.csv`; see [Data/README.md](Data/README.md).

Contact: Ismail Guennouni, ismail.guennouni@iwr.uni-heidelberg.de. Manuscript authors are listed in [CITATION.txt](CITATION.txt).

The fixed submission snapshot is the `submission-v4` release. Later updates to the repository may differ.

## Run the analyses

1. Check the package:

   ```sh
   python3 verify_package.py
   ```

2. Prepare R 4.5.1 and the package versions in `renv.lock`, following [environment/README.md](environment/README.md). Python, Pandoc, and LaTeX are also required.
3. Keep the two authorised input files in a directory outside this package.
4. Run from the package directory, replacing the example paths:

   ```sh
   python3 run_analysis.py \
     --data-dir /path/to/authorised-inputs \
     --work-dir /path/to/new-private-run \
     --mode fresh --workers 8
   ```

The working directory must be new and outside the package. The command copies the code and inputs there, runs the analyses, and renders the main PDF, main Word document, and supplement PDF. Generated files may contain participant information and must remain private.

Successful completion creates `RUN_STATUS.json` with `completed: true`. The documents are `VTC_RTG_analysis.pdf`, `VTC_RTG_analysis.docx`, and `supplement.pdf`. Numerical results and saved analysis state are in `output/`; command logs are `stage_*.stdout.log` and `stage_*.stderr.log`.

To regenerate documents from a completed private run, use `Rscript render.R both` inside that working directory. Rendering checks that the source and input files match the saved analysis state.

## Optional run modes

| Mode | Purpose |
|---|---|
| `fresh` | Fit models from the two input files, recalculate analyses, and render documents. |
| `reuse-fits` | Check and reuse the manuscript's fitted models, then recalculate subsequent analyses and render. |
| `archived` | Render from the manuscript's existing result archives, without rerunning the bootstrap. |

The latter two modes require an authorised original analysis archive and the additional argument `--private-archive /path/to/original-analysis`. Those archives are not distributed here.

For a short pipeline check, add `--bootstrap 80 --no-render`. This reduced run cannot support the manuscript's statistical conclusions; document rendering requires the default 10,000 draws.

## Analysis settings

Primary analyses include eligible participants across administration modes. Investor hidden-state comparisons include community controls and participants with borderline personality disorder; trustee comparisons also include the affective-symptom sample. Sample definitions and group coding follow the manuscript.

The bootstrap varies transition counts and parameters while holding the fitted state structure and emission distributions fixed. It does not include uncertainty from selecting the number of states. Each analysis pools eight batches of 1,250 successful draws before calculating intervals and corrected tests.

| Analysis | Seeds |
|---|---|
| Primary five-state trustee model | 11 through 18 |
| Four-state check | 15 through 7015, in steps of 1000 |
| Six-state check | 17 through 7017, in steps of 1000 |
| All modes, age/gender adjusted | 21011 through 28011, in steps of 1000 |
| Online, age/gender adjusted | 31011 through 38011, in steps of 1000 |
| All online participants | 41011 through 48011, in steps of 1000 |
| Online with complete demographics, unadjusted | 51011 through 58011, in steps of 1000 |

`--workers` changes the number of concurrent batches, not their seeds. Model-search settings and seeds are in `analysis_main.Rmd`. Numerical optimizers can differ across software versions and platforms.

## File guide

- `run_analysis.py` and `scripts/recompute.R`: entry point and computation order.
- `analysis_main.Rmd` and `supplement.Rmd`: analyses and manuscript content; `VTC_RTG_analysis.Rmd` and `render.R` coordinate rendering.
- `scripts/`: model, bootstrap, and sensitivity-analysis functions. Original `run_trustee_*.R` scripts support source-file checks used by the manuscript; use the main entry point to run this package.
- `data_dictionary.csv`: required input fields and coding, without participant records.
- `output_map.csv`: code locations producing each main and supplementary table or figure.
- `renv.lock` and `environment/`: software versions and setup instructions.
- `RELEASE_MANIFEST.json` and `verify_package.py`: distributed file list and integrity checks.

## Licence

See [LICENSE_STATUS.txt](LICENSE_STATUS.txt) for the current code and third-party licensing status.
