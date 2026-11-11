# EpiState-T
Reference implementation for longitudinal epigenetic measurements and a simulated reader of historical records. Authors: Jade Sunghyeok Choi and Jung In Choi. Software version **1.1.0**.

## Start on macOS with VS Code
1. Unzip the source ZIP. Open the extracted folder in VS Code using **File > Open Folder**.
2. Choose **Terminal > New Terminal**. In a Python 3.12 environment, run:
```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python demo.py
python check_demo.py
```
The first example uses intentionally synthetic population measurements (0.20, 0.55, 0.21). It reports a sampled return for the population but cannot establish a same-cell return or a daily rate. The separate model example reads a known virtual checkpoint; it is not evidence that biological cells do so.

## Use the actual fibroblast example
Download and extract the separate **EpiState-T_v1.1.0_Code_and_Data.zip** distributed with the corresponding preprint. Put the extracted data-release folder next to this source folder, then:
```sh
python demo.py --data-release ../EpiState-T_v1.1.0
```
If you saved the data elsewhere, pass its full path. The specified folder must contain `VERSION` and `addons/source/opening_example.json`. The command prints the observed **treated and mock** trajectories for O2/cg24947637 from Gill et al. [GSE165179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165179); its original probe is flagged and does not pass the sensitivity masks. A second illustrative site, cg06036945, passes both filters and also shows an observed treated but not mock excursion. The program labels the toy-model output separately. No bulk assay data is stored in this source tree. The data-and-code release contains the full reproducible study, study datasets or reconstruction scripts, output tables, checks and figure sources.

## Layout
- `model/tet_readh.py`: executable cell, neighbor, retention and Read(H) simulation from Paper 1.
- `query/src/`: reference tensor query rules, pandas comparators and sentence output. Use `demo.py` as the entry point in this source-only edition.
- `analysis_scripts/`: corresponding study-analysis source scripts. For a full empirical rerun, use the separate, matched data-and-code release with its exact directory layout and instructions.
- `CITATION.cff`, `VERSION`: citation and revision metadata for the study.
- `SOURCES_AND_RIGHTS.md`: source attribution and rights status.

This tree does not assign itself a DOI or GitHub address. The authors can supply these after deposit.
