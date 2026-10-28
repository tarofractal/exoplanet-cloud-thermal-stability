# Exoplanet cloud thermal stability

Research code associated with **Taro M. Shefferson-Nagata and Kazumasa Ohno (2026)**, *Radiatively Controlled Thermal Stability of High-altitude Clouds in Exoplanetary Atmospheres*, The Astrophysical Journal, 1005, 210.

Paper: https://doi.org/10.3847/1538-4357/ae7d27

**Status: partial recovery of research scripts.** Two historical source files have been recovered. This repository is not yet a complete release of the code and input data used for the published figures. See [recovery status](docs/RECOVERY_STATUS.md) before interpreting outputs.

## Contents

| File | Purpose |
| --- | --- |
| `scripts/qabs_loggrid.py` | Mie absorption efficiencies for eight condensates, wavelength integrals, optional Planck means, CSV/PDF output. |
| `scripts/thermal_decoupling.py` | Particle energy balance including gas collisions and the pressure at which the fractional temperature difference reaches 5%. |
| `archive/` | The two recovered source files, preserved without edits. |
| `data/README.md` | Required input filenames, units, and missing data. |
| `tests/test_numerics.py` | Synthetic numerical checks; these do not reproduce the published figures. |

The runnable Qabs script uses `miepython`, the library cited in the paper. Its recovered source uses PyMieScatt. This migration and the small execution-related changes are documented separately from the historical source.

## Installation

Use Python 3.11 or later in a virtual environment:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

## Run

The original optical-constant files and Qabs CSV inputs have **not** been recovered. Supply the files described in [data/README.md](data/README.md). Do not substitute the synthetic test data for research inputs.

```sh
python scripts/qabs_loggrid.py --data_dir data/optical_constants --save_dir outputs/qabs --no_plot
python scripts/thermal_decoupling.py --data-dir outputs/qabs --output-dir outputs/decoupling
```

The Qabs calculation creates `<Material>_logR700_Qabs.csv` and an integrals CSV for each material. Omit `--no_plot` to create the absorption-efficiency plots. The decoupling script creates `decoupling_full.pdf` using the historical settings described in the recovery notes.

```sh
python -m unittest discover -s tests -v
```

Checks cover the Rayleigh limit, the nonabsorbing-sphere limit, Planck-mean normalization, consistent units between the two scripts, the grey-particle equilibrium limit, and missing-input failure reporting. Passing them does not establish agreement with the final paper.

Validation on 2026-10-28: all seven checks passed under Python 3.12.14 with the versions pinned in `requirements.txt`. Both executable scripts also passed Python syntax checks.

## Citation and credits

Please cite the paper when using the scientific work; `CITATION.cff` contains its bibliographic information. The paper credits both authors with developing the numerical code. The optical-constant compilation and underlying laboratory measurements must also be cited when those data are used.

No redistribution license has been selected for the recovered source code. The article's publication license does not automatically determine the license of these separate source files.
