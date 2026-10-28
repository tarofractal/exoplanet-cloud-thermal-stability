# Scientific inputs

The original numerical input files have not been recovered. No physical dataset is bundled or fabricated.

## Optical constants

Put these files in `data/optical_constants/`, or provide a directory with `--data_dir`:

| Material | Expected filename |
| --- | --- |
| Al2O3 | `Al2O3.dat` |
| Fe | `Fe.dat` |
| KCl | `KCl.dat` |
| Mg2SiO4 | `Mg2SiO4_amorph_sol-gel.dat` |
| MnS | `MnS.dat` |
| Na2S | `Na2S.dat` |
| SiO2 | `SiO2_alpha.dat` |
| ZnS | `ZnS.dat` |

Each file is expected to contain at least three numeric columns: wavelength in micrometres, real refractive index n, and absorption index k. Lines beginning with `#` are comments.

The paper cites the compilation of Kitzmann & Heng (2018), *Optical properties of potential condensates in exoplanetary atmospheres*, MNRAS, 475, 94, https://doi.org/10.1093/mnras/stx3141. The exact files, processing history, original laboratory citations, and redistribution permission still need to be matched to the recovered code.

The published UV extension for Mg2SiO4 and Al2O3 is not implemented in the recovered snapshot. See `docs/RECOVERY_STATUS.md`.

## Qabs tables

The decoupling code expects the eight `<Material>_logR700_Qabs.csv` files in one directory. Each table has `Wavelength (μm)` and radius-dependent columns such as `Qabs_0.1μm`. The Qabs script produces that format.

Tests generate temporary, synthetic optical constants solely for numerical checks. Their results are not observational data or the inputs used in the publication.
