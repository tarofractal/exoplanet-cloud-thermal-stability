# Recovery status

This repository was assembled on 2026-10-28 from two preserved research-code uploads. It is a recovery checkpoint, not a declaration that all final publication code has been recovered.

## Preserved source

| Archived file | Original upload date (UTC) | SHA-256 |
| --- | --- | --- |
| `archive/qabs_loggrid_20260428.txt` | 2026-04-28 | `be68556c2581a028b255d6a921140859ac5422ef661ebae584d8607c996c9396` |
| `archive/thermal_decoupling_20260406.py` | 2026-04-06 | `2e36af55d01747eac31c4d40897a5fdccd85febc7c0a2375ddadc08d7d9d5abc` |

The archived Qabs source contains a notebook-only `!pip` command before a `from __future__` import and uses PyMieScatt. It is retained as a text file, exactly as recovered. The archived decoupling source is retained as Python, also without edits.

## Packaging changes

- The executable Qabs copy removes the notebook installation command. Dependencies are installed through `requirements.txt`.
- Its Mie backend is migrated to `miepython.efficiencies`. Diameter and wavelength are both passed in micrometres; absorption efficiency is calculated as `Qext - Qsca`. The absorption sign convention is explicit.
- Its command line returns a nonzero status if any material fails, instead of reporting success after missing inputs.
- The executable decoupling copy accepts explicit input/output directories. It does not recursively search notebook or unrelated filesystem directories for input CSVs.
- The recovered scientific grids, interpolation, energy-balance expressions, and plotting parameters are otherwise retained. None of the missing final figure scripts has been invented or described as recovered.

## Known differences from the published paper

| Topic | Recovered source | Published paper |
| --- | --- | --- |
| Mie backend | PyMieScatt in the preserved source; miepython in the runnable copy | miepython, section 2.2 |
| UV range for Mg2SiO4 and Al2O3 | The grid clips to the measured data range | Constant n and k are extended from 0.2 down to 0.1 micrometres |
| Decoupling stellar temperature | 5000 K | Figure 1 uses 6000 K |
| Second decoupling panel | Decoupling pressure versus particle radius | Figure 1 shows particle temperature versus pressure |

Consequently, running these scripts unchanged will not recreate the final Figure 1 or establish reproduction of the paper. The short-wavelength treatment also requires reconciliation with the actual final input files before scientific use.

## Still missing

- The exact optical-constant files used in the study and their redistribution terms.
- The corresponding Qabs CSV products.
- The final revised decoupling and blackbody-overlay plotting scripts.
- The full particle-temperature, condensation/critical-temperature, and WASP-17b vertical-profile scripts and their inputs.
- Final figure-output comparisons against the publication.

The paper supplies scientific methods and constraints, but it does not contain the original Python source. Numerical checks here use explicitly synthetic inputs and do not replace recovery of these missing materials.
