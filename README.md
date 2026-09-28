# CSTR pH controller benchmark

Compare controller compensation of the nonlinear local process gain in a lab-scale strong-acid/strong-base CSTR.

## Open the notebook

[Main notebook](Lab_Scale_CSTR_Simulation_2_corrected.ipynb) · [Independent schedule experiment](CSTR_Schedule_Cause_Separation.ipynb)

The main notebook contains embedded figures and analysis helpers, so uploading just the notebook to Google Colab does not require a local `robustness` package. In Colab, use **File → Upload notebook**. To run the companion command-line analyses, download the complete repository.

For local execution with Python 3.11 or 3.12:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m jupyterlab
```

Run the main notebook in order through Section 6, then run the **User-controller benchmark** near the end. Edit `USER_GAIN` to provide a scalar, gain scheduler, or callable local gain. The default is explicitly a copy of C4b, not a new proposed controller. Full tuning sweeps and dynamic simulations can take substantial time.

## What the benchmark measures

The candidate is compared with the three existing controllers having the lowest proportional loop-gain max/min ratio on 20,001 pH points between 5 and 9. Absolute and normalized curves separate tuning strength from compensation shape. A ±10% band is editable and is not a stability criterion. Static compensation does not establish tracking performance or robustness. Existing numerical validation did not establish convergence or robustness over the tested uncertainty set.

## Included evidence

- `output/`: saved summary images (PNG/PDF), numerical tables, protocols and reports.
- `output/original_schedule_ablation/design_manifest.csv`: frozen controller schedules.
- `robustness/`: scripts to reproduce compensation figures and numerical validation.
- `images/notebook/`: exported images from saved notebook cell outputs, also visible inside the notebooks.
- `ASSETS.md`: image inventory with previews and links.
- `MANIFEST.sha256`: checksums for packaged files (excluding itself and Git metadata).

Saved outputs and embedded figure snapshots are historical evidence; this packaging update does not claim to rerun the full simulation study. Executing analysis scripts can overwrite saved outputs. The companion notebook uses its own simplified simulator, while the robustness scripts use the main notebook's simulator. See [validation methodology](robustness/README.md).

The distribution contains notebook-related code and results. No reuse license has been assigned; the owner can add one before broader distribution.
