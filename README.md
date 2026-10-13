# CSTR pH controller benchmark

Compare controller compensation of the nonlinear local process gain in a lab-scale strong-acid/strong-base CSTR.

## See the loop-gain flatness versus operating pH graph

**[Open the loop-gain-versus-operating-pH graph](output/compensation_vs_tracking/loop_gain_vs_pH.png)** — view the saved figure directly on GitHub; no notebook execution is required.

![Proportional loop gain versus operating pH, with flatness F for each schedule](output/compensation_vs_tracking/loop_gain_vs_pH.png)

The horizontal axis is operating pH; the curves show local proportional loop gain, `L_P(pH) = Kc(pH) Kp,local(pH)`. Each legend reports the flatness ratio `F = max(L_P)/min(L_P)` over pH 5–9. A flatter curve indicates more complete compensation of process-gain variation; ideal compensation has `F = 1`. A curve normalized to its reference value is unity everywhere under ideal compensation. The absolute loop gain need not equal one.

In the **[main notebook](Lab_Scale_CSTR_Simulation_2_corrected.ipynb)**, look for **Section 6: Central result: loop-gain flatness and nonlinear gain compensation** for the executable analysis. The **Packaged results and figure gallery** near the end displays the saved graph under **Loop-gain flatness versus operating pH — saved graph**. The [Section 6 summary figure](output/loop_gain_flatness/loop_gain_flatness.png) also shows process gain, controller gain and proportional loop gain.

These saved figures quantify static compensation. The notebook's **Nonlinear compensation versus setpoint tracking** analysis relates it to tracking error under the declared model and tuning conditions.

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
