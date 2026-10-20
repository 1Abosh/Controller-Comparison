"""Render the six selected controllers on the notebook's standardised scenario.

Run from the repository root: python controller_responses.py
Requires the packages in requirements.txt. No parameter search is performed.
"""
from pathlib import Path
import ast
import contextlib
import io
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "robustness"))
from run_validation import NS, CELLS, DESIGNS, DEV, pars


def main():
    definitions = ast.Module(
        body=[node for node in ast.parse("".join(CELLS[20]["source"])).body
              if isinstance(node, ast.FunctionDef)], type_ignores=[])
    exec(compile(definitions, "notebook_imc_tuning", "exec"), NS)
    kp, tau, theta = np.mean(pars, axis=0)
    kc, ti, _ = NS["imc_pid_rivera"](kp, tau, theta, 5 * theta)
    controllers = [
        ("C1", "Conventional PI (fixed)", DESIGNS["fixed/C1"]),
        ("C2", "GS-PI (3 references)", DESIGNS["reference_PI/A_original3"]),
        ("C2b", "GS-PI (5 references)", DESIGNS["reference_PI/E_original5_retuned"]),
        ("C3", "IMC PI (fixed)", dict(Kc=kc, Ti=ti, Td=0.0)),
        ("C4", "GS-IMC (3 references)", DESIGNS["IMC_PID/A_original3"]),
        ("C4b", "GS-IMC (5 references)", DESIGNS["IMC_PID/E_original5_retuned"]),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True, sharey=True)
    for ax, (cid, label, params) in zip(axes.flat, controllers):
        # Use the original simulator and its original sampling/delay conventions.
        with contextlib.redirect_stdout(io.StringIO()):
            result = NS["run_closed_loop_sim"](
                **params, **DEV["standardised"], plot=False)
        assert np.isfinite(result["pH"]).all(), cid
        ax.plot(result["t"] / 1000, result["SP"], "k--", lw=1.3, label="Setpoint")
        ax.plot(result["t"] / 1000, result["pH"], color="#146c94", lw=1.2,
                label="Reactor pH")
        for index, event in enumerate(DEV["standardised"]["dist_events"]):
            ax.axvline(event["time"] / 1000, color="#b66b23", ls=":", lw=1,
                       label="Acid-feed change" if index == 0 else None)
        ax.set_title(f"{cid} — {label}", fontsize=12)
        ax.set_xlim(0, 170)
        ax.grid(alpha=0.2)
        ax.legend(loc="upper right", fontsize=8)
        print(f"{cid}: {len(result['t'])} samples, finite response", flush=True)
    for ax in axes[:, 0]:
        ax.set_ylabel("pH")
    for ax in axes[-1, :]:
        ax.set_xlabel("Time (ks)")
    fig.suptitle("Six-controller comparison: setpoint tracking and disturbances", fontsize=16)
    fig.text(0.5, 0.015,
             "Setpoint: 7 → 8 → 7 → 6   |   Acid-feed factors: 1.05, 0.95, 1.05 at 20, 65, 120 ks",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    fig.savefig(ROOT / "controller_responses.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
