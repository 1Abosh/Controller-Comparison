# Model and controller code extracted verbatim from the supplied local notebook.
# Cell numbers are zero-based. Plotting/export code is in regenerate_figures.py.

# ORIGINAL NOTEBOOK CELL 2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ============================================================
# PROCESS CONSTANTS (lab-scale CSTR)
# ============================================================
Qa = 0.005 * 1000   # mL/s, acid volumetric flow
Ca = 0.01            # mol/L, acid concentration
Cb = 0.01            # mol/L, base concentration
Kw = 1e-14           # water equilibrium constant
V  = 10 * 1000       # mL, reactor volume

print("--- System Inputs ---")
print(f"Acid Volumetric Flow (Qa): {Qa} mL/s")
print(f"Acid Concentration (Ca): {Ca} mol/L")
print(f"Base Concentration (Cb): {Cb} mol/L")
print(f"Water Equilibrium Constant (Kw): {Kw}")
print(f"Reactor Volume (V): {V} mL")
print("---------------------")

# ORIGINAL NOTEBOOK CELL 3
# ============================================================
# CHARGE-BALANCE <-> pH CONVERSIONS
# ============================================================
# x = [H+] - [OH-] is the conserved reaction-invariant charge balance used as the
# state variable throughout (Kalafatis et al., 1992; McAvoy et al.).

def x_to_pH(x):
    """Convert charge-balance state x to pH (vectorised)."""
    H = (x + np.sqrt(x**2 + 4 * Kw)) / 2
    return -np.log10(H)

def pH_to_x(pH):
    """Convert a target pH to the corresponding charge-balance value x."""
    H = 10.0 ** (-pH)
    return H - Kw / H

def residual_Qb(Qb, Qa_val, Ca_val, Cb_val, x_target):
    """Steady-state charge balance residual, used with fsolve to find the base
    flow Qb that holds the reactor at a target x (equivalently, a target pH)."""
    return (Qa_val * Ca_val - Qb * Cb_val) / (Qa_val + Qb) - x_target

def solve_Qb_for_pH(pH_target, guess, Ca_val=None):
    """Steady-state base flow Qb0 required to hold the reactor at pH_target."""
    Ca_val = Ca if Ca_val is None else Ca_val
    x_target = pH_to_x(pH_target)
    return fsolve(residual_Qb, guess, args=(Qa, Ca_val, Cb, x_target))[0]

def default_F_max(pH_ref=9.0):
    """A generous upper bound on the base flow actuator, computed from the physics
    (the flow that would drive the reactor to pH_ref) rather than a magic number."""
    x_ref = pH_to_x(pH_ref)
    Qb_ref = fsolve(residual_Qb, 1e-3, args=(Qa, Ca, Cb, x_ref))[0]
    return max(2 * Qb_ref, 0.1)

# ORIGINAL NOTEBOOK CELL 4
# ============================================================
# STEP-TEST ROUTINE (process identification)
# ============================================================
# Runs a two-stage step test around a target operating pH: first bring the reactor from
# acid-only to the target pH, then apply a small step in Qb and fit a first-order-plus-
# dead-time (FOPDT) model to the response (Kp, tau, theta) using the 63.2%/28.3% method.

def _time_to_percent_change(pc, initial_val, final_val, time_arr, val_arr, step_t):
    """Time after step_t at which val_arr reaches fraction pc of its total change."""
    total_change = final_val - initial_val
    if abs(total_change) < 1e-6:
        return np.nan

    target = initial_val + pc * total_change
    post = time_arr >= step_t
    t_resp, v_resp = time_arr[post], val_arr[post]
    if len(t_resp) == 0:
        return np.nan

    cross = np.where(v_resp >= target)[0] if total_change > 0 else np.where(v_resp <= target)[0]
    if len(cross) == 0:
        return np.nan

    idx = cross[0]
    if idx == 0:
        return t_resp[idx] - step_t

    t1, v1 = t_resp[idx - 1], v_resp[idx - 1]
    t2, v2 = t_resp[idx], v_resp[idx]
    if v2 == v1:
        return t2 - step_t
    return (t1 + (target - v1) * (t2 - t1) / (v2 - v1)) - step_t


def run_step_test(pH_target, step_pct, Qb0_guess, region_label="",
                   t_step=2000, t_step2=20000, tend=30000, n=600, plot=True):
    """
    Identify FOPDT parameters (Kp, tau, theta) around pH_target.

    pH_target : operating point to identify around
    step_pct  : fractional step applied to Qb0 at t_step2 (e.g. 0.0125 for +1.25%)
    Qb0_guess : initial guess for fsolve when finding Qb0
    """
    Qb0 = solve_Qb_for_pH(pH_target, Qb0_guess)
    Qb_new = (1 + step_pct) * Qb0

    def dxdt(t, x):
        if t < t_step:
            Qb = 0.0
        elif t < t_step2:
            Qb = Qb0
        else:
            Qb = Qb_new
        Fout = Qa + Qb
        dx = (Qa * Ca / V) + (Qb * (-Cb) / V) - (Fout / V) * x[0]
        return [dx]

    t_eval = np.linspace(0, tend, n)
    sol = solve_ivp(dxdt, (0, tend), [Ca], t_eval=t_eval)
    pH = x_to_pH(sol.y[0])

    # Fit FOPDT parameters from the response to the second step
    end_settling = t_step2 + (tend - t_step2) * 0.75
    idx_ss = sol.t >= end_settling
    pH_ss = np.mean(pH[idx_ss]) if np.any(idx_ss) else pH[-1]
    delta_Y = pH_ss - pH_target

    t63 = _time_to_percent_change(0.632, pH_target, pH_ss, sol.t, pH, t_step2)
    t28 = _time_to_percent_change(0.283, pH_target, pH_ss, sol.t, pH, t_step2)

    Kp = delta_Y / (Qb_new - Qb0)
    tau = 1.5 * (t63 - t28)
    theta = abs(t63 - tau)

    if plot:
        fig, ax1 = plt.subplots(figsize=(11, 4.5))
        ax1.plot(sol.t, pH, color='tab:blue', label='pH')
        ax1.axvline(t_step, color='tab:red', ls='--', alpha=0.7, label=f'Base flow on (t={t_step}s)')
        ax1.axvline(t_step2, color='tab:orange', ls='--', alpha=0.7,
                     label=f'Qb step {step_pct*100:+.3f}% (t={t_step2}s)')
        ax1.axhline(pH_target, color='purple', ls=':', alpha=0.7, label=f'Target pH = {pH_target}')
        ax1.set_xlabel("Time (s)"); ax1.set_ylabel("pH")
        ax1.grid(alpha=0.3)
        ax1.set_title(f"Step test — {region_label} region" if region_label else "Step test")
        ax1.legend(loc='best', fontsize=9)
        plt.tight_layout(); plt.show()

    print(f"--- {region_label} region (around pH {pH_target}) ---")
    print(f"Qb0 = {Qb0:.5f} mL/s  ->  Qb_new = {Qb_new:.5f} mL/s  ({step_pct*100:+.3f}%)")
    print(f"Process gain Kp = {Kp:.4f} pH/(mL/s)")
    print(f"Time constant tau = {tau:.2f} s")
    print(f"Dead time theta = {theta:.2f} s\n")

    return Kp, tau, theta

# ORIGINAL NOTEBOOK CELL 5
# ============================================================
# CLOSED-LOOP PERFORMANCE METRICS
# ============================================================

PH_SETTLING_BAND = 0.05   # settling band in absolute pH units (not a % of SP -- pH is
                           # logarithmic with an arbitrary origin, so a %-of-SP band gives a
                           # tighter criterion for acidic setpoints than neutral ones purely
                           # as an artefact of the setpoint's numeric value)

def calculate_performance_metrics(t_eval, SP_trace, pH_arr):
    """IAE / ISE / ITAE / overshoot / settling time (+/- PH_SETTLING_BAND pH) for every
    distinct, real setpoint change in SP_trace. Segments where the setpoint does not
    change (SP_old == SP_new) are skipped: they are not a setpoint step and would
    otherwise dilute the per-controller group means downstream."""
    dt = t_eval[1] - t_eval[0]
    sp_change_indices = [0] + list(np.where(np.diff(SP_trace) != 0)[0] + 1)

    results = []
    step_num = 0
    for i in range(len(sp_change_indices)):
        start_idx = sp_change_indices[i]
        end_idx = sp_change_indices[i + 1] if (i + 1) < len(sp_change_indices) else len(t_eval)

        SP_old = SP_trace[start_idx - 1] if start_idx > 0 else SP_trace[0]
        SP_new = SP_trace[start_idx]

        if SP_old == SP_new:
            continue

        if end_idx > start_idx:
            t_slice = t_eval[start_idx:end_idx]
            SP_slice = SP_trace[start_idx:end_idx]
            pH_slice = pH_arr[start_idx:end_idx]
        else:
            t_slice = t_eval[start_idx:start_idx + 1]
            SP_slice = SP_trace[start_idx:start_idx + 1]
            pH_slice = pH_arr[start_idx:start_idx + 1]
        if len(t_slice) == 0:
            continue

        error_slice = SP_slice - pH_slice
        iae = np.sum(np.abs(error_slice)) * dt
        ise = np.sum(error_slice**2) * dt
        t_shifted = t_slice - t_slice[0]
        itae = np.sum(t_shifted * np.abs(error_slice)) * dt

        overshoot = 0.0
        if SP_new > SP_old:
            dev = pH_slice - SP_new
            dev = dev[dev > 0]
            if len(dev) > 0:
                overshoot = np.max(dev)
        elif SP_new < SP_old:
            dev = SP_new - pH_slice
            dev = dev[dev > 0]
            if len(dev) > 0:
                overshoot = np.max(dev)

        settling_time = np.inf
        tol = PH_SETTLING_BAND
        upper, lower = SP_new + tol, SP_new - tol
        out_of_band = (pH_slice < lower) | (pH_slice > upper)
        out_idx = np.where(out_of_band)[0]
        if len(out_idx) == 0:
            settling_time = 0.0
        else:
            last = out_idx[-1]
            if last == len(t_slice) - 1:
                settling_time = np.inf
            elif np.all(~out_of_band[last + 1:]):
                settling_time = t_slice[last + 1] - t_slice[0]

        step_num += 1
        results.append({
            'Step': step_num, 'Time': t_slice[0],
            'SP_Transition': f"{SP_old:.2f} -> {SP_new:.2f}",
            'IAE': iae, 'ISE': ise, 'ITAE': itae,
            'Overshoot': overshoot, 'Settling_Time': settling_time,
        })
    return results


def metrics_table(metrics_list):
    """Tidy pandas table of a metrics_list, ready to paste into the report."""
    df = pd.DataFrame(metrics_list)
    if df.empty:
        return df
    df = df.rename(columns={
        'Step': 'Step', 'Time': 'Time (s)', 'SP_Transition': 'SP transition',
        'IAE': 'IAE', 'ISE': 'ISE', 'ITAE': 'ITAE',
        'Overshoot': 'Overshoot (pH)', 'Settling_Time': 'Settling time (s)',
    })
    return df.round(3)


# ORIGINAL NOTEBOOK CELL 6
# ============================================================
# PIECEWISE-LINEAR GAIN SCHEDULER
# ============================================================

class PiecewiseLinearGainController:
    """Interpolates a scheduled gain (Kc, Ti or Td) linearly between (PV, gain)
    breakpoints. Used both for the IMC-scheduled controller (few, region-based
    breakpoints) and the pure gain-scheduled controller (many breakpoints tracing
    out a continuous curve)."""

    def __init__(self, breakpoints):
        self.breakpoints = sorted(breakpoints, key=lambda p: p[0])
        if len(self.breakpoints) < 2:
            raise ValueError("At least two breakpoints are required for interpolation.")

    def calculate_gain(self, current_pv):
        bp = self.breakpoints
        if current_pv <= bp[0][0]:
            return bp[0][1]
        if current_pv >= bp[-1][0]:
            return bp[-1][1]
        for (x1, y1), (x2, y2) in zip(bp[:-1], bp[1:]):
            if x1 <= current_pv <= x2:
                if x2 == x1:
                    return y1
                return y1 + (current_pv - x1) * (y2 - y1) / (x2 - x1)


def print_gain_schedule_details(name, breakpoints):
    print(f"--- {name} schedule ({len(breakpoints)} breakpoints) ---")
    bp = sorted(breakpoints, key=lambda p: p[0])
    for i in range(min(len(bp) - 1, 8)):
        x1, y1 = bp[i]; x2, y2 = bp[i + 1]
        if x2 == x1:
            print(f"  pH = {x1:.2f}: constant gain = {y1:.4f}")
            continue
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
        print(f"  pH {x1:.2f} to {x2:.2f}: gain = {m:.4f}*pH + {b:.4f}")
    if len(bp) - 1 > 8:
        print(f"  ... ({len(bp) - 1 - 8} more segments)")

# ORIGINAL NOTEBOOK CELL 7
# ============================================================
# CLOSED-LOOP SIMULATOR (used by every controller in this notebook)
# ============================================================
# Kc, Ti, Td may each be a constant (fixed-gain controller) or a
# PiecewiseLinearGainController (gain-scheduled on measured pH) -- the same function
# handles the conventional, IMC, IMC-scheduled and pure gain-scheduled controllers.

def run_closed_loop_sim(Kc, Ti, Td=0.0, sp_events=None, dist_events=None,
                         initial_pH=None, theta_dead=None,
                         F_min=0.0, F_max=None, tend=100000, n_points=None,
                         schedule_filter_tau=0.0,
                         title="Closed-Loop pH Control", show_mv=False, plot=True):
    """
    schedule_filter_tau : if > 0, the pH used to *look up* scheduled gains is a
        first-order-filtered version of the measurement (a standard practical fix for
        gain-scheduling on a fast/noisy PV: filter it before indexing the schedule, don't
        let the raw measurement -- and hence the schedule -- chatter). The unfiltered
        measurement is still used for the control error itself.
    """
    if not sp_events:
        raise ValueError("sp_events must be a non-empty list, e.g. [{'time': 0, 'SP': 7.0}]")
    dist_events = dist_events or []
    theta_dead = theta if theta_dead is None else theta_dead
    F_max = default_F_max() if F_max is None else F_max

    # Keep the controller sample time dt well below the fastest dynamics anywhere in the
    # operating range (down to ~18-40 s near equivalence). A coarse dt here doesn't just
    # lose accuracy -- for the more aggressive gain-scheduled controllers it can produce a
    # simulation artifact that looks like instability but is really just under-sampling.
    if n_points is None:
        n_points = max(800, int(tend / 20))

    t_eval = np.linspace(0, tend, n_points)
    n = len(t_eval)
    dt = t_eval[1] - t_eval[0]
    theta_steps = max(1, int(theta_dead / dt))

    def sp_at(t):
        cur = sp_events[0]['SP']
        for ev in sp_events:
            if t >= ev['time']:
                cur = ev['SP']
        return cur

    def ca_at(t):
        cur = Ca
        for ev in dist_events:
            if t >= ev['time']:
                cur = Ca * ev['Ca_factor']
        return cur

    def scheduled(val, pH):
        return val.calculate_gain(pH) if isinstance(val, PiecewiseLinearGainController) else val

    initial_pH = sp_events[0]['SP'] if initial_pH is None else initial_pH
    x0_val = pH_to_x(initial_pH)

    F = np.zeros(n); pH_arr = np.zeros(n); e = np.zeros(n)
    pH_arr[0] = initial_pH
    e[0] = sp_at(0) - pH_arr[0]
    F[0] = solve_Qb_for_pH(initial_pH, guess=0.0)
    x0 = [x0_val]
    MV_buffer = np.full(theta_steps, F[0])

    t_span_start = 0.0
    pH_sched_filt = pH_arr[0]
    for k in range(1, n):
        current_SP = sp_at(t_eval[k])
        if schedule_filter_tau > 0:
            alpha = dt / (schedule_filter_tau + dt)
            pH_sched_filt = pH_sched_filt + alpha * (pH_arr[k - 1] - pH_sched_filt)
            current_pH = pH_sched_filt
        else:
            current_pH = pH_arr[k - 1]
        current_Kc = scheduled(Kc, current_pH)
        current_Ti = scheduled(Ti, current_pH)
        current_Td = scheduled(Td, current_pH)

        e[k] = current_SP - pH_arr[k - 1]
        de_val = e[k] - e[k - 1]
        pv_prev1 = pH_arr[k - 2] if k >= 2 else pH_arr[0]
        pv_prev2 = pH_arr[k - 3] if k >= 3 else pH_arr[0]
        de_val2 = pH_arr[k - 1] - 2 * pv_prev1 + pv_prev2

        ti_term = (dt / current_Ti) * e[k] if current_Ti > 1e-6 else 0.0
        td_term = (current_Td / dt) * de_val2 if dt > 1e-6 else 0.0

        dF = current_Kc * (de_val + ti_term - td_term)
        F[k] = np.clip(F[k - 1] + dF, F_min, F_max)

        Qb_delayed = MV_buffer[-1]
        MV_buffer = np.roll(MV_buffer, 1)
        MV_buffer[0] = F[k]

        current_Ca = ca_at(t_eval[k])

        def dxdt_func(t, x_vec, Qb_in, Ca_val):
            Fout = Qa + Qb_in
            dx = (Qa * Ca_val / V) + (Qb_in * (-Cb) / V) - (Fout / V) * x_vec[0]
            return [dx]

        sol = solve_ivp(dxdt_func, (t_span_start, t_span_start + dt), x0,
                         args=(Qb_delayed, current_Ca), method='RK45')
        x0 = [sol.y[0, -1]]
        t_span_start += dt
        pH_arr[k] = x_to_pH(x0[0])

    SP_trace = np.array([sp_at(t) for t in t_eval])
    metrics_list = calculate_performance_metrics(t_eval, SP_trace, pH_arr)

    if plot:
        fig, ax1 = plt.subplots(figsize=(11, 4.8))
        ax1.plot(t_eval, pH_arr, color='tab:blue', label='pH (PV)')
        ax1.plot(t_eval, SP_trace, color='tab:red', ls='--', label='Setpoint (SP)')
        for dv in dist_events:
            label = dv.get('label', f"Ca x{dv['Ca_factor']}")
            ax1.axvline(dv['time'], color='tab:green', ls=':', alpha=0.7, label=f"{label} (t={dv['time']}s)")
        ax1.set_xlabel("Time (s)"); ax1.set_ylabel("pH")
        ax1.grid(alpha=0.3)
        handles, labels = ax1.get_legend_handles_labels()
        uniq = dict(zip(labels, handles))
        ax1.legend(uniq.values(), uniq.keys(), loc='lower right', fontsize=8)
        ax1.set_title(title)
        plt.tight_layout(); plt.show()

        if show_mv:
            fig2, ax2 = plt.subplots(figsize=(11, 3.2))
            ax2.step(t_eval, F, color='tab:green', where='post', label='Base flow (MV)')
            ax2.set_xlabel("Time (s)"); ax2.set_ylabel("Flow (mL/s)")
            ax2.grid(alpha=0.3); ax2.legend(loc='best', fontsize=8)
            plt.tight_layout(); plt.show()

    df = metrics_table(metrics_list)
    print(f"--- {title}: performance per setpoint step ---")
    if not df.empty:
        try:
            from IPython.display import display
            display(df)
        except Exception:
            print(df.to_string(index=False))

    return {'t': t_eval, 'pH': pH_arr, 'F': F, 'SP': SP_trace,
            'metrics': metrics_list, 'metrics_df': df}

# ORIGINAL NOTEBOOK CELL 11
Kp_ac, tau_ac, theta_ac = run_step_test(pH_target=4, step_pct=0.0125, Qb0_guess=0.05,
                                         region_label="Acidic")

# ORIGINAL NOTEBOOK CELL 12
Kp_neu, tau_neu, theta_neu = run_step_test(pH_target=6, step_pct=-0.0010, Qb0_guess=5.0,
                                            region_label="Near-neutral")

# ORIGINAL NOTEBOOK CELL 13
Kp_alk, tau_alk, theta_alk = run_step_test(pH_target=9, step_pct=0.0500, Qb0_guess=5.0,
                                            region_label="Alkaline")

# ORIGINAL NOTEBOOK CELL 15
Kp = (Kp_ac + Kp_neu + Kp_alk) / 3
tau = (tau_ac + tau_neu + tau_alk) / 3
theta = (theta_ac + theta_neu + theta_alk) / 3

print("--- Individual region parameters ---")
print(f"Acidic:       Kp={Kp_ac:.2f}, tau={tau_ac:.2f} s, theta={theta_ac:.2f} s")
print(f"Near-neutral: Kp={Kp_neu:.2f}, tau={tau_neu:.2f} s, theta={theta_neu:.2f} s")
print(f"Alkaline:     Kp={Kp_alk:.2f}, tau={tau_alk:.2f} s, theta={theta_alk:.2f} s")
print(f"\n--- Averaged process parameters ---")
print(f"Kp = {Kp:.2f}, tau = {tau:.2f} s, theta = {theta:.2f} s")

# Ziegler-Nichols open-loop tuning (Table 2, Modelling and Designing Controllers for pH Process)
Kc = 0.9 * tau / (Kp * theta)
tau_i = 3.33 * theta
print(f"\nZiegler-Nichols PI tuning: Kc = {Kc:.6f}, Ti = {tau_i:.2f} s")

# ORIGINAL NOTEBOOK CELL 20
# Single source of truth for the Rivera/Morari/Skogestad IMC-PID formula, used by both
# the fixed IMC-PI tuning directly below and the regional IMC tuning further down, so the
# formula cannot drift between the two call sites again.
def imc_pid_rivera(Kp, tau, theta, lam):
    """Conventional IMC-PID (Rivera, Morari & Skogestad, 1986).
    First-order Pade approximation of the dead time. Returns (Kc, Ti, Td)."""
    Kc = (2 * tau + theta) / (Kp * (2 * lam + theta))
    Ti = tau + theta / 2
    Td = tau * theta / (2 * tau + theta)
    return Kc, Ti, Td


# ORIGINAL NOTEBOOK CELL 26
def calculate_bode_data(Kp_val, tau_val, theta_val, omega_values):
    """G(s) = Kp / (tau*s + 1) * Pade(1, theta) approximation of dead time."""
    mag_term1 = Kp_val / np.sqrt(1 + (tau_val * omega_values) ** 2)
    phase_term1 = -np.arctan(tau_val * omega_values)
    phase_term2 = -2 * np.arctan(theta_val * omega_values / 2)
    return mag_term1, phase_term1 + phase_term2   # magnitude of the Pade term is 1


def find_ultimate_gain_period(Kp_val, tau_val, theta_val, omega_values):
    """Ku, Pu at the frequency where the open-loop phase crosses -180 degrees."""
    theta_val = max(1e-3, theta_val)
    magnitude, phase = calculate_bode_data(Kp_val, tau_val, theta_val, omega_values)
    crossings = np.where(np.diff(np.sign(phase + np.pi)))[0]
    if len(crossings) == 0:
        return np.nan, np.nan, magnitude, phase
    idx = crossings[0]
    w1, p1 = omega_values[idx], phase[idx]
    w2, p2 = omega_values[idx + 1], phase[idx + 1]
    omega_u = w1 + (w2 - w1) * (-np.pi - p1) / (p2 - p1)
    mag_u = np.interp(omega_u, omega_values, magnitude)
    Ku = 1 / mag_u
    Pu = 2 * np.pi / omega_u
    return Ku, Pu, magnitude, phase


omega = np.logspace(-3, 1, 500)
Ku, Pu, magnitude_avg, phase_avg = find_ultimate_gain_period(Kp, tau, theta, omega)
print(f"Averaged process: Ku = {Ku:.4f}, Pu = {Pu:.2f} s")

# ORIGINAL NOTEBOOK CELL 27
regions = {'Acidic': (Kp_ac, tau_ac, theta_ac), 'Near-neutral': (Kp_neu, tau_neu, theta_neu),
           'Alkaline': (Kp_alk, tau_alk, theta_alk)}

bode_results = {}
for name, (Kp_r, tau_r, theta_r) in regions.items():
    bode_results[name] = find_ultimate_gain_period(Kp_r, tau_r, theta_r, omega)

Ku_ac_region, Pu_ac_region = bode_results['Acidic'][:2]
Ku_neu_region, Pu_neu_region = bode_results['Near-neutral'][:2]
Ku_alk_region, Pu_alk_region = bode_results['Alkaline'][:2]

fig, (ax_mag, ax_ph) = plt.subplots(2, 1, figsize=(12, 8))
colors = {'Acidic': 'tab:blue', 'Near-neutral': 'tab:green', 'Alkaline': 'tab:red'}
for name, (Kp_r, tau_r, theta_r) in regions.items():
    mag, ph = calculate_bode_data(Kp_r, tau_r, max(1e-3, theta_r), omega)
    ax_mag.semilogx(omega, 20 * np.log10(mag), label=f'{name} region', color=colors[name])
    ax_ph.semilogx(omega, ph * 180 / np.pi, label=f'{name} region', color=colors[name])

ax_mag.set_title('Bode Plot for pH Process in Different Regions')
ax_mag.set_ylabel('Magnitude (dB)'); ax_mag.grid(True, which='both'); ax_mag.legend()
ax_ph.set_xlabel('Frequency (rad/s)'); ax_ph.set_ylabel('Phase (degrees)')
ax_ph.grid(True, which='both'); ax_ph.legend()
plt.tight_layout(); plt.show()

print("--- Ultimate gain (Ku) and ultimate period (Pu) per region ---")
for name in regions:
    Ku_r, Pu_r = bode_results[name][:2]
    print(f"{name}: Ku = {Ku_r:.4f}, Pu = {Pu_r:.2f} s")

# ORIGINAL NOTEBOOK CELL 48
# ============================================================
# ANALYTIC LOCAL PROCESS GAIN dpH/dQb
# ============================================================
# Differentiate the steady-state charge-balance relation used everywhere else in this
# notebook, x_ss(Qb) = (Qa*Ca - Qb*Cb)/(Qa+Qb), directly, instead of taking a finite
# difference on a fixed grid that is too coarse near equivalence.
#
#   dx/dQb = -Qa*(Ca+Cb) / (Qa+Qb)^2                          (quotient rule on x_ss(Qb))
#   dpH/dx = -(1/(H*ln(10))) * (1 + x/sqrt(x^2+4*Kw)) / 2      (chain rule through x_to_pH)
#
# dpH/dQb = dpH/dx * dx/dQb is then exact, in closed form, for any Qb -- no grid at all.

def dx_dQb_analytic(Qb):
    return -Qa * (Ca + Cb) / (Qa + Qb) ** 2

def dpH_dx_analytic(x):
    H = (x + np.sqrt(x ** 2 + 4 * Kw)) / 2
    dH_dx = (1 + x / np.sqrt(x ** 2 + 4 * Kw)) / 2
    return -1.0 / (H * np.log(10)) * dH_dx

def local_gain_analytic(pH_query, Qb_guess=5.0):
    """Exact local steady-state gain dpH/dQb at pH_query (analytic differentiation,
    no grid, no np.gradient)."""
    Qb0 = solve_Qb_for_pH(pH_query, Qb_guess)
    x0 = (Qa * Ca - Qb0 * Cb) / (Qa + Qb0)
    return dpH_dx_analytic(x0) * dx_dQb_analytic(Qb0)

# --- Validate against a central difference on a refined grid (independent of the analytic
# derivation above) -- h = 1e-7 mL/s balances truncation error against floating-point
# round-off; both are visible if h is swept away from this value, see the printed table. ---
_CENTRAL_DIFF_H = 1e-7   # mL/s

def _central_diff_gain(pH_query, Qb_guess=5.0):
    Qb0 = solve_Qb_for_pH(pH_query, Qb_guess)
    x_fwd = (Qa * Ca - (Qb0 + _CENTRAL_DIFF_H) * Cb) / (Qa + Qb0 + _CENTRAL_DIFF_H)
    x_bwd = (Qa * Ca - (Qb0 - _CENTRAL_DIFF_H) * Cb) / (Qa + Qb0 - _CENTRAL_DIFF_H)
    return (x_to_pH(x_fwd) - x_to_pH(x_bwd)) / (2 * _CENTRAL_DIFF_H)

print("--- Analytic local gain dpH/dQb vs. central-difference validation ---")
print(f"{'pH':>4s}  {'analytic':>14s}  {'central diff':>14s}  {'rel. diff':>10s}")
local_gain_by_pH = {}
for _pH in [4, 6, 7, 9]:
    g_analytic = local_gain_analytic(_pH)
    g_cd = _central_diff_gain(_pH)
    local_gain_by_pH[_pH] = g_analytic
    rel = abs(g_analytic - g_cd) / abs(g_cd)
    print(f"{_pH:4d}  {g_analytic:14.6f}  {g_cd:14.6f}  {rel:10.2e}")
    if rel > 1e-6:
        raise AssertionError(f"Analytic local gain at pH={_pH} disagrees with central "
                              f"difference beyond 6 significant figures (rel. diff {rel:.2e}).")
print("\nAnalytic derivative validated to >=6 significant figures at all four points.")


# ORIGINAL NOTEBOOK CELL 52
# ============================================================
# DESIGN RULES (Section 2) -- reused identically by all six controllers below
# ============================================================

def pi_reference(Ku_val, Pu_val, f):
    """Ziegler-Nichols-family reference tuning from the ultimate gain/period. Returns
    (Kc, Ti, Td); Td is always 0 -- PI only, matching every other reference-tuned/fixed
    controller in this notebook."""
    return f * Ku_val, 2.0 * Pu_val, 0.0

REFERENCE_TUNING_F_DEFAULT = 0.2
IMC_LAMBDA_K_DEFAULT = 2.25   # == LAMBDA_MULTIPLIER above, named again here so this
                              # section's rule statement is self-contained


# ORIGINAL NOTEBOOK CELL 54
# ============================================================
# SIX CONTROLLERS (Section 3) -- builder functions, parametrised by f or k
# ============================================================
# region_map_5ref: which identified region's (tau, theta) each of the five C2b/C4b reference
# points borrows -- see the assumption above. `regions` (Frequency response section) already
# holds the three identified (Kp, tau, theta) triples; only tau/theta are used here.
REF_PH_5 = [4, 6, 7, 8, 9]
region_map_5ref = {
    4: ('Acidic', tau_ac, theta_ac),
    6: ('Near-neutral', tau_neu, theta_neu),
    7: ('Near-neutral', tau_neu, theta_neu),
    8: ('Near-neutral', tau_neu, theta_neu),
    9: ('Alkaline', tau_alk, theta_alk),
}

def build_C1(f):
    """Conventional PI, reference tuning, fixed at the near-neutral reference point."""
    Kc1, Ti1, Td1 = pi_reference(Ku_neu_region, Pu_neu_region, f)
    return dict(Kc=Kc1, Ti=Ti1, Td=Td1)

def build_C2(f):
    """GS-PI, 3 references (acidic/near-neutral/alkaline identification points)."""
    Kc_ac, Ti_ac, _ = pi_reference(Ku_ac_region, Pu_ac_region, f)
    Kc_neu, Ti_neu, _ = pi_reference(Ku_neu_region, Pu_neu_region, f)
    Kc_alk, Ti_alk, _ = pi_reference(Ku_alk_region, Pu_alk_region, f)
    kc_sched = PiecewiseLinearGainController([(4, Kc_ac), (6, Kc_neu), (9, Kc_alk)])
    ti_sched = PiecewiseLinearGainController([(4, Ti_ac), (6, Ti_neu), (9, Ti_alk)])
    return dict(Kc=kc_sched, Ti=ti_sched, Td=0.0, schedule_filter_tau=schedule_filter_tau)

def build_C3(k):
    """IMC PI, fixed, on the averaged FOPDT model."""
    Kc3, Ti3, _ = imc_pid_rivera(Kp, tau, theta, k * theta)
    return dict(Kc=Kc3, Ti=Ti3, Td=0.0)   # PI only, matching the fixed IMC PI above

def build_C4(k):
    """GS-IMC, 3 references (acidic/near-neutral/alkaline identification points)."""
    Kc_ac, Ti_ac, Td_ac = imc_pid_rivera(Kp_ac, tau_ac, theta_ac, k * theta_ac)
    Kc_neu, Ti_neu, Td_neu = imc_pid_rivera(Kp_neu, tau_neu, theta_neu, k * theta_neu)
    Kc_alk, Ti_alk, Td_alk = imc_pid_rivera(Kp_alk, tau_alk, theta_alk, k * theta_alk)
    kc_sched = PiecewiseLinearGainController([(4, Kc_ac), (6, Kc_neu), (9, Kc_alk)])
    ti_sched = PiecewiseLinearGainController([(4, Ti_ac), (6, Ti_neu), (9, Ti_alk)])
    td_sched = PiecewiseLinearGainController([(4, Td_ac), (6, Td_neu), (9, Td_alk)])
    return dict(Kc=kc_sched, Ti=ti_sched, Td=td_sched, schedule_filter_tau=schedule_filter_tau)

def build_C2b(f):
    """GS-PI, 5 references, on analytic local gains (Section 1) with nearest-region tau/theta."""
    kc_pts, ti_pts = [], []
    for ph in REF_PH_5:
        _, tau_r, theta_r = region_map_5ref[ph]
        Kp_r = local_gain_analytic(ph)
        Ku_r, Pu_r, _, _ = find_ultimate_gain_period(Kp_r, tau_r, theta_r, omega)
        Kc_r, Ti_r, _ = pi_reference(Ku_r, Pu_r, f)
        kc_pts.append((ph, Kc_r)); ti_pts.append((ph, Ti_r))
    kc_sched = PiecewiseLinearGainController(kc_pts)
    ti_sched = PiecewiseLinearGainController(ti_pts)
    return dict(Kc=kc_sched, Ti=ti_sched, Td=0.0, schedule_filter_tau=schedule_filter_tau)

def build_C4b(k):
    """GS-IMC, 5 references, on analytic local gains (Section 1) with nearest-region tau/theta."""
    kc_pts, ti_pts, td_pts = [], [], []
    for ph in REF_PH_5:
        _, tau_r, theta_r = region_map_5ref[ph]
        Kp_r = local_gain_analytic(ph)
        Kc_r, Ti_r, Td_r = imc_pid_rivera(Kp_r, tau_r, theta_r, k * theta_r)
        kc_pts.append((ph, Kc_r)); ti_pts.append((ph, Ti_r)); td_pts.append((ph, Td_r))
    kc_sched = PiecewiseLinearGainController(kc_pts)
    ti_sched = PiecewiseLinearGainController(ti_pts)
    td_sched = PiecewiseLinearGainController(td_pts)
    return dict(Kc=kc_sched, Ti=ti_sched, Td=td_sched, schedule_filter_tau=schedule_filter_tau)

CONTROLLER_BUILDERS = {
    'C1': ('f', build_C1), 'C2': ('f', build_C2), 'C2b': ('f', build_C2b),
    'C3': ('k', build_C3), 'C4': ('k', build_C4), 'C4b': ('k', build_C4b),
}


# ORIGINAL NOTEBOOK CELL 56
standardised_scenario = dict(
    sp_events=[
        {'time': 0,      'SP': 7.0},
        {'time': 60000,  'SP': 8.0},
        {'time': 80000,  'SP': 7.0},
        {'time': 100000, 'SP': 6.0},
    ],
    dist_events=[
        {'time': 20000,  'Ca_factor': 1.05},
        {'time': 65000,  'Ca_factor': 0.95},
        {'time': 120000, 'Ca_factor': 1.05},
    ],
    initial_pH=7.0, tend=170000,
)

wide_scenario = dict(
    sp_events=[
        {'time': 0,      'SP': 4},
        {'time': 5000,   'SP': 5},
        {'time': 30000,  'SP': 6},
        {'time': 60000,  'SP': 7},
        {'time': 80000,  'SP': 8},
        {'time': 100000, 'SP': 9},
        {'time': 120000, 'SP': 10},
    ],
    initial_pH=4, tend=150000,
)

regulatory_scenario = dict(
    sp_events=[{'time': 0, 'SP': 7.0}],
    dist_events=[{'time': 20000, 'Ca_factor': 1.05}],
    initial_pH=7.0, tend=60000,
)

def whole_trace_iae(res):
    """IAE over the entire simulated trace, regardless of whether the setpoint ever
    changes -- same integral as inside calculate_performance_metrics (Shared setup), just
    not gated on a setpoint-change segment. Used for the regulatory scenario above."""
    t_eval, SP_trace, pH_arr = res['t'], res['SP'], res['pH']
    dt = t_eval[1] - t_eval[0]
    return float(np.sum(np.abs(SP_trace - pH_arr)) * dt)

def mean_step_iae(res):
    """Mean IAE across a scenario's setpoint steps (standardised/wide -- both have real
    setpoint changes, so calculate_performance_metrics returns one row per step)."""
    return float(np.mean([m['IAE'] for m in res['metrics']]))

def steps_not_settled(res):
    return int(sum(1 for m in res['metrics'] if m['Settling_Time'] == np.inf))


