# ============================================================
# CLOSED-LOOP SIMULATOR (used by every controller in this notebook)
# ============================================================
# Kc, Ti, Td may each be a constant (fixed-gain controller) or a
# PiecewiseLinearGainController (gain-scheduled on measured pH) -- the same function
# handles the conventional, IMC, IMC-scheduled and pure gain-scheduled controllers.

def audited_sim(Kc, Ti, Td=0.0, sp_events=None, dist_events=None,
                         initial_pH=None, theta_dead=None,
                         F_min=0.0, F_max=None, tend=100000, n_points=None,
                         schedule_filter_tau=0.0,
                         title="Closed-Loop pH Control", show_mv=False, plot=True, noise=None, delay_mode="legacy"):
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
    measured = np.zeros(n)
    measured[0] = initial_pH + noise[0]
    e[0] = sp_at(0) - measured[0]
    F[0] = solve_Qb_for_pH(initial_pH, guess=0.0)
    x0 = [x0_val]
    MV_buffer = np.full(theta_steps, F[0])

    t_span_start = 0.0
    pH_sched_filt = measured[0]
    for k in range(1, n):
        current_SP = sp_at(t_eval[k])
        if schedule_filter_tau > 0:
            alpha = dt / (schedule_filter_tau + dt)
            pH_sched_filt = pH_sched_filt + alpha * (measured[k - 1] - pH_sched_filt)
            current_pH = pH_sched_filt
        else:
            current_pH = measured[k - 1]
        current_Kc = scheduled(Kc, current_pH)
        current_Ti = scheduled(Ti, current_pH)
        current_Td = scheduled(Td, current_pH)

        e[k] = current_SP - measured[k - 1]
        de_val = e[k] - e[k - 1]
        pv_prev1 = measured[k - 2] if k >= 2 else measured[0]
        pv_prev2 = measured[k - 3] if k >= 3 else measured[0]
        de_val2 = measured[k - 1] - 2 * pv_prev1 + pv_prev2

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

        if delay_mode == 'fractional':
            ratio = theta_dead / dt
            lag = int(np.floor(ratio))
            fraction = ratio - lag
            older = F[max(0, k-lag-1)]
            newer = F[max(0, k-lag)]
            mid = advance(dxdt_func, x0, older, current_Ca, fraction*dt)
            sol = advance(dxdt_func, [mid.y[0,-1]], newer, current_Ca, (1-fraction)*dt)
        else:
            sol = advance(dxdt_func, x0, Qb_delayed, current_Ca, dt)
        x0 = [sol.y[0, -1]]
        t_span_start += dt
        pH_arr[k] = x_to_pH(x0[0])
        measured[k] = pH_arr[k] + noise[k]

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