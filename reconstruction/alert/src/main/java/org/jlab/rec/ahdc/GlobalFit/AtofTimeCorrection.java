package org.jlab.rec.ahdc.GlobalFit;

/**
 * The ATOF time correction of the time switch (plan D2): the bank time of a wedge row is brought back to the engine time before
 * any T0, then the calibrated T0 and the gain-group walk are applied, as the pinned Python does:
 * <pre>
 *   t = t_bank + dep_t0 + (run - dep_icpt) * dep_slope - t0[ch] - W_g(x)           (atof_cal.py:10, 202-212: t_c = t_dep + t0_dep + dt_dep - t0_c - dt_c;
 *                                                                                   tcal.py:292-315 corrected(): raw - t_ch - walk_eval)
 *   ToT = rint(energy / TOT2ENERGY_WEDGE / ATT)                                     (atof_cal.py:41-42, 174 tot_of; ATOFHit.convertTotToEnergy wedge)
 *   x   = ln(max(ToT, 1) / p50[ch]),  g = group[ch]                                 (tcal.py:181-184 _xvar uw_g; 102-106 group_of)
 *   W_g = np.interp(x, walkX[g], walkW[g]), linear, clamped at the end nodes        (tcal.py:292-315 walk_eval)
 * </pre>
 * Bar rows (component 10) keep the bank time (tcal.py:292-315 corrected: "deployed time elsewhere"). The start time is not
 * subtracted here (as in the Python dumps). Pure function of the row and the constants; allocation-free except for
 * {@code ev.aTime} when it is null or shorter than the number of ATOF rows.
 */
public final class AtofTimeCorrection {

	/** ToT to energy of a wedge and its attenuation factor (atof_cal.py:41-42). */
	public static final double TOT2ENERGY_WEDGE = 1.956 * 2.0 / 1000, ATT = Math.exp(10.0 / 1600.0);

	private AtofTimeCorrection() {}

	/** fills {@code ev.aTime} (rows 0 .. nAtof-1) from {@code ev.tBank}, {@code ev.aEnergy}, {@code ev.sector/layer/comp/isBar} for run number {@code run}. */
	public static void correct(GFitEvent ev, SwitchConstants.TimeChannels tc, int run) {
		int n = ev.nAtof();
		if (ev.aTime == null || ev.aTime.length < n) ev.aTime = new double[n];
		for (int i = 0; i < n; i++) {
			if (ev.isBar[i]) {
				ev.aTime[i] = ev.tBank[i];
				continue;
			}
			int ch = SwitchConstants.channel(ev.sector[i], ev.layer[i], ev.comp[i]);
			double tot = Math.rint(ev.aEnergy[i] / TOT2ENERGY_WEDGE / ATT);
			double x = Math.log(Math.max(tot, 1.0) / tc.p50[ch]);
			ev.aTime[i] = ev.tBank[i] + tc.depT0[ch] + (run - tc.depIcpt[ch]) * tc.depSlope[ch] - tc.t0[ch]
				- interp(x, tc.walkX[tc.group[ch]], tc.walkW[tc.group[ch]]);
		}
	}

	/** np.interp for ascending nodes: end values outside, linear between (NaN stays NaN). */
	static double interp(double x, double[] xp, double[] fp) {
		int n = xp.length;
		if (x != x) return x;
		if (x <= xp[0]) return fp[0];
		if (x >= xp[n - 1]) return fp[n - 1];
		int j = 0;
		while (xp[j + 1] <= x) j++;
		return (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]) * (x - xp[j]) + fp[j];
	}
}
