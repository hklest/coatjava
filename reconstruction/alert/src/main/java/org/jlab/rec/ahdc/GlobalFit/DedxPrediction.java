package org.jlab.rec.ahdc.GlobalFit;

/**
 * dE/dx prediction of a track on the ln p grid (P4 D6, porting spec STEP 2; pinned Python est/dedx.py): for each grid momentum the
 * geometric mean over the kept hits of the Landau most probable deposit and of beta*gamma at the hit (residual range of each hit),
 * the toy tau and sigma at the track's hit count, and the inversion of the prediction curve. Proton only (mass M, z2 = 1, family 0).
 * One instance per thread: it owns the scratch arrays, {@link #pred0}, {@link #sgc} and the slope of the last {@link #invert}.
 */
public final class DedxPrediction {

	/** proton mass [GeV], the hypothesis of the fit (mscale = 1). */
	public static final double MASS = 0.938272088;
	/** reference chord [mm] of the per-hit deposit (DX:50 CHORD), gas density [g/cm3] (DX:51 RHO_GAS), electron mass [MeV] (DX:52 ME_MEV). */
	public static final double CHORD = 4.0, RHO_GAS = 5.015e-4, ME_MEV = 0.51099895;
	/** mean excitation energies of He and CO2 [MeV] (DX:53 I_HE, I_CO2) and the He mass fraction of the gas (prop.py:124 _W_HE_GAS). */
	public static final double I_HE = 41.8e-6, I_CO2 = 85.0e-6;
	public static final double W_HE_GAS = 0.8 * 1.664e-4 / (0.8 * 1.664e-4 + 0.2 * 1.84e-3);
	/** grid of ln p (p in GeV): np.linspace(ln 0.08, ln 1.6, 900) (DX:560 LNP_GRID). */
	public static final int NG = 900;
	public static final double P_LO = 0.08, P_HI = 1.6;
	public static final double[] LNP_GRID = grid();
	/** a hit that ends below this kinetic energy [MeV] (stopped) makes the whole grid point NaN (DX:104 {@code th < 0.05}). */
	public static final double T_STOP = 0.05;
	/** largest hit count of the toy tables (DX:55 NMAX). */
	public static final int NMAX = 16;
	/** the I [MeV] and zoa of the gas of the proton tables: DX:58-63 gas_constants. */
	public static final double GAS_LN_I = (W_HE_GAS * 0.49968 * Math.log(I_HE) + (1 - W_HE_GAS) * 0.49989 * Math.log(I_CO2))
		/ (W_HE_GAS * 0.49968 + (1 - W_HE_GAS) * 0.49989);

	private static double[] grid() {
		double a = Math.log(P_LO), b = Math.log(P_HI), step = (b - a) / (NG - 1);
		double[] g = new double[NG];
		for (int k = 0; k < NG; k++) g[k] = a + k * step;
		g[NG - 1] = b;
		return g;
	}

	private final StoppingTables st;
	private final double zoa, gasI;
	/** prediction ln(gm Delta_p [keV]) + toy tau (DX:563-606: 0 + lndp + tau) and toy sigma per grid point; NaN where a hit stops. */
	public final double[] pred0 = new double[NG], sgc = new double[NG];
	/** response curve a_run + g pred0 of the last {@link #response}. */
	public final double[] predt = new double[NG];
	/** d curve / d ln p at the solution of the last {@link #invert} (NaN if it failed). */
	public double slope = Double.NaN;
	private final double[] cc = new double[NG], gg = new double[NG];

	public DedxPrediction(StoppingTables st) {
		this.st = st;
		this.zoa = st.zoa[StoppingTables.GAS];
		this.gasI = Math.exp(GAS_LN_I);
	}

	/** most probable energy loss [MeV] of a hit of chord {@link #CHORD} at speed beta (DX:66-76 _landau_mpv, z2 = 1, no density effect). */
	public static double landauMpv(double beta, double zoa, double gasI) {
		double b2 = beta * beta;
		double g2 = 1.0 / (1.0 - b2);
		double xi = 0.153537 * zoa * RHO_GAS * (CHORD / 10.0) / b2;
		return xi * (Math.log(2 * ME_MEV * b2 * g2 / gasI) + Math.log(xi / gasI) + 0.200 - b2);
	}

	/** zoa of the gas in the stopping tables. */
	public double zoa() { return zoa; }

	/** I of the gas [MeV]. */
	public double gasI() { return gasI; }

	/** kinetic energy [MeV] after the grammage x [g/cm2] from tvert [MeV]; 0 if stopped (DX:81-88 _t_hit with mscale 1). */
	private double tHit(double r0, double x) {
		double rr = r0 - x;
		return rr <= 0 ? 0.0 : st.tFromRange(StoppingTables.GAS, rr);
	}

	/**
	 * Fills {@link #pred0} and {@link #sgc} from the grammage x[0..n-1] of the kept hits (all hits of H, not only the used ones),
	 * the track's used-hit count nUsed (toy row, clipped to 1..16) and the half whose toy table is used (DX:94-118 _curves,
	 * 257-264 toy_lookup, 563-606). A track with n = 0 is NaN everywhere.
	 */
	public void compute(double[] x, int n, int nUsed, int half, SwitchConstants.DedxConstants c) {
		int nn = Math.min(Math.max(nUsed, 1), NMAX);
		double[] lb = c.lnBg;
		double[] tau = c.toyTau[half][nn], sig = c.toySig[half][nn];
		for (int j = 0; j < NG; j++) {
			double p = Math.exp(LNP_GRID[j]);
			double tv = (Math.sqrt(p * p + MASS * MASS) - MASS) * 1e3;
			double r0 = st.rangeG(StoppingTables.GAS, tv);
			double sd = 0.0, sb = 0.0;
			boolean ok = n > 0;
			for (int i = 0; i < n && ok; i++) {
				double th = tHit(r0, x[i]);
				if (th < T_STOP) {
					ok = false;
					break;
				}
				double e = th * 1e-3 + MASS;
				double pp = Math.sqrt(e * e - MASS * MASS);
				sd += Math.log(landauMpv(pp / e, zoa, gasI) * 1e3);
				sb += Math.log(pp / MASS);
			}
			if (!ok) {
				pred0[j] = Double.NaN;
				sgc[j] = Double.NaN;
				continue;
			}
			double lndp = sd / n, lnbg = sb / n;
			double xx = Math.min(Math.max(lnbg, lb[0]), lb[lb.length - 1]);
			int k = Math.min(Math.max(lowerBound(lb, lb.length, xx) - 1, 0), lb.length - 2);
			double f = (xx - lb[k]) / (lb[k + 1] - lb[k]);
			pred0[j] = 0.0 + lndp + (tau[k] * (1 - f) + tau[k + 1] * f);
			sgc[j] = sig[k] * (1 - f) + sig[k + 1] * f;
		}
	}

	/** response curve a_run + g pred0 into {@link #predt} (CR:278-279, CD:97-105). */
	public double[] response(SwitchConstants.DedxConstants c) {
		for (int j = 0; j < NG; j++) predt[j] = c.aRun + c.g * pred0[j];
		return predt;
	}

	/** numpy searchsorted, side left, of key in a[0..n-1] (first index with a[i] >= key; NaN key and NaN data as numpy's comparison). */
	static int lowerBound(double[] a, int n, double key) {
		int lo = 0, hi = n;
		while (lo < hi) {
			int mid = lo + ((hi - lo) >> 1);
			if (a[mid] < key || (Double.isNaN(key) && !Double.isNaN(a[mid]))) lo = mid + 1;
			else hi = mid;
		}
		return lo;
	}

	/**
	 * curve value at ln p by linear interpolation on the grid; NaN if lnp is NaN or outside [g0, g899] or the curve is NaN there
	 * (DX:598-605 at_lnp).
	 */
	public static double atLnp(double[] curve, double lnp) {
		if (!(lnp >= LNP_GRID[0] && lnp <= LNP_GRID[NG - 1])) return Double.NaN;
		double x = Math.min(Math.max(lnp, LNP_GRID[0]), LNP_GRID[NG - 1]);
		int j = Math.min(Math.max(lowerBound(LNP_GRID, NG, x) - 1, 0), NG - 2);
		double f = (x - LNP_GRID[j]) / (LNP_GRID[j + 1] - LNP_GRID[j]);
		return curve[j] * (1 - f) + curve[j + 1] * f;
	}

	/**
	 * ln p where the decreasing curve equals lnE; NaN if lnE is not finite or outside the range of the part after the maximum of the
	 * finite curve points (Bragg side excluded). Sets {@link #slope} to d curve / d ln p there (DX:574-595 invert).
	 */
	public double invert(double lnE, double[] curve) {
		slope = Double.NaN;
		if (Double.isNaN(lnE) || Double.isInfinite(lnE)) return Double.NaN;
		int m = 0;
		for (int j = 0; j < NG; j++) {
			double v = curve[j];
			if (!Double.isNaN(v) && !Double.isInfinite(v)) {
				cc[m] = v;
				gg[m] = LNP_GRID[j];
				m++;
			}
		}
		if (m < 3) return Double.NaN;
		int i0 = 0;
		for (int k = 1; k < m; k++) if (cc[k] > cc[i0]) i0 = k;
		int len = m - i0;
		if (len < 3 || lnE > cc[i0] || lnE < cc[m - 1]) return Double.NaN;
		// searchsorted(-cc, -y) on cc[i0..m-1]
		int lo = 0, hi = len;
		double key = -lnE;
		while (lo < hi) {
			int mid = lo + ((hi - lo) >> 1);
			if (-cc[i0 + mid] < key) lo = mid + 1;
			else hi = mid;
		}
		int j = Math.min(Math.max(lo, 1), len - 1);
		double c0 = cc[i0 + j - 1], c1 = cc[i0 + j], g0 = gg[i0 + j - 1], g1 = gg[i0 + j];
		double f = (lnE - c0) / (c1 - c0);
		slope = (c1 - c0) / (g1 - g0);
		return g0 + f * (g1 - g0);
	}
}
