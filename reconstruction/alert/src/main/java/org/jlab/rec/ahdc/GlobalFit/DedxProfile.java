package org.jlab.rec.ahdc.GlobalFit;

/**
 * Post-fit dE/dx profile of a track (P4 D6, porting spec STEP 3; pinned Python est/dedx.py, comb_dedx.py, dedx_run.py): the track
 * statistic ln E, the response curve a + g pred0, the dE/dx-alone ln p, and the one-dimensional profile in ln p around the base fit
 * (chi2 = (ln p - ln p_base)^2 / v - 2 ln f_toy + 2 ln sigma), minimised on the grid inside a window with parabolic refinement. The
 * outputs {@link #pDedx} and {@link #sigPDedx} go to the bank (state and covariance untouched). One instance per thread; the
 * outputs are the public fields of the last {@link #evaluate}. No allocation per call.
 */
public final class DedxProfile {

	/** window of the profile around the base ln p (strict), minimum number of grid points and of used hits (DX:608 win = 1.2, nmin = 3; CD:25 NMIN). */
	public static final double WINDOW = 1.2;
	public static final int NMIN = 3;
	/** outlier floor of the toy pdf and the width of its uniform part (DX:267 floor = 0.02, zwin = 20.0). */
	public static final double FLOOR = 0.02, ZWIN = 20.0;
	/** log-likelihood where not finite (DX:633 {@code -50.0}). */
	public static final double LL_BAD = -50.0;
	/** bank flag of a track with fewer than 3 dE/dx hits whose base fit is kept (Python flag 64, DR:28 FLAG_FEWHITS; Java bit 512, D6). */
	public static final int FLAG_FEWHITS = 512;

	/** hit-count class: -1 below 3 hits, 0 up to 6, 1 for 7, 2 from 8 (DR:31-32 ncls). */
	public static int ncls(int n) { return n < 3 ? -1 : n <= 6 ? 0 : n == 7 ? 1 : 2; }

	/** the prediction (curves of the last evaluated track). */
	public final DedxPrediction pred;
	/** ln E of the used hits and their number (DX:517-526 track_stat); NaN for none. */
	public double lnE;
	public int n;
	/** dE/dx alone (CD:97-105): ln p of the inversion (NaN for n &lt; 3) and its error (ln p, from the unmasked inversion). */
	public double dedxLnp, dedxSigLnp;
	/** class variance used (1.0 for n &lt; 3), profile minimum and error in ln p (NaN where undefined). */
	public double classVar, profLnp, profErr;
	/** outputs of the bank: combined p and sigma_p (base ones when n &lt; 3), NaN if the base failed or the profile is not finite. */
	public double pDedx, sigPDedx;
	/** 0 if p_dedx is valid, else 1; flags {@link #FLAG_FEWHITS} or 0. */
	public int status, flags;
	/** the prompt class that was given (DR:35-39), echoed. */
	public boolean prompt;

	private final double[] c2 = new double[DedxPrediction.NG];

	public DedxProfile(StoppingTables st) { pred = new DedxPrediction(st); }

	/**
	 * Evaluates the dE/dx term of one track. x, lnu, use hold the kept hits (the on-track table H) at indices 0..nKept-1: grammage
	 * [g/cm2], ln u and the use flag of the calibrated hit. {@code baseP}, {@code baseSigP} and {@code baseStatus} are the fit the
	 * profile refines (default fit for ge10_dedx), {@code prompt} the wedge class computed by the caller, {@code half} the half that
	 * uses the constants. Results in the public fields.
	 */
	public void evaluate(double[] x, double[] lnu, boolean[] use, int nKept, double baseP, double baseSigP, int baseStatus, boolean prompt,
		int half, SwitchConstants.DedxConstants c) {
		this.prompt = prompt;
		// track statistic (DX:517-526)
		int nu = 0;
		double s = 0.0;
		for (int i = 0; i < nKept; i++)
			if (use[i] && !Double.isNaN(lnu[i]) && !Double.isInfinite(lnu[i])) {
				nu++;
				s += Math.exp(-lnu[i]);
			}
		n = nu;
		lnE = nu > 0 ? -Math.log(s / nu) : Double.NaN;

		pred.compute(x, nKept, nu, half, c);
		double[] predt = pred.response(c), sgc = pred.sgc;
		double ld = pred.invert(lnE, predt);
		dedxSigLnp = Math.abs(DedxPrediction.atLnp(sgc, ld) / pred.slope);
		dedxLnp = nu >= NMIN ? ld : Double.NaN;

		int nc = ncls(nu);
		double v = 1.0;
		if (nc >= 0) {
			v = c.classVar[half][prompt ? 1 : 0][nc];
			if (Double.isNaN(v) || Double.isInfinite(v)) v = 1.0;
		}
		classVar = v;
		boolean ok = baseStatus == 0 && !Double.isNaN(baseP) && !Double.isInfinite(baseP) && baseP > 0;
		double lnpB = ok ? Math.log(baseP) : Double.NaN;
		profile(lnpB, v, predt, sgc, half, c);

		boolean few = nu < NMIN && ok;
		double p = Double.isNaN(profLnp) || Double.isInfinite(profLnp) ? Double.NaN : Math.exp(profLnp);
		double sig = p * profErr;
		if (few) {
			p = baseP;
			sig = baseSigP;
		}
		boolean good = ok && !Double.isNaN(p) && !Double.isInfinite(p) && p > 0;
		status = good ? 0 : 1;
		pDedx = good ? p : Double.NaN;
		sigPDedx = sig;
		flags = few ? FLAG_FEWHITS : 0;
	}

	/** DX:608-646 combine_profile for one track: sets {@link #profLnp} and {@link #profErr}. */
	private void profile(double lnpZt, double vZt, double[] predt, double[] sgc, int half, SwitchConstants.DedxConstants c) {
		profLnp = Double.NaN;
		profErr = Double.NaN;
		if (Double.isNaN(lnpZt) || Double.isInfinite(lnpZt) || Double.isNaN(vZt) || Double.isInfinite(vZt)) return;
		double[] g = DedxPrediction.LNP_GRID;
		int lo0 = -1, hi0 = -1;
		for (int k = 0; k < DedxPrediction.NG; k++)
			if (Math.abs(g[k] - lnpZt) < WINDOW) {
				if (lo0 < 0) lo0 = k;
				hi0 = k;
			}
		int m = lo0 < 0 ? 0 : hi0 - lo0 + 1;
		if (m < 3) return;
		for (int k = 0; k < m; k++) {
			double d = g[lo0 + k] - lnpZt;
			c2[k] = d * d / vZt;
		}
		if (Double.isNaN(lnE) || Double.isInfinite(lnE) || n < NMIN) {
			profLnp = lnpZt;
			profErr = Math.sqrt(vZt);
			return;
		}
		int nn = Math.min(Math.max(n, 1), DedxPrediction.NMAX);
		double[] pdf = c.pdf[half][nn];
		for (int k = 0; k < m; k++) {
			double z = (lnE - predt[lo0 + k]) / sgc[lo0 + k];
			double ll = logPdf(z, pdf, c.zgrid) - Math.log(sgc[lo0 + k]);
			if (Double.isNaN(ll) || Double.isInfinite(ll)) ll = LL_BAD;
			c2[k] -= 2 * ll;
		}
		int j = 0;
		for (int k = 1; k < m; k++) if (c2[k] < c2[j]) j = k;
		if (j > 0 && j < m - 1) {
			double y0 = c2[j - 1], y1 = c2[j], y2 = c2[j + 1];
			double den = y0 - 2 * y1 + y2;
			double sh = den > 0 ? 0.5 * (y0 - y2) / den : 0.0;
			profLnp = g[lo0 + j] + sh * (g[1] - g[0]);
			int lo = j;
			while (lo > 0 && c2[lo] - c2[j] < 1.0) lo--;
			int hi = j;
			while (hi < m - 1 && c2[hi] - c2[j] < 1.0) hi++;
			profErr = 0.5 * (g[lo0 + hi] - g[lo0 + lo]);
		} else {
			profLnp = g[lo0 + j];
		}
	}

	/** ln of the standardised track pdf with the outlier floor (DX:267-273 toy_logpdf); z outside the grid or NaN gives the floor only. */
	private static double logPdf(double z, double[] pdf, double[] zg) {
		int nz = zg.length;
		int k = Math.min(Math.max(DedxPrediction.lowerBound(zg, nz, z) - 1, 0), nz - 2);
		boolean inside = z >= zg[0] && z < zg[nz - 1];
		double f = inside ? pdf[k] : 0.0;
		return Math.log((1 - FLOOR) * f + FLOOR / ZWIN);
	}
}
