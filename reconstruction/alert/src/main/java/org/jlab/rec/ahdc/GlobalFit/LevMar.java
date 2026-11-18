package org.jlab.rec.ahdc.GlobalFit;

/**
 * Levenberg-Marquardt minimisation of {@link TrackCost} on a state in place: port of fit._lm (fit.py:451-498) with the
 * in-bounds test fit._inb (fit.py:445-448). Damping starts at 1e-3, is scaled by max(|H_jj|, 1e-10), shrinks by 10 on an
 * accepted step (floor 1e-9) and grows by 10 otherwise; a step leaving the bounds or with a non-finite solution counts as
 * rejected (give up above 1e10), and no step lowering the cost above 1e8 is reported as converged. Immutable and thread
 * safe; the scratch lives in the workspace, the outcome in {@code ws.lmConv} and {@code ws.lmIt}.
 */
public final class LevMar {

	private static final int NPAR = GFitWorkspace.NPAR;

	private final TrackCost cost;
	private final double tolCost, tolStep;

	public LevMar(TrackCost cost, GFitConfig cfg) {
		this.cost = cost;
		this.tolCost = cfg.tolCost();
		this.tolStep = cfg.tolStep();
	}

	/** fit._inb: 0.05 &lt; theta &lt; pi - 0.05, 0.03 &lt; p &lt; 3 GeV, |d0| &lt; 20 mm (false for NaN). */
	public static boolean inb(double[] u) {
		double p = Math.exp(u[2]);
		return 0.05 < u[4] && u[4] < Math.PI - 0.05 && 0.03 < p && p < 3.0 && Math.abs(u[0]) < 20.0;
	}

	private static boolean allFinite(double[] a, int n) {
		for (int i = 0; i < n; i++) if (Double.isNaN(a[i]) || Double.isInfinite(a[i])) return false;
		return true;
	}

	/** the leading npar x npar block (stride NPAR) is finite. */
	private static boolean blockFinite(double[] a, int npar) {
		for (int j = 0; j < npar; j++)
			for (int m = 0; m < npar; m++) {
				double x = a[j * NPAR + m];
				if (Double.isNaN(x) || Double.isInfinite(x)) return false;
			}
		return true;
	}

	/**
	 * Minimises over u[0..npar) in place and returns the final cost; ws.lmConv and ws.lmIt hold the converged flag and the number
	 * of iterations. th0f is the frozen Highland theta0 of the kink prior. A singular damped matrix (zero pivot) throws, as numpy does.
	 */
	public double run(double[] u, int imodel, int npar, double ze, double[] th0f, int maxit, GFitWorkspace ws) {
		double F = cost.cost(u, imodel, npar, true, ze, th0f, ws);
		double lam = 1e-3;
		boolean conv = false;
		int it = 0;
		double[] gs = ws.lmGs, hs = ws.lmHs, a = ws.lmA, dx = ws.lmDx, un = ws.lmUn, uu = ws.lmUu;
		for (int j = 0; j < npar; j++) gs[j] = ws.g[j];
		for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) hs[j * NPAR + m] = ws.hc[j * NPAR + m];
		while (it < maxit) {
			it++;
			if (!(blockFinite(hs, npar) && allFinite(gs, npar))) break;
			for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) a[j * NPAR + m] = hs[j * NPAR + m];
			for (int j = 0; j < npar; j++) {
				double d = Math.abs(hs[j * NPAR + j]);
				a[j * NPAR + j] += lam * (1e-10 > d ? 1e-10 : d);
			}
			for (int j = 0; j < npar; j++) dx[j] = -gs[j];
			if (!SmallLinAlg.solve(a, NPAR, npar, dx, ws.piv)) throw new ArithmeticException("singular matrix in LM step");
			if (!allFinite(dx, npar)) {
				lam *= 10.0;
				if (lam > 1e10) break;
				continue;
			}
			System.arraycopy(u, 0, un, 0, NPAR);
			for (int j = 0; j < npar; j++) un[j] = u[j] + dx[j];
			if (!inb(un)) {
				lam *= 10.0;
				if (lam > 1e10) break;
				continue;
			}
			System.arraycopy(un, 0, uu, 0, NPAR);
			double Fn = cost.cost(uu, imodel, npar, false, ze, th0f, ws);
			if (Fn <= F) {
				double dF = F - Fn;
				for (int j = 0; j < npar; j++) u[j] = un[j];
				F = cost.cost(u, imodel, npar, true, ze, th0f, ws);
				for (int j = 0; j < npar; j++) gs[j] = ws.g[j];
				for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) hs[j * NPAR + m] = ws.hc[j * NPAR + m];
				double l2 = lam * 0.1;
				lam = 1e-9 > l2 ? 1e-9 : l2;
				double st = 0.0;
				for (int j = 0; j < npar; j++) {
					double h = hs[j * NPAR + j];
					double v = Math.abs(dx[j]) * Math.sqrt((1e-30 > h ? 1e-30 : h) * 0.5);
					if (v > st) st = v;
				}
				if (dF < tolCost && st < tolStep) {
					conv = true;
					break;
				}
			} else {
				lam *= 10.0;
				if (lam > 1e8) {     // no step lowers the cost: a minimum to numerical precision
					conv = true;
					break;
				}
			}
		}
		ws.lmConv = conv;
		ws.lmIt = it;
		return F;
	}
}
