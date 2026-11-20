package org.jlab.rec.ahdc.GlobalFit;

/**
 * Prediction, cost, gradient and both Hessians of one track at a fixed state: ports fit._predict (fit.py:328-387) for the
 * helix model (imodel 0) and the propagation model (imodel 1, forward differences of the whole propagation; the kink
 * columns restart from the nominal snapshots) and fit._cost (fit.py:391-442), with the hit-list setup of
 * fit._eval_track (fit.py:881-906). All scratch lives in the {@link GFitWorkspace}; the ATOF term comes from {@link AtofTerm}.
 */
public final class TrackCost {

	private static final int NPAR = GFitWorkspace.NPAR, NC = GFitWorkspace.NC, NK = GFitWorkspace.NK;

	private final Propagator prop;
	private final CellGeometry geo;
	private final GFitConfig cfg;
	private final MaterialStack stk;
	private final double bz, beamX, beamY, sigZ, sigD0, wMin;
	private final boolean zOn, d0On;
	/** kink columns restart from the snapshots; the equivalence test switches it off to run full propagations. */
	boolean useRestart = true;

	public TrackCost(Propagator prop, CellGeometry geo, GFitConfig cfg) {
		this.prop = prop;
		this.stk = prop.stack();
		this.geo = geo;
		this.cfg = cfg;
		bz = cfg.bz(); beamX = cfg.beamX(); beamY = cfg.beamY(); sigZ = cfg.sigZ(); sigD0 = cfg.sigD0(); wMin = cfg.wMin();
		zOn = cfg.zOn(); d0On = cfg.d0On();
	}

	/** hit list of the track's B0 hits (fit._fit_one 603-619, _eval_track 881-894); all kinds 0, warm arc lengths -1. */
	public void loadHits(GFitTrack t, GFitWorkspace ws) {
		int nb = t.nHits();
		ws.nh = nb;
		for (int k = 0; k < ws.cap; k++) {
			ws.lg[k] = 0; ws.lli[k] = 0; ws.lw[k] = 0; ws.lism[k] = 0; ws.lcls[k] = 0; ws.lkind[k] = 0; ws.lrow[k] = -1;
			ws.lm[k] = 0.0; ws.lsig[k] = 1.0; ws.lh[k] = 1.0; ws.sarc[k] = -1.0;
		}
		int[] layoff = geo.layoff;
		for (int k = 0; k < nb; k++) {
			int g = t.hg[k];
			ws.lg[k] = g; ws.lm[k] = t.hm[k]; ws.lism[k] = t.hism[k]; ws.lsig[k] = t.hsig[k]; ws.lh[k] = t.hh[k];
			ws.lcls[k] = t.hcls[k]; ws.lrow[k] = t.hrow[k];
			for (int li = 0; li < 8; li++)
				if (g >= layoff[li] && g < layoff[li + 1]) {
					ws.lli[k] = li;
					ws.lw[k] = g - layoff[li];
				}
		}
	}

	/** signed doca, PCA point and flag of every hit of the table in ws.out into d[k], x[3k..3k+2], pf[k] (fit._prop_run, fit.py:322-324). */
	private void pcaAll(GFitWorkspace ws, int nn, double[] d, double[] x, int[] pf) {
		double[] r = ws.pcaR;
		for (int k = 0; k < ws.nh; k++) {
			int f = HitTerm.pca1(ws.out, nn, geo.wo, geo.we, ws.lg[k], ws.pcaG, ws.pcaD, r);
			d[k] = r[0]; x[3 * k] = r[1]; x[3 * k + 1] = r[2]; x[3 * k + 2] = r[3]; pf[k] = f;
		}
	}

	/**
	 * fit._predict: ws.D (signed d), ws.S (cell margin; computed when want or wantS, otherwise left as it was), ws.JD / ws.JS
	 * (nh x NPAR) with want, ws.PF; the nominal Highland theta0 into th0out; with ws.atofOn and imodel 1 the 16 intercept records
	 * into ws.ia. Returns the trajectory flags (imodel 1) or 0.
	 */
	public int predict(double[] u, int imodel, int npar, boolean want, boolean wantS, double[] th0out, GFitWorkspace ws) {
		int nh = ws.nh;
		if (imodel == 0) {
			double[] o = ws.o9, o2 = ws.o9b;
			for (int k = 0; k < nh; k++) {
				int g = ws.lg[k];
				double sd = HelixModel.hhit(u, geo.wo, geo.wu, g, bz, beamX, beamY, o, ws.sarc[k], want, ws.jx5, ws.j53, ws.sArc1);
				double sa = ws.sArc1[0];
				boolean fd = !(Double.isNaN(sd) || Double.isInfinite(sd));
				ws.D[k] = sd;
				ws.PF[k] = fd ? 0 : 1;
				if (!(Double.isNaN(sa) || Double.isInfinite(sa))) ws.sarc[k] = sa;
				if (want) for (int j = 0; j < npar; j++) ws.JD[k * NPAR + j] = (j < 5 && fd) ? ws.jx5[j] : 0.0;
				if (wantS) {
					if (!fd) {
						ws.S[k] = Double.NaN;
						if (want) for (int j = 0; j < NPAR; j++) ws.JS[k * NPAR + j] = 0.0;
						continue;
					}
					HelixModel.pt(u, sa, bz, beamX, beamY, o);
					geo.eval(ws.lli[k], ws.lw[k], o[0], o[1], o[2], ws.cell4);
					ws.S[k] = ws.cell4[0];
					if (want) {
						double gx = ws.cell4[1], gy = ws.cell4[2], gz = ws.cell4[3];
						for (int j = 0; j < npar; j++) ws.JS[k * NPAR + j] = 0.0;
						for (int j = 0; j < 5; j++) {
							double hj = (j == 0 || j == 3) ? 1e-4 : 1e-6;
							double[] up = ws.upert;
							System.arraycopy(u, 0, up, 0, NPAR);
							up[j] += hj;
							double s2 = HelixModel.closest(up, geo.wo, geo.wu, g, bz, beamX, beamY, o2, sa);
							if (!(Double.isNaN(s2) || Double.isInfinite(s2))) {
								HelixModel.pt(up, s2, bz, beamX, beamY, o2);
								ws.JS[k * NPAR + j] = (gx * (o2[0] - o[0]) + gy * (o2[1] - o[1]) + gz * (o2[2] - o[2])) / hj;
							}
						}
					}
				}
			}
			return 0;
		}
		boolean saveOld = ws.saveSnapshots;
		ws.saveSnapshots = want;
		int nn = prop.run(u, ws);
		ws.saveSnapshots = false;
		int tflag = ws.flags;
		System.arraycopy(ws.theta0, 0, th0out, 0, GFitWorkspace.NPL);
		pcaAll(ws, nn, ws.D, ws.nomX, ws.PF);
		if (ws.atofOn) Propagator.atofIntercepts(ws.out, nn, ws.regnode, ws.aux, stk.rOut, ws.ia, 0);
		double[] gg = ws.gradG;
		if (wantS || want) {
			for (int k = 0; k < nh; k++) {
				geo.eval(ws.lli[k], ws.lw[k], ws.nomX[3 * k], ws.nomX[3 * k + 1], ws.nomX[3 * k + 2], ws.cell4);
				ws.S[k] = ws.cell4[0];
				gg[3 * k] = ws.cell4[1]; gg[3 * k + 1] = ws.cell4[2]; gg[3 * k + 2] = ws.cell4[3];
			}
		}
		if (want) {
			double[] up = ws.upert;
			int[] pf1 = ws.PF1;
			for (int j = 0; j < npar; j++) {
				double hj = j < 5 ? GFitConfig.stepU(j) : GFitConfig.STEP_KINK;
				System.arraycopy(u, 0, up, 0, NPAR);
				up[j] += hj;
				int nnj = j < 5 || !useRestart ? prop.run(up, ws) : prop.runRestart(up, (j - 5) / 2, ws);
				if (ws.atofOn) Propagator.atofIntercepts(ws.out, nnj, ws.regnode, ws.aux, stk.rOut, ws.ia, (1 + j) * 20);
				pcaAll(ws, nnj, ws.pertD, ws.pertX, pf1);
				for (int k = 0; k < nh; k++) {
					ws.JD[k * NPAR + j] = (ws.pertD[k] - ws.D[k]) / hj;
					ws.JS[k * NPAR + j] = (gg[3 * k] * (ws.pertX[3 * k] - ws.nomX[3 * k]) + gg[3 * k + 1] * (ws.pertX[3 * k + 1] - ws.nomX[3 * k + 1])
						+ gg[3 * k + 2] * (ws.pertX[3 * k + 2] - ws.nomX[3 * k + 2])) / hj;
				}
			}
		}
		ws.saveSnapshots = saveOld;
		return tflag;
	}

	/**
	 * fit._cost (term dbox, ATOF term not added here): returns F = hit + z prior + d0 prior + kink prior; sets ws.tflag, ws.outH,
	 * the term fields and, with want, ws.g and the symmetric ws.hc (clipped Gauss-Newton) and ws.he (exact) over npar parameters.
	 * th0f = frozen Highland theta0 per plane (used by the kink prior when th0_mode is 0).
	 */
	public double cost(double[] u, int imodel, int npar, boolean want, double ze, double[] th0f, GFitWorkspace ws) {
		boolean wantS = false;   // term == 2 or add_hits: neither is ported
		double[] g = ws.g, hc = ws.hc, he = ws.he;
		ws.tflag = predict(u, imodel, npar, want, wantS, ws.th0, ws);
		if (want) {
			for (int j = 0; j < npar; j++) {
				g[j] = 0.0;
				for (int m = 0; m < npar; m++) {
					hc[j * NPAR + m] = 0.0;
					he[j * NPAR + m] = 0.0;
				}
			}
		}
		double F = 0.0;
		double[] b3 = ws.box3;
		int nh = ws.nh;
		for (int k = 0; k < nh; k++) {
			double d = ws.D[k];
			HitTerm.dbox(d, ws.lm[k], ws.lsig[k], ws.lh[k], false, cfg, b3);
			double l = -0.5 * b3[0], ld = -0.5 * b3[1], ldd = -0.5 * b3[2];
			double ls = 0.0, lss = 0.0, lds = 0.0;
			ws.outH[3 * k] = l; ws.outH[3 * k + 1] = d; ws.outH[3 * k + 2] = ws.S[k];
			F -= 2.0 * l;
			if (want && finite(d) && allFinite(ws.JD, k * NPAR, npar)) {
				double a = Math.max(-ldd, wMin), c = Math.max(-lss, 0.0), b = -lds;
				if (b * b > a * c) b = Math.copySign(Math.sqrt(a * c), b);
				for (int j = 0; j < npar; j++) {
					double gd = ws.JD[k * NPAR + j], gs = 0.0;
					g[j] += -2.0 * (ld * gd + ls * gs);
					for (int m = j; m < npar; m++) {
						double hd = ws.JD[k * NPAR + m], hs = 0.0;
						hc[j * NPAR + m] += 2.0 * (a * gd * hd + c * gs * hs + b * (gd * hs + gs * hd));
						he[j * NPAR + m] += 2.0 * (-ldd * gd * hd - lss * gs * hs - lds * (gd * hs + gs * hd));
					}
				}
			}
		}
		ws.fHit = F;
		ws.fZ = 0.0; ws.fD0 = 0.0; ws.fK = 0.0;
		if (zOn) {
			double r = (u[3] - ze) / sigZ;
			F += r * r;
			ws.fZ = r * r;
			if (want) {
				g[3] += 2.0 * r / sigZ;
				hc[3 * NPAR + 3] += 2.0 / (sigZ * sigZ);
				he[3 * NPAR + 3] += 2.0 / (sigZ * sigZ);
			}
		}
		if (d0On) {
			double r = u[0] / sigD0;
			F += r * r;
			ws.fD0 = r * r;
			if (want) {
				g[0] += 2.0 * r / sigD0;
				hc[0] += 2.0 / (sigD0 * sigD0);
				he[0] += 2.0 / (sigD0 * sigD0);
			}
		}
		ws.fAtof = 0.0;
		if (imodel == 1 && ws.atofOn) {      // ATOF position term (fit.py:428-429): value from the nominal record, derivatives through the FD records
			ws.fAtof = AtofTerm.termCost(u, ws.ia, want, npar, ze, ws.atof, g, hc, he);
			F += ws.fAtof;
		}
		if (npar > 5) {
			double[] th = cfg.th0Mode() == 1 ? ws.th0 : th0f;
			for (int j = 0; j < NK; j++) {
				int pl = j / 2;
				double sk = th[pl] > 0 ? cfg.fMs(pl) * th[pl] : 1e-4;
				double r = u[5 + j] / sk;
				F += r * r;
				ws.fK += r * r;
				if (want) {
					g[5 + j] += 2.0 * r / sk;
					hc[(5 + j) * NPAR + 5 + j] += 2.0 / (sk * sk);
					he[(5 + j) * NPAR + 5 + j] += 2.0 / (sk * sk);
				}
			}
		}
		if (want) {
			for (int j = 0; j < npar; j++)
				for (int m = 0; m < j; m++) {
					hc[j * NPAR + m] = hc[m * NPAR + j];
					he[j * NPAR + m] = he[m * NPAR + j];
				}
		}
		return F;
	}

	private static boolean finite(double x) { return !(Double.isNaN(x) || Double.isInfinite(x)); }

	private static boolean allFinite(double[] a, int off, int n) {
		for (int i = 0; i < n; i++) if (!finite(a[off + i])) return false;
		return true;
	}
}
