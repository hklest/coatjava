package org.jlab.rec.ahdc.GlobalFit;

/**
 * Per-hit dE/dx geometry and calibrated hit quantities of one track at a frozen fit state (pinned est/dedx.py: DX:277-285
 * _cellpath, DX:366-442 hit_table, DX:121-134 _hit_shape with DX:66-76 _landau_mpv and DX:72-79 _t_hit, DX:90-91 _grammage,
 * DX:446-455 design, DX:517-526 track_stat; dedx_cal.py:99-101 hit_use, DC:116-183 calibrate; proton hypothesis, mscale = 1).
 *
 * <p>{@link #compute} propagates the state, and for every original hit with a bank row (kind 0, row >= 0): the PCA arc
 * length s (HitTerm.pca1 slot 10), the in-cell path (0.04 mm steps over s +- 10 mm, cell margin > 0), the grammage X from the
 * vertex to the hit, and with the calibration constants ln u = ln(max(ADC - offset, 1e-3)) - ln G - (ln Delta_p - gm). Hits
 * with no PCA, a failed interpolation or a path <= 0.05 mm are dropped (not kept). Results are in the public fields, per hit
 * indexed like the result's hit list (original hits come first, the same index as GFitTrack.hg); entries of hits that are not kept hold
 * NaN (path, s, X, lnu) and false (kept, use). One instance per thread, never shared; nothing is allocated per call.
 *
 * <p>Constants that D2's holder feeds as plain arguments: {@code offset} = dx_offset (8, additive ADC offsets by layer index,
 * subtracted), {@code f3} = dx_f3[half] (10: gains of layers 1..7, doca per mm, ped per 50, path exponent), {@code f3Norm} =
 * dx_f3_norm[half], with half = tag_event % 2 of the track.
 */
public final class DedxHits {

	/** vertex-to-hit scan: half range and step of the cell path (DX:277), minimum path of a kept hit (DX:441). */
	static final double HALF = 10.0, DS = 0.04, MIN_PATH = 0.05;
	/** reference chord of the per-hit deposit (DX:43), median pedestal of the design row (DX:42), density of the gas (DX:44), electron mass MeV (DX:45). */
	static final double CHORD = 4.0, PED_REF = 262.0, RHO_GAS = 5.015e-4, ME_MEV = 0.51099895;
	private static final double M = Propagator.MASS;
	private static final int NC = GFitWorkspace.NC;
	private static final int C_P = 7;
	private static final double I_HE = 41.8e-6, I_CO2 = 85.0e-6;

	private final Propagator prop;
	private final CellGeometry geo;
	private final StoppingTables tab;
	private final GFitWorkspace ws = new GFitWorkspace(0);
	private final double[] gs = new double[GFitWorkspace.NMAX], ds = new double[GFitWorkspace.NMAX], r12 = new double[12], c4 = new double[4];
	private final double zoa, iGas;

	/** per track: status (0 ok, 1 fit failed, 2 propagation failed), vertex momentum, hits on track, used hits, mean ln Delta_p of the used hits, ln E (NaN with no used hit). */
	public int status, nOn, nUsed;
	public double pTraj0, gm, lnE;
	/** per hit (index of the result's hit list): kept in the on-track table H, path [mm], PCA arc length, grammage X [g/cm2], ln u, use flag. */
	public final boolean[] kept, use;
	public final double[] path, s, X, lnu;
	private final double[] lndp;

	public DedxHits(Propagator prop, CellGeometry geo, StoppingTables tab, int maxHits) {
		this.prop = prop;
		this.geo = geo;
		this.tab = tab;
		int n = maxHits + GFitConfig.NADD;
		kept = new boolean[n]; use = new boolean[n];
		path = new double[n]; s = new double[n]; X = new double[n]; lnu = new double[n]; lndp = new double[n];
		zoa = tab.zoa[StoppingTables.GAS];
		double w = StoppingTables.W_HE_GAS;  // DX:58-63 gas_constants
		double lnI = (w * 0.49968 * Math.log(I_HE) + (1 - w) * 0.49989 * Math.log(I_CO2)) / (w * 0.49968 + (1 - w) * 0.49989);
		iGas = Math.exp(lnI);
	}

	/**
	 * Hit quantities of the track whose fit state and hit list are in {@code r} (u = d0, phi0, ln p, z0, theta, 10 kinks; nh,
	 * hitKind, hitRow, hitWire of the final hit set; only kind 0 hits with row >= 0 are used). hAdc, hPed and hDoca [mm] are
	 * per hit, indexed like r (entries of unused hits are ignored); offset, f3 and f3Norm as described in the class comment.
	 * Fills the public fields and returns the status (0 ok, 1 r.status != 0, 2 propagation failed or fewer than 3 nodes).
	 */
	public int compute(GFitResult r, double[] hAdc, double[] hPed, double[] hDoca, double[] offset, double[] f3, double f3Norm) {
		final int nh = r.nh;
		nOn = 0; nUsed = 0; pTraj0 = Double.NaN; gm = 0.0; lnE = Double.NaN;
		for (int k = 0; k < nh; k++) {
			kept[k] = false; use[k] = false;
			path[k] = Double.NaN; s[k] = Double.NaN; X[k] = Double.NaN; lnu[k] = Double.NaN; lndp[k] = Double.NaN;
		}
		if (r.status != 0) return status = 1;
		int nn = prop.run(r.u, ws);  // kin = est_kinks + [0, 0]: plane 5 has none (DX:380)
		if (nn < 3) return status = 2;
		status = 0;
		final double[] a = ws.out;
		final double p0 = a[C_P];
		pTraj0 = p0;
		final double tv = (Math.sqrt(p0 * p0 + M * M) - M) * 1e3;  // Phys.tkin, DX:147
		final int[] layoff = geo.layoff;
		for (int k = 0; k < nh; k++) {
			if (r.hitKind[k] != 0 || r.hitRow[k] < 0) continue;
			int g = r.hitWire[k];
			HitTerm.pca1(a, nn, geo.wo, geo.we, g, gs, ds, r12);
			double sp = r12[10];
			if (Double.isNaN(sp) || Double.isInfinite(sp)) continue;  // DX:405-407
			int li = geo.layerOf(g), w = g - layoff[li];
			double pa = cellPath(a, nn, sp, li, w);
			int i = seg(a, nn, sp, 0);
			if (i < 0 || pa <= MIN_PATH) continue;  // at_s not ok or path <= 0.05 (DX:408-411)
			double d = a[(i + 1) * NC] - a[i * NC];
			double tau = (sp - a[i * NC]) / d;
			double pr = a[i * NC + C_P] * (1 - tau) + a[(i + 1) * NC + C_P] * tau;
			double th = (Math.sqrt(pr * pr + M * M) - M) * 1e3;
			kept[k] = true;
			path[k] = pa;
			s[k] = sp;
			X[k] = (tab.rangeG(StoppingTables.GAS, tv) - tab.rangeG(StoppingTables.GAS, Math.max(th, 1e-3))) / 1.0;  // DX:90-91
			nOn++;
		}
		// calibrated hit quantities (DC:116-173)
		int cnt = 0;
		double sum = 0.0;
		for (int k = 0; k < nh; k++) {
			if (!kept[k]) continue;
			double l = hitLndp(X[k], tv);
			lndp[k] = l;
			double adc = hAdc[k] - offset[geo.layerOf(r.hitWire[k])];
			use[k] = adc > 0 && path[k] > MIN_PATH && !Double.isNaN(l);
			if (use[k]) {
				cnt++;
				sum += l;
			}
		}
		gm = sum / Math.max(cnt, 1);
		double inv = 0.0;
		int n = 0;
		for (int k = 0; k < nh; k++) {
			if (!kept[k]) continue;
			int li = geo.layerOf(r.hitWire[k]);
			double adc = hAdc[k] - offset[li];
			double lnadc = Math.log(Math.max(adc, 1e-3));
			double des = 0.0;  // design row . f3, columns in order (DX:446-455)
			for (int c = 1; c < 8; c++) des += (li == c ? 1.0 : 0.0) * f3[c - 1];
			des += hDoca[k] * f3[7];
			des += (hPed[k] - PED_REF) / 50.0 * f3[8];
			des += Math.log(Math.max(path[k], 0.2) / CHORD) * f3[9];
			double lnG = des - f3Norm;
			lnu[k] = lnadc - lnG - (lndp[k] - gm);
			if (use[k] && !Double.isNaN(lnu[k]) && !Double.isInfinite(lnu[k])) {
				inv += Math.exp(-lnu[k]);
				n++;
			}
		}
		nUsed = n;
		lnE = n == 0 ? Double.NaN : -Math.log(inv / n);  // DX:517-526
		return status;
	}

	/** first segment i with nonzero length and a[i].s <= sp <= a[i+1].s, searched from segment i0 on; -1 if none (prop._at_s, prop.py:653-667). */
	private static int seg(double[] a, int nn, double sp, int i0) {
		for (int i = i0; i < nn - 1; i++) {
			double d = a[(i + 1) * NC] - a[i * NC];
			if (d < 1e-12) continue;
			if (a[i * NC] <= sp && sp <= a[(i + 1) * NC]) return i;
		}
		return -1;
	}

	/** DX:277-285 _cellpath: path of the trajectory within the cell of wire (li, w) over s0 +- 10 mm in steps of 0.04 mm. */
	private double cellPath(double[] a, int nn, double s0, int li, int w) {
		int n = (int) (2 * HALF / DS);
		double tot = 0.0;
		int i0 = 0;  // s_k increases, so a segment before the last match cannot match again
		for (int k = 0; k < n; k++) {
			double sk = s0 - HALF + (k + 0.5) * DS;
			int i = seg(a, nn, sk, i0);
			if (i < 0) continue;
			i0 = i;
			double d = a[(i + 1) * NC] - a[i * NC];
			double tau = (sk - a[i * NC]) / d;
			double t2 = tau * tau, t3 = tau * tau * tau;
			double h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + tau, h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
			int p = i * NC, q = (i + 1) * NC;
			double x = h00 * a[p + 1] + h10 * d * a[p + 4] + h01 * a[q + 1] + h11 * d * a[q + 4];
			double y = h00 * a[p + 2] + h10 * d * a[p + 5] + h01 * a[q + 2] + h11 * d * a[q + 5];
			double z = h00 * a[p + 3] + h10 * d * a[p + 6] + h01 * a[q + 3] + h11 * d * a[q + 6];
			geo.eval(li, w, x, y, z, c4);
			if (c4[0] > 0) tot += DS;
		}
		return tot;
	}

	/** DX:121-134 _hit_shape for one hit: ln Delta_p [keV] at vertex kinetic energy tv and grammage x; NaN if the proton stops within it (T < 0.05 MeV). */
	double hitLndp(double x, double tv) {
		double r0 = tab.rangeG(StoppingTables.GAS, tv) / 1.0;  // DX:78-85 _t_hit
		double rr = r0 - x;
		double th = rr <= 0 ? 0.0 : tab.tFromRange(StoppingTables.GAS, rr * 1.0) / 1.0;
		if (th < 0.05) return Double.NaN;
		double e = th * 1e-3 + M;
		double pp = Math.sqrt(e * e - M * M);
		double beta = pp / e;
		double b2 = beta * beta;  // DX:66-76 _landau_mpv, z2 = 1
		double g2 = 1.0 / (1.0 - b2);
		double xi = 0.153537 * 1.0 * zoa * RHO_GAS * (CHORD / 10.0) / b2;
		double dp = xi * (Math.log(2 * ME_MEV * b2 * g2 / iGas) + Math.log(xi / iGas) + 0.200 - b2);
		return Math.log(dp * 1e3);
	}
}
