package org.jlab.rec.ahdc.GlobalFit;

/**
 * Immutable configuration of the global track fit: the fields of the pinned Python {@code Config} (fit.py:53-97) that the
 * ge5_vz_zt variant uses (propagator model, dbox hit term, ATOF term on, fallback off), their packing into the numba
 * arrays {@code _cfarr}/{@code _ciarr} (fit.py:1038-1054, values in fixture unit/config.txt), the grid constants and
 * finite-difference steps of fit.py, the recovery windows of the ge9 freeze file, and the search mode.
 * {@link #pythonDefault()} is the Python configuration; the {@code with*} methods return modified copies. The FAST constants
 * (rule 6, c95) are those of the reference's search measurement (md5
 * 6dd26bbd8cde71c5e6ffe1314cccc51a).
 */
public final class GFitConfig {

	/** Start search: FULL = the Python sequence, FAST = the seeded fit with fallback to FULL. */
	public enum Search { FULL, FAST }

	/** Starting points: seed only (ci[12] = 2, Python "b0") or grid plus seed (ci[12] = 3, "grid+b0"). */
	public enum Starts {
		B0(2), GRID_B0(3);
		public final int code;
		Starts(int code) { this.code = code; }
	}

	/** speed of light factor of the field scale (prop.py SOL), cfg_grid[8]. */
	public static final double SOL = 1.00414;
	/** grid of the start search (cfg_grid: NG, SMAX, DS, NPL, NK, NPAR, NADD, NMISS). */
	public static final int NG = 131;
	public static final double SMAX = 129.75;
	public static final double DS = SMAX / (NG - 1);
	public static final int NPL = 6, NK = 10, NPAR = 15, NADD = 8, NMISS = 40;
	/** forward-difference steps of u0..u4 (cfg_steps_u) and of the kinks. */
	private static final double[] STEPS_U = {1e-3, 1e-5, 1e-6, 1e-3, 1e-5};
	public static final double STEP_KINK = 1e-5;

	// Python Config fields
	private double bz = 5.0 * SOL;
	private double beamX = -1.05, beamY = 0.0;
	private double sigZ = 13.8, sigD0 = 0.3, dz = 0.0;
	private double boxEdge = 0.3, boxEps = 0.02, dboxClip = 2.0;
	private double delta = 0.020, wMin = 1e-2;
	private double[] fMs = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	private double tolCost = 1e-4, tolStep = 1e-3;
	private double bzHelix = 5.0 * SOL;
	private double b0Pscale = 1.0;
	private boolean zOn = true, d0On = true, prefit = true, atofOn = true, atofFallback = false;
	private int nFlip = 2, maxit = 40, npass = 4, th0Mode = 0;
	private Starts starts = Starts.GRID_B0;
	/** hit term: 1 = dbox (ci[1]), the only one ported. */
	private final int hitTerm = 1;
	// recovery windows (ge9_recover.json, fixture cfg rec_windows / rec_pool_w / rec_t_stock)
	private double recTLo = -60.0, recTHi = 240.0, recTotMin = 0.0, recAdcMin = 0.0, recWfMax = 6.0, recKRes = 3.0;
	private double recQ = 1.466250985505345;
	private int recRounds = 2, recPoolW = 3;
	/** hit sigma scale of recovery rows (Config.errscale, fit.py:85: sigma = errscale x kf_error) and the drop-worst retry of recover.Windows.drop_worst. */
	private double errScale = 2.0;
	private boolean recDropWorst = true;
	private double recTStockLo = -60.0, recTStockHi = 240.0;
	// search mode; FAST constants copied from search_measure.json (see the class javadoc)
	private Search search = Search.FULL;
	/** Ladder rule of the FAST search, chosen by measurement: 6 = the full search for every track (no seeded rule passed gate P on both runs). */
	private int fastRule = 6;
	/** 95th percentile of cost per hit of the full search on 23014 (status 0), the threshold of rules 2, 3 and 5. */
	private double c95 = 0.351160986443719;
	/** recovery refits start from the baseline state; false, as the measurement recorded in search_measure.json. */
	private boolean refitSeeded = false;
	/** ATOF time correction (SwitchConstants.TimeChannels, AtofTimeCorrection) and the dE/dx term on or off; both off in the Python default. Not part of cf()/ci(). */
	private boolean timeOn = false, dedxOn = false;

	private GFitConfig() {}

	private GFitConfig copy() {
		GFitConfig c = new GFitConfig();
		c.bz = bz; c.beamX = beamX; c.beamY = beamY; c.sigZ = sigZ; c.sigD0 = sigD0; c.dz = dz;
		c.boxEdge = boxEdge; c.boxEps = boxEps; c.dboxClip = dboxClip; c.delta = delta; c.wMin = wMin;
		c.fMs = fMs.clone(); c.tolCost = tolCost; c.tolStep = tolStep; c.bzHelix = bzHelix; c.b0Pscale = b0Pscale;
		c.zOn = zOn; c.d0On = d0On; c.prefit = prefit; c.atofOn = atofOn; c.atofFallback = atofFallback;
		c.nFlip = nFlip; c.maxit = maxit; c.npass = npass; c.th0Mode = th0Mode; c.starts = starts;
		c.recTLo = recTLo; c.recTHi = recTHi; c.recTotMin = recTotMin; c.recAdcMin = recAdcMin; c.recWfMax = recWfMax;
		c.recKRes = recKRes; c.errScale = errScale; c.recDropWorst = recDropWorst; c.recQ = recQ; c.recRounds = recRounds; c.recPoolW = recPoolW;
		c.recTStockLo = recTStockLo; c.recTStockHi = recTStockHi;
		c.search = search; c.fastRule = fastRule; c.c95 = c95; c.refitSeeded = refitSeeded;
		c.timeOn = timeOn; c.dedxOn = dedxOn;
		return c;
	}

	/** the Python Config of variant ge5_vz_zt (full search). */
	public static GFitConfig pythonDefault() { return new GFitConfig(); }

	public GFitConfig withSearch(Search s) { GFitConfig c = copy(); c.search = s; return c; }

	/** ladder rule 1-6 used by the FAST search (tests exercise the other rules). */
	public GFitConfig withFastRule(int rule) { GFitConfig c = copy(); c.fastRule = rule; return c; }

	public GFitConfig withNFlip(int n) { GFitConfig c = copy(); c.nFlip = n; return c; }

	public GFitConfig withStarts(Starts s) { GFitConfig c = copy(); c.starts = s; return c; }

	/** beam position (x, y) in mm, AHDC frame: the CCDB value of the run in place of the Python default (-1.05, 0). */
	public GFitConfig withBeam(double x, double y) { GFitConfig c = copy(); c.beamX = x; c.beamY = y; return c; }

	/** electron-vertex z prior on or off (off for an event without an electron). */
	public GFitConfig withZPrior(boolean on) { GFitConfig c = copy(); c.zOn = on; return c; }

	/** ATOF time correction switch (time_offsets, walk and T0 tables, SwitchConstants) and dE/dx term switch. */
	public GFitConfig withSwitches(boolean time, boolean dedx) { GFitConfig c = copy(); c.timeOn = time; c.dedxOn = dedx; return c; }

	public double bz() { return bz; }
	public double beamX() { return beamX; }
	public double beamY() { return beamY; }
	public double sigZ() { return sigZ; }
	public double sigD0() { return sigD0; }
	public double dz() { return dz; }
	public double boxEdge() { return boxEdge; }
	public double boxEps() { return boxEps; }
	public double dboxClip() { return dboxClip; }
	public double delta() { return delta; }
	public double wMin() { return wMin; }
	public double fMs(int plane) { return fMs[plane]; }
	public double tolCost() { return tolCost; }
	public double tolStep() { return tolStep; }
	public double bzHelix() { return bzHelix; }
	public double b0Pscale() { return b0Pscale; }
	public boolean zOn() { return zOn; }
	public boolean d0On() { return d0On; }
	public boolean prefit() { return prefit; }
	public boolean atofOn() { return atofOn; }
	public boolean atofFallback() { return atofFallback; }
	public int nFlip() { return nFlip; }
	public int maxit() { return maxit; }
	public int npass() { return npass; }
	public int th0Mode() { return th0Mode; }
	public Starts starts() { return starts; }
	public Search search() { return search; }
	public int fastRule() { return fastRule; }
	public double c95() { return c95; }
	public boolean refitSeeded() { return refitSeeded; }
	public boolean timeOn() { return timeOn; }
	public boolean dedxOn() { return dedxOn; }
	public double recTLo() { return recTLo; }
	public double recTHi() { return recTHi; }
	public double recTotMin() { return recTotMin; }
	public double recAdcMin() { return recAdcMin; }
	public double recWfMax() { return recWfMax; }
	public double recKRes() { return recKRes; }
	public double recQ() { return recQ; }
	public double errScale() { return errScale; }
	public boolean recDropWorst() { return recDropWorst; }
	public int recRounds() { return recRounds; }
	public int recPoolW() { return recPoolW; }
	public double recTStockLo() { return recTStockLo; }
	public double recTStockHi() { return recTStockHi; }

	/** forward-difference step of u0..u4 (cfg_steps_u). */
	public static double stepU(int k) { return STEPS_U[k]; }

	/** the float array of fit._cfarr (fit.py:1038-1045), 24 entries. */
	public double[] cf() {
		double[] cf = new double[24];
		cf[0] = bz; cf[1] = beamX; cf[2] = beamY; cf[3] = 1.0; cf[4] = sigZ; cf[5] = sigD0; cf[6] = dz;
		cf[7] = boxEdge; cf[8] = boxEps; cf[9] = dboxClip; cf[10] = delta; cf[11] = wMin;
		System.arraycopy(fMs, 0, cf, 12, 6);
		cf[18] = tolCost; cf[19] = tolStep; cf[20] = bzHelix; cf[21] = b0Pscale;
		return cf;
	}

	/** the integer array of fit._ciarr (fit.py:1047-1054), 16 entries. */
	public int[] ci() {
		int[] ci = new int[16];
		ci[0] = 1; ci[1] = hitTerm; ci[2] = zOn ? 1 : 0; ci[3] = d0On ? 1 : 0; ci[4] = 0; ci[5] = nFlip;
		ci[6] = 0; ci[7] = 0; ci[8] = maxit; ci[9] = npass; ci[10] = th0Mode; ci[11] = prefit ? 1 : 0;
		ci[12] = starts.code; ci[13] = atofOn ? 1 : 0; ci[14] = atofFallback ? 1 : 0;
		return ci;
	}
}
