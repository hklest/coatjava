package org.jlab.rec.ahdc.GlobalFit;

import java.util.HashMap;
import java.util.Map;

/**
 * Test provider of {@link GFitConstants} from fixture unit/atof_tables.txt.gz: run 23014 serves per half the table fitted on
 * the other half (the Python 23014 fit, atof.load_par with paths [h0, h1]), every other run the pooled table; the beam is
 * (-1.05, 0) mm (plan D8). {@code bankTimes} selects the t_deployed_entry_* time keys (engine), else t_final_entry_* (faithful).
 */
public final class GFitFixtureConstants implements GFitConstants {

	public static final int RUN_PER_HALF = 23014;

	private final AtofPositionTable[] half = new AtofPositionTable[2];
	private final AtofPositionTable pooled;

	public GFitFixtureConstants(boolean bankTimes) {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/atof_tables.txt.gz");
		String calset = bankTimes ? "deployed" : "final";
		AtofPositionTable[] t = {table(fx, "h0", calset), table(fx, "h1", calset)};
		half[0] = t[1];
		half[1] = t[0];
		pooled = table(fx, "hall", calset);
	}

	private static AtofPositionTable table(GFitFixtures.Fx fx, String name, String calset) {
		return AtofPositionTable.fromKeys(k -> fx.d0(name + "." + k), fx.d(name + ".bz_recal"), calset);
	}

	private final Map<String, GFitFixtures.Fx> consts = new HashMap<>();

	/** switch/const_&lt;run&gt;.txt.gz of the run: 23014 -&gt; const_23014, 23015 -&gt; const_23015A. */
	private GFitFixtures.Fx consts(int run) {
		String name = run == 23014 ? "23014" : run == 23015 ? "23015A" : null;
		if (name == null) throw new IllegalArgumentException("no switch constants for run " + run);
		return consts.computeIfAbsent(name, n -> GFitFixtures.resource("switch/const_" + n + ".txt.gz"));
	}

	/** time_keys row h, the row of the half that uses it (stored indexed by the using half, never swapped). */
	@Override
	public SwitchConstants.TimeKeys timeKeys(int run, int h) { return new SwitchConstants.TimeKeys(consts(run).row("time_keys", h)); }

	@Override
	public SwitchConstants.TimeChannels timeChannels(int run) {
		GFitFixtures.Fx fx = consts(run);
		double[][] wx = new double[SwitchConstants.NGROUP][], ww = new double[SwitchConstants.NGROUP][];
		for (int g = 0; g < SwitchConstants.NGROUP; g++) {
			wx[g] = fx.row("walk_x", g);
			ww[g] = fx.row("walk_w", g);
		}
		return new SwitchConstants.TimeChannels(fx.d("t0"), fx.d("p50"), fx.i("group"), wx, ww, fx.d("dep_t0"), fx.d("dep_slope"), fx.d("dep_icpt"));
	}

	private final Map<String, GFitFixtures.Fx> dedxFx = new HashMap<>();

	private GFitFixtures.Fx dedxResource(String name) {
		return dedxFx.computeIfAbsent(name, n -> GFitFixtures.resource("switch/dedx_const_" + n + ".txt.gz"));
	}

	private static double[][] shape2(double[] a, int n0, int n1) {
		double[][] o = new double[n0][n1];
		for (int i = 0; i < n0; i++) System.arraycopy(a, i * n1, o[i], 0, n1);
		return o;
	}

	private static double[][][] shape3(double[] a, int n0, int n1, int n2) {
		double[][][] o = new double[n0][n1][n2];
		for (int i = 0; i < n0; i++)
			for (int j = 0; j < n1; j++) System.arraycopy(a, (i * n1 + j) * n2, o[i][j], 0, n2);
		return o;
	}

	/** dE/dx constants with the all class variance (the CCDB table's set). */
	@Override
	public SwitchConstants.DedxConstants dedx(int run) { return dedx(run, true); }

	/**
	 * dE/dx constants of run 23014 (switch/dedx_const_23014) or 23015 (dedx_const_23015A: its own a, g, ln c, F3 norm and class
	 * variances, every other array those of 23014, dx_transfer_from); {@code classVarAll} selects dx_class_var_all (ge10_all),
	 * else dx_class_var_dedx (ge10_dedx).
	 */
	public SwitchConstants.DedxConstants dedx(int run, boolean classVarAll) {
		if (run != 23014 && run != 23015) throw new IllegalArgumentException("no dedx constants for run " + run);
		GFitFixtures.Fx base = dedxResource("23014"), own = run == 23014 ? base : dedxResource("23015A");
		String cv = classVarAll ? "dx_class_var_all" : "dx_class_var_dedx";
		return new SwitchConstants.DedxConstants(own.d0("dx_a"), own.d0("dx_g"), own.d0("dx_lnc_ref"), own.d0("dx_lnc"), base.d("dx_offset"),
			shape2(base.d("dx_f3"), 2, 10), own.d("dx_f3_norm"), base.d("dx_sig_inst"), shape2(base.d("dx_eff_par"), 2, 3), base.d("dx_toy_bg"),
			shape3(base.d("dx_toy_tau"), 2, 17, 36), shape3(base.d("dx_toy_sig"), 2, 17, 36), base.d("dx_zgrid"),
			shape3(base.d("dx_pdf"), 2, 17, 320), shape3(own.d(cv), 2, 2, 3));
	}

	@Override
	public double[] beamXY(int run) { return new double[] {-1.05, 0.0}; }

	@Override
	public AtofPositionTable atofPosition(int run, int h) { return run == RUN_PER_HALF ? half[h] : pooled; }
}
