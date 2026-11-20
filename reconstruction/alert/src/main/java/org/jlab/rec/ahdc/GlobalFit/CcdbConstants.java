package org.jlab.rec.ahdc.GlobalFit;

import java.util.ArrayList;
import java.util.List;
import org.jlab.utils.groups.IndexedTable;

/**
 * {@link GFitConstants} read from CCDB tables. The tables come through a {@link TableSource}, so the engine
 * adapter passes {@code ConstantsManager::getConstants}; this class is the only one of the package that knows
 * {@link IndexedTable}. The beam position is 10 x (/geometry/beam/position + /geometry/alert/ahdc/beam_offset), cm to mm,
 * index (0, 0, 0). Both halves get the pooled ATOF table of /calibration/alert/atof/gfit_position and gfit_bar_z with the
 * bank time keys ({@code t_deployed_entry_*}, stored as {@code t_off ... t_fbg}).
 * <p>
 * A table the source does not have (null) is reported by {@link #missing}, and by {@link MissingTableException} naming the
 * tables when a value is asked for anyway; the engine adapter calls {@code missing(run)} and skips the event.
 */
public final class CcdbConstants implements GFitConstants {

	public static final String BEAM = "/geometry/beam/position", BEAM_OFFSET = "/geometry/alert/ahdc/beam_offset",
		ATOF_POSITION = "/calibration/alert/atof/gfit_position", ATOF_BAR_Z = "/calibration/alert/atof/gfit_bar_z";

	private static final String[] BEAM_TABLES = {BEAM, BEAM_OFFSET}, ATOF_TABLES = {ATOF_POSITION, ATOF_BAR_Z};

	private static final String[] TIME_TABLES = {TIME_KEYS, TIME_CHANNEL, TIME_WALK, TIME_OFFSETS};

	/** run and table name to the table of that run, null when CCDB has none. */
	@FunctionalInterface
	public interface TableSource {
		IndexedTable table(int run, String name);
	}

	/** a table of the fit is not in CCDB for the run; {@link #tables()} names them. */
	public static final class MissingTableException extends RuntimeException {
		private static final long serialVersionUID = 1L;
		private final transient List<String> tables;

		MissingTableException(int run, List<String> tables) {
			super("CCDB table(s) missing for run " + run + ": " + String.join(", ", tables));
			this.tables = tables;
		}

		public List<String> tables() { return tables; }
	}

	private final TableSource source;
	private final boolean timeOn, dedxOn;

	/** the default fit: switches off, {@link #missing} covers the four tables of the position term only. */
	public CcdbConstants(TableSource source) {
		this(source, false, false);
	}

	/** {@link #missing} also covers the tables of the switched-on terms of {@code cfg}. */
	public CcdbConstants(TableSource source, GFitConfig cfg) {
		this(source, cfg.timeOn(), cfg.dedxOn());
	}

	private CcdbConstants(TableSource source, boolean timeOn, boolean dedxOn) {
		this.source = source;
		this.timeOn = timeOn;
		this.dedxOn = dedxOn;
	}

	/** names of the tables of the fit (the four of the position term, plus those of the switched-on terms) that the source does not have for {@code run}; empty when complete. */
	public List<String> missing(int run) {
		List<String> m = new ArrayList<>();
		for (String n : BEAM_TABLES) if (source.table(run, n) == null) m.add(n);
		for (String n : ATOF_TABLES) if (source.table(run, n) == null) m.add(n);
		if (timeOn) for (String n : TIME_TABLES) if (source.table(run, n) == null) m.add(n);
		if (dedxOn) for (String n : DEDX_TABLES) if (source.table(run, n) == null) m.add(n);
		return m;
	}

	private IndexedTable[] tables(int run, String[] names) {
		IndexedTable[] t = new IndexedTable[names.length];
		List<String> m = new ArrayList<>();
		for (int i = 0; i < names.length; i++) {
			t[i] = source.table(run, names[i]);
			if (t[i] == null) m.add(names[i]);
		}
		if (!m.isEmpty()) throw new MissingTableException(run, m);
		return t;
	}

	@Override
	public double[] beamXY(int run) {
		IndexedTable[] t = tables(run, BEAM_TABLES);
		return new double[] {10.0 * (t[0].getDoubleValue("x_offset", 0, 0, 0) + t[1].getDoubleValue("x_offset", 0, 0, 0)),
			10.0 * (t[0].getDoubleValue("y_offset", 0, 0, 0) + t[1].getDoubleValue("y_offset", 0, 0, 0))};
	}

	@Override
	public AtofPositionTable atofPosition(int run, int half) {
		IndexedTable[] t = tables(run, ATOF_TABLES);
		double[] bz = new double[2 * 2 * AtofPositionTable.NBAR];      // bz_recal flat 2 x 60 x 2: [0, b, 0] offset, [1, b, 0] scale
		for (int b = 0; b < AtofPositionTable.NBAR; b++) {
			bz[2 * b] = t[1].getDoubleValue("offset", b / 4, b % 4, 10);
			bz[2 * AtofPositionTable.NBAR + 2 * b] = t[1].getDoubleValue("scale", b / 4, b % 4, 10);
		}
		return AtofPositionTable.fromKeys(k -> t[0].getDoubleValue(column(k), 0, 0, 0), bz, "deployed");
	}

	/** column of gfit_position of a key of {@link AtofPositionTable#fromKeys}. */
	private static String column(String key) {
		if (key.equals("bz_off_wW")) return "bz_off";
		if (key.startsWith("t_deployed_entry_")) return "t_" + key.substring("t_deployed_entry_".length());
		return key;
	}

	private static final String[] TIME_COLUMNS = {"off", "s1", "s2", "f2", "mlate", "slate", "flate", "fbg"};

	/** gfit_time_keys row (0, 0, half), columns off ... fbg: the row of the half that uses it. */
	@Override
	public SwitchConstants.TimeKeys timeKeys(int run, int half) {
		IndexedTable t = tables(run, new String[] {TIME_KEYS})[0];
		double[] v = new double[TIME_COLUMNS.length];
		for (int i = 0; i < v.length; i++) v[i] = t.getDoubleValue(TIME_COLUMNS[i], 0, 0, half);
		return new SwitchConstants.TimeKeys(v);
	}

	/**
	 * gfit_time_channel (sector 0-14, layer 0-3, component 0-9; t0, p50), gfit_time_walk (0, 0, group; x0..x11, w0..w11) and
	 * time_offsets (sector, layer, component, order 0; t0, extra1 = slope, extra2 = intercept). The gain group follows from
	 * p50 (tcal.group_of: 0 below 1200, 2 above 2800, else 1).
	 */
	@Override
	public SwitchConstants.TimeChannels timeChannels(int run) {
		IndexedTable[] t = tables(run, new String[] {TIME_CHANNEL, TIME_WALK, TIME_OFFSETS});
		int n = SwitchConstants.NCH;
		double[] t0 = new double[n], p50 = new double[n], dt0 = new double[n], slope = new double[n], icpt = new double[n];
		int[] group = new int[n];
		for (int s = 0; s < 15; s++)
			for (int l = 0; l < 4; l++)
				for (int c = 0; c < 10; c++) {
					int ch = SwitchConstants.channel(s, l, c);
					t0[ch] = t[0].getDoubleValue("t0", s, l, c);
					p50[ch] = t[0].getDoubleValue("p50", s, l, c);
					group[ch] = SwitchConstants.groupOf(p50[ch]);
					dt0[ch] = t[2].getDoubleValue("t0", s, l, c, 0);
					slope[ch] = t[2].getDoubleValue("extra1", s, l, c, 0);
					icpt[ch] = t[2].getDoubleValue("extra2", s, l, c, 0);
				}
		double[][] wx = new double[SwitchConstants.NGROUP][SwitchConstants.NNODE], ww = new double[SwitchConstants.NGROUP][SwitchConstants.NNODE];
		for (int g = 0; g < SwitchConstants.NGROUP; g++)
			for (int k = 0; k < SwitchConstants.NNODE; k++) {
				wx[g][k] = t[1].getDoubleValue("x" + k, 0, 0, g);
				ww[g][k] = t[1].getDoubleValue("w" + k, 0, 0, g);
			}
		return new SwitchConstants.TimeChannels(t0, p50, group, wx, ww, dt0, slope, icpt);
	}

	private static final String[] F3_COLUMNS = {"gain_21", "gain_22", "gain_31", "gain_32", "gain_41", "gain_42", "gain_51", "doca_per_mm", "ped_per_50", "path_exp"},
		EFF_COLUMNS = {"eff_lnthr", "eff_width", "eff_plateau"};

	/**
	 * the dE/dx constants of the five {@link GFitConstants#DEDX_TABLES}: gfit_dedx row (0, 0, half) (a, g, lnc_ref, lnc
	 * and the run-wide values are read from half 0; sig_inst, eff_*, F3 and f3_norm per half, the half that uses them, never swapped),
	 * gfit_dedx_offset (0, li, 0), gfit_dedx_toy (half, n, k; bg, tau, sig), gfit_dedx_pdf (half, n, iz) and gfit_dedx_classvar
	 * (half, prompt, class), which holds the class_var_all set. The z grid is the fixed linspace(-6, 10, 321); n = 0 rows are zeros.
	 */
	@Override
	public SwitchConstants.DedxConstants dedx(int run) {
		IndexedTable[] t = tables(run, DEDX_TABLES);
		int nn = SwitchConstants.DX_NN, nbg = SwitchConstants.DX_NBG, nz = SwitchConstants.DX_NZ;
		double[] offset = new double[SwitchConstants.DX_NLAYER], f3Norm = new double[2], sigInst = new double[2], toyBg = new double[nbg], zgrid = new double[nz];
		double[][] f3 = new double[2][SwitchConstants.DX_NF3], effPar = new double[2][EFF_COLUMNS.length];
		double[][][] tau = new double[2][nn][nbg], sig = new double[2][nn][nbg], pdf = new double[2][nn][nz - 1], cv = new double[2][2][3];
		for (int li = 0; li < offset.length; li++) offset[li] = t[1].getDoubleValue("offset", 0, li, 0);
		for (int h = 0; h < 2; h++) {
			f3Norm[h] = t[0].getDoubleValue("f3_norm", 0, 0, h);
			sigInst[h] = t[0].getDoubleValue("sig_inst", 0, 0, h);
			for (int i = 0; i < EFF_COLUMNS.length; i++) effPar[h][i] = t[0].getDoubleValue(EFF_COLUMNS[i], 0, 0, h);
			for (int i = 0; i < F3_COLUMNS.length; i++) f3[h][i] = t[0].getDoubleValue(F3_COLUMNS[i], 0, 0, h);
			for (int n = 0; n < nn; n++) {
				for (int k = 0; k < nbg; k++) {
					tau[h][n][k] = t[2].getDoubleValue("tau", h, n, k);
					sig[h][n][k] = t[2].getDoubleValue("sig", h, n, k);
				}
				for (int iz = 0; iz < nz - 1; iz++) pdf[h][n][iz] = t[3].getDoubleValue("pdf", h, n, iz);
			}
			for (int w = 0; w < 2; w++)
				for (int c = 0; c < 3; c++) cv[h][w][c] = t[4].getDoubleValue("var", h, w, c);
		}
		for (int k = 0; k < nbg; k++) toyBg[k] = t[2].getDoubleValue("bg", 0, nn - 1, k);
		for (int k = 0; k < nz; k++) zgrid[k] = k == nz - 1 ? 10.0 : -6.0 + k * 0.05;
		return new SwitchConstants.DedxConstants(t[0].getDoubleValue("a", 0, 0, 0), t[0].getDoubleValue("g", 0, 0, 0),
			t[0].getDoubleValue("lnc_ref", 0, 0, 0), t[0].getDoubleValue("lnc", 0, 0, 0), offset, f3, f3Norm, sigInst, effPar, toyBg, tau, sig, zgrid, pdf, cv);
	}
}
