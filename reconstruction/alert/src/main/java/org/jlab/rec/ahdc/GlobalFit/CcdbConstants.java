package org.jlab.rec.ahdc.GlobalFit;

import java.util.ArrayList;
import java.util.List;
import org.jlab.utils.groups.IndexedTable;

/**
 * {@link GFitConstants} read from CCDB tables (plan D7, D8). The tables come through a {@link TableSource}, so the engine
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

	public CcdbConstants(TableSource source) {
		this.source = source;
	}

	/** names of the tables of the fit (all four) that the source does not have for {@code run}; empty when complete. */
	public List<String> missing(int run) {
		List<String> m = new ArrayList<>();
		for (String n : BEAM_TABLES) if (source.table(run, n) == null) m.add(n);
		for (String n : ATOF_TABLES) if (source.table(run, n) == null) m.add(n);
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
}
