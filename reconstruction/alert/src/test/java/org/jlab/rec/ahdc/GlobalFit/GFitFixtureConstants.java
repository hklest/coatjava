package org.jlab.rec.ahdc.GlobalFit;

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

	@Override
	public double[] beamXY(int run) { return new double[] {-1.05, 0.0}; }

	@Override
	public AtofPositionTable atofPosition(int run, int h) { return run == RUN_PER_HALF ? half[h] : pooled; }
}
