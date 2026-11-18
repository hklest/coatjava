package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import org.junit.Test;

public class GFitConstantsTest {

	private static final String[] FOUR = {"/geometry/beam/position", "/geometry/alert/ahdc/beam_offset",
		"/calibration/alert/atof/gfit_position", "/calibration/alert/atof/gfit_bar_z"};
	private static final String[] TIME = {"/calibration/alert/ahdc/time_offsets", "/calibration/alert/ahdc/time_to_distance_wire"};

	@Test
	public void requiredTables() {
		GFitConfig cfg = GFitConfig.pythonDefault();
		Map<String, Integer> off = GFitConstants.requiredTables(cfg, false);
		assertEquals(4, off.size());
		for (String t : FOUR) assertEquals(Integer.valueOf(3), off.get(t));
		for (String t : TIME) assertFalse(off.containsKey(t));
		assertEquals(off, GFitConstants.requiredTables(cfg));
		Map<String, Integer> on = GFitConstants.requiredTables(cfg, true);
		assertEquals(6, on.size());
		for (String t : FOUR) assertTrue(on.containsKey(t));
		for (String t : TIME) assertEquals(Integer.valueOf(3), on.get(t));
	}

	@Test
	public void fixtureProviderMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/atof_tables.txt.gz");
		GFitFixtureConstants c = new GFitFixtureConstants(false);
		String[] names = {"h0", "h1"};
		for (int h = 0; h < 2; h++) {
			AtofPositionTable t = AtofPositionTable.fromKeys(k -> fx.d0(names[1 - 0 * 0] + "." + k), fx.d("h1.bz_recal"), "final");
			assertEquals(34, t.par().length);
			AtofPositionTable p = c.atofPosition(23014, h);
			assertArrayEquals("par_23014 row " + h, fx.row("par_23014", h), p.par(), 0.0);
			assertArrayEquals("offset " + h, slice(fx.d("recal_23014"), h, 0), p.recalOffset(), 0.0);
			assertArrayEquals("scale " + h, slice(fx.d("recal_23014"), h, 1), p.recalScale(), 0.0);
			AtofPositionTable q = c.atofPosition(23015, h);
			assertArrayEquals("par_23015A row " + h, fx.row("par_23015A", h), q.par(), 0.0);
			assertArrayEquals(slice(fx.d("recal_23015A"), h, 0), q.recalOffset(), 0.0);
			assertArrayEquals(slice(fx.d("recal_23015A"), h, 1), q.recalScale(), 0.0);
		}
		assertArrayEquals(new double[] {-1.05, 0.0}, c.beamXY(23014), 0.0);
		assertArrayEquals(new double[] {-1.05, 0.0}, c.beamXY(23015), 0.0);
	}

	@Test
	public void bankTimeKeys() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/atof_tables.txt.gz");
		GFitFixtureConstants fin = new GFitFixtureConstants(false), bank = new GFitFixtureConstants(true);
		String[] keys = {"off", "s1", "s2", "f2", "mlate", "slate", "flate", "fbg"};
		for (int run : new int[] {23014, 23015}) {
			for (int h = 0; h < 2; h++) {
				double[] f = fin.atofPosition(run, h).par(), b = bank.atofPosition(run, h).par();
				// 23014 half h uses the table of half 1-h; the pooled run uses hall
				String tab = run == 23014 ? "h" + (1 - h) : "hall";
				for (int i = 0; i < AtofPositionTable.NP; i++) {
					if (i >= 16 && i <= 23) assertEquals(fx.d0(tab + ".t_deployed_entry_" + keys[i - 16]), b[i], 0.0);
					else assertEquals("slot " + i, f[i], b[i], 0.0);
				}
				int differ = 0;
				for (int i = 16; i <= 23; i++) if (f[i] != b[i]) differ++;
				assertTrue(differ > 0);
			}
		}
	}

	/** row {@code h} of a 2 x 2 x 60 array, component {@code k}. */
	private static double[] slice(double[] a, int h, int k) {
		double[] o = new double[60];
		System.arraycopy(a, (h * 2 + k) * 60, o, 0, 60);
		return o;
	}
}
