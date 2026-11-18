package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import org.jlab.utils.groups.IndexedTable;
import org.junit.Test;

public class GFitTimeCorrectionTest {

	/** Time constants of the run read through in-memory CCDB tables (CcdbConstants), rows of switch/time_rows_&lt;rows&gt;.txt.gz. */
	private static void check(int run, String consts, String rows) {
		Map<String, IndexedTable> m = GFitSwitchConstantsTest.timeTables(new GFitFixtureConstants(true), run);
		CcdbConstants c = new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault().withSwitches(true, false));
		SwitchConstants.TimeChannels tc = c.timeChannels(run);
		GFitFixtures.Fx fx = GFitFixtures.resource("switch/time_rows_" + rows + ".txt.gz");
		int n = fx.d("t_bank").length;
		GFitEvent ev = new GFitEvent();
		ev.ax = new double[n];
		ev.tBank = fx.d("t_bank");
		ev.aEnergy = fx.d("energy");
		ev.sector = fx.i("sector");
		ev.layer = fx.i("layer");
		ev.comp = fx.i("comp");
		ev.isBar = new boolean[n];
		int wedges = 0, bars = 0;
		for (int i = 0; i < n; i++) {
			ev.isBar[i] = ev.comp[i] == 10;
			if (ev.isBar[i]) bars++; else wedges++;
		}
		// the test resource holds the first 600 wedge and 60 bar rows of the 2,000 + 200 draw of $FX/switch/time_rows_<run>.txt
		assertEquals(600, wedges);
		assertEquals(60, bars);
		for (long r : fx.l("run")) assertEquals(run, r);
		AtofTimeCorrection.correct(ev, tc, run);
		double[] want = fx.d("t_expect"), bank = fx.d("t_bank");
		double shift = 0.0;
		for (int i = 0; i < n; i++) {
			assertEquals("row " + i, want[i], ev.aTime[i], 1e-9);
			if (ev.isBar[i]) assertEquals("bar row " + i, bank[i], ev.aTime[i], 0.0);
			else shift = Math.max(shift, Math.abs(ev.aTime[i] - bank[i]));
		}
		assertTrue("the correction moves wedge times (" + consts + ")", shift > 0.01);
		double[] again = ev.aTime;
		AtofTimeCorrection.correct(ev, tc, run);
		assertTrue("aTime reused", again == ev.aTime);
	}

	@Test
	public void timeCorrectionMatchesPython() { check(23015, "const_23015A", "23015A"); }

	@Test
	public void timeCorrectionMatches23014() { check(23014, "const_23014", "23014"); }
}
