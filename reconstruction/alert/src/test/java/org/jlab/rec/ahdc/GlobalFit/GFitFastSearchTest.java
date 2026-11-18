package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * FAST search of GlobalTrackFitter on t60: under the measured rule 6 it equals FULL bit for bit and meets the tiers;
 * fallbacks of the seeded fit equal the FULL result.
 */
public class GFitFastSearchTest {

	private static GFitWorkspace ws;
	private static GFitResult[] full, seedB;
	private static final GFitConfig FULL = GFitConfig.pythonDefault();

	@BeforeClass
	public static void load() {
		GFitFullSearchTest.load();
		ws = GFitFullSearchTest.workspace();
		full = GFitFullSearchTest.runAll(GFitFullSearchTest.fitter(FULL), ws);
		seedB = GFitFullSearchTest.runAll(GFitFullSearchTest.fitter(FULL.withStarts(GFitConfig.Starts.B0).withNFlip(0)), ws);
	}

	private static GFitResult[] fast(int rule) {
		GFitConfig c = FULL.withSearch(GFitConfig.Search.FAST).withFastRule(rule);
		return GFitFullSearchTest.runAll(GFitFullSearchTest.fitter(c), ws);
	}

	@Test
	public void fastMatchesComposite() {
		assertEquals("measured rule", 6, FULL.fastRule());
		GFitResult[] f = fast(FULL.fastRule());
		for (int k = 0; k < f.length; k++) assertNull("track " + k, GFitFullSearchTest.diff(full[k], f[k]));
		GFitFullSearchTest.checkTiers("fast rule " + FULL.fastRule(), f);
	}

	@Test
	public void fastFallsBack() {
		GFitResult out = new GFitResult(GFitFullSearchTest.maxHits), ref = new GFitResult(GFitFullSearchTest.maxHits);
		GFitTrack t = GFitFullSearchTest.set.tracks[0];
		GFitEvent e = GFitFullSearchTest.set.events[0];
		double[] pbad = {Double.NaN, Double.POSITIVE_INFINITY, 0.03, -1.0};
		for (int rule : new int[] {6, 1, 3}) {
			GlobalTrackFitter ft = GFitFullSearchTest.fitter(FULL.withSearch(GFitConfig.Search.FAST).withFastRule(rule));
			for (double p : pbad) {
				GFitTrack tb = new GFitTrack(t.trackid, t.hg, t.hm, t.hsig, t.hh, t.hism, t.hcls, t.hrow, t.seedPhi, p, t.seedTheta, t.seedZ);
				GFitFullSearchTest.fitOne(ft, ws, tb, e, out);
				GFitFullSearchTest.fitOne(GFitFullSearchTest.fitter(FULL), ws, tb, e, ref);
				assertNull("rule " + rule + " seed p " + p, GFitFullSearchTest.diff(ref, out));
			}
		}
		// seeded fit with a status != 0 (rule 1: G fires): the result equals the full search
		int n1 = 0;
		GFitResult[] f1 = fast(1);
		for (int k = 0; k < f1.length; k++) {
			if (seedB[k].status != 0) {
				n1++;
				assertNull("track " + k + " seeded status " + seedB[k].status, GFitFullSearchTest.diff(full[k], f1[k]));
			}
		}
		int n1s = 0;
		for (int k = 0; k < f1.length; k++) if (seedB[k].status == 1) n1s++;
		System.out.printf("fallback: %d tracks with seeded status != 0 (%d with status 1) fell back to the full search bit for bit%n", n1, n1s);
		assertTrue("no t60 track has a failing seeded fit", n1 > 0);
		assertTrue("no t60 track has a seeded fit with status 1", n1s > 0);
	}
}
