package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

/** RecoveryPool.fromRawRows: track rows excluded, time = let - t0 - start, negative time gives distance 0, frozen windows. */
public class GFitRecoveryPoolTest {

	@Test
	public void poolFromRawRows() {
		GFitConfig cfg = GFitConfig.pythonDefault();
		double lo = cfg.recTLo(), hi = cfg.recTHi();
		GFitEvent ev = new GFitEvent();
		ev.startTime = 10.0;
		int n = 12;
		// t0 = 5 for every wire, so time = let - 15
		double[] let = {115, 115, 115, 5, 15 + hi + 1, 15 + lo - 1, 215, 215, 215, 15 + hi, 15 + lo, 15};
		ev.rLet = let;
		ev.rLi = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3};
		ev.rWire = new int[] {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		ev.eg = new int[] {40, 51, 62, 73, 84, 95, 106, 117, 128, 139, 140, 141};
		ev.rTot = new double[n]; ev.rAdc = new int[n]; ev.rWf = new int[n]; ev.rPed = new double[n];
		java.util.Arrays.fill(ev.rTot, cfg.recTotMin() + 1);
		java.util.Arrays.fill(ev.rAdc, (int) cfg.recAdcMin() + 1);
		for (int i = 0; i < n; i++) ev.rPed[i] = 100 + i;
		ev.rWf[6] = (int) cfg.recWfMax() + 1;   // dropped: waveform type
		ev.rTot[7] = cfg.recTotMin() - 1;       // dropped: ToT
		ev.rAdc[8] = (int) cfg.recAdcMin() - 1; // dropped: ADC
		boolean[] inTrack = new boolean[n];
		inTrack[0] = true;                      // excluded: in a track

		int k = RecoveryPool.fromRawRows(ev, cfg, inTrack, (li, w) -> 5.0, (li, w, t) -> 100 * li + w + 0.5 * t);
		// kept: 1, 2, 3 (t = -10, distance 0), 9, 10, 11 (t = hi, lo, 0); dropped: 0 track, 4 and 5 time, 6 wf, 7 tot, 8 adc
		assertEquals(6, k);
		assertEquals(6, ev.nPool());
		assertArrayEquals(new int[] {1, 2, 3, 9, 10, 11}, ev.pEv);
		assertArrayEquals(new int[] {51, 62, 73, 139, 140, 141}, ev.pG);
		assertArrayEquals(new int[] {1, 2, 3, 1, 2, 3}, ev.pLi);
		assertArrayEquals(new int[] {1, 1, 1, 1, 1, 1}, ev.pSrc);
		assertArrayEquals(new double[] {100, 100, -10, hi, lo, 0}, ev.pTime, 0);
		assertArrayEquals(new double[] {155, 256, 0, 113 + 0.5 * hi, 0, 315}, ev.pDoca, 1e-12);
		assertEquals(101, ev.pPed[0], 0);
	}
}
