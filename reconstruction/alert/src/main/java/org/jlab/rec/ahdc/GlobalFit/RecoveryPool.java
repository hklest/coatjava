package org.jlab.rec.ahdc.GlobalFit;

/**
 * Recovery pool from the raw AHDC::adc rows (decision D11): the rows of the event in no track's hit list, turned into time
 * and distance as HitReader does (time = leading edge - t0 - start, distance 0 for a negative time), then cut by the frozen
 * recovery windows (GFitConfig rec*). The raw cuts of HitReader are not applied. Python reference: recover.derive_rows and
 * pass_windows.
 */
public final class RecoveryPool {

	private RecoveryPool() {}

	/** Wire time offset t0 [ns]; layerIdx is 0-7 in the order of CellGeometry.layoff, wire the index within that layer. */
	public interface WireTime {
		double t0(int layerIdx, int wire);
	}

	/** Time-to-distance relation [mm] for a time &ge; 0; same layer and wire convention as {@link WireTime}. */
	public interface TimeToDistance {
		double doca(int layerIdx, int wire, double time);
	}

	/**
	 * Fills the pool fields (pEv, pG, pLi, pSrc, pAdc, pWf, pTime, pDoca, pLet, pTot, pPed) of ev from its raw rows.
	 * A row is kept when inTrack[row] is false, time in [recTLo, recTHi], ToT &ge; recTotMin, ADC &ge; recAdcMin and
	 * waveform type &le; recWfMax. pEv is the row index, pG the global wire ev.eg[row], pLi the layer index, pSrc 1 (derived).
	 *
	 * @param inTrack per raw row: true when the row is in a track's hit list
	 * @return number of rows kept
	 */
	public static int fromRawRows(GFitEvent ev, GFitConfig cfg, boolean[] inTrack, WireTime t0, TimeToDistance t2d) {
		int n = ev.rLet.length, k = 0;
		int[] ie = new int[n], ig = new int[n], il = new int[n], ia = new int[n], iw = new int[n];
		double[] tm = new double[n], dc = new double[n], le = new double[n], to = new double[n], pe = new double[n];
		for (int r = 0; r < n; r++) {
			if (inTrack[r]) continue;
			double t = ev.rLet[r] - t0.t0(ev.rLi[r], ev.rWire[r]) - ev.startTime;
			if (t < cfg.recTLo() || t > cfg.recTHi() || ev.rTot[r] < cfg.recTotMin() || ev.rAdc[r] < cfg.recAdcMin()
					|| ev.rWf[r] > cfg.recWfMax()) continue;
			ie[k] = r; ig[k] = ev.eg[r]; il[k] = ev.rLi[r]; ia[k] = ev.rAdc[r]; iw[k] = ev.rWf[r];
			tm[k] = t; dc[k] = t < 0 ? 0.0 : t2d.doca(ev.rLi[r], ev.rWire[r], t);
			le[k] = ev.rLet[r]; to[k] = ev.rTot[r]; pe[k] = ev.rPed[r];
			k++;
		}
		ev.pEv = java.util.Arrays.copyOf(ie, k); ev.pG = java.util.Arrays.copyOf(ig, k); ev.pLi = java.util.Arrays.copyOf(il, k);
		ev.pAdc = java.util.Arrays.copyOf(ia, k); ev.pWf = java.util.Arrays.copyOf(iw, k);
		ev.pSrc = new int[k]; java.util.Arrays.fill(ev.pSrc, 1);
		ev.pTime = java.util.Arrays.copyOf(tm, k); ev.pDoca = java.util.Arrays.copyOf(dc, k); ev.pLet = java.util.Arrays.copyOf(le, k);
		ev.pTot = java.util.Arrays.copyOf(to, k); ev.pPed = java.util.Arrays.copyOf(pe, k);
		return k;
	}
}
