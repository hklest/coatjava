package org.jlab.rec.ahdc.GlobalFit;

/**
 * Inputs shared by the tracks of one event: electron vertex ze, table half, start time, the ATOF rows, the AHDC wires of
 * the event and the recovery pool (rows with trackid &lt; 0), as laid out by fit.load_tracks TC (fit.py:1024-1033) and the
 * t500 fixtures. A plain holder: the replay loader or the engine adapter fills the public fields once.
 */
public final class GFitEvent {
	/** electron vertex z [mm] (TC[:,0]). */
	public double ze;
	/** table half 0/1 = tag_event % 2 (TC[:,5]). */
	public int half;
	/** start time [ns]. */
	public double startTime;
	/** key of the run (calibration set selection); 0 when unused. */
	public int runKey;

	/** ATOF rows: position, bar flag, sector, layer, component, bank time, final time. */
	public double[] ax, ay, az, tBank, tFinal;
	public boolean[] isBar;
	public int[] sector, layer, comp;

	/** ATOF candidates as the fit consumes them (atof.attach): row in the event ATOF rows, bar-z corrected z, phi, association weight, bar flag, time - start (final, bank), sector, component. */
	public int[] cRow, cBar, cSec, cComp;
	public double[] cz, cphi, ct, cTime, cTimeBank;
	/** no candidate at all (ATOF_NONE) and the seed intercept record (20). */
	public int fnone;
	public double[] atofSeed;

	/** global wire index of every AHDC row of the event, in row order. */
	public int[] eg;

	/** recovery pool: row in event, wire, layer index, time, doca, source (0 stock, 1 derived), let, tot, adc, wf, ped. */
	public int[] pEv, pG, pLi, pSrc, pAdc, pWf;
	public double[] pTime, pDoca, pLet, pTot, pPed;

	public int nAtof() { return ax == null ? 0 : ax.length; }

	public int nPool() { return pG == null ? 0 : pG.length; }
}
