package org.jlab.rec.ahdc.GlobalFit;

/**
 * Inputs of one track fit as primitive arrays: the B0 hits and the seed of one row of fit.load_tracks (fit.py:971-1035;
 * hits hg/hm/hism/hsig/hh/hcls/hrow, seed = TC[:, 1:5]). Immutable by convention; the arrays are not copied.
 */
public final class GFitTrack {
	public final int trackid;
	/** global wire index, doca [mm], sigma [mm], half pitch [mm], smear class (0 ok, 1 early, 2 very early), layer class, row in event. */
	public final int[] hg;
	public final double[] hm, hsig, hh;
	public final int[] hism, hcls, hrow;
	/** seed of the B0 track: phi [rad], p [GeV], theta [rad], z [mm]. */
	public final double seedPhi, seedP, seedTheta, seedZ;

	public GFitTrack(int trackid, int[] hg, double[] hm, double[] hsig, double[] hh, int[] hism, int[] hcls, int[] hrow,
			double seedPhi, double seedP, double seedTheta, double seedZ) {
		this.trackid = trackid;
		this.hg = hg; this.hm = hm; this.hsig = hsig; this.hh = hh; this.hism = hism; this.hcls = hcls; this.hrow = hrow;
		this.seedPhi = seedPhi; this.seedP = seedP; this.seedTheta = seedTheta; this.seedZ = seedZ;
	}

	public int nHits() { return hg.length; }
}
