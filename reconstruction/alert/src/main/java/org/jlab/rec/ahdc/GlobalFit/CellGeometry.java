package org.jlab.rec.ahdc.GlobalFit;

/**
 * Wire table and drift-cell polygons of the AHDC. Ports fit.wire_tables (fit.py:937-951: wu, layer offsets, mean wire
 * radius lrad, half pitch) and cell.cell_geom as built by common/geom.py:89-107 (polygon vertices from the nominal layer
 * radius and the azimuth of each wire end, cell.py:36-41), and cell._eval (cell.py:105-141): signed distance of a point to
 * the cell polygon at z and its gradient. Wires are ordered by (layer index, wire); row of wire w of layer li is layoff[li] + w.
 * Immutable and thread safe.
 */
public final class CellGeometry {

	/** end plate z of the wire origin and of the wire end (cell.py:34) */
	static final double Z0 = -188.0, Z1 = 162.5;
	/** nominal layer radii (cell.py:36) and radial half width of the cell (cell.py:38), mm */
	private static final double[] RAD = {32., 38., 42., 48., 52., 58., 62., 68.};
	private static final double HW = 2.0;
	private static final int NLAYER = 8;
	/** vertex k: radial offset and azimuth offset in units of pi/N (cell.py / geom.py:100) */
	private static final double[] DR = {HW, HW, -HW, -HW, -HW, HW};
	private static final double[] DA = {0, 1, 1, 0, -1, -1};

	public final double[] wo, we, wu, lrad, half, cellarr;
	public final int[] layoff;

	private CellGeometry(double[] wo, double[] we, double[] wu, int[] layoff, double[] lrad, double[] half, double[] cellarr) {
		this.wo = wo;
		this.we = we;
		this.wu = wu;
		this.layoff = layoff;
		this.lrad = lrad;
		this.half = half;
		this.cellarr = cellarr;
	}

	public int nWires() { return wu.length / 3; }

	public int layerOf(int g) {
		int li = 0;
		while (g >= layoff[li + 1]) li++;
		return li;
	}

	/** builds wu, lrad, half and cellarr from the wire end points (n x 3 flat each) and the 9 layer offsets. */
	public static CellGeometry fromWireEnds(double[] wo, double[] we, int[] layoff) {
		int n = wo.length / 3;
		double[] wu = new double[3 * n];
		for (int g = 0; g < n; g++) {
			double dx = we[3 * g] - wo[3 * g], dy = we[3 * g + 1] - wo[3 * g + 1], dz = we[3 * g + 2] - wo[3 * g + 2];
			double nr = Math.sqrt(dx * dx + dy * dy + dz * dz);
			wu[3 * g] = dx / nr;
			wu[3 * g + 1] = dy / nr;
			wu[3 * g + 2] = dz / nr;
		}
		double[] lrad = new double[NLAYER];
		for (int k = 0; k < NLAYER; k++) {
			double sum = 0;
			for (int g = layoff[k]; g < layoff[k + 1]; g++) sum += Math.hypot(wo[3 * g], wo[3 * g + 1]);
			lrad[k] = sum / (layoff[k + 1] - layoff[k]);
		}
		double[] half = new double[n];
		double[] cell = new double[n * 24];
		int[] lo = layoff.clone();
		for (int li = 0; li < NLAYER; li++) {
			int nw = lo[li + 1] - lo[li];
			double d = Math.PI / nw;
			for (int g = lo[li]; g < lo[li + 1]; g++) {
				half[g] = Math.PI * lrad[li] / nw;
				for (int k = 0; k < 6; k++) {
					double dr = DR[k], da = DA[k] * d;
					for (int p = 0; p < 2; p++) {
						double[] pt = p == 0 ? wo : we;
						double a = Math.atan2(pt[3 * g + 1], pt[3 * g]) + da;
						cell[g * 24 + 4 * k + 2 * p] = (RAD[li] + dr) * Math.cos(a);
						cell[g * 24 + 4 * k + 2 * p + 1] = (RAD[li] + dr) * Math.sin(a);
					}
				}
			}
		}
		return new CellGeometry(wo.clone(), we.clone(), wu, lo, lrad, half, cell);
	}

	/**
	 * (s, ds/dx, ds/dy, ds/dz) of the xy signed distance (positive inside) to the cell polygon of wire w of layer li at
	 * (x, y, z) into r4. Port of cell._eval (cell.py:105-141).
	 */
	public void eval(int li, int w, double x, double y, double z, double[] r4) {
		double[] a = cellarr;
		int c = (layoff[li] + w) * 24;
		double t = (z - Z0) / (Z1 - Z0);
		double px = a[c + 20] + (a[c + 22] - a[c + 20]) * t;
		double py = a[c + 21] + (a[c + 23] - a[c + 21]) * t;
		double best = 1e300;
		int bi = 0;
		double bu = 0.0;
		boolean inside = false;
		for (int i = 0; i < 6; i++) {
			double ax = px, ay = py;
			px = a[c + 4 * i] + (a[c + 4 * i + 2] - a[c + 4 * i]) * t;
			py = a[c + 4 * i + 1] + (a[c + 4 * i + 3] - a[c + 4 * i + 1]) * t;
			double ex = px - ax, ey = py - ay;
			double u = ((x - ax) * ex + (y - ay) * ey) / (ex * ex + ey * ey);
			u = Math.min(1.0, Math.max(0.0, u));
			double dx = x - (ax + u * ex), dy = y - (ay + u * ey);
			double d2 = dx * dx + dy * dy;
			if (d2 < best) {
				best = d2;
				bi = i;
				bu = u;
			}
			if ((ay > y) != (py > y) && x < ax + (y - ay) / (py - ay) * ex) inside = !inside;
		}
		int i0 = (bi + 5) % 6;
		int c0 = c + 4 * i0, c1 = c + 4 * bi;
		double cx = a[c0] + (a[c0 + 2] - a[c0]) * t;
		double cy = a[c0 + 1] + (a[c0 + 3] - a[c0 + 1]) * t;
		cx += bu * (a[c1] + (a[c1 + 2] - a[c1]) * t - cx);
		cy += bu * (a[c1 + 1] + (a[c1 + 3] - a[c1 + 1]) * t - cy);
		double d = Math.sqrt(best);
		double sg = inside ? 1.0 : -1.0;
		double dd = Math.max(d, 1e-12);
		double gx = sg * (x - cx) / dd, gy = sg * (y - cy) / dd;
		double vzx = ((1 - bu) * (a[c0 + 2] - a[c0]) + bu * (a[c1 + 2] - a[c1])) / (Z1 - Z0);
		double vzy = ((1 - bu) * (a[c0 + 3] - a[c0 + 1]) + bu * (a[c1 + 3] - a[c1 + 1])) / (Z1 - Z0);
		r4[0] = sg * d;
		r4[1] = gx;
		r4[2] = gy;
		r4[3] = -(gx * vzx + gy * vzy);
	}
}
