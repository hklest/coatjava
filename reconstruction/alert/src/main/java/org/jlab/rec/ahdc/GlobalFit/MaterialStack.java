package org.jlab.rec.ahdc.GlobalFit;

/**
 * Cylindrical region stack about the AHDC axis for the propagator: port of prop.make_stack("data_best") with its default
 * variants (pinned prop.py:187-242): target H2, 60 um Kapton wall, HECO2 container, He clear space, AHDC gas chunks with the
 * smeared wire shell between r = 30 and 70 mm, the ATOF bar and wedge. The mean layer radii that place the scattering planes
 * are constants equal to the mean wire radii of the pinned wire CSV (layer_radius of the fixture unit/stopping.txt), because
 * the port does not read the CSV. Immutable and thread safe. Scattering plane index of region k is {@code pla[k]} (-1: none);
 * the thin-scatterer segment is {@code seg[k]}.
 */
public final class MaterialStack {

	public static final int NPL = 6;

	// mean wire radii of layers 21 22 31 32 41 42 (mm), from the fixture unit/stopping.txt layer_radius (prop._layer_radii, prop.py:178-184)
	private static final double R21 = 37.47442305218791, R22 = 41.41756974404924;
	private static final double R31 = 47.2729205170592, R32 = 51.21398752227948;
	private static final double R41 = 57.05704187533805, R42 = 60.98490803619638;

	private static MaterialStack cached;

	public final int nreg;
	public final double[] rOut, rho, xk, zoarho, planeR;
	public final int[] mat, nst, seg, pla;
	public final int npl = NPL;

	/** the data_best stack on the proton tables (cached, immutable). */
	public static synchronized MaterialStack dataBest() {
		if (cached == null) cached = new MaterialStack(StoppingTables.proton());
		return cached;
	}

	private final double[][] regs = new double[32][];  // r_in, r_out, mat, rho, seg, plane_after
	private int n;

	private void add(double a, double b, int mat, double rho, int seg, int pl) {
		regs[n++] = new double[] {a, b, mat, rho, seg, pl};
	}

	private MaterialStack(StoppingTables t) {
		double wallUm = 60.0;
		double[] pr = {3.0 + wallUm * 1e-3, 0.5 * (4.0 + 29.99), 0.5 * (R21 + R22), 0.5 * (R31 + R32), 0.5 * (R41 + R42), 78.5};
		double[] cuts = {0.5 * (pr[2] + pr[3]), 0.5 * (pr[3] + pr[4])};
		double g = 5.015e-4;
		add(0.0, 3.0, StoppingTables.H2, 4.69e-4, 0, -1);
		add(3.0, pr[0], StoppingTables.KAPTON, 1.42, 0, 0);
		add(pr[0], 4.0, StoppingTables.HECO2, 4.87e-4, 1, -1);
		add(4.0, pr[1], StoppingTables.HE, 1.664e-4, 1, 1);
		add(pr[1], 29.99, StoppingTables.HE, 1.664e-4, 1, -1);
		// AHDC gas chunks: [29.99, mid23] seg 2 (plane at pr[2]); [mid23, mid34] seg 3; [mid34, 77] seg 4
		double[][] bounds = {{29.99, cuts[0], pr[2], 2}, {cuts[0], cuts[1], pr[3], 3}, {cuts[1], 77.0, pr[4], 4}};
		for (double[] bd : bounds) {
			double a = bd[0], b = bd[1], plr = bd[2];
			int sg = (int) bd[3];
			double[][] pieces = {{a, plr}, {plr, b}};
			for (int k = 0; k < 2; k++) {
				double u = pieces[k][0], v = pieces[k][1];
				// split at the wire shell edge 70 mm
				for (int w = 0; w < 2; w++) {
					double lo = w == 0 ? u : Math.max(u, 70.0), hi = w == 0 ? Math.min(v, 70.0) : v;
					if (hi > lo + 1e-9) {
						boolean pl = (k == 0 && hi == v) || (k == 0 && hi == Math.min(v, 70.0) && v <= 70.0);
						add(lo, hi, w == 0 ? StoppingTables.GASW : StoppingTables.GAS,
							w == 0 ? (g + StoppingTables.RHO_WIRE_EQ) : g, sg, pl ? sg : -1);
					}
				}
			}
		}
		add(77.0, 78.5, StoppingTables.VT, 1.032, 5, 5);
		add(78.5, 80.0, StoppingTables.VT, 1.032, 5, -1);
		add(80.0, 100.0, StoppingTables.VT, 1.032, 5, -1);
		nreg = n;
		rOut = new double[nreg];
		rho = new double[nreg];
		mat = new int[nreg];
		seg = new int[nreg];
		pla = new int[nreg];
		nst = new int[nreg];
		xk = new double[nreg];
		zoarho = new double[nreg];
		for (int k = 0; k < nreg; k++) {
			double[] r = regs[k];
			rOut[k] = r[1];
			rho[k] = r[3];
			mat[k] = (int) r[2];
			seg[k] = (int) r[4];
			pla[k] = (int) r[5];
			double th = r[1] - r[0];
			boolean solid = mat[k] == StoppingTables.KAPTON || mat[k] == StoppingTables.VT || mat[k] == StoppingTables.AL;
			// solid_per_mm 2, thick_step 1 mm, gas_step 4 mm (prop.py:237); scale 1 leaves the counts unchanged
			nst[k] = (solid && th <= 2.0) ? Math.max(1, (int) Math.ceil(th * 2.0))
				: Math.max(1, (int) Math.ceil(th / (solid ? 1.0 : 4.0)));
			xk[k] = 0.1 * rho[k] / t.x0g[mat[k]];
			zoarho[k] = rho[k] * t.zoa[mat[k]];
		}
		planeR = pr;
	}
}
