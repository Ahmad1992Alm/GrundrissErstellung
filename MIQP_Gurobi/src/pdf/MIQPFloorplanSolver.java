package pdf;
//////////

////////
//////////import java.util.*;
//////////import org.ojalgo.optimisation.ExpressionsBasedModel;
//////////import org.ojalgo.optimisation.Optimisation;
//////////import org.ojalgo.optimisation.Expression;
//////////import org.ojalgo.optimisation.Variable;
//////////
///////////**
////////// * Einfache MIQP-Variante des Grundriss-Solvers. Verwendet das freie
////////// * ojAlgo-Optimierungspaket und modelliert jeden Raum als Rechteck mit
////////// * ganzzahligen Koordinaten (x,y,w,h). L-förmige Räume werden aktuell als
////////// * Rechteck angenähert.
////////// */
//////////public final class MIQPFloorplanSolver {
//////////
//////////	private MIQPFloorplanSolver() {
//////////	}
//////////
//////////	public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {
//////////		var G = p.G;
//////////		int n = G.nodes.size();
//////////		if (n == 0)
//////////			throw new IllegalArgumentException("Graph hat keine Knoten");
//////////
//////////		ExpressionsBasedModel model = new ExpressionsBasedModel();
//////////		Variable[] x = new Variable[n];
//////////		Variable[] y = new Variable[n];
//////////		Variable[] w = new Variable[n];
//////////		Variable[] h = new Variable[n];
//////////
//////////		for (int i = 0; i < n; i++) {
//////////			var id = G.nodes.get(i);
//////////			var ro = p.perRoomOptions.get(id);
//////////			int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
//////////			int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
//////////			x[i] = model.addVariable("x_" + i).lower(0).upper(p.hullWCells).integer(true);
//////////			y[i] = model.addVariable("y_" + i).lower(0).upper(p.hullHCells).integer(true);
//////////			w[i] = model.addVariable("w_" + i).lower(minDim).upper(maxDim).integer(true);
//////////			h[i] = model.addVariable("h_" + i).lower(minDim).upper(maxDim).integer(true);
//////////
//////////			// Hüllbegrenzung
//////////			model.addExpression("bndX_" + i).upper(p.hullWCells).set(x[i], 1).set(w[i], 1);
//////////			model.addExpression("bndY_" + i).upper(p.hullHCells).set(y[i], 1).set(h[i], 1);
//////////
//////////			// Fläche exakt einhalten (Quadratische Gleichung)
////////////			Integer area = p.areaSqmByRoom.get(id);
////////////			if (area != null && p.enforceAreas) {
////////////				Expression areaExpr = model.addExpression("area_" + i).level(area);
////////////				areaExpr.set(w[i], h[i]);
////////////			}
//////////
//////////			Integer areaSqm = p.areaSqmByRoom.get(id);
//////////			if (areaSqm != null && p.enforceAreas) {
//////////				double minAsp = (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect;
//////////				double maxAsp = (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect;
//////////				minDim = (ro != null && ro.minDimCells != null) ? ro.minDimCells : p.minDim;
//////////				maxDim = (ro != null && ro.maxDimCells != null) ? ro.maxDimCells : p.maxDim;
//////////
//////////				var pairs = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim, p.areaTolCells, minAsp, maxAsp);
//////////				if (pairs.isEmpty())
//////////					throw new IllegalStateException("Keine (w,h)-Paare für " + id);
//////////
//////////				// Binärwahl
//////////				Variable[] z = new Variable[pairs.size()];
//////////				for (int k = 0; k < pairs.size(); k++)
//////////					z[k] = model.addVariable("z_" + i + "_" + k).binary();
//////////
//////////				// genau ein Paar
//////////				var sel = model.addExpression("sel_" + i).level(1);
//////////				for (var zk : z)
//////////					sel.set(zk, 1);
//////////
//////////				// w = Sum z_k * w_k, h = Sum z_k * h_k
//////////				var wPick = model.addExpression("wPick_" + i).level(0).set(w[i], 1);
//////////				var hPick = model.addExpression("hPick_" + i).level(0).set(h[i], 1);
//////////				for (int k = 0; k < pairs.size(); k++) {
//////////					wPick.set(z[k], -pairs.get(k)[0]);
//////////					hPick.set(z[k], -pairs.get(k)[1]);
//////////				}
//////////			}
//////////
//////////		}
//////////
////////////		int M = p.hullWCells + p.hullHCells; // Big-M
//////////		int M = p.hullWCells + p.hullHCells + p.maxDim + 10;
//////////
//////////		// Adjazenz und Nicht-Nachbarn
//////////		for (int i = 0; i < n; i++) {
//////////			for (int j = i + 1; j < n; j++) {
//////////				String idI = G.nodes.get(i);
//////////				String idJ = G.nodes.get(j);
//////////				boolean neighbors = G.areNeighbors(idI, idJ);
//////////				if (neighbors) {
//////////					Variable L = model.addVariable("L_" + i + "_" + j).binary();
//////////					Variable R = model.addVariable("R_" + i + "_" + j).binary();
//////////					Variable T = model.addVariable("T_" + i + "_" + j).binary();
//////////					Variable B = model.addVariable("B_" + i + "_" + j).binary();
////////////					model.addExpression("ori_" + i + "_" + j).level(1).set(L, 1).set(R, 1).set(T, 1).set(B, 1);
//////////					model.addExpression("ori_" + i + "_" + j).lower(1).set(L, 1).set(R, 1).set(T, 1).set(B, 1);
//////////
//////////					// L: i rechts an j links
//////////					model.addExpression("Lx1_" + i + "_" + j).upper(M).set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L,
//////////							M);
//////////					model.addExpression("Lx2_" + i + "_" + j).upper(M).set(x[j], 1).set(x[i], -1).set(w[i], -1).set(L,
//////////							M);
//////////					model.addExpression("Ly1_" + i + "_" + j).upper(M - p.minContactCells).set(y[i], 1).set(y[j], -1)
//////////							.set(h[j], -1).set(L, M);
//////////					model.addExpression("Ly2_" + i + "_" + j).upper(M - p.minContactCells).set(y[j], 1).set(y[i], -1)
//////////							.set(h[i], -1).set(L, M);
//////////
//////////					// R: i links an j rechts
//////////					model.addExpression("Rx1_" + i + "_" + j).upper(M).set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R,
//////////							M);
//////////					model.addExpression("Rx2_" + i + "_" + j).upper(M).set(x[i], 1).set(x[j], -1).set(w[j], -1).set(R,
//////////							M);
//////////					model.addExpression("Ry1_" + i + "_" + j).upper(M - p.minContactCells).set(y[i], 1).set(y[j], -1)
//////////							.set(h[j], -1).set(R, M);
//////////					model.addExpression("Ry2_" + i + "_" + j).upper(M - p.minContactCells).set(y[j], 1).set(y[i], -1)
//////////							.set(h[i], -1).set(R, M);
//////////
//////////					// T: i unten an j oben
//////////					model.addExpression("Tx1_" + i + "_" + j).upper(M).set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T,
//////////							M);
//////////					model.addExpression("Tx2_" + i + "_" + j).upper(M).set(y[j], 1).set(y[i], -1).set(h[i], -1).set(T,
//////////							M);
//////////					model.addExpression("Tx3_" + i + "_" + j).upper(M - p.minContactCells).set(x[i], 1).set(x[j], -1)
//////////							.set(w[j], -1).set(T, M);
//////////					model.addExpression("Tx4_" + i + "_" + j).upper(M - p.minContactCells).set(x[j], 1).set(x[i], -1)
//////////							.set(w[i], -1).set(T, M);
//////////
//////////					// B: i oben an j unten
//////////					model.addExpression("Bx1_" + i + "_" + j).upper(M).set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B,
//////////							M);
//////////					model.addExpression("Bx2_" + i + "_" + j).upper(M).set(y[i], 1).set(y[j], -1).set(h[j], -1).set(B,
//////////							M);
//////////					model.addExpression("Bx3_" + i + "_" + j).upper(M - p.minContactCells).set(x[i], 1).set(x[j], -1)
//////////							.set(w[j], -1).set(B, M);
//////////					model.addExpression("Bx4_" + i + "_" + j).upper(M - p.minContactCells).set(x[j], 1).set(x[i], -1)
//////////							.set(w[i], -1).set(B, M);
//////////				} else if (p.forbidUnwantedContacts) {
//////////					Variable left = model.addVariable("left_" + i + "_" + j).binary();
//////////					Variable right = model.addVariable("right_" + i + "_" + j).binary();
//////////					Variable above = model.addVariable("above_" + i + "_" + j).binary();
//////////					Variable below = model.addVariable("below_" + i + "_" + j).binary();
//////////					model.addExpression("sep_" + i + "_" + j).level(1).set(left, 1).set(right, 1).set(above, 1)
//////////							.set(below, 1);
//////////					model.addExpression("sepL_" + i + "_" + j).upper(M).set(x[i], 1).set(w[i], 1).set(x[j], -1)
//////////							.set(left, M);
//////////					model.addExpression("sepR_" + i + "_" + j).upper(M).set(x[j], 1).set(w[j], 1).set(x[i], -1)
//////////							.set(right, M);
//////////					model.addExpression("sepA_" + i + "_" + j).upper(M).set(y[i], 1).set(h[i], 1).set(y[j], -1)
//////////							.set(above, M);
//////////					model.addExpression("sepB_" + i + "_" + j).upper(M).set(y[j], 1).set(h[j], 1).set(y[i], -1)
//////////							.set(below, M);
//////////				}
//////////			}
//////////		}
//////////
//////////		// Bounding Box und Zielfunktion
//////////		Variable maxX = model.addVariable("maxX").lower(0).upper(p.hullWCells).integer(true);
//////////		Variable maxY = model.addVariable("maxY").lower(0).upper(p.hullHCells).integer(true);
//////////		for (int i = 0; i < n; i++) {
//////////			model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
//////////			model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
//////////		}
//////////
////////////		Expression obj = model.addExpression("OBJ").weight(1.0);
////////////		obj.set(maxX, maxX);
////////////		obj.set(maxY, maxY);
//////////		Expression obj = model.addExpression("OBJ").weight(1.0);
//////////		obj.set(maxX, 1.0);
//////////		obj.set(maxY, 1.0);
//////////
//////////		for (int i = 0; i < n; i++) {
//////////			Expression asp = model.addExpression("asp_" + i).weight(0.1);
//////////			asp.set(w[i], w[i]);
//////////			asp.set(h[i], h[i]);
//////////			asp.set(w[i], h[i], -2);
//////////		}
//////////
////////////		// OPTIONAL: "möglichst quadratisch" als L1-Abweichung nur für Räume, die es
////////////		// brauchen
////////////		// (wenn du eine Room-Kategorie hast; sonst einfach weglassen)
////////////		for (int i = 0; i < n; i++) {
////////////			// Beispiel: L1-Strafe |w - h| mit Hilfsvariable d >= |w-h|
////////////			Variable d = model.addVariable("d_" + i).lower(0); // continuous reicht
////////////			model.addExpression("d1_" + i).lower(0).set(d, 1).set(w[i], -1).set(h[i], +1);
////////////			model.addExpression("d2_" + i).lower(0).set(d, 1).set(w[i], +1).set(h[i], -1);
////////////			obj.set(d, 0.1); // Gewicht anpassen oder für "ANY"-Räume weglassen
////////////		}
//////////
//////////		// harte/softe Limits für den Solver
//////////		Optimisation.Options opts = new Optimisation.Options();
//////////		opts.time_abort = 8_000L; // harter Abbruch nach 8s (Millisekunden)
//////////		opts.time_suffice = 3_000L; // „gut genug“ nach 3s (soft stop)
//////////		opts.iterations_abort = 5_000_000; // optional
//////////		opts.iterations_suffice = 1_000_000; // optional
//////////
//////////		// Falls deine Version ein MIP-Gap hat, heißt es je nach Version
//////////		// unterschiedlich:
//////////		// -> diese Helfer setzen es nur, wenn es existiert (sonst no-op):
//////////		setIfExists(opts, "mip_gap", 0.03); // ältere Schreibweise?
//////////		setIfExists(opts, "mipGap", 0.03); // neuere Schreibweise?
//////////		// (optional, je nach Größe)
//////////		// opts.iterations_suffice = 1_000_000;
//////////		// opts.iterations_abort = 5_000_000;
//////////
//////////		Optimisation.Result res = model.minimise();
//////////		if (!res.getState().isFeasible()) {
//////////			return List.of();
//////////		}
//////////
//////////		FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
//////////		sol.maxX = maxX.getValue().intValue();
//////////		sol.maxY = maxY.getValue().intValue();
//////////		for (int i = 0; i < n; i++) {
//////////			FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
//////////			rp.x0 = x[i].getValue().intValue();
//////////			rp.y0 = y[i].getValue().intValue();
//////////			rp.wH = rp.tX = w[i].getValue().intValue();
//////////			rp.tY = rp.hV = h[i].getValue().intValue();
//////////			rp.wBB = rp.wH;
//////////			rp.hBB = rp.tY;
//////////			sol.piecesByNode.put(G.nodes.get(i), rp);
//////////			System.out.println("Sol: ");
//////////
//////////		}
//////////
//////////		return List.of(sol);
//////////	}
//////////
//////////	private static void setIfExists(Object obj, String field, double value) {
//////////		try {
//////////			var f = obj.getClass().getField(field);
//////////			f.setAccessible(true);
//////////			if (f.getType() == double.class)
//////////				f.setDouble(obj, value);
//////////			else if (f.getType() == Double.class)
//////////				f.set(obj, value);
//////////		} catch (NoSuchFieldException | IllegalAccessException ignore) {
//////////		}
//////////	}
//////////
////////////	static List<int[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim, int tolCells,
////////////			double minAsp, double maxAsp) {
////////////		int target = (int) Math.round(areaSqm / (cellM * cellM));
////////////		List<int[]> pairs = new ArrayList<>();
////////////		for (int W = minDim; W <= maxDim; W++) {
////////////			for (int H = minDim; H <= maxDim; H++) {
////////////				int cells = W * H;
////////////				if (Math.abs(cells - target) <= tolCells) {
////////////					double ar = (double) W / (double) H;
////////////					if (ar >= minAsp && ar <= maxAsp)
////////////						pairs.add(new int[] { W, H });
////////////				}
////////////			}
////////////		}
////////////		return pairs;
////////////	}
//////////	static List<int[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim, int tolCells,
//////////			double minAsp, double maxAsp) {
//////////		int target = (int) Math.round(areaSqm / (cellM * cellM));
//////////		List<int[]> pairs = new ArrayList<>();
//////////		for (int W = minDim; W <= maxDim; W++) {
//////////			for (int H = minDim; H <= maxDim; H++) {
//////////				int cells = W * H;
//////////				if (Math.abs(cells - target) <= tolCells) {
//////////					double ar = (double) W / (double) H;
//////////					if (ar >= minAsp && ar <= maxAsp)
//////////						pairs.add(new int[] { W, H });
//////////				}
//////////			}
//////////		}
//////////		// nach Nähe zur Zielzellenzahl, dann nach "quadratisch"
//////////		pairs.sort(Comparator.comparingInt((int[] wh) -> Math.abs(wh[0] * wh[1] - target))
//////////				.thenComparingInt(wh -> Math.abs(wh[0] - wh[1])));
//////////		int K = Math.min(12, pairs.size()); // <= 12 behalten
//////////		return pairs.subList(0, K);
//////////	}
//////////
//////////}
////////
////////import java.util.*;
////////
////////import org.ojalgo.optimisation.ExpressionsBasedModel;
////////import org.ojalgo.optimisation.Optimisation;
////////import org.ojalgo.optimisation.Expression;
////////import org.ojalgo.optimisation.Variable;
////////
/////////**
//////// * Einfache MIQP-Variante des Grundriss-Solvers. Verwendet das freie
//////// * ojAlgo-Optimierungspaket und modelliert jeden Raum als Rechteck mit
//////// * ganzzahligen Koordinaten (x,y,w,h). L-förmige Räume werden aktuell als
//////// * Rechteck angenähert.
//////// */
////////public final class MIQPFloorplanSolver {
////////
////////	private MIQPFloorplanSolver() {
////////	}
////////
////////	public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {
////////		// Eingabeprüfung
////////		var G = p.G;
////////		int n = G.nodes.size();
////////		if (n == 0)
////////			throw new IllegalArgumentException("Graph hat keine Knoten");
////////		if (p.minDim > p.maxDim)
////////			throw new IllegalArgumentException("minDim > maxDim");
////////		if (p.hullWCells <= 0 || p.hullHCells <= 0)
////////			throw new IllegalArgumentException("Ungültige Bounding-Box");
////////		if (p.areaTolCells < 0)
////////			throw new IllegalArgumentException("areaTolCells < 0"); // FIX: Robustheit
////////		if (p.minAspect > p.maxAspect)
////////			throw new IllegalArgumentException("minAspect > maxAspect"); // FIX
////////
////////		// Modellinitialisierung
////////		ExpressionsBasedModel model = new ExpressionsBasedModel();
////////
////////		// HINWEIS: Die alte Debug-Zeile war vermutlich inkompatibel mit
////////		// ojAlgo-Versionen.
////////		// Wenn Debugging nötig ist, bitte entsprechend der verwendeten Version setzen.
////////		// model.options.debug(true);
////////
////////		Variable[] x = new Variable[n];
////////		Variable[] y = new Variable[n];
////////		Variable[] w = new Variable[n];
////////		Variable[] h = new Variable[n];
////////
////////		// FIX: Achsenspezifische, straffere Big-Ms
//////////		final int Mx = p.hullWCells + p.maxDim;
//////////		final int My = p.hullHCells + p.maxDim;
////////		final int Mx = p.hullWCells ;
////////		final int My = p.hullHCells ;
////////		
////////		// Variablen für jeden Raum
////////		for (int i = 0; i < n; i++) {
////////			var id = G.nodes.get(i);
////////			var ro = p.perRoomOptions.get(id);
////////
////////			int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
////////			int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
////////			if (minDim > maxDim)
////////				throw new IllegalArgumentException("minDim > maxDim für Raum " + id);
////////
////////			if (ro != null && ro.minAspect != null && ro.maxAspect != null && ro.minAspect > ro.maxAspect) {
////////				throw new IllegalArgumentException("minAspect > maxAspect für Raum " + id); // FIX
////////			}
////////
////////			x[i] = model.addVariable("x_" + i).lower(0).upper(p.hullWCells).integer(true);
////////			y[i] = model.addVariable("y_" + i).lower(0).upper(p.hullHCells).integer(true);
////////			w[i] = model.addVariable("w_" + i).lower(minDim).upper(maxDim).integer(true);
////////			h[i] = model.addVariable("h_" + i).lower(minDim).upper(maxDim).integer(true);
////////
////////			// Hüllbegrenzung: x + w ≤ hullW, y + h ≤ hullH
////////			model.addExpression("bndX_" + i).upper(p.hullWCells).set(x[i], 1).set(w[i], 1);
////////			model.addExpression("bndY_" + i).upper(p.hullHCells).set(y[i], 1).set(h[i], 1);
////////
////////			// Flächenanforderungen
////////			Integer areaSqm = p.areaSqmByRoom.get(id);
////////			if (areaSqm != null && p.enforceAreas) {
////////				double minAsp = (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect;
////////				double maxAsp = (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect;
////////				var pairs = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim, p.areaTolCells, minAsp, maxAsp);
////////				if (pairs.isEmpty())
////////					throw new IllegalStateException("Keine (w,h)-Paare für " + id);
////////
////////				Variable[] z = new Variable[pairs.size()];
////////				for (int k = 0; k < pairs.size(); k++) {
////////					z[k] = model.addVariable("z_" + i + "_" + k).binary();
////////				}
////////				// Genau ein Paar auswählen
////////				var sel = model.addExpression("sel_" + i).level(1);
////////				for (var zk : z)
////////					sel.set(zk, 1);
////////				// w = Sum z_k * w_k, h = Sum z_k * h_k
////////				var wPick = model.addExpression("wPick_" + i).level(0).set(w[i], 1);
////////				var hPick = model.addExpression("hPick_" + i).level(0).set(h[i], 1);
////////				for (int k = 0; k < pairs.size(); k++) {
////////					wPick.set(z[k], -pairs.get(k)[0]);
////////					hPick.set(z[k], -pairs.get(k)[1]);
////////				}
////////			}
////////		}
////////
////////		// Adjazenz und Nicht-Überlappung
////////		for (int i = 0; i < n; i++) {
////////			for (int j = i + 1; j < n; j++) {
////////				String idI = G.nodes.get(i);
////////				String idJ = G.nodes.get(j);
////////				boolean neighbors = G.areNeighbors(idI, idJ);
////////
////////				// Binäre Variablen für Nicht-Überlappung (immer benötigt)
////////				Variable left = model.addVariable("left_" + i + "_" + j).binary();
////////				Variable right = model.addVariable("right_" + i + "_" + j).binary();
////////				Variable above = model.addVariable("above_" + i + "_" + j).binary();
////////				Variable below = model.addVariable("below_" + i + "_" + j).binary();
////////
////////				// Nicht-Überlappung: Mindestens eine Richtung
////////				model.addExpression("sep_" + i + "_" + j).lower(1).set(left, 1).set(right, 1).set(above, 1).set(below,
////////						1);
////////
////////				// FIX: Achsenspezifische Big-Ms
////////				model.addExpression("sepL_" + i + "_" + j).upper(Mx).set(x[i], 1).set(w[i], 1).set(x[j], -1).set(left,
////////						Mx); // x_i + w_i - x_j ≤ Mx*(1-left)
////////				model.addExpression("sepR_" + i + "_" + j).upper(Mx).set(x[j], 1).set(w[j], 1).set(x[i], -1).set(right,
////////						Mx); // x_j + w_j - x_i ≤ Mx*(1-right)
////////				model.addExpression("sepA_" + i + "_" + j).upper(My).set(y[i], 1).set(h[i], 1).set(y[j], -1).set(above,
////////						My); // y_i + h_i - y_j ≤ My*(1-above)
////////				model.addExpression("sepB_" + i + "_" + j).upper(My).set(y[j], 1).set(h[j], 1).set(y[i], -1).set(below,
////////						My); // y_j + h_j - y_i ≤ My*(1-below)
////////
////////				if (neighbors) {
////////					// Orientierung nur für Nachbarn
////////					Variable L = model.addVariable("L_" + i + "_" + j).binary();
////////					Variable R = model.addVariable("R_" + i + "_" + j).binary();
////////					Variable T = model.addVariable("T_" + i + "_" + j).binary();
////////					Variable B = model.addVariable("B_" + i + "_" + j).binary();
////////
////////					// Mindestens eine Orientierung
////////					model.addExpression("ori_" + i + "_" + j).lower(1).set(L, 1).set(R, 1).set(T, 1).set(B, 1);
////////
////////					// FIX: Konsistenz zwischen Separations- und Orientierungs-Binärs (>= reicht)
////////					model.addExpression("link_L_" + i + "_" + j).lower(0).set(left, 1).set(L, -1);
////////					model.addExpression("link_R_" + i + "_" + j).lower(0).set(right, 1).set(R, -1);
////////					model.addExpression("link_T_" + i + "_" + j).lower(0).set(above, 1).set(T, -1);
////////					model.addExpression("link_B_" + i + "_" + j).lower(0).set(below, 1).set(B, -1);
////////
////////					// --- FIX: Adjazenz-GLEICHHEITEN korrekt via zwei Ungleichungen ---
////////
////////					// L: i links an j rechts → x_i + w_i = x_j (wenn L=1)
////////					model.addExpression("Lx_le_" + i + "_" + j).upper(Mx).set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L,
////////							Mx); // ≤ Mx*(1-L)
////////					model.addExpression("Lx_ge_" + i + "_" + j).upper(Mx).set(x[j], 1).set(x[i], -1).set(w[i], -1)
////////							.set(L, Mx); // -(...) ≤ Mx*(1-L)
////////
////////					// R: i rechts an j links → x_j + w_j = x_i (wenn R=1)
////////					model.addExpression("Rx_le_" + i + "_" + j).upper(Mx).set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R,
////////							Mx);
////////					model.addExpression("Rx_ge_" + i + "_" + j).upper(Mx).set(x[i], 1).set(x[j], -1).set(w[j], -1)
////////							.set(R, Mx);
////////
////////					// T: i oben an j unten → y_i + h_i = y_j (wenn T=1)
////////					model.addExpression("Ty_le_" + i + "_" + j).upper(My).set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T,
////////							My);
////////					model.addExpression("Ty_ge_" + i + "_" + j).upper(My).set(y[j], 1).set(y[i], -1).set(h[i], -1)
////////							.set(T, My);
////////
////////					// B: i unten an j oben → y_j + h_j = y_i (wenn B=1)
////////					model.addExpression("By_le_" + i + "_" + j).upper(My).set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B,
////////							My);
////////					model.addExpression("By_ge_" + i + "_" + j).upper(My).set(y[i], 1).set(y[j], -1).set(h[j], -1)
////////							.set(B, My);
////////
////////					// --- FIX: Überlappung orthogonaler Achse + Mindest-Kontakt ---
////////
////////					// Y-Überlappung für L/R
////////					Variable overlapY = model.addVariable("overlapY_" + i + "_" + j).lower(0)
////////							.upper(Math.min(h[i].getUpperLimit().intValue(), h[j].getUpperLimit().intValue()));
////////
////////					// Für L=1: y_i - y_j ≤ h_j und y_j - y_i ≤ h_i
////////					model.addExpression("Ly_le1_" + i + "_" + j).upper(My).set(y[i], 1).set(y[j], -1).set(h[j], -1)
////////							.set(L, My);
////////					model.addExpression("Ly_le2_" + i + "_" + j).upper(My).set(y[j], 1).set(y[i], -1).set(h[i], -1)
////////							.set(L, My);
////////					// Für R=1 analog
////////					model.addExpression("Ry_le1_" + i + "_" + j).upper(My).set(y[i], 1).set(y[j], -1).set(h[j], -1)
////////							.set(R, My);
////////					model.addExpression("Ry_le2_" + i + "_" + j).upper(My).set(y[j], 1).set(y[i], -1).set(h[i], -1)
////////							.set(R, My);
////////
////////					// Bindung overlapY an tatsächliche Schnittmenge (aktiv wenn L=1 oder R=1)
////////					model.addExpression("overlapY_leA_" + i + "_" + j).upper(My).set(overlapY, 1).set(y[i], 1)
////////							.set(y[j], -1).set(h[j], -1).set(L, My).set(R, My);
////////					model.addExpression("overlapY_leB_" + i + "_" + j).upper(My).set(overlapY, 1).set(y[j], 1)
////////							.set(y[i], -1).set(h[i], -1).set(L, My).set(R, My);
////////
////////					// Mindest-Kontakt: overlapY ≥ minContact * (L + R) (linear, da Konstante *
////////					// Binär)
////////					if (p.minContactCells > 0) {
////////						model.addExpression("overlapY_min_" + i + "_" + j).lower(0).set(overlapY, 1)
////////								.set(L, -p.minContactCells).set(R, -p.minContactCells);
////////					}
////////
////////					// X-Überlappung für T/B
////////					Variable overlapX = model.addVariable("overlapX_" + i + "_" + j).lower(0)
////////							.upper(Math.min(w[i].getUpperLimit().intValue(), w[j].getUpperLimit().intValue()));
////////
////////					// Für T=1: x_i - x_j ≤ w_j und x_j - x_i ≤ w_i
////////					model.addExpression("Tx_le1_" + i + "_" + j).upper(Mx).set(x[i], 1).set(x[j], -1).set(w[j], -1)
////////							.set(T, Mx);
////////					model.addExpression("Tx_le2_" + i + "_" + j).upper(Mx).set(x[j], 1).set(x[i], -1).set(w[i], -1)
////////							.set(T, Mx);
////////					// Für B=1 analog
////////					model.addExpression("Bx_le1_" + i + "_" + j).upper(Mx).set(x[i], 1).set(x[j], -1).set(w[j], -1)
////////							.set(B, Mx);
////////					model.addExpression("Bx_le2_" + i + "_" + j).upper(Mx).set(x[j], 1).set(x[i], -1).set(w[i], -1)
////////							.set(B, Mx);
////////
////////					// Bindung overlapX an tatsächliche Schnittmenge (aktiv wenn T=1 oder B=1)
////////					model.addExpression("overlapX_leA_" + i + "_" + j).upper(Mx).set(overlapX, 1).set(x[i], 1)
////////							.set(x[j], -1).set(w[j], -1).set(T, Mx).set(B, Mx);
////////					model.addExpression("overlapX_leB_" + i + "_" + j).upper(Mx).set(overlapX, 1).set(x[j], 1)
////////							.set(x[i], -1).set(w[i], -1).set(T, Mx).set(B, Mx);
////////
////////					// Mindest-Kontakt: overlapX ≥ minContact * (T + B)
////////					if (p.minContactCells > 0) {
////////						model.addExpression("overlapX_min_" + i + "_" + j).lower(0).set(overlapX, 1)
////////								.set(T, -p.minContactCells).set(B, -p.minContactCells);
////////					}
////////				}
////////			}
////////		}
////////
////////		// Bounding-Box und Zielfunktion
////////		// FIX: maxX/maxY müssen nicht ganzzahlig sein → entspannteres Modell
////////		Variable maxX = model.addVariable("maxX").lower(0).upper(p.hullWCells);
////////		Variable maxY = model.addVariable("maxY").lower(0).upper(p.hullHCells);
////////		for (int i = 0; i < n; i++) {
////////			model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
////////			model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
////////		}
////////
////////		Expression obj = model.addExpression("OBJ").weight(1.0);
////////		obj.set(maxX, maxX, 1.0).set(maxY, maxY, 1.0); // maxX^2 + maxY^2
////////
////////		for (int i = 0; i < n; i++) {
////////			Expression asp = model.addExpression("asp_" + i).weight(0.1);
////////			asp.set(w[i], w[i], 1.0).set(h[i], h[i], 1.0).set(w[i], h[i], -2.0); // (w-h)^2
////////		}
////////
////////		// Optimierung
////////		Optimisation.Result res = model.minimise();
////////		if (!res.getState().isFeasible()) {
////////			System.out.println("Keine feasible Lösung gefunden");
////////			return List.of();
////////		}
////////
////////		// Lösungsextraktion
////////		FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
////////		sol.maxX = maxX.getValue().intValue();
////////		sol.maxY = maxY.getValue().intValue();
////////		for (int i = 0; i < n; i++) {
////////			FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
////////			rp.x0 = x[i].getValue().intValue();
////////			rp.y0 = y[i].getValue().intValue();
////////			rp.wH = rp.tX = w[i].getValue().intValue();
////////			rp.tY = rp.hV = h[i].getValue().intValue();
////////			rp.wBB = rp.wH;
////////			rp.hBB = rp.tY;
////////			sol.piecesByNode.put(G.nodes.get(i), rp);
////////		}
////////		return List.of(sol);
////////	}
////////
////////	static List<int[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim, int tolCells,
////////			double minAsp, double maxAsp) {
////////		if (areaSqm < 0 || cellM <= 0)
////////			throw new IllegalArgumentException("Ungültige Flächenparameter");
////////		if (tolCells < 0)
////////			throw new IllegalArgumentException("tolCells < 0"); // FIX
////////		if (minDim > maxDim)
////////			throw new IllegalArgumentException("minDim > maxDim");
////////		if (minAsp > maxAsp)
////////			throw new IllegalArgumentException("minAsp > maxAsp"); // FIX
////////
////////		int target = (int) Math.round(areaSqm / (cellM * cellM));
////////		List<int[]> pairs = new ArrayList<>();
////////		for (int W = minDim; W <= maxDim; W++) {
////////			for (int H = minDim; H <= maxDim; H++) {
////////				int cells = W * H;
////////				if (Math.abs(cells - target) <= tolCells) {
////////					double ar = (double) W / H;
////////					if (ar >= minAsp && ar <= maxAsp) {
////////						pairs.add(new int[] { W, H });
////////					}
////////				}
////////			}
////////		}
////////		return pairs;
////////	}
////////}
//////import java.util.*;
//////
//////import org.ojalgo.optimisation.ExpressionsBasedModel;
//////import org.ojalgo.optimisation.Optimisation;
//////import org.ojalgo.optimisation.Expression;
//////import org.ojalgo.optimisation.Variable;
//////
///////**
////// * Einfache MIQP-Variante des Grundriss-Solvers. Verwendet ojAlgo und modelliert jeden
////// * Raum als Rechteck mit ganzzahligen Koordinaten (x,y) und Maßen (w,h).
////// * L-förmige Räume werden aktuell als Rechteck angenähert.
////// */
//////public final class MIQPFloorplanSolver {
//////
//////    private MIQPFloorplanSolver() {}
//////
//////    public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {
//////
//////        // ---------- Eingabeprüfung ----------
//////        var G = p.G;
//////        int n = G.nodes.size();
//////        if (n == 0) throw new IllegalArgumentException("Graph hat keine Knoten");
//////        if (p.minDim > p.maxDim) throw new IllegalArgumentException("minDim > maxDim");
//////        if (p.hullWCells <= 0 || p.hullHCells <= 0) throw new IllegalArgumentException("Ungültige Bounding-Box");
//////        if (p.areaTolCells < 0) throw new IllegalArgumentException("areaTolCells < 0");
//////        if (p.minAspect > p.maxAspect) throw new IllegalArgumentException("minAspect > maxAspect");
//////        if (p.borderCells < 0) throw new IllegalArgumentException("borderCells < 0");
//////
//////        // ---------- Modell ----------
//////        ExpressionsBasedModel model = new ExpressionsBasedModel();
//////        // Tipp: Debug/Timelimit nur versionsspezifisch setzen; hier weggelassen.
//////
//////        Variable[] x = new Variable[n];
//////        Variable[] y = new Variable[n];
//////        Variable[] w = new Variable[n];
//////        Variable[] h = new Variable[n];
//////
//////        // FIX: Achsenspezifische straffe Big-M-Werte.
//////        // Wegen x + w ≤ hullWUse und y + h ≤ hullHUse reichen Mx= hullWUse, My= hullHUse
//////        final int hullWUse = Math.max(0, p.hullWCells - 2 * p.borderCells); // nutzbare Breite
//////        final int hullHUse = Math.max(0, p.hullHCells - 2 * p.borderCells); // nutzbare Höhe
//////        if (hullWUse <= 0 || hullHUse <= 0) throw new IllegalArgumentException("Nutzbare Hülle ≤ 0 (borderCells zu groß?)");
//////
//////        final int Mx = hullWUse; // FIX
//////        final int My = hullHUse; // FIX
//////
//////        // ---------- Variablen pro Raum ----------
//////        for (int i = 0; i < n; i++) {
//////            var id = G.nodes.get(i);
//////            var ro = p.perRoomOptions.get(id);
//////
//////            int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
//////            int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
//////            if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim für Raum " + id);
//////
//////            Double minAspLocal = (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect;
//////            Double maxAspLocal = (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect;
//////            if (minAspLocal > maxAspLocal) throw new IllegalArgumentException("minAspect > maxAspect für Raum " + id);
//////
//////            // FIX: x/y beginnen am inneren Rand (borderCells)
//////            x[i] = model.addVariable("x_" + i).lower(p.borderCells).upper(p.hullWCells - p.borderCells).integer(true);
//////            y[i] = model.addVariable("y_" + i).lower(p.borderCells).upper(p.hullHCells - p.borderCells).integer(true);
//////
//////            // FIX: Wenn Fläche mit z-Paaren erzwungen wird, müssen w/h nicht integer sein
//////            Integer areaSqm = p.areaSqmByRoom.get(id);
//////            boolean discreteWH = !(p.enforceAreas && areaSqm != null);
//////
//////            w[i] = model.addVariable("w_" + i).lower(minDim).upper(maxDim).integer(discreteWH);
//////            h[i] = model.addVariable("h_" + i).lower(minDim).upper(maxDim).integer(discreteWH);
//////
//////            // Hüllbegrenzung innerhalb des inneren Rechtecks:
//////            // (x - border) + w ≤ hullWUse  -> x + w ≤ hullWCells - border
//////            model.addExpression("bndX_" + i).upper(p.hullWCells - p.borderCells)
//////                 .set(x[i], 1).set(w[i], 1);
//////            model.addExpression("bndY_" + i).upper(p.hullHCells - p.borderCells)
//////                 .set(y[i], 1).set(h[i], 1);
//////
//////            // Flächen-/Aspektanforderungen via diskreter (w,h)-Auswahl (optional)
//////            if (areaSqm != null && p.enforceAreas) {
//////                var pairs = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim,
//////                                                p.areaTolCells, minAspLocal, maxAspLocal);
//////                if (pairs.isEmpty()) throw new IllegalStateException("Keine (w,h)-Paare für " + id);
//////
//////                Variable[] z = new Variable[pairs.size()];
//////                for (int k = 0; k < pairs.size(); k++) {
//////                    z[k] = model.addVariable("z_" + i + "_" + k).binary();
//////                }
//////                // genau ein Paar
//////                var sel = model.addExpression("sel_" + i).level(1);
//////                for (var zk : z) sel.set(zk, 1);
//////
//////                // w = Sum z_k * w_k; h = Sum z_k * h_k
//////                var wPick = model.addExpression("wPick_" + i).level(0).set(w[i], 1);
//////                var hPick = model.addExpression("hPick_" + i).level(0).set(h[i], 1);
//////                for (int k = 0; k < pairs.size(); k++) {
//////                    wPick.set(z[k], -pairs.get(k)[0]);
//////                    hPick.set(z[k], -pairs.get(k)[1]);
//////                }
//////            }
//////        }
//////
//////        // ---------- Paarweise Constraints ----------
//////        for (int i = 0; i < n; i++) {
//////            for (int j = i + 1; j < n; j++) {
//////                String idI = G.nodes.get(i);
//////                String idJ = G.nodes.get(j);
//////                boolean neighbors = G.areNeighbors(idI, idJ);
//////
//////                // Binärvariablen für Nicht-Überlappung
//////                Variable left  = model.addVariable("left_"  + i + "_" + j).binary();
//////                Variable right = model.addVariable("right_" + i + "_" + j).binary();
//////                Variable above = model.addVariable("above_" + i + "_" + j).binary();
//////                Variable below = model.addVariable("below_" + i + "_" + j).binary();
//////
//////                // Mindestens eine Richtung aktiv
//////                model.addExpression("sep_" + i + "_" + j).lower(1)
//////                     .set(left, 1).set(right, 1).set(above, 1).set(below, 1);
//////
//////                // FIX: Bei Nicht-Nachbarn optional Kontaktverbot (mind. 1 Zelle Abstand)
//////                int gapX = (!neighbors && p.forbidUnwantedContacts) ? 1 : 0;
//////                int gapY = (!neighbors && p.forbidUnwantedContacts) ? 1 : 0;
//////
//////                // Trennung mit (ggf.) Abstand:
//////                // x_i + w_i - x_j + Mx*left  ≤  Mx - gapX     (aktiv wenn left=1)
//////                model.addExpression("sepL_" + i + "_" + j).upper(Mx - gapX)
//////                     .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(left, Mx);
//////                model.addExpression("sepR_" + i + "_" + j).upper(Mx - gapX)
//////                     .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(right, Mx);
//////                model.addExpression("sepA_" + i + "_" + j).upper(My - gapY)
//////                     .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(above, My);
//////                model.addExpression("sepB_" + i + "_" + j).upper(My - gapY)
//////                     .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(below, My);
//////
//////                // ---------- Adjazenz nur für Nachbarn ----------
//////                if (neighbors) {
//////                    Variable L = model.addVariable("L_" + i + "_" + j).binary();
//////                    Variable R = model.addVariable("R_" + i + "_" + j).binary();
//////                    Variable T = model.addVariable("T_" + i + "_" + j).binary();
//////                    Variable B = model.addVariable("B_" + i + "_" + j).binary();
//////
//////                    // Mindestens/Genau eine Orientierung (abhängig von Option)
//////                    if (p.adjAtLeastOneSide) {
//////                        model.addExpression("ori_" + i + "_" + j).lower(1).set(L,1).set(R,1).set(T,1).set(B,1);
//////                    } else {
//////                        model.addExpression("ori_" + i + "_" + j).level(1).set(L,1).set(R,1).set(T,1).set(B,1);
//////                    }
//////
//////                    // FIX: Konsistenz zu Separations-Flags (L ⇒ left usw.)
//////                    model.addExpression("link_L_" + i + "_" + j).lower(0).set(left, 1).set(L, -1);
//////                    model.addExpression("link_R_" + i + "_" + j).lower(0).set(right, 1).set(R, -1);
//////                    model.addExpression("link_T_" + i + "_" + j).lower(0).set(above, 1).set(T, -1);
//////                    model.addExpression("link_B_" + i + "_" + j).lower(0).set(below, 1).set(B, -1);
//////
//////                    // --- FIX: Gleichheiten über 2 Ungleichungen (aktiv, wenn Binär=1) ---
//////
//////                    // L: x_i + w_i = x_j (wenn L=1)
//////                    model.addExpression("Lx_le_" + i + "_" + j).upper(Mx)
//////                         .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L, Mx);
//////                    model.addExpression("Lx_ge_" + i + "_" + j).upper(Mx)
//////                         .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(L, Mx);
//////
//////                    // R: x_j + w_j = x_i (wenn R=1)
//////                    model.addExpression("Rx_le_" + i + "_" + j).upper(Mx)
//////                         .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R, Mx);
//////                    model.addExpression("Rx_ge_" + i + "_" + j).upper(Mx)
//////                         .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(R, Mx);
//////
//////                    // T: y_i + h_i = y_j (wenn T=1)
//////                    model.addExpression("Ty_le_" + i + "_" + j).upper(My)
//////                         .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T, My);
//////                    model.addExpression("Ty_ge_" + i + "_" + j).upper(My)
//////                         .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(T, My);
//////
//////                    // B: y_j + h_j = y_i (wenn B=1)
//////                    model.addExpression("By_le_" + i + "_" + j).upper(My)
//////                         .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B, My);
//////                    model.addExpression("By_ge_" + i + "_" + j).upper(My)
//////                         .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(B, My);
//////
//////                    // --- FIX: Überlappung in der orthogonalen Achse + Mindestkontakt ---
//////
//////                    // Y-Überlappung (für L oder R)
//////                    Variable overlapY = model.addVariable("overlapY_" + i + "_" + j)
//////                            .lower(0)
//////                            .upper(Math.min(h[i].getUpperLimit().intValue(), h[j].getUpperLimit().intValue()));
//////
//////                    // Für L/R aktiv: |y_i - y_j| ≤ {h_j, h_i}
//////                    model.addExpression("Ly_le1_" + i + "_" + j).upper(My)
//////                         .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(L, My).set(R, My);
//////                    model.addExpression("Ly_le2_" + i + "_" + j).upper(My)
//////                         .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(L, My).set(R, My);
//////
//////                    // overlapY ≤ y_j + h_j - y_i   und   overlapY ≤ y_i + h_i - y_j
//////                    model.addExpression("overlapY_leA_" + i + "_" + j).upper(My)
//////                         .set(overlapY, 1).set(y[i], 1).set(y[j], -1).set(h[j], -1)
//////                         .set(L, My).set(R, My);
//////                    model.addExpression("overlapY_leB_" + i + "_" + j).upper(My)
//////                         .set(overlapY, 1).set(y[j], 1).set(y[i], -1).set(h[i], -1)
//////                         .set(L, My).set(R, My);
//////
//////                    // Mindestkontakt: overlapY ≥ minContact * (L + R)
//////                    if (p.minContactCells > 0) {
//////                        model.addExpression("overlapY_min_" + i + "_" + j).lower(0)
//////                             .set(overlapY, 1)
//////                             .set(L, -p.minContactCells)
//////                             .set(R, -p.minContactCells);
//////                    }
//////
//////                    // X-Überlappung (für T oder B)
//////                    Variable overlapX = model.addVariable("overlapX_" + i + "_" + j)
//////                            .lower(0)
//////                            .upper(Math.min(w[i].getUpperLimit().intValue(), w[j].getUpperLimit().intValue()));
//////
//////                    // Für T/B aktiv: |x_i - x_j| ≤ {w_j, w_i}
//////                    model.addExpression("Tx_le1_" + i + "_" + j).upper(Mx)
//////                         .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(T, Mx).set(B, Mx);
//////                    model.addExpression("Tx_le2_" + i + "_" + j).upper(Mx)
//////                         .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(T, Mx).set(B, Mx);
//////
//////                    // overlapX ≤ x_j + w_j - x_i   und   overlapX ≤ x_i + w_i - x_j
//////                    model.addExpression("overlapX_leA_" + i + "_" + j).upper(Mx)
//////                         .set(overlapX, 1).set(x[i], 1).set(x[j], -1).set(w[j], -1)
//////                         .set(T, Mx).set(B, Mx);
//////                    model.addExpression("overlapX_leB_" + i + "_" + j).upper(Mx)
//////                         .set(overlapX, 1).set(x[j], 1).set(x[i], -1).set(w[i], -1)
//////                         .set(T, Mx).set(B, Mx);
//////
//////                    // Mindestkontakt: overlapX ≥ minContact * (T + B)
//////                    if (p.minContactCells > 0) {
//////                        model.addExpression("overlapX_min_" + i + "_" + j).lower(0)
//////                             .set(overlapX, 1)
//////                             .set(T, -p.minContactCells)
//////                             .set(B, -p.minContactCells);
//////                    }
//////                }
//////            }
//////        }
//////
//////        // ---------- Bounding-Box & Ziel ----------
//////        // FIX: maxX/maxY brauchen nicht integer zu sein (entspannt die LP-Relaxation)
//////        Variable maxX = model.addVariable("maxX").lower(0).upper(p.hullWCells);
//////        Variable maxY = model.addVariable("maxY").lower(0).upper(p.hullHCells);
//////
//////        for (int i = 0; i < n; i++) {
//////            model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
//////            model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
//////        }
//////
//////        // Optional: Hülle vollständig füllen (innen, d. h. bis äußere Hülle)
//////        if (p.forceFillHull) {
//////            model.addExpression("fillX").level(p.hullWCells - p.borderCells); // maxX = W - border
//////            model.getExpression("fillX").set(maxX, 1);
//////            model.addExpression("fillY").level(p.hullHCells - p.borderCells); // maxY = H - border
//////            model.getExpression("fillY").set(maxY, 1);
//////        }
//////
//////        // Zielfunktion: quadratisch wie gehabt (MIQP)
//////        Expression obj = model.addExpression("OBJ").weight(1.0);
//////        obj.set(maxX, maxX, 1.0).set(maxY, maxY, 1.0); // maxX^2 + maxY^2
//////
//////        // Weiche Aspect-Präferenz (w≈h) – ggf. entfernen, falls zu schwer
//////        for (int i = 0; i < n; i++) {
//////            Expression asp = model.addExpression("asp_" + i).weight(0.1);
//////            asp.set(w[i], w[i], 1.0).set(h[i], h[i], 1.0).set(w[i], h[i], -2.0); // (w-h)^2
//////        }
//////
//////        // ---------- Optimierung ----------
//////        Optimisation.Result res = model.minimise();
//////        if (!res.getState().isFeasible()) {
//////            System.out.println("Keine feasible Lösung gefunden – Status: " + res.getState());
//////            return List.of();
//////        }
//////
//////        // ---------- Lösung extrahieren ----------
//////        FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
//////        sol.maxX = maxX.getValue().intValue();
//////        sol.maxY = maxY.getValue().intValue();
//////
//////        for (int i = 0; i < n; i++) {
//////            FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
//////            rp.x0 = x[i].getValue().intValue();
//////            rp.y0 = y[i].getValue().intValue();
//////            rp.wH = rp.tX = w[i].getValue().intValue();
//////            rp.tY = rp.hV = h[i].getValue().intValue();
//////            rp.wBB = rp.wH;
//////            rp.hBB = rp.tY;
//////            sol.piecesByNode.put(G.nodes.get(i), rp);
//////        }
//////        return List.of(sol);
//////    }
//////
//////    // ---------- Hilfsfunktion: erlaubte (w,h)-Paare für Flächenziel ----------
//////    static List<int[]> allowedPairsForArea(int areaSqm,
//////                                           double cellM,
//////                                           int minDim, int maxDim,
//////                                           int tolCells,
//////                                           double minAsp, double maxAsp) {
//////        if (areaSqm < 0 || cellM <= 0) throw new IllegalArgumentException("Ungültige Flächenparameter");
//////        if (tolCells < 0) throw new IllegalArgumentException("tolCells < 0");
//////        if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim");
//////        if (minAsp > maxAsp) throw new IllegalArgumentException("minAsp > maxAsp");
//////
//////        int target = (int) Math.round(areaSqm / (cellM * cellM)); // Zielzellen
//////        List<int[]> pairs = new ArrayList<>();
//////        for (int W = minDim; W <= maxDim; W++) {
//////            for (int H = minDim; H <= maxDim; H++) {
//////                int cells = W * H;
//////                if (Math.abs(cells - target) <= tolCells) {
//////                    double ar = (double) W / H;
//////                    if (ar >= minAsp && ar <= maxAsp) {
//////                        pairs.add(new int[]{W, H});
//////                    }
//////                }
//////            }
//////        }
//////        return pairs;
//////    }
//////}
////import java.util.*;
////
////import org.ojalgo.optimisation.ExpressionsBasedModel;
////import org.ojalgo.optimisation.Optimisation;
////import org.ojalgo.optimisation.Expression;
////import org.ojalgo.optimisation.Variable;
////
/////**
//// * MILP-Variante des Grundriss-Solvers mit straffen Big-Ms und
//// * minimalem Binärsatz. Rechteckräume mit (x,y,w,h).
//// * L-förmige Räume werden aktuell als Rechteck angenähert.
//// */
////public final class MIQPFloorplanSolver {
////
////    private MIQPFloorplanSolver() {}
////
////    public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {
////
////        // --- Eingaben prüfen ---
////        var G = p.G;
////        int n = G.nodes.size();
////        if (n == 0) throw new IllegalArgumentException("Graph hat keine Knoten");
////        if (p.minDim > p.maxDim) throw new IllegalArgumentException("minDim > maxDim");
////        if (p.hullWCells <= 0 || p.hullHCells <= 0) throw new IllegalArgumentException("Ungültige Bounding-Box");
////        if (p.borderCells < 0) throw new IllegalArgumentException("borderCells < 0");
////        if (p.areaTolCells < 0) throw new IllegalArgumentException("areaTolCells < 0");
////        if (p.minAspect > p.maxAspect) throw new IllegalArgumentException("minAspect > maxAspect");
////
////        // --- Nutzbare Innenhülle & straffe Big-Ms ---
////        final int hullWUse = Math.max(0, p.hullWCells - 2 * p.borderCells);
////        final int hullHUse = Math.max(0, p.hullHCells - 2 * p.borderCells);
////        if (hullWUse <= 0 || hullHUse <= 0) {
////            throw new IllegalArgumentException("Nutzbare Hülle ≤ 0 (borderCells zu groß?)");
////        }
////        final int Mx = hullWUse; // straff
////        final int My = hullHUse; // straff
////
////        // --- Modell & Variablen ---
////        ExpressionsBasedModel model = new ExpressionsBasedModel();
////
////        Variable[] x = new Variable[n];
////        Variable[] y = new Variable[n];
////        Variable[] w = new Variable[n];
////        Variable[] h = new Variable[n];
////
////        for (int i = 0; i < n; i++) {
////            String id = G.nodes.get(i);
////            var ro = p.perRoomOptions.get(id);
////
////            int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
////            int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
////            if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim für Raum " + id);
////
////            double minAspLocal = (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect;
////            double maxAspLocal = (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect;
////            if (minAspLocal > maxAspLocal) throw new IllegalArgumentException("minAspect > maxAspect für Raum " + id);
////
////            // Koordinaten am Innenrand
////            x[i] = model.addVariable("x_" + i).lower(p.borderCells).upper(p.hullWCells - p.borderCells).integer(true);
////            y[i] = model.addVariable("y_" + i).lower(p.borderCells).upper(p.hullHCells - p.borderCells).integer(true);
////
////            // w/h: bei enforceAreas (mit Ziel-Fläche) kontinuierlich – z-Variablen diskretisieren ohnehin
////            Integer areaSqm = p.areaSqmByRoom.get(id);
////            boolean makeWHInteger = !(p.enforceAreas && areaSqm != null);
////
////            w[i] = model.addVariable("w_" + i).lower(minDim).upper(maxDim).integer(makeWHInteger);
////            h[i] = model.addVariable("h_" + i).lower(minDim).upper(maxDim).integer(makeWHInteger);
////
////            // Innenhülle: x + w ≤ hullWCells - border, y + h ≤ hullHCells - border
////            Expression bndX = model.addExpression("bndX_" + i).upper(p.hullWCells - p.borderCells);
////            bndX.set(x[i], 1).set(w[i], 1);
////            Expression bndY = model.addExpression("bndY_" + i).upper(p.hullHCells - p.borderCells);
////            bndY.set(y[i], 1).set(h[i], 1);
////
////            // Diskrete (w,h)-Auswahl bei enforceAreas
////            if (areaSqm != null && p.enforceAreas) {
////                var pairs = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim,
////                        p.areaTolCells, minAspLocal, maxAspLocal);
////                if (pairs.isEmpty()) throw new IllegalStateException("Keine (w,h)-Paare für " + id);
////
////                Variable[] z = new Variable[pairs.size()];
////                for (int k = 0; k < pairs.size(); k++) z[k] = model.addVariable("z_" + i + "_" + k).binary();
////
////                // genau ein Paar
////                Expression sel = model.addExpression("sel_" + i).level(1);
////                for (var zk : z) sel.set(zk, 1);
////
////                // w = Sum z_k * w_k; h = Sum z_k * h_k
////                Expression wPick = model.addExpression("wPick_" + i).level(0).set(w[i], 1);
////                Expression hPick = model.addExpression("hPick_" + i).level(0).set(h[i], 1);
////                for (int k = 0; k < pairs.size(); k++) {
////                    wPick.set(z[k], -pairs.get(k)[0]);
////                    hPick.set(z[k], -pairs.get(k)[1]);
////                }
////            }
////        }
////
////        // --- Paarweise Constraints ---
////        for (int i = 0; i < n; i++) {
////            for (int j = i + 1; j < n; j++) {
////                boolean neighbors = G.areNeighbors(G.nodes.get(i), G.nodes.get(j));
////
////                if (neighbors) {
////                    // Nur Adjazenz-Orientierung (L/R/T/B), KEINE sep-Binärs (spart viel Branching)
////                    Variable L = model.addVariable("L_" + i + "_" + j).binary();
////                    Variable R = model.addVariable("R_" + i + "_" + j).binary();
////                    Variable T = model.addVariable("T_" + i + "_" + j).binary();
////                    Variable B = model.addVariable("B_" + i + "_" + j).binary();
////
////                    // Mindestens oder genau 1 Orientierung
////                    if (p.adjAtLeastOneSide) {
////                        model.addExpression("ori_" + i + "_" + j).lower(1).set(L,1).set(R,1).set(T,1).set(B,1);
////                    } else {
////                        model.addExpression("ori_" + i + "_" + j).level(1).set(L,1).set(R,1).set(T,1).set(B,1);
////                    }
////
////                    // L: x_i + w_i = x_j (wenn L=1) – zwei Ungleichungen
////                    model.addExpression("Lx_le_" + i + "_" + j).upper(Mx)
////                         .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L, Mx);
////                    model.addExpression("Lx_ge_" + i + "_" + j).upper(Mx)
////                         .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(L, Mx);
////
////                    // R: x_j + w_j = x_i (wenn R=1)
////                    model.addExpression("Rx_le_" + i + "_" + j).upper(Mx)
////                         .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R, Mx);
////                    model.addExpression("Rx_ge_" + i + "_" + j).upper(Mx)
////                         .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(R, Mx);
////
////                    // T: y_i + h_i = y_j (wenn T=1)
////                    model.addExpression("Ty_le_" + i + "_" + j).upper(My)
////                         .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T, My);
////                    model.addExpression("Ty_ge_" + i + "_" + j).upper(My)
////                         .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(T, My);
////
////                    // B: y_j + h_j = y_i (wenn B=1)
////                    model.addExpression("By_le_" + i + "_" + j).upper(My)
////                         .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B, My);
////                    model.addExpression("By_ge_" + i + "_" + j).upper(My)
////                         .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(B, My);
////
////                    // Orthogonale Überlappung: für L/R → Y-Überlappung, für T/B → X-Überlappung
////                    // (keine Mindestkontaktlänge, weil minContactCells=0 schnell machen soll)
////                    // Y-Überlappung aktiv, wenn L oder R
////                    model.addExpression("Yov1_" + i + "_" + j).upper(My)
////                         .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(L, My).set(R, My); // y_i - y_j ≤ h_j
////                    model.addExpression("Yov2_" + i + "_" + j).upper(My)
////                         .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(L, My).set(R, My); // y_j - y_i ≤ h_i
////
////                    // X-Überlappung aktiv, wenn T oder B
////                    model.addExpression("Xov1_" + i + "_" + j).upper(Mx)
////                         .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(T, Mx).set(B, Mx); // x_i - x_j ≤ w_j
////                    model.addExpression("Xov2_" + i + "_" + j).upper(Mx)
////                         .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(T, Mx).set(B, Mx); // x_j - x_i ≤ w_i
////                    
////                    
////                    
////                    
////                    
////                    
////                    
////                    
////                    
////                    
////                 // Gemeinsame Kantenlänge modellieren:
////                    Variable overlapY = model.addVariable("overlapY_" + i + "_" + j)
////                            .lower(0)
////                            .upper(Math.min(h[i].getUpperLimit().intValue(),
////                                            h[j].getUpperLimit().intValue()));
////
////                    Variable overlapX = model.addVariable("overlapX_" + i + "_" + j)
////                            .lower(0)
////                            .upper(Math.min(w[i].getUpperLimit().intValue(),
////                                            w[j].getUpperLimit().intValue()));
////
////                    // --- Y-Überlappung, aktiv wenn L oder R (horizontale Adjazenz) ---
////                    // overlapY ≤ y_j + h_j - y_i
////                    model.addExpression("ovY_leA_" + i + "_" + j).upper(My)
////                         .set(overlapY, 1).set(y[i], 1).set(y[j], -1).set(h[j], -1)
////                         .set(L, My).set(R, My);
////                    // overlapY ≤ y_i + h_i - y_j
////                    model.addExpression("ovY_leB_" + i + "_" + j).upper(My)
////                         .set(overlapY, 1).set(y[j], 1).set(y[i], -1).set(h[i], -1)
////                         .set(L, My).set(R, My);
////
////                    // Mindestkontakt in Y: overlapY ≥ minContactCells * (L + R)
////                    if (p.minContactCells > 0) {
////                        model.addExpression("ovY_min_" + i + "_" + j).lower(0)
////                             .set(overlapY, 1)
////                             .set(L, -p.minContactCells)
////                             .set(R, -p.minContactCells);
////                    }
////
////                    // --- X-Überlappung, aktiv wenn T oder B (vertikale Adjazenz) ---
////                    // overlapX ≤ x_j + w_j - x_i
////                    model.addExpression("ovX_leA_" + i + "_" + j).upper(Mx)
////                         .set(overlapX, 1).set(x[i], 1).set(x[j], -1).set(w[j], -1)
////                         .set(T, Mx).set(B, Mx);
////                    // overlapX ≤ x_i + w_i - x_j
////                    model.addExpression("ovX_leB_" + i + "_" + j).upper(Mx)
////                         .set(overlapX, 1).set(x[j], 1).set(x[i], -1).set(w[i], -1)
////                         .set(T, Mx).set(B, Mx);
////
////                    // Mindestkontakt in X: overlapX ≥ minContactCells * (T + B)
////                    if (p.minContactCells > 0) {
////                        model.addExpression("ovX_min_" + i + "_" + j).lower(0)
////                             .set(overlapX, 1)
////                             .set(T, -p.minContactCells)
////                             .set(B, -p.minContactCells);
////                    }
////                    
////                    
////                    
////                    
////                    
////                    
////                    
////                    
////
////                } else {
////                    // Nicht-Nachbarn: klassische Nicht-Überlappung mit 4 sep-Binärs
////                    Variable left  = model.addVariable("left_"  + i + "_" + j).binary();
////                    Variable right = model.addVariable("right_" + i + "_" + j).binary();
////                    Variable above = model.addVariable("above_" + i + "_" + j).binary();
////                    Variable below = model.addVariable("below_" + i + "_" + j).binary();
////
////                    // Mindestens eine Richtung aktiv
////                    model.addExpression("sep_" + i + "_" + j).lower(1)
////                         .set(left, 1).set(right, 1).set(above, 1).set(below, 1);
////
////                    // ggf. 1-Zellen-Abstand (Berührungsverbot)
////                    int gapX = p.forbidUnwantedContacts ? 1 : 0;
////                    int gapY = p.forbidUnwantedContacts ? 1 : 0;
////
////                    model.addExpression("sepL_" + i + "_" + j).upper(Mx - gapX)
////                         .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(left, Mx);   // x_i + w_i ≤ x_j - gapX
////                    model.addExpression("sepR_" + i + "_" + j).upper(Mx - gapX)
////                         .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(right, Mx);  // x_j + w_j ≤ x_i - gapX
////                    model.addExpression("sepA_" + i + "_" + j).upper(My - gapY)
////                         .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(above, My);  // y_i + h_i ≤ y_j - gapY
////                    model.addExpression("sepB_" + i + "_" + j).upper(My - gapY)
////                         .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(below, My);  // y_j + h_j ≤ y_i - gapY
////                }
////            }
////        }
////
////        // --- Bounding-Box (rechts/oben) & Ziel ---
////        Variable maxX = model.addVariable("maxX").lower(0).upper(p.hullWCells);
////        Variable maxY = model.addVariable("maxY").lower(0).upper(p.hullHCells);
////
////        for (int i = 0; i < n; i++) {
////            model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
////            model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
////        }
////
////        if (p.forceFillHull) {
////            Expression fillX = model.addExpression("fillX").level(p.hullWCells - p.borderCells);
////            fillX.set(maxX, 1);
////            Expression fillY = model.addExpression("fillY").level(p.hullHCells - p.borderCells);
////            fillY.set(maxY, 1);
////        }
////
////        // **LINEARE** Zielfunktion (MILP statt MIQP): minimiert Ausdehnung
////        Expression obj = model.addExpression("OBJ").weight(1.0);
////        obj.set(maxX, 1.0).set(maxY, 1.0);
////
////        // --- Lösen ---
////        Optimisation.Result res = model.minimise();
////        if (!res.getState().isFeasible()) {
////            System.out.println("Keine feasible Lösung gefunden – Status: " + res.getState());
////            return List.of();
////        }
////
////        // --- Lösung extrahieren ---
////        FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
////        sol.maxX = maxX.getValue().intValue();
////        sol.maxY = maxY.getValue().intValue();
////
////        for (int i = 0; i < n; i++) {
////            FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
////            rp.x0 = x[i].getValue().intValue();
////            rp.y0 = y[i].getValue().intValue();
////            rp.wH = rp.tX = w[i].getValue().intValue();
////            rp.tY = rp.hV = h[i].getValue().intValue();
////            rp.wBB = rp.wH;
////            rp.hBB = rp.tY;
////            sol.piecesByNode.put(G.nodes.get(i), rp);
////        }
////        return List.of(sol);
////    }
////
////    // Erlaubte (w,h)-Paare für Flächenziel (in Zellen, via cellMeters)
////    static List<int[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim,
////                                           int tolCells, double minAsp, double maxAsp) {
////        if (areaSqm < 0 || cellM <= 0) throw new IllegalArgumentException("Ungültige Flächenparameter");
////        if (tolCells < 0) throw new IllegalArgumentException("tolCells < 0");
////        if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim");
////        if (minAsp > maxAsp) throw new IllegalArgumentException("minAsp > maxAsp");
////
////        int target = (int) Math.round(areaSqm / (cellM * cellM));
////        List<int[]> pairs = new ArrayList<>();
////        for (int W = minDim; W <= maxDim; W++) {
////            for (int H = minDim; H <= maxDim; H++) {
////                int cells = W * H;
////                if (Math.abs(cells - target) <= tolCells) {
////                    double ar = (double) W / H;
////                    if (ar >= minAsp && ar <= maxAsp) {
////                        pairs.add(new int[]{W, H});
////                    }
////                }
////            }
////        }
////        return pairs;
////    }
////}
//import java.util.*;
//import org.ojalgo.optimisation.*;
//
//public final class MIQPFloorplanSolver {
//
//    private MIQPFloorplanSolver() {}
//
//    // wie viele unterschiedliche gute Lösungen sollen geliefert werden?
//    private static final int K_BEST = 5; // bei Bedarf 1 setzen
//
//    public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {
//
//        // --- Vorab-Checks ---
//        var G = p.G;
//        int n = G.nodes.size();
//        if (n == 0) throw new IllegalArgumentException("Graph hat keine Knoten");
//        if (p.minDim > p.maxDim) throw new IllegalArgumentException("minDim > maxDim");
//        if (p.hullWCells <= 0 || p.hullHCells <= 0) throw new IllegalArgumentException("Ungültige Hülle");
//        if (p.borderCells < 0) throw new IllegalArgumentException("borderCells < 0");
//        if (p.areaTolCells < 0) throw new IllegalArgumentException("areaTolCells < 0");
//        if (p.minAspect > p.maxAspect) throw new IllegalArgumentException("minAspect > maxAspect");
//
//        // Innenhülle + straffe Big-Ms
//        final int hullWUse = Math.max(0, p.hullWCells - 2 * p.borderCells);
//        final int hullHUse = Math.max(0, p.hullHCells - 2 * p.borderCells);
//        if (hullWUse <= 0 || hullHUse <= 0)
//            throw new IllegalArgumentException("Nutzbare Hülle ≤ 0 (borderCells zu groß?)");
//        final int Mx = hullWUse, My = hullHUse;
//
//        // Liste der No-Good-Muster (über alle Binärvariablen); jedes Muster ist ein boolean[]
//        List<boolean[]> noGoods = new ArrayList<>();
//        List<FloorplanCrossModalSolver.Solution> out = new ArrayList<>();
//
//        for (int ksol = 0; ksol < K_BEST; ksol++) {
//            // --- Modell neu aufbauen (deterministische Reihenfolge!) ---
//            ExpressionsBasedModel model = new ExpressionsBasedModel();
//
//            Variable[] x = new Variable[n];
//            Variable[] y = new Variable[n];
//            Variable[] w = new Variable[n];
//            Variable[] h = new Variable[n];
//
//            // wir sammeln alle Binärvariablen in genau der Reihenfolge, in der wir sie anlegen
//            final List<Variable> binVars = new ArrayList<>();
//
//            // Overlap-Variablen (nur für Nachbarn) sammeln, um sie in Stufe 2 zu gewichten
//            final List<Variable> overlaps = new ArrayList<>();
//
//            // ---- Räume & Bounds ----
//            for (int i = 0; i < n; i++) {
//                String id = G.nodes.get(i);
//                var ro = p.perRoomOptions.get(id);
//
//                int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
//                int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
//                if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim für Raum " + id);
//
//                double minAspLocal = (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect;
//                double maxAspLocal = (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect;
//                if (minAspLocal > maxAspLocal) throw new IllegalArgumentException("minAspect > maxAspect für Raum " + id);
//
//                x[i] = model.addVariable("x_" + i).lower(p.borderCells).upper(p.hullWCells - p.borderCells).integer(true);
//                y[i] = model.addVariable("y_" + i).lower(p.borderCells).upper(p.hullHCells - p.borderCells).integer(true);
//
//                Integer areaSqm = p.areaSqmByRoom.get(id);
//                boolean whInteger = !(p.enforceAreas && areaSqm != null); // bei enforceAreas→z: w/h können kontinuierlich
//                w[i] = model.addVariable("w_" + i).lower(minDim).upper(maxDim).integer(whInteger);
//                h[i] = model.addVariable("h_" + i).lower(minDim).upper(maxDim).integer(whInteger);
//
//                model.addExpression("bndX_" + i).upper(p.hullWCells - p.borderCells).set(x[i], 1).set(w[i], 1);
//                model.addExpression("bndY_" + i).upper(p.hullHCells - p.borderCells).set(y[i], 1).set(h[i], 1);
//
//                // diskrete (w,h)-Paare bei enforceAreas
//                if (areaSqm != null && p.enforceAreas) {
//                    var pairs = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim,
//                            p.areaTolCells, minAspLocal, maxAspLocal);
//                    if (pairs.isEmpty()) return List.of(); // infeasible: kein (w,h)
//
//                    Variable[] z = new Variable[pairs.size()];
//                    for (int kk = 0; kk < pairs.size(); kk++) {
//                        z[kk] = model.addVariable("z_" + i + "_" + kk).binary();
//                        binVars.add(z[kk]);
//                    }
//                    var sel = model.addExpression("sel_" + i).level(1);
//                    for (var zk : z) sel.set(zk, 1);
//                    var wPick = model.addExpression("wPick_" + i).level(0).set(w[i], 1);
//                    var hPick = model.addExpression("hPick_" + i).level(0).set(h[i], 1);
//                    for (int kk = 0; kk < pairs.size(); kk++) {
//                        wPick.set(z[kk], -pairs.get(kk)[0]);
//                        hPick.set(z[kk], -pairs.get(kk)[1]);
//                    }
//                }
//            }
//
//            // ---- Paarweise Constraints ----
//            for (int i = 0; i < n; i++) {
//                for (int j = i + 1; j < n; j++) {
//                    boolean neighbors = G.areNeighbors(G.nodes.get(i), G.nodes.get(j));
//
//                    if (neighbors) {
//                        // Orientierung (L/R/T/B) und Kanten-Gleichheit
//                        Variable L = model.addVariable("L_" + i + "_" + j).binary();
//                        Variable R = model.addVariable("R_" + i + "_" + j).binary();
//                        Variable T = model.addVariable("T_" + i + "_" + j).binary();
//                        Variable B = model.addVariable("B_" + i + "_" + j).binary();
//                        binVars.add(L); binVars.add(R); binVars.add(T); binVars.add(B);
//
//                        // mind. eine Orientierung
//                        model.addExpression("ori_" + i + "_" + j).lower(1).set(L,1).set(R,1).set(T,1).set(B,1);
//
//                        // Gleichheit über je zwei Ungleichungen (aktiv bei Binär=1)
//                        model.addExpression("Lx_le_" + i + "_" + j).upper(Mx)
//                                .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L, Mx);
//                        model.addExpression("Lx_ge_" + i + "_" + j).upper(Mx)
//                                .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(L, Mx);
//
//                        model.addExpression("Rx_le_" + i + "_" + j).upper(Mx)
//                                .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R, Mx);
//                        model.addExpression("Rx_ge_" + i + "_" + j).upper(Mx)
//                                .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(R, Mx);
//
//                        model.addExpression("Ty_le_" + i + "_" + j).upper(My)
//                                .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T, My);
//                        model.addExpression("Ty_ge_" + i + "_" + j).upper(My)
//                                .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(T, My);
//
//                        model.addExpression("By_le_" + i + "_" + j).upper(My)
//                                .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B, My);
//                        model.addExpression("By_ge_" + i + "_" + j).upper(My)
//                                .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(B, My);
//
//                        // Overlap-Variablen: für horizontale Adjazenz Y-Overlap, für vertikale X-Overlap
//                        Variable ovY = model.addVariable("ovY_" + i + "_" + j)
//                                .lower(0).upper(Math.min(h[i].getUpperLimit().intValue(), h[j].getUpperLimit().intValue()));
//                        Variable ovX = model.addVariable("ovX_" + i + "_" + j)
//                                .lower(0).upper(Math.min(w[i].getUpperLimit().intValue(), w[j].getUpperLimit().intValue()));
//
//                        // Bindung von ovY (aktiv, wenn L oder R)
//                        model.addExpression("ovY_leA_" + i + "_" + j).upper(My)
//                                .set(ovY, 1).set(y[i], 1).set(y[j], -1).set(h[j], -1).set(L, My).set(R, My);
//                        model.addExpression("ovY_leB_" + i + "_" + j).upper(My)
//                                .set(ovY, 1).set(y[j], 1).set(y[i], -1).set(h[i], -1).set(L, My).set(R, My);
//                        if (p.minContactCells > 0) {
//                            model.addExpression("ovY_min_" + i + "_" + j).lower(0)
//                                    .set(ovY, 1).set(L, -p.minContactCells).set(R, -p.minContactCells);
//                        }
//
//                        // Bindung von ovX (aktiv, wenn T oder B)
//                        model.addExpression("ovX_leA_" + i + "_" + j).upper(Mx)
//                                .set(ovX, 1).set(x[i], 1).set(x[j], -1).set(w[j], -1).set(T, Mx).set(B, Mx);
//                        model.addExpression("ovX_leB_" + i + "_" + j).upper(Mx)
//                                .set(ovX, 1).set(x[j], 1).set(x[i], -1).set(w[i], -1).set(T, Mx).set(B, Mx);
//                        if (p.minContactCells > 0) {
//                            model.addExpression("ovX_min_" + i + "_" + j).lower(0)
//                                    .set(ovX, 1).set(T, -p.minContactCells).set(B, -p.minContactCells);
//                        }
//
//                        // Für Stufe 2 summieren wir beide Kontaktlängen (nur eine der Orientierungen ist aktiv)
//                        overlaps.add(ovX);
//                        overlaps.add(ovY);
//
//                    } else {
//                        // Nicht-Nachbarn: reine Nicht-Überlappung (ggf. mit 1-Zellen-Abstand)
//                        Variable left  = model.addVariable("left_"  + i + "_" + j).binary();
//                        Variable right = model.addVariable("right_" + i + "_" + j).binary();
//                        Variable above = model.addVariable("above_" + i + "_" + j).binary();
//                        Variable below = model.addVariable("below_" + i + "_" + j).binary();
//                        binVars.add(left); binVars.add(right); binVars.add(above); binVars.add(below);
//
//                        model.addExpression("sep_" + i + "_" + j).lower(1)
//                                .set(left,1).set(right,1).set(above,1).set(below,1);
//
//                        int gapX = p.forbidUnwantedContacts ? 1 : 0;
//                        int gapY = p.forbidUnwantedContacts ? 1 : 0;
//
//                        model.addExpression("sepL_" + i + "_" + j).upper(Mx - gapX)
//                                .set(x[i],1).set(w[i],1).set(x[j],-1).set(left, Mx);
//                        model.addExpression("sepR_" + i + "_" + j).upper(Mx - gapX)
//                                .set(x[j],1).set(w[j],1).set(x[i],-1).set(right, Mx);
//                        model.addExpression("sepA_" + i + "_" + j).upper(My - gapY)
//                                .set(y[i],1).set(h[i],1).set(y[j],-1).set(above, My);
//                        model.addExpression("sepB_" + i + "_" + j).upper(My - gapY)
//                                .set(y[j],1).set(h[j],1).set(y[i],-1).set(below, My);
//                    }
//                }
//            }
//
//            // --- Bounding-Box rechts/oben ---
//            Variable maxX = model.addVariable("maxX").lower(0).upper(p.hullWCells);
//            Variable maxY = model.addVariable("maxY").lower(0).upper(p.hullHCells);
//            for (int i = 0; i < n; i++) {
//                model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
//                model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
//            }
//            if (p.forceFillHull) {
//                model.addExpression("fillX").level(p.hullWCells - p.borderCells).set(maxX, 1);
//                model.addExpression("fillY").level(p.hullHCells - p.borderCells).set(maxY, 1);
//            }
//
//            // --- Stufe 1: kompakteste Lösung ---
//            Expression obj1 = model.addExpression("OBJ1").weight(1.0);
//            obj1.set(maxX, 1.0).set(maxY, 1.0);
//
//            // No-Good-Cuts aus vorherigen Lösungen hinzufügen
//            for (int c = 0; c < noGoods.size(); c++) {
//                boolean[] pat = noGoods.get(c);
//                Expression ng = model.addExpression("nogood_" + c).lower(1);
//                for (int bi = 0; bi < pat.length; bi++) {
//                    if (pat[bi]) ng.set(binVars.get(bi), -1);  // (1 - b)
//                    else         ng.set(binVars.get(bi),  1);  // (b)
//                }
//                // Die Konstante +1 fehlt noch: ojAlgo-Expressions ohne freie Konstante,
//                // trick: addiere alle (1 - b) als (-b) und zähle die Anzahl der 'true' auf die rechte Seite.
//                // Einfacher Workaround: füge eine Hilfsvariable c0 mit fixem Level 1 hinzu.
//            }
//            // (kleiner Hack für die Konstante 1:)
//            Variable cst = model.addVariable("cst").lower(1).upper(1); // konstante 1
//            // Jetzt jede No-Good-Expression um +cst ergänzen:
//            for (int c = 0; c < noGoods.size(); c++) {
//                Expression ng = model.getExpression("nogood_" + c);
//                if (ng != null) ng.set(cst, 1);
//            }
//
//            Optimisation.Result r1 = model.minimise();
//            if (!r1.getState().isFeasible()) break; // keine weitere Lösung
//
//            double opt1 = maxX.getValue().doubleValue() + maxY.getValue().doubleValue();
//
//            // --- Stufe 2: bei gleicher Kompaktheit Kontaktlänge maximieren ---
//            // Fixiere Stufe-1-Optimum:
//            model.getExpression("OBJ1").weight(0.0);
//            Expression fix1 = model.addExpression("FIX1").level(opt1);
//            fix1.set(maxX, 1.0).set(maxY, 1.0);
//
//            // neues Ziel: - Summe Overlaps (minimieren ⇒ maximiert Summe Overlaps)
//            Expression obj2 = model.addExpression("OBJ2").weight(1.0);
//            for (Variable ov : overlaps) obj2.set(ov, -1.0);
//
//            Optimisation.Result r2 = model.minimise();
//            if (!r2.getState().isFeasible()) {
//                // sollte nicht passieren – dann nimm r1
//            }
//
//            // --- Lösung extrahieren ---
//            FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
//            sol.maxX = maxX.getValue().intValue();
//            sol.maxY = maxY.getValue().intValue();
//            for (int i = 0; i < n; i++) {
//                FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
//                rp.x0 = x[i].getValue().intValue();
//                rp.y0 = y[i].getValue().intValue();
//                rp.wH = rp.tX = w[i].getValue().intValue();
//                rp.tY = rp.hV = h[i].getValue().intValue();
//                rp.wBB = rp.wH;
//                rp.hBB = rp.tY;
//                sol.piecesByNode.put(G.nodes.get(i), rp);
//            }
//            out.add(sol);
//
//            // --- No-Good-Muster aufzeichnen (über alle Binärvariablen) ---
//            boolean[] pat = new boolean[binVars.size()];
//            for (int bi = 0; bi < binVars.size(); bi++) {
//                double v = binVars.get(bi).getValue().doubleValue();
//                pat[bi] = v >= 0.5;
//            }
//            noGoods.add(pat);
//        }
//
//        return out.isEmpty() ? List.of() : out;
//    }
//
//    static List<int[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim,
//                                           int tolCells, double minAsp, double maxAsp) {
//        if (areaSqm < 0 || cellM <= 0) throw new IllegalArgumentException("Ungültige Flächenparameter");
//        if (tolCells < 0) throw new IllegalArgumentException("tolCells < 0");
//        if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim");
//        if (minAsp > maxAsp) throw new IllegalArgumentException("minAsp > maxAsp");
//        int target = (int) Math.round(areaSqm / (cellM * cellM));
//        List<int[]> pairs = new ArrayList<>();
//        for (int W = minDim; W <= maxDim; W++) {
//            for (int H = minDim; H <= maxDim; H++) {
//                int cells = W * H;
//                if (Math.abs(cells - target) <= tolCells) {
//                    double ar = (double) W / H;
//                    if (ar >= minAsp && ar <= maxAsp) pairs.add(new int[]{W, H});
//                }
//            }
//        }
//        return pairs;
//    }
//}

import java.util.*;
import org.ojalgo.optimisation.ExpressionsBasedModel;
import org.ojalgo.optimisation.Optimisation;
import org.ojalgo.optimisation.Expression;
import org.ojalgo.optimisation.Variable;

/**
 * MILP-Variante des Grundriss-Solvers (rechteckige Räume). - Nicht-Überlappung
 * für alle Paare (auch Nachbarn) - Nachbarn: bündige Adjazenz mit
 * Mindest-Kontaktlänge (optional) - Straffe Big-Ms, Border-Rand, optionale
 * Flächendiskretisierung (enforceAreas)
 */
public final class MIQPFloorplanSolver {

	private MIQPFloorplanSolver() {
	}

	public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {

		// --- Eingaben prüfen ---
		var G = p.G;
		int n = G.nodes.size();
		if (n == 0)
			throw new IllegalArgumentException("Graph hat keine Knoten");
		if (p.minDim > p.maxDim)
			throw new IllegalArgumentException("minDim > maxDim");
		if (p.hullWCells <= 0 || p.hullHCells <= 0)
			throw new IllegalArgumentException("Ungültige Bounding-Box");
		if (p.borderCells < 0)
			throw new IllegalArgumentException("borderCells < 0");
		if (p.areaTolCells < 0)
			throw new IllegalArgumentException("areaTolCells < 0");
		if (p.minAspect > p.maxAspect)
			throw new IllegalArgumentException("minAspect > maxAspect");

		// --- Innenhülle & straffe Big-Ms ---
		final int hullWUse = Math.max(0, p.hullWCells - 2 * p.borderCells);
		final int hullHUse = Math.max(0, p.hullHCells - 2 * p.borderCells);
		if (hullWUse <= 0 || hullHUse <= 0) {
			throw new IllegalArgumentException("Nutzbare Hülle ≤ 0 (borderCells zu groß?)");
		}
		final int Mx = hullWUse;
		final int My = hullHUse;

		// --- Modell & Variablen ---
		ExpressionsBasedModel model = new ExpressionsBasedModel();

		Variable[] x = new Variable[n];
		Variable[] y = new Variable[n];
		Variable[] w = new Variable[n];
		Variable[] h = new Variable[n];

		for (int i = 0; i < n; i++) {
			String id = G.nodes.get(i);
			var ro = p.perRoomOptions.get(id);

			int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
			int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
			if (minDim > maxDim)
				throw new IllegalArgumentException("minDim > maxDim für Raum " + id);

			double minAspLocal = (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect;
			double maxAspLocal = (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect;
			if (minAspLocal > maxAspLocal)
				throw new IllegalArgumentException("minAspect > maxAspect für Raum " + id);

			// Koordinaten am Innenrand
			x[i] = model.addVariable("x_" + i).lower(p.borderCells).upper(p.hullWCells - p.borderCells).integer(true);
			y[i] = model.addVariable("y_" + i).lower(p.borderCells).upper(p.hullHCells - p.borderCells).integer(true);

			// w/h: bei enforceAreas mit Ziel-Fläche können kontinuierlich sein (z
			// diskretisiert)
			Integer areaSqm = p.areaSqmByRoom.get(id);
			boolean makeWHInteger = !(p.enforceAreas && areaSqm != null);
			w[i] = model.addVariable("w_" + i).lower(minDim).upper(maxDim).integer(makeWHInteger);
			h[i] = model.addVariable("h_" + i).lower(minDim).upper(maxDim).integer(makeWHInteger);

			// Innenhülle
			model.addExpression("bndX_" + i).upper(p.hullWCells - p.borderCells).set(x[i], 1).set(w[i], 1);
			model.addExpression("bndY_" + i).upper(p.hullHCells - p.borderCells).set(y[i], 1).set(h[i], 1);

			// Diskrete (w,h)-Paare bei enforceAreas
			if (areaSqm != null && p.enforceAreas) {
				var pairs = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim, p.areaTolCells, minAspLocal,
						maxAspLocal);
				if (pairs.isEmpty())
					throw new IllegalStateException("Keine (w,h)-Paare für " + id);

				Variable[] z = new Variable[pairs.size()];
				for (int kk = 0; kk < pairs.size(); kk++)
					z[kk] = model.addVariable("z_" + i + "_" + kk).binary();

				// genau ein Paar
				Expression sel = model.addExpression("sel_" + i).level(1);
				for (var zk : z)
					sel.set(zk, 1);

				// w = Sum z_k * w_k; h = Sum z_k * h_k
				Expression wPick = model.addExpression("wPick_" + i).level(0).set(w[i], 1);
				Expression hPick = model.addExpression("hPick_" + i).level(0).set(h[i], 1);
				for (int kk = 0; kk < pairs.size(); kk++) {
					wPick.set(z[kk], -pairs.get(kk)[0]);
					hPick.set(z[kk], -pairs.get(kk)[1]);
				}
			}
		}

		// --- Paarweise Constraints ---
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				String idI = G.nodes.get(i);
				String idJ = G.nodes.get(j);
				boolean neighbors = G.areNeighbors(idI, idJ);

				// Immer: Nicht-Überlappung (4 Binärs)
				Variable left = model.addVariable("left_" + i + "_" + j).binary();
				Variable right = model.addVariable("right_" + i + "_" + j).binary();
				Variable above = model.addVariable("above_" + i + "_" + j).binary();
				Variable below = model.addVariable("below_" + i + "_" + j).binary();

				// Mindestens eine Richtung
				model.addExpression("sep_" + i + "_" + j).lower(1).set(left, 1).set(right, 1).set(above, 1).set(below,
						1);

				// Abstand nur für Nicht-Nachbarn
				int gapX = (!neighbors && p.forbidUnwantedContacts) ? 1 : 0;
				int gapY = (!neighbors && p.forbidUnwantedContacts) ? 1 : 0;

				// Trennungen (auch für Nachbarn; dort gap=0 erlaubt bündiges Anlegen)
				model.addExpression("sepL_" + i + "_" + j).upper(Mx - gapX).set(x[i], 1).set(w[i], 1).set(x[j], -1)
						.set(left, Mx);
				model.addExpression("sepR_" + i + "_" + j).upper(Mx - gapX).set(x[j], 1).set(w[j], 1).set(x[i], -1)
						.set(right, Mx);
				model.addExpression("sepA_" + i + "_" + j).upper(My - gapY).set(y[i], 1).set(h[i], 1).set(y[j], -1)
						.set(above, My);
				model.addExpression("sepB_" + i + "_" + j).upper(My - gapY).set(y[j], 1).set(h[j], 1).set(y[i], -1)
						.set(below, My);

				if (neighbors) {
					// Orientierung
					Variable L = model.addVariable("L_" + i + "_" + j).binary();
					Variable R = model.addVariable("R_" + i + "_" + j).binary();
					Variable T = model.addVariable("T_" + i + "_" + j).binary();
					Variable B = model.addVariable("B_" + i + "_" + j).binary();

					if (p.adjAtLeastOneSide)
						model.addExpression("ori_" + i + "_" + j).lower(1).set(L, 1).set(R, 1).set(T, 1).set(B, 1);
					else
						model.addExpression("ori_" + i + "_" + j).level(1).set(L, 1).set(R, 1).set(T, 1).set(B, 1);

					// Harte Kopplung: Trenn-Binärs == Orientierungs-Binärs
					model.addExpression("eq_left_" + i + "_" + j).level(0).set(left, 1).set(L, -1);
					model.addExpression("eq_right_" + i + "_" + j).level(0).set(right, 1).set(R, -1);
					model.addExpression("eq_above_" + i + "_" + j).level(0).set(above, 1).set(T, -1);
					model.addExpression("eq_below_" + i + "_" + j).level(0).set(below, 1).set(B, -1);

					// Gleichheiten (aktiv, wenn jeweilige Orientierung=1)
					model.addExpression("Lx_le_" + i + "_" + j).upper(Mx).set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L,
							Mx);
					model.addExpression("Lx_ge_" + i + "_" + j).upper(Mx).set(x[j], 1).set(x[i], -1).set(w[i], -1)
							.set(L, Mx);

					model.addExpression("Rx_le_" + i + "_" + j).upper(Mx).set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R,
							Mx);
					model.addExpression("Rx_ge_" + i + "_" + j).upper(Mx).set(x[i], 1).set(x[j], -1).set(w[j], -1)
							.set(R, Mx);

					model.addExpression("Ty_le_" + i + "_" + j).upper(My).set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T,
							My);
					model.addExpression("Ty_ge_" + i + "_" + j).upper(My).set(y[j], 1).set(y[i], -1).set(h[i], -1)
							.set(T, My);

					model.addExpression("By_le_" + i + "_" + j).upper(My).set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B,
							My);
					model.addExpression("By_ge_" + i + "_" + j).upper(My).set(y[i], 1).set(y[j], -1).set(h[j], -1)
							.set(B, My);

					// Mindest-Kontaktlänge (optional)
					if (p.minContactCells > 0) {
						Variable overlapY = model.addVariable("overlapY_" + i + "_" + j).lower(0)
								.upper(Math.min(h[i].getUpperLimit().intValue(), h[j].getUpperLimit().intValue()));
						Variable overlapX = model.addVariable("overlapX_" + i + "_" + j).lower(0)
								.upper(Math.min(w[i].getUpperLimit().intValue(), w[j].getUpperLimit().intValue()));

						// Y-Overlap aktiv bei L oder R
						model.addExpression("ovY_leA_" + i + "_" + j).upper(My).set(overlapY, 1).set(y[i], 1)
								.set(y[j], -1).set(h[j], -1).set(L, My).set(R, My);
						model.addExpression("ovY_leB_" + i + "_" + j).upper(My).set(overlapY, 1).set(y[j], 1)
								.set(y[i], -1).set(h[i], -1).set(L, My).set(R, My);
						model.addExpression("ovY_min_" + i + "_" + j).lower(0).set(overlapY, 1)
								.set(L, -p.minContactCells).set(R, -p.minContactCells);

						// X-Overlap aktiv bei T oder B
						model.addExpression("ovX_leA_" + i + "_" + j).upper(Mx).set(overlapX, 1).set(x[i], 1)
								.set(x[j], -1).set(w[j], -1).set(T, Mx).set(B, Mx);
						model.addExpression("ovX_leB_" + i + "_" + j).upper(Mx).set(overlapX, 1).set(x[j], 1)
								.set(x[i], -1).set(w[i], -1).set(T, Mx).set(B, Mx);
						model.addExpression("ovX_min_" + i + "_" + j).lower(0).set(overlapX, 1)
								.set(T, -p.minContactCells).set(B, -p.minContactCells);
					}
				} else { // ---- NICHT-NACHBARN: keine gemeinsame Kante, Ecken erlaubt ----
							// Trenn-Binärs wie gehabt
					 left = model.addVariable("left_" + i + "_" + j).binary();
					 right = model.addVariable("right_" + i + "_" + j).binary();
					 above = model.addVariable("above_" + i + "_" + j).binary();
					 below = model.addVariable("below_" + i + "_" + j).binary();

					// Mindestens eine Richtung aktiv
					model.addExpression("sep_" + i + "_" + j).lower(1).set(left, 1).set(right, 1).set(above, 1)
							.set(below, 1);

					// Aggregatoren: H=1, wenn (left oder right); V=1, wenn (above oder below)
					Variable H = model.addVariable("H_" + i + "_" + j).binary();
					Variable V = model.addVariable("V_" + i + "_" + j).binary();
					model.addExpression("H_ge_left_" + i + "_" + j).lower(0).set(H, 1).set(left, -1);
					model.addExpression("H_ge_right_" + i + "_" + j).lower(0).set(H, 1).set(right, -1);
					model.addExpression("H_le_sum_" + i + "_" + j).upper(0).set(H, 1).set(left, -1).set(right, -1);

					model.addExpression("V_ge_above_" + i + "_" + j).lower(0).set(V, 1).set(above, -1);
					model.addExpression("V_ge_below_" + i + "_" + j).lower(0).set(V, 1).set(below, -1);
					model.addExpression("V_le_sum_" + i + "_" + j).upper(0).set(V, 1).set(above, -1).set(below, -1);

					// Trennungen (klassisches Big-M), aber:
					// Wenn *keine* vertikale Trennung (V=0) und horizontal getrennt (left=1 oder
					// right=1),
					// dann erzwinge *mindestens 1 Zelle* horizontalen Abstand (=> keine gemeinsame
					// vertikale Kante).
					model.addExpression("sepL_" + i + "_" + j).upper(Mx - 1).set(x[i], 1).set(w[i], 1).set(x[j], -1)
							.set(left, Mx) // Big-M
							.set(V, -1); // -V verschiebt RHS: bei V=0 -> ≤ -1, bei V=1 -> ≤ 0

					model.addExpression("sepR_" + i + "_" + j).upper(Mx - 1).set(x[j], 1).set(w[j], 1).set(x[i], -1)
							.set(right, Mx).set(V, -1);

					// Symmetrisch: Wenn *keine* horizontale Trennung (H=0) und vertikal getrennt,
					// dann erzwinge *mindestens 1 Zelle* vertikalen Abstand (=> keine gemeinsame
					// horizontale Kante).
					model.addExpression("sepA_" + i + "_" + j).upper(My - 1).set(y[i], 1).set(h[i], 1).set(y[j], -1)
							.set(above, My).set(H, -1);

					model.addExpression("sepB_" + i + "_" + j).upper(My - 1).set(y[j], 1).set(h[j], 1).set(y[i], -1)
							.set(below, My).set(H, -1);
				}

			}
		}

		// --- Bounding-Box & Ziel (linear) ---
		Variable maxX = model.addVariable("maxX").lower(0).upper(p.hullWCells);
		Variable maxY = model.addVariable("maxY").lower(0).upper(p.hullHCells);
		for (int i = 0; i < n; i++) {
			model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
			model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
		}
		if (p.forceFillHull) {
			model.addExpression("fillX").level(p.hullWCells - p.borderCells).set(maxX, 1);
			model.addExpression("fillY").level(p.hullHCells - p.borderCells).set(maxY, 1);
		}

		Expression obj = model.addExpression("OBJ").weight(1.0);
		obj.set(maxX, 1.0).set(maxY, 1.0); // MILP: kompakte Hülle

		// --- Lösen ---
		Optimisation.Result res = model.minimise();
		if (!res.getState().isFeasible()) {
			System.out.println("Keine feasible Lösung gefunden – Status: " + res.getState());
			return List.of();
		}

		// --- Lösung extrahieren (runden statt truncaten) ---
		FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
		sol.maxX = (int) Math.round(maxX.getValue().doubleValue());
		sol.maxY = (int) Math.round(maxY.getValue().doubleValue());
		for (int i = 0; i < n; i++) {
			FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
			rp.x0 = (int) Math.round(x[i].getValue().doubleValue());
			rp.y0 = (int) Math.round(y[i].getValue().doubleValue());
			rp.wH = rp.tX = (int) Math.round(w[i].getValue().doubleValue());
			rp.tY = rp.hV = (int) Math.round(h[i].getValue().doubleValue());
			rp.wBB = rp.wH;
			rp.hBB = rp.tY;
			sol.piecesByNode.put(G.nodes.get(i), rp);
		}

		return List.of(sol);
	}

	// Erlaubte (w,h)-Paare zur Flächenerzwingung
	static List<int[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim, int tolCells,
			double minAsp, double maxAsp) {
		if (areaSqm < 0 || cellM <= 0)
			throw new IllegalArgumentException("Ungültige Flächenparameter");
		if (tolCells < 0)
			throw new IllegalArgumentException("tolCells < 0");
		if (minDim > maxDim)
			throw new IllegalArgumentException("minDim > maxDim");
		if (minAsp > maxAsp)
			throw new IllegalArgumentException("minAsp > maxAsp");
		int target = (int) Math.round(areaSqm / (cellM * cellM));
		List<int[]> pairs = new ArrayList<>();
		for (int W = minDim; W <= maxDim; W++) {
			for (int H = minDim; H <= maxDim; H++) {
				int cells = W * H;
				if (Math.abs(cells - target) <= tolCells) {
					double ar = (double) W / H;
					if (ar >= minAsp && ar <= maxAsp) {
						pairs.add(new int[] { W, H });
					}
				}
			}
		}
		return pairs;
	}
}
