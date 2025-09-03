////package pdf;
////
////import com.google.ortools.Loader;
////import com.google.ortools.sat.*;
////import java.util.*;
////
////public final class FloorplanCrossModalSolver {
////
////	// ---------- Hilfsfunktionen für LinearExpr ----------
////	private static LinearExpr S(LinearArgument... terms) {
////		return LinearExpr.newBuilder().addSum(terms).build();
////	}
////
////	private static LinearExpr Splus(long constant, LinearArgument... terms) {
////		LinearExprBuilder b = LinearExpr.newBuilder();
////		b.add(constant);
////		for (LinearArgument t : terms)
////			b.add(t);
////		return b.build();
////	}
////
////	// ---------- Variablen pro Raum (L-Form als 2 Balken, gemeinsamer Ursprung)
////	// --------
////	static final class VarsL {
////		IntVar x0, y0; // gemeinsame Ecke
////		// Horizontaler Balken (wH x tY)
////		IntVar wH, tY, rightH, topH;
////		// Vertikaler Balken (tX x hV)
////		IntVar tX, hV, rightV, topV;
////		// Bounding Box
////		IntVar rightBB, topBB;
////		// Fläche
////		IntVar p1, p2, p3, areaCells;
////	}
////
////	// ---------- Lösung ----------
////	public static final class RoomPieces {
////		public int x0, y0, wH, tY, tX, hV, wBB, hBB;
////	}
////
////	public static final class Solution {
////		public final Map<String, RoomPieces> piecesByNode = new LinkedHashMap<>();
////		public int maxX, maxY;
////
////		@Override
////		public String toString() {
////			StringBuilder sb = new StringBuilder();
////			sb.append("maxX=").append(maxX).append(", maxY=").append(maxY).append('\n');
////			for (var e : piecesByNode.entrySet()) {
////				var p = e.getValue();
////				int areaCells = p.wH * p.tY + p.tX * p.hV - p.tX * p.tY;
////				double areaSqm = areaCells * 1.0; // echte m² skaliert der Viewer
////				sb.append(String.format("%s -> x0=%d y0=%d H:(%d×%d) V:(%d×%d) BB:(%d×%d) cells=%d\n", e.getKey(), p.x0,
////						p.y0, p.wH, p.tY, p.tX, p.hV, p.wBB, p.hBB, areaCells));
////			}
////			return sb.toString();
////		}
////	}
////
////	// ---------- Callback zum Einsammeln ----------
//////	static final class MultipleSolutionsCallback extends CpSolverSolutionCallback {
//////		private final List<Solution> solutions = new ArrayList<>();
//////		private final VarsL[] V;
//////		private final CrossModalMapper.AdjacencyGraph G;
//////		private final IntVar maxX, maxY;
//////
//////		MultipleSolutionsCallback(VarsL[] V, CrossModalMapper.AdjacencyGraph G, IntVar maxX, IntVar maxY) {
//////			this.V = V;
//////			this.G = G;
//////			this.maxX = maxX;
//////			this.maxY = maxY;
//////		}
//////
//////		public List<Solution> getSolutions() {
//////			return solutions;
//////		}
//////
//////		@Override
//////		public void onSolutionCallback() {
//////			Solution sol = new Solution();
//////			sol.maxX = (int) value(maxX);
//////			sol.maxY = (int) value(maxY);
//////			for (int i = 0; i < V.length; i++) {
//////				var vi = V[i];
//////				var rp = new RoomPieces();
//////				rp.x0 = (int) value(vi.x0);
//////				rp.y0 = (int) value(vi.y0);
//////				rp.wH = (int) value(vi.wH);
//////				rp.tY = (int) value(vi.tY);
//////				rp.tX = (int) value(vi.tX);
//////				rp.hV = (int) value(vi.hV);
//////				int rightBB = (int) value(vi.rightBB);
//////				int topBB = (int) value(vi.topBB);
//////				rp.wBB = rightBB - rp.x0;
//////				rp.hBB = topBB - rp.y0;
//////				sol.piecesByNode.put(G.nodes.get(i), rp);
//////			}
//////			solutions.add(sol);
//////		}
//////	}
////
////	// ---------- Utilities ----------
////	private static BoolVar ltBool(CpModel m, IntVar a, IntVar b, String name) {
////		BoolVar z = m.newBoolVar(name);
////		m.addLessOrEqual(Splus(1, a), b).onlyEnforceIf(z); // a+1 <= b ⇒ a < b
////		m.addGreaterOrEqual(a, b).onlyEnforceIf(z.not()); // a >= b ⇒ ¬(a<b)
////		return z;
////	}
////
////	private static void noOverlap(CpModel m, IntVar ax, IntVar ay, IntVar aw, IntVar ah, IntVar bx, IntVar by,
////			IntVar bw, IntVar bh, int sep, String tag) {
////		BoolVar left = m.newBoolVar("left_" + tag);
////		BoolVar right = m.newBoolVar("right_" + tag);
////		BoolVar above = m.newBoolVar("above_" + tag);
////		BoolVar below = m.newBoolVar("below_" + tag);
////		m.addEquality(LinearExpr.sum(new LinearArgument[] { left, right, above, below }), LinearExpr.constant(1));
////		// ax+aw <= bx - sep (A links von B)
////		m.addLessOrEqual(S(ax, aw), Splus(-sep, bx)).onlyEnforceIf(left);
////		// bx+bw <= ax - sep
////		m.addLessOrEqual(S(bx, bw), Splus(-sep, ax)).onlyEnforceIf(right);
////		// ay+ah <= by - sep
////		m.addLessOrEqual(S(ay, ah), Splus(-sep, by)).onlyEnforceIf(above);
////		// by+bh <= ay - sep
////		m.addLessOrEqual(S(by, bh), Splus(-sep, ay)).onlyEnforceIf(below);
////	}
////
////	private static void contactAtLeastOneSide(CpModel m, IntVar ax, IntVar ay, IntVar aw, IntVar ah, IntVar aright,
////			IntVar atop, IntVar bx, IntVar by, IntVar bw, IntVar bh, IntVar bright, IntVar btop, int minContactCells,
////			List<BoolVar> collectors, String tag) {
////		// Vier Richtungen A↔B, je ein Bool; wir fügen die Bools dem collectors hinzu
////		BoolVar L = m.newBoolVar("L_" + tag); // A rechts an B links
////		BoolVar R = m.newBoolVar("R_" + tag); // A links an B rechts
////		BoolVar T = m.newBoolVar("T_" + tag); // A unten an B oben
////		BoolVar B = m.newBoolVar("B_" + tag); // A oben an B unten
////
////		// L: aright == bx & y-Intervalle überlappen mind. minContact
////		m.addEquality(aright, bx).onlyEnforceIf(L);
////		m.addLessOrEqual(ay, Splus(-minContactCells, by, bh)).onlyEnforceIf(L);
////		m.addLessOrEqual(by, Splus(-minContactCells, ay, ah)).onlyEnforceIf(L);
////
////		// R: bright == ax
////		m.addEquality(bright, ax).onlyEnforceIf(R);
////		m.addLessOrEqual(ay, Splus(-minContactCells, by, bh)).onlyEnforceIf(R);
////		m.addLessOrEqual(by, Splus(-minContactCells, ay, ah)).onlyEnforceIf(R);
////
////		// T: atop == by & x-Überlappung >= minContact
////		m.addEquality(atop, by).onlyEnforceIf(T);
////		m.addLessOrEqual(ax, Splus(-minContactCells, bx, bw)).onlyEnforceIf(T);
////		m.addLessOrEqual(bx, Splus(-minContactCells, ax, aw)).onlyEnforceIf(T);
////
////		// B: btop == ay
////		m.addEquality(btop, ay).onlyEnforceIf(B);
////		m.addLessOrEqual(ax, Splus(-minContactCells, bx, bw)).onlyEnforceIf(B);
////		m.addLessOrEqual(bx, Splus(-minContactCells, ax, aw)).onlyEnforceIf(B);
////
////		collectors.add(L);
////		collectors.add(R);
////		collectors.add(T);
////		collectors.add(B);
////	}
////
////	// ---------- Hauptmethode --------------------------------------------------
////	public static List<Solution> layout(CrossModalMapper.MappedParams p) {
////		Loader.loadNativeLibraries();
////		CpModel model = new CpModel();
////
////		var G = p.G;
////		if (G.nodes.isEmpty())
////			throw new IllegalArgumentException("Graph hat keine Knoten.");
////		int n = G.nodes.size();
////
////		VarsL[] V = new VarsL[n];
////
////		// Variablen je Raum
////		for (int i = 0; i < n; i++) {
////			V[i] = new VarsL();
////			var ro = p.perRoomOptions.get(G.nodes.get(i));
////			int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
////			int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
////			int minT = ro != null && ro.minThicknessCells != null ? ro.minThicknessCells : 1;
////			int maxT = ro != null && ro.maxThicknessCells != null ? ro.maxThicknessCells : Math.min(4, maxDim);
////
////			V[i].x0 = model.newIntVar(0, p.hullWCells, "x0_" + i);
////			V[i].y0 = model.newIntVar(0, p.hullHCells, "y0_" + i);
////
////			V[i].wH = model.newIntVar(minDim, maxDim, "wH_" + i);
////			V[i].tY = model.newIntVar(minT, maxT, "tY_" + i);
////			V[i].tX = model.newIntVar(minT, maxT, "tX_" + i);
////			V[i].hV = model.newIntVar(minDim, maxDim, "hV_" + i);
////
////			V[i].rightH = model.newIntVar(0, p.hullWCells, "rightH_" + i);
////			V[i].topH = model.newIntVar(0, p.hullHCells, "topH_" + i);
////			V[i].rightV = model.newIntVar(0, p.hullWCells, "rightV_" + i);
////			V[i].topV = model.newIntVar(0, p.hullHCells, "topV_" + i);
////			model.addEquality(V[i].rightH, S(V[i].x0, V[i].wH));
////			model.addEquality(V[i].topH, S(V[i].y0, V[i].tY));
////			model.addEquality(V[i].rightV, S(V[i].x0, V[i].tX));
////			model.addEquality(V[i].topV, S(V[i].y0, V[i].hV));
////
////			V[i].rightBB = model.newIntVar(0, p.hullWCells, "rightBB_" + i);
////			V[i].topBB = model.newIntVar(0, p.hullHCells, "topBB_" + i);
////			model.addMaxEquality(V[i].rightBB, new IntVar[] { V[i].rightH, V[i].rightV });
////			model.addMaxEquality(V[i].topBB, new IntVar[] { V[i].topH, V[i].topV });
////
////			// Hülle einhalten
////			model.addGreaterOrEqual(V[i].x0, p.borderCells);
////			model.addGreaterOrEqual(V[i].y0, p.borderCells);
////
////			model.addLessOrEqual(V[i].rightBB, LinearExpr.constant(p.hullWCells - p.borderCells));
////			model.addLessOrEqual(V[i].topBB, LinearExpr.constant(p.hullHCells - p.borderCells));
////
////			// Fläche: areaCells = wH*tY + tX*hV - tX*tY
////			V[i].p1 = model.newIntVar(0, p.hullWCells * p.hullHCells, "p1_" + i);
////			V[i].p2 = model.newIntVar(0, p.hullWCells * p.hullHCells, "p2_" + i);
////			V[i].p3 = model.newIntVar(0, p.hullWCells * p.hullHCells, "p3_" + i);
////			V[i].areaCells = model.newIntVar(0, p.hullWCells * p.hullHCells, "area_" + i);
////			model.addMultiplicationEquality(V[i].p1, V[i].wH, V[i].tY);
////			model.addMultiplicationEquality(V[i].p2, V[i].tX, V[i].hV);
////			model.addMultiplicationEquality(V[i].p3, V[i].tX, V[i].tY);
////			LinearExprBuilder ab = LinearExpr.newBuilder();
////			ab.add(V[i].p1);
////			ab.add(V[i].p2);
////			ab.addTerm(V[i].p3, -1);
////			model.addEquality(V[i].areaCells, ab.build());
////
////			// Shape-Restriktionen
////			CrossModalSpec.Shape shape = ro != null && ro.shape != null ? ro.shape : CrossModalSpec.Shape.RECT_OR_L;
////			if (shape == CrossModalSpec.Shape.RECT) {
////				model.addEquality(V[i].tX, V[i].wH);
////				model.addEquality(V[i].tY, V[i].hV);
////			} else if (shape == CrossModalSpec.Shape.L) {
////				BoolVar bx = ltBool(model, V[i].tX, V[i].wH, "tX_lt_wH_" + i);
////				BoolVar by = ltBool(model, V[i].tY, V[i].hV, "tY_lt_hV_" + i);
////				model.addEquality(LinearExpr.sum(new LinearArgument[] { bx, by }), LinearExpr.constant(2)); // beide
////																											// echt
////																											// kleiner
//////				model.addGreaterOrEqual(V[i].wH, Splus(1, V[i].tX));
//////				model.addGreaterOrEqual(V[i].hV, Splus(1, V[i].tY));
////
////			} // RECT_OR_L: keine Zusatzbindung
////		}
////		// Flächenvorgaben
////
////		// maxX/maxY (BB) optional auf Hülle festklemmen
////		IntVar maxX = model.newIntVar(0, p.hullWCells, "maxX");
////		IntVar maxY = model.newIntVar(0, p.hullHCells, "maxY");
////		model.addMaxEquality(maxX, Arrays.stream(V).map(v -> v.rightBB).toArray(IntVar[]::new));
////		model.addMaxEquality(maxY, Arrays.stream(V).map(v -> v.topBB).toArray(IntVar[]::new));
////		if (p.forceFillHull) {
////			model.addEquality(maxX, p.hullWCells);
////			model.addEquality(maxY, p.hullHCells);
////		}
////
////		// Anker (optional): ersten Raum oben links fixieren
////		model.addEquality(V[0].x0, p.borderCells);
////		model.addEquality(V[0].y0, p.borderCells);
////
////		// Paare
////		List<int[]> pairs = new ArrayList<>();
////		for (int i = 0; i < n; i++)
////			for (int j = i + 1; j < n; j++)
////				pairs.add(new int[] { i, j });
////
////		// Adjazenzen
////		for (int[] pr : pairs) {
////			int i = pr[0], j = pr[1];
////			String idI = G.nodes.get(i), idJ = G.nodes.get(j);
////			if (G.areNeighbors(idI, idJ)) {
////				List<BoolVar> contactBools = new ArrayList<>();
////				// (i.H) vs (j.H)
////				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[i].rightH, V[i].topH, V[j].x0,
////						V[j].y0, V[j].wH, V[j].tY, V[j].rightH, V[j].topH, p.minContactCells, contactBools,
////						"iH_jH_" + i + "_" + j);
////				// (i.H) vs (j.V)
////				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[i].rightH, V[i].topH, V[j].x0,
////						V[j].y0, V[j].tX, V[j].hV, V[j].rightV, V[j].topV, p.minContactCells, contactBools,
////						"iH_jV_" + i + "_" + j);
////				// (i.V) vs (j.H)
////				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[i].rightV, V[i].topV, V[j].x0,
////						V[j].y0, V[j].wH, V[j].tY, V[j].rightH, V[j].topH, p.minContactCells, contactBools,
////						"iV_jH_" + i + "_" + j);
////				// (i.V) vs (j.V)
////				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[i].rightV, V[i].topV, V[j].x0,
////						V[j].y0, V[j].tX, V[j].hV, V[j].rightV, V[j].topV, p.minContactCells, contactBools,
////						"iV_jV_" + i + "_" + j);
////
////				// Insgesamt mindestens eine Seitenberührung zwischen i und j
////				if (p.adjAtLeastOneSide) {
////					mAddGreaterOrEqualSum(model, contactBools, 1);
////				} else {
////					// selten sinnvoll mit L-Formen, aber möglich: genau 1 Seite
////					model.addEquality(LinearExpr.sum(contactBools.toArray(new LinearArgument[0])),
////							LinearExpr.constant(1));
////				}
////			}
////		}
////
////		// Nicht-Nachbarn trennen (alle Stück-Kombinationen)
////		int sep = p.forbidUnwantedContacts ? 1 : 0;
////		for (int[] pr : pairs) {
////			int i = pr[0], j = pr[1];
////			String idI = G.nodes.get(i), idJ = G.nodes.get(j);
////			if (!G.areNeighbors(idI, idJ)) {
////				noOverlap(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[j].x0, V[j].y0, V[j].wH, V[j].tY, sep,
////						"iH_jH_" + i + "_" + j);
////				noOverlap(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[j].x0, V[j].y0, V[j].tX, V[j].hV, sep,
////						"iH_jV_" + i + "_" + j);
////				noOverlap(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[j].x0, V[j].y0, V[j].wH, V[j].tY, sep,
////						"iV_jH_" + i + "_" + j);
////				noOverlap(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[j].x0, V[j].y0, V[j].tX, V[j].hV, sep,
////						"iV_jV_" + i + "_" + j);
////			}
////		}
////
////		// Solver
////		CpSolver solver = new CpSolver();
////		solver.getParameters().setEnumerateAllSolutions(true);
//////		solver.getParameters().setEnumerateAllSolutions(false); // nur 1 Lösung
////		solver.getParameters().setNumSearchWorkers(1);
////		solver.getParameters().setMaxTimeInSeconds(30.0);
//////		solver.getParameters().setSolutionPoolSize(10); // stoppe nach 20 Lösungen
////
//////		MultipleSolutionsCallback cb = new MultipleSolutionsCallback(V, G, maxX, maxY);
////		// Callback mit Limits
////		int MAX_SOLUTIONS = 20;
////		double MAX_WALL_TIME = 10.0; // Sekunden, optional
////		MultipleSolutionsCallback cb =
////		    new MultipleSolutionsCallback(V, G, maxX, maxY, MAX_SOLUTIONS, MAX_WALL_TIME);
////
////
////		CpSolverStatus status = solver.solve(model, cb);
////		if (status != CpSolverStatus.OPTIMAL && status != CpSolverStatus.FEASIBLE)
////			throw new RuntimeException("Kein Layout gefunden: " + status);
////		return cb.getSolutions();
////	}
////
////	private static void mAddGreaterOrEqualSum(CpModel m, List<BoolVar> bs, int k) {
////		LinearArgument[] arr = bs.toArray(new LinearArgument[0]);
////		m.addGreaterOrEqual(LinearExpr.sum(arr), LinearExpr.constant(k));
////	}
////	// in FloorplanCrossModalSolver
////
////	static final class MultipleSolutionsCallback extends CpSolverSolutionCallback {
////		private final List<Solution> solutions = new ArrayList<>();
////		private final VarsL[] V;
////		private final CrossModalMapper.AdjacencyGraph G;
////		private final IntVar maxX, maxY;
////
////		// NEU:
////		private final int maxSolutions; // z. B. 20
////		private final double maxWallTimeSeconds; // optional: Sicherheitsbremse (z. B. 10.0)
////
////		MultipleSolutionsCallback(VarsL[] V, CrossModalMapper.AdjacencyGraph G, IntVar maxX, IntVar maxY,
////				int maxSolutions, double maxWallTimeSeconds) {
////			this.V = V;
////			this.G = G;
////			this.maxX = maxX;
////			this.maxY = maxY;
////			this.maxSolutions = maxSolutions;
////			this.maxWallTimeSeconds = maxWallTimeSeconds;
////			System.out.println("a");
////		}
////
////		public List<Solution> getSolutions() {
////			return solutions;
////		}
////
////		@Override
////		public void onSolutionCallback() {
////			// Lösung einsammeln
////			Solution sol = new Solution();
////			sol.maxX = (int) value(maxX);
////			sol.maxY = (int) value(maxY);
////			for (int i = 0; i < V.length; i++) {
////				var vi = V[i];
////				var rp = new RoomPieces();
////				rp.x0 = (int) value(vi.x0);
////				rp.y0 = (int) value(vi.y0);
////				rp.wH = (int) value(vi.wH);
////				rp.tY = (int) value(vi.tY);
////				rp.tX = (int) value(vi.tX);
////				rp.hV = (int) value(vi.hV);
////				int rightBB = (int) value(vi.rightBB);
////				int topBB = (int) value(vi.topBB);
////				rp.wBB = rightBB - rp.x0;
////				rp.hBB = topBB - rp.y0;
////				sol.piecesByNode.put(G.nodes.get(i), rp);
////			}
////			solutions.add(sol);
////
////			// Ausgabe (optional)
////			System.out.println("Solution #" + solutions.size() + " @ " + wallTime() + "s");
////
////			// HARTES LIMIT: Anzahl
////			if (solutions.size() >= maxSolutions) {
////				System.out.println("Reached maxSolutions=" + maxSolutions + ". Stopping search.");
////				stopSearch(); // <<< das beendet den CP-SAT sofort
////				return;
////			}
////
////			// OPTIONAL: Zeitbremse (falls Enumeration aus dem Ruder läuft)
////			if (maxWallTimeSeconds > 0 && wallTime() >= maxWallTimeSeconds) {
////				System.out.println(
////						"Reached wall time " + wallTime() + "s (limit " + maxWallTimeSeconds + "s). Stopping search.");
////				stopSearch();
////			}
////		}
////	}
////
////}
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
///////**************************************////////////////////////////////////////*************************************/////////////
//package pdf;
//
//import com.google.ortools.Loader;
//import com.google.ortools.sat.*;
//import com.google.ortools.sat.TableConstraint; // wichtig
//
//import java.util.*;
//
///**
// * CP-SAT Floorplan-Solver mit Google OR-Tools. Rechteckige Räume (x,y,w,h) im
// * Raster, harte Flächenvorgaben via Allowed-Tuples. - Nicht-Überlappung
// * reifiziert (ohne Big-M) - Nachbarn: bündige Adjazenz (L/R/T/B) +
// * Mindestkontakt - Nicht-Nachbarn: keine gemeinsame Kante (optional per
// * forbidUnwantedContacts), Ecken erlaubt - Ziel: minimiere maxX + maxY
// */
//public final class FloorplanCrossModalSolver {
//
//	/** Einzelnes (ggf. L-förmig darstellbares) Raumstück – hier als Rechteck. */
//	public static final class RoomPieces {
//		public int x0, y0; // linke obere Ecke (Rasterzellen)
//		public int wH, hV; // Breite/Höhe der Bounding-Box
//		public int tX, tY; // für L-Form (hier = w/h)
//		public int wBB, hBB; // Bounding-Box (hier = w/h)
//	}
//
//	/** komplette Lösung: maxX/maxY und Zuordnung Raum->Pieces */
//	public static final class Solution {
//		public int maxX, maxY;
//		public final Map<String, RoomPieces> piecesByNode = new LinkedHashMap<>();
//	}
//
//	private FloorplanCrossModalSolver() {
//	}
//
//	public static List<Solution> layout(CrossModalMapper.MappedParams p) {
//		Loader.loadNativeLibraries();
//
//		// ---- Eingabeprüfung ----
//		var G = p.G;
//		int n = G.nodes.size();
//		if (n == 0)
//			throw new IllegalArgumentException("Graph hat keine Knoten");
//		if (p.minDim > p.maxDim)
//			throw new IllegalArgumentException("minDim > maxDim");
//		if (p.hullWCells <= 0 || p.hullHCells <= 0)
//			throw new IllegalArgumentException("Ungültige Bounding-Box");
//
//		final int border = Math.max(0, p.borderCells);
//		final int hullW = p.hullWCells;
//		final int hullH = p.hullHCells;
//		final int hullWUse = Math.max(0, hullW - 2 * border);
//		final int hullHUse = Math.max(0, hullH - 2 * border);
//		if (hullWUse <= 0 || hullHUse <= 0) {
//			throw new IllegalArgumentException("Nutzbare Hülle ≤ 0 (borderCells zu groß?)");
//		}
//
//		final CpModel model = new CpModel();
//
//		// ---- Variablen je Raum ----
//		IntVar[] x = new IntVar[n];
//		IntVar[] y = new IntVar[n];
//		IntVar[] w = new IntVar[n];
//		IntVar[] h = new IntVar[n];
//
//		for (int i = 0; i < n; i++) {
//			String id = G.nodes.get(i);
//			var ro = p.perRoomOptions.get(id);
//
//			int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
//			int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
//			if (minDim > maxDim)
//				throw new IllegalArgumentException("minDim > maxDim für Raum " + id);
//
//			x[i] = model.newIntVar(border, hullW - border, "x_" + i);
//			y[i] = model.newIntVar(border, hullH - border, "y_" + i);
//			w[i] = model.newIntVar(minDim, maxDim, "w_" + i);
//			h[i] = model.newIntVar(minDim, maxDim, "h_" + i);
//
//			// Innenhülle: x+w ≤ hullW-border, y+h ≤ hullH-border
//			model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[i], w[i] }), hullW - border);
//			model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y[i], h[i] }), hullH - border);
//
//			// Fläche (optional): Allowed (w,h)-Paare in Zellen mit Aspectfenster/Toleranz
//			Integer areaSqm = p.areaSqmByRoom.get(id);
//			if (p.enforceAreas && areaSqm != null) {
//				List<long[]> tuples = allowedPairsForArea(areaSqm, p.cellMeters, minDim, maxDim, p.areaTolCells,
//						(ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect,
//						(ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect);
//				if (tuples.isEmpty()) {
//					throw new IllegalStateException(
//							"Keine (w,h)-Paare für Raum " + id + " (Fläche/Aspekt/Toleranz/Dims prüfen).");
//				}
////				model.addAllowedAssignments(new IntVar[] { w[i], h[i] }, tuples);
//				TableConstraint tab = model.addAllowedAssignments(new IntVar[] { w[i], h[i] });
//				for (long[] t : tuples) {
//				    tab.addTuple(t); // z.B. new long[]{W, H}
//				}
//
//			}
//		}
//
//		// ---- Paarweise Constraints ----
//		for (int i = 0; i < n; i++) {
//			for (int j = i + 1; j < n; j++) {
//				String idI = G.nodes.get(i);
//				String idJ = G.nodes.get(j);
//				boolean neighbors = G.areNeighbors(idI, idJ);
//
//				// Basis: 4 Richtungs-Bools + reifizierte Trennung
//				BoolVar left = model.newBoolVar("left_" + i + "_" + j);
//				BoolVar right = model.newBoolVar("right_" + i + "_" + j);
//				BoolVar above = model.newBoolVar("above_" + i + "_" + j);
//				BoolVar below = model.newBoolVar("below_" + i + "_" + j);
//
//				// Mindestens eine Richtung aktiv
//				model.addBoolOr(new Literal[] { left, right, above, below });
//
//				// Reifizierte Trennungen (bündig erlaubt)
//				model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[i], w[i] }), x[j]).onlyEnforceIf(left);
//				model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[j], w[j] }), x[i]).onlyEnforceIf(right);
//				model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y[i], h[i] }), y[j]).onlyEnforceIf(above);
//				model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y[j], h[j] }), y[i]).onlyEnforceIf(below);
//
//				if (neighbors) {
//					// --- Nachbarn: Adjazenz-Orientierung + bündig + Mindestkontakt ---
//					BoolVar L = model.newBoolVar("L_" + i + "_" + j);
//					BoolVar R = model.newBoolVar("R_" + i + "_" + j);
//					BoolVar T = model.newBoolVar("T_" + i + "_" + j);
//					BoolVar B = model.newBoolVar("B_" + i + "_" + j);
//
//					// Mindestens eine Orientierung oder genau eine – per Option
//					model.addBoolOr(new Literal[] { L, R, T, B });
//					if (!p.adjAtLeastOneSide) {
//						model.addAtMostOne(new Literal[] { L, R, T, B });
//						// zusammen mit addBoolOr ⇒ exakt eine
//					}
//
//					// Kopplung sep-Bools <-> Orientierung (verhindert Widerspruch)
//					model.addEquality(L, left);
//					model.addEquality(R, right);
//					model.addEquality(T, above);
//					model.addEquality(B, below);
//
//					// Bündig-Gleichheit
//					model.addEquality(LinearExpr.sum(new IntVar[] { x[i], w[i] }), x[j]).onlyEnforceIf(L);
//					model.addEquality(LinearExpr.sum(new IntVar[] { x[j], w[j] }), x[i]).onlyEnforceIf(R);
//					model.addEquality(LinearExpr.sum(new IntVar[] { y[i], h[i] }), y[j]).onlyEnforceIf(T);
//					model.addEquality(LinearExpr.sum(new IntVar[] { y[j], h[j] }), y[i]).onlyEnforceIf(B);
//
//					// Mindestkontakt (optional)
//					int m = Math.max(0, p.minContactCells);
//					if (m > 0) {
//						// ovY aktiv bei L/R
//						IntVar ovY = model.newIntVar(0, hullHUse, "ovY_" + i + "_" + j);
//						// ovY ≤ y_j + h_j - y_i ⇒ ovY + y_i ≤ y_j + h_j
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovY, y[i] }),
//								LinearExpr.sum(new IntVar[] { y[j], h[j] })).onlyEnforceIf(L);
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovY, y[j] }),
//								LinearExpr.sum(new IntVar[] { y[i], h[i] })).onlyEnforceIf(L);
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovY, y[i] }),
//								LinearExpr.sum(new IntVar[] { y[j], h[j] })).onlyEnforceIf(R);
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovY, y[j] }),
//								LinearExpr.sum(new IntVar[] { y[i], h[i] })).onlyEnforceIf(R);
//						model.addGreaterOrEqual(ovY, m).onlyEnforceIf(L);
//						model.addGreaterOrEqual(ovY, m).onlyEnforceIf(R);
//
//						// ovX aktiv bei T/B
//						IntVar ovX = model.newIntVar(0, hullWUse, "ovX_" + i + "_" + j);
//						// ovX ≤ x_j + w_j - x_i ⇒ ovX + x_i ≤ x_j + w_j
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovX, x[i] }),
//								LinearExpr.sum(new IntVar[] { x[j], w[j] })).onlyEnforceIf(T);
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovX, x[j] }),
//								LinearExpr.sum(new IntVar[] { x[i], w[i] })).onlyEnforceIf(T);
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovX, x[i] }),
//								LinearExpr.sum(new IntVar[] { x[j], w[j] })).onlyEnforceIf(B);
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { ovX, x[j] }),
//								LinearExpr.sum(new IntVar[] { x[i], w[i] })).onlyEnforceIf(B);
//						model.addGreaterOrEqual(ovX, m).onlyEnforceIf(T);
//						model.addGreaterOrEqual(ovX, m).onlyEnforceIf(B);
//					}
//
//				} else {
//					// --- Nicht-Nachbarn ---
//					// Keine gemeinsame Kante, Ecken erlaubt – nur wenn gewünscht:
//					if (p.forbidUnwantedContacts) {
//						// Aggregatoren H = (left ∨ right), V = (above ∨ below)
//						BoolVar H = model.newBoolVar("H_" + i + "_" + j);
//						BoolVar V = model.newBoolVar("V_" + i + "_" + j);
//						model.addBoolOr(new Literal[] { left, right, H.not() });
//						model.addImplication(left, H);
//						model.addImplication(right, H);
//						model.addBoolOr(new Literal[] { above, below, V.not() });
//						model.addImplication(above, V);
//						model.addImplication(below, V);
//
//						// wenn nur horizontal getrennt (H=1 & V=0) ⇒ mind. 1 Zelle horizontaler Abstand
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[i], w[i] }),
//								LinearExpr.sum(new IntVar[] { x[j] }).build())
//								.onlyEnforceIf(new Literal[] { left, V.not() });
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[j], w[j] }),
//								LinearExpr.sum(new IntVar[] { x[i] }).build())
//								.onlyEnforceIf(new Literal[] { right, V.not() });
//						// … mit +1 Spalt: x_i + w_i + 1 ≤ x_j bzw. x_j + w_j + 1 ≤ x_i
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[i], w[i] }),
//								LinearExpr.newBuilder().add(x[j]).add(-1).build())
//								.onlyEnforceIf(new Literal[] { left, V.not() });
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[j], w[j] }),
//								LinearExpr.newBuilder().add(x[i]).add(-1).build())
//								.onlyEnforceIf(new Literal[] { right, V.not() });
//
//						// wenn nur vertikal getrennt (V=1 & H=0) ⇒ mind. 1 Zelle vertikaler Abstand
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y[i], h[i] }),
//								LinearExpr.newBuilder().add(y[j]).add(-1).build())
//								.onlyEnforceIf(new Literal[] { above, H.not() });
//						model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y[j], h[j] }),
//								LinearExpr.newBuilder().add(y[i]).add(-1).build())
//								.onlyEnforceIf(new Literal[] { below, H.not() });
//					}
//					// sonst: nur Basis-Trennung (bündig erlaubt) reicht – keine Kante wird NICHT
//					// verboten
//				}
//			}
//		}
//
//		// ---- Bounding-Box & Ziel ----
//		IntVar maxX = model.newIntVar(0, hullW, "maxX");
//		IntVar maxY = model.newIntVar(0, hullH, "maxY");
//		for (int i = 0; i < n; i++) {
//			model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x[i], w[i] }), maxX);
//			model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y[i], h[i] }), maxY);
//		}
//		if (p.forceFillHull) {
//			model.addEquality(maxX, hullW - border);
//			model.addEquality(maxY, hullH - border);
//		}
//
//		model.minimize(LinearExpr.sum(new IntVar[] { maxX, maxY }));
//
//		// ---- Lösen ----
//		CpSolver solver = new CpSolver();
//		// Optional: Limits setzen
//		// solver.getParameters().setMaxTimeInSeconds(30.0);
//		// solver.getParameters().setNumSearchWorkers(Math.max(1,
//		// Runtime.getRuntime().availableProcessors()-1));
//		CpSolverStatus st = solver.solve(model);
//		if (!(st == CpSolverStatus.OPTIMAL || st == CpSolverStatus.FEASIBLE)) {
//			System.out.println("Keine feasible Lösung gefunden – Status: " + st);
//			return List.of();
//		}
//
//		// ---- Lösung extrahieren ----
//		Solution sol = new Solution();
//		sol.maxX = (int) solver.value(maxX);
//		sol.maxY = (int) solver.value(maxY);
//
//		for (int i = 0; i < n; i++) {
//			RoomPieces rp = new RoomPieces();
//			rp.x0 = (int) solver.value(x[i]);
//			rp.y0 = (int) solver.value(y[i]);
//			int wi = (int) solver.value(w[i]);
//			int hi = (int) solver.value(h[i]);
//			rp.wH = rp.tX = wi;
//			rp.tY = rp.hV = hi;
//			rp.wBB = wi;
//			rp.hBB = hi;
//			sol.piecesByNode.put(G.nodes.get(i), rp);
//		}
//
//		return List.of(sol);
//	}
//
//	/**
//	 * Erlaubte (w,h)-Paare in Zellen für eine Ziel-Fläche in m² inkl. Toleranz und
//	 * Aspect-Fenster.
//	 */
//	static List<long[]> allowedPairsForArea(int areaSqm, double cellM, int minDim, int maxDim, int tolCells,
//			double minAsp, double maxAsp) {
//		if (areaSqm < 0 || cellM <= 0)
//			throw new IllegalArgumentException("Ungültige Flächenparameter");
//		if (tolCells < 0)
//			throw new IllegalArgumentException("tolCells < 0");
//		if (minDim > maxDim)
//			throw new IllegalArgumentException("minDim > maxDim");
//		if (minAsp > maxAsp)
//			throw new IllegalArgumentException("minAsp > maxAsp");
//
//		int target = (int) Math.round(areaSqm / (cellM * cellM)); // Zielzellen
//		List<long[]> tuples = new ArrayList<>();
//		for (int W = minDim; W <= maxDim; W++) {
//			for (int H = minDim; H <= maxDim; H++) {
//				int cells = W * H;
//				if (Math.abs(cells - target) <= tolCells) {
//					double ar = (double) W / H;
//					if (ar >= minAsp && ar <= maxAsp) {
//						tuples.add(new long[] { W, H });
//					}
//				}
//			}
//		}
//		return tuples;
//	}
//}
package pdf;

import com.google.ortools.Loader;
import com.google.ortools.sat.*;

import java.util.*;

/**
 * CP-SAT Floorplan-Solver mit Google OR-Tools.
 * Rechteckige Räume (x,y,w,h) im Raster, harte Flächenvorgaben via Allowed-Tuples.
 * - Nicht-Überlappung reifiziert (ohne Big-M)
 * - Nachbarn: bündige Adjazenz (L/R/T/B) + Mindestkontakt
 * - Nicht-Nachbarn: keine gemeinsame Kante (optional per forbidUnwantedContacts), Ecken erlaubt
 * - Ziel: gewichtete Summe ->  (W_PRIMARY)*(maxX+maxY) + (W_SQUARE)*sum|w-h| für ausgewählte Räume
 */
public final class FloorplanCrossModalSolver {

    /** Einzelnes (ggf. L-förmig darstellbares) Raumstück – hier als Rechteck. */
    public static final class RoomPieces {
        public int x0, y0;     // linke obere Ecke (Rasterzellen)
        public int wH, hV;     // Breite/Höhe der Bounding-Box
        public int tX, tY;     // für L-Form (hier = w/h)
        public int wBB, hBB;   // Bounding-Box (hier = w/h)
    }

    /** komplette Lösung: maxX/maxY und Zuordnung Raum->Pieces */
    public static final class Solution {
        public int maxX, maxY;
        public final Map<String, RoomPieces> piecesByNode = new LinkedHashMap<>();
    }

    private FloorplanCrossModalSolver() {}

    public static List<Solution> layout(CrossModalMapper.MappedParams p) {
        Loader.loadNativeLibraries();

        // ---- Eingabeprüfung ----
        var G = p.G;
        int n = G.nodes.size();
        if (n == 0) throw new IllegalArgumentException("Graph hat keine Knoten");
        if (p.minDim > p.maxDim) throw new IllegalArgumentException("minDim > maxDim");
        if (p.hullWCells <= 0 || p.hullHCells <= 0) throw new IllegalArgumentException("Ungültige Bounding-Box");

        final int border = Math.max(0, p.borderCells);
        final int hullW = p.hullWCells;
        final int hullH = p.hullHCells;
        final int hullWUse = Math.max(0, hullW - 2 * border);
        final int hullHUse = Math.max(0, hullH - 2 * border);
        if (hullWUse <= 0 || hullHUse <= 0) {
            throw new IllegalArgumentException("Nutzbare Hülle ≤ 0 (borderCells zu groß?)");
        }

        // ---- Quadratik-Präferenz: HIER anpassen, welche Räume „quadratisch“ werden sollen ----
        final Set<String> preferSquare = new HashSet<>(List.of("A","B")); // z.B. List.of("A","C")

        // ---- Gewichte für die Zielfunktion ----
        final long W_PRIMARY = 10_000L; // groß: Kompaktheit hat Vorrang
        final long W_SQUARE  = 1L;      // klein: Quadratnähe als Tie-Break

        final CpModel model = new CpModel();

        // ---- Variablen je Raum ----
        IntVar[] x = new IntVar[n];
        IntVar[] y = new IntVar[n];
        IntVar[] w = new IntVar[n];
        IntVar[] h = new IntVar[n];

        // |w-h|-Variablen für „quadratische“ Räume
        List<IntVar> sqDiffs = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            String id = G.nodes.get(i);
            var ro = p.perRoomOptions.get(id);

            int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
            int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
            if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim für Raum " + id);

            x[i] = model.newIntVar(border, hullW - border, "x_" + i);
            y[i] = model.newIntVar(border, hullH - border, "y_" + i);
            w[i] = model.newIntVar(minDim, maxDim, "w_" + i);
            h[i] = model.newIntVar(minDim, maxDim, "h_" + i);

            // Innenhülle: x+w ≤ hullW-border, y+h ≤ hullH-border
            model.addLessOrEqual(LinearExpr.sum(new IntVar[]{x[i], w[i]}), hullW - border);
            model.addLessOrEqual(LinearExpr.sum(new IntVar[]{y[i], h[i]}), hullH - border);

            // Fläche (optional): Allowed (w,h)-Paare in Zellen mit Aspectfenster/Toleranz
            Integer areaSqm = p.areaSqmByRoom.get(id);
            if (p.enforceAreas && areaSqm != null) {
                List<long[]> tuples = allowedPairsForArea(
                        areaSqm, p.cellMeters, minDim, maxDim,
                        p.areaTolCells,
                        (ro != null && ro.minAspect != null) ? ro.minAspect : p.minAspect,
                        (ro != null && ro.maxAspect != null) ? ro.maxAspect : p.maxAspect
                );
                if (tuples.isEmpty()) {
                    throw new IllegalStateException("Keine (w,h)-Paare für Raum " + id
                            + " (Fläche/Aspekt/Toleranz/Dims prüfen).");
                }
                // Kompatibel mit älteren OR-Tools: TableConstraint + addTuple
                TableConstraint tab = model.addAllowedAssignments(new IntVar[]{w[i], h[i]});
                for (long[] t : tuples) tab.addTuple(t);
            }

            // Quadratnähe für ausgewählte Räume: d = |w-h|
            if (preferSquare.contains(id)) {
                int maxDiff = Math.max(0, (ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim)
                                          - (ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim));
                if (maxDiff == 0) maxDiff = p.maxDim - p.minDim;
                IntVar d = model.newIntVar(0, Math.max(0, maxDiff), "sqdiff_" + id);
                // d = |w - h|
                model.addAbsEquality(d, LinearExpr.newBuilder().add(w[i]).addTerm(h[i], -1).build());
                sqDiffs.add(d);
            }
        }

        // ---- Paarweise Constraints ----
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                String idI = G.nodes.get(i);
                String idJ = G.nodes.get(j);
                boolean neighbors = G.areNeighbors(idI, idJ);

                // Basis: 4 Richtungs-Bools + reifizierte Trennung (bündig erlaubt)
                BoolVar left  = model.newBoolVar("left_"  + i + "_" + j);
                BoolVar right = model.newBoolVar("right_" + i + "_" + j);
                BoolVar above = model.newBoolVar("above_" + i + "_" + j);
                BoolVar below = model.newBoolVar("below_" + i + "_" + j);

                // Mindestens eine Richtung aktiv
                model.addBoolOr(new Literal[]{left, right, above, below});

                // Reifizierte Trennungen (keine Big-Ms)
                model.addLessOrEqual(LinearExpr.sum(new IntVar[]{x[i], w[i]}), x[j]).onlyEnforceIf(left);
                model.addLessOrEqual(LinearExpr.sum(new IntVar[]{x[j], w[j]}), x[i]).onlyEnforceIf(right);
                model.addLessOrEqual(LinearExpr.sum(new IntVar[]{y[i], h[i]}), y[j]).onlyEnforceIf(above);
                model.addLessOrEqual(LinearExpr.sum(new IntVar[]{y[j], h[j]}), y[i]).onlyEnforceIf(below);

                if (neighbors) {
                    // --- Nachbarn: Adjazenz-Orientierung + bündig + Mindestkontakt ---
                    BoolVar L = model.newBoolVar("L_" + i + "_" + j);
                    BoolVar R = model.newBoolVar("R_" + i + "_" + j);
                    BoolVar T = model.newBoolVar("T_" + i + "_" + j);
                    BoolVar B = model.newBoolVar("B_" + i + "_" + j);

                    model.addBoolOr(new Literal[]{L, R, T, B});
                    if (!p.adjAtLeastOneSide) {
                        model.addAtMostOne(new Literal[]{L, R, T, B}); // zusammen mit BoolOr => genau 1
                    }

                    // Kopplung sep-Bools <-> Orientierung
                    model.addEquality(L, left);
                    model.addEquality(R, right);
                    model.addEquality(T, above);
                    model.addEquality(B, below);

                    // Bündig-Gleichheit
                    model.addEquality(LinearExpr.sum(new IntVar[]{x[i], w[i]}), x[j]).onlyEnforceIf(L);
                    model.addEquality(LinearExpr.sum(new IntVar[]{x[j], w[j]}), x[i]).onlyEnforceIf(R);
                    model.addEquality(LinearExpr.sum(new IntVar[]{y[i], h[i]}), y[j]).onlyEnforceIf(T);
                    model.addEquality(LinearExpr.sum(new IntVar[]{y[j], h[j]}), y[i]).onlyEnforceIf(B);

                    // Mindestkontakt (optional)
                    int m = Math.max(0, p.minContactCells);
                    if (m > 0) {
                        IntVar ovY = model.newIntVar(0, hullHUse, "ovY_" + i + "_" + j);
                        // ovY ≤ y_j + h_j - y_i  &  ovY ≤ y_i + h_i - y_j   (aktiv bei L/R)
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovY, y[i]}),
                                LinearExpr.sum(new IntVar[]{y[j], h[j]})).onlyEnforceIf(L);
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovY, y[j]}),
                                LinearExpr.sum(new IntVar[]{y[i], h[i]})).onlyEnforceIf(L);
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovY, y[i]}),
                                LinearExpr.sum(new IntVar[]{y[j], h[j]})).onlyEnforceIf(R);
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovY, y[j]}),
                                LinearExpr.sum(new IntVar[]{y[i], h[i]})).onlyEnforceIf(R);
                        model.addGreaterOrEqual(ovY, m).onlyEnforceIf(L);
                        model.addGreaterOrEqual(ovY, m).onlyEnforceIf(R);

                        IntVar ovX = model.newIntVar(0, hullWUse, "ovX_" + i + "_" + j);
                        // ovX ≤ x_j + w_j - x_i  &  ovX ≤ x_i + w_i - x_j   (aktiv bei T/B)
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovX, x[i]}),
                                LinearExpr.sum(new IntVar[]{x[j], w[j]})).onlyEnforceIf(T);
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovX, x[j]}),
                                LinearExpr.sum(new IntVar[]{x[i], w[i]})).onlyEnforceIf(T);
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovX, x[i]}),
                                LinearExpr.sum(new IntVar[]{x[j], w[j]})).onlyEnforceIf(B);
                        model.addLessOrEqual(LinearExpr.sum(new IntVar[]{ovX, x[j]}),
                                LinearExpr.sum(new IntVar[]{x[i], w[i]})).onlyEnforceIf(B);
                        model.addGreaterOrEqual(ovX, m).onlyEnforceIf(T);
                        model.addGreaterOrEqual(ovX, m).onlyEnforceIf(B);
                    }

                } else {
                    // --- Nicht-Nachbarn ---
                    if (p.forbidUnwantedContacts) {
                        // Aggregatoren H = (left ∨ right), V = (above ∨ below)
                        BoolVar H = model.newBoolVar("H_" + i + "_" + j);
                        BoolVar V = model.newBoolVar("V_" + i + "_" + j);
                        model.addBoolOr(new Literal[]{left, right, H.not()});
                        model.addImplication(left, H);
                        model.addImplication(right, H);
                        model.addBoolOr(new Literal[]{above, below, V.not()});
                        model.addImplication(above, V);
                        model.addImplication(below, V);

                        // Wenn nur horizontal getrennt (H=1 & V=0): mind. 1 Zelle horizontaler Abstand
                        model.addLessOrEqual(
                                LinearExpr.sum(new IntVar[]{x[i], w[i]}),
                                LinearExpr.newBuilder().add(x[j]).add(-1).build()
                        ).onlyEnforceIf(new Literal[]{left, V.not()});
                        model.addLessOrEqual(
                                LinearExpr.sum(new IntVar[]{x[j], w[j]}),
                                LinearExpr.newBuilder().add(x[i]).add(-1).build()
                        ).onlyEnforceIf(new Literal[]{right, V.not()});

                        // Wenn nur vertikal getrennt (V=1 & H=0): mind. 1 Zelle vertikaler Abstand
                        model.addLessOrEqual(
                                LinearExpr.sum(new IntVar[]{y[i], h[i]}),
                                LinearExpr.newBuilder().add(y[j]).add(-1).build()
                        ).onlyEnforceIf(new Literal[]{above, H.not()});
                        model.addLessOrEqual(
                                LinearExpr.sum(new IntVar[]{y[j], h[j]}),
                                LinearExpr.newBuilder().add(y[i]).add(-1).build()
                        ).onlyEnforceIf(new Literal[]{below, H.not()});
                    }
                    // sonst: nur Basis-Trennung (bündig erlaubt)
                }
            }
        }

        // ---- Bounding-Box (rechts/oben) ----
        IntVar maxX = model.newIntVar(0, hullW, "maxX");
        IntVar maxY = model.newIntVar(0, hullH, "maxY");
        for (int i = 0; i < n; i++) {
            model.addLessOrEqual(LinearExpr.sum(new IntVar[]{x[i], w[i]}), maxX);
            model.addLessOrEqual(LinearExpr.sum(new IntVar[]{y[i], h[i]}), maxY);
        }
        if (p.forceFillHull) {
            model.addEquality(maxX, hullW - border);
            model.addEquality(maxY, hullH - border);
        }

        // ---- Gewichtete Zielfunktion ----
//        LinearExpr.Builder obj = LinearExpr.newBuilder();
//        obj.addTerm(maxX, W_PRIMARY);
//        obj.addTerm(maxY, W_PRIMARY);
//        for (IntVar d : sqDiffs) obj.addTerm(d, W_SQUARE);
//        model.minimize(obj.build());
        List<LinearArgument> terms = new ArrayList<>();
        terms.add(LinearExpr.term(maxX, W_PRIMARY));
        terms.add(LinearExpr.term(maxY, W_PRIMARY));
        for (IntVar d : sqDiffs) terms.add(LinearExpr.term(d, W_SQUARE));
        model.minimize(LinearExpr.sum(terms.toArray(new LinearArgument[0])));

        // ---- Lösen ----
        CpSolver solver = new CpSolver();
        // Beispiel-Parameter:
        // solver.getParameters().setMaxTimeInSeconds(30.0);
        // solver.getParameters().setNumSearchWorkers(Math.max(1, Runtime.getRuntime().availableProcessors()-1));
        CpSolverStatus st = solver.solve(model);
        if (!(st == CpSolverStatus.OPTIMAL || st == CpSolverStatus.FEASIBLE)) {
            System.out.println("Keine feasible Lösung gefunden – Status: " + st);
            return List.of();
        }

        // ---- Lösung extrahieren ----
        Solution sol = new Solution();
        sol.maxX = (int) solver.value(maxX);
        sol.maxY = (int) solver.value(maxY);

        for (int i = 0; i < n; i++) {
            RoomPieces rp = new RoomPieces();
            rp.x0 = (int) solver.value(x[i]);
            rp.y0 = (int) solver.value(y[i]);
            int wi = (int) solver.value(w[i]);
            int hi = (int) solver.value(h[i]);
            rp.wH = rp.tX = wi;
            rp.tY = rp.hV = hi;
            rp.wBB = wi;
            rp.hBB = hi;
            sol.piecesByNode.put(G.nodes.get(i), rp);
        }

        return List.of(sol);
    }

    /** Erlaubte (w,h)-Paare in Zellen für eine Ziel-Fläche in m² inkl. Toleranz und Aspect-Fenster. */
    static List<long[]> allowedPairsForArea(int areaSqm,
                                            double cellM,
                                            int minDim, int maxDim,
                                            int tolCells,
                                            double minAsp, double maxAsp) {
        if (areaSqm < 0 || cellM <= 0) throw new IllegalArgumentException("Ungültige Flächenparameter");
        if (tolCells < 0) throw new IllegalArgumentException("tolCells < 0");
        if (minDim > maxDim) throw new IllegalArgumentException("minDim > maxDim");
        if (minAsp > maxAsp) throw new IllegalArgumentException("minAsp > maxAsp");

        int target = (int) Math.round(areaSqm / (cellM * cellM)); // Zielzellen
        List<long[]> tuples = new ArrayList<>();
        for (int W = minDim; W <= maxDim; W++) {
            for (int H = minDim; H <= maxDim; H++) {
                int cells = W * H;
                if (Math.abs(cells - target) <= tolCells) {
                    double ar = (double) W / H;
                    if (ar >= minAsp && ar <= maxAsp) {
                        tuples.add(new long[]{W, H});
                    }
                }
            }
        }
        return tuples;
    }
}
