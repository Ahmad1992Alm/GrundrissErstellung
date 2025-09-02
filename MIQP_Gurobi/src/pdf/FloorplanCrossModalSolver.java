package pdf;

import com.google.ortools.Loader;
import com.google.ortools.sat.*;
import java.util.*;

public final class FloorplanCrossModalSolver {

	// ---------- Hilfsfunktionen für LinearExpr ----------
	private static LinearExpr S(LinearArgument... terms) {
		return LinearExpr.newBuilder().addSum(terms).build();
	}

	private static LinearExpr Splus(long constant, LinearArgument... terms) {
		LinearExprBuilder b = LinearExpr.newBuilder();
		b.add(constant);
		for (LinearArgument t : terms)
			b.add(t);
		return b.build();
	}

	// ---------- Variablen pro Raum (L-Form als 2 Balken, gemeinsamer Ursprung)
	// --------
	static final class VarsL {
		IntVar x0, y0; // gemeinsame Ecke
		// Horizontaler Balken (wH x tY)
		IntVar wH, tY, rightH, topH;
		// Vertikaler Balken (tX x hV)
		IntVar tX, hV, rightV, topV;
		// Bounding Box
		IntVar rightBB, topBB;
		// Fläche
		IntVar p1, p2, p3, areaCells;
	}

	// ---------- Lösung ----------
	public static final class RoomPieces {
		public int x0, y0, wH, tY, tX, hV, wBB, hBB;
	}

	public static final class Solution {
		public final Map<String, RoomPieces> piecesByNode = new LinkedHashMap<>();
		public int maxX, maxY;

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append("maxX=").append(maxX).append(", maxY=").append(maxY).append('\n');
			for (var e : piecesByNode.entrySet()) {
				var p = e.getValue();
				int areaCells = p.wH * p.tY + p.tX * p.hV - p.tX * p.tY;
				double areaSqm = areaCells * 1.0; // echte m² skaliert der Viewer
				sb.append(String.format("%s -> x0=%d y0=%d H:(%d×%d) V:(%d×%d) BB:(%d×%d) cells=%d\n", e.getKey(), p.x0,
						p.y0, p.wH, p.tY, p.tX, p.hV, p.wBB, p.hBB, areaCells));
			}
			return sb.toString();
		}
	}

	// ---------- Callback zum Einsammeln ----------
	static final class MultipleSolutionsCallback extends CpSolverSolutionCallback {
		private final List<Solution> solutions = new ArrayList<>();
		private final VarsL[] V;
		private final CrossModalMapper.AdjacencyGraph G;
		private final IntVar maxX, maxY;

		MultipleSolutionsCallback(VarsL[] V, CrossModalMapper.AdjacencyGraph G, IntVar maxX, IntVar maxY) {
			this.V = V;
			this.G = G;
			this.maxX = maxX;
			this.maxY = maxY;
		}

		public List<Solution> getSolutions() {
			return solutions;
		}

		@Override
		public void onSolutionCallback() {
			Solution sol = new Solution();
			sol.maxX = (int) value(maxX);
			sol.maxY = (int) value(maxY);
			for (int i = 0; i < V.length; i++) {
				var vi = V[i];
				var rp = new RoomPieces();
				rp.x0 = (int) value(vi.x0);
				rp.y0 = (int) value(vi.y0);
				rp.wH = (int) value(vi.wH);
				rp.tY = (int) value(vi.tY);
				rp.tX = (int) value(vi.tX);
				rp.hV = (int) value(vi.hV);
				int rightBB = (int) value(vi.rightBB);
				int topBB = (int) value(vi.topBB);
				rp.wBB = rightBB - rp.x0;
				rp.hBB = topBB - rp.y0;
				sol.piecesByNode.put(G.nodes.get(i), rp);
			}
			solutions.add(sol);
		}
	}

	// ---------- Utilities ----------
	private static BoolVar ltBool(CpModel m, IntVar a, IntVar b, String name) {
		BoolVar z = m.newBoolVar(name);
		m.addLessOrEqual(Splus(1, a), b).onlyEnforceIf(z); // a+1 <= b ⇒ a < b
		m.addGreaterOrEqual(a, b).onlyEnforceIf(z.not()); // a >= b ⇒ ¬(a<b)
		return z;
	}

	private static void noOverlap(CpModel m, IntVar ax, IntVar ay, IntVar aw, IntVar ah, IntVar bx, IntVar by,
			IntVar bw, IntVar bh, int sep, String tag) {
		BoolVar left = m.newBoolVar("left_" + tag);
		BoolVar right = m.newBoolVar("right_" + tag);
		BoolVar above = m.newBoolVar("above_" + tag);
		BoolVar below = m.newBoolVar("below_" + tag);
		m.addEquality(LinearExpr.sum(new LinearArgument[] { left, right, above, below }), LinearExpr.constant(1));
		// ax+aw <= bx - sep (A links von B)
		m.addLessOrEqual(S(ax, aw), Splus(-sep, bx)).onlyEnforceIf(left);
		// bx+bw <= ax - sep
		m.addLessOrEqual(S(bx, bw), Splus(-sep, ax)).onlyEnforceIf(right);
		// ay+ah <= by - sep
		m.addLessOrEqual(S(ay, ah), Splus(-sep, by)).onlyEnforceIf(above);
		// by+bh <= ay - sep
		m.addLessOrEqual(S(by, bh), Splus(-sep, ay)).onlyEnforceIf(below);
	}

	private static void contactAtLeastOneSide(CpModel m, IntVar ax, IntVar ay, IntVar aw, IntVar ah, IntVar aright,
			IntVar atop, IntVar bx, IntVar by, IntVar bw, IntVar bh, IntVar bright, IntVar btop, int minContactCells,
			List<BoolVar> collectors, String tag) {
		// Vier Richtungen A↔B, je ein Bool; wir fügen die Bools dem collectors hinzu
		BoolVar L = m.newBoolVar("L_" + tag); // A rechts an B links
		BoolVar R = m.newBoolVar("R_" + tag); // A links an B rechts
		BoolVar T = m.newBoolVar("T_" + tag); // A unten an B oben
		BoolVar B = m.newBoolVar("B_" + tag); // A oben an B unten

		// L: aright == bx & y-Intervalle überlappen mind. minContact
		m.addEquality(aright, bx).onlyEnforceIf(L);
		m.addLessOrEqual(ay, Splus(-minContactCells, by, bh)).onlyEnforceIf(L);
		m.addLessOrEqual(by, Splus(-minContactCells, ay, ah)).onlyEnforceIf(L);

		// R: bright == ax
		m.addEquality(bright, ax).onlyEnforceIf(R);
		m.addLessOrEqual(ay, Splus(-minContactCells, by, bh)).onlyEnforceIf(R);
		m.addLessOrEqual(by, Splus(-minContactCells, ay, ah)).onlyEnforceIf(R);

		// T: atop == by & x-Überlappung >= minContact
		m.addEquality(atop, by).onlyEnforceIf(T);
		m.addLessOrEqual(ax, Splus(-minContactCells, bx, bw)).onlyEnforceIf(T);
		m.addLessOrEqual(bx, Splus(-minContactCells, ax, aw)).onlyEnforceIf(T);

		// B: btop == ay
		m.addEquality(btop, ay).onlyEnforceIf(B);
		m.addLessOrEqual(ax, Splus(-minContactCells, bx, bw)).onlyEnforceIf(B);
		m.addLessOrEqual(bx, Splus(-minContactCells, ax, aw)).onlyEnforceIf(B);

		collectors.add(L);
		collectors.add(R);
		collectors.add(T);
		collectors.add(B);
	}

	// ---------- Hauptmethode --------------------------------------------------
	public static List<Solution> layout(CrossModalMapper.MappedParams p) {
		Loader.loadNativeLibraries();
		CpModel model = new CpModel();

		var G = p.G;
		if (G.nodes.isEmpty())
			throw new IllegalArgumentException("Graph hat keine Knoten.");
		int n = G.nodes.size();

		VarsL[] V = new VarsL[n];

		// Variablen je Raum
		for (int i = 0; i < n; i++) {
			V[i] = new VarsL();
			var ro = p.perRoomOptions.get(G.nodes.get(i));
			int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
			int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
			int minT = ro != null && ro.minThicknessCells != null ? ro.minThicknessCells : 1;
			int maxT = ro != null && ro.maxThicknessCells != null ? ro.maxThicknessCells : Math.min(4, maxDim);

			V[i].x0 = model.newIntVar(0, p.hullWCells, "x0_" + i);
			V[i].y0 = model.newIntVar(0, p.hullHCells, "y0_" + i);

			V[i].wH = model.newIntVar(minDim, maxDim, "wH_" + i);
			V[i].tY = model.newIntVar(minT, maxT, "tY_" + i);
			V[i].tX = model.newIntVar(minT, maxT, "tX_" + i);
			V[i].hV = model.newIntVar(minDim, maxDim, "hV_" + i);

			V[i].rightH = model.newIntVar(0, p.hullWCells, "rightH_" + i);
			V[i].topH = model.newIntVar(0, p.hullHCells, "topH_" + i);
			V[i].rightV = model.newIntVar(0, p.hullWCells, "rightV_" + i);
			V[i].topV = model.newIntVar(0, p.hullHCells, "topV_" + i);
			model.addEquality(V[i].rightH, S(V[i].x0, V[i].wH));
			model.addEquality(V[i].topH, S(V[i].y0, V[i].tY));
			model.addEquality(V[i].rightV, S(V[i].x0, V[i].tX));
			model.addEquality(V[i].topV, S(V[i].y0, V[i].hV));

			V[i].rightBB = model.newIntVar(0, p.hullWCells, "rightBB_" + i);
			V[i].topBB = model.newIntVar(0, p.hullHCells, "topBB_" + i);
			model.addMaxEquality(V[i].rightBB, new IntVar[] { V[i].rightH, V[i].rightV });
			model.addMaxEquality(V[i].topBB, new IntVar[] { V[i].topH, V[i].topV });

			// Hülle einhalten
			model.addGreaterOrEqual(V[i].x0, p.borderCells);
			model.addGreaterOrEqual(V[i].y0, p.borderCells);
			
		
			
			model.addLessOrEqual(V[i].rightBB, LinearExpr.constant(p.hullWCells - p.borderCells));
			model.addLessOrEqual(V[i].topBB, LinearExpr.constant(p.hullHCells - p.borderCells));

			// Fläche: areaCells = wH*tY + tX*hV - tX*tY
			V[i].p1 = model.newIntVar(0, p.hullWCells * p.hullHCells, "p1_" + i);
			V[i].p2 = model.newIntVar(0, p.hullWCells * p.hullHCells, "p2_" + i);
			V[i].p3 = model.newIntVar(0, p.hullWCells * p.hullHCells, "p3_" + i);
			V[i].areaCells = model.newIntVar(0, p.hullWCells * p.hullHCells, "area_" + i);
			model.addMultiplicationEquality(V[i].p1, V[i].wH, V[i].tY);
			model.addMultiplicationEquality(V[i].p2, V[i].tX, V[i].hV);
			model.addMultiplicationEquality(V[i].p3, V[i].tX, V[i].tY);
			LinearExprBuilder ab = LinearExpr.newBuilder();
			ab.add(V[i].p1);
			ab.add(V[i].p2);
			ab.addTerm(V[i].p3, -1);
			model.addEquality(V[i].areaCells, ab.build());

			// Shape-Restriktionen
			CrossModalSpec.Shape shape = ro != null && ro.shape != null ? ro.shape : CrossModalSpec.Shape.RECT_OR_L;
			if (shape == CrossModalSpec.Shape.RECT) {
				model.addEquality(V[i].tX, V[i].wH);
				model.addEquality(V[i].tY, V[i].hV);
			} else if (shape == CrossModalSpec.Shape.L) {
				BoolVar bx = ltBool(model, V[i].tX, V[i].wH, "tX_lt_wH_" + i);
				BoolVar by = ltBool(model, V[i].tY, V[i].hV, "tY_lt_hV_" + i);
				model.addEquality(LinearExpr.sum(new LinearArgument[] { bx, by }), LinearExpr.constant(2)); // beide
																											// echt
																											// kleiner
//				model.addGreaterOrEqual(V[i].wH, Splus(1, V[i].tX));
//				model.addGreaterOrEqual(V[i].hV, Splus(1, V[i].tY));

			} // RECT_OR_L: keine Zusatzbindung
		}
		// Flächenvorgaben

		// maxX/maxY (BB) optional auf Hülle festklemmen
		IntVar maxX = model.newIntVar(0, p.hullWCells, "maxX");
		IntVar maxY = model.newIntVar(0, p.hullHCells, "maxY");
		model.addMaxEquality(maxX, Arrays.stream(V).map(v -> v.rightBB).toArray(IntVar[]::new));
		model.addMaxEquality(maxY, Arrays.stream(V).map(v -> v.topBB).toArray(IntVar[]::new));
		if (p.forceFillHull) {
			model.addEquality(maxX, p.hullWCells);
			model.addEquality(maxY, p.hullHCells);
		}

		// Anker (optional): ersten Raum oben links fixieren
		model.addEquality(V[0].x0, p.borderCells);
		model.addEquality(V[0].y0, p.borderCells);

		// Paare
		List<int[]> pairs = new ArrayList<>();
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				pairs.add(new int[] { i, j });

		// Adjazenzen
		for (int[] pr : pairs) {
			int i = pr[0], j = pr[1];
			String idI = G.nodes.get(i), idJ = G.nodes.get(j);
			if (G.areNeighbors(idI, idJ)) {
				List<BoolVar> contactBools = new ArrayList<>();
				// (i.H) vs (j.H)
				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[i].rightH, V[i].topH, V[j].x0,
						V[j].y0, V[j].wH, V[j].tY, V[j].rightH, V[j].topH, p.minContactCells, contactBools,
						"iH_jH_" + i + "_" + j);
				// (i.H) vs (j.V)
				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[i].rightH, V[i].topH, V[j].x0,
						V[j].y0, V[j].tX, V[j].hV, V[j].rightV, V[j].topV, p.minContactCells, contactBools,
						"iH_jV_" + i + "_" + j);
				// (i.V) vs (j.H)
				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[i].rightV, V[i].topV, V[j].x0,
						V[j].y0, V[j].wH, V[j].tY, V[j].rightH, V[j].topH, p.minContactCells, contactBools,
						"iV_jH_" + i + "_" + j);
				// (i.V) vs (j.V)
				contactAtLeastOneSide(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[i].rightV, V[i].topV, V[j].x0,
						V[j].y0, V[j].tX, V[j].hV, V[j].rightV, V[j].topV, p.minContactCells, contactBools,
						"iV_jV_" + i + "_" + j);

				// Insgesamt mindestens eine Seitenberührung zwischen i und j
				if (p.adjAtLeastOneSide) {
					mAddGreaterOrEqualSum(model, contactBools, 1);
				} else {
					// selten sinnvoll mit L-Formen, aber möglich: genau 1 Seite
					model.addEquality(LinearExpr.sum(contactBools.toArray(new LinearArgument[0])),
							LinearExpr.constant(1));
				}
			}
		}

		// Nicht-Nachbarn trennen (alle Stück-Kombinationen)
		int sep = p.forbidUnwantedContacts ? 1 : 0;
		for (int[] pr : pairs) {
			int i = pr[0], j = pr[1];
			String idI = G.nodes.get(i), idJ = G.nodes.get(j);
			if (!G.areNeighbors(idI, idJ)) {
				noOverlap(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[j].x0, V[j].y0, V[j].wH, V[j].tY, sep,
						"iH_jH_" + i + "_" + j);
				noOverlap(model, V[i].x0, V[i].y0, V[i].wH, V[i].tY, V[j].x0, V[j].y0, V[j].tX, V[j].hV, sep,
						"iH_jV_" + i + "_" + j);
				noOverlap(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[j].x0, V[j].y0, V[j].wH, V[j].tY, sep,
						"iV_jH_" + i + "_" + j);
				noOverlap(model, V[i].x0, V[i].y0, V[i].tX, V[i].hV, V[j].x0, V[j].y0, V[j].tX, V[j].hV, sep,
						"iV_jV_" + i + "_" + j);
			}
		}

		// Solver
		CpSolver solver = new CpSolver();
		solver.getParameters().setEnumerateAllSolutions(true);
		solver.getParameters().setNumSearchWorkers(1);
		solver.getParameters().setMaxTimeInSeconds(30.0);

		MultipleSolutionsCallback cb = new MultipleSolutionsCallback(V, G, maxX, maxY);
		CpSolverStatus status = solver.solve(model, cb);
		if (status != CpSolverStatus.OPTIMAL && status != CpSolverStatus.FEASIBLE)
			throw new RuntimeException("Kein Layout gefunden: " + status);
		return cb.getSolutions();
	}

	private static void mAddGreaterOrEqualSum(CpModel m, List<BoolVar> bs, int k) {
		LinearArgument[] arr = bs.toArray(new LinearArgument[0]);
		m.addGreaterOrEqual(LinearExpr.sum(arr), LinearExpr.constant(k));
	}
}