//package pdf;
//
//import com.google.ortools.Loader;
//import com.google.ortools.sat.*;
//
//public final class OrtoolsFloorplanSolver {
//	public static void main(String[] args) {
//		Loader.loadNativeLibraries();
//		CpModel model = new CpModel();
//
//		int hullW = 15, hullH = 14, border = 0;
//		int minDim = 2, maxDim = 10;
//		double cellM = 0.5;
//
//		// Beispiel: zwei Räume mit Fläche 12 m² (48 Zellen)
//		int targetCells = (int) Math.round(12.0 / (cellM * cellM)); // 48
//
//		IntVar x1 = model.newIntVar(border, hullW - border, "x1");
//		IntVar y1 = model.newIntVar(border, hullH - border, "y1");
//		IntVar w1 = model.newIntVar(minDim, maxDim, "w1");
//		IntVar h1 = model.newIntVar(minDim, maxDim, "h1");
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x1, w1 }), hullW - border);
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y1, h1 }), hullH - border);
//		// Fläche via Table:
//		java.util.List<long[]> t1 = java.util.List.of(new long[] { 6, 8 }, new long[] { 8, 6 });
//		model.addAllowedAssignments(new IntVar[] { w1, h1 }, t1);
//
//		IntVar x2 = model.newIntVar(border, hullW - border, "x2");
//		IntVar y2 = model.newIntVar(border, hullH - border, "y2");
//		IntVar w2 = model.newIntVar(minDim, maxDim, "w2");
//		IntVar h2 = model.newIntVar(minDim, maxDim, "h2");
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x2, w2 }), hullW - border);
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y2, h2 }), hullH - border);
//		java.util.List<long[]> t2 = java.util.List.of(new long[] { 4, 12 }, new long[] { 6, 8 }, new long[] { 8, 6 });
//		model.addAllowedAssignments(new IntVar[] { w2, h2 }, t2);
//
//		// Nicht-Nachbarn: keine Kante, Ecke erlaubt
//		BoolVar left = model.newBoolVar("left_1_2");
//		BoolVar right = model.newBoolVar("right_1_2");
//		BoolVar above = model.newBoolVar("above_1_2");
//		BoolVar below = model.newBoolVar("below_1_2");
//		model.addBoolOr(new Literal[] { left, right, above, below });
//		// Reifizierte Trennungen:
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x1, w1 }), x2).onlyEnforceIf(left);
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x2, w2 }), x1).onlyEnforceIf(right);
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y1, h1 }), y2).onlyEnforceIf(above);
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y2, h2 }), y1).onlyEnforceIf(below);
//		// H/V-Aggregatoren + 1-Zellen-Gap wenn nur eine Achse trennt:
//		BoolVar H = model.newBoolVar("H_1_2");
//		BoolVar V = model.newBoolVar("V_1_2");
//		model.addBoolOr(new Literal[] { left, right, H.not() });
//		model.addImplication(left, H);
//		model.addImplication(right, H);
//		model.addBoolOr(new Literal[] { above, below, V.not() });
//		model.addImplication(above, V);
//		model.addImplication(below, V);
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x1, w1 }), LinearExpr.newBuilder().add(x2).add(-1).build())
//				.onlyEnforceIf(new Literal[] { left, V.not() });
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { x2, w2 }), LinearExpr.newBuilder().add(x1).add(-1).build())
//				.onlyEnforceIf(new Literal[] { right, V.not() });
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y1, h1 }), LinearExpr.newBuilder().add(y2).add(-1).build())
//				.onlyEnforceIf(new Literal[] { above, H.not() });
//		model.addLessOrEqual(LinearExpr.sum(new IntVar[] { y2, h2 }), LinearExpr.newBuilder().add(y1).add(-1).build())
//				.onlyEnforceIf(new Literal[] { below, H.not() });
//
//		// Zielfunktion (optional): minimiere maxX+maxY
//		IntVar maxX = model.newIntVar(0, hullW, "maxX");
//		IntVar maxY = model.newIntVar(0, hullH, "maxY");
//		model.addGreaterOrEqual(maxX, LinearExpr.sum(new IntVar[] { x1, w1 }));
//		model.addGreaterOrEqual(maxX, LinearExpr.sum(new IntVar[] { x2, w2 }));
//		model.addGreaterOrEqual(maxY, LinearExpr.sum(new IntVar[] { y1, h1 }));
//		model.addGreaterOrEqual(maxY, LinearExpr.sum(new IntVar[] { y2, h2 }));
//		model.minimize(LinearExpr.sum(new IntVar[] { maxX, maxY }));
//
//		CpSolver solver = new CpSolver();
//		CpSolverStatus st = solver.solve(model);
//		System.out.println(st);
//		if (st == CpSolverStatus.OPTIMAL || st == CpSolverStatus.FEASIBLE) {
//			System.out.printf("R1 x=%d y=%d w=%d h=%d%n", solver.value(x1), solver.value(y1), solver.value(w1),
//					solver.value(h1));
//			System.out.printf("R2 x=%d y=%d w=%d h=%d%n", solver.value(x2), solver.value(y2), solver.value(w2),
//					solver.value(h2));
//		}
//	}
//}
