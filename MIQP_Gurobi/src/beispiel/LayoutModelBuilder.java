
package beispiel;

import com.google.ortools.sat.*;
import java.util.*;

// ===== Datenklassen (package-sichtbar, nicht public) =====
enum WallSide {
	ANY, LEFT, RIGHT, TOP, BOTTOM, NONE
}

class RoomOptions {
	Double aspectRatioMin; // z.B. 1.2 (null = aus)
	Double aspectRatioMax; // z.B. 2.0 (null = aus)
	WallSide mustTouch = WallSide.NONE;
	Double wallToleranceM; // z.B. 0.2 (null = exakt)
}

//class Adjacency {
//	String other;
//	int weight = 1; // alpha_ij >= 1
//	boolean requireDoor = false;
//}

// class Room {
//	String name;
//	// Zielmaße (weich), Meter
//	double wStarM, hStarM;
//	// harte Bounds (Meter)
//	double wMinM, wMaxM, hMinM, hMaxM;
//	List<Adjacency> neighbors = new ArrayList<>();
//	RoomOptions options = new RoomOptions();
//}

/**
 * ===== Öffentliche Builder-Klasse (Dateiname: LayoutModelBuilder.java) =====
 */
public class LayoutModelBuilder {
	// === Skalierung ===
	public static final int SCALE = 10; // 1 Einheit = 0.1 m

	private final java.util.Set<String> requiredAdj = new java.util.HashSet<>();

	// === Hülle ===
	final double WM, HM; // Meter
	final int W, H; // skaliert (Int)

	// === Optionaler Korridor ===
	boolean useVerticalCorridor = false;
	double corridorWidthM = 1.2;
	int corridorWidth;

	// Türbreite (für requireDoor)
	double doorMinM = 0.9;
	int doorMin;

	// Zielfunktions-Gewicht
	int beta = 2;

	// Eingabe
	final List<Room> rooms;

	// OR-Tools
	final CpModel model = new CpModel();

	// Variablen
	final Map<String, IntVar> x = new HashMap<>();
	final Map<String, IntVar> y = new HashMap<>();
	final Map<String, IntVar> w = new HashMap<>();
	final Map<String, IntVar> h = new HashMap<>();
	final Map<String, IntVar> cx = new HashMap<>();
	final Map<String, IntVar> cy = new HashMap<>();

	final Map<String, BoolVar> s = new HashMap<>(); // i links von j
	final Map<String, BoolVar> u = new HashMap<>(); // i unter j

	private String undirectedKey(String a, String b) {
		return (a.compareTo(b) < 0) ? (a + "||" + b) : (b + "||" + a);
	}

	public void requireAdjacency(String a, String b) {
		requiredAdj.add(undirectedKey(a, b));
	}

	// Verbotene Zonen (Cut-outs)
	static class Rect {
		final int x, y, w, h;

		Rect(int x, int y, int w, int h) {
			this.x = x;
			this.y = y;
			this.w = w;
			this.h = h;
		}
	}

	final List<Rect> forbidden = new ArrayList<>();

	// === Helfer für Keys und LinearExpr ===
	private String key(String a, String b) {
		return a + "||" + b;
	}

	private static LinearExpr lin(IntVar v) {
		return LinearExpr.term(v, 1);
	}

	private static LinearExpr sum(LinearExpr a, LinearExpr b) {
		return LinearExpr.sum(new LinearExpr[] { a, b });
	}

	private static LinearExpr sum(LinearExpr a, LinearExpr b, LinearExpr c) {
		return LinearExpr.sum(new LinearExpr[] { a, b, c });
	}

	private static LinearExpr sum(LinearExpr a, LinearExpr b, LinearExpr c, LinearExpr d) {
		return LinearExpr.sum(new LinearExpr[] { a, b, c, d });
	}

	// === API für Cut-outs in Metern ===
	public void addForbiddenRectMeters(double xM, double yM, double wM, double hM) {
		forbidden.add(new Rect((int) Math.round(xM * SCALE), (int) Math.round(yM * SCALE), (int) Math.round(wM * SCALE),
				(int) Math.round(hM * SCALE)));
	}

	public LayoutModelBuilder(double wMeters, double hMeters, List<Room> rooms) {
		this.WM = wMeters;
		this.HM = hMeters;
		this.W = (int) Math.round(WM * SCALE);
		this.H = (int) Math.round(HM * SCALE);
		this.rooms = rooms;
		this.corridorWidth = (int) Math.round(corridorWidthM * SCALE);
		this.doorMin = (int) Math.round(doorMinM * SCALE);
	}

	public void setVerticalCorridor(boolean use) {
		this.useVerticalCorridor = use;
	}

	// ================= Build =================
	public CpModel build() {
		defineVarsAndHull();
		addForbiddenZones(); // Weg A: Cut-outs
		addNonOverlap();
		addCenters();
		addAspectRatios();
		addWallRequirements();
		if (useVerticalCorridor)
			addVerticalCorridor();
		addObjective();
		return model;
	}

	// ===== Variablen & Hülle =====
	private void defineVarsAndHull() {
		for (Room r : rooms) {
			IntVar xi = model.newIntVar(0, W, "x_" + r.name);
			IntVar yi = model.newIntVar(0, H, "y_" + r.name);
			IntVar wi = model.newIntVar((int) Math.round(r.wMinM * SCALE), (int) Math.round(r.wMaxM * SCALE),
					"w_" + r.name);
			IntVar hi = model.newIntVar((int) Math.round(r.hMinM * SCALE), (int) Math.round(r.hMaxM * SCALE),
					"h_" + r.name);
			x.put(r.name, xi);
			y.put(r.name, yi);
			w.put(r.name, wi);
			h.put(r.name, hi);

			// x + w <= W ; y + h <= H
			model.addLessOrEqual(sum(lin(xi), lin(wi)), LinearExpr.constant(W));
			model.addLessOrEqual(sum(lin(yi), lin(hi)), LinearExpr.constant(H));
		}

	}

	// ===== Weg A: verbotene Zonen =====
	private void addForbiddenZones() {
		for (Room r : rooms) {
			IntVar xi = x.get(r.name), yi = y.get(r.name), wi = w.get(r.name), hi = h.get(r.name);
			for (int k = 0; k < forbidden.size(); k++) {
				Rect f = forbidden.get(k);

				BoolVar L = model.newBoolVar("F_L_" + r.name + "_" + k);
				BoolVar Rr = model.newBoolVar("F_R_" + r.name + "_" + k);
				BoolVar B = model.newBoolVar("F_B_" + r.name + "_" + k);
				BoolVar T = model.newBoolVar("F_T_" + r.name + "_" + k);

				// L+R+B+T >= 1
				LinearExpr atLeastOne = sum(sum(lin(L), lin(Rr)), sum(lin(B), lin(T)));
				model.addLinearConstraint(atLeastOne, 1, 4);

				// links: x + w <= f.x
				model.addLessOrEqual(sum(lin(xi), lin(wi)), LinearExpr.constant(f.x)).onlyEnforceIf(L);
				// rechts: f.x + f.w <= x
				model.addLessOrEqual(LinearExpr.constant(f.x + f.w), lin(xi)).onlyEnforceIf(Rr);
				// unten: y + h <= f.y
				model.addLessOrEqual(sum(lin(yi), lin(hi)), LinearExpr.constant(f.y)).onlyEnforceIf(B);
				// oben: f.y + f.h <= y
				model.addLessOrEqual(LinearExpr.constant(f.y + f.h), lin(yi)).onlyEnforceIf(T);
			}
		}
	}

	private void addNonOverlap() {
		final int Mx = W; // Big-M in X (groß genug: >= W)
		final int My = H; // Big-M in Y (groß genug: >= H)

		for (int i = 0; i < rooms.size(); i++) {
			for (int j = i + 1; j < rooms.size(); j++) {
				Room ri = rooms.get(i), rj = rooms.get(j);

				IntVar xi = x.get(ri.name), yi = y.get(ri.name), wi = w.get(ri.name), hi = h.get(ri.name);
				IntVar xj = x.get(rj.name), yj = y.get(rj.name), wj = w.get(rj.name), hj = h.get(rj.name);

				// Vier Disjunktions-Booleans:
				// L_ij: i komplett links von j
				// R_ij: i komplett rechts von j
				// B_ij: i komplett unter j
				// T_ij: i komplett über j
				BoolVar L = model.newBoolVar("L_" + ri.name + "_" + rj.name);
				BoolVar R = model.newBoolVar("R_" + ri.name + "_" + rj.name);
				BoolVar B = model.newBoolVar("B_" + ri.name + "_" + rj.name);
				BoolVar T = model.newBoolVar("T_" + ri.name + "_" + rj.name);

				// Mindestens eine Richtung aktiv (>=1). Wenn du "genau eine" willst, setz 1..1.
				model.addLinearConstraint(LinearExpr.sum(new LinearExpr[] { LinearExpr.term(L, 1),
						LinearExpr.term(R, 1), LinearExpr.term(B, 1), LinearExpr.term(T, 1) }), 1, 4);

				// Links: x_i + w_i <= x_j + Mx*(1 - L)
				model.addLessOrEqual(
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xi, 1), LinearExpr.term(wi, 1) }),
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xj, 1), LinearExpr.term(L, -Mx),
								LinearExpr.constant(Mx) }));

				// Rechts: x_j + w_j <= x_i + Mx*(1 - R)
				model.addLessOrEqual(
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xj, 1), LinearExpr.term(wj, 1) }),
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xi, 1), LinearExpr.term(R, -Mx),
								LinearExpr.constant(Mx) }));

				// Unten: y_i + h_i <= y_j + My*(1 - B)
				model.addLessOrEqual(
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yi, 1), LinearExpr.term(hi, 1) }),
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yj, 1), LinearExpr.term(B, -My),
								LinearExpr.constant(My) }));

				// Oben: y_j + h_j <= y_i + My*(1 - T)
				model.addLessOrEqual(
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yj, 1), LinearExpr.term(hj, 1) }),
						LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yi, 1), LinearExpr.term(T, -My),
								LinearExpr.constant(My) }));

				
				
				boolean mustTouch = requiredAdj.contains(undirectedKey(ri.name, rj.name));
				if (mustTouch) {
				  // Genau EINE relative Lage (stabil)
				  model.addLinearConstraint(
				    LinearExpr.sum(new LinearExpr[]{
				      LinearExpr.term(L,1), LinearExpr.term(R,1),
				      LinearExpr.term(B,1), LinearExpr.term(T,1)
				    }),
				    1, 1
				  );

				  // Kanten-Kontakt erzwingen je nach Lage:
				  // L: x_i + w_i = x_j
				  model.addEquality(
				    LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(xi,1), LinearExpr.term(wi,1) }),
				    LinearExpr.term(xj,1)
				  ).onlyEnforceIf(L);

				  // R: x_j + w_j = x_i
				  model.addEquality(
				    LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(xj,1), LinearExpr.term(wj,1) }),
				    LinearExpr.term(xi,1)
				  ).onlyEnforceIf(R);

				  // B: y_i + h_i = y_j
				  model.addEquality(
				    LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(yi,1), LinearExpr.term(hi,1) }),
				    LinearExpr.term(yj,1)
				  ).onlyEnforceIf(B);

				  // T: y_j + h_j = y_i
				  model.addEquality(
				    LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(yj,1), LinearExpr.term(hj,1) }),
				    LinearExpr.term(yi,1)
				  ).onlyEnforceIf(T);

				  // UND: sinnvolle Tür-Überlappung aktivieren (mind. doorMin)
				  // (nutzt deinen bestehenden Tür-Block; Stelle sicher, dass edgeHasDoor(...) true ist)
				}

				
				
				// (Optional, aber empfehlenswert) Kopplung X vs Y, um "schräge" Freischaltungen
				// zu vermeiden:
				// Wenn horizontal getrennt (L oder R), darf vertikale Trennung theoretisch 0
				// sein (kein Zwang).
				// Wenn vertikal getrennt (B oder T), darf horizontale Trennung 0 sein.
				// D.h. KEINE zusätzliche Kopplung nötig — die Big-M-Formeln reichen.

//	      // TÜR-BEDINGUNG (falls gefordert):
//	      if (edgeHasDoor(ri, rj)) {
//	        // Tür bei horizontaler Nachbarschaft (L oder R): vertikale Überlappung >= doorMin
//	        // Wir modellieren das separat mit zwei Implikationen:
//
//	        // Wenn L=1: y_i <= y_j + h_j - doorMin  UND  y_j <= y_i + h_i - doorMin
//	        model.addLessOrEqual(
//	          LinearExpr.term(yi, 1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(yj,1), LinearExpr.term(hj,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(L,  My) }) // wenn L=0, relax um My
//	        );
//	        model.addLessOrEqual(
//	          LinearExpr.term(yj, 1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(yi,1), LinearExpr.term(hi,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(L,  My) })
//	        );
//
//	        // Wenn R=1: gleiche Bedingung
//	        model.addLessOrEqual(
//	          LinearExpr.term(yi, 1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(yj,1), LinearExpr.term(hj,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(R,  My) })
//	        );
//	        model.addLessOrEqual(
//	          LinearExpr.term(yj, 1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(yi,1), LinearExpr.term(hi,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(R,  My) })
//	        );
//
//	        // Tür bei vertikaler Nachbarschaft (B oder T): horizontale Überlappung >= doorMin
//	        model.addLessOrEqual(
//	          LinearExpr.term(xi,1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(xj,1), LinearExpr.term(wj,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(B,  Mx) })
//	        );
//	        model.addLessOrEqual(
//	          LinearExpr.term(xj,1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(xi,1), LinearExpr.term(wi,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(B,  Mx) })
//	        );
//
//	        model.addLessOrEqual(
//	          LinearExpr.term(xi,1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(xj,1), LinearExpr.term(wj,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(T,  Mx) })
//	        );
//	        model.addLessOrEqual(
//	          LinearExpr.term(xj,1),
//	          LinearExpr.sum(new LinearExpr[]{ LinearExpr.term(xi,1), LinearExpr.term(wi,1), LinearExpr.constant(-doorMin),
//	                                           LinearExpr.term(T,  Mx) })
//	        );
//	      }
				// Türbreite, falls gefordert
				if (edgeHasDoor(ri, rj)) {

					// Wenn L=1 ODER R=1 => vertikale Überlappung >= doorMin
					// y_i <= y_j + h_j - doorMin + My*(1 - L)
					model.addLessOrEqual(LinearExpr.term(yi, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yj, 1), LinearExpr.term(hj, 1),
									LinearExpr.constant(-doorMin),
									// + My*(1 - L) = + My - My*L
									LinearExpr.constant(My), LinearExpr.term(L, -My) }));
					// y_j <= y_i + h_i - doorMin + My*(1 - L)
					model.addLessOrEqual(LinearExpr.term(yj, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yi, 1), LinearExpr.term(hi, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(My), LinearExpr.term(L, -My) }));
					// Gleiches für R:
					model.addLessOrEqual(LinearExpr.term(yi, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yj, 1), LinearExpr.term(hj, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(My), LinearExpr.term(R, -My) }));
					model.addLessOrEqual(LinearExpr.term(yj, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(yi, 1), LinearExpr.term(hi, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(My), LinearExpr.term(R, -My) }));

					// Wenn B=1 ODER T=1 => horizontale Überlappung >= doorMin
					// x_i <= x_j + w_j - doorMin + Mx*(1 - B)
					model.addLessOrEqual(LinearExpr.term(xi, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xj, 1), LinearExpr.term(wj, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(Mx), LinearExpr.term(B, -Mx) }));
					model.addLessOrEqual(LinearExpr.term(xj, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xi, 1), LinearExpr.term(wi, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(Mx), LinearExpr.term(B, -Mx) }));
					// Gleiches für T:
					model.addLessOrEqual(LinearExpr.term(xi, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xj, 1), LinearExpr.term(wj, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(Mx), LinearExpr.term(T, -Mx) }));
					model.addLessOrEqual(LinearExpr.term(xj, 1),
							LinearExpr.sum(new LinearExpr[] { LinearExpr.term(xi, 1), LinearExpr.term(wi, 1),
									LinearExpr.constant(-doorMin), LinearExpr.constant(Mx), LinearExpr.term(T, -Mx) }));
				}

			}
		}
	}

//  // ===== Nicht-Überlappung & Türbreite =====
//  private void addNonOverlap(){
//    for (int i=0;i<rooms.size();i++){
//      for (int j=i+1;j<rooms.size();j++){
//        Room ri = rooms.get(i), rj = rooms.get(j);
//        IntVar xi = x.get(ri.name), yi = y.get(ri.name), wi = w.get(ri.name), hi = h.get(ri.name);
//        IntVar xj = x.get(rj.name), yj = y.get(rj.name), wj = w.get(rj.name), hj = h.get(rj.name);
//
//        BoolVar sij = model.newBoolVar("s_"+ri.name+"_"+rj.name); // i links von j
//        BoolVar uij = model.newBoolVar("u_"+ri.name+"_"+rj.name); // i unter j
//        s.put(key(ri.name,rj.name), sij);
//        u.put(key(ri.name,rj.name), uij);
//
//        // Links/Rechts:
//        // (sij=1) -> x_i + w_i <= x_j + W
//        model.addLessOrEqual( sum(lin(xi), lin(wi)), sum(lin(xj), LinearExpr.constant(W)) ).onlyEnforceIf(sij);
//        // allgemein: x_j + w_j <= x_i + W*sij
//        model.addLessOrEqual( sum(lin(xj), lin(wj)), sum(lin(xi), LinearExpr.term(sij, W)) );
//
//        // Oben/Unten:
//        model.addLessOrEqual( sum(lin(yi), lin(hi)), sum(lin(yj), LinearExpr.constant(H)) ).onlyEnforceIf(uij);
//        model.addLessOrEqual( sum(lin(yj), lin(hj)), sum(lin(yi), LinearExpr.term(uij, H)) );
//
//        // s_ij + u_ij >= 1
//        model.addLinearConstraint( sum(lin(sij), lin(uij)), 1, 2 );
//
//        // Türbreite, falls gefordert
//        if (edgeHasDoor(ri, rj)){
//          // wenn i links von j → vertikale Überlappung >= doorMin
//          model.addLessOrEqual( lin(yi), sum(lin(yj), lin(hj), LinearExpr.constant(-doorMin)) ).onlyEnforceIf(sij);
//          model.addLessOrEqual( lin(yj), sum(lin(yi), lin(hi), LinearExpr.constant(-doorMin)) ).onlyEnforceIf(sij);
//          // wenn i unter j → horizontale Überlappung >= doorMin
//          model.addLessOrEqual( lin(xi), sum(lin(xj), lin(wj), LinearExpr.constant(-doorMin)) ).onlyEnforceIf(uij);
//          model.addLessOrEqual( lin(xj), sum(lin(xi), lin(wi), LinearExpr.constant(-doorMin)) ).onlyEnforceIf(uij);
//        }
//      }
//    }
//  }

	private boolean edgeHasDoor(Room a, Room b) {
		for (Adjacency adj : a.neighbors)
			if (adj.other.equals(b.name) && adj.requireDoor)
				return true;
		for (Adjacency adj : b.neighbors)
			if (adj.other.equals(a.name) && adj.requireDoor)
				return true;
		return false;
	}

	// ===== Zentren (2*cx = 2*x + w; 2*cy = 2*y + h) =====
	private void addCenters() {
		for (Room r : rooms) {
			IntVar cxi = model.newIntVar(0, W, "cx_" + r.name);
			IntVar cyi = model.newIntVar(0, H, "cy_" + r.name);
			cx.put(r.name, cxi);
			cy.put(r.name, cyi);

			LinearExpr cxEq = LinearExpr.newBuilder().addTerm(cxi, 2).addTerm(x.get(r.name), -2)
					.addTerm(w.get(r.name), -1).build();
			LinearExpr cyEq = LinearExpr.newBuilder().addTerm(cyi, 2).addTerm(y.get(r.name), -2)
					.addTerm(h.get(r.name), -1).build();
			model.addEquality(cxEq, 0);
			model.addEquality(cyEq, 0);
		}
	}

	// ===== Harte Seitenverhältnisse =====
	private void addAspectRatios() {
		for (Room r : rooms) {
			if (r.options == null)
				continue;
			if (r.options.aspectRatioMin != null) {
				int K = 10, AR_MIN = (int) Math.round(r.options.aspectRatioMin * K);
				model.addGreaterOrEqual(LinearExpr.term(w.get(r.name), K), LinearExpr.term(h.get(r.name), AR_MIN));
			}
			if (r.options.aspectRatioMax != null) {
				int K = 10, AR_MAX = (int) Math.round(r.options.aspectRatioMax * K);
				model.addLessOrEqual(LinearExpr.term(w.get(r.name), K), LinearExpr.term(h.get(r.name), AR_MAX));
			}
		}
	}

	// ===== Außenwand-Pflichten =====
	private void addWallRequirements() {
		for (Room r : rooms) {
			if (r.options == null || r.options.mustTouch == WallSide.NONE)
				continue;
			IntVar xi = x.get(r.name), yi = y.get(r.name), wi = w.get(r.name), hi = h.get(r.name);
			int delta = (r.options.wallToleranceM == null) ? 0 : (int) Math.round(r.options.wallToleranceM * SCALE);

			switch (r.options.mustTouch) {
			case LEFT:
				if (delta == 0)
					model.addEquality(lin(xi), LinearExpr.constant(0));
				else
					model.addLessOrEqual(lin(xi), LinearExpr.constant(delta));
				break;
			case RIGHT:
				if (delta == 0)
					model.addEquality(sum(lin(xi), lin(wi)), LinearExpr.constant(W));
				else
					model.addGreaterOrEqual(sum(lin(xi), lin(wi)), LinearExpr.constant(W - delta));
				break;
			case BOTTOM:
				if (delta == 0)
					model.addEquality(lin(yi), LinearExpr.constant(0));
				else
					model.addLessOrEqual(lin(yi), LinearExpr.constant(delta));
				break;
			case TOP:
				if (delta == 0)
					model.addEquality(sum(lin(yi), lin(hi)), LinearExpr.constant(H));
				else
					model.addGreaterOrEqual(sum(lin(yi), lin(hi)), LinearExpr.constant(H - delta));
				break;
			case ANY:
				BoolVar bL = model.newBoolVar("bL_" + r.name);
				BoolVar bR = model.newBoolVar("bR_" + r.name);
				BoolVar bB = model.newBoolVar("bB_" + r.name);
				BoolVar bT = model.newBoolVar("bT_" + r.name);
				model.addLinearConstraint(sum(sum(lin(bL), lin(bR)), sum(lin(bB), lin(bT))), 1, 4);
				if (delta == 0) {
					model.addEquality(lin(xi), LinearExpr.constant(0)).onlyEnforceIf(bL);
					model.addEquality(sum(lin(xi), lin(wi)), LinearExpr.constant(W)).onlyEnforceIf(bR);
					model.addEquality(lin(yi), LinearExpr.constant(0)).onlyEnforceIf(bB);
					model.addEquality(sum(lin(yi), lin(hi)), LinearExpr.constant(H)).onlyEnforceIf(bT);
				} else {
					model.addLessOrEqual(lin(xi), LinearExpr.constant(delta)).onlyEnforceIf(bL);
					model.addGreaterOrEqual(sum(lin(xi), lin(wi)), LinearExpr.constant(W - delta)).onlyEnforceIf(bR);
					model.addLessOrEqual(lin(yi), LinearExpr.constant(delta)).onlyEnforceIf(bB);
					model.addGreaterOrEqual(sum(lin(yi), lin(hi)), LinearExpr.constant(H - delta)).onlyEnforceIf(bT);
				}
				break;
			}
		}
	}

	// ===== Optional: vertikaler Korridor =====
	private void addVerticalCorridor() {
		IntVar xC = model.newIntVar(0, W, "x_C");
		LinearExpr wC = LinearExpr.constant(corridorWidth);
		// xC + wC <= W
		model.addLessOrEqual(sum(lin(xC), wC), LinearExpr.constant(W));

		for (Room r : rooms) {
			BoolVar leftOf = model.newBoolVar("leftOfC_" + r.name);
			BoolVar rightOf = model.newBoolVar("rightOfC_" + r.name);
			model.addLinearConstraint(sum(lin(leftOf), lin(rightOf)), 1, 1);

			// links: x_r + w_r <= xC
			model.addLessOrEqual(sum(lin(x.get(r.name)), lin(h.get(r.name)) /* <- FALSCH: sollte w !!! */ ), lin(xC))
					.onlyEnforceIf(leftOf);
			// ^^^ ACHTUNG: korrigieren wir unten im Hinweis!

			// rechts: xC + wC <= x_r
			model.addLessOrEqual(sum(lin(xC), wC), lin(x.get(r.name))).onlyEnforceIf(rightOf);
		}
	}

	// ===== Zielfunktion =====
//  private void addObjective(){
//    LinearExpr.Builder obj = LinearExpr.newBuilder();
//
//    // |w - w*| + |h - h*|
//    for (Room r : rooms){
//      int wStar = (int)Math.round(r.wStarM * SCALE);
//      int hStar = (int)Math.round(r.hStarM * SCALE);
//
//      IntVar wdiff = model.newIntVar(-W, W, "wdiff_"+r.name);
//      IntVar hdiff = model.newIntVar(-H, H, "hdiff_"+r.name);
//
//      model.addEquality( lin(wdiff), sum(lin(w.get(r.name)), LinearExpr.constant(-wStar)) );
//      model.addEquality( lin(hdiff), sum(lin(h.get(r.name)), LinearExpr.constant(-hStar)) );
//
//      IntVar dw = model.newIntVar(0, W, "absw_"+r.name);
//      IntVar dh = model.newIntVar(0, H, "absh_"+r.name);
//
//      model.addGreaterOrEqual( lin(dw), lin(wdiff) );
//      model.addGreaterOrEqual( lin(dw), LinearExpr.term(wdiff, -1) );
//      model.addGreaterOrEqual( lin(dh), lin(hdiff) );
//      model.addGreaterOrEqual( lin(dh), LinearExpr.term(hdiff, -1) );
//
//      obj.addTerm(dw, beta);
//      obj.addTerm(dh, beta);
//    }
//
//    // ∑ alpha_ij (|Δcx| + |Δcy|)
//    for (int i=0;i<rooms.size();i++){
//      Room ri = rooms.get(i);
//      for (Adjacency adj : ri.neighbors){
//        int j = indexOf(adj.other);
//        if (j <= i) continue;
//        int a = Math.max(1, adj.weight);
//
//        IntVar sdx = model.newIntVar(-W, W, "sdx_"+ri.name+"_"+adj.other);
//        IntVar sdy = model.newIntVar(-H, H, "sdy_"+ri.name+"_"+adj.other);
//        IntVar dx  = model.newIntVar(0, W,  "dx_"+ri.name+"_"+adj.other);
//        IntVar dy  = model.newIntVar(0, H,  "dy_"+ri.name+"_"+adj.other);
//
//        model.addEquality( lin(sdx), sum(lin(cx.get(ri.name)), LinearExpr.term(cx.get(adj.other), -1)) );
//        model.addEquality( lin(sdy), sum(lin(cy.get(ri.name)), LinearExpr.term(cy.get(adj.other), -1)) );
//
//        model.addGreaterOrEqual( lin(dx), lin(sdx) );
//        model.addGreaterOrEqual( lin(dx), LinearExpr.term(sdx, -1) );
//        model.addGreaterOrEqual( lin(dy), lin(sdy) );
//        model.addGreaterOrEqual( lin(dy), LinearExpr.term(sdy, -1) );
//
//        obj.addTerm(dx, a);
//        obj.addTerm(dy, a);
//      }
//    }
//
//    model.minimize(obj.build());
//  }
	private void addObjective() {
		java.util.List<LinearExpr> terms = new java.util.ArrayList<>();

		// (i) |w - w*| + |h - h*|
		for (Room r : rooms) {
			int wStar = (int) Math.round(r.wStarM * SCALE);
			int hStar = (int) Math.round(r.hStarM * SCALE);

			IntVar wdiff = model.newIntVar(-W, W, "wdiff_" + r.name);
			IntVar hdiff = model.newIntVar(-H, H, "hdiff_" + r.name);

			model.addEquality(LinearExpr.term(wdiff, 1), LinearExpr
					.sum(new LinearExpr[] { LinearExpr.term(w.get(r.name), 1), LinearExpr.constant(-wStar) }));
			model.addEquality(LinearExpr.term(hdiff, 1), LinearExpr
					.sum(new LinearExpr[] { LinearExpr.term(h.get(r.name), 1), LinearExpr.constant(-hStar) }));

			IntVar dw = model.newIntVar(0, W, "absw_" + r.name);
			IntVar dh = model.newIntVar(0, H, "absh_" + r.name);

			model.addGreaterOrEqual(LinearExpr.term(dw, 1), LinearExpr.term(wdiff, 1));
			model.addGreaterOrEqual(LinearExpr.term(dw, 1), LinearExpr.term(wdiff, -1));
			model.addGreaterOrEqual(LinearExpr.term(dh, 1), LinearExpr.term(hdiff, 1));
			model.addGreaterOrEqual(LinearExpr.term(dh, 1), LinearExpr.term(hdiff, -1));

			// Gewichtung via coeff in term(...)
			terms.add(LinearExpr.term(dw, beta));
			terms.add(LinearExpr.term(dh, beta));
		}

		// (ii) ∑ alpha_ij (|Δcx| + |Δcy|)
		for (int i = 0; i < rooms.size(); i++) {
			Room ri = rooms.get(i);
			for (Adjacency adj : ri.neighbors) {
				int j = indexOf(adj.other);
				if (j <= i)
					continue;
				int a = Math.max(1, adj.weight);

				IntVar sdx = model.newIntVar(-W, W, "sdx_" + ri.name + "_" + adj.other);
				IntVar sdy = model.newIntVar(-H, H, "sdy_" + ri.name + "_" + adj.other);
				IntVar dx = model.newIntVar(0, W, "dx_" + ri.name + "_" + adj.other);
				IntVar dy = model.newIntVar(0, H, "dy_" + ri.name + "_" + adj.other);

				model.addEquality(LinearExpr.term(sdx, 1), LinearExpr.sum(new LinearExpr[] {
						LinearExpr.term(cx.get(ri.name), 1), LinearExpr.term(cx.get(adj.other), -1) }));
				model.addEquality(LinearExpr.term(sdy, 1), LinearExpr.sum(new LinearExpr[] {
						LinearExpr.term(cy.get(ri.name), 1), LinearExpr.term(cy.get(adj.other), -1) }));

				model.addGreaterOrEqual(LinearExpr.term(dx, 1), LinearExpr.term(sdx, 1));
				model.addGreaterOrEqual(LinearExpr.term(dx, 1), LinearExpr.term(sdx, -1));
				model.addGreaterOrEqual(LinearExpr.term(dy, 1), LinearExpr.term(sdy, 1));
				model.addGreaterOrEqual(LinearExpr.term(dy, 1), LinearExpr.term(sdy, -1));

				terms.add(LinearExpr.term(dx, a));
				terms.add(LinearExpr.term(dy, a));
			}
		}

		// Ziel setzen
		LinearExpr objective = terms.isEmpty() ? LinearExpr.constant(0)
				: LinearExpr.sum(terms.toArray(new LinearExpr[0]));
		model.minimize(objective);
	}

	private int indexOf(String name) {
		for (int k = 0; k < rooms.size(); k++)
			if (rooms.get(k).name.equals(name))
				return k;
		return -1;
	}

	// ===== Getter (praktisch fürs Auslesen) =====
	public Map<String, IntVar> getX() {
		return x;
	}

	public Map<String, IntVar> getY() {
		return y;
	}

	public Map<String, IntVar> getW() {
		return w;
	}

	public Map<String, IntVar> getH() {
		return h;
	}

	public Map<String, IntVar> getCX() {
		return cx;
	}

	public Map<String, IntVar> getCY() {
		return cy;
	}
}
