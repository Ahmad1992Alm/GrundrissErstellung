package pdf;

import java.util.*;

public final class CrossModalSpec {

	// --- Raumform --------------------------------------------------------------
	public enum Shape {
		RECT, L, RECT_OR_L
	}

	// --- Boundary --------------------------------------------------------------
	public static final class Boundary {
		public double widthMeters; // z. B. 7.5
		public double heightMeters; // z. B. 7.0
		public double cellMeters = 0.5; // Rastergröße
		public int borderCells = 0; // innenliegender Rand

		public Boundary() {
		}

		public Boundary(double w, double h, double cell) {
			this.widthMeters = w;
			this.heightMeters = h;
			this.cellMeters = cell;
		}

		public Boundary setBorderCells(int b) {
			this.borderCells = b;
			return this;
		}
	}

	// --- Room -----------------------------------------------------------------
	public static final class Room {
		public String id; // Pflicht
		public String type; // optional
		public Double areaSqm; // optional Ziel-Fläche
		public Double minAspect; // (nur für reine Rechteckräume relevant)
		public Double maxAspect;
		public Integer minDimCells; // min Breite/Höhe in Zellen (rechteckig)
		public Integer maxDimCells; // max Breite/Höhe in Zellen (rechteckig)

		// L-Form Parameter
		public Shape shape = Shape.RECT_OR_L;
		public Integer minThicknessCells = 1; // min Dicke tX/tY
		public Integer maxThicknessCells = 4; // max Dicke tX/tY

		public Room() {
		}

		public Room(String id) {
			this.id = id;
		}

		public Room id(String v) {
			this.id = v;
			return this;
		}

		public Room type(String v) {
			this.type = v;
			return this;
		}

		public Room area(double v) {
			this.areaSqm = v;
			return this;
		}

		public Room aspect(double min, double max) {
			this.minAspect = min;
			this.maxAspect = max;
			return this;
		}

		public Room dims(Integer min, Integer max) {
			this.minDimCells = min;
			this.maxDimCells = max;
			return this;
		}

		public Room shape(Shape s) {
			this.shape = s;
			return this;
		}

		public Room thickness(Integer minT, Integer maxT) {
			this.minThicknessCells = minT;
			this.maxThicknessCells = maxT;
			return this;
		}
	}

	// --- Edges ----------------------------------------------------------------
	public enum Relation {
		ADJACENT, NOT_ADJACENT
	}

	public static final class Edge {
		public String u;
		public String v;
		public Relation rel = Relation.ADJACENT;
		public Integer minContactCells;

		public Edge() {
		}

		public Edge(String u, String v) {
			this.u = u;
			this.v = v;
		}

		public Edge rel(Relation r) {
			this.rel = r;
			return this;
		}

		public Edge minContact(int c) {
			this.minContactCells = c;
			return this;
		}
	}

	// --- ModelOptions ---------------------------------------------------------
	public static final class ModelOptions {
		/** Auswahl des zugrunde liegenden Solvers. */
		public enum Solver {
			CP_SAT, MIQP
		}

		/** Welcher Solver soll verwendet werden? */
		public Solver solver = Solver.CP_SAT;

		public int minDimCells = 1;
		public int maxDimCells = 10;
//		public boolean enforceAreas = false;
		public boolean enforceAreas = true;
		public int areaTolCells = 1;
		public double minAspect = 0.2;
		public double maxAspect = 5.0;
//		public int minContactCells = 0;
		public int minContactCells = 1;
		public boolean forbidUnwantedContacts = true; // Nicht-Nachbarn dürfen sich nicht berühren
//		public boolean forceFillHull = false; // maxX==W && maxY==H
		public boolean forceFillHull = true; // maxX==W && maxY==H
		public boolean adjAtLeastOneSide = true; // statt "genau 1 Seite"
	}

	// --- Spec -----------------------------------------------------------------
	public static final class Spec {
		public final List<Room> rooms = new ArrayList<>();
		public final List<Edge> edges = new ArrayList<>();
		public Boundary boundary = new Boundary(7.5, 7.0, 0.5);
		public ModelOptions options = new ModelOptions();

		public Spec addRoom(Room r) {
			rooms.add(r);
			return this;
		}

		public Spec addEdge(Edge e) {
			edges.add(e);
			return this;
		}

		public Spec boundary(Boundary b) {
			this.boundary = b;
			return this;
		}

		public Spec options(ModelOptions o) {
			this.options = o;
			return this;
		}
	}
}
