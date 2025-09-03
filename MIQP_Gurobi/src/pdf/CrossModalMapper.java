package pdf;

import java.util.*;
//
///**
// * Mappt eine {@link CrossModalSpec.Spec} in ein Parameter-Bündel, das der
// * Solver konsumiert.
// */
//public final class CrossModalMapper {
//	// ---- Graph für den Solver -------------------------------------------------
//	public static class AdjacencyGraph {
//		public final List<String> nodes = new ArrayList<>();
//		public final Map<String, Set<String>> nbrs = new HashMap<>();
//
//		public void addNode(String id) {
//			if (!nbrs.containsKey(id)) {
//				nbrs.put(id, new HashSet<>());
//				nodes.add(id);
//			}
//		}
//
//		public void addUndirectedEdge(String a, String b) {
//			if (!a.equals(b)) {
//				addNode(a);
//				addNode(b);
//				nbrs.get(a).add(b);
//				nbrs.get(b).add(a);
//			}
//		}
//
//		public boolean areNeighbors(String a, String b) {
//			return nbrs.getOrDefault(a, Set.of()).contains(b);
//		}
//	}
//
//	// ---- Per-Raum Overrides, inkl. Shape & Thickness -------------------------
//	public static final class RoomLocalOptions {
//		public final Integer minDimCells;
//		public final Integer maxDimCells;
//		public final Double minAspect;
//		public final Double maxAspect;
//		public final CrossModalSpec.Shape shape;
//		public final Integer minThicknessCells;
//		public final Integer maxThicknessCells;
//
//		public RoomLocalOptions(Integer minDim, Integer maxDim, Double minAsp, Double maxAsp,
//				CrossModalSpec.Shape shape, Integer minT, Integer maxT) {
//			this.minDimCells = minDim;
//			this.maxDimCells = maxDim;
//			this.minAspect = minAsp;
//			this.maxAspect = maxAsp;
//			this.shape = shape;
//			this.minThicknessCells = minT;
//			this.maxThicknessCells = maxT;
//		}
//	}
//
//	// ---- MappedParams: alles, was der Solver braucht -------------------------
//	public static final class MappedParams {
//		public final CrossModalSpec.ModelOptions.Solver solver;
//		public final int hullWCells, hullHCells, borderCells;
//		public final double cellMeters;
//		public final int minDim, maxDim;
//		public final boolean enforceAreas;
//		public final int areaTolCells;
//		public final double minAspect, maxAspect;
//		public final int minContactCells;
//		public final boolean forbidUnwantedContacts;
//		public final boolean forceFillHull;
//		public final boolean adjAtLeastOneSide;
//		public final Map<String, Integer> areaSqmByRoom;
//		public final Map<String, RoomLocalOptions> perRoomOptions;
//		public final AdjacencyGraph G;
//
//		MappedParams(CrossModalSpec.ModelOptions.Solver solver, int w, int h, int border, double cell, int minDim,
//				int maxDim, boolean enforceAreas, int areaTol, double minAsp, double maxAsp, int minContact,
//				boolean forbid, boolean forceFill, boolean adjAtLeast, Map<String, Integer> areaMap,
//				Map<String, RoomLocalOptions> roomOpts, AdjacencyGraph g) {
//			this.solver = solver;
//			this.hullWCells = w;
//			this.hullHCells = h;
//			this.borderCells = border;
//			this.cellMeters = cell;
//			this.minDim = minDim;
//			this.maxDim = maxDim;
//			this.enforceAreas = enforceAreas;
//			this.areaTolCells = areaTol;
//			this.minAspect = minAsp;
//			this.maxAspect = maxAsp;
//			this.minContactCells = minContact;
//			this.forbidUnwantedContacts = forbid;
//			this.forceFillHull = forceFill;
//			this.adjAtLeastOneSide = adjAtLeast;
//			this.areaSqmByRoom = Collections.unmodifiableMap(areaMap);
//			this.perRoomOptions = Collections.unmodifiableMap(roomOpts);
//			this.G = g;
//		}
//	}
//
//	private static int toCells(double meters, double cellMeters) {
//		return (int) Math.round(meters / cellMeters);
//	}
//
//	public static MappedParams map(CrossModalSpec.Spec spec) {
//		Objects.requireNonNull(spec, "spec");
//		var b = spec.boundary;
//		var o = spec.options;
//		int W = toCells(b.widthMeters, b.cellMeters);
//		int H = toCells(b.heightMeters, b.cellMeters);
//
//		Map<String, Integer> areaMap = new LinkedHashMap<>();
//		Map<String, RoomLocalOptions> perRoom = new HashMap<>();
//
//		var G = new AdjacencyGraph();
////		for (var r : spec.rooms) {
////			G.addNode(r.id);
////			if (r.areaSqm != null)
////				areaMap.put(r.id, (int) Math.round(r.areaSqm));
////			perRoom.put(r.id, new RoomLocalOptions(r.minDimCells, r.maxDimCells, r.minAspect, r.maxAspect, r.shape,
////					r.minThicknessCells, r.maxThicknessCells));
////		}
//		for (var r : spec.rooms) {
//			G.addNode(r.id);
////			if (r.areaSqm != null) {
////				double cellArea = b.cellMeters * b.cellMeters;
////				int areaCells = (int) Math.round(r.areaSqm / cellArea);
////				areaMap.put(r.id, areaCells);
////			}
//			if (r.areaSqm != null) {
//			    areaMap.put(r.id, (int) Math.round(r.areaSqm)); // m² speichern
//			}
//			perRoom.put(r.id, new RoomLocalOptions(r.minDimCells, r.maxDimCells, r.minAspect, r.maxAspect, r.shape,
//					r.minThicknessCells, r.maxThicknessCells));
//		}
//		for (var e : spec.edges) {
//			if (e.rel == CrossModalSpec.Relation.ADJACENT)
//				G.addUndirectedEdge(e.u, e.v);
//		}
//
//        return new MappedParams(o.solver, W, H, b.borderCells, b.cellMeters, o.minDimCells, o.maxDimCells,
//                o.enforceAreas, o.areaTolCells, o.minAspect, o.maxAspect, o.minContactCells,
//                o.forbidUnwantedContacts, o.forceFillHull, o.adjAtLeastOneSide, areaMap, perRoom, G);
//}
//
//}


import java.util.Objects;

/**
 * Mappt eine {@link CrossModalSpec.Spec} in ein Parameter-Bündel, das der Solver konsumiert.
 * - Speichert Raumflächen (falls vorhanden) in m² (Integer, gerundet).
 * - Rechnet Boundary-Maße in Zellen um (auf Basis von cellMeters).
 * - Baut einen ungerichteten Adjazenzgraphen.
 *
 * Hinweise:
 * - Der MIQPFloorplanSolver erwartet areaSqmByRoom in m² (Integer); die Umrechnung in Zellen
 *   erfolgt solverseitig (unter Einbezug von cellMeters und areaTolCells).
 */
public final class CrossModalMapper {

    // ---- Graph für den Solver -------------------------------------------------
    public static class AdjacencyGraph {
        public final List<String> nodes = new ArrayList<>();
        public final Map<String, Set<String>> nbrs = new HashMap<>();

        public void addNode(String id) {
            Objects.requireNonNull(id, "node id");
            if (!nbrs.containsKey(id)) {
                nbrs.put(id, new HashSet<>());
                nodes.add(id);
            }
        }

        public void addUndirectedEdge(String a, String b) {
            Objects.requireNonNull(a, "edge u");
            Objects.requireNonNull(b, "edge v");
            if (!a.equals(b)) {
                addNode(a);
                addNode(b);
                nbrs.get(a).add(b);
                nbrs.get(b).add(a);
            }
        }

        public boolean areNeighbors(String a, String b) {
            return nbrs.getOrDefault(a, Set.of()).contains(b);
        }
    }

    // ---- Per-Raum Overrides, inkl. Shape & Thickness -------------------------
    public static final class RoomLocalOptions {
        public final Integer minDimCells;
        public final Integer maxDimCells;
        public final Double minAspect;
        public final Double maxAspect;
        public final CrossModalSpec.Shape shape;
        public final Integer minThicknessCells;
        public final Integer maxThicknessCells;

        public RoomLocalOptions(Integer minDim, Integer maxDim, Double minAsp, Double maxAsp,
                                CrossModalSpec.Shape shape, Integer minT, Integer maxT) {
            this.minDimCells = minDim;
            this.maxDimCells = maxDim;
            this.minAspect = minAsp;
            this.maxAspect = maxAsp;
            this.shape = shape;
            this.minThicknessCells = minT;
            this.maxThicknessCells = maxT;
        }
    }

    // ---- MappedParams: alles, was der Solver braucht -------------------------
    public static final class MappedParams {
        public final CrossModalSpec.ModelOptions.Solver solver;
        public final int hullWCells, hullHCells, borderCells;
        public final double cellMeters;
        public final int minDim, maxDim;
        public final boolean enforceAreas;
        public final int areaTolCells;
        public final double minAspect, maxAspect;
        public final int minContactCells;
        public final boolean forbidUnwantedContacts;
        public final boolean forceFillHull;
        public final boolean adjAtLeastOneSide;
        public final Map<String, Integer> areaSqmByRoom;           // m², Integer (gerundet)
        public final Map<String, RoomLocalOptions> perRoomOptions;
        public final AdjacencyGraph G;

        MappedParams(CrossModalSpec.ModelOptions.Solver solver, int w, int h, int border, double cell,
                     int minDim, int maxDim, boolean enforceAreas, int areaTol, double minAsp, double maxAsp,
                     int minContact, boolean forbid, boolean forceFill, boolean adjAtLeast,
                     Map<String, Integer> areaMap, Map<String, RoomLocalOptions> roomOpts, AdjacencyGraph g) {
            this.solver = solver;
            this.hullWCells = w;
            this.hullHCells = h;
            this.borderCells = border;
            this.cellMeters = cell;
            this.minDim = minDim;
            this.maxDim = maxDim;
            this.enforceAreas = enforceAreas;
//            this.enforceAreas = false;
            this.areaTolCells = areaTol;
            this.minAspect = minAsp;
            this.maxAspect = maxAsp;
            this.minContactCells = minContact;
//            this.minContactCells = 0;
            this.forbidUnwantedContacts = forbid;
            this.forceFillHull = forceFill;
            this.adjAtLeastOneSide = adjAtLeast;
            this.areaSqmByRoom = Collections.unmodifiableMap(areaMap);
            this.perRoomOptions = Collections.unmodifiableMap(roomOpts);
            this.G = g;
        }
    }

    // --- Hilfen & Validierung --------------------------------------------------

    private static int toCells(double meters, double cellMeters) {
        if (!(cellMeters > 0.0)) {
            throw new IllegalArgumentException("cellMeters muss > 0 sein");
        }
        double raw = meters / cellMeters;
        int cells = (int) Math.round(raw);
        if (cells <= 0) {
            throw new IllegalArgumentException("Boundary in Zellen ≤ 0 (Eingabe prüfen: meters=" + meters + ", cellMeters=" + cellMeters + ")");
        }
        return cells;
    }

    private static void validateGlobalOptions(CrossModalSpec.Spec spec) {
        var b = spec.boundary;
        var o = spec.options;

        if (b == null) throw new IllegalArgumentException("boundary fehlt");
        if (!(b.widthMeters > 0) || !(b.heightMeters > 0))
            throw new IllegalArgumentException("boundary width/height müssen > 0 sein");
        if (!(b.cellMeters > 0))
            throw new IllegalArgumentException("cellMeters muss > 0 sein");
        if (o == null) throw new IllegalArgumentException("options fehlen");

        if (o.minDimCells <= 0) throw new IllegalArgumentException("minDimCells muss > 0 sein");
        if (o.maxDimCells <= 0) throw new IllegalArgumentException("maxDimCells muss > 0 sein");
        if (o.minDimCells > o.maxDimCells) throw new IllegalArgumentException("minDimCells > maxDimCells");

        if (o.areaTolCells < 0) throw new IllegalArgumentException("areaTolCells < 0");
        if (o.minContactCells < 0) throw new IllegalArgumentException("minContactCells < 0");

        if (!(o.minAspect > 0.0) || !(o.maxAspect > 0.0))
            throw new IllegalArgumentException("minAspect/maxAspect müssen > 0 sein");
        if (o.minAspect > o.maxAspect)
            throw new IllegalArgumentException("minAspect > maxAspect");
    }

    private static void validateRoomOverride(String id, RoomLocalOptions r) {
        if (r == null) return;
        if (r.minDimCells != null && r.maxDimCells != null && r.minDimCells > r.maxDimCells) {
            throw new IllegalArgumentException("minDimCells > maxDimCells für Raum " + id);
        }
        if (r.minAspect != null && r.maxAspect != null && r.minAspect > r.maxAspect) {
            throw new IllegalArgumentException("minAspect > maxAspect für Raum " + id);
        }
        if (r.minThicknessCells != null && r.maxThicknessCells != null && r.minThicknessCells > r.maxThicknessCells) {
            throw new IllegalArgumentException("minThicknessCells > maxThicknessCells für Raum " + id);
        }
    }

    // --- Hauptmapping ----------------------------------------------------------

    public static MappedParams map(CrossModalSpec.Spec spec) {
        Objects.requireNonNull(spec, "spec");
        validateGlobalOptions(spec);

        var b = spec.boundary;
        var o = spec.options;

        final int W = toCells(b.widthMeters,  b.cellMeters);
        final int H = toCells(b.heightMeters, b.cellMeters);

        // Räume & Area-Map
        Map<String, Integer> areaMap = new LinkedHashMap<>();
        Map<String, RoomLocalOptions> perRoom = new HashMap<>();

        var G = new AdjacencyGraph();
        for (var r : spec.rooms) {
            Objects.requireNonNull(r, "room");
            if (r.id == null || r.id.isBlank()) throw new IllegalArgumentException("Raum mit leerer/null ID");

            G.addNode(r.id);

            // Fläche: in m² (Integer gerundet) speichern – Solver wandelt selbst in Zellen um
            if (r.areaSqm != null) {
                long rounded = Math.round(r.areaSqm);
                if (rounded < 0) throw new IllegalArgumentException("Negative areaSqm bei Raum " + r.id);
                areaMap.put(r.id, (int) rounded);
            }

            var local = new RoomLocalOptions(
                    r.minDimCells, r.maxDimCells,
                    r.minAspect, r.maxAspect,
                    r.shape,
                    r.minThicknessCells, r.maxThicknessCells
            );
            validateRoomOverride(r.id, local);
            perRoom.put(r.id, local);
        }

        // Adjazenzkanten aus Spec
        for (var e : spec.edges) {
            Objects.requireNonNull(e, "edge");
            if (e.rel == CrossModalSpec.Relation.ADJACENT) {
                // Kante nur hinzufügen, wenn beide IDs bekannt sind
                if (e.u == null || e.v == null) {
                    throw new IllegalArgumentException("Edge mit null-Endpoint");
                }
                G.addUndirectedEdge(e.u, e.v);
            }
        }

        return new MappedParams(
                o.solver,
                W, H, b.borderCells,
                b.cellMeters,
                o.minDimCells, o.maxDimCells,
                o.enforceAreas, o.areaTolCells,
                o.minAspect, o.maxAspect,
                o.minContactCells,
                o.forbidUnwantedContacts,
                o.forceFillHull,
                o.adjAtLeastOneSide,
                areaMap, perRoom, G
        );
    }
}
