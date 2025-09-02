package pdf;

import java.util.*;

/**
 * Mappt eine {@link CrossModalSpec.Spec} in ein Parameter-Bündel, das der
 * Solver konsumiert.
 */
public final class CrossModalMapper {
	// ---- Graph für den Solver -------------------------------------------------
	public static class AdjacencyGraph {
		public final List<String> nodes = new ArrayList<>();
		public final Map<String, Set<String>> nbrs = new HashMap<>();

		public void addNode(String id) {
			if (!nbrs.containsKey(id)) {
				nbrs.put(id, new HashSet<>());
				nodes.add(id);
			}
		}

		public void addUndirectedEdge(String a, String b) {
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
                public final Map<String, Integer> areaSqmByRoom;
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
			this.areaTolCells = areaTol;
			this.minAspect = minAsp;
			this.maxAspect = maxAsp;
			this.minContactCells = minContact;
			this.forbidUnwantedContacts = forbid;
			this.forceFillHull = forceFill;
			this.adjAtLeastOneSide = adjAtLeast;
			this.areaSqmByRoom = Collections.unmodifiableMap(areaMap);
			this.perRoomOptions = Collections.unmodifiableMap(roomOpts);
			this.G = g;
		}
	}

	private static int toCells(double meters, double cellMeters) {
		return (int) Math.round(meters / cellMeters);
	}

	public static MappedParams map(CrossModalSpec.Spec spec) {
		Objects.requireNonNull(spec, "spec");
		var b = spec.boundary;
		var o = spec.options;
		int W = toCells(b.widthMeters, b.cellMeters);
		int H = toCells(b.heightMeters, b.cellMeters);

		Map<String, Integer> areaMap = new LinkedHashMap<>();
		Map<String, RoomLocalOptions> perRoom = new HashMap<>();

		var G = new AdjacencyGraph();
		for (var r : spec.rooms) {
			G.addNode(r.id);
			if (r.areaSqm != null)
				areaMap.put(r.id, (int) Math.round(r.areaSqm));
			perRoom.put(r.id, new RoomLocalOptions(r.minDimCells, r.maxDimCells, r.minAspect, r.maxAspect, r.shape,
					r.minThicknessCells, r.maxThicknessCells));
		}
		for (var e : spec.edges) {
			if (e.rel == CrossModalSpec.Relation.ADJACENT)
				G.addUndirectedEdge(e.u, e.v);
		}

                return new MappedParams(o.solver, W, H, b.borderCells, b.cellMeters, o.minDimCells, o.maxDimCells,
                                o.enforceAreas, o.areaTolCells, o.minAspect, o.maxAspect, o.minContactCells,
                                o.forbidUnwantedContacts, o.forceFillHull, o.adjAtLeastOneSide, areaMap, perRoom, G);
        }

}