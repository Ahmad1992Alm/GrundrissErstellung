package pdf;

import java.util.*;
import org.ojalgo.optimisation.ExpressionsBasedModel;
import org.ojalgo.optimisation.Optimisation;
import org.ojalgo.optimisation.Expression;
import org.ojalgo.optimisation.Variable;
/**
 * Einfache MIQP-Variante des Grundriss-Solvers. Verwendet das freie
 * ojAlgo-Optimierungspaket und modelliert jeden Raum als Rechteck mit
 * ganzzahligen Koordinaten (x,y,w,h). L-förmige Räume werden aktuell als
 * Rechteck angenähert.
 */
public final class MIQPFloorplanSolver {

    private MIQPFloorplanSolver() {
    }

    public static List<FloorplanCrossModalSolver.Solution> layout(CrossModalMapper.MappedParams p) {
        var G = p.G;
        int n = G.nodes.size();
        if (n == 0) throw new IllegalArgumentException("Graph hat keine Knoten");

        ExpressionsBasedModel model = new ExpressionsBasedModel();
        Variable[] x = new Variable[n];
        Variable[] y = new Variable[n];
        Variable[] w = new Variable[n];
        Variable[] h = new Variable[n];

        for (int i = 0; i < n; i++) {
            var id = G.nodes.get(i);
            var ro = p.perRoomOptions.get(id);
            int minDim = ro != null && ro.minDimCells != null ? ro.minDimCells : p.minDim;
            int maxDim = ro != null && ro.maxDimCells != null ? ro.maxDimCells : p.maxDim;
            x[i] = Variable.make("x_" + i).integer(true).lower(0).upper(p.hullWCells);
            y[i] = Variable.make("y_" + i).integer(true).lower(0).upper(p.hullHCells);
            w[i] = Variable.make("w_" + i).integer(true).lower(minDim).upper(maxDim);
            h[i] = Variable.make("h_" + i).integer(true).lower(minDim).upper(maxDim);
            model.addVariable(x[i]);
            model.addVariable(y[i]);
            model.addVariable(w[i]);
            model.addVariable(h[i]);

            // Hüllbegrenzung
            model.addExpression("bndX_" + i).upper(p.hullWCells).set(x[i], 1).set(w[i], 1);
            model.addExpression("bndY_" + i).upper(p.hullHCells).set(y[i], 1).set(h[i], 1);

            // Fläche exakt einhalten (Quadratische Gleichung)
            Integer area = p.areaSqmByRoom.get(id);
            if (area != null && p.enforceAreas) {
                Expression areaExpr = model.addExpression("area_" + i).level(area);
                areaExpr.set(w[i], h[i]);
            }
        }

        int M = p.hullWCells + p.hullHCells; // Big-M

        // Adjazenz und Nicht-Nachbarn
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                String idI = G.nodes.get(i);
                String idJ = G.nodes.get(j);
                boolean neighbors = G.areNeighbors(idI, idJ);
                if (neighbors) {
                    Variable L = Variable.binary("L_" + i + "_" + j);
                    Variable R = Variable.binary("R_" + i + "_" + j);
                    Variable T = Variable.binary("T_" + i + "_" + j);
                    Variable B = Variable.binary("B_" + i + "_" + j);
                    model.addVariable(L);
                    model.addVariable(R);
                    model.addVariable(T);
                    model.addVariable(B);
                    model.addExpression("ori_" + i + "_" + j).level(1).set(L, 1).set(R, 1).set(T, 1).set(B, 1);

                    // L: i rechts an j links
                    model.addExpression("Lx1_" + i + "_" + j).upper(M)
                            .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(L, M);
                    model.addExpression("Lx2_" + i + "_" + j).upper(M)
                            .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(L, M);
                    model.addExpression("Ly1_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(L, M);
                    model.addExpression("Ly2_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(L, M);

                    // R: i links an j rechts
                    model.addExpression("Rx1_" + i + "_" + j).upper(M)
                            .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(R, M);
                    model.addExpression("Rx2_" + i + "_" + j).upper(M)
                            .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(R, M);
                    model.addExpression("Ry1_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(R, M);
                    model.addExpression("Ry2_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(R, M);

                    // T: i unten an j oben
                    model.addExpression("Tx1_" + i + "_" + j).upper(M)
                            .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(T, M);
                    model.addExpression("Tx2_" + i + "_" + j).upper(M)
                            .set(y[j], 1).set(y[i], -1).set(h[i], -1).set(T, M);
                    model.addExpression("Tx3_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(T, M);
                    model.addExpression("Tx4_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(T, M);

                    // B: i oben an j unten
                    model.addExpression("Bx1_" + i + "_" + j).upper(M)
                            .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(B, M);
                    model.addExpression("Bx2_" + i + "_" + j).upper(M)
                            .set(y[i], 1).set(y[j], -1).set(h[j], -1).set(B, M);
                    model.addExpression("Bx3_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(x[i], 1).set(x[j], -1).set(w[j], -1).set(B, M);
                    model.addExpression("Bx4_" + i + "_" + j).upper(M - p.minContactCells)
                            .set(x[j], 1).set(x[i], -1).set(w[i], -1).set(B, M);
                } else if (p.forbidUnwantedContacts) {
                    Variable left = Variable.binary("left_" + i + "_" + j);
                    Variable right = Variable.binary("right_" + i + "_" + j);
                    Variable above = Variable.binary("above_" + i + "_" + j);
                    Variable below = Variable.binary("below_" + i + "_" + j);
                    model.addVariable(left);
                    model.addVariable(right);
                    model.addVariable(above);
                    model.addVariable(below);
                    model.addExpression("sep_" + i + "_" + j).level(1).set(left, 1).set(right, 1).set(above, 1).set(below, 1);
                    model.addExpression("sepL_" + i + "_" + j).upper(M)
                            .set(x[i], 1).set(w[i], 1).set(x[j], -1).set(left, M);
                    model.addExpression("sepR_" + i + "_" + j).upper(M)
                            .set(x[j], 1).set(w[j], 1).set(x[i], -1).set(right, M);
                    model.addExpression("sepA_" + i + "_" + j).upper(M)
                            .set(y[i], 1).set(h[i], 1).set(y[j], -1).set(above, M);
                    model.addExpression("sepB_" + i + "_" + j).upper(M)
                            .set(y[j], 1).set(h[j], 1).set(y[i], -1).set(below, M);
                }
            }
        }

        // Bounding Box und Zielfunktion
        Variable maxX = Variable.make("maxX").integer(true).lower(0).upper(p.hullWCells);
        Variable maxY = Variable.make("maxY").integer(true).lower(0).upper(p.hullHCells);
        model.addVariable(maxX);
        model.addVariable(maxY);
        for (int i = 0; i < n; i++) {
            model.addExpression("maxX_" + i).lower(0).set(maxX, 1).set(x[i], -1).set(w[i], -1);
            model.addExpression("maxY_" + i).lower(0).set(maxY, 1).set(y[i], -1).set(h[i], -1);
        }

        org.ojalgo.optimisation.Expression obj = model.addExpression("OBJ").weight(1.0);
        obj.set(maxX, maxX);
        obj.set(maxY, maxY);
        for (int i = 0; i < n; i++) {
            Expression asp = model.addExpression("asp_" + i).weight(0.1);
            asp.set(w[i], w[i]);
            asp.set(h[i], h[i]);
            asp.set(w[i], h[i], -2);
        }

        Optimisation.Result res = model.minimise();
        if (!res.getState().isFeasible()) {
            return List.of();
        }

        FloorplanCrossModalSolver.Solution sol = new FloorplanCrossModalSolver.Solution();
        sol.maxX = (int) Math.round(res.get(maxX));
        sol.maxY = (int) Math.round(res.get(maxY));
        for (int i = 0; i < n; i++) {
            FloorplanCrossModalSolver.RoomPieces rp = new FloorplanCrossModalSolver.RoomPieces();
            rp.x0 = (int) Math.round(res.get(x[i]));
            rp.y0 = (int) Math.round(res.get(y[i]));
            rp.wH = rp.tX = (int) Math.round(res.get(w[i]));
            rp.tY = rp.hV = (int) Math.round(res.get(h[i]));
            rp.wBB = rp.wH;
            rp.hBB = rp.tY;
            sol.piecesByNode.put(G.nodes.get(i), rp);
        }

        return List.of(sol);
    }
}
