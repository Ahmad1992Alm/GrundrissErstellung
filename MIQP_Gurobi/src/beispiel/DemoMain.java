package beispiel;

import com.google.ortools.Loader;
import com.google.ortools.sat.*;

import java.util.*;

public class DemoMain {
  public static void main(String[] args) {
    Loader.loadNativeLibraries();

    Room k = new Room();
    k.name="Kueche"; k.wStarM=3.0; k.hStarM=3.0; k.wMinM=2.0; k.wMaxM=6.0; k.hMinM=2.5; k.hMaxM=5.0;

    Room e = new Room();
    e.name="Essen"; e.wStarM=3.5; e.hStarM=3.0; e.wMinM=2.0; e.wMaxM=6.0; e.hMinM=2.5; e.hMaxM=5.0;

    Room w = new Room();
    w.name="Wohnen"; w.wStarM=4.5; w.hStarM=3.5; w.wMinM=2.0; w.wMaxM=6.0; w.hMinM=2.5; w.hMaxM=5.0;
    w.options.mustTouch = WallSide.TOP;
    w.options.aspectRatioMin = 1.2;
    w.options.aspectRatioMax = 2.0;

    Adjacency ke = new Adjacency(); ke.other="Essen";  ke.weight=10; ke.requireDoor=true;
    Adjacency kw = new Adjacency(); kw.other="Wohnen"; kw.weight=3;
    Adjacency ew = new Adjacency(); ew.other="Wohnen"; ew.weight=4;
    k.neighbors.add(ke); k.neighbors.add(kw); e.neighbors.add(ew);

    List<Room> rooms = Arrays.asList(k,e,w);

    LayoutModelBuilder b = new LayoutModelBuilder(10.0, 8.0, rooms);

    // Beispiel-Cutouts (Meter anpassen):
    b.addForbiddenRectMeters(3.5, 6.0, 2.0, 1.0);
    b.addForbiddenRectMeters(6.5, 0.0, 1.5, 1.0);

    CpModel model = b.build();

    CpSolver solver = new CpSolver();
    solver.getParameters().setMaxTimeInSeconds(15.0);
    CpSolverStatus st = solver.solve(model);
    System.out.println("Status: " + st + "  Obj: " + solver.objectiveValue());
    if (st == CpSolverStatus.OPTIMAL || st == CpSolverStatus.FEASIBLE) {
      for (Room r : rooms) {
        double xi = solver.value(b.getX().get(r.name))/LayoutModelBuilder.SCALE;
        double yi = solver.value(b.getY().get(r.name))/LayoutModelBuilder.SCALE;
        double wi = solver.value(b.getW().get(r.name))/LayoutModelBuilder.SCALE;
        double hi = solver.value(b.getH().get(r.name))/LayoutModelBuilder.SCALE;
        System.out.printf("%-8s x=%.2f  y=%.2f  w=%.2f  h=%.2f%n", r.name, xi, yi, wi, hi);
      }
    }
  }
}
