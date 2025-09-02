package beispiel.fx;

import com.google.ortools.Loader;
import com.google.ortools.sat.CpModel;
import com.google.ortools.sat.CpSolver;
import com.google.ortools.sat.CpSolverStatus;

import beispiel.Adjacency;
import beispiel.FloorplanViewer;
import beispiel.LayoutModelBuilder;
import beispiel.Room;
import javafx.application.Application;

import java.util.*;

public class GraphToFloorplanMain {
  public static void main(String[] args) {
    Loader.loadNativeLibraries();

    // ------- dein Graph: A neben B, B neben C (aber A nicht zwingend neben C)
    Map<String, List<String>> graph = new HashMap<>();
    graph.put("Kueche", Arrays.asList("Essen"));
    graph.put("Essen", Arrays.asList("Kueche", "Wohnen"));
    graph.put("Wohnen", Arrays.asList("Kueche","Essen")); // Achte: kein vollständiges Klicken nötig

    // ------- Builder + Räume aus Graph
    // 3) Builder + Räume aus Graph erzeugen
    double hullW = 10.0, hullH = 8.0; // Hülle in Metern
    LayoutModelBuilder tmp = new LayoutModelBuilder(hullW, hullH, new ArrayList<>());
    List<Room> rooms = roomsFromGraph(graph,
        /*w* ,h* */ 3.0, 3.0,
        /*wMin,wMax*/ 2.0, 6.0,
        /*hMin,hMax*/ 2.5, 5.0,
        tmp);
    // neuen Builder mit den echten Räumen erstellen (tmp nur, um requireAdjacency() zu setzen)
    LayoutModelBuilder b = new LayoutModelBuilder(hullW, hullH, rooms);
    // ggf. Cut-outs:
    // b.addForbiddenRectMeters(6.5, 0.0, 1.5, 1.0);

    // 4) Modell bauen & lösen
    CpModel model = b.build();
    CpSolver solver = new CpSolver();
    solver.getParameters().setMaxTimeInSeconds(20.0);
    CpSolverStatus st = solver.solve(model);
    System.out.println("Status: " + st + "  Obj: " + solver.objectiveValue());

    // 5) Lösung einsammeln → PlacedRect in Metern
    List<PlacedRect> placed = new ArrayList<>();
    double S = LayoutModelBuilder.SCALE; // z.B. 10 => 0.1 m Raster
    for (Room r : rooms) {
      double xi = solver.value(b.getX().get(r.name)) / S;
      double yi = solver.value(b.getY().get(r.name)) / S;
      double wi = solver.value(b.getW().get(r.name)) / S;
      double hi = solver.value(b.getH().get(r.name)) / S;
      System.out.printf("%-8s x=%.2f y=%.2f w=%.2f h=%.2f%n", r.name, xi, yi, wi, hi);
      placed.add(new PlacedRect(r.name, xi, yi, wi, hi));
    }

    // 6) JavaFX-Viewer füttern & starten (einmal pro JVM!)
    FloorplanViewer.setData(
        hullW, hullH,
        placed,
        /* cutouts: */ Collections.emptyList(),  // oder aus b.getForbidden() umrechnen
        /* Zoom px/m: */ 80.0,
        "Grundriss aus Graph"
    );
    Application.launch(FloorplanViewer.class);
  }

  // Trick: wir müssen builder.requireAdjacency(...) während des Room-Baus aufrufen.
  // Deshalb erzeugen wir den Builder zuerst separat:
//Räume + Pflicht-Adjazenzen aus der Map aufbauen
 static List<Room> roomsFromGraph(
     Map<String, List<String>> graph,
     double wStar, double hStar,
     double wMin, double wMax,
     double hMin, double hMax,
     LayoutModelBuilder builderForAdj
 ) {
   Map<String, Room> byName = new HashMap<>();
   for (String name : graph.keySet()) {
     Room r = new Room();
     r.name = name;
     r.wStarM = wStar; r.hStarM = hStar;
     r.wMinM = wMin;  r.wMaxM = wMax;
     r.hMinM = hMin;  r.hMaxM = hMax;
     byName.put(name, r);
   }
   for (var e : graph.entrySet()) {
     String a = e.getKey();
     for (String b : e.getValue()) {
       if (!byName.containsKey(b)) continue;
       // Zielfunktionsgewicht + Türpflicht
       Adjacency adj = new Adjacency();
       adj.other = b;
       adj.weight = 20;
       adj.requireDoor = true;
       byName.get(a).neighbors.add(adj);
       // Kantenkontakt als Pflicht (ungerichtet)
       if (builderForAdj != null) builderForAdj.requireAdjacency(a, b);
     }
   }
   return new ArrayList<>(byName.values());
 }
}
