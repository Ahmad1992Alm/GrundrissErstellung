package beispiel.fx;

import javafx.application.Application;
import com.google.ortools.Loader;
import com.google.ortools.sat.*;

import beispiel.FloorplanViewer;

import java.util.*;

public class DemoShowFX {
  public static void main(String[] args) {
    Loader.loadNativeLibraries();

    // ... hier baust/löst du dein Modell wie gehabt ...
    // Angenommen, du hast:
    // LayoutModelBuilder b = new LayoutModelBuilder(10.0, 8.0, rooms);
    // (evtl. Cut-outs: b.addForbiddenRectMeters(...);)
    // CpModel model = b.build();
    // CpSolver solver = new CpSolver();
    // solver.getParameters().setMaxTimeInSeconds(15.0);
    // CpSolverStatus st = solver.solve(model);

    // Für die Demo nehmen wir einfach die Werte aus deiner letzten Lösung:
    double hullWm = 10.0, hullHm = 8.0;

    List<PlacedRect> roomRects = new ArrayList<>();
    roomRects.add(new PlacedRect("Küche",  6.0, 2.0, 2.0, 3.0));
    roomRects.add(new PlacedRect("Essen",  6.0, 2.0, 2.0, 3.0));
    roomRects.add(new PlacedRect("Wohnen", 5.0, 5.0, 4.0, 2.0));
//    roomRects.add(new PlacedRect("schlaf", 2.0, 5.0, 4.0, 2.0));

    // Wenn du Cut-outs hast, trage sie hier ein (Meter):
    List<PlacedRect> cutouts = new ArrayList<>();
    // cutouts.add(new PlacedRect("Cutout1", 3.5, 6.0, 2.0, 1.0));
    // cutouts.add(new PlacedRect("Cutout2", 6.5, 0.0, 1.5, 1.0));

    // Daten in den Viewer
    FloorplanViewer.setData(hullWm, hullHm, roomRects, cutouts,
                            80.0,                     // px pro Meter (Zoom)
                            "Floorplan (JavaFX)");

    // JavaFX starten
    Application.launch(FloorplanViewer.class);
  }
}
