package pdf;

//import javafx.application.Application;
//import javafx.geometry.Insets;
//import javafx.scene.Scene;
//import javafx.scene.Group;
//import javafx.scene.control.Button;
//import javafx.scene.control.Label;
//import javafx.scene.layout.BorderPane;
//import javafx.scene.layout.HBox;
//import javafx.scene.paint.Color;
//import javafx.scene.shape.Rectangle;
//import javafx.scene.text.Text;
//import javafx.stage.Stage;
//
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Map;
//
//public class CrossModalJavaFxDemo extends Application {
//
//	private static final int CELL_PX = 60;
//	private static final int PAD = 20;
//	private List<FloorplanCrossModalSolver.Solution> solutions = new ArrayList<>();
//	private int index = 0;
//	private CrossModalMapper.MappedParams params;
//
//	public static void main(String[] args) {
//		launch(args);
//	}
//
//	@Override
//	public void start(Stage stage) {
//
//		CrossModalSpec.Spec spec = buildSpec();
//		this.params = CrossModalMapper.map(spec);
//		try {
//			if (params.solver == CrossModalSpec.ModelOptions.Solver.MIQP) {
//				this.solutions = MIQPFloorplanSolver.layout(params);
//			} else {
//				this.solutions = FloorplanCrossModalSolver.layout(params);
//			}
//		} catch (RuntimeException e) {
//			throw e; // zeige Solver-Status statt generischem Fallback
//		}
//
//		if (solutions.isEmpty())
//			throw new IllegalStateException("Keine Lösung gefunden.");
//
//		BorderPane root = new BorderPane();
//		Group canvas = new Group();
//		Label title = new Label();
//		title.setPadding(new Insets(10));
//		Label info = new Label();
//		info.setPadding(new Insets(10));
//		Button prev = new Button("◀ Vorherige");
//		Button next = new Button("Nächste ▶");
//		HBox controls = new HBox(10, prev, next, info);
//		controls.setPadding(new Insets(10));
//		prev.setOnAction(e -> {
//			index = (index - 1 + solutions.size()) % solutions.size();
//			redraw(canvas, title, info);
//		});
//		next.setOnAction(e -> {
//			index = (index + 1) % solutions.size();
//			redraw(canvas, title, info);
//		});
//
//		root.setTop(title);
//		root.setCenter(canvas);
//		root.setBottom(controls);
//		redraw(canvas, title, info);
//
//		int sceneW = PAD * 2 + params.hullWCells * CELL_PX + 1;
//		int sceneH = PAD * 2 + params.hullHCells * CELL_PX + 1 + 60;
//		Scene scene = new Scene(root, sceneW, sceneH, Color.WHITE);
//		scene.setOnKeyPressed(ke -> {
//			switch (ke.getCode()) {
//			case LEFT -> prev.fire();
//			case RIGHT -> next.fire();
//			}
//		});
//		stage.setTitle("Cross-Modal Floorplan Viewer – " + params.hullWCells + "×" + params.hullHCells
//				+ " Zellen, 1 Raster = " + params.cellMeters + " m");
//		stage.setScene(scene);
//		stage.show();
//	}
//
//	private void redraw(Group canvas, Label title, Label info) {
//		// Zeichenfläche zurücksetzen
//		canvas.getChildren().clear();
//
//		// Hülle zeichnen
//		Rectangle hull = new Rectangle(PAD, PAD, params.hullWCells * CELL_PX, params.hullHCells * CELL_PX);
//		hull.setFill(null);
//		hull.setStroke(Color.BLACK);
//		hull.setStrokeWidth(3.0);
//		canvas.getChildren().add(hull);
//
//		// Farbpalette
//		Color[] palette = new Color[] { Color.web("#8ecae6"), Color.web("#ffb703"), Color.web("#219ebc"),
//				Color.web("#fb8500"), Color.web("#90be6d"), Color.web("#f28482"), Color.web("#bdb2ff"),
//				Color.web("#ffd6a5") };
//
//		// aktuelle Lösung
//		FloorplanCrossModalSolver.Solution sol = solutions.get(index);
//
//		int k = 0;
//		for (Map.Entry<String, FloorplanCrossModalSolver.RoomPieces> e : sol.piecesByNode.entrySet()) {
//			String id = e.getKey();
//			FloorplanCrossModalSolver.RoomPieces rp = e.getValue(); // <- hier ist rp definiert
//
//			double X = PAD + rp.x0 * CELL_PX;
//			double Y = PAD + rp.y0 * CELL_PX;
//
//			// Prüfen, ob es in Wahrheit ein Rechteck ist (Balken decken sich vollständig)
//			boolean isRect = (rp.wH == rp.tX) && (rp.hV == rp.tY);
//
//			Color c = palette[k++ % palette.length];
//
//			if (isRect) {
//				// Ein einziges Rechteck (Bounding-Box entspricht wH×hV)
//				Rectangle rect = new Rectangle(X, Y, rp.wH * CELL_PX, rp.hV * CELL_PX);
//				rect.setFill(c.deriveColor(0, 1, 1, 0.75));
//				rect.setStroke(Color.BLACK);
//				rect.setStrokeWidth(1.5);
//
//				int areaCells = rp.wH * rp.hV;
//				double areaSqm = areaCells * params.cellMeters * params.cellMeters;
//				String labelText = id + " — " + String.format("%.1f m²", areaSqm) + " (" + rp.wH + "×" + rp.hV + ")";
//
//				Text label = new Text(X + 8, Y + 18, labelText);
//
//				canvas.getChildren().addAll(rect, label);
//			} else {
//				// L-Form: zwei Balken zeichnen (H und V), die sich an der Ecke überlappen
//				Rectangle rH = new Rectangle(X, Y, rp.wH * CELL_PX, rp.tY * CELL_PX);
//				Rectangle rV = new Rectangle(X, Y, rp.tX * CELL_PX, rp.hV * CELL_PX);
//
//				rH.setFill(c.deriveColor(0, 1, 1, 0.75));
//				rV.setFill(c.deriveColor(0, 1, 1, 0.75));
//				rH.setStroke(Color.BLACK);
//				rV.setStroke(Color.BLACK);
//				rH.setStrokeWidth(1.5);
//				rV.setStrokeWidth(1.5);
//
//				int areaCells = rp.wH * rp.tY + rp.tX * rp.hV - rp.tX * rp.tY; // Inklusions-Exklusionsfläche
//				double areaSqm = areaCells * params.cellMeters * params.cellMeters;
//				String labelText = id + " — " + String.format("%.1f m²", areaSqm) + "  H:(" + rp.wH + "×" + rp.tY
//						+ ") V:(" + rp.tX + "×" + rp.hV + ")";
//
//				Text label = new Text(X + 8, Y + 18, labelText);
//
//				canvas.getChildren().addAll(rH, rV, label);
//			}
//		}
//
//		// Kopf-/Infotext aktualisieren
//		title.setText("Lösung " + (index + 1) + "/" + solutions.size());
//		info.setText("maxX=" + sol.maxX + ", maxY=" + sol.maxY + "  |  Räume: " + sol.piecesByNode.size());
//	}
//
//	private static CrossModalSpec.Spec buildSpec() {
//		var spec = new CrossModalSpec.Spec().boundary(new CrossModalSpec.Boundary(7.5, 7.0, 0.5).setBorderCells(0));
////        var spec = new CrossModalSpec.Spec().boundary(new CrossModalSpec.Boundary(10.0, 10.0, 0.5).setBorderCells(0));
//
//		// Beispielräume (A–F) – A–D mit Flächen, E/F ohne
//		spec.addRoom(new CrossModalSpec.Room("A").area(15).shape(CrossModalSpec.Shape.RECT));
//		spec.addRoom(new CrossModalSpec.Room("B").area(15).shape(CrossModalSpec.Shape.RECT));
//		spec.addRoom(new CrossModalSpec.Room("C").area(15).shape(CrossModalSpec.Shape.RECT));
//		spec.addRoom(new CrossModalSpec.Room("D").area(15).shape(CrossModalSpec.Shape.RECT));
////		spec.addRoom(new CrossModalSpec.Room("E").shape(CrossModalSpec.Shape.RECT));
////		spec.addRoom(new CrossModalSpec.Room("F").shape(CrossModalSpec.Shape.RECT));
////		spec.addRoom(new CrossModalSpec.Room("g").shape(CrossModalSpec.Shape.RECT_OR_L));
//
//		// Adjazenzen (wie in deinem Beispiel)
//		spec.addEdge(new CrossModalSpec.Edge("A", "B"));
////		spec.addEdge(new CrossModalSpec.Edge("B", "D"));
//		spec.addEdge(new CrossModalSpec.Edge("D", "C"));
//////		spec.addEdge(new CrossModalSpec.Edge("C", "A"));
////		spec.addEdge(new CrossModalSpec.Edge("E", "A"));
//////		spec.addEdge(new CrossModalSpec.Edge("E", "C"));
////		spec.addEdge(new CrossModalSpec.Edge("E", "D"));
////		spec.addEdge(new CrossModalSpec.Edge("E", "F"));
//////		spec.addEdge(new CrossModalSpec.Edge("F", "E"));
////		spec.addEdge(new CrossModalSpec.Edge("F", "A"));
////		spec.addEdge(new CrossModalSpec.Edge("F", "B"));
////		spec.addEdge(new CrossModalSpec.Edge("F", "C"));
////		spec.addEdge(new CrossModalSpec.Edge("F", "D"));
//
//		// Lockeres Default-Setup, um Feasibility zu erhöhen
//		spec.options.forbidUnwantedContacts = true; // Nicht-Nachbarn dürfen nicht berühren
//		spec.options.forceFillHull = false; // Hülle nicht erzwingen (später aktivieren)
//		spec.options.adjAtLeastOneSide = true; // "mind. 1 Seite" statt "genau 1"
//		spec.options.areaTolCells = 8; // Toleranz etwas großzügiger
//		spec.options.solver = CrossModalSpec.ModelOptions.Solver.MIQP; // neuen MIQP-Solver verwenden
//
//		return spec;
//	}
//}


import javafx.application.Application;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.Group;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressIndicator;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.Rectangle;
import javafx.scene.text.Text;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class CrossModalJavaFxDemo extends Application {

    private static final int CELL_PX = 60;
    private static final int PAD = 20;

    private volatile List<FloorplanCrossModalSolver.Solution> solutions = new ArrayList<>();
    private int index = 0;
    private CrossModalMapper.MappedParams params;

    // UI
    private Label title = new Label();
    private Label info = new Label();
    private Label status = new Label();
    private Group canvas = new Group();
    private Button prev = new Button("◀ Vorherige");
    private Button next = new Button("Nächste ▶");
    private Button cancel = new Button("Abbrechen");
    private Button retryMilp = new Button("Als MILP neu versuchen");
    private ProgressIndicator spinner = new ProgressIndicator();

    private Task<List<FloorplanCrossModalSolver.Solution>> currentTask;

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    public void start(Stage stage) {
        // Spec bauen und Params mappen (das ist schnell)
        CrossModalSpec.Spec spec = buildSpec();
        this.params = CrossModalMapper.map(spec);

        // UI zusammenbauen (sofort anzeigen)
        BorderPane root = new BorderPane();
        title.setPadding(new Insets(10));
        info.setPadding(new Insets(10));
        status.setPadding(new Insets(10));
        spinner.setPrefSize(24, 24);
        spinner.setVisible(false);

        HBox controls = new HBox(10, prev, next, cancel, retryMilp, info, spinner);
        controls.setPadding(new Insets(10));

        prev.setOnAction(e -> { index = (index - 1 + solutions.size()) % solutions.size(); redraw(); });
        next.setOnAction(e -> { index = (index + 1) % solutions.size(); redraw(); });
        cancel.setOnAction(e -> { if (currentTask != null) currentTask.cancel(true); });
        retryMilp.setOnAction(e -> retryAsMilp());

        prev.setDisable(true);
        next.setDisable(true);
        cancel.setDisable(true);
        retryMilp.setDisable(true);

        root.setTop(title);
        root.setCenter(canvas);
        root.setBottom(controls);
        root.setLeft(status);

        int sceneW = PAD * 2 + params.hullWCells * CELL_PX + 1;
        int sceneH = PAD * 2 + params.hullHCells * CELL_PX + 1 + 60;
        Scene scene = new Scene(root, sceneW, sceneH, Color.WHITE);
        scene.setOnKeyPressed(ke -> { switch (ke.getCode()) {
            case LEFT -> prev.fire();
            case RIGHT -> next.fire();
            default -> {}
        }});
        stage.setTitle("Cross-Modal Floorplan Viewer – " + params.hullWCells + "×" + params.hullHCells
                + " Zellen, 1 Raster = " + params.cellMeters + " m");
        stage.setScene(scene);
        stage.show();

        // Solver im Hintergrund starten
        startSolveAsync();
    }

    private void startSolveAsync() {
        setBusy(true, "Löse Modell (" + params.G.nodes.size() + " Räume) …");
        currentTask = new Task<>() {
            @Override
            protected List<FloorplanCrossModalSolver.Solution> call() {
                // WICHTIG: Keine UI-Zugriffe hier!
                if (params.solver == CrossModalSpec.ModelOptions.Solver.MIQP) {
                    return MIQPFloorplanSolver.layout(params);
                } else {
                    return FloorplanCrossModalSolver.layout(params);
                }
            }
        };

        currentTask.setOnSucceeded(ev -> {
            solutions = currentTask.getValue();
            if (solutions == null || solutions.isEmpty()) {
                setBusy(false, "Keine Lösung gefunden.");
                prev.setDisable(true);
                next.setDisable(true);
                retryMilp.setDisable(false);
                return;
            }
            index = 0;
            setBusy(false, "Fertig.");
            prev.setDisable(solutions.size() <= 1);
            next.setDisable(solutions.size() <= 1);
            retryMilp.setDisable(false);
            redraw();
        });

        currentTask.setOnFailed(ev -> {
            Throwable ex = currentTask.getException();
            setBusy(false, "Fehler: " + (ex != null ? ex.getMessage() : "unbekannt"));
            retryMilp.setDisable(false);
        });

        currentTask.setOnCancelled(ev -> {
            setBusy(false, "Abgebrochen.");
            retryMilp.setDisable(false);
        });

        Thread t = new Thread(currentTask, "solver-thread");
        t.setDaemon(true);
        t.start();
    }

    private void retryAsMilp() {
        // „Als MILP neu versuchen“: schnellere Heuristik — wir schalten hier ein paar
        // Optionen in der Spec/Params um (linearere Zielsetzung etc.).
        status.setText("Starte Neuversuch als MILP (vereinfachte Zielfunktion) …");

        // Einfache Heuristik: CP-SAT verwenden (falls vorhanden), sonst MIQP behalten.
        // Alternativ könntest du in MIQPFloorplanSolver eine lineare Zielfunktion anbieten.
        // Hier demonstrativ auf CP_SAT umschalten:
        var spec = buildSpec();
        spec.options.solver = CrossModalSpec.ModelOptions.Solver.CP_SAT;
        this.params = CrossModalMapper.map(spec);

        startSolveAsync();
    }

    private void setBusy(boolean busy, String message) {
        status.setText(message);
        spinner.setVisible(busy);
        cancel.setDisable(!busy);
        prev.setDisable(busy || solutions.size() <= 1);
        next.setDisable(busy || solutions.size() <= 1);
        retryMilp.setDisable(busy);
    }

    private void redraw() {
        // Zeichenfläche zurücksetzen
        canvas.getChildren().clear();

        // Hülle zeichnen
        Rectangle hull = new Rectangle(PAD, PAD, params.hullWCells * CELL_PX, params.hullHCells * CELL_PX);
        hull.setFill(null);
        hull.setStroke(Color.BLACK);
        hull.setStrokeWidth(3.0);
        canvas.getChildren().add(hull);

        // Farbpalette
        Color[] palette = new Color[] {
                Color.web("#8ecae6"), Color.web("#ffb703"), Color.web("#219ebc"),
                Color.web("#fb8500"), Color.web("#90be6d"), Color.web("#f28482"),
                Color.web("#bdb2ff"), Color.web("#ffd6a5")
        };

        // aktuelle Lösung
        FloorplanCrossModalSolver.Solution sol = solutions.get(index);

        int k = 0;
        for (Map.Entry<String, FloorplanCrossModalSolver.RoomPieces> e : sol.piecesByNode.entrySet()) {
            String id = e.getKey();
            FloorplanCrossModalSolver.RoomPieces rp = e.getValue();

            double X = PAD + rp.x0 * CELL_PX;
            double Y = PAD + rp.y0 * CELL_PX;

            boolean isRect = (rp.wH == rp.tX) && (rp.hV == rp.tY);
            Color c = palette[k++ % palette.length];

            if (isRect) {
                Rectangle rect = new Rectangle(X, Y, rp.wH * CELL_PX, rp.hV * CELL_PX);
                rect.setFill(c.deriveColor(0, 1, 1, 0.75));
                rect.setStroke(Color.BLACK);
                rect.setStrokeWidth(1.5);

                int areaCells = rp.wH * rp.hV;
                double areaSqm = areaCells * params.cellMeters * params.cellMeters;
                String labelText = id + " — " + String.format("%.1f m²", areaSqm) + " (" + rp.wH + "×" + rp.hV + ")";
                Text label = new Text(X + 8, Y + 18, labelText);

                canvas.getChildren().addAll(rect, label);
            } else {
                Rectangle rH = new Rectangle(X, Y, rp.wH * CELL_PX, rp.tY * CELL_PX);
                Rectangle rV = new Rectangle(X, Y, rp.tX * CELL_PX, rp.hV * CELL_PX);

                rH.setFill(c.deriveColor(0, 1, 1, 0.75));
                rV.setFill(c.deriveColor(0, 1, 1, 0.75));
                rH.setStroke(Color.BLACK);
                rV.setStroke(Color.BLACK);
                rH.setStrokeWidth(1.5);
                rV.setStrokeWidth(1.5);

                int areaCells = rp.wH * rp.tY + rp.tX * rp.hV - rp.tX * rp.tY; // Inklusions-Exklusionsfläche
                double areaSqm = areaCells * params.cellMeters * params.cellMeters;
                String labelText = id + " — " + String.format("%.1f m²", areaSqm)
                        + "  H:(" + rp.wH + "×" + rp.tY + ") V:(" + rp.tX + "×" + rp.hV + ")";
                Text label = new Text(X + 8, Y + 18, labelText);

                canvas.getChildren().addAll(rH, rV, label);
            }
        }

        // Kopf-/Infotext aktualisieren
        title.setText("Lösung " + (index + 1) + "/" + solutions.size());
        info.setText("maxX=" + sol.maxX + ", maxY=" + sol.maxY + "  |  Räume: " + sol.piecesByNode.size());
    }

    private static CrossModalSpec.Spec buildSpec() {
        var spec = new CrossModalSpec.Spec()
                .boundary(new CrossModalSpec.Boundary(6.5, 6.5, 0.5).setBorderCells(0));

        // Beispielräume (A–F) – A–D mit Flächen, E/F ohne
        spec.addRoom(new CrossModalSpec.Room("A").area(6).shape(CrossModalSpec.Shape.RECT));
        spec.addRoom(new CrossModalSpec.Room("B").area(6).shape(CrossModalSpec.Shape.RECT));
        spec.addRoom(new CrossModalSpec.Room("C").area(6).shape(CrossModalSpec.Shape.RECT));
        spec.addRoom(new CrossModalSpec.Room("D").area(6).shape(CrossModalSpec.Shape.RECT));
        spec.addRoom(new CrossModalSpec.Room("E").area(4).shape(CrossModalSpec.Shape.RECT));
        spec.addRoom(new CrossModalSpec.Room("F").area(3).shape(CrossModalSpec.Shape.RECT));
        spec.addRoom(new CrossModalSpec.Room("g").area(4).shape(CrossModalSpec.Shape.RECT));

        // Adjazenzen
        spec.addEdge(new CrossModalSpec.Edge("A", "B"));
        spec.addEdge(new CrossModalSpec.Edge("A", "F"));
        spec.addEdge(new CrossModalSpec.Edge("A", "E"));
        spec.addEdge(new CrossModalSpec.Edge("B", "A"));
        spec.addEdge(new CrossModalSpec.Edge("B", "F"));
        spec.addEdge(new CrossModalSpec.Edge("C", "D"));
        spec.addEdge(new CrossModalSpec.Edge("C", "F"));
        spec.addEdge(new CrossModalSpec.Edge("D", "C"));
        spec.addEdge(new CrossModalSpec.Edge("D", "F"));
        spec.addEdge(new CrossModalSpec.Edge("D", "E"));
        spec.addEdge(new CrossModalSpec.Edge("E", "A"));
        spec.addEdge(new CrossModalSpec.Edge("E", "D"));
        spec.addEdge(new CrossModalSpec.Edge("E", "F"));
        spec.addEdge(new CrossModalSpec.Edge("F", "A"));
        spec.addEdge(new CrossModalSpec.Edge("F", "B"));
        spec.addEdge(new CrossModalSpec.Edge("F", "C"));
        spec.addEdge(new CrossModalSpec.Edge("F", "D"));
        spec.addEdge(new CrossModalSpec.Edge("g", "B"));
        spec.addEdge(new CrossModalSpec.Edge("g", "F"));
        spec.addEdge(new CrossModalSpec.Edge("g", "C"));

        // Optionen (etwas „entschärft“)
        spec.options.forbidUnwantedContacts = true; // Nicht-Nachbarn dürfen sich nicht berühren
        spec.options.forceFillHull = true;         // Hülle nicht erzwingen
//        spec.options.forceFillHull = false;         // Hülle nicht erzwingen
        spec.options.adjAtLeastOneSide = false;      // mind. 1 Seite
        spec.options.areaTolCells = 12;              // großzügigere Toleranz
        spec.options.solver = CrossModalSpec.ModelOptions.Solver.CP_SAT; // MIQP

        return spec;
    }
}
