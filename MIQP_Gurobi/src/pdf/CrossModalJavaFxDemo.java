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
import java.util.HashMap;
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
//////*************************////////////////////***********************/////////////////****************
	// Topologie zur aktuellen Lösung (wird vor jedem redraw() aktualisiert)
	private Map<Long, FloorplanCrossModalSolver.Ori> topology = new HashMap<>();

	// Grabbable Overlays für Kanten
	private final List<javafx.scene.shape.Rectangle> seamOverlays = new ArrayList<>();

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

//		prev.setOnAction(e -> {
//			index = (index - 1 + solutions.size()) % solutions.size();
//			redraw();
//		});
//		next.setOnAction(e -> {
//			index = (index + 1) % solutions.size();
//			redraw();
//		});
		/// ****************/////********
		prev.setOnAction(e -> {
			if (solutions == null || solutions.isEmpty())
				return;
			index = (index - 1 + solutions.size()) % solutions.size();
			topology = FloorplanCrossModalSolver.inferTopology(solutions.get(index), params.G.nodes);
			redraw();
		});
		next.setOnAction(e -> {
			if (solutions == null || solutions.isEmpty())
				return;
			index = (index + 1) % solutions.size();
			topology = FloorplanCrossModalSolver.inferTopology(solutions.get(index), params.G.nodes);
			redraw();
		});

		cancel.setOnAction(e -> {
			if (currentTask != null)
				currentTask.cancel(true);
		});
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
		scene.setOnKeyPressed(ke -> {
			switch (ke.getCode()) {
			case LEFT -> prev.fire();
			case RIGHT -> next.fire();
			default -> {
			}
			}
		});
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
//			solutions = currentTask.getValue();
			//// *********************//////////////******************************
			//// *********************//////////////******************************
			//// *********************//////////////******************************
			solutions = new ArrayList<>(currentTask.getValue()); // -> mutierbar

			
			if (solutions == null || solutions.isEmpty()) {
				setBusy(false, "Keine Lösung gefunden.");
				prev.setDisable(true);
				next.setDisable(true);
				retryMilp.setDisable(false);
				return;
			}
			index = 0;
			//// *********************//////////////******************************
			//// *********************//////////////******************************
			//// *********************//////////////******************************
			this.topology = FloorplanCrossModalSolver.inferTopology(solutions.get(0), params.G.nodes);

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
			this.topology = new HashMap<>();

		});

		currentTask.setOnCancelled(ev -> {
			setBusy(false, "Abgebrochen.");
			retryMilp.setDisable(false);
			this.topology = new HashMap<>();

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
		// Alternativ könntest du in MIQPFloorplanSolver eine lineare Zielfunktion
		// anbieten.
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
		Color[] palette = new Color[] { Color.web("#8ecae6"), Color.web("#ffb703"), Color.web("#219ebc"),
				Color.web("#fb8500"), Color.web("#90be6d"), Color.web("#f28482"), Color.web("#bdb2ff"),
				Color.web("#ffd6a5") };

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
				String labelText = id + " — " + String.format("%.1f m²", areaSqm) + "  H:(" + rp.wH + "×" + rp.tY
						+ ") V:(" + rp.tX + "×" + rp.hV + ")";
				Text label = new Text(X + 8, Y + 18, labelText);

				canvas.getChildren().addAll(rH, rV, label);
			}
		}

		// Kopf-/Infotext aktualisieren
		title.setText("Lösung " + (index + 1) + "/" + solutions.size());
		info.setText("maxX=" + sol.maxX + ", maxY=" + sol.maxY + "  |  Räume: " + sol.piecesByNode.size());

		//// *****************////////////////***************//////////////////
		//// *****************////////////////***************//////////////////
		//// *****************////////////////***************//////////////////
		rebuildSeamOverlays(); 

	}

	private static CrossModalSpec.Spec buildSpec() {
		var spec = new CrossModalSpec.Spec().boundary(new CrossModalSpec.Boundary(6.5, 6.5, 0.5).setBorderCells(0));

		// Beispielräume (A–F) – A–D mit Flächen, E/F ohne
		spec.addRoom(new CrossModalSpec.Room("A").area(4).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("B").area(4).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("C").area(4).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("D").area(4).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("E").area(4).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("F").area(3).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("g").area(2).shape(CrossModalSpec.Shape.RECT));
		spec.addRoom(new CrossModalSpec.Room("l").area(4).shape(CrossModalSpec.Shape.L));

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
		spec.addEdge(new CrossModalSpec.Edge("g", "l"));
		spec.addEdge(new CrossModalSpec.Edge("C", "l"));
		spec.addEdge(new CrossModalSpec.Edge("D", "l"));
		spec.addEdge(new CrossModalSpec.Edge("E", "l"));

		// Optionen (etwas „entschärft“)
		spec.options.forbidUnwantedContacts = true; // Nicht-Nachbarn dürfen sich nicht berühren
//		spec.options.forceFillHull = true; // Hülle nicht erzwingen
        spec.options.forceFillHull = true;         // Hülle nicht erzwingen
		spec.options.adjAtLeastOneSide = true; // mind. 1 Seite
		spec.options.areaTolCells = 5; // großzügigere Toleranz
		spec.options.solver = CrossModalSpec.ModelOptions.Solver.CP_SAT; // MIQP

		return spec;
	}

	/////// ********************************************/////////////////////******************************************///////////
	/////// ********************************************/////////////////////******************************************///////////
	/////// ********************************************/////////////////////******************************************///////////
	/////// ********************************************/////////////////////******************************************///////////
	/////// ********************************************/////////////////////******************************************///////////
	/////// ********************************************/////////////////////******************************************///////////

	private static long pairKey(int i, int j) {
		return (((long) i) << 32) | (j & 0xffffffffL);
	}

	private static int clamp(int v, int lo, int hi) {
		return Math.max(lo, Math.min(hi, v));
	}

	// Szene → Gitterzelle (X/Y), robust gegen Layout-Shifts:
	private int sceneToCellX(double sceneX, double sceneY) {
		var p = canvas.sceneToLocal(sceneX, sceneY);
		return clamp((int) Math.round((p.getX() - PAD) / CELL_PX), 0, params.hullWCells);
	}

	private int sceneToCellY(double sceneX, double sceneY) {
		var p = canvas.sceneToLocal(sceneX, sceneY);
		return clamp((int) Math.round((p.getY() - PAD) / CELL_PX), 0, params.hullHCells);
	}

//	private void rebuildSeamOverlays() {
//	    // 1) Alte Overlays aus der Szene entfernen und Liste leeren
//	    for (var r : seamOverlays) {
//	        canvas.getChildren().remove(r);
//	    }
//	    seamOverlays.clear();
//
//	    // 2) Guard
//	    if (solutions == null || solutions.isEmpty()) return;
//	    final var sol0 = solutions.get(index);
//
//	    // 3) Für jedes Nachbarpaar ein Handle auf die gemeinsame Kante
//	    for (int i = 0; i < params.G.nodes.size(); i++) {
//	        final String idI = params.G.nodes.get(i);
//	        final var ri = sol0.piecesByNode.get(idI);
//	        final int xi = ri.x0, yi = ri.y0, wi = ri.wH, hi = ri.hV;
//
//	        for (int j = i + 1; j < params.G.nodes.size(); j++) {
//	            final int ii = i;
//	            final int jj = j;
//
//	            final String idJ = params.G.nodes.get(j);
//	            if (!params.G.areNeighbors(idI, idJ)) continue;
//
//	            final var rj = sol0.piecesByNode.get(idJ);
//	            final int xj = rj.x0, yj = rj.y0, wj = rj.wH, hj = rj.hV;
//
//	            final var ori = topology.get(pairKey(ii, jj));
//	            if (ori == null) continue;
//
//	            switch (ori) {
//	                case L -> { // i links, j rechts: vertikale Kante bei S = xj
//	                    int S = xj;
//	                    int yTop = Math.max(yi, yj);
//	                    int yBot = Math.min(yi + hi, yj + hj);
//	                    int hCells = Math.max(0, yBot - yTop);
//	                    if (hCells <= 0) break;
//
//	                    double rx = PAD + S * CELL_PX - 3;      // 6px breit
//	                    double ry = PAD + yTop * CELL_PX;
//	                    double rw = 6;
//	                    double rh = hCells * CELL_PX;
//
//	                    var handle = new Rectangle(rx, ry, rw, rh);
//	                    handle.setFill(Color.color(1, 0, 0, 0.18)); // << sichtbar zum Testen
//	                    handle.setStroke(Color.color(1, 0, 0, 0.35));
//	                    handle.setPickOnBounds(true);
//
//	                    handle.setOnMouseReleased(ev -> {
//	                        int newX = sceneToCellX(ev.getSceneX(), ev.getSceneY());
//	                        var edited = FloorplanCrossModalSolver.localEditSolve(
//	                                params, sol0, topology,
//	                                new FloorplanCrossModalSolver.MoveVerticalSeam(ii, jj, newX)
//	                        );
//	                        if (edited != null) {
//	                        	List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
//	                        	copy.set(index, edited);
//	                        	solutions = copy; // Austausch der (nun mutierbaren) Liste
//
//	                        	topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
//	                        	redraw();
//
//	                        }
//	                    });
//
//	                    canvas.getChildren().add(handle);
//	                    seamOverlays.add(handle);
//	                }
//	                case R -> { // i rechts, j links: vertikale Kante bei S = xi
//	                    int S = xi;
//	                    int yTop = Math.max(yi, yj);
//	                    int yBot = Math.min(yi + hi, yj + hj);
//	                    int hCells = Math.max(0, yBot - yTop);
//	                    if (hCells <= 0) break;
//
//	                    double rx = PAD + S * CELL_PX - 3;
//	                    double ry = PAD + yTop * CELL_PX;
//	                    double rw = 6;
//	                    double rh = hCells * CELL_PX;
//
//	                    var handle = new Rectangle(rx, ry, rw, rh);
//	                    handle.setFill(Color.color(1, 0, 0, 0.18));
//	                    handle.setStroke(Color.color(1, 0, 0, 0.35));
//	                    handle.setPickOnBounds(true);
//
//	                    handle.setOnMouseReleased(ev -> {
//	                        int newX = sceneToCellX(ev.getSceneX(), ev.getSceneY());
//	                        // R: (x_j + w_j = x_i) -> Reihenfolge umdrehen
//	                        var edited = FloorplanCrossModalSolver.localEditSolve(
//	                                params, sol0, topology,
//	                                new FloorplanCrossModalSolver.MoveVerticalSeam(jj, ii, newX)
//	                        );
//	                        if (edited != null) {
//	                        	List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
//	                        	copy.set(index, edited);
//	                        	solutions = copy; // Austausch der (nun mutierbaren) Liste
//
//	                        	topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
//	                        	redraw();
//
//	                        }
//	                    });
//
//	                    canvas.getChildren().add(handle);
//	                    seamOverlays.add(handle);
//	                }
//	                case T -> { // i oben, j unten: horizontale Kante bei T = yj
//	                    int T = yj;
//	                    int xL = Math.max(xi, xj);
//	                    int xR = Math.min(xi + wi, xj + wj);
//	                    int wCells = Math.max(0, xR - xL);
//	                    if (wCells <= 0) break;
//
//	                    double rx = PAD + xL * CELL_PX;
//	                    double ry = PAD + T * CELL_PX - 3;
//	                    double rw = wCells * CELL_PX;
//	                    double rh = 6;
//
//	                    var handle = new Rectangle(rx, ry, rw, rh);
//	                    handle.setFill(Color.color(0, 0, 1, 0.18));
//	                    handle.setStroke(Color.color(0, 0, 1, 0.35));
//	                    handle.setPickOnBounds(true);
//
//	                    handle.setOnMouseReleased(ev -> {
//	                        int newY = sceneToCellY(ev.getSceneX(), ev.getSceneY());
//	                        var edited = FloorplanCrossModalSolver.localEditSolve(
//	                                params, sol0, topology,
//	                                new FloorplanCrossModalSolver.MoveHorizontalSeam(ii, jj, newY)
//	                        );
//	                        if (edited != null) {
//	                        	List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
//	                        	copy.set(index, edited);
//	                        	solutions = copy; // Austausch der (nun mutierbaren) Liste
//
//	                        	topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
//	                        	redraw();
//
//	                        }
//	                    });
//
//	                    canvas.getChildren().add(handle);
//	                    seamOverlays.add(handle);
//	                }
//	                case B -> { // i unten, j oben: horizontale Kante bei T = yi
//	                    int T = yi;
//	                    int xL = Math.max(xi, xj);
//	                    int xR = Math.min(xi + wi, xj + wj);
//	                    int wCells = Math.max(0, xR - xL);
//	                    if (wCells <= 0) break;
//
//	                    double rx = PAD + xL * CELL_PX;
//	                    double ry = PAD + T * CELL_PX - 3;
//	                    double rw = wCells * CELL_PX;
//	                    double rh = 6;
//
//	                    var handle = new Rectangle(rx, ry, rw, rh);
//	                    handle.setFill(Color.color(0, 0, 1, 0.18));
//	                    handle.setStroke(Color.color(0, 0, 1, 0.35));
//	                    handle.setPickOnBounds(true);
//
//	                    handle.setOnMouseReleased(ev -> {
//	                        int newY = sceneToCellY(ev.getSceneX(), ev.getSceneY());
//	                        // B: (y_j + h_j = y_i) -> Reihenfolge umdrehen
//	                        var edited = FloorplanCrossModalSolver.localEditSolve(
//	                                params, sol0, topology,
//	                                new FloorplanCrossModalSolver.MoveHorizontalSeam(jj, ii, newY)
//	                        );
//	                        if (edited != null) {
//	                        	List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
//	                        	copy.set(index, edited);
//	                        	solutions = copy; // Austausch der (nun mutierbaren) Liste
//
//	                        	topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
//	                        	redraw();
//
//	                        }
//	                    });
//
//	                    canvas.getChildren().add(handle);
//	                    seamOverlays.add(handle);
//	                }
//	            }
//	        }
//	    }
//	}
	private void rebuildSeamOverlays() {
	    // 1) Alte Overlays weg
	    for (var r : seamOverlays) canvas.getChildren().remove(r);
	    seamOverlays.clear();

	    if (solutions == null || solutions.isEmpty()) return;
	    final var sol0 = solutions.get(index);

	    // ---- A) Gemeinsame (Nachbarschafts-)Kanten: L/R/T/B ----
	    for (int i = 0; i < params.G.nodes.size(); i++) {
	        final String idI = params.G.nodes.get(i);
	        final var ri = sol0.piecesByNode.get(idI);
	        final int xi = ri.x0, yi = ri.y0, wi = ri.wH, hi = ri.hV;

	        for (int j = i + 1; j < params.G.nodes.size(); j++) {
	            final int ii = i, jj = j;

	            final String idJ = params.G.nodes.get(j);
	            if (!params.G.areNeighbors(idI, idJ)) continue;

	            final var rj = sol0.piecesByNode.get(idJ);
	            final int xj = rj.x0, yj = rj.y0, wj = rj.wH, hj = rj.hV;

	            final var ori = topology.get(pairKey(ii, jj));
	            if (ori == null) continue;

	            switch (ori) {
	                case L -> { // vertikale gemeinsame Kante bei S = xj
	                    int S = xj;
	                    int yTop = Math.max(yi, yj);
	                    int yBot = Math.min(yi + hi, yj + hj);
	                    int hCells = Math.max(0, yBot - yTop);
	                    if (hCells <= 0) break;

	                    double rx = PAD + S * CELL_PX - 3, ry = PAD + yTop * CELL_PX;
	                    Rectangle handle = new Rectangle(rx, ry, 6, hCells * CELL_PX);
	                    handle.setFill(Color.color(1, 0, 0, 0.15));
	                    handle.setStroke(Color.color(1, 0, 0, 0.35));
	                    handle.setPickOnBounds(true);

	                    handle.setOnMousePressed(ev -> {
	                        handle.setFill(Color.color(1, 0, 0, 0.35)); // highlight
	                    });
	                    handle.setOnMouseReleased(ev -> {
	                        int newX = sceneToCellX(ev.getSceneX(), ev.getSceneY());
	                        var edited = FloorplanCrossModalSolver.localEditSolve(
	                                params, sol0, topology,
	                                new FloorplanCrossModalSolver.MoveVerticalSeam(ii, jj, newX)
	                        );
	                        handle.setFill(Color.color(1, 0, 0, 0.15)); // zurück
	                        if (edited != null) {
	                            // mutable replace
	                            List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                            copy.set(index, edited);
	                            solutions = copy;
	                            topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                            redraw();
	                        }
	                    });

	                    canvas.getChildren().add(handle);
	                    seamOverlays.add(handle);
	                }
	                case R -> { // vertikale gemeinsame Kante bei S = xi
	                    int S = xi;
	                    int yTop = Math.max(yi, yj);
	                    int yBot = Math.min(yi + hi, yj + hj);
	                    int hCells = Math.max(0, yBot - yTop);
	                    if (hCells <= 0) break;

	                    double rx = PAD + S * CELL_PX - 3, ry = PAD + yTop * CELL_PX;
	                    Rectangle handle = new Rectangle(rx, ry, 6, hCells * CELL_PX);
	                    handle.setFill(Color.color(1, 0, 0, 0.15));
	                    handle.setStroke(Color.color(1, 0, 0, 0.35));
	                    handle.setPickOnBounds(true);

	                    handle.setOnMousePressed(ev -> handle.setFill(Color.color(1, 0, 0, 0.35)));
	                    handle.setOnMouseReleased(ev -> {
	                        int newX = sceneToCellX(ev.getSceneX(), ev.getSceneY());
	                        var edited = FloorplanCrossModalSolver.localEditSolve(
	                                params, sol0, topology,
	                                new FloorplanCrossModalSolver.MoveVerticalSeam(jj, ii, newX) // R ⇒ (j,i)
	                        );
	                        handle.setFill(Color.color(1, 0, 0, 0.15));
	                        if (edited != null) {
	                            List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                            copy.set(index, edited);
	                            solutions = copy;
	                            topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                            redraw();
	                        }
	                    });

	                    canvas.getChildren().add(handle);
	                    seamOverlays.add(handle);
	                }
	                case T -> { // horizontale gemeinsame Kante bei T = yj
	                    int T = yj;
	                    int xL = Math.max(xi, xj);
	                    int xR = Math.min(xi + wi, xj + wj);
	                    int wCells = Math.max(0, xR - xL);
	                    if (wCells <= 0) break;

	                    double rx = PAD + xL * CELL_PX, ry = PAD + T * CELL_PX - 3;
	                    Rectangle handle = new Rectangle(rx, ry, wCells * CELL_PX, 6);
	                    handle.setFill(Color.color(1, 0, 0, 0.15));
	                    handle.setStroke(Color.color(1, 0, 0, 0.35));
	                    handle.setPickOnBounds(true);

	                    handle.setOnMousePressed(ev -> handle.setFill(Color.color(1, 0, 0, 0.35)));
	                    handle.setOnMouseReleased(ev -> {
	                        int newY = sceneToCellY(ev.getSceneX(), ev.getSceneY());
	                        var edited = FloorplanCrossModalSolver.localEditSolve(
	                                params, sol0, topology,
	                                new FloorplanCrossModalSolver.MoveHorizontalSeam(ii, jj, newY)
	                        );
	                        handle.setFill(Color.color(1, 0, 0, 0.15));
	                        if (edited != null) {
	                            List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                            copy.set(index, edited);
	                            solutions = copy;
	                            topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                            redraw();
	                        }
	                    });

	                    canvas.getChildren().add(handle);
	                    seamOverlays.add(handle);
	                }
	                case B -> { // horizontale gemeinsame Kante bei T = yi
	                    int T = yi;
	                    int xL = Math.max(xi, xj);
	                    int xR = Math.min(xi + wi, xj + wj);
	                    int wCells = Math.max(0, xR - xL);
	                    if (wCells <= 0) break;

	                    double rx = PAD + xL * CELL_PX, ry = PAD + T * CELL_PX - 3;
	                    Rectangle handle = new Rectangle(rx, ry, wCells * CELL_PX, 6);
	                    handle.setFill(Color.color(1, 0, 0, 0.15));
	                    handle.setStroke(Color.color(1, 0, 0, 0.35));
	                    handle.setPickOnBounds(true);

	                    handle.setOnMousePressed(ev -> handle.setFill(Color.color(1, 0, 0, 0.35)));
	                    handle.setOnMouseReleased(ev -> {
	                        int newY = sceneToCellY(ev.getSceneX(), ev.getSceneY());
	                        var edited = FloorplanCrossModalSolver.localEditSolve(
	                                params, sol0, topology,
	                                new FloorplanCrossModalSolver.MoveHorizontalSeam(jj, ii, newY) // B ⇒ (j,i)
	                        );
	                        handle.setFill(Color.color(1, 0, 0, 0.15));
	                        if (edited != null) {
	                            List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                            copy.set(index, edited);
	                            solutions = copy;
	                            topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                            redraw();
	                        }
	                    });

	                    canvas.getChildren().add(handle);
	                    seamOverlays.add(handle);
	                }
	            }
	        }
	    }

	    // ---- B) Einzelkanten pro Raum: Left / Right / Top / Bottom ----
	    for (int i = 0; i < params.G.nodes.size(); i++) {
	        final int ii = i;
	        final String id = params.G.nodes.get(i);
	        final var r = sol0.piecesByNode.get(id);
	        final int xi = r.x0, yi = r.y0, wi = r.wH, hi = r.hV;

	        // LEFT edge at X = xi
	        {
	            double rx = PAD + xi * CELL_PX - 2, ry = PAD + yi * CELL_PX;
	            Rectangle h = new Rectangle(rx, ry, 4, hi * CELL_PX);
	            h.setFill(Color.color(0, 0, 0, 0.08));
	            h.setStroke(null);
	            h.setPickOnBounds(true);
	            h.setOnMousePressed(ev -> h.setFill(Color.color(1, 0, 0, 0.35)));
	            h.setOnMouseReleased(ev -> {
	                int newX = sceneToCellX(ev.getSceneX(), ev.getSceneY());
	                var edited = FloorplanCrossModalSolver.localEditSolve(
	                        params, sol0, topology,
	                        new FloorplanCrossModalSolver.MoveLeftEdge(ii, newX)
	                );
	                h.setFill(Color.color(0, 0, 0, 0.08));
	                if (edited != null) {
	                    List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                    copy.set(index, edited);
	                    solutions = copy;
	                    topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                    redraw();
	                }
	            });
	            canvas.getChildren().add(h);
	            seamOverlays.add(h);
	        }
	        // RIGHT edge at X = xi + wi
	        {
	            double rx = PAD + (xi + wi) * CELL_PX - 2, ry = PAD + yi * CELL_PX;
	            Rectangle h = new Rectangle(rx, ry, 4, hi * CELL_PX);
	            h.setFill(Color.color(0, 0, 0, 0.08));
	            h.setPickOnBounds(true);
	            h.setOnMousePressed(ev -> h.setFill(Color.color(1, 0, 0, 0.35)));
	            h.setOnMouseReleased(ev -> {
	                int newX = sceneToCellX(ev.getSceneX(), ev.getSceneY());
	                var edited = FloorplanCrossModalSolver.localEditSolve(
	                        params, sol0, topology,
	                        new FloorplanCrossModalSolver.MoveRightEdge(ii, newX)
	                );
	                h.setFill(Color.color(0, 0, 0, 0.08));
	                if (edited != null) {
	                    List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                    copy.set(index, edited);
	                    solutions = copy;
	                    topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                    redraw();
	                }
	            });
	            canvas.getChildren().add(h);
	            seamOverlays.add(h);
	        }
	        // TOP edge at Y = yi
	        {
	            double rx = PAD + xi * CELL_PX, ry = PAD + yi * CELL_PX - 2;
	            Rectangle h = new Rectangle(rx, ry, wi * CELL_PX, 4);
	            h.setFill(Color.color(0, 0, 0, 0.08));
	            h.setPickOnBounds(true);
	            h.setOnMousePressed(ev -> h.setFill(Color.color(1, 0, 0, 0.35)));
	            h.setOnMouseReleased(ev -> {
	                int newY = sceneToCellY(ev.getSceneX(), ev.getSceneY());
	                var edited = FloorplanCrossModalSolver.localEditSolve(
	                        params, sol0, topology,
	                        new FloorplanCrossModalSolver.MoveTopEdge(ii, newY)
	                );
	                h.setFill(Color.color(0, 0, 0, 0.08));
	                if (edited != null) {
	                    List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                    copy.set(index, edited);
	                    solutions = copy;
	                    topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                    redraw();
	                }
	            });
	            canvas.getChildren().add(h);
	            seamOverlays.add(h);
	        }
	        // BOTTOM edge at Y = yi + hi
	        {
	            double rx = PAD + xi * CELL_PX, ry = PAD + (yi + hi) * CELL_PX - 2;
	            Rectangle h = new Rectangle(rx, ry, wi * CELL_PX, 4);
	            h.setFill(Color.color(0, 0, 0, 0.08));
	            h.setPickOnBounds(true);
	            h.setOnMousePressed(ev -> h.setFill(Color.color(1, 0, 0, 0.35)));
	            h.setOnMouseReleased(ev -> {
	                int newY = sceneToCellY(ev.getSceneX(), ev.getSceneY());
	                var edited = FloorplanCrossModalSolver.localEditSolve(
	                        params, sol0, topology,
	                        new FloorplanCrossModalSolver.MoveBottomEdge(ii, newY)
	                );
	                h.setFill(Color.color(0, 0, 0, 0.08));
	                if (edited != null) {
	                    List<FloorplanCrossModalSolver.Solution> copy = new ArrayList<>(solutions);
	                    copy.set(index, edited);
	                    solutions = copy;
	                    topology = FloorplanCrossModalSolver.inferTopology(edited, params.G.nodes);
	                    redraw();
	                }
	            });
	            canvas.getChildren().add(h);
	            seamOverlays.add(h);
	        }
	    }
	}


}
