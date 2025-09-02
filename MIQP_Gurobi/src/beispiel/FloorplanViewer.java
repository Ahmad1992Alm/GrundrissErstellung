package beispiel;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.Group;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;

import beispiel.fx.PlacedRect;

public class FloorplanViewer extends Application {

  // ---- statische Übergabe (einfach & robust) ----
  private static double hullWidthM;   // Hüllenbreite in Metern
  private static double hullHeightM;  // Hüllenhöhe in Metern
  private static final List<PlacedRect> rooms = new ArrayList<>();
  private static final List<PlacedRect> cutouts = new ArrayList<>();
  private static String title = "Floorplan";

  // Pixel pro Meter (Zoom)
  private static double PPM = 80.0; // 1 m = 80 px (anpassbar)

  // öffentlicher Setup-Entry:
  public static void setData(double hullWm, double hullHm,
                             List<PlacedRect> roomRects,
                             List<PlacedRect> cutoutRects,
                             double pixelsPerMeter,
                             String windowTitle) {
    hullWidthM = hullWm;
    hullHeightM = hullHm;
    PPM = pixelsPerMeter > 0 ? pixelsPerMeter : PPM;

    rooms.clear();
    if (roomRects != null) rooms.addAll(roomRects);

    cutouts.clear();
    if (cutoutRects != null) cutouts.addAll(cutoutRects);

    if (windowTitle != null) title = windowTitle;
  }

  @Override
  public void start(Stage stage) {
    // Canvasgröße mit Rand
    double padding = 40;
    double widthPx  = hullWidthM  * PPM + 2*padding;
    double heightPx = hullHeightM * PPM + 2*padding;

    Canvas canvas = new Canvas(widthPx, heightPx);
    GraphicsContext g = canvas.getGraphicsContext2D();

    // Hintergrund
    g.setFill(Color.web("#f7f7fb"));
    g.fillRect(0, 0, widthPx, heightPx);

    // Koordinatensystem: Modell hat Ursprung unten-links.
    // Hilfsfunktionen:
    java.util.function.DoubleUnaryOperator X = (xm) -> padding + xm * PPM;
    java.util.function.DoubleUnaryOperator Y = (ym) -> padding + (hullHeightM - ym) * PPM;
    java.util.function.DoubleUnaryOperator H = (hm) -> hm * PPM;
    java.util.function.DoubleUnaryOperator W = (wm) -> wm * PPM;

    // Grid
    drawGrid(g, padding, widthPx, heightPx, hullWidthM, hullHeightM);

    // Hülle
    g.setStroke(Color.BLACK);
    g.setLineWidth(2);
    g.strokeRect(X.applyAsDouble(0), Y.applyAsDouble(hullHeightM), W.applyAsDouble(hullWidthM), H.applyAsDouble(hullHeightM));

    // Cut-outs (verbotene Zonen)
    if (!cutouts.isEmpty()) {
      g.setGlobalAlpha(0.35);
      g.setFill(Color.GRAY);
      for (PlacedRect c : cutouts) {
        double x = X.applyAsDouble(c.x);
        double y = Y.applyAsDouble(c.y + c.h);
        double w = W.applyAsDouble(c.w);
        double h = H.applyAsDouble(c.h);
        g.fillRect(x, y, w, h);
        // Rahmen
        g.setGlobalAlpha(1.0);
        g.setStroke(Color.DARKGRAY);
        g.setLineDashes(6);
        g.strokeRect(x, y, w, h);
        g.setLineDashes(null);
        g.setGlobalAlpha(0.35);
      }
      g.setGlobalAlpha(1.0);
    }

    // Räume
    Font labelFont = Font.font("Arial", 14);
    int colorIdx = 0;
    Color[] palette = new Color[]{
      Color.web("#5B8FF9"), Color.web("#61DDAA"), Color.web("#65789B"),
      Color.web("#F6BD16"), Color.web("#7262fd"), Color.web("#78D3F8"),
      Color.web("#9661BC"), Color.web("#F6903D"), Color.web("#1E9493")
    };

    for (PlacedRect r : rooms) {
      Color fill = palette[colorIdx % palette.length];
      colorIdx++;

      double x = X.applyAsDouble(r.x);
      double y = Y.applyAsDouble(r.y + r.h);
      double w = W.applyAsDouble(r.w);
      double h = H.applyAsDouble(r.h);

      // Füllung
      g.setFill(fill.deriveColor(0,1,1,0.25));
      g.fillRect(x, y, w, h);
      // Rahmen
      g.setStroke(fill.darker());
      g.setLineWidth(2);
      g.strokeRect(x, y, w, h);

      // Label
      g.setFill(Color.BLACK);
      g.setFont(labelFont);
      String label = r.name + String.format("  (%.1fm × %.1fm)", r.w, r.h);
      g.fillText(label, x + 6, y + 18);

      // Maße (optional)
      g.setFill(Color.DARKGRAY);
      g.fillText(String.format("x=%.1f,y=%.1f", r.x, r.y), x + 6, y + 36);
    }

    Group root = new Group(canvas);
    Scene scene = new Scene(root, widthPx, heightPx, Color.WHITE);
    stage.setTitle(title);
    stage.setScene(scene);
    stage.show();
  }

  private void drawGrid(GraphicsContext g, double pad, double widthPx, double heightPx,
                        double hullWm, double hullHm) {
    g.setStroke(Color.rgb(0,0,0,0.07));
    g.setLineWidth(1);

    // vertikale Linien pro 1 m
    for (int xm = 0; xm <= (int)Math.round(hullWm); xm++) {
      double x = pad + xm * PPM;
      g.strokeLine(x, pad, x, pad + hullHm*PPM);
    }
    // horizontale Linien pro 1 m
    for (int ym = 0; ym <= (int)Math.round(hullHm); ym++) {
      double y = pad + ym * PPM;
      g.strokeLine(pad, y, pad + hullWm*PPM, y);
    }

    // Achsen-Beschriftung
    g.setFill(Color.GRAY);
    g.fillText("x (m)", pad + hullWm*PPM - 30, pad + hullHm*PPM + 20);
    g.fillText("y (m)", pad - 25, pad + 12);
  }
}
