package beispiel0309;

import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FloorplanMIQP {
	public static void main(String[] args) {
		Loader.loadNativeLibraries();
		MPSolver solver = MPSolver.createSolver("SCIP"); // Für MIQP
		if (solver == null) {
			System.out.println("Solver nicht erstellt.");
			return;
		}

		// Beispiel: Räume (aus Paper: Rechtecke mit x,y,w,d)
		int numRooms = 5; // Anzahl Räume
		double domainWidth = 10.0; // Gebäudebreite
		double domainDepth = 10.0; // Gebäudetiefe
		double M = domainWidth + domainDepth; // Big-M für Constraints

		// Variablen für jedes Rechteck (pro Raum)
		MPVariable[] x = new MPVariable[numRooms];
		MPVariable[] y = new MPVariable[numRooms];
		MPVariable[] w = new MPVariable[numRooms];
		MPVariable[] d = new MPVariable[numRooms];
		for (int i = 0; i < numRooms; i++) {
			x[i] = solver.makeNumVar(0, domainWidth, "x" + i);
			y[i] = solver.makeNumVar(0, domainDepth, "y" + i);
			w[i] = solver.makeNumVar(1, domainWidth, "w" + i); // Min/Max Größe
			d[i] = solver.makeNumVar(1, domainDepth, "d" + i);
		}

		// Inside-Constraint (Paper Eq. 1)
		for (int i = 0; i < numRooms; i++) {
//			solver.makeConstraint(0, domainWidth).setCoefficient(x[i] + w[i], 1); // xi + wi <= w
//			// Ähnlich für y, d
			
			  MPConstraint inX = solver.makeConstraint(-MPSolver.infinity(), domainWidth);
			  inX.setCoefficient(x[i], 1.0);
			  inX.setCoefficient(w[i], 1.0);
		}

		// Non-Overlap-Constraint (Paper Eq. 2) mit binären Variablen
		for (int i = 0; i < numRooms; i++) {
			for (int j = i + 1; j < numRooms; j++) {
				MPVariable[] sigma = new MPVariable[4]; // front, back, left, right
				for (int k = 0; k < 4; k++)
					sigma[k] = solver.makeBoolVar("sigma" + i + j + k);
				// Constraints: xi - wj >= xj - M*(1-sigmaR), usw.
				// Sum sigma >=1
			}
		}

		// Adjacency aus Graph (deine Map)
		Map<String, List<String>> adjacencyMap = new HashMap<>(); // Dein Graph
		// Fülle Map...
		// Für jedes Paar: Füge Cadj-Constraints (Paper Eq. 10)

		// Objective: Min Ecover + Esize (Paper Eq. 3/6/12) - quadratic
		MPObjective objective = solver.objective();
		// Für quadratic: objective.setQuadraticCoefficient(w[i], d[i], -1); // Für
		// Area-Term (wi*di)
		objective.setMinimization();

		// Löse
		MPSolver.ResultStatus result = solver.solve();
		if (result == MPSolver.ResultStatus.OPTIMAL) {
			// Extrahiere Lösung: Drucke x,y,w,d
		}
	}
}
