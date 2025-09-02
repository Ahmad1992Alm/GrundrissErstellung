package beispiel;

import java.util.ArrayList;
import java.util.List;

public class Room {
	public String name;
	// Zielmaße (weich), Meter
	public double wStarM, hStarM;
	// harte Bounds (Meter)
	public double wMinM, wMaxM, hMinM, hMaxM;
	public List<Adjacency> neighbors = new ArrayList<>();
	public RoomOptions options = new RoomOptions();
}
