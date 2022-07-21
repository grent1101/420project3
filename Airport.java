package cmsc420_s22;

// YOU SHOULD NOT MODIFY THIS FILE

/**
 * Airport object. An airport is represented by its: - 3-letter IATA airport
 * code (e.g., "LAX") - city (e.g., "Los Angeles") - coordinates (x,y), given as
 * longitude (-180...180) and latitude (-90...+90)
 */
public class Airport implements LabeledPoint2D, Comparable<Airport> {

	// Private data (The x,y coordinates are provided by the Point2D class)

	private final String code; // IATA code
	private final String city; // city
	private final Point2D location; // coordinates

	/**
	 * Basic constructor.
	 * 
	 * @param code IATA code (String)
	 * @param city airport's city (String)
	 * @param x    x-coordinate [longitude] (double)
	 * @param y    y-coordinate [latitude] (double)
	 */
	public Airport(String code, String city, double x, double y) {
		this.code = code;
		this.city = city;
		this.location = new Point2D(x, y);
	}

	/**
	 * Minimal constructor. Constructs from just coordinates. This is used to refer
	 * to an arbitrary location in space.
	 * 
	 * @param x x-coordinate [longitude] (double)
	 * @param y y-coordinate [latitude] (double)
	 */
	public Airport(double x, double y) {
		this.code = "???";
		this.city = "Not-yet-initialized";
		this.location = new Point2D(x, y);
	}

	// Standard functions - Getters and toString
	public String getCode() {
		return code;
	}

	public String getCity() {
		return city;
	}

	public String toString() {
		return code + ":" + location;
	}

	// Required by LabeledPoint2D interface

	public String getLabel() {
		return code;
	}

	public double getX() {
		return location.getX();
	}

	public double getY() {
		return location.getY();
	}

	public double get(int i) {
		return location.get(i);
	}

	public Point2D getPoint2D() {
		return location;
	}

	public boolean equals(Object o) {
		if (this == o)
			return true;
		if (o instanceof Airport) {
			Airport ap = (Airport) o;
			if (code.equals(ap.code) && city.equals(ap.city) && location.equals(ap.location))
				return true;
		}
		return false;
	}

	/**
	 * Comparison function uses code for sorting.
	 * 
	 * @param ap Other airport
	 * @return Result of comparing codes
	 */
	public int compareTo(Airport ap) {
		return code.compareTo(ap.code);
	}

	/**
	 * Generate string according to desired format.
	 * 
	 * @param format Formatting code
	 */
	public String getString(String format) {
		if (format.compareTo("default") == 0 || format.compareTo("") == 0) {
			return toString();
		} else if (format.compareTo("code-only") == 0) {
			return code;
		} else if (format.compareTo("full") == 0) {
			return code + ", " + city + ", (" + location.getX() + "," + location.getY() + ")";
		} else if (format.compareTo("attributes") == 0) {
			return city + ", (" + location.getX() + "," + location.getY() + ")";
		} else
			return "Invalid format!";
	}
	
	/**
	 * Returns a hash code for this Airport.
	 */
	@Override
	public int hashCode() {
	    int hash = 7;
	    hash = 31 * hash + (code == null ? 0 : code.hashCode());
	    hash = 31 * hash + (city == null ? 0 : city.hashCode());
	    hash = 31 * hash + (location == null ? 0 : location.hashCode());
	    return hash;
	}

}
