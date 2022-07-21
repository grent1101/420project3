package cmsc420_s22;

// YOU SHOULD NOT MODIFY THIS FILE. If you want to add additional functionality,
// create a new class/file (e.g., MyPoint2D.java) and modify that.

/**
 * A 2-dimensional point. We represent a point as a 2-element array of type
 * double.
 *
 * The functionality is a superset of Java's builtin Point2D, and some of the
 * code here has been taken from the Java OpenJDK source code.
 */

public class Point2D {
	public final static int DIM = 2; // spatial dimension
	private double[] coord; // coordinates

	/**
	 * Default constructor. Generates (0,0)
	 * 
	 */
	public Point2D() {
		coord = new double[DIM];
		for (int i = 0; i < DIM; i++)
			coord[i] = 0;
	}

	/**
	 * Construct from coordinates.
	 * 
	 * @param coord The array of coordinates.
	 */
	public Point2D(double x, double y) {
		coord = new double[DIM];
		coord[0] = x;
		coord[1] = y;
	}

	/**
	 * Construct from a 2-element coordinate array.
	 * 
	 * @param coord The array of coordinates.
	 */
	public Point2D(double[] coord) {
		assert (coord.length == DIM);
		this.coord = new double[DIM];
		for (int i = 0; i < DIM; i++)
			this.coord[i] = coord[i];
	}

	/**
	 * Copy constructor.
	 * 
	 * @param pt The point to copy.
	 */
	public Point2D(Point2D pt) {
		coord = new double[DIM];
		for (int i = 0; i < DIM; i++)
			coord[i] = pt.get(i);
	}

	/**
	 * Get the dimension.
	 * 
	 * @return The dimension
	 */
	public static int getDimension() {
		return DIM;
	}

	/**
	 * Get the x-coordinate.
	 * 
	 * @return The x-coordinate of the point.
	 */
	public double getX() {
		return coord[0];
	}

	/**
	 * Get the y-coordinate.
	 * 
	 * @return The y-coordinate of the point.
	 */
	public double getY() {
		return coord[1];
	}

	/**
	 * Get the i-th coordinate (i=0 for x, i=1 for y).
	 * 
	 * @param i The coordinate index 0 or 1.
	 * @return The i-th coordinate of the point.
	 */
	public double get(int i) {
		return coord[i];
	}

	/**
	 * Set both coordinates
	 * 
	 * @param x The new x coordinate.
	 * @param y The new y coordinate.
	 */
	public void setLocation(double x, double y) {
		coord[0] = x;
		coord[1] = y;
	}

	/**
	 * Set coordinates from a given point.
	 * 
	 * @param pt The point to copy from.
	 */
	public void setLocation(Point2D pt) {
		for (int i = 0; i < DIM; i++) {
			coord[i] = pt.coord[i];
		}
	}

	/**
	 * Set the i-th coordinate (i=0 for x, i=1 for y).
	 * 
	 * @param i The coordinate index 0 or 1.
	 * @param x The new coordinate value.
	 */
	public void set(int i, double x) {
		coord[i] = x;
	}

	/**
	 * Squared Euclidean distance between two points.
	 * 
	 * @param x1 The x coordinate of the first point.
	 * @param y1 The y coordinate of the first point.
	 * @param x2 The x coordinate of the second point.
	 * @param y2 The y coordinate of the second point.
	 * @return The squared distance between the two points.
	 */
	public static double distanceSq(double x1, double y1, double x2, double y2) {
		return Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2);
	}

	/**
	 * Euclidean distance between two points.
	 * 
	 * @param x1 The x coordinate of the first point.
	 * @param y1 The y coordinate of the first point.
	 * @param x2 The x coordinate of the second point.
	 * @param y2 The y coordinate of the second point.
	 * @return The distance between the two points.
	 */
	public static double distance(double x1, double y1, double x2, double y2) {
		return Math.sqrt(distanceSq(x1, y1, x2, y2));
	}

	/**
	 * Squared Euclidean distance to another point.
	 * 
	 * @param px The x coordinate of the other point.
	 * @param py The y coordinate of the other point.
	 * @return The squared distance between the two points.
	 */
	public double distanceSq(double px, double py) {
		return Math.pow(coord[0] - px, 2) + Math.pow(coord[1] - py, 2);
	}

	/**
	 * Euclidean distance to another point.
	 * 
	 * @param px The x coordinate of the other point.
	 * @param py The y coordinate of the other point.
	 * @return The squared distance between the two points.
	 */
	public double distance(double px, double py) {
		return Math.sqrt(distanceSq(px, py));
	}

	/**
	 * Squared Euclidean distance to another point.
	 * 
	 * @param pt The other point.
	 * @return The squared distance between the two points.
	 */
	public double distanceSq(Point2D pt) {
		if (pt == null)
			return Double.POSITIVE_INFINITY;
		else
			return Math.pow(coord[0] - pt.coord[0], 2) + Math.pow(coord[1] - pt.coord[1], 2);
	}

	/**
	 * Determines if two points are equal.
	 * 
	 * @param pt The other point
	 * @return True if equal
	 */
	public boolean equals(Object obj) {
		if (obj instanceof Point2D) {
			Point2D p2d = (Point2D) obj;
			return (getX() == p2d.getX()) && (getY() == p2d.getY());
		}
		return super.equals(obj);
	}

	/**
	 * Euclidean distance to point.
	 * 
	 * @param pt The other point
	 * @return The distance
	 */
	public double distance(Point2D pt) {
		if (pt == null)
			return Double.POSITIVE_INFINITY;
		else {
			double sum = 0;
			for (int i = 0; i < DIM; i++) {
				sum += Math.pow(pt.coord[i] - coord[i], 2);
			}
			return (double) Math.sqrt(sum);
		}
	}

	/**
	 * String representation.
	 * 
	 * @return String representation of the point.
	 */
	public String toString() {
		return "(" + coord[0] + "," + coord[1] + ")";
	}

	/**
	 * Returns the hashcode for this Point2D.
	 * 
	 * @return A hash code for this Point2D.
	 */
	public int hashCode() {
		long bits = java.lang.Double.doubleToLongBits(getX());
		bits ^= java.lang.Double.doubleToLongBits(getY()) * 31;
		return (((int) bits) ^ ((int) (bits >> 32)));
	}

}
