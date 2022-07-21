package cmsc420_s22;

// YOU SHOULD NOT MODIFY THIS FILE. If you want to add additional functionality,
// create a new class/file (e.g., MyRectangle2D.java) and modify that.

/**
 * A simple 2-dimensional rectangle.
 */

public class Rectangle2D {
	final static int DIM = 2; // spatial dimension
	Point2D low; // lower-left corner
	Point2D high; // upper-right corner

	/**
	 * Construct an "empty" rectangle. This rectangle can serve as the starting
	 * point to a process of expansion by adding points.
	 */
	public Rectangle2D() {
		this.low  = new Point2D(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY);
		this.high = new Point2D(Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY);
	}

	/**
	 * Construct from any two extreme corner points.
	 */
	public Rectangle2D(Point2D c1, Point2D c2) {
		this.low  = new Point2D(Math.min(c1.getX(), c2.getX()), Math.min(c1.getY(), c2.getY()));
		this.high = new Point2D(Math.max(c1.getX(), c2.getX()), Math.max(c1.getY(), c2.getY()));
	}

	/**
	 * Copy constructor
	 */
	public Rectangle2D(Rectangle2D r) {
		low = new Point2D(r.low);
		high = new Point2D(r.high);
	}

	public String toString() {
		return "[" + low + "," + high + "]";
	}

	/**
	 * Getters
	 */

	public Point2D getLow() {
		return low;
	}

	public Point2D getHigh() {
		return high;
	}

	/**
	 * Get width along dimension i.
	 */
	public double getWidth(int i) {
		return high.get(i) - low.get(i);
	}

	/**
	 * Get center point.
	 */
	public Point2D getCenter() {
		return new Point2D((low.getX() + high.getX())/2, (low.getY() + high.getY())/2);
	}

	/**
	 * Check whether we contain a given point.
	 */
	public boolean contains(Point2D q) {
		for (int i = 0; i < DIM; i++) {
			if (q.get(i) < low.get(i) || q.get(i) > high.get(i))
				return false;
		}
		return true;
	}

	/**
	 * Check whether we contain another rectangle.
	 */
	public boolean contains(Rectangle2D c) {
		for (int i = 0; i < DIM; i++) {
			if (c.low.get(i) < low.get(i) || c.high.get(i) > high.get(i))
				return false;
		}
		return true;
	}

	/**
	 * Check whether we are disjoint from another rectangle. Rectangles are closed,
	 * so if their boundaries overlap, they are not disjoint.
	 */
	public boolean disjointFrom(Rectangle2D c) {
		for (int i = 0; i < DIM; i++) {
			if (c.high.get(i) < low.get(i) || c.low.get(i) > high.get(i))
				return true;
		}
		return false;
	}

	/**
	 * Compute the squared Euclidean distance to a point. Returns zero if the point
	 * is contained within this rectangle.
	 */
	public double distanceSq(Point2D pt) {
		double sum = 0; // sum of squared coordinate distances
		for (int i = 0; i < DIM; i++) {
			double coord = pt.get(i); // q's i-th coordinate
			double lc = low.get(i); // low's i-th coordinate
			double hc = high.get(i); // high's i-th coordinate
			if (coord < lc) { // to the left of the rectangle
				sum += Math.pow((lc - coord), 2);
			} else if (coord > hc) {
				sum += Math.pow((coord - hc), 2);
			}
		}
		return sum;
	}

	/**
	 * Compute a minimum rectangle that encloses the union of two
	 * rectangles. To do this, we start with one rectangle and expand
	 * it to include the opposing corners of the other.
	 */
	public static Rectangle2D union(Rectangle2D r1, Rectangle2D r2) {
		Rectangle2D result = new Rectangle2D(r1);
		result.expand(r2.low);
		result.expand(r2.high);
		return result;
	}

	/**
	 * Compute the left part of a rectangle that is split by a 
	 * line orthogonal to the cutting dimension at the given
	 * cutting value. It is assumed that the cut passes through
	 * the rectangle. (No error checking!)
	 */
	public Rectangle2D leftPart(int cutDim, double cutVal) {
		Rectangle2D result = new Rectangle2D(this);
		result.high.set(cutDim, cutVal);
		return result;
	}

	/**
	 * Compute the right part of a rectangle that is split by a 
	 * line orthogonal to the cutting dimension at the given
	 * cutting value. It is assumed that the cut passes through
	 * the rectangle. (No error checking!)
	 */
	public Rectangle2D rightPart(int cutDim, double cutVal) {
		Rectangle2D result = new Rectangle2D(this);
		result.low.set(cutDim, cutVal);
		return result;
	}

	/**
	 * Expand the rectangle to contain the given point.
	 */
	public void expand(Point2D pt) {
		for (int i = 0; i < DIM; i++) {
			double coord = pt.get(i);
			if (coord < low.get(i))
				low.set(i, coord);
			if (coord > high.get(i))
				high.set(i, coord);
		}
	}

}
