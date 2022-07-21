package cmsc420_s22;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * A height-balanced kd-tree. This height-balanced kd-tree is an extended
 * version of a kd-tree, which is balanced by height in the sense that for each
 * internal node u, the difference in heights between its left and right
 * subtrees is bounded by a parameter, called the height-difference.
 *
 * Generic Elements: The tree is parameterized by an LPoint, which implements
 * the LabeledPoint2D interface.
 * 
 * In addition to its left and right children, each internal node stores a
 * cutting dimension (either 0 for x or 1 for y), a cutting value, and the
 * height of the subtree rooted here. Each external node stores just a single
 * point.
 */

public class HBkdTree<LPoint extends LabeledPoint2D> {

	// =================================================================
	// KDNode
	//
	// Each node contains a point and cutting dimension and links to its
	// two children. To help with height balancing, we store the height
	// of the subtree rooted at this node.
	// =================================================================

	private class KDNode {

		private LPoint point; // the associated point
		private int cutDim; // cutting dimension (0 == x, 1 == y)
		private int height; // subtree height
		private KDNode left, right; // children

		public KDNode(LPoint point, int cutDim, KDNode left, KDNode right) { // standard constructor
			this.point = point;
			this.cutDim = cutDim;
			this.left = left;
			this.right = right;
			updateStats();
		}

		public KDNode(LPoint point, int cutDim) { // leaf constructor
			this.point = point;
			this.cutDim = cutDim;
			height = 0;
			left = right = null;
		}

		void updateStats() { // update subtree height
			height = 1 + Math.max(getHeight(left), getHeight(right));
		}

		boolean onLeft(LPoint pt) { // in the left subtree? (for Labeled points)
			return pt.get(cutDim) < point.get(cutDim);
		}

		boolean onLeft(Point2D pt) { // in the left subtree? (for points)
			return pt.get(cutDim) < point.get(cutDim);
		}

		boolean onSplit(Point2D pt) { // on the splitting line? (for points)
			return pt.get(cutDim) == point.get(cutDim);
		}

		Rectangle2D leftPart(Rectangle2D cell) { // left part of cell
			return cell.leftPart(cutDim, point.get(cutDim));
		}

		Rectangle2D rightPart(Rectangle2D cell) { // right part of cell
			return cell.rightPart(cutDim, point.get(cutDim));
		}

		public String toString() { // string representation
			String cut = (cutDim == 0 ? "x" : "y");
			return "(" + cut + "=" + point.get(cutDim) + " ht=" + height + ") " + point.toString();
		}
	}

	// -----------------------------------------------------------------
	// Utilities
	// -----------------------------------------------------------------
	int getHeight(KDNode p) { // get subtree height
		return (p == null ? -1 : p.height);
	}

	double distanceSq(Rectangle2D cell, LPoint pt) { // distance cell to LPoint
		return cell.distanceSq(pt.getPoint2D());
	}

	double distanceSq(Point2D pt1, LPoint pt2) { // distance Point2D to LPoint
		if (pt2 == null)
			return Double.POSITIVE_INFINITY;
		else
			return pt1.distanceSq(pt2.getPoint2D());
	}

	// -----------------------------------------------------------------
	// Recursive helpers for main functions
	// -----------------------------------------------------------------

	/**
	 * Finds a point in the node's subtree.
	 */
	LPoint find(KDNode p, Point2D pt) { // find point in subtree
		if (p == null) {
			return null;
		} else if (p.point.getPoint2D().equals(pt)) {
			return p.point;
		} else if (p.onLeft(pt)) {
			return find(p.left, pt);
		} else {
			if (p.onSplit(pt)) {
				LPoint lResult = find(p.left, pt);
				if (lResult != null) {
					System.err.println("Warning: Find " + pt + " returns left-side point");
				}
			}
			return find(p.right, pt);
		}
	}

	/**
	 * Insert a point in the node's subtree.
	 */
	KDNode insert(LPoint pt, KDNode p, Rectangle2D cell) throws Exception { // insert point into subtree
		if (p == null) {
			return new KDNode(pt, getCutDim(cell));
		} else if (pt.getPoint2D().equals(p.point.getPoint2D())) {
			throw new Exception("Attempt to insert a duplicate point");
		} else if (p.onLeft(pt)) { // insert on appropriate side
			p.left = insert(pt, p.left, p.leftPart(cell));
		} else {
			p.right = insert(pt, p.right, p.rightPart(cell));
		}
		return rebalance(p, cell);
	}

	/**
	 * Delete a point from node's subtree.
	 */
	KDNode delete(Point2D pt, KDNode p, Rectangle2D cell) throws Exception {
		if (p == null) { // fell out of tree?
			throw new Exception("Attempt to delete a nonexistent point");
		} else if (pt.equals(p.point.getPoint2D())) { // found it
			if (p.right != null) { // can replace from right
				p.point = findMin(p.right, p.cutDim); // find and copy replacement
				p.right = delete(p.point.getPoint2D(), p.right, p.rightPart(cell)); // delete from right
			} else if (p.left != null) { // can replace from left
				p.point = findMin(p.left, p.cutDim); // find and copy replacement
				p.right = delete(p.point.getPoint2D(), p.left, p.leftPart(cell)); // delete left but move to right!!
				p.left = null; // left subtree is now empty
			} else { // deleted point in leaf
				p = null; // remove this leaf
			}
		} else if (p.onLeft(pt)) {
			p.left = delete(pt, p.left, p.leftPart(cell)); // delete from left subtree
		} else { // delete from right subtree
			p.right = delete(pt, p.right, p.rightPart(cell));
		}
		return rebalance(p, cell);
	}

	/**
	 * Find min node in subtree along coordinate i.
	 */
	LPoint findMin(KDNode p, int i) {
		if (p == null) { // fell out of tree?
			return null;
		} else if (p.cutDim == i) { // cutting dimension matches i?
			if (p.left == null) { // no left child?
				return p.point; // use this point
			} else {
				return findMin(p.left, i); // get min from left subtree
			}
		} else { // check both sides and this point as well
			return min(i, p.point, min(i, findMin(p.left, i), findMin(p.right, i)));
		}
	}

	/**
	 * Return the minimum non-null point w.r.t. coordinate i.
	 */
	LPoint min(int i, LPoint pt1, LPoint pt2) {
		if (pt1 == null) {
			return pt2;
		} else if (pt2 == null) {
			return pt1;
		} else if (pt1.get(i) < pt2.get(i)) {
			return pt1;
		} else if (pt2.get(i) < pt1.get(i)) {
			return pt2;
		} else {
			int j = 1 - i; // swap coordinate
			return (pt1.get(j) < pt2.get(j) ? pt1 : pt2);
		}
	}

	/**
	 * Builds a preorder list of node contents.
	 */
	ArrayList<String> getPreorderList(KDNode p) { // list entries in preorder
		ArrayList<String> list = new ArrayList<String>();
		if (p == null) {
			list.add("[]"); // null node indicator
		} else {
			list.add(p.toString()); // add this node
			list.addAll(getPreorderList(p.left)); // process left
			list.addAll(getPreorderList(p.right)); // process right
		}
		return list;
	}

	/**
	 * Returns list of points lying with a rectangle.
	 */
	ArrayList<LPoint> orthogRangeReport(Rectangle2D query, KDNode p, Rectangle2D cell) { // orthogonal range query
		ArrayList<LPoint> list = new ArrayList<LPoint>();
		if (p == null) {
			return list;
		} else if (query.disjointFrom(cell)) { // no overlap with query region?
		} else {
			if (query.contains(p.point.getPoint2D())) {
				list.add(p.point); // add this node's point
			}
			list.addAll(orthogRangeReport(query, p.left, p.leftPart(cell)));
			list.addAll(orthogRangeReport(query, p.right, p.rightPart(cell)));
		}
		return list;
	}

	/**
	 * Nearest neighbor. Returns the closest point to center. The best point seen so
	 * far is best. If we fall out of the tree, we return best. If the point stored
	 * here is closer than best, we update best. Otherwise, we search both subtrees
	 * in order according to which side is closer.
	 */

	LPoint nearestNeighbor(Point2D q, KDNode p, Rectangle2D cell, LPoint best) {
		if (p != null) {
			if (distanceSq(q, p.point) < distanceSq(q, best)) { // better than best?
				best = p.point; // new best
			}
			Rectangle2D leftCell = p.leftPart(cell); // left child's cell
			Rectangle2D rightCell = p.rightPart(cell); // right child's cell
			int cd = p.cutDim; // cutting dimension
			if (q.get(cd) < p.point.get(cd)) { // q is closer to left
				best = nearestNeighbor(q, p.left, leftCell, best);
				if (rightCell.distanceSq(q) < distanceSq(q, best)) { // worth visiting right?
					best = nearestNeighbor(q, p.right, rightCell, best);
				}
			} else { // q is closer to right
				best = nearestNeighbor(q, p.right, rightCell, best);
				if (leftCell.distanceSq(q) < distanceSq(q, best)) { // worth visiting left?
					best = nearestNeighbor(q, p.left, leftCell, best);
				}
			}
		}
		return best;
	}

	// -----------------------------------------------------------------
	// Other utilities
	// -----------------------------------------------------------------

	String debugPrint(KDNode p, String prefix) { // print for debugging
		String left = new String();
		String right = new String();
		if (p.left != null)
			left = System.lineSeparator() + debugPrint(p.left, prefix + "| ");
		if (p.right != null)
			right = debugPrint(p.right, prefix + "| ") + System.lineSeparator();
		return right + prefix + p.toString() + left;
	}

	// -----------------------------------------------------------------
	// Rebalancing utilities
	// -----------------------------------------------------------------

	/**
	 * Rebalance subtree by updating height, checking height condition, and
	 * rebuilding the subtree if needed.
	 */
	KDNode rebalance(KDNode p, Rectangle2D cell) {
		if (p == null)
			return null;
		p.updateStats(); // update p's height and size
		if (Math.abs(getHeight(p.left) - getHeight(p.right)) > maxHeightDifference) {
			return rebuild(p, cell);
		} else {
			return p;
		}
	}

	/**
	 * Rebuild the subtree rooted at this node. We first compute a sorted list of
	 * external nodes and then rebuild the subtree recursively.
	 */
	KDNode rebuild(KDNode p, Rectangle2D cell) {
		List<LPoint> list = buildPointList(p); // collect all the points
		return buildSubtree(list, cell);
	}

	/**
	 * Builds a list of external nodes in inorder.
	 */
	List<LPoint> buildPointList(KDNode p) {
		ArrayList<LPoint> list = new ArrayList<LPoint>();
		if (p != null) {
			list.addAll(buildPointList(p.left)); // add left
			list.add(p.point); // add this point
			list.addAll(buildPointList(p.right)); // add right
		}
		return list;
	}

	/**
	 * Recursively rebuilds the subtree rooted at this node given a sorted list of
	 * external nodes and their weights. If the list is empty, we return a null
	 * tree. Otherwise, it computes the median and splits about this node. It
	 * recursively builds subtrees for these sublists and joins them under a common
	 * node that splits at the median.
	 */

	KDNode buildSubtree(List<LPoint> list, Rectangle2D cell) {
		if (list.size() == 0) {
			return null; // copy this node
		} else {
			// Rectangle2D cell = getBoundingBox(list); // bounding box for points in list
			int cutDim = getCutDim(cell); // get the cutting dimension
			if (cutDim == 0) { // sort by the cutting dimension
				Collections.sort(list, new ByXThenY());
			} else {
				Collections.sort(list, new ByYThenX());
			}
			int size = list.size();
			int m = size / 2; // median index
			LPoint medPt = list.get(m); // median point
			// split and recurse
			Rectangle2D leftCell = cell.leftPart(cutDim, medPt.get(cutDim));
			Rectangle2D rightCell = cell.rightPart(cutDim, medPt.get(cutDim));
			KDNode left = buildSubtree(list.subList(0, m), leftCell); // sublist from [0..m-1]
			KDNode right = buildSubtree(list.subList(m + 1, size), rightCell); // sublist from [m..size-1]
			return new KDNode(medPt, cutDim, left, right); // join under internal node
		}
	}

	/**
	 * Utilities for sorting and splitting.
	 */

	private static int getCutDim(Rectangle2D cell) { // get the wider of the sides
		// System.out.println("Getting cutting dimension of cell of width=" +
		// cell.getWidth(0) + " and height=" + cell.getWidth(1));
		return (cell.getWidth(0) >= cell.getWidth(1) ? 0 : 1);
	}

	private class ByXThenY implements Comparator<LPoint> { // lexicographic (x,y)
		public int compare(LPoint pt1, LPoint pt2) {
			double x1 = pt1.getX();
			double x2 = pt2.getX();
			if (x1 < x2)
				return -1;
			else if (x1 > x2)
				return +1;
			else {
				double y1 = pt1.getY();
				double y2 = pt2.getY();
				if (y1 < y2)
					return -1;
				else if (y1 > y2)
					return +1;
				else {
					System.err.println("Duplicate coordinates!");
					return 0;
				}
			}
		}
	}

	private class ByYThenX implements Comparator<LPoint> { // lexicographic (y,x)
		public int compare(LPoint pt1, LPoint pt2) {
			double y1 = pt1.getY();
			double y2 = pt2.getY();
			if (y1 < y2)
				return -1;
			else if (y1 > y2)
				return +1;
			else {
				double x1 = pt1.getX();
				double x2 = pt2.getX();
				if (x1 < x2)
					return -1;
				else if (x1 > x2)
					return +1;
				else {
					System.err.println("Duplicate coordinates!");
					return 0;
				}
			}
		}
	}

	// -----------------------------------------------------------------
	// Private data
	// -----------------------------------------------------------------

	private KDNode root; // root of the tree
	private int nPoints; // number of points in the tree
	private int maxHeightDifference; // allowed height difference
	private Rectangle2D bbox; // the bounding box

	// -----------------------------------------------------------------
	// Public members
	// -----------------------------------------------------------------

	/**
	 * Creates an empty tree.
	 */
	public HBkdTree(int maxHeightDifference, Rectangle2D bbox) {
		root = null;
		nPoints = 0;
		this.bbox = new Rectangle2D(bbox);
		this.maxHeightDifference = maxHeightDifference;
	}

	/**
	 * Number of entries in the dictionary.
	 */
	public int size() {
		return nPoints;
	}

	/**
	 * Find an point in the tree.
	 */
	public LPoint find(Point2D pt) {
		return find(root, pt);
	}

	/**
	 * Insert a point.
	 */
	public void insert(LPoint pt) throws Exception {
		if (!bbox.contains(pt.getPoint2D())) {
			throw new Exception("Attempt to insert a point outside bounding box");
		} else {
			root = insert(pt, root, bbox); // insert the point
		}
		nPoints += 1; // one more point
	}

	/**
	 * Delete a point. Note that the point being deleted does not need to match
	 * fully. It suffices that it has enough information to satisfy the comparator.
	 */
	public void delete(Point2D pt) throws Exception {
		root = delete(pt, root, bbox); // delete the point
		nPoints -= 1; // one fewer point
	}

	/**
	 * Get a list of entries in preorder
	 */
	public ArrayList<String> getPreorderList() {
		return getPreorderList(root);
	}

	/**
	 * Rectangular range query.
	 */
	public ArrayList<LPoint> orthogRangeReport(Rectangle2D query) {
		ArrayList<LPoint> result = orthogRangeReport(query, root, bbox);
		return result;
	}

	/**
	 * Nearest neighbor query.
	 */
	public LPoint nearestNeighbor(Point2D center) {
		LPoint best = null;
		LPoint result = nearestNeighbor(center, root, bbox, best);
		return result;
	}

	/**
	 * Remove all items, resulting in an empty tree
	 */
	public void clear() {
		root = null;
		nPoints = 0;
	}

	/**
	 * Delete a point. Note that the point being deleted does not need to match
	 * fully. It suffices that it has enough information to satisfy the comparator.
	 */
	public void setHeightDifference(int newDiff) throws Exception {
		if (newDiff < 1) {
			throw new Exception("Height difference must be at least 1");
		} else {
			maxHeightDifference = newDiff;
		}
	}

	// -----------------------------------------------------------------
	// Debugging utilities
	// -----------------------------------------------------------------

	/**
	 * Print the tree for debugging purposes
	 */
	String debugPrint(String prefix) {
		if (root == null)
			return new String(prefix + "[Empty]");
		else
			return debugPrint(root, prefix);
	}

}
