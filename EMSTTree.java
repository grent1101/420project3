package cmsc420_s22;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

/**
 * Implements a variant of Prim's algorithm for constructing a Euclidean minimum
 * spanning tree. Points are added by the operation addPoint. The tree is built
 * by the operation buildEMST. The tree can be listed by listTree.
 */

public class EMSTree<LPoint extends LabeledPoint2D> {

	// -------------------------------------------------------------------------
	// Private members
	// -------------------------------------------------------------------------
	
	private int maxHeightDifference = 2;
	private int quakeMaxLevels = 10;
	private ArrayList<LPoint> pointList; // list of points
	private HashMap<LPoint, ArrayList<LPoint>> dependents; // list of dependent points
	private ArrayList<Pair<LPoint>> edgeList; // partial EMST as a list of edges
	private HashSet<LPoint> inEMST; // set of points that are in tree so far
	private HBkdTree<LPoint> kdTree; // kd-tree for the points
	private QuakeHeap<Double, Pair<LPoint>> heap; // priority queue for Prim's

	// -------------------------------------------------------------------------
	// General utilities. We also include utilities testing for duplicate points
	// and distances, since this may result in variable valid results.
	// -------------------------------------------------------------------------

	/**
	 * Get squared distance between two LPoints, or infinity if either is null.
	 */
	double distanceSq(LPoint pt1, LPoint pt2) {
		if (pt1 == null || pt2 == null)
			return Double.POSITIVE_INFINITY;
		else
			return pt1.getPoint2D().distanceSq(pt2.getPoint2D());
	}

	// -------------------------------------------------------------------------
	// Utilities used by buildEMST
	// -------------------------------------------------------------------------

	/**
	 * Initialize everything before building the EMST
	 */
	void initializeEMST(LPoint start) throws Exception {
		edgeList.clear(); // clear the edge list
		inEMST.clear(); // clear the EMST set
		heap.clear(); // clear the heap
		Iterator<Entry<LPoint, ArrayList<LPoint>>> iter = dependents.entrySet().iterator();
		while (iter.hasNext()) { // clear all the dependents lists
			iter.next().getValue().clear();
		}
		kdTree.clear(); // clear the spatial map
		for (LPoint pt : pointList) { // insert all but start in kd-tree
			if (!start.equals(pt)) {
				kdTree.insert(pt); // insert into kd-tree
			}
		}
		inEMST.add(start); // add start to EMST
	}

	/**
	 * Adds a single edge to the EMST and updates the nearest neighbors.
	 */
	String addEdge(Pair<LPoint> edge) throws Exception {
		LPoint pt2 = edge.getSecond();
		String result = "add: " + edge; // include this in the result
		edgeList.add(edge); // add edge to the EMST
		inEMST.add(pt2); // flag pt2 as being in EMST
		kdTree.delete(pt2.getPoint2D()); // remove pt2 from the spatial map
		ArrayList<LPoint> dep2 = dependents.get(pt2); // get pt2's dependent points
		dep2.add(pt2); // a trick to force pt2 to compute its NN
		ArrayList<String> nnList = new ArrayList<String>();
		for (LPoint pt3 : dep2) { // update dependents
			LPoint nn3 = kdTree.nearestNeighbor(pt3.getPoint2D()); // get pt3's new nearest
			if (nn3 == null)
				break; // out of points? -- we're done
			nnList.add(addNearNeighbor(pt3, nn3));
		}
		Collections.sort(nnList); // sort nnList
		result += " new-nn:";
		for (String item : nnList) { // copy to result
			result += " " + item;
		}
		return result;
	}

	/**
	 * Adds a new nearest-neighbor pair (pt->nn).
	 */
	String addNearNeighbor(LPoint pt, LPoint nn) {
		double dist = distanceSq(pt, nn); // distance to nearest neighbor
		Pair<LPoint> pair = new Pair<LPoint>(pt, nn); // new nearest neighbor pair
		heap.insert(dist, pair); // add to priority queue
		dependents.get(nn).add(pt); // add pt3 to nn3's dependents list
		return "(" + pt.getLabel() + "->" + nn.getLabel() + ")"; // return summary
	}

	/**
	 * Cleanup after building the EMST
	 */
	void cleanUp() {
		inEMST.clear(); // clear the EMST set
		kdTree.clear(); // clear the spatial map
		heap.clear(); // clear the heap
	}

	// -------------------------------------------------------------------------
	// Public members
	// -------------------------------------------------------------------------

	/**
	 * Create an empty spanning tree
	 */
	public EMSTree(Rectangle2D bbox) {
		pointList = new ArrayList<LPoint>();
		dependents = new HashMap<LPoint, ArrayList<LPoint>>();
		edgeList = new ArrayList<Pair<LPoint>>();
		inEMST = new HashSet<LPoint>();
		kdTree = new HBkdTree<LPoint>(maxHeightDifference, bbox);
		heap = new QuakeHeap<Double, Pair<LPoint>>(quakeMaxLevels);
	}

	/**
	 * Adds a point to the point set.
	 */
	public void addPoint(LPoint pt) {
		pointList.add(pt); // add to the point list
		dependents.put(pt, new ArrayList<LPoint>()); // create a new dependents list
	}

	/**
	 * Builds the EMST by a variant of Prim's algorithm for constructing a Euclidean
	 * minimum spanning tree. Points are added by the operation addPoint. The tree
	 * is built by the operation buildEMST. The tree can be listed by listTree.
	 * 
	 * The algorithm works as follows for a given point set P. Given a start point
	 * start, we create a kd-tree containing all the points of P except start. Using
	 * the kd-tree we then compute the nearest neighbor in P to start, call it
	 * startNN. We initialize a priority queue with the edge (start, startNN) whose
	 * priority is the distance between them. Initially, only start is in the EMST.
	 * 
	 * We then repeat the following process until all the points are added to the
	 * EMST. We assume that we have computed a partial spanning tree T containing i
	 * points in the tree and i-1 edges. We maintain a priority queue, Q which
	 * contains one entry for each of the i points. Each entry consists of a pair
	 * (p, q), where p is a point of T, and q is its nearest neighbor that is not in
	 * the tree. The associated priority is the distance between p and q.
	 * 
	 * At each iteration, the algorithm extracts the entry (p,q) from Q with the
	 * smallest distance. We add the edge (p,q) to the tree. We remove q from the
	 * kd-tree. Then for each point r in q's dependency list, we update r's nearest
	 * neighbor and add these pairs to the priority queue.
	 */
	public ArrayList<String> buildEMST(LPoint start) throws Exception {
		ArrayList<String> summary = new ArrayList<String>();
		initializeEMST(start); // initialize everything
		// -------------------------------------------------------------
		// Add the initial near-neighbor pair (start -> start's NN)
		// -------------------------------------------------------------
		LPoint nn = kdTree.nearestNeighbor(start.getPoint2D()); // get start's NN
		if (nn == null)
			return summary; // out of points? -- we're done
		addNearNeighbor(start, nn); // add nearest neighbor pair
		summary.add("new-nn: (" + start.getLabel() + "->" + nn.getLabel() + ")"); // include in result

		// -------------------------------------------------------------
		// Repeated add edges until the kd-tree is empty.
		// -------------------------------------------------------------
		while (kdTree.size() != 0) { // points remain to be added?
			Pair<LPoint> edge = heap.extractMin(); // extract best from priority queue
			if (edge == null) { // there should be something there
				System.err.println("Internal error - Heap is empty");
				throw new Error("Empty heap");
			}
			LPoint pt2 = edge.getSecond(); // get destination end point
			if (!inEMST.contains(pt2)) { // end point not in EMST?
				summary.add(addEdge(edge)); // add the edge to the EMST and update summary
			}
		}
		cleanUp(); // clean up after building the tree
		return summary; // return the final summary
	}

	/**
	 * Returns a copy of the list of the edges in the current tree
	 */
	public ArrayList<String> listEMST() {
		ArrayList<String> list = new ArrayList<String>();
		for (int i = 0; i < edgeList.size(); i++) {
			list.add("(" + edgeList.get(i).getFirst().getLabel() + "," + edgeList.get(i).getSecond().getLabel() + ")");
		}
		return list;
	}

	/**
	 * Clear everything.
	 */
	public void clear() {
		pointList.clear();
		edgeList.clear();
		kdTree.clear();
		heap.clear();
		inEMST.clear();
		dependents.clear();
	}

	/**
	 * Number of entries.
	 */
	public int size() {
		return pointList.size();
	}

}
