package cmsc420_s22;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;

public class QuakeHeap<Key extends Comparable<Key>, Value> {

	// -----------------------------------------------------------------
	// Node
	// -----------------------------------------------------------------

	class Node {
		Key key; // key (for sorting)
		Value value; // value (application dependent)
		int level; // level in the tree (leaf = level 0)
		Node left; // children
		Node right;
		Node parent; // parent
		Node leftmostLeaf; // leaf node along left chain
		Node topLeftNode; // highest node along left ascending chain

		/**
		 * Basic constructor.
		 */
		Node(Key x, Value v, int level, Node left, Node right, Node leftmostLeaf, Node topLeftNode) {
			this.key = x;
			this.value = v;
			this.level = level;
			this.left = left;
			this.right = right;
			this.parent = null;
			this.leftmostLeaf = leftmostLeaf;
			this.topLeftNode = topLeftNode;
		}
		
		Key getKey() {
			return leftmostLeaf.key;
		}
		
		Value getValue() {
			return leftmostLeaf.value;
		}

	}

	/**
	 * Node comparator used for sorting root nodes.
	 */

	private class ByKey implements Comparator<Node> {
		public int compare(Node u, Node v) {
			return u.getKey().compareTo(v.getKey());
		}
	}

	/**
	 * Compares nodes based on keys. If either node is null,
	 * we treat the other as being smaller.
	 */
	private int compare(Node u, Node v) {
		if (u == null) return +1;
		if (v == null) return -1;
		return u.getKey().compareTo(v.getKey());
	}

	// -----------------------------------------------------------------
	// Locator - Used to locate a previously inserted item
	// The Node reference always points to a leaf node.
	// -----------------------------------------------------------------

	public class Locator {
		private final Node u;

		private Locator(Node u) { // basic constructor
			this.u = u;
		}

		private Node get() { // get the associated node
			return u;
		}
	}

	// -----------------------------------------------------------------
	// Private members
	// -----------------------------------------------------------------

	private final double defaultQuakeRatio = 0.75f; // default quake ratio
	private double quakeRatio; // ratio used for triggering quake
	private int nEntries; // number of entries (for size)
	private int nLevels; // number of levels
	private LinkedList<Node>[] roots; // list of roots per level
	private int[] nodeCt; // number of nodes per level

	// -----------------------------------------------------------------
	// Local utilities
	// -----------------------------------------------------------------

	/**
	 * Add a new node as root.
	 */
	void makeRoot(Node u) {
		u.parent = null; // null out parent link
		roots[u.level].add(u); // add node at u's level
	}

	/**
	 * Create a new leaf node.
	 * 
	 * Creates a new leaf node at level 0. The leftmostLeaf and topLeftNode 
	 * entries are set only if doing fast decrease-key.
	 */
	Node newLeafNode(Key x, Value v) {
		Node u = new Node(x, v, 0, null, null, null, null);
		u.leftmostLeaf = u.topLeftNode = u;
		return u;
	}

	/**
	 * Create a trivial single-node tree
	 */
	Node trivialTree(Key x, Value v) { // create a trivial single-node tree
		Node u = newLeafNode(x, v); // create new leaf node
		nodeCt[0] += 1; // increment node count
		makeRoot(u); // make it a root
		return u;
	}
	
	/**
	 * Create a new internal node.
	 * 
	 * Creates a new internal node at the specified level that is a parent to
	 * the two given nodes.
	 */
	Node newInternalNode(Node u, Node v) {
		Node w;
		if (compare(u, v) <= 0) { // u's key is smaller?
			w = new Node(u.key, u.value, u.level+1, u, v, null, null);
		} else {
			w = new Node(v.key, v.value, v.level+1, v, u, null, null);
		}
		w.key = null;
		w.value = null;
		w.leftmostLeaf = w.left.leftmostLeaf;
		w.leftmostLeaf.topLeftNode = w;
		return w;
	}

	/**
	 * Link trees rooted at u and v (from same level) together to form a new tree
	 * (one level higher). By convention, the smaller key is stored in the left
	 * subtree. The tasks of removing u and v from the list of roots and adding the
	 * new root w is handled in the calling function.
	 * 
	 * For fast decrease-key, we to not set the key-value pair, but we do set the
	 * leftmostLeaf pointer, and we update the associated topLeftNode. Otherwise, we
	 * copy the key-value pair from the left child.
	 */
	Node link(Node u, Node v) { // link u and v into new tree
		assert (u.level == v.level); // nodes must be at the same level
		Node w = newInternalNode(u, v); // join under new internal node
		nodeCt[w.level] += 1; // increment node count
		u.parent = v.parent = w; // w is the new parent
		return w;
	}

	/**
	 * Cut u's right child, making it a new root.
	 */
	void cut(Node u) { // cut off u's right child
		Node v = u.right;
		if (v != null) {
			u.right = null; // cut off v
			makeRoot(v); // ... and make it a root
		}
	}

	/**
	 * Search the roots of all the trees and return a reference to the one having
	 * the smallest key value.
	 */
	Node findRootWithSmallestKey() {
		Node min = null;
		for (int lev = 0; lev < nLevels; lev++) { // process all levels
			for (Node u : roots[lev]) {
				if (min == null || compare(u, min) < 0) {
					min = u;
				}
			}
		}
		return min;
	}

	/**
	 * Delete the leftmost path for a root node u. (By our convention, all of these
	 * nodes have the same key.) We unlink each right child using the cut operation
	 * and decrement the node count at each level. We leave the task of removing u
	 * as a root to the calling function.
	 */
	void deleteLeftPath(Node u) {
		while (u != null) { // repeat until falling out of the tree
			cut(u); // cut off u's right child
			nodeCt[u.level] -= 1; // u is gone from this level
			u = u.left; // go to the left child
		}
	}

	/**
	 * Merge all pairs of trees at the same level. We work bottom-up because merging
	 * two trees creates a tree one level higher, which can then be merged with
	 * others. Note that we stop merging at nLevels-2, since we cannot create nodes
	 * above that level.
	 * 
	 * Alert: We sort the roots by key value. This is not part of the QuakeHeap
	 * algorithm. It is done for the sake of having deterministic behavior.
	 */
	void mergeTrees() {
		for (int lev = 0; lev < nLevels - 1; lev++) { // process levels bottom-up
			Collections.sort(roots[lev], new ByKey()); // sort roots by key
			while (roots[lev].size() >= 2) { // at least two trees?
				Node u = roots[lev].remove(); // remove two trees
				Node v = roots[lev].remove();
				Node w = link(u, v); // ... and merge them
				makeRoot(w); // ... and make this a root
			}
		}
	}

	/**
	 * Clear all nodes strictly above level top. Working top-down, we remove all
	 * root nodes and make their two children into roots.
	 */
	void clearAllAboveLevel(int top) { // clear all nodes above top
		for (int lev = nLevels - 1; lev > top; lev--) { // process all levels
			while (roots[lev].size() > 0) {
				Node u = roots[lev].remove(); // remove root u from this level
				unRoot(u);
			}
			nodeCt[lev] = 0; // zero out the node count
		}
	}

	/**
	 * Un-make u as a root. It is assumed that u is a root. This makes u's children
	 * into roots (assuming they are non-null). If we are using fast decrease-key, 
	 * u's leftmost leaf now points to u's left child.
	 */
	void unRoot(Node u) {
		if (u.left != null)
			makeRoot(u.left); // make u's children roots
		if (u.right != null)
			makeRoot(u.right);
		Node ll = u.leftmostLeaf; // u's leftmost leaf
		ll.topLeftNode = u.left; // update top-left node to u			
	}

	/**
	 * Flatten the structure if needed. This checks each level and if the one above
	 * as more the 3/4 of what we have here, we clear everything above the current
	 * level.
	 * 
	 * For testing purposes, if the quake takes place, we return the level at which
	 * the quake occurs, and -1 otherwise.
	 */
	int quake() { // flatten if needed
		for (int lev = 0; lev < nLevels - 1; lev++) { // process all levels
			if (nodeCt[lev + 1] > quakeRatio * nodeCt[lev]) { // too many nodes above?
				clearAllAboveLevel(lev); // clear all above this level
				return lev;
			}
		}
		return -1; // no quake
	}

	/**
	 * Get a list of the nodes in preorder of a single subtree.
	 */
	ArrayList<String> getPreorderList(Node u) {
		ArrayList<String> list = new ArrayList<String>();
		if (u == null) {
			list.add("[null]");
		} else if (u.level > 0) {
			list.add("(" + u.getKey() + ")");
			list.addAll(getPreorderList(u.left));
			list.addAll(getPreorderList(u.right));
		} else {
			list.add("[" + u.getKey() + " " + u.getValue() + "]");
		}
		return list;
	}

	// -----------------------------------------------------------------
	// Public members
	// -----------------------------------------------------------------

	/**
	 * Default constructor.
	 */
	@SuppressWarnings("unchecked")
	public QuakeHeap(int nLevels) {
		this.nLevels = nLevels;
		quakeRatio = defaultQuakeRatio;
		roots = new LinkedList[nLevels];
		nodeCt = new int[nLevels];
		for (int i = 0; i < nLevels; i++) {
			roots[i] = new LinkedList<Node>();
			nodeCt[i] = 0;
		}
		nEntries = 0;
	}

	/**
	 * Return the size (number of keys) in the heap.
	 */
	public int size() {
		return nEntries;
	}

	/**
	 * Clear the entire structure.
	 */
	public void clear() {
		for (int lev = 0; lev < nLevels; lev++) {
			roots[lev].clear();
			nodeCt[lev] = 0;
		}
		nEntries = 0;
	}

	/**
	 * Set a new quake ratio. (Should be between 1/2 and 1. A value of 1/2 causes
	 * the structure to be completely flat, and a value of 1 never flattens.)
	 */
	public void setQuakeRatio(double newRatio) throws Exception {
		if (newRatio < 0.5 || newRatio > 1.0) {
			throw new Exception("Quake ratio is outside valid bounds");
		} else {
			quakeRatio = newRatio;
		}
	}

	/**
	 * Set number of levels. If the number of levels has fallen, we clear all nodes
	 * above the new number of levels. Next, we save the current state of the
	 * levels. We then allocate new arrays of the desired size, and copy the saved
	 * contents here.
	 */
	@SuppressWarnings("unchecked")
	public void setNLevels(int nl) throws Exception {
		if (nl < 1) { // need at least one level
			throw new Exception("Attempt to set an invalid number of levels");
		}
		if (nl < nLevels) { // clear out higher levels
			clearAllAboveLevel(nl - 1);
		}
		LinkedList<Node>[] saveRoots = new LinkedList[nLevels];
		int[] saveNodeCt = new int[nLevels];
		for (int i = 0; i < nLevels; i++) { // save the current root and counts
			saveRoots[i] = roots[i];
			saveNodeCt[i] = nodeCt[i];
		}
		roots = new LinkedList[nl]; // allocate new arrays
		nodeCt = new int[nl];
		for (int i = 0; i < nl; i++) { // copy the old contents over
			if (i < nLevels) {
				roots[i] = saveRoots[i];
				nodeCt[i] = saveNodeCt[i];
			} else {
				roots[i] = new LinkedList<Node>();
				nodeCt[i] = 0;
			}
		}
		nLevels = nl; // update the number of levels
	}

	/**
	 * Insert key-value pair.
	 */
	public Locator insert(Key x, Value v) {
		Node u = trivialTree(x, v); // create a one-node tree storing x
		nEntries += 1; // one more entry
		return new Locator(u); // return a reference to it
	}

	/**
	 * Decrease key for item at location r to y.
	 * 
	 * We implement two versions. In the simple version, we walk up the left-child
	 * path, updating keys. A cut is performed at the parent of the last node in
	 * this path. In the faster version, we change only leftmost leaf key and then
	 * use the top-left node link to find where to perform the cut.
	 */
	public void decreaseKey(Locator r, Key newKey) throws Exception {
		Node u = r.get(); // leaf node to be changed
		Node cutPt = u; // node where cut will be applied
		if (newKey.compareTo(u.key) > 0) { // weight shouldn't increase
			throw new Exception("Invalid key for decrease-key");
		}
		u.key = newKey; // update the leaf key value
		cutPt = u.topLeftNode.parent;
		if (cutPt != null) { // is there a node to cut?
			cut(cutPt); // cut u's subtree off from parent
		}
	}

	/**
	 * Get the minimum key from the heap. In addition to returning the minimum, this
	 * also applies the merging and quaking part of the reorganization.
	 */
	public Key getMinKey() throws Exception {
		if (size() == 0) {
			throw new Exception("Empty heap");
		}
		Node u = findRootWithSmallestKey(); // find the min root
		mergeTrees(); // merge trees
		return u.getKey();
	}

	/**
	 * Get the maximum level of for the entry specified by locator r. This is
	 * defined to be the length of the longest chain of left-child links leading to
	 * the leaf node referenced by r.
	 */
	public int getMaxLevel(Locator r) {
		Node u = r.get(); // leaf node to be changed
		while (u.parent != null && u == u.parent.left) { // climb up left links
			u = u.parent;
		}
		return u.level;
	}

	/**
	 * Extract the minimum item from the heap.
	 */
	public Value extractMin() throws Exception {
		if (size() == 0) {
			throw new Exception("Empty heap");
		}
		Node u = findRootWithSmallestKey(); // find the min root
		Value result = u.getValue(); // final return result
		if (result == null) {
			System.err.println("Null result from extractMin");
		}
		deleteLeftPath(u); // delete entire left path
		roots[u.level].remove(u); // remove u from this level
		mergeTrees(); // merge trees
		quake(); // perform the quake operation
		nEntries -= 1; // one fewer entry
		return result;
	}

	/**
	 * Get a list of entries in preorder.
	 */
	public ArrayList<String> listHeap() {
		ArrayList<String> list = new ArrayList<String>();
		for (int lev = 0; lev < nLevels; lev++) {
			if (nodeCt[lev] > 0) {
				list.add("{lev: " + lev + " nodeCt: " + nodeCt[lev] + "}");
			}
			if (roots[lev].size() > 0) { // has at least one root?
				Collections.sort(roots[lev], new ByKey()); // sort roots by key
				for (Node u : roots[lev]) {
					list.addAll(getPreorderList(u));
				}
			}
		}
		return list;
	}

}
