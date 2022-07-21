package cmsc420_s22;

// YOU SHOULD NOT MODIFY THIS FILE

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 * Command handler. This inputs a single command line, processes the command (by
 * invoking the appropriate method(s) on the EMST.
 */

public class CommandHandler {

	private boolean initialized; // have we initialized the structure yet?
	private HashMap<String, Airport> airports; // airport codes seen so far
	private EMSTree<Airport> emsTree; // the EMST

	/**
	 * Initialize command handler
	 */
	public CommandHandler() {
		airports = new HashMap<String, Airport>();
		emsTree = null;
	}

	/**
	 * Process a single command and return the string output. Each line begins with
	 * a command (e.g., find, insert, delete) followed by a list of arguments. The
	 * arguments are separated by colons (":").
	 * 
	 * @param inputLine The input line with the command and parameters.
	 * @return A short summary of the command's execution/result.
	 */

	public String processCommand(String inputLine) throws Exception {
		String output = new String(); // for storing summary output
		Scanner line = new Scanner(inputLine);
		try {
			line.useDelimiter(":"); // use ":" to separate arguments
			String cmd = (line.hasNext() ? line.next() : ""); // next command
			// -----------------------------------------------------
			// INITIALIZE
			// - this command must come first in the input
			// - sets the bounding box
			// -----------------------------------------------------
			if (cmd.compareTo("initialize") == 0) {
				double xMin = line.nextDouble(); // bounding box
				double xMax = line.nextDouble();
				double yMin = line.nextDouble();
				double yMax = line.nextDouble();
				if (xMin > xMax || yMin > yMax) {
					throw new Exception("Error - invalid bounding box dimensions");
				}
				if (initialized) {
					throw new Exception("Error - Second attempt to initialize");
				} else {
					Rectangle2D bbox = new Rectangle2D(new Point2D(xMin, yMin), new Point2D(xMax, yMax));
					emsTree = new EMSTree<Airport>(bbox); // create a new tree
					output += "initialize: bounding-box = " + bbox
							+ System.lineSeparator();
					initialized = true;
				}
			}
			// -----------------------------------------------------
			// COMMENT string
			// - comment line for the output
			// -----------------------------------------------------
			else if (cmd.compareTo("comment") == 0) {
				String comment = line.next(); // read the comment
				output += "[" + comment + "]" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// ADD-POINT code city x y
			// -----------------------------------------------------
			else if (cmd.compareTo("add-point") == 0) {
				confirmInitialized(); // confirm that we are initialized
				String code = line.next(); // get parameters
				String city = line.next();
				double x = line.nextDouble();
				double y = line.nextDouble();
				Airport ap = new Airport(code, city, x, y); // create airport object
				output += "add-point(" + code + "): ";
				Airport ap2 = airports.get(code);
				if (ap2 != null) { // code already exists?
					throw new Exception("Insertion of duplicate airport code");
				}
				emsTree.addPoint(ap); // insert into the tree's point set
				airports.put(code, ap); // insert into dictionary
				output += "successful {" + ap.getString("attributes") + "}" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CLEAR - Clear all the points added
			// -----------------------------------------------------
			else if (cmd.compareTo("clear") == 0) {
				confirmInitialized(); // confirm that we are initialized
				emsTree.clear(); // get the emst's cost
				airports.clear(); // clear the airports map
				output += "clear: successful" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// SIZE
			// -----------------------------------------------------
			else if (cmd.compareTo("size") == 0) {
				confirmInitialized(); // confirm that we are initialized
				int size = emsTree.size(); // get the tree's current size
				output += "size: " + size + System.lineSeparator();
			}
			// -----------------------------------------------------
			// BUILD-EMST code 
			// - build EMST starting at the given point
			// -----------------------------------------------------
			else if (cmd.compareTo("build-emst") == 0) {
				confirmInitialized(); // confirm that we are initialized
				String code = line.next(); // get start node's name			
				output += "build-emst(" + code + "): ";
				Airport ap = airports.get(code); // look up the airport
				if (ap == null) { // no such airport?
					throw new Exception("Start point is does not exist");
				}			
				ArrayList<String> summary = emsTree.buildEMST(ap); // build the EMST
				output += "successful" + System.lineSeparator();
				for (String item : summary) {
					output += " " + item + System.lineSeparator();
				}
			}
			// -----------------------------------------------------
			// LIST-EMST - list the edges of the EMST
			// -----------------------------------------------------
			else if (cmd.compareTo("list-emst") == 0) {
				confirmInitialized(); // confirm that we are initialized
				ArrayList<String> list = emsTree.listEMST();
				if (list == null) throw new Exception("Error - list returned a null result");
				output += "list-emst:";
				for (String edge : list) { // output the list
					output += " " + edge;
				}
				output += System.lineSeparator();
			}
			// 
			// -----------------------------------------------------
			// Invalid command or empty
			// -----------------------------------------------------
			else {
				if (cmd.compareTo("") == 0)
					System.err.println("Error: Empty command line (Ignored)");
				else
					System.err.println("Error: Invalid command - \"" + cmd + "\" (Ignored)");
			}
		} catch (Exception e) { // exception thrown?
			output += "failure due to exception: \"" + e.getMessage() + "\"" + System.lineSeparator();
		} catch (Error e) { // error occurred?
			System.err.print("Operation failed due to error: " + e.getMessage());
			e.printStackTrace(System.err);
		} finally { // always executed
			line.close(); // close the input scanner
		}
		return output; // return summary output
	}

	/**
	 * Confirm that the data structure has been initialized, or throw an exception.
	 */
	void confirmInitialized() throws Exception {
		if (!initialized) {
			throw new Exception("Error: First command must be 'initialize'.");
		}
	}
}
