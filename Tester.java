package cmsc420_s22;

// YOU SHOULD NOT MODIFY THIS FILE, EXCEPT TO ALTER THE INPUT/OUTPUT SOURCES

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Scanner;

/**
 * Testing program for Tour data structure
 */

public class Tester {

	// --------------------------------------------------------------------------------------------
	// Uncomment these to read from standard input and output
	// private static final boolean USE_STD_IO = true;
	// private static String inputFileName = "";
	// private static String outputFileName = "";
	// --------------------------------------------------------------------------------------------
	// Uncomment these to read from files
	private static final boolean USE_STD_IO = false;
	private static String inputFileName = "tests/test01-input.txt";
	private static String outputFileName = "tests/test01-output.txt";
	// --------------------------------------------------------------------------------------------

	public static void main(String[] args) {

		// configure to read from file rather than standard input/output
		if (!USE_STD_IO) {
			try {
				System.setIn(new FileInputStream(inputFileName));
				System.setOut(new PrintStream(outputFileName));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}

		try {
			Scanner scanner = new Scanner(System.in); // input scanner
			CommandHandler commandHandler = new CommandHandler(); // initialize command handler
			while (scanner.hasNextLine()) {
				String line = scanner.nextLine(); // input next line
				String output = commandHandler.processCommand(line); // process this command
				System.out.print(output); // output summary
			}
			scanner.close();
		} catch (Exception e) {
			System.err.println("Unexpected error: " + e.getMessage());
			e.printStackTrace(System.err);
		}
	}
}
