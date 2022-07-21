package cmsc420_s22;

/**
 * A generic ordered pair of objects of type X.
 *
 * Stolen from http://www.java2s.com/example/java-src/pkg/javafx/util/pair-7999d.html and
 * modified until unrecognizable.
 */

public class Pair<X>  {

    private X first; // the key
    private X second; // the value

    /**
     * Standard getters
     */
    public X getFirst() {
        return first;
    }

    public X getSecond() {
        return second;
    }

    /**
     * Constructor
     */
    public Pair(X first, X second) {
        this.first = first;
        this.second = second;
    }

    /**
     * String representation
     */
    @Override
    public String toString() {
        return "(" + first + "--" + second + ")";
    }

    /**
     * Equality test, based on comparing components
     */
    @Override
    public boolean equals(Object o) {
        if (this == o)
            return true;
        if (o instanceof Pair) {
            @SuppressWarnings("rawtypes")
			Pair pair = (Pair) o;
            if (first != null ? !first.equals(pair.first) : pair.first != null)
                return false;
            if (second != null ? !second.equals(pair.second) : pair.second != null)
                return false;
            return true;
        }
        return false;
    }
}
