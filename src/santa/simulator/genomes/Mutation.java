
package santa.simulator.genomes;

import java.util.List;
import java.util.ArrayList;

import org.apache.commons.lang3.Range;


/**
 * @author rambaut
 *         Date: Apr 26, 2005
 *         Time: 10:28:35 AM
 *         
 */
public class Mutation  implements Comparable<Mutation>  {
    protected Mutation() { state = 0; }

    private Mutation(int position, byte state) {
        this.position = position;
        this.state = state;
    }

    /*
     * position is 0-based -- koen.
     */
    public int position;
    public final byte state;

    public static Mutation getMutation(int position, byte state) {
        return new Mutation(position, state);
    }

	public boolean apply(Genome genome) {
		return(genome.substitute(position, state));
	}
	
	public Range<Integer> apply(Range<Integer> r) { return(r); }

    public int compareTo(Mutation other) {
        return other.position - position;
    }

    public boolean equals(Object other) {
        return ((Mutation) other).position == position;
    }


	public int length() {
		return 1;
	}

	public int getPosition() {
		return position;
	}
	
	/**
	 * create a list of nucleotides changed by this mutation.
	 */
	public List<StateChange> getChanges(Genome genome, int[] featureSiteTable) {
		
		List<StateChange> l = new ArrayList<StateChange>();
		
		try {
			if (featureSiteTable[this.position] != -1) {
				byte oldState = genome.getNucleotide(this.position);
				StateChange c = new StateChange(featureSiteTable[this.position], oldState, this.state);
				l.add(c);
			}
		} catch (java.lang.ArrayIndexOutOfBoundsException e) {
			System.out.println("WARNING! The expected index error has happened: ");
			System.out.println(this.position);
			System.out.println(genome.getLength());	
			System.out.println(featureSiteTable[this.position]);				
			System.exit(0);
		}		
		
		return (l);
	}

}
