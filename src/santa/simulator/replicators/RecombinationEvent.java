package santa.simulator.replicators;

import java.util.Arrays;
import java.util.List;
import santa.simulator.Virus;
import santa.simulator.fitness.FitnessFunction;
import santa.simulator.genomes.*;
import santa.simulator.mutators.Mutator;
import java.util.SortedSet;
import santa.simulator.replicators.RecombinantTracker;

/**
 * @author Phillip Swanepoel
 *         Date: Mar 20, 2021
 *         Time: 9:40:33 AM
 */
public class RecombinationEvent {

    public RecombinationEvent(Genome recombinant, String recombinantSequence, List<Genome> parents, List<String> parentalSequences, SortedSet<Integer> breakpoints, int gen) {
        //to get parental sequences:
        //parents.get(0).getSequence().getNucleotides());
        //parents.get(1).getSequence().getNucleotides());
        this.recombinant = recombinant;
        this.recombinantSequence = recombinantSequence;
        this.parents = parents;
        this.parentalSequences = parentalSequences;
        this.breakpoints = breakpoints;
        this.gen = gen;
    }

    public Genome getRecombinant() {
        return recombinant;
    }

    public String getRecombinantSequence() {
        return recombinantSequence;
    }  

    public List<Genome> getParents() {
        return parents;
    }

    public List<String> getParentalSequences() {
        return parentalSequences;
    }

    public SortedSet<Integer> getBreakpoints() {
        return breakpoints;
    }

    public int getGeneration() {
        return gen;
    }

    private final Genome recombinant;
    private final String recombinantSequence;
    private final List<Genome> parents;
    private final List<String> parentalSequences;
    private final SortedSet<Integer> breakpoints;
    private final int gen;
    
}
