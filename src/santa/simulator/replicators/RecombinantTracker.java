package santa.simulator.replicators;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import santa.simulator.Virus;
import santa.simulator.fitness.FitnessFunction;
import santa.simulator.genomes.*;
import santa.simulator.mutators.Mutator;
import santa.simulator.replicators.RecombinationEvent;

import java.util.SortedSet;

/**
 * @author rambaut
 *         Date: Apr 27, 2005
 *         Time: 9:40:33 AM
 */
public class RecombinantTracker {
    
    public RecombinantTracker() {

        this.recombinationCount = recombinationCount;
        this.recombinationList = recombinationList;
    }  

    public int getRecombinationCount() {
        return recombinationCount;
    }

    public List<RecombinationEvent> getrecombinationList() {
        return recombinationList;
    }    
    
    private static RecombinantTracker instance = null;
    public static int recombinationCount = 0;
    public static List<RecombinationEvent> recombinationList = new ArrayList<>();

    public static synchronized RecombinantTracker getInstance() {
    if (instance == null) {
        instance = new RecombinantTracker();
    }
    return instance;

    }
}
