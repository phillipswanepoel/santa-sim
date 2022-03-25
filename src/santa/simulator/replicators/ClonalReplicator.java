package santa.simulator.replicators;

import santa.simulator.Virus;
import santa.simulator.fitness.FitnessFunction;
import santa.simulator.genomes.*;
import santa.simulator.mutators.Mutator;
import java.util.List;
import java.util.ArrayList;
import java.util.SortedSet;

/**
 * @author rambaut
 *         Date: Apr 27, 2005
 *         Time: 9:40:33 AM
 */
public class ClonalReplicator implements Replicator {

    public ClonalReplicator() {
        // nothing to do
    }

	public int getParentCount() {
		return 1;
	}

    public void replicate(Virus virus, Virus[] parents, Mutator mutator, FitnessFunction fitnessFunction, GenePool genePool, int generation) {

        Genome parentGenome = parents[0].getGenome();

        SortedSet<Mutation> mutations = mutator.mutate(parentGenome);

        Genome genome = genePool.duplicateGenome(parentGenome, mutations, fitnessFunction);

        virus.setGenome(genome);
        virus.setParent(parents[0]);

        //Adding indel events to virus indel list
        int pos = 0;
        int size = 0; 

        for (Mutation x: mutations) {
            boolean isDeletion = x instanceof Deletion;
            boolean isInsertion = x instanceof Insertion;

            if (isDeletion || isInsertion) {             
                pos = x.getPosition();
                size = x.length();  

                List<Integer> indels = new ArrayList<Integer>();
                indels.add(pos);
                indels.add(size);
                virus.getGenome().getSequence().addIndelEvent(indels);              
            }               
        }

    }
}
