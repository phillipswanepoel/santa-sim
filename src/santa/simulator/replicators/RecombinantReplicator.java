package santa.simulator.replicators;

import java.util.Arrays;
import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.LinkedHashSet; 
import java.util.List;
import java.util.ArrayList;
import java.util.stream.Collectors;
import java.util.logging.*;

import org.apache.commons.math3.distribution.BinomialDistribution;

import santa.simulator.EventLogger;
import santa.simulator.Random;
import santa.simulator.Virus;
import santa.simulator.SimulationEpoch;
import santa.simulator.fitness.FitnessFunction;
import santa.simulator.genomes.GenePool;
import santa.simulator.genomes.Genome;
import santa.simulator.genomes.GenomeDescription;
import santa.simulator.genomes.Mutation;
import santa.simulator.genomes.Insertion;
import santa.simulator.genomes.Deletion;
import santa.simulator.genomes.Sequence;
import santa.simulator.mutators.Mutator;
import santa.simulator.replicators.RecombinationEvent;
import santa.simulator.samplers.SamplingSchedule;

/**
 * @author rambaut
 *         Date: Apr 27, 2005
 *         Time: 9:40:33 AM
 */

/**
So you would need to insert a logger here. You would need to create the logger itself also. 

In the class mentioned above the variable breakPoints contains the basepair position of breakpoints. 

The variable parents contains the two genomes that will recombine. contains the parents. 
So you should get the identifier of these parents and insert it into your logger, alongside the current generation. I think there are functions for this purpose already made. 

Similarly you should log the sequence id that is created after recombination. This would be the genome variable, which is created using createGenome function. 
Recombination occurs before this function is called, and random mutations are added afterwards usingng the applyMutations fucntion. So your logger should be there after that.
*/

public class RecombinantReplicator implements Replicator {

    public RecombinantReplicator(double dualInfectionProbability, double recombinationProbability) {
        this.dualInfectionProbability = dualInfectionProbability;
        this.recombinationProbability = recombinationProbability;
    }

	public int getParentCount() {
		return 2;
	}

    public void replicate(Virus virus, Virus[] vparents, Mutator mutator, FitnessFunction fitnessFunction, GenePool genePool, int generation) {
		Logger logger = Logger.getLogger("santa.simulator.replicators");


        if (Random.nextUniform(0.0, 1.0) < dualInfectionProbability * recombinationProbability) {        	    	

            // dual infection and recombination
			List<Genome> parents = Arrays.asList(vparents).stream().map(Virus::getGenome).collect(Collectors.toList());

			//Storing parental sequences to store in recombination event
			String parent1seq = parents.get(0).getSequence().getNucleotides();
			String parent2seq = parents.get(1).getSequence().getNucleotides();
			List<String> parentSeqs = Arrays.asList(parent1seq, parent2seq);
			

			// sort the parents by increasing genome length
			parents.sort((p1, p2) -> p1.getLength() - p2.getLength());

			// get minimum length of the parents
			int length = parents.stream().map(g -> g.getLength()).reduce(Integer::min).get() - 1 ;

			// pick number of breakpoints
			//BinomialDistribution binomialDeviate = new BinomialDistribution(Random.randomData.getRandomGenerator(), length, recombinationProbability);
			//int nbreaks = binomialDeviate.sample() + 1;
			int nbreaks;
			if(Math.random() < 0.5) {
				nbreaks = 1;
			}
			else {
				nbreaks = 2;
			}
			
			// Then draw the positions.
			// Don't repeat a breakpoint, and only break at codon boundaries.
			// Change: to ensure a breakpoint is actually selected:
			// dont discard if bp%3!=0, rather add so that bp%3==0 always.
			SortedSet<Integer> breakPoints = new TreeSet<Integer>();
			for (int i = 0; i < nbreaks; i++) {
				int bp = Random.nextInt(1, length-3);
				if (bp % 3 == 0) {
					breakPoints.add(bp);
				}
				else {
					bp = bp + (3-(bp % 3));
					breakPoints.add(bp);
				}

			}

			logger.finest("recombination: " + breakPoints.size() + "@" + breakPoints);			

			// create the recombinant genome description
			GenomeDescription[] gd_parents = parents.stream().map(Genome::getDescription).toArray(GenomeDescription[]::new);

			GenomeDescription recombinantGenome = GenomeDescription.recombine(gd_parents, breakPoints);
			
			Sequence recombinantSequence = getRecombinantSequence(parents, breakPoints);	
			
			Genome genome = genePool.createGenome(recombinantSequence, recombinantGenome);
			
	        SortedSet<Mutation> mutations = mutator.mutate(genome);

	        genome.setFrequency(1);

	        genome.applyMutations(mutations);

	        // we can't just update some of the fitness so recompute...
	        fitnessFunction.computeLogFitness(genome);

            virus.setGenome(genome);
            virus.setParent(vparents[0]);                

            //Loggin stuff			
	        String fitnessStr = parents.stream().map(Genome::getLogFitness).map(Object::toString).collect(Collectors.joining(", "));
	        EventLogger.log("Recombination: (" + fitnessStr + ") -> " + genome.getLogFitness());

	        //Only sample if in last 80% of generations
            if (generation >= ((int) SamplingSchedule.getSampFreq()*0.20))
            {
            	
            	//Getting recombination event list from both parents and concatenating them for recombinant
	            LinkedHashSet<Integer> recombinations0 = vparents[0].getRecombinationList();
	            LinkedHashSet<Integer> recombinations1 = vparents[1].getRecombinationList();
	            LinkedHashSet<Integer> newRecombinationList = new LinkedHashSet<Integer>();


	            newRecombinationList.addAll(recombinations0);
	            newRecombinationList.addAll(recombinations1);
	            //System.out.println("^^^");
	            //System.out.println(newRecombinationList);
	            virus.setRecombinationList(newRecombinationList); 

	            //Creating recombination event to store recombination information
				RecombinationEvent rec = new RecombinationEvent(genome, genome.getSequence().getNucleotides(), parents, parentSeqs, breakPoints, generation);	
				
				//System.out.println("********************************");
				int len = RecombinantTracker.recombinationList.size();
				//System.out.println(len);

				if (len > 0) {
					
					RecombinationEvent ev = RecombinantTracker.recombinationList.get(len-1);
					//System.out.println(ev.getRecombinant());
					//System.out.println(ev.getParents());
					//System.out.println(ev.getBreakpoints());
					//System.out.println(nbreaks);
					//System.out.println(ev.getGeneration());
					//System.out.println("********************************");
				}
				

				RecombinantTracker.recombinationList.add(rec);
				virus.addRecombinationEvent(len);			

				//System.out.println(virus.getRecombinationList());
				//System.out.println("");
				//System.out.println("");
			}

        } else {

            // single infection - no recombination...
            Genome parentGenome = vparents[0].getGenome();
                     
            LinkedHashSet<Integer> recombinations = vparents[0].getRecombinationList();          

            SortedSet<Mutation> mutations = mutator.mutate(parentGenome);
              
            Genome genome = genePool.duplicateGenome(parentGenome, mutations, fitnessFunction);

            virus.setGenome(genome);
            virus.setParent(vparents[0]);
            virus.setRecombinationList(recombinations); 
            
        }

    }


	/**
	 * Create a recombined nucleotide sequence from two parents.
	 *
	 * Given a pair of parent genomes and a set of breakpoints, create
	 * a new sequence that is a combination of fragments from both
	 * parents.  breakPoints describes the positions at which we
	 * switch from one template to the other.  

	 * 'len' is the length of the recombined sequence.  If
	 * 'breakPoints' is empty, this routine simply copies the sequence
	 * from first genome in 'parents'.
	 **/
    private static Sequence getRecombinantSequence(List<Genome> parents, SortedSet<Integer> breakPoints) {
		return parents.stream().map(Genome::getSequence).reduce((s1, s2) -> s1.recombineWith(s2, breakPoints)).get();
	}


    private final double dualInfectionProbability;
    private final double recombinationProbability;
}
