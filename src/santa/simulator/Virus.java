
package santa.simulator;

import santa.simulator.genomes.Genome;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * @author rambaut
 *         Date: Apr 22, 2005
 *         Time: 2:23:14 PM
 */
public class Virus {

    public Virus() {
    }

    public Virus(Genome genome, Virus parent) {
        this.genome = genome;
        this.parent = parent;        
    }
 

    public Genome getGenome() {
        return genome;
    }

    public Virus getParent() {
        return parent;
    }

	public double getLogFitness() {
		return genome.getLogFitness();
	}

	public double getFitness() {
		return genome.getFitness();
	}

    public List<Integer> getRecombinationList() {
        return recombinationList;
    }


    public void setGenome(Genome genome) {
        this.genome = genome;
    }

    public void setParent(Virus parent) {
        this.parent = parent;
    }

    public void addRecombinationEvent(int recombIndex) {
        this.recombinationList.add(recombIndex);
    }

    public void setRecombinationList(List<Integer> recombinationList) {
        List<Integer> newList = new ArrayList<Integer>(recombinationList);
        this.recombinationList = newList;
    }

    public int getOffspringCount() {
        return offspringCount;
    }

    public void setOffspringCount(int offspringCount) {
        this.offspringCount = offspringCount;
    }

    private Genome genome = null;
    private Virus parent = null;
    private int offspringCount = 0;
    private List<Integer> recombinationList = new ArrayList<>();   

}
