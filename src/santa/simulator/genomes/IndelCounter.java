package santa.simulator.genomes;

/**
 * @file   Indel.java
 * @author cswarth
 * @date   Thu Jun 11 16:12:43 PDT 2015
 * 
 * 
 */


/**
 * Counts indel events as they happen, allowing each indel to have a unique identifier, useful
 * later know which viruses share which indel events
 */

public enum IndelCounter  {
	INSTANCE(0);
	
	private int counter;
	
	private IndelCounter(int startingValue) {
		this.counter = startingValue;
	}
	
	public IndelCounter getInstance() {
		return INSTANCE;
	}
	
	public int getCount() {
		return counter;
	}
	
	public void addCount() {
		this.counter++;
	}
	
}
