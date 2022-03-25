/*
 * Created on Jul 17, 2006
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package santa.simulator.genomes;

import java.util.Arrays;
import java.util.SortedSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author kdforc0
 *
 * A genomic mutable sequence.
 */
public final class SimpleSequence implements Sequence {

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */

	/*
	 * Internal representation: every byte corresponds to a nucleotide.
	 */
	private byte states[];
	
	private List<List<Integer>> indelList = new ArrayList<>();
	
	/**
	 * Create a new sequence with given nucleotide length. The new sequence
	 * is initialized to all A's.
	 */
	public SimpleSequence(int length) {
		states = new byte[length];
	}

	public SimpleSequence() {
		
	}

	/**
	 * Create a new sequence with the given nucleotides.
	 */
	public SimpleSequence(String nucleotides) {
		states = new byte[nucleotides.length()];

		for (int i = 0; i < nucleotides.length(); ++i) {
			setNucleotide(i,  Nucleotide.parse(nucleotides.charAt(i)));
		}
	}


	public SimpleSequence(byte[] states) {
		// it is possible for {@code states} to be null at this point.
		if (states == null)
			states = new byte[0];
		this.states = states;
	}

	public SimpleSequence(Sequence other) {
		states = new byte[other.getLength()];

		for (List<Integer> del: other.getIndelList()) {
			this.indelList.add(del);
		}		

		copyNucleotides(0, other, 0, other.getLength());
	}

	public SimpleSequence(SimpleSequence other) {
		states = new byte[other.getLength()];
		
		for (List<Integer> del: other.getIndelList()) {
			this.indelList.add(del);
		}

		copyNucleotides(0, other, 0, other.getLength());
	}


	public SimpleSequence(Sequence other, int start, int length) {
		states = new byte[length];

		for (List<Integer> del: other.getIndelList()) {
			this.indelList.add(del);
		}			

		copyNucleotides(0, other, start, start + length);
	}

	public SimpleSequence(SimpleSequence other, int start, int length) {
		states = new byte[length];

		for (List<Integer> del: other.getIndelList()) {
			this.indelList.add(del);
		}			

		copyNucleotides(0, other, start, length);
	}

	/* Override hashCode() and equals() so Sequence instamces can be
	   used as keys in hash collections,
	   e.g. Population::initialize().
	*/
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(states);
		return result;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SimpleSequence other = (SimpleSequence) obj;
		if (!Arrays.equals(states, other.states))
			return false;
		return true;
	}

	public Sequence getSubSequence(int start, int length) {
		return new SimpleSequence(this, start, length);
	}

	/* (non-Javadoc)
		 * @see santa.simulator.genomes.Sequence#getLength()
		 */
	public int getLength() {
		return states.length;
	}

	/* (non-Javadoc)
		 * @see santa.simulator.genomes.Sequence#getAminoacidsLength()
		 */
	public int getAminoAcidsLength() {
		return getLength() / 3;
	}

	/* (non-Javadoc)
		 * @see santa.simulator.genomes.Sequence#getNucleotide(int)
		 */
	public byte getNucleotide(int i) {
		return states[i];
	}

	public void setNucleotide(int i, byte state) {
		states[i] = state;
	}

	protected void copyNucleotides(int start, SimpleSequence source,
	                               int sourceStart, int sourceCount) {
		System.arraycopy(source.states, sourceStart, states, start, sourceCount);
	}

	protected void copyNucleotides(int start, Sequence source,
	                               int sourceStart, int sourceStop) {
		for (int i = 0; i < sourceStop - sourceStart; ++i) {
			setNucleotide(start + i, source.getNucleotide(sourceStart + i));
		}
	}

	private int countChanges(int pos) {
		int cumChanges = 0;
		int tempPos = pos;

		for (List<Integer> indels : this.indelList) {

			if (indels.get(0) <= tempPos) {
				cumChanges = cumChanges + indels.get(1);
				tempPos = tempPos - indels.get(1);
			}
		}

		return cumChanges;
	}
	
	private List<List<Integer>> updateIndels(List<List<Integer>> indels, int pos, int size) {
		List<List<Integer>> updated = new ArrayList<List<Integer>>(); 
		
		for (List<Integer> indel: indels) {
			int old_pos = indel.get(0);
			if (pos <= old_pos) {
				List<Integer> newIndel = new ArrayList<Integer>();					
				newIndel.add(old_pos + size);
				newIndel.add(indel.get(1));
				newIndel.add(indel.get(2));
				newIndel.add(indel.get(3));
				
				updated.add(newIndel);
			}	
			else {
				updated.add(indel);
			}
		}
		
		return updated;
	}
	
	static List<List<Integer>> mergeOverlaps(List<List<Integer>> indels) {
		List<List<Integer>> merged = new ArrayList<List<Integer>>(indels); 
		boolean done = false;
		//if indel A goes from 1->3, aka A = [1, 3]
		//then for indel B to have a seamless overlap, needs to be:
		//B = [1, X], [2, X], [3, X] or [4, X]
		//aka if indel A = [X, Y] B = [ (X,X+1,X+2,..,X+Y+1), Z]	
		while (!done) {
			int i_counter = 0;
			int j_counter = 0;
			boolean modified = false;
			outerloop:			
			for (int i = 0; i < merged.size()-1; i++) {
				i_counter = i;
				for (int j = i+1; j < merged.size(); j++) {
					j_counter = j;
					List<Integer> indelA = merged.get(i);
					List<Integer> indelB = merged.get(j);
					
					int posA = indelA.get(0);
					int sizeA = Math.abs(indelA.get(1));
					int idA = indelA.get(3);
					//int signA = (int) Math.signum(indelA.get(1));	
					
					int posB = indelB.get(0);
					int sizeB = Math.abs(indelB.get(1));
					int idB = indelB.get(3);
					//int signB = (int) Math.signum(indelB.get(1));		
					
					if (posB >= posA && posB <= posA+sizeA) {
						//there is an overlap, can merge
						List<Integer> mergedIndel = new ArrayList<Integer>();					
						mergedIndel.add(posA);
						mergedIndel.add(sizeA + sizeB);
						mergedIndel.add(indelA.get(2));
						int modifier = (idA + idB)/2;
						mergedIndel.add(modifier);
						
						merged.remove(i);
						merged.remove(j-1);
						merged.add(mergedIndel);								
						
						modified = true;
						break outerloop;
						
					} else if (posA >= posB && posA <= posB+sizeB) {
						List<Integer> mergedIndel = new ArrayList<Integer>();					
						mergedIndel.add(posB);
						mergedIndel.add(sizeB + sizeA);
						mergedIndel.add(indelB.get(2));
						int modifier = (idB + idB)/2;
						mergedIndel.add(modifier);
						
						merged.remove(i);
						merged.remove(j-1);
						merged.add(mergedIndel);						
						
						modified = true;
						break outerloop;
					} 
				}			
			}		
			if (i_counter >= merged.size()-2 &&  j_counter >= merged.size()-1) {
				done = true;
			}						
		}
		
		
		return merged;
	}
	
	public boolean deleteSubSequence(int pos, int count) {
		// assert that insert preserved frame...
		assert (states.length % 3) == 0: "States.length not divisible by 3, pos = " + pos;
		assert (count % 3) == 0: "count not divisible by 3, pos = " + pos;

		byte newstates[] = new byte[states.length - count];
		System.arraycopy(states, 0, newstates, 0, pos);
		System.arraycopy(states, pos+count, newstates, pos, states.length-pos-count);
		states = newstates;		
		assert (states.length % 3) == 0: "States.length changed, but now not divisible by 3, pos = " + pos;	

		List<Integer> indel = new ArrayList<Integer>();	
		indel.add(pos);
		indel.add(-count);
		indel.add(countChanges(pos));
		
		//adding unique identifier to indel event
		IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
		indelcount.addCount();		
		indel.add(indelcount.getCount());

		this.indelList.add(indel);			
		
		return(true);
	}
	
	
	public boolean insertSequence(int start, SimpleSequence source) {
		// allocate more space and copy the old contents
		// assert that insert preserved frame...
		// assert (states.length % 3) == 0;
		// assert (source.getLength() % 3) == 0;

		byte newstates[] = Arrays.copyOf(states, states.length + source.getLength());
		System.arraycopy(source.states, 0, newstates, start, source.getLength());
		System.arraycopy(states, start, newstates, start+source.getLength(), states.length-start);
		states = newstates;
		// assert (states.length % 3) == 0;
		List<Integer> indel = new ArrayList<Integer>();	
		int pos = start;
		int size = (source.getLength());			
		this.indelList = updateIndels(this.indelList, pos, size);	
		
		indel.add(pos);
		indel.add(size);
		indel.add(countChanges(start));
		
		//adding unique identifier to indel event
		IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
		indelcount.addCount();		
		indel.add(indelcount.getCount());
		
		this.indelList.add(indel);				
		this.indelList = mergeOverlaps(this.indelList);	
		
		return(true);
	}

	/* (non-Javadoc)
		 * @see santa.simulator.genomes.Sequence#getAminoAcid(int)
		 */
	public byte getAminoAcid(int i) {
		int aa_i = i * 3;

		return AminoAcid.STANDARD_GENETIC_CODE[getNucleotide(aa_i)]
				[getNucleotide(aa_i + 1)]
				[getNucleotide(aa_i + 2)];
	}

	/* (non-Javadoc)
		 * @see santa.simulator.genomes.Sequence#getNucleotides()
		 */
	public String getNucleotides() {
		StringBuffer result = new StringBuffer(getLength());
		result.setLength(getLength());

		for (int i = 0; i < getLength(); ++i) {
			result.setCharAt(i, Nucleotide.asChar(getNucleotide(i)));
		}

		return result.toString();
	}

	/* (non-Javadoc)
	 * @see santa.simulator.genomes.Sequence#getNucleotideStates()
	 */
	public byte[] getNucleotideStates() {
		byte[] result = new byte[getLength()];

		for (int i = 0; i < getLength(); ++i) {
			result[i] = getNucleotide(i);
		}

		return result;
	}

	/* (non-Javadoc)
		 * @see santa.simulator.genomes.Sequence#getAminoacids()
		 */
	public String getAminoAcids() {
		StringBuffer result = new StringBuffer(getAminoAcidsLength());
		result.setLength(getAminoAcidsLength());

		for (int i = 0; i < getAminoAcidsLength(); ++i) {
			result.setCharAt(i, AminoAcid.asChar(getAminoAcid(i)));
		}

		return result.toString();
	}

	/* (non-Javadoc)
	 * @see santa.simulator.genomes.Sequence#getAminoAcidStates()
	 */
	public byte[] getAminoAcidStates() {
		byte[] result = new byte[getAminoAcidsLength()];

		for (int i = 0; i < getAminoAcidsLength(); ++i) {
			result[i] = getAminoAcid(i);
		}

		return result;
	}

	public int getLength(SequenceAlphabet alphabet) {
		return getLength() / alphabet.getTokenSize();
	}

	public byte getState(SequenceAlphabet alphabet, int i) {
		if (alphabet == SequenceAlphabet.NUCLEOTIDES)
			return getNucleotide(i);
		else
			return getAminoAcid(i);
	}

	public String getStateString(SequenceAlphabet alphabet) {
		if (alphabet == SequenceAlphabet.NUCLEOTIDES)
			return getNucleotides();
		else
			return getAminoAcids();
	}

	public byte[] getStates(SequenceAlphabet alphabet) {
		if (alphabet == SequenceAlphabet.NUCLEOTIDES)
			return getNucleotideStates();
		else
			return getAminoAcidStates();
	}


	public Sequence recombineWith(Sequence other, SortedSet<Integer> breakPoints) {
		SimpleSequence[] parents = {this, (SimpleSequence) other};
		return SimpleSequence.getRecombinantSequence(parents, breakPoints);
	}
	
	/**
	 * calculate how long the product of recombination will be.
	 *
	 **/
	static int getRecombinantLength(SimpleSequence[] parents, SortedSet<Integer> breakPoints)  {
		assert(parents.length == 2);
		assert(parents[0].getLength() <= parents[1].getLength());

		/**
		 * NOTE: for non-homologous, isolocus recombination, the
		 * length of the product genome is completely determined by
		 * the length of the parent genomes and the modulus of the
		 * number of crossings.  Homologous recombination will have a
		 * similar shortcut.
		 **/
		int fastlen = parents[breakPoints.size() % 2 == 0 ? 0 : 1].getLength();

		return fastlen;
	}
	
	/**
	 * Factory method to create a recombined nucleotide sequence from two parents.
	 *
	 * Given a pair of parent sequences and a set of breakpoints,
	 * create a new sequence that is a combination of fragments from
	 * both parents.  breakPoints describes the positions at which we
	 * switch from one template to the other.
	 *
	 * If 'breakPoints' is empty, this routine simply copies the
	 * sequence from first genome in 'parents'.
	 *
	 * NOTE: This routine implements non-homologuous recombination.
	 * It uses only a single vector to specify isoloci in each parent
	 * where recombination should occur.  Homologous recombination
	 * would need two breakpoint vectors identifying the points of
	 * homology in each parent.
	 *
	 * PHILLIP CHANGE:
	 * Now also inherits a recombined indel list
	 **/


	//class used to sort indel events by position first, then size
	//Sorts in ascending order
	static class LengthComparator implements Comparator<List<Integer>> {
        public int compare(final List<Integer> list1, final List<Integer> list2) {        
            if (list1.get(3) == list2.get(3)) {
                return Integer.compare(list1.get(3), list2.get(3));
            }
            else {
                return Integer.compare(list1.get(3), list2.get(3));
            }        
        }
    }
	
	static class LengthComparator2 implements Comparator<List<Integer>> {
        public int compare(final List<Integer> list1, final List<Integer> list2) {        
            if (list1.get(0) == list2.get(0)) {
                return Integer.compare(list1.get(0), list2.get(0));
            }
            else {
                return Integer.compare(list1.get(0), list2.get(0));
            }        
        }
    }

    static List<List<Integer>> sortIndels(List<List<Integer>> insertions) {
        List<List<Integer>> withoutOverlaps = insertions;      
        Collections.sort(withoutOverlaps, new LengthComparator());
        return withoutOverlaps;
    }
    
    static List<List<Integer>> sortIndelsByPosition(List<List<Integer>> insertions) {
        List<List<Integer>> withoutOverlaps = insertions;      
        Collections.sort(withoutOverlaps, new LengthComparator2());
        return withoutOverlaps;
    }

    //returns an updated indelList, which is created by combining the lists from both parents in a "special" way
    //the merger ensures events can still be "unpacked" in reverse order without disrupting the coordinate system
    /*
    static List<List<Integer>> getNewIndelList(SimpleSequence[] parents, SortedSet<Integer> breakPoints) {
    	List<List<Integer>> indels1 = new ArrayList<List<Integer>>(parents[0].getIndelList());
		List<List<Integer>> indels2 = new ArrayList<List<Integer>>(parents[1].getIndelList());
		List<List<Integer>> final_list = new ArrayList<List<Integer>>();
		
		return final_list;
    }
    */
    
    static List<List<Integer>> selectIndels(List<List<Integer>> oldList, List<List<Integer>> currentList, int start, int end, int shift_diff) {    	
    	List<List<Integer>> selectedIndels = new ArrayList<>();
    	
    	if (currentList.size() > 0) {
    		int counter = 0;

			List<Integer> startIndel = new ArrayList<Integer>(currentList.get(0));

			while (counter < currentList.size()) {
				
				startIndel = currentList.get(counter);					
				int old_pos = oldList.get(counter).get(0);				
				int updated_pos = startIndel.get(0);
				int size = Math.abs(startIndel.get(1));
				int sign = (int) Math.signum(startIndel.get(1));
				int unique_id = startIndel.get(3);
				
				if (start <= updated_pos && (updated_pos + size) <= end) {
					//Indel cleanly within block
					List<Integer> newIndel = new ArrayList<Integer>();
					newIndel.add(updated_pos - shift_diff);
					newIndel.add(size*sign);
					newIndel.add(startIndel.get(2) - shift_diff);					
					newIndel.add(unique_id);		
					
					selectedIndels.add(newIndel);

				} else if ( ((start <= updated_pos) && (updated_pos < end)) && ((updated_pos + size) > end) ) {
					//Indel divided by end breakpoint						
					List<Integer> dividedIndel = new ArrayList<Integer>();
					int new_size = (end - updated_pos);				
					dividedIndel.add(updated_pos - shift_diff);
					dividedIndel.add(new_size*sign);
					dividedIndel.add(startIndel.get(2) - shift_diff);
					dividedIndel.add(unique_id);	
					/*
					if (shift_diff != 0) {
						//if size has been changed, we need to treat this as if it is a new indel
						IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
						indelcount.addCount();		
						unique_id = indelcount.getCount();	
					}
					*/
					selectedIndels.add(dividedIndel);

				} else if ( (updated_pos < start) && ((updated_pos + size) > start) && ((updated_pos + size) <= end)) {
					//Indel divided by starting breakpoint						
					List<Integer> dividedIndel = new ArrayList<Integer>();
					//dividedIndel.add(old_pos + (updated_pos - start));
					dividedIndel.add(start);
					dividedIndel.add((updated_pos + size - start)*sign);
					dividedIndel.add(startIndel.get(2) - shift_diff);
					dividedIndel.add(unique_id);
					/*
					if (shift_diff != 0) {
						//if size has been changed, we need to treat this as if it is a new indel
						IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
						indelcount.addCount();		
						unique_id = indelcount.getCount();	
					}
					*/
					selectedIndels.add(dividedIndel);
				} 

				counter++;							
			}

		}
    	
    	return selectedIndels;
    }
	
	static List<List<Integer>> updatePositions(List<List<Integer>> indels) {
		//updates positions and translations so that they are accurate for current time, not time of specific indel event
		//updated = positions are updated to time of recombination		
		List<List<Integer>> updated_List = new ArrayList<List<Integer>>(); 
		
		if (indels.size() >= 2) {

			for (int i = 0; i < indels.size()-1; i++) {
				List<Integer> baseIndel = new ArrayList<Integer>(indels.get(i));
				int base_pos = baseIndel.get(0);
				int newTranslation = 0;

				for (int j = i + 1; j < indels.size(); j++) {
					List<Integer> newerIndel = new ArrayList<Integer>(indels.get(j));

					//newer indel lies before base indel, and thus base indel's position is changed at time of recombination
					if (newerIndel.get(0) <= base_pos) {
						//adding size of indel
						int size = newerIndel.get(1);
						newTranslation = newTranslation + size;
						//base_pos = base_pos + size;
					}
				}

				List<Integer> updated_indel = new ArrayList<Integer>();
				//Calculating new position and translation
				//New translation such that position+translation is equal for old and new indel		
				
				updated_indel.add(baseIndel.get(0) + newTranslation);
				updated_indel.add(baseIndel.get(1));
				updated_indel.add(baseIndel.get(2) + newTranslation);
				updated_indel.add(baseIndel.get(3));

				updated_List.add(updated_indel);
				newTranslation = 0;
			}

			//adding final event
			List<Integer> finalIndel = new ArrayList<Integer>(indels.get(indels.size()-1));			
			updated_List.add(finalIndel);			
		

		} else if (indels.size() == 1) {
			updated_List.add(indels.get(0));			
		}		
		
		return updated_List;
		
	}
	
	static int calcFrameShift(int pos, List<List<Integer>> sorted_indels, boolean hard) {
		int frameshift = 0;
		int bp_pos = pos;	
		
		for (List<Integer> indel: sorted_indels) {	
			int updated_pos = indel.get(0);
			int sign = (int) Math.signum(indel.get(1));
			int size = Math.abs(indel.get(1));		
			
			if (updated_pos < bp_pos && updated_pos >= 0) {
				if (updated_pos + size >= bp_pos) {	
					if (hard) {
						frameshift = frameshift + ((size)*sign);
						//frameshift = frameshift + (bp_pos - updated_pos)*sign;						
						bp_pos = bp_pos + (size)*sign;
					}
					else {
						frameshift = frameshift + (bp_pos - updated_pos)*sign;
						//frameshift = frameshift + ((size)*sign);
						//bp_pos = bp_pos + (size)*sign;
					}					
				} else {
					if (hard) {
						frameshift = frameshift + (size)*sign;
						bp_pos = bp_pos + (size)*sign;
					}
					else {
						frameshift = frameshift + (size)*sign;
						//bp_pos = bp_pos + (size)*sign;
					}
					
					//bp_pos = bp_pos + (size)*sign;
				}					
			}	
		}	
		
		return frameshift;
	}
	
	//calcShiftDiff(bp, indels1, indels2, updatedIndels1, updatedIndels2)
	static int calcHomoBreakPoint(int pos, List<List<Integer>> indels1_sorted, List<List<Integer>> indels2_sorted) {
		//calculates difference in frameshift between two sequences (caused by indels) up to a certain nucleotide position	
		System.out.println("---------------------------------------");		
		System.out.println(indels1_sorted);
		System.out.println(indels2_sorted);	
		
		int frameshift1 = calcFrameShift(pos, indels1_sorted, false);
		int frameshift2 = calcFrameShift(Math.max(pos - frameshift1, 0), indels2_sorted, true);		
				
		return pos + frameshift2 - frameshift1;
	}
	
	
	static private class indelCollection {
		private List<List<Integer>> indels1;
		private List<List<Integer>> indels2;
		private List<List<Integer>> updatedIndels1;
		private List<List<Integer>> updatedIndels2;
		private List<List<Integer>> sortedIndels1;
		private List<List<Integer>> sortedIndels2;
		private List<Integer> breakPointsList;
		private List<Integer> homologousBreakPointsList;
		int[] parentLengths;
		
		public indelCollection(List<List<Integer>> parent1indels, List<List<Integer>> parent2indels, SortedSet<Integer> breakPoints, int[] parentLengths) {	
			this.indels1 = new ArrayList<List<Integer>>(parent1indels);		
			this.indels2 = new ArrayList<List<Integer>>(parent2indels);	
			//this.updatedIndels1 =  new ArrayList<List<Integer>>(updatePositions(indels1));
			//this.updatedIndels2 =  new ArrayList<List<Integer>>(updatePositions(indels2));	
			this.updatedIndels1 = new ArrayList<List<Integer>>((indels1));
			this.updatedIndels2 = new ArrayList<List<Integer>>((indels2));			
		
			this.sortedIndels1 = sortUpdated(this.updatedIndels1);
			this.sortedIndels2 = sortUpdated(this.updatedIndels2);
		
			this.breakPointsList = new ArrayList<Integer>(breakPoints);
			this.parentLengths = parentLengths;
		}			
	
		public List<List<Integer>> getIndels1() {
			return indels1;
		}
		public List<List<Integer>> getIndels2() {
			return indels1;
		}
		public List<List<Integer>> getUpdatedIndels1() {
			return updatedIndels1;
		}
		public List<List<Integer>> getUpdatedIndels2() {
			return updatedIndels2;
		}
		public List<Integer> getBreakPoints() {
			return breakPointsList;
		}
		
		private List<List<Integer>> sortUpdated(List<List<Integer>> indels) {
			List<List<Integer>> sorted_indels = new ArrayList<List<Integer>>(indels); 		
			sorted_indels = sortIndelsByPosition(sorted_indels);
			
			return sorted_indels;
		}
		
		private int refineBP(int bp, int homobp) {
			//need to be smarter about how indels are dealth with that go across the breakpoint
			//if this indel is shared between both parents, need to make homo breakpoint such that it is preserved
			int refined_bp = homobp;	
			int nucleotides_inherited = 0;			
			int pos = -1;
			int og_size = 0;
			for (List<Integer> indel : updatedIndels1) {
				int updated_pos = indel.get(0);
				int sign = (int) Math.signum(indel.get(1));
				int size = Math.abs(indel.get(1));	
				
				if (updated_pos < bp && updated_pos+size >= bp && sign > 0) {					
					nucleotides_inherited = (bp - updated_pos);
					og_size = size*sign;
					pos = updated_pos;
				}
			}
			//Now check if this indel is shared with other parent:
			for (List<Integer> indel : updatedIndels2) {
				int updated_pos = indel.get(0);						
				int size = Math.abs(indel.get(1));
				int sign = (int) Math.signum(indel.get(1));
				if (updated_pos == pos || pos != -1 && updated_pos <= homobp && updated_pos+size > homobp && sign > 0) {
					//we have a match, update homobp accordingly
					if (nucleotides_inherited <= updated_pos+size-homobp) {
						refined_bp = refined_bp + nucleotides_inherited;
					}					
				}
			}
			
			return refined_bp;
		}
		
		public List<Integer> getHomologousBreakPoints() {
			List<Integer> homologousBreakPoints = new ArrayList<Integer>();
			
			for (int bp :breakPointsList) {
				
				int homo_bp = calcHomoBreakPoint(bp, sortedIndels1, sortedIndels2);				
				homo_bp = refineBP(bp, homo_bp);				
				homo_bp = Math.min(homo_bp, this.parentLengths[1]);		
				if (homo_bp < 0) {
					homo_bp = 0;
				}
				System.out.println("Refined homo bp: " + (homo_bp));
				homologousBreakPoints.add(homo_bp);				
			}		
			
			this.homologousBreakPointsList = homologousBreakPoints;
			return homologousBreakPoints;
		}
		
		public List<List<Integer>> getNewIndelList() {		
			List<List<Integer>> updated_indels1 = new ArrayList<List<Integer>>(updatedIndels1);
			List<List<Integer>> updated_indels2 =  new ArrayList<List<Integer>>(updatedIndels2);
			//Now use the updated lists to merge them according to breakpoints
			//Need to add consideration for end of sequence as breakpoint
			List<Integer> breakPoints_with_end = new ArrayList<Integer>(breakPointsList);
			List<Integer> homoBreakPoints_with_end = new ArrayList<Integer>(homologousBreakPointsList);
			
			int recombinant_len = parentLengths[1] - (homologousBreakPointsList.get(0) - breakPointsList.get(0));			
			//Length of recombinant is as follows:
			//parents are sorted such that parent1 < parent2
			//if only 1 breakpoint, len == parent2
			//if 2 breakpoints, len == parent1
			if (breakPointsList.size() == 1) {
				breakPoints_with_end.add(recombinant_len);
				homoBreakPoints_with_end.add(recombinant_len);
			} else {
				breakPoints_with_end.add(this.parentLengths[0]);
				homoBreakPoints_with_end.add(this.parentLengths[0]);
			}				
			//Looping through breakpoints and adding indels that fall between consecutive ones	
			//selectIndels(List<List<Integer>> oldList, List<List<Integer>> currentList, int start, int end, int shift_diff)		
			updated_indels1 =  selectIndels(indels1, updated_indels1, 0, breakPoints_with_end.get(0), 0);
			updated_indels2 =  selectIndels(indels2, updated_indels2, homoBreakPoints_with_end.get(0),  parentLengths[1],  homoBreakPoints_with_end.get(0) - breakPoints_with_end.get(0));
			
			List<List<Integer>> merged_indels =  new ArrayList<List<Integer>>();	
			merged_indels.addAll(updated_indels1);
			merged_indels.addAll(updated_indels2);
			
			merged_indels = sortIndelsByPosition(merged_indels);
			
			return merged_indels;

		}
		
	}
	

    static SimpleSequence getRecombinantSequence(SimpleSequence[] parents, SortedSet<Integer> breakPoints) {
	 	assert(parents.length == 2);
		assert(parents[0].getLength() <= parents[1].getLength());
				
		List<List<Integer>> indels1 = new ArrayList<List<Integer>>(parents[0].getIndelList());
		List<List<Integer>> indels2 = new ArrayList<List<Integer>>(parents[1].getIndelList());
		indelCollection indels = new indelCollection(indels1, indels2, breakPoints, new int[]{parents[0].getLength(), parents[1].getLength()});
		
		List<Integer> listBreakPoints = new ArrayList<Integer>(indels.getBreakPoints());	
		List<Integer> homologousBreakPoints = new ArrayList<Integer>(indels.getHomologousBreakPoints());
		//List<Integer> homologousBreakPoints = 
		//homologousBreakPoints = calcHomologousBreakpoints();
		
		int recombinant_len = parents[1].getLength() - (homologousBreakPoints.get(0) - listBreakPoints.get(0));

		int lastBreakPoint = 0;		// previous recombination location
		int currentSeq = 0;			// index of currently selected parent
		//int newlen = getRecombinantLength(parents, breakPoints);
		SimpleSequence product = new SimpleSequence(recombinant_len);

		byte[] dest = product.states;	// where to put the product
		SimpleSequence seq = parents[currentSeq];
		
		String parent1 = parents[0].getNucleotides();
		String parent2 = parents[1].getNucleotides();				
	
		System.out.println("Breakpoint = " + listBreakPoints.get(0));
		System.out.println(parent1.substring(0, listBreakPoints.get(0)) + "|" + parent1.substring(listBreakPoints.get(0)));
		System.out.println(parent2.substring(0, homologousBreakPoints.get(0)) + "|" + parent2.substring(homologousBreakPoints.get(0)));	
		
		int counter = 0;
		for (int nextBreakPoint : breakPoints) {				
			int homologousNextBreakPoint = homologousBreakPoints.get(counter);	
			
			System.arraycopy(seq.states, lastBreakPoint, 
							 dest, lastBreakPoint, nextBreakPoint-lastBreakPoint);			
			
			lastBreakPoint = homologousNextBreakPoint;
			homologousBreakPoints.add(lastBreakPoint);
			currentSeq = 1 - currentSeq;
			seq = parents[currentSeq];	
			counter++;
		}
		
		
		
		System.arraycopy(seq.states, lastBreakPoint, 
						 dest, breakPoints.last(), parents[1].getLength()-lastBreakPoint);
		
		String recombinant = product.getNucleotides();		
		System.out.println(recombinant.substring(0, breakPoints.last()) + "|" + recombinant.substring(breakPoints.last()));		
		
		List<List<Integer>> newIndels = indels.getNewIndelList();
		product.setIndelList(newIndels);
		
		System.out.println("======");
		System.out.println(newIndels);		
		
		return(product);
	}

	public List<List<Integer>> getIndelList() {
        return indelList;
    }

    public void setIndelList(List<List<Integer>> indels) {
    	List<List<Integer>> newList = new ArrayList<List<Integer>>(indels);
        this.indelList = newList;
    }

    public void addIndelEvent(List<Integer> indelEvent) {
        this.indelList.add(indelEvent);
    }

    public void removeIndelEvent(List<Integer> indelEvent) {
    	List<List<Integer>> newIndelList = new ArrayList<List<Integer>>(this.indelList);
        newIndelList.remove(indelEvent);
        this.indelList = new ArrayList<List<Integer>>(newIndelList);
    }
}
