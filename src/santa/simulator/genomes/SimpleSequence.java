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
			if ((del.get(0) > start) && (del.get(0) < start+length)) {
				this.indelList.add(del);
			}			
		}			

		copyNucleotides(0, other, start, start + length);
	}

	public SimpleSequence(SimpleSequence other, int start, int length) {
		states = new byte[length];

		for (List<Integer> del: other.getIndelList()) {
			if ((del.get(0) > start) && (del.get(0) < start+length)) {
				this.indelList.add(del);
			}	
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
		
		List<List<Integer>> indels_reversed = new ArrayList<List<Integer>>(this.indelList);		
		Collections.reverse(indels_reversed);	
		
		for (List<Integer> indels : indels_reversed) {
			
			if (indels.get(1) > 0) {
				if (indels.get(0) + indels.get(1) <= tempPos) {
					cumChanges = cumChanges + indels.get(1);
					tempPos = tempPos - indels.get(1);
				} else if (indels.get(0) < tempPos && indels.get(0) + indels.get(1) > tempPos) {
					int change = tempPos - indels.get(0);
					cumChanges = cumChanges + change;
					tempPos = tempPos - change;				
				}
			} else {
				if (indels.get(0) <= tempPos) {
					cumChanges = cumChanges + indels.get(1);
					
				} else {
					
				}
			}
			
		}

		return cumChanges;
	}
	
	private List<Integer> makeIndel(int p, int s, int t, int i) {
		List<Integer> newIndel = new ArrayList<Integer>();					
		newIndel.add(p);
		newIndel.add(s);
		newIndel.add(t);
		newIndel.add(i);
		
		return newIndel;
	}
	
	private List<List<Integer>> updateIndels(List<List<Integer>> indels, int pos, int size, List<Integer> new_indel) {
		List<List<Integer>> updated = new ArrayList<List<Integer>>(); 		
		int size_mod = 0;
		
		for (List<Integer> indel: indels) {
			int old_pos = indel.get(0);
			int old_size = indel.get(1);
			if (size > 0) {
				//if new indel is an insertion
				if (pos <= old_pos) {		
					updated.add(makeIndel(old_pos + size, indel.get(1), indel.get(2) + size, indel.get(3)));
					
				} else if (pos > old_pos && pos < old_pos + old_size) {
					//indel has been split into two parts by newer indel
					updated.add(makeIndel(old_pos, pos - old_pos, indel.get(2), indel.get(3)));		
					updated.add(makeIndel(pos + size, indel.get(1) - (pos - old_pos), indel.get(2) + size + (pos - old_pos), indel.get(3)));	
				}
				else {
					updated.add(indel);
				}
			} else {		
				//new indel is a deletion
				if (old_size < 0) {
					//if old indel is a deletion
					if (pos == old_pos) {
						//combine deletions into one larger one
						size_mod += old_size;
						
					} else if (pos - size <= old_pos) {
						//if before, just shift position	
						if (pos - size == old_pos) {
							size_mod += old_size;
						} else {
							updated.add(makeIndel(old_pos + size, indel.get(1), indel.get(2) + size, indel.get(3)));
						}						
						
					} else if (pos < old_pos && pos - size > old_pos) {
						//new deletion overlaps old one, no need to shift
						size_mod += old_size;
						
					} else {
						//if not, can just add, since deletions can't overlap with other deletions in any way
						updated.add(indel);
					}
				} else {
					//if old indel is an insertion, overlaps can happen, deleting parts of the insertion
					if (pos - size <= old_pos) {
						//no overlap just shift position						
						updated.add(makeIndel(old_pos + size, indel.get(1), indel.get(2) + size, indel.get(3)));
						
					} else if (pos >= old_pos + old_size) {
						//no overlap, no position shift needed
						updated.add(makeIndel(old_pos , old_size, indel.get(2), indel.get(3)));
						
					} else if (pos <= old_pos && pos - size > old_pos ) {
						//overlap, deletion overlaps with start of insertion
						int size_remaining = old_pos + old_size - (pos - size);						
						
						if (size_remaining > 0) {
							//if some of the insertion has survived, need to add this remaining fragment							
							updated.add(makeIndel(pos, size_remaining, indel.get(2) + (pos - old_pos), indel.get(3)));		
							size_mod += old_size - size_remaining;
						} else {							
							size_mod += old_size;
						}
						
						
					} else if (pos > old_pos && pos < old_pos + old_size) {
						//overlap, deletion is "inside" of insertion, could either delete from end, or be completely inside and split insertion into two parts
						if (pos - size >= old_pos + old_size) {
							//indel not split, end fragment is just deleted
							updated.add(makeIndel(old_pos, pos - old_pos, indel.get(2), indel.get(3)));
							size_mod += old_pos + old_size - pos;
						} else {
							//indel has been split into two parts
							updated.add(makeIndel(old_pos, old_size + size, indel.get(2), indel.get(3)));	
							//possible add second indel here?
							size_mod = -size;
						}
					} 
				}
			}
			
		}
						
		if ( size + size_mod != 0) {
			updated.add(makeIndel(pos, size + size_mod, new_indel.get(2), new_indel.get(3)));		
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
					int signA = (int) Math.signum(indelA.get(1));	
					
					int posB = indelB.get(0);
					int sizeB = Math.abs(indelB.get(1));
					int idB = indelB.get(3);
					int signB = (int) Math.signum(indelB.get(1));		
					//if (idA == idB && posB >= posA && posB <= posA+sizeA) 
					if (posA == posB && signA == signB && signA < 0) {
						//two deletions can be made one
						List<Integer> mergedIndel = new ArrayList<Integer>();					
						mergedIndel.add(posA);
						mergedIndel.add((sizeA + sizeB)*signA);
						mergedIndel.add(indelA.get(2));
						int modifier = (idA + idB)/2;
						mergedIndel.add(modifier);
						
						merged.remove(i);
						merged.remove(j-1);
						merged.add(mergedIndel);								
						
						modified = true;
						break outerloop;							
						
						
					} else if (idA == idB && posB >= posA && posB <= posA+sizeA && signA == signB && signA > 0) {
						//there is an overlap, can merge
						List<Integer> mergedIndel = new ArrayList<Integer>();					
						mergedIndel.add(posA);
						mergedIndel.add((sizeA + sizeB)*signA);
						mergedIndel.add(indelA.get(2));
						int modifier = (idA + idB)/2;
						mergedIndel.add(modifier);
						
						merged.remove(i);
						merged.remove(j-1);
						merged.add(mergedIndel);								
						
						modified = true;
						break outerloop;
						
					//else if (idA == idB && posA >= posB && posA <= posB+sizeB)
					} else if (idA == idB && posA >= posB && posA <= posB+sizeB && signA == signB && signA > 0) {
						List<Integer> mergedIndel = new ArrayList<Integer>();					
						mergedIndel.add(posB);
						mergedIndel.add((sizeB + sizeA)*signA);
						mergedIndel.add(indelB.get(2));
						int modifier = (idB + idB)/2;
						mergedIndel.add(modifier);
						
						merged.remove(i);
						merged.remove(j-1);
						merged.add(mergedIndel);						
						
						modified = true;
						break outerloop;
						
					} else if (idA == idB && posA == posB && sizeA == sizeB && signA == signB && signA < 0) {
						List<Integer> mergedIndel = new ArrayList<Integer>();					
						mergedIndel.add(posA);
						mergedIndel.add((sizeA)*signA);
						mergedIndel.add(indelB.get(2));						
						mergedIndel.add(indelB.get(3));
						
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
		
		this.indelList = sortIndelsByPosition(this.indelList);		
		List<Integer> indel = new ArrayList<Integer>();	
		indel.add(pos);
		indel.add(-count);
		int translation = countChanges(pos);
		indel.add(translation);
		
		//adding unique identifier to indel event
		IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
		indelcount.addCount();		
		indel.add(indelcount.getCount());
	
		this.indelList = updateIndels(this.indelList, pos, -count, indel);	
				
		if (pos - translation < 0) {
			System.out.println("OH POESSSSSS");
			System.out.println(getNucleotides());
			System.out.println(this.indelList);
			System.out.println(indel);
			System.exit(0);
		}	
			
		return(true);
	}
	
	public int countNetChange(List<List<Integer>> indels) {
		int sum = 0;
		for (List<Integer> indel: indels) {
			sum += indel.get(1);
		}
		
		return sum;
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
			
		this.indelList = sortIndelsByPosition(this.indelList);
		
		indel.add(pos);
		indel.add(size);
		int translation = countChanges(pos);
		indel.add(translation);		
		
		
		//adding unique identifier to indel event
		IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
		indelcount.addCount();		
		indel.add(indelcount.getCount());			
			
		this.indelList = updateIndels(this.indelList, pos, size, indel);		
			
		//this.indelList.add(indel);	
		//this.indelList = mergeOverlaps(this.indelList);	
		
		if (pos - translation < 0) {
			System.out.println("OH POES");
			System.out.println(getNucleotides());
			System.out.println(this.indelList);
			System.out.println(indel);
			System.exit(0);
		}		
		
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
                return Integer.compare(list1.get(1), list2.get(1));
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
    
    static int countTotalChange(List<List<Integer>> indels) {
    	int total = 0;
    	for (List<Integer> indel: indels) {
    		total += indel.get(1);
    	}
    	return total;
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
	

	

	
	
	static private class indelCollection {
		private List<List<Integer>> indels1;
		private List<List<Integer>> indels2;
		private List<List<Integer>> updatedIndels1;
		private List<List<Integer>> updatedIndels2;
		private List<List<Integer>> sortedIndels1;
		private List<List<Integer>> sortedIndels2;
		private List<Integer> breakPointsList;
		private List<Integer> unrefinedHomologousBreaks;
		private List<Integer> homologousBreakPointsList;
		int[] parentLengths;		
		private int del_shift;
		
		public indelCollection(List<List<Integer>> parent1indels, List<List<Integer>> parent2indels, SortedSet<Integer> breakPoints, int[] parentLengths) {	
			this.del_shift = 1;
			this.unrefinedHomologousBreaks = new ArrayList<Integer>();
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
		
		private List<List<Integer>> selectIndels(List<List<Integer>> oldList, List<List<Integer>> currentList, int start, int end, int shift_diff, int type) {    	
	    	List<List<Integer>> selectedIndels = new ArrayList<>();    
	    	
	    	if (currentList.size() > 0) {
	    		int counter = 0;

				List<Integer> startIndel = new ArrayList<Integer>(currentList.get(0));
				List<Integer> lastIndel = new ArrayList<Integer>(startIndel);

				while (counter < currentList.size()) {
					
					startIndel = currentList.get(counter);					
					int old_pos = oldList.get(counter).get(0);				
					int updated_pos = startIndel.get(0);
					int size = Math.abs(startIndel.get(1));
					int sign = (int) Math.signum(startIndel.get(1));
					int unique_id = startIndel.get(3);		
					int unrefined_bp = this.unrefinedHomologousBreaks.get(0);	
					
					
					if ((start <= updated_pos || (sign < 0 && (type== 1 || type == 2) && updated_pos == unrefined_bp)) && ((sign > 0 && updated_pos + size <= end) || 
													(sign < 0 && ((updated_pos < end && (type==0 || type==2)) || (updated_pos <= end && (type==1 || type ==3))))))  {
						//Indel cleanly within block
						if (sign > 0 || updated_pos > start) {
							List<Integer> newIndel = new ArrayList<Integer>();
							newIndel.add(updated_pos - shift_diff);
							newIndel.add(size*sign);
							newIndel.add(startIndel.get(2) - shift_diff);					
							newIndel.add(unique_id);	
							newIndel.add(-3);	
							
							selectedIndels.add(newIndel);
						} else if (sign < 0 && (updated_pos == start || updated_pos == unrefined_bp)) {							
							
							if (del_shift < 1 && (type==1 || type == 2)) {
								if (del_shift < 0) {
									List<Integer> newIndel = new ArrayList<Integer>();
									newIndel.add(updated_pos - shift_diff );
									newIndel.add(del_shift);
									newIndel.add(startIndel.get(2) - shift_diff);					
									newIndel.add(unique_id);	
									newIndel.add(-4);
									
									selectedIndels.add(newIndel);
								}								
							} else {
								List<Integer> newIndel = new ArrayList<Integer>();
								newIndel.add(updated_pos - shift_diff);
								newIndel.add(size*sign);
								newIndel.add(startIndel.get(2) - shift_diff);					
								newIndel.add(unique_id);	
								newIndel.add(-3);	
								
								selectedIndels.add(newIndel);
							}
						}
						

					} else if ( ((start <= updated_pos) && (updated_pos < end)) && ((updated_pos + size) > end)) {
						//Indel divided by end breakpoint						
						List<Integer> dividedIndel = new ArrayList<Integer>();
						int new_size = (end - updated_pos);				
						dividedIndel.add(updated_pos - shift_diff);
						dividedIndel.add(new_size*sign);
						dividedIndel.add(startIndel.get(2) - shift_diff);
						dividedIndel.add(unique_id);	
						dividedIndel.add(-2);	
						/*
						IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
						indelcount.addCount();		
						unique_id = indelcount.getCount();	
						dividedIndel.add(unique_id);						
						*/					
						selectedIndels.add(dividedIndel);

					} else if ((updated_pos < start) && ((updated_pos + size*sign) > start) && ((updated_pos + size) <= end)) {
						//Indel divided by starting breakpoint						
						List<Integer> dividedIndel = new ArrayList<Integer>();
						//dividedIndel.add(old_pos + (updated_pos - start));					
						int change = start - updated_pos;
						dividedIndel.add(updated_pos - shift_diff + change);
						dividedIndel.add((updated_pos + size - start)*sign);
						dividedIndel.add(startIndel.get(2) - shift_diff + change);	
						dividedIndel.add(unique_id);	
						
//						IndelCounter indelcount = IndelCounter.INSTANCE.getInstance();
//						indelcount.addCount();		
//						unique_id = indelcount.getCount();	
//						dividedIndel.add(unique_id);
						
						dividedIndel.add(-1);	
						selectedIndels.add(dividedIndel);
					} 
					
					lastIndel = new ArrayList<Integer>(startIndel);
					counter++;					
							    	
				}

			}   
	    	
	    	return selectedIndels;
	    }
		
		private int calcFrameShift(int pos, List<List<Integer>> sorted_indels, boolean hard) {
			
			int frameshift = 0;
			int bp_pos = pos;	
			
			
			for (List<Integer> indel: sorted_indels) {	
				int updated_pos = indel.get(0);
				int sign = (int) Math.signum(indel.get(1));
				int size = Math.abs(indel.get(1));		
				
				if (updated_pos < bp_pos && updated_pos >= 0) {				
					if (hard) {
						if (sign < 0) {
							if (updated_pos + size < bp_pos) {
								frameshift = frameshift + size*sign;											
								bp_pos = bp_pos + size*sign;
							} else {
								this.del_shift = (updated_pos + size - bp_pos)*sign;								
								frameshift = frameshift + (bp_pos - updated_pos)*sign;		
								bp_pos = bp_pos + (bp_pos - updated_pos)*sign;									
							}								
							
						} else {
							frameshift = frameshift + size*sign;											
							bp_pos = bp_pos + size*sign;																		
						}	
						
					}
					else {
						if (sign < 0) {
							if (updated_pos + size >= bp_pos) {
								//frameshift = frameshift + (bp_pos - updated_pos)*sign;
								frameshift = frameshift + size*sign;
							} else {
								frameshift = frameshift + size*sign;
							}	
						} else {
							if (updated_pos + size >= bp_pos) {
								frameshift = frameshift + (bp_pos - updated_pos)*sign;
							} else {
								frameshift = frameshift + size*sign;
							}						
						}					
					}				
					
				}	
			}
			
			return frameshift;
		}
		
		//calcShiftDiff(bp, indels1, indels2, updatedIndels1, updatedIndels2)
		private int calcHomoBreakPoint(int pos, List<List<Integer>> indels1_sorted, List<List<Integer>> indels2_sorted) {
			//calculates difference in frameshift between two sequences (caused by indels) up to a certain nucleotide position	
			
			int frameshift1 = calcFrameShift(pos, indels1_sorted, false);
			int frameshift2 = calcFrameShift(Math.max(pos - frameshift1, 0), indels2_sorted, true);			
							
			return pos + frameshift2 - frameshift1;
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
				if ((updated_pos == pos || pos != -1) && updated_pos <= homobp && updated_pos+size > homobp && sign > 0) {
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
			
			for (int bp: breakPointsList) {
				
				int homo_bp = calcHomoBreakPoint(bp, sortedIndels1, sortedIndels2);						
				this.unrefinedHomologousBreaks.add(homo_bp);
				homo_bp = refineBP(bp, homo_bp);	
				
				
				homo_bp = Math.min(homo_bp, this.parentLengths[1]);		
				if (homo_bp < 0) {
					homo_bp = 0;
				}				
				homologousBreakPoints.add(homo_bp);				
			}		
			
			this.homologousBreakPointsList = homologousBreakPoints;
			return homologousBreakPoints;
		}
		
		public List<List<Integer>> getNewIndelList() {		
			List<List<Integer>> updated_indels1 = new ArrayList<List<Integer>>(updatedIndels1);
			List<List<Integer>> updated_indels2 =  new ArrayList<List<Integer>>(updatedIndels2);
			List<List<Integer>> updated_indels3 =  new ArrayList<List<Integer>>(updatedIndels1);
			//Now use the updated lists to merge them according to breakpoints
			//Need to add consideration for end of sequence as breakpoint
			List<Integer> breakPoints_with_end = new ArrayList<Integer>(breakPointsList);
			List<Integer> homoBreakPoints_with_end = new ArrayList<Integer>(homologousBreakPointsList);			
			
			//Length of recombinant is as follows:
			//parents are sorted such that parent1 < parent2
			//listBreakPoints.get(0) + parents[1].getLength() - homologousBreakPoints.get(0)
			int recombinant_len = 0;
			int bp_size = homologousBreakPointsList.size();
			if (bp_size == 1) {
				//recombinant_len = parentLengths[1] - (homologousBreakPointsList.get(0) - breakPointsList.get(0));
				recombinant_len = breakPointsList.get(0) + parentLengths[1] - homologousBreakPointsList.get(0);
			} else {
				recombinant_len = parentLengths[0] - (homologousBreakPointsList.get(1) - breakPointsList.get(1));
			}			
			
			breakPoints_with_end.add(recombinant_len);
			homoBreakPoints_with_end.add(recombinant_len);
			List<List<Integer>> merged_indels =  new ArrayList<List<Integer>>();		
			//Looping through breakpoints and adding indels that fall between consecutive ones	
			//selectIndels(List<List<Integer>> oldList, List<List<Integer>> currentList, int start, int end, int shift_diff)
			
			//For very start of genome and very end, deletions should been handled somewhat differently
			//included if == start only start = 0. included if end == end only if end == genomelength
			if (bp_size == 1) {
				updated_indels1 =  selectIndels(indels1, updated_indels1, 0, breakPoints_with_end.get(0), 0, 0);
				updated_indels2 =  selectIndels(indels2, updated_indels2, homoBreakPoints_with_end.get(0),  
						parentLengths[1],  homoBreakPoints_with_end.get(0) - breakPoints_with_end.get(0), 1);	
				merged_indels.addAll(updated_indels1);
				merged_indels.addAll(updated_indels2);				
			} else {
				updated_indels1 =  selectIndels(indels1, updated_indels1, 0, breakPoints_with_end.get(0), 0, 0);
				updated_indels2 =  selectIndels(indels2, updated_indels2, homoBreakPoints_with_end.get(0),  
						homoBreakPoints_with_end.get(1),  homoBreakPoints_with_end.get(0) - breakPoints_with_end.get(0), 2);
				
				int length_diff = (homologousBreakPointsList.get(1) - homologousBreakPointsList.get(0)) - 
						(breakPointsList.get(1) - breakPointsList.get(0));
				
				updated_indels3 = selectIndels(indels1, updated_indels3, breakPoints_with_end.get(1), parentLengths[0], 
						length_diff, 3);		
					
				merged_indels.addAll(updated_indels1);
				merged_indels.addAll(updated_indels2);
				merged_indels.addAll(updated_indels3);
			}
			
			this.del_shift = 1;
				
			
			merged_indels = sortIndelsByPosition(merged_indels);	
			merged_indels = mergeOverlaps(merged_indels);	
			
			return merged_indels;

		}
		
	}
	
    static SimpleSequence getRecombinantSeq(SimpleSequence[] parents, SortedSet<Integer> breakPoints) {
	 	assert(parents.length == 2);
		//assert(parents[0].getLength() <= parents[1].getLength());
	 	int max_length = Math.max(parents[0].getLength(), parents[1].getLength());
				
		List<List<Integer>> indels1 = new ArrayList<List<Integer>>(parents[0].getIndelList());
		List<List<Integer>> indels2 = new ArrayList<List<Integer>>(parents[1].getIndelList());
		indelCollection indels = new indelCollection(indels1, indels2, breakPoints, new int[]{parents[0].getLength(), parents[1].getLength()});
		
		List<Integer> listBreakPoints = new ArrayList<Integer>(indels.getBreakPoints());	
		List<Integer> homologousBreakPoints = new ArrayList<Integer>(indels.getHomologousBreakPoints());
		
		//checking if homo_bp is at end of genome, then we don't have to do anything (except maybe chop off a bit?)
		if (homologousBreakPoints.get(0) >= parents[1].getLength()) {			
			homologous_breakpoints.add(homologousBreakPoints.get(0));
			
			//take parents[0], but only until last breakpoint			
			return parents[0];
		}
		
		
		//List<Integer> homologousBreakPoints = 
		//homologousBreakPoints = calcHomologousBreakpoints();
		int recombinant_len = 0;
		int bp_size = homologousBreakPoints.size();
		if (bp_size == 1) {			
			recombinant_len = listBreakPoints.get(0) + parents[1].getLength() - homologousBreakPoints.get(0);					
			homologous_breakpoints.add(homologousBreakPoints.get(0));
		} else {
			recombinant_len = parents[0].getLength() + 
					((homologousBreakPoints.get(1) - homologousBreakPoints.get(0)) -  (listBreakPoints.get(1) - listBreakPoints.get(0)));			
		}
		
				
		int lastBreakPoint = 0;		// previous recombination location
		int currentSeq = 0;			// index of currently selected parent
		//int newlen = getRecombinantLength(parents, breakPoints);
		SimpleSequence product = new SimpleSequence(recombinant_len);

		byte[] dest = product.states;	// where to put the product
		SimpleSequence seq = parents[currentSeq];			
		
		/*
		System.out.println("!");
		System.out.println(recombinant_len);
		System.out.println(listBreakPoints);
		System.out.println(homologousBreakPoints);
		System.out.println("-----------------");
		System.out.println(parents[0].getLength());
		System.out.println(parents[1].getLength());		
		System.out.println("");
		*/
		
		int lastHomologous = 0;
		int counter = 0;
		for (int nextBreakPoint : breakPoints) {				
			int homologousNextBreakPoint = homologousBreakPoints.get(counter);			
			if (counter == 1) {		
				//System.arraycopy(seq.states, lastHomologous, 
						 //dest, lastBreakPoint, homologousNextBreakPoint-lastHomologous);
				System.arraycopy(seq.states, lastHomologous, 
						 dest, lastBreakPoint, homologousNextBreakPoint-lastHomologous);
			} else {
				//ERROR HAPPENS HERE
				System.arraycopy(seq.states, lastBreakPoint, 
						 dest, lastBreakPoint, nextBreakPoint-lastBreakPoint);
			}			
									
			//lastBreakPoint = homologousNextBreakPoint;
			lastBreakPoint = nextBreakPoint;
			lastHomologous = homologousNextBreakPoint;
			homologousBreakPoints.add(lastBreakPoint);
			currentSeq = 1 - currentSeq;
			seq = parents[currentSeq];	
			counter++;
		}
		
		if (counter == 2) {		
			System.arraycopy(seq.states, breakPoints.last(), 
					 dest, breakPoints.first() + (homologousBreakPoints.get(1) - homologousBreakPoints.get(0)), 
					 seq.getLength()-breakPoints.last());
		} else {			
			System.arraycopy(seq.states, homologousBreakPoints.get(0), 
					 dest, breakPoints.last(), seq.getLength()-homologousBreakPoints.get(0));
		}
		
		
		//String recombinant = product.getNucleotides();				
		
		List<List<Integer>> newIndels = indels.getNewIndelList();
		product.setIndelList(newIndels);			
		
		return(product);
	}	
	
static int last_homo;
	
    static SimpleSequence getRecombinantSequence(SimpleSequence[] parents, SortedSet<Integer> breakPoints) {
    	homologous_breakpoints.clear();
    	normal_breakpoints.clear();
    	
    	
	 	if (breakPoints.size() == 1) {	 
	 		normal_breakpoints.add(breakPoints.first());
	 		return getRecombinantSeq(parents, breakPoints);
	 		
	 	} else {		
	 		
	 		SortedSet<Integer> new_breaks = new TreeSet<Integer>();
	 		int first_break = breakPoints.first();
	 		new_breaks.add(first_break);
	 		normal_breakpoints.add(first_break);
	 		
	 		SimpleSequence seq1 = getRecombinantSeq(parents, new_breaks);
	 		
	 		new_breaks.clear();	 	
	 		
	 		int last_bp = breakPoints.last();	 			 				 		
	 		
	 		//sometimes homology leads to small events being ignored basically, since first homologous breakpoint will lie further
	 		//than the second non-homologous breakpoint. This is a non-sensical event, caused since non-homologous breakpoints are
	 		//chosen without respecting homology.
	 		
	 		//Can fix by just moving things to the right a bit
	 		int first_homo = homologous_breakpoints.first();
	 	
	 		if (first_homo >= last_bp) {
	 			last_bp += first_homo - breakPoints.first();	 			
	 		}
	 		
	 		last_bp = Math.min(last_bp, seq1.getLength());
	 		new_breaks.add(last_bp);
	 		normal_breakpoints.add(last_bp);
	 		
	 		SimpleSequence[] parents_rev = new SimpleSequence[] {seq1, parents[0]};
	 		SimpleSequence final_seq = getRecombinantSeq(parents_rev, new_breaks);	
	 		
	 		//int final_bp = Math.min(last_homo, final_seq.getLength());	 		
	 			 			
	 		return final_seq;	 		
	 		
	 	}
	}	
	
static SortedSet<Integer> homologous_breakpoints = new TreeSet<Integer>();
static SortedSet<Integer> normal_breakpoints = new TreeSet<Integer>();
 
 	public SortedSet<Integer> get_homologous_breakpoints(){
 		return homologous_breakpoints;
 	}
 	
 	public SortedSet<Integer> get_normal_breakpoints(){
 		return normal_breakpoints;
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
