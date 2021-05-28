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

		List<List<Integer>> reversedIndels = new ArrayList<List<Integer>>(this.indelList);
		Collections.reverse(reversedIndels);

		int tempPos = pos;

		for (List<Integer> indels : this.indelList) {

			if (indels.get(0) < tempPos) {
				cumChanges = cumChanges + indels.get(1);
				tempPos = tempPos - indels.get(1);
			}
		}

		return cumChanges;
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
		indel.add(start);
		indel.add(source.getLength());
		indel.add(countChanges(start));

		this.indelList.add(indel);

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
	class LengthComparator implements Comparator<List<Integer>> {
        public int compare(final List<Integer> list1, final List<Integer> list2) {        
            if (list1.get(0) == list2.get(0)) {
                return Integer.compare(list1.get(0), list2.get(0));
            }
            else {
                return Integer.compare(list1.get(0), list2.get(0));
            }        
        }
    }

    public List<List<Integer>> sortIndels(List<List<Integer>> insertions) {
        List<List<Integer>> withoutOverlaps = insertions;      
        Collections.sort(withoutOverlaps, new LengthComparator());
        return withoutOverlaps;
    }

    //returns an updated indelList, which is created by combining the lists from both parents in a "special" way
    //the merger ensures events can still be "unpacked" in reverse order without disrupting the coordinate system

	static List<List<Integer>> getNewIndelList(SimpleSequence[] parents, SortedSet<Integer> breakPoints) {

		List<List<Integer>> indels1 = new ArrayList<List<Integer>>(parents[0].getIndelList());
		List<List<Integer>> indels2 = new ArrayList<List<Integer>>(parents[1].getIndelList());
		List<List<Integer>> updated_indels1 = new ArrayList<List<Integer>>();
		List<List<Integer>> updated_indels2 = new ArrayList<List<Integer>>();

		List<List<Integer>> recombinant_indel_list = new ArrayList<List<Integer>>();

		//System.out.println("NEW INDEL LIST BEING CREATED: ");
		//System.out.println(breakPoints);
		//System.out.println(indels1);
		//System.out.println(indels2);
		//System.out.println("");		
		
		if (indels1.size() == 0 && indels2.size() == 0) {
			//System.out.println("-------------------------------------");
			return indels1;
		}
		//generating updated lists for both parents
		//updated = positions are updated to time of recombination
		List<List<Integer>> currentList;
		List<List<Integer>> updated_currentList;
		for (int h = 0; h < 2; h++) {

			if (h==0) {
				currentList = indels1;
				updated_currentList = updated_indels1;
			}
			else {
				currentList = indels2;
				updated_currentList = updated_indels2;
			}

			if (currentList.size() >= 2) {

				for (int i = 0; i < currentList.size()-1; i++) {
					List<Integer> baseIndel = new ArrayList<Integer>(currentList.get(i));
					int newTranslation = 0;

					for (int j = i + 1; j < currentList.size(); j++) {
						List<Integer> newerIndel = new ArrayList<Integer>(currentList.get(j));

						//newer indel lies before base indel, and thus base indel's position is changed at time of recombination
						if (newerIndel.get(0)-newerIndel.get(2) < baseIndel.get(0)-baseIndel.get(2)) {
							//adding size of indel
							newTranslation = newTranslation + newerIndel.get(1);
						}

					}

					List<Integer> updated_indel = new ArrayList<Integer>();
					//Calculating new position and translation
					//New translation such that position+translation is equal for old and new indel

					//*NOTE!!!: should this be +newTranslation or -newTranslation???
					int newpos = baseIndel.get(0)+newTranslation;
					updated_indel.add(newpos);					
					updated_indel.add(baseIndel.get(1));
					updated_indel.add(baseIndel.get(2)-newTranslation);

					updated_currentList.add(updated_indel);
					newTranslation = 0;
				}

				//adding final event
				List<Integer> finalIndel = new ArrayList<Integer>(currentList.get(currentList.size()-1));
				updated_currentList.add(finalIndel);

			} else if (currentList.size() == 1) {
				updated_currentList.add(currentList.get(0));
			}

		}

		//Now sort both updated indel lists
		updated_indels1 = new SimpleSequence().sortIndels(updated_indels1);
		updated_indels2 = new SimpleSequence().sortIndels(updated_indels2);
		//recombinant_indel_list
		//System.out.println(updated_indels1);
		//System.out.println(updated_indels2);
		//System.out.println("");

		int currentStart = 0;
		int currentEnd = 0;
		boolean whichList = true;

		//Now use the updated lists to merge them according to breakpoints
		//Need to add consideration for end of sequence as breakpoint
		SortedSet<Integer> breakPoints_with_end = new TreeSet<Integer>(breakPoints);
		if (breakPoints_with_end.size() == 1 ) {

			breakPoints_with_end.add(parents[1].getLength());
		} else {

			breakPoints_with_end.add(parents[0].getLength());
		}
		
		//Looping through breakpoints and adding indels that fall between consecutive ones
		//Also compare the accumulated frameshifts between the two sequences. If different, they will be misaligned after breakpoint.
		//Thus need to add an indel event to compensate.

		//Comparing frameshifts between breakpoints and adding indel to adjust for it.	
		int startp = 0;
		int endp = 0;
		int c = 0;	

		//int[] frameShift = new int[]{0,0};
		int[] frameShift = new int[]{0,0};			
	
		for (int bpoint: breakPoints) {

			//int temp1 = frameShift[c];
			//int temp2 = frameShift[1-c];		
			//frameShift[c] = temp2;
			//frameShift[1-c] = temp1;
			
			endp = bpoint;

			for (List<Integer> indel1: updated_indels1) {				
				if (indel1.get(0) < endp && indel1.get(0) >= startp) {
					if (indel1.get(0) + indel1.get(1) >= endp) {
						int sign1 = (int) Math.signum(indel1.get(1));
						frameShift[0] = frameShift[0] + ((endp - indel1.get(0))*sign1);
					} else {
						frameShift[0] = frameShift[0] + indel1.get(1);
					}					
				}
			}
			for (List<Integer> indel2: updated_indels2) {
				if (indel2.get(0) < endp && indel2.get(0) >= startp) {
					if (indel2.get(0) + indel2.get(1) >= endp) {
						int sign2 = (int) Math.signum(indel2.get(1));					
						frameShift[1] = frameShift[1] + ((endp - indel2.get(0))*sign2);
					} else {
						frameShift[1] = frameShift[1] + indel2.get(1);
					}					
				}				
			}

			int shift_diff = frameShift[1-c] - frameShift[c];			

			//System.out.println(shift_diff);

			if (shift_diff > 0) {
				List<Integer> shiftIndel = new ArrayList<Integer>();
				//shiftIndel.add(bpoint);
				shiftIndel.add(bpoint-shift_diff);
				shiftIndel.add(shift_diff);
				shiftIndel.add(frameShift[c]);
				recombinant_indel_list.add(shiftIndel);		
			} else if (shift_diff < 0) {
				List<Integer> shiftIndel = new ArrayList<Integer>();
				shiftIndel.add(bpoint);
				shiftIndel.add(shift_diff);
				shiftIndel.add(frameShift[c]);
				recombinant_indel_list.add(shiftIndel);
			}

			c = 1-c;
			startp = endp;			
		}		
		

		//Looping through breakpoints and adding indels that fall between consecutive ones,		
		for (int bp : breakPoints_with_end) {
			currentEnd = bp;
			int counter = 0;	

			//swapping between parental sequences 
			if (whichList) {
				currentList = updated_indels1;
			}
			else {
				currentList = updated_indels2;
			}

			if (currentList.size() > 0) {

				List<Integer> startIndel = new ArrayList<Integer>(currentList.get(0));

				while (startIndel.get(0) < currentEnd && counter < currentList.size()) {
					startIndel = currentList.get(counter);					
					int pos = startIndel.get(0);
					int size = Math.abs(startIndel.get(1));
					int sign = (int) Math.signum(startIndel.get(1));
					
					if (currentStart <= pos && (pos + size) < currentEnd) {
						//Indel cleanly within block
						recombinant_indel_list.add(startIndel);

					} else if ( ((currentStart <= pos) && (pos < currentEnd)) && ((pos + size) >= currentEnd) ) {
						//Indel divided by end breakpoint
						List<Integer> dividedIndel = new ArrayList<Integer>();
						dividedIndel.add(pos);
						dividedIndel.add((currentEnd - pos)*sign);
						dividedIndel.add(startIndel.get(2));

						recombinant_indel_list.add(dividedIndel);

					} else if ( (pos < currentStart) && ((pos + size) > currentStart) && ((pos + size) < currentEnd)) {
						//Indel divided by starting breakpoint
						List<Integer> dividedIndel = new ArrayList<Integer>();
						dividedIndel.add(currentStart);
						dividedIndel.add((pos + size - currentStart)*sign);
						dividedIndel.add(startIndel.get(2));

						recombinant_indel_list.add(dividedIndel);

					} else {

					}

					counter++;							
				}

			}	

			whichList = !whichList;
			currentStart = bp;
		}		
		
		//System.out.println("FINAL LIST: ");
		//System.out.println(recombinant_indel_list);

		recombinant_indel_list = new SimpleSequence().sortIndels(recombinant_indel_list);

		//System.out.println(recombinant_indel_list);
		//System.out.println("-------------------------------------");

		return recombinant_indel_list;

	}


    static SimpleSequence getRecombinantSequence(SimpleSequence[] parents, SortedSet<Integer> breakPoints) {
	 	assert(parents.length == 2);
		assert(parents[0].getLength() <= parents[1].getLength());

		int lastBreakPoint = 0;		// previous recombination location
		int currentSeq = 0;			// index of currently selected parent
		int newlen = getRecombinantLength(parents, breakPoints);
		SimpleSequence product = new SimpleSequence(newlen);

		byte[] dest = product.states;	// where to put the product
		SimpleSequence seq = parents[currentSeq];
		for (int nextBreakPoint : breakPoints) {
			System.arraycopy(seq.states, lastBreakPoint, 
							 dest, lastBreakPoint, nextBreakPoint-lastBreakPoint);
			
			lastBreakPoint = nextBreakPoint;
			currentSeq = 1 - currentSeq;
			seq = parents[currentSeq];
		}
		int nextBreakPoint =  seq.getLength();
		System.arraycopy(seq.states, lastBreakPoint, 
						 dest, lastBreakPoint, nextBreakPoint-lastBreakPoint);

		List<List<Integer>> newIndels = new ArrayList<List<Integer>>(getNewIndelList(parents, breakPoints));
		product.setIndelList(newIndels);

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
