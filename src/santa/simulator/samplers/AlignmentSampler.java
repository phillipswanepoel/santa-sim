package santa.simulator.samplers;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.HashMap;

import santa.simulator.Random;
import santa.simulator.Virus;
import santa.simulator.genomes.AminoAcid;
import santa.simulator.genomes.Feature;
import santa.simulator.genomes.Nucleotide;
import santa.simulator.population.Population;
import santa.simulator.replicators.RecombinantTracker;
import santa.simulator.replicators.RecombinationEvent;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: AlignmentSampler.java,v 1.6 2006/07/18 07:37:47 kdforc0 Exp $
 * Eddited by Abbas Jariani Sep 2013
 */
public class AlignmentSampler implements Sampler {
    public enum Format {
        FASTA,
        NEXUS,
        XML
    };

    private final Feature feature;
    private final Set<Integer> sites;
    //Abbas: final modifier was removed from sampleSize
    private int sampleSize;
    private Format format;
    private String label;
    private String fileName;
    private PrintStream destination;
    private Map<Integer,Integer> schedule;
    private int replicate;
    private boolean consensus;

    /**
     * Construct an alignment sampler
     * @param sampleSize  amount of sequences to sample at regular intervals
     * @param consensus   write the consensus sequence of the sample rather than writing the sample ?
     * @param schedule    amount of sequences to sample at irregular intervals
     * @param format      format
     * @param label       label with possible %g, %s, %t, and %f variables
     * @param fileName    name of the file to write the samples
     */
    public AlignmentSampler(Feature feature, Set<Integer> sites, int sampleSize, boolean consensus,
                            Map<Integer,Integer> schedule, Format format, String label, String fileName) {
        this.format = format;
        this.fileName = fileName;

        if (label == null) {
            this.label = "virus_%g_%s";
        } else {
            this.label = label;
        }

        this.feature = feature;
        this.sites = sites;
        this.fileName = fileName;
        this.sampleSize = sampleSize;
        this.consensus = consensus;
        this.schedule = schedule;
    }

    public void initialize(int replicate) {

        this.replicate = replicate;
        String fName = substituteVariables(fileName, 0, 0, 0.0);

        try {
            destination = new PrintStream(fName);
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Could not open file for writing: " + fName);
        }

        if (format == Format.NEXUS) {
            destination.println("#NEXUS");
            destination.println();
            destination.println("BEGIN DATA;");

            int nchar = sites.size();
            int samplesize = computeSampleSize();

            destination.println("\tDIMENSIONS NTAX=" + samplesize + " NCHAR=" + nchar + ";");
            if (feature.getFeatureType() == Feature.Type.AMINO_ACID) {
                destination.println("\tFORMAT DATATYPE=PROTEIN GAP=-;");
            } else {
                destination.println("\tFORMAT DATATYPE=NUCLEOTIDE GAP=-;");
            }

            destination.println("\tMATRIX");
        }
    }

    private String substituteVariables(String name, int generation, int seq, double fitness) {
        String result = name.replaceAll("%r", String.valueOf(replicate+1));
        result = result.replaceAll("%g", String.valueOf(generation));
        result = result.replaceAll("%s", String.valueOf(seq));
        result = result.replaceAll("%f", String.valueOf(fitness));
        return result;
    }

    private int computeSampleSize() {
        if (schedule == null)
            return 0; // do not know really, should be samplesize * number of samples
        else {
            int result = 0;

            for (Integer item : schedule.values()) {
                result += item;
            }

            return result;
        }
    }

    public void sample(int generation, Population population) {
        Virus[] sample = getSample(generation, population);

        if (sample != null) {
            if (format == Format.NEXUS) {
                writeNexusFormat(generation, sample);
            } else if (format == Format.XML) {
                writeXMLFormat(generation, sample);
            } else if (format == Format.FASTA) {
                writeFastaFormat(generation, sample);
            }
        }
    }

    protected Virus[] getSample(int generation, Population population) {

        /**
         * Abbas: if the population size is smaller than sampleSize we get error
         * so the following minimum is applied
         * also final modifier was removed from sampleSize
         */
    	if (schedule == null) {
            List<Virus> viruses = population.getCurrentGeneration();
            if (sampleSize>viruses.size())
            {
            	System.out.println("warning: sampleSize of alignmentSampler was shrunk because of small population size to "+ viruses.size());
            }	
        	sampleSize = Math.min(sampleSize,viruses.size());

            Object[] tmp = Random.nextSample(viruses, sampleSize);
            Virus[] sample = new Virus[tmp.length];
            System.arraycopy(tmp, 0, sample, 0, tmp.length);
            return sample;
        } else {
            if (schedule.containsKey(generation)) {
                int count = schedule.get(generation);
                List<Virus> viruses = population.getCurrentGeneration();
                Object[] tmp = Random.nextSample(viruses, count);
                Virus[] sample = new Virus[tmp.length];
                System.arraycopy(tmp, 0, sample, 0, tmp.length);
                return sample;
            } else
                return null;
        }
    }

    private void writeNexusFormat(int generation, Virus[] sample) {
        if (consensus) {
            String l = substituteVariables(label, generation, 0, 0.0);

            destination.print(l + "\t");
            destination.println(computeConsensus(sample));
        } else {
            int i = 1;
            for (Virus virus : sample) {
                String l = substituteVariables(label, generation, i, virus.getFitness());

                destination.print(l + "\t");

                byte[] states = virus.getGenome().getStates(feature);
                if (feature.getFeatureType() == Feature.Type.AMINO_ACID) {
                    for (int site : sites) {
                        destination.print(AminoAcid.asChar(states[site]));
                    }
                    destination.println();
                } else {
                    for (int site : sites) {
                        destination.print(Nucleotide.asChar(states[site]));
                    }
                    destination.println();

                }
                i++;
            }
        }
    }

    private String computeConsensus(Virus[] sample) {
        String result = "";

        byte[][] states = new byte[sample.length][];
        int j = 0;
        for (Virus virus : sample) {
            states[j] = virus.getGenome().getStates(feature);
            j++;
        }

        int freqs[] = new int[feature.getAlphabet().getStateCount()];
        for (int site : sites) {
            int maxfreq = 0;
            byte maxS = Nucleotide.A;

            for (int i = 0; i < states.length; i++) {
                byte s = states[site][i];
                freqs[s]++;
                if (freqs[s] > maxfreq) {
                    maxfreq = freqs[s];
                    maxS = s;
                }
            }

            /* for now just take the maximum one */
            if (maxfreq < sample.length/2)
                result += "N";
            else
                result += Nucleotide.asChar(maxS);
        }

        return result;
    }

    class LengthComparator implements Comparator<List<Integer>> {
        public int compare(final List<Integer> list1, final List<Integer> list2) {        
            if (list1.get(0) == list2.get(0)) {
                return Integer.compare(list1.get(1), list2.get(1));
            }
            else {
                return Integer.compare(list1.get(0), list2.get(0));
            }        
        }
    }

    private List<List<Integer>> sortIndels(List<List<Integer>> insertions) {
        List<List<Integer>> withoutOverlaps = insertions;      
        Collections.sort(withoutOverlaps, new LengthComparator());
        return withoutOverlaps;
    }

    //adds size gaps to string starting at position pos
    private String addGaps(String s, int pos, int size) {    
        String gapChars = "-".repeat(size);
        s = s.substring(0, pos) + gapChars + s.substring(pos, s.length());
   
        return s;
    }
    
    //the idea is have a collection of insert_events
    //these will be then done in reverse order positionally (from end of sequence to beginning)
    //all viruses in affected_viruses will have subseq inserted, rest will have gap characters inserted
    private class Insert_Event {
        private int position;
        private int size;   
        private int translation;     
        private Set<Integer> affected_viruses = new HashSet<Integer>(); //idea is to list all viruses where seq inserted, rest gaps
        private String subseq;

        public Insert_Event(int pos, int s, int translation, String subs) {
            this.position = pos;
            this.size = s;  
            this.translation = translation;          
            this.subseq = subs;
        }

        public Insert_Event(int pos, int s, int translation, int virus_index, String subs) {
            this.position = pos;
            this.size = s;           
            this.subseq = subs;
            this.translation = translation;   
            this.affected_viruses.add(virus_index);
        }

        public void add_virus(int virus_index) {
            this.affected_viruses.add(virus_index);
        }

        public List<Integer> getIns() {
            List<Integer> ins = new ArrayList<>();
            ins.add(this.position);
            ins.add(this.size);
            ins.add(this.translation);
            return ins;
        }

        public Set<Integer> getViruses() {
            return this.affected_viruses;
        }

        public String getSubseq() {
            return this.subseq;
        }
    }

    
    //calculates the final position of insert, once all remaining insertions and deletions have happened
    private int calc_insert_pos(List<List<Integer>> indels, int insert_index) {

        int initial_pos = indels.get(insert_index).get(0);
        int finalpos = initial_pos;

        if (insert_index >= indels.size()-1) {
            return finalpos;
        }
        //add impact 
        for (int k = insert_index+1; k < indels.size(); k++) {
            if (indels.get(k).get(0) < initial_pos) {
                finalpos = finalpos - indels.get(k).get(1);
            }
        }

        return finalpos;
    }

    private List<Insert_Event> insertion_events = new ArrayList<>();
    private Set<Integer> seenInsertions = new HashSet<Integer>();
    private HashMap<Integer, Insert_Event> insertion_map = new HashMap<Integer, Insert_Event>();    

    //returns string with gap characters added and insertions removed
    private String remove_indels(String s, List<List<Integer>> indels, int virus_index) {   
    	
        for (List<Integer> indel : indels) {
            String subseq = "";
            
            int ini_size = s.length();
            int sizy = indel.get(1);

            if (indel.get(1) < 0) {
            	
            	if (indel.get(0) > s.length()) {            		
            		System.out.println("deletion out of bounds: ");
            		System.out.println(indel);
            		System.out.println(s.length());
            	}
            	
                s = addGaps(s, indel.get(0), -indel.get(1));
            }
            else {
            	
            	if (indel.get(0) >= s.length()) {            		
            		System.out.println("insertion out of bounds: ");
            		System.out.println(indel);
            		System.out.println(s.length());
            	}
            	
                subseq = s.substring(indel.get(0), indel.get(0)+indel.get(1));
                s = s.substring(0, indel.get(0)) + s.substring(indel.get(0)+indel.get(1), s.length());
              
                int indel_id = indel.get(3);

                if (seenInsertions.contains(indel_id)) {
                    insertion_map.get(indel_id).add_virus(virus_index);
                }
                else {
                    seenInsertions.add(indel_id);                   
                    //insertion_map.put(mod_indel,  new Insert_Event(indel.get(0), indel.get(1), indel.get(2), virus_index, subseq));
                    insertion_map.put(indel_id,  new Insert_Event(indel.get(0), indel.get(1), indel.get(2), virus_index, subseq));
                }

            }
            
            int final_size = s.length();
            
            if (final_size != ini_size - sizy) {
            	System.out.println("POESOP");
            	System.out.println(ini_size);
            	System.out.println(sizy);
            	System.out.println(final_size);
            	
            }            
        }
        
        return s;
    }
 

    private List<String> add_insertions(List<String> genomeStringsList) {

        List<String> inserted_strings = new ArrayList<String>(genomeStringsList);
        Set<Integer> keys = insertion_map.keySet();       
        List<List<Integer>> cum_changes = new ArrayList<List<Integer>>();       
        
        HashMap<Integer, List<List<Integer>>> gaps_added = new HashMap<Integer, List<List<Integer>>>();  

        for (Integer key : keys) {
            //System.out.println(key);
            //System.out.println(insertion_map.get(key).getSubseq());
            //System.out.println(insertion_map.get(key).getViruses());
            //System.out.println("");
        	
            Insert_Event temp = insertion_map.get(key);

            List<Integer> indel = new ArrayList<Integer>(temp.getIns());
            int pos = indel.get(0);
            int size = indel.get(1);           
            String subseq = temp.getSubseq();
            Set<Integer> affected = temp.getViruses();   

            for (int p = 0; p < inserted_strings.size(); p++) {
            	//THIS LOOP ALMOST CERTAINLY SOURCE OF SOME ISSUES
            	//HOW CAN YOU USE THE SAME TRANSLATION FOR ALL (WHEN RECOMBINATION IS A THING)
            	//rather than using temp from insertion_map, need to get pos+trans from original stringlist
            	int cum = 0;
            	
            	for (List<Integer> indy : cum_changes) {
            		int temp_pos = indy.get(0);
            		int temp_size = indy.get(1);
            		
            		if (temp_pos < pos) {
            			cum = cum + temp_size;
            		}
            	}
                if (affected.contains(p)) {
                	
                    String old_string = inserted_strings.get(p); 
                    String new_string = old_string.substring(0, pos + cum) + subseq + old_string.substring(pos + cum, old_string.length());

                    inserted_strings.set(p, new_string);
                    
                }
                else {     
                    String old_string = inserted_strings.get(p);
                    
                    if (pos- + cum > old_string.length()) { 
                    	System.out.println(pos+ cum);
                		System.out.println(pos);                		
                		System.out.println(cum); 
                		System.out.println("");   
                	}
                    
                    String new_string = addGaps(old_string, pos + cum, size);

                    inserted_strings.set(p, new_string);
                }
            }     
            
            cum_changes.add(indel);

        }

        return inserted_strings;
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
    
    static List<List<Integer>> sortIndelsByPosition(List<List<Integer>> insertions) {
        List<List<Integer>> withoutOverlaps = insertions;      
        Collections.sort(withoutOverlaps, new LengthComparator2());
        return withoutOverlaps;
    }
    
    private void writeFastaFormat(int generation, Virus[] sample) {
        Set<Integer> seenEvents = new HashSet<Integer>();
        List<String> genomeStringsList = new ArrayList<String>();  
        List<String> finalStringsList = new ArrayList<String>(); 
        List<List<Integer>> list_eventsList = new ArrayList<List<Integer>>(); 
        Integer endLength = 0;    

        if (consensus) {
            String l = substituteVariables(label, generation, 0, 0.0);

            destination.println(">" + l);
            destination.println(computeConsensus(sample));
        } else {
            int i = 1;   
            
            List<String> nameList = new ArrayList<String>();

            int counter = 0;
            //Writing each virus to fasta
            //Inserting deletions
            for (Virus virus : sample) {   

                //1. reverse chrono order
                //2. remove deletions AND insertions (save insertion sequence)
                //3. re-introduce insertions (insert seqs or gaps)
                //4. voila

                LinkedHashSet<Integer> eventList = virus.getRecombinationList(); 
                List<Integer> events_to_add = new ArrayList<>(eventList);               
                list_eventsList.add(events_to_add);

                String l = substituteVariables(label, generation, i, virus.getFitness());
                nameList.add(l);

                String genomeString = virus.getGenome().getSequence().getNucleotides();
                               
                List<List<Integer>> indelsList = new ArrayList<List<Integer>>(virus.getGenome().getSequence().getIndelList());                
                indelsList = sortIndelsByPosition(indelsList);
                System.out.println("> " + i);
                System.out.println(indelsList);    
                System.out.println("");           

                //System.out.println("GENOME: " + Integer.toString(counter+1)); 
                //System.out.println(indelsList);              
                Collections.reverse(indelsList); 
                
                genomeString = remove_indels(genomeString, indelsList, counter);
                  
                genomeStringsList.add(genomeString);
                //destination.println(">" + l + "_" + virus.getGenome() + "_" + eventList);
                //destination.println(genomeString);

                for (Integer event : eventList) {
                    seenEvents.add(event);
                    //add a counter here, number of occurences of events
                    //Then if an event appears in like 40%+ seqs, create another event that is the inverse.
                }

                counter++;
                i++;
            }              
            
            int temp_i = genomeStringsList.get(0).length();
            for (String k : genomeStringsList) {   
            	
            	if (k.length() != temp_i) {
            		System.out.println("OH POESA, all strings not equal length after removing indels!");
            	}
            	temp_i = k.length();    
            }
              
            finalStringsList = new ArrayList<String>(add_insertions(genomeStringsList));

            int ini_len = finalStringsList.get(0).length();
        
            for (String seq: finalStringsList) {

                if (seq.length() != ini_len) {
                    System.out.println("WARNING, NOT ALL SAME LENGTH. Misalignment very likely");
                    break;
                    //System.exit(1);
                }         

            }  

            for (int h = 0; h < finalStringsList.size(); h++) {
                //destination.println(">" + nameList.get(h) + "_" + sample[h].getGenome() + "_" + list_eventsList.get(h));
                //destination.println(">" + (h+1) + "_" + sample[h].getGenome());
                destination.println(">" + (h+1));
                destination.println(finalStringsList.get(h));                
            }


            //System.out.println("------------------------");

        }       
        

        //First remove all events that are present in all sequences
        //Remove from both recList and list_eventsList

        //Finding events present in ALL sequences
        List<Integer> events_to_remove = new ArrayList<Integer>(list_eventsList.get(0)); 

        for (int a = 1; a < list_eventsList.size(); a++) {            
            events_to_remove.retainAll(list_eventsList.get(a));
        }     

        //Removing those events
        for (List<Integer> events : list_eventsList) {
            events.removeAll(events_to_remove);
        } 
        seenEvents.removeAll(events_to_remove);
      
              
        //Printing recombination events that are seen in sample
        try {

            //System.out.println("LENGTH: ");
            endLength = finalStringsList.get(0).length();
            //System.out.println(endLength);

            PrintStream recombPrinter = new PrintStream("recombination_events.txt");
            recombPrinter.println("EventNum*Breakpoints*Generation*Recombinant*Parents");
            List<RecombinationEvent> recList = RecombinantTracker.recombinationList;

            for (int j = 0; j < recList.size(); j++) {
                if (seenEvents.contains(j)) {

                    RecombinationEvent event = recList.get(j);

                    //Modifying single breakpoints by adding start/end of genome.
                    //This is to make analysis at later stages easier.
                    SortedSet<Integer> breakpoints = new TreeSet<Integer>(event.getBreakpoints());                   
                    if (breakpoints.size() == 1) {
                        if (breakpoints.first() < endLength/2) {
                            breakpoints.add(0);
                        }
                        else {
                            breakpoints.add(endLength);
                        }
                    }

                    //Before writing breakpoints, modify them by adding endpoint.    
                    recombPrinter.println(j + "*" + breakpoints + "*" + 
                    event.getGeneration() + "*" + event.getRecombinant() + "*" + 
                    event.getParents());
                }                
            }
            recombPrinter.close();

        }
        catch(FileNotFoundException ex) {
            System.out.println("RECOMBINATION_EVENTS FILE NOT FOUND!");
        }

        //Printing all sampled sequence names together with associated recombination events
        try {

            PrintStream recombPrinter2 = new PrintStream("sequence_events_map.txt");
            recombPrinter2.println("Sequence*Events");        
            
            for (int h = 0; h < finalStringsList.size(); h++) {            
                recombPrinter2.println((h+1) + "*" + list_eventsList.get(h));                            
            }

            recombPrinter2.close();

        }
        catch(FileNotFoundException ex) {
            System.out.println("RECOMBINATION_EVENTS FILE NOT FOUND!");
        }

    }

    private void writeXMLFormat(int generation, Virus[] sample) {
        destination.println();
        destination.println("<sequences>");
        destination.println("<!-- Generation = " + generation + " -->");
        int i = 1;
        for (Virus virus : sample) {
            String l = substituteVariables(label, generation, i, virus.getFitness());
            destination.println("\t<sequence label=\"" + l + "\">");
            destination.println("\t\t" + virus.getGenome().getSequence().getNucleotides());
            destination.println("\t</sequence>");
            i++;
        }
        destination.println("</sequences>");
    }

    public void cleanUp() {
        if (format == Format.NEXUS) {
            destination.println("\t;");
            destination.println("END;");
        }

        destination.close();
        destination = null;
    }
}
