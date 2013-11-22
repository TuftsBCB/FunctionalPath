import java.io.*;
import java.util.*;

public class IO {


	public static void loadCoveredProteins(String filename, HashMap<String, HashMap<String, Integer>> indices, HashMap<String, ArrayList<GOTerm>> collections) {
		
		try {
		
            BufferedReader in = new BufferedReader(new FileReader(filename));
            
            String thisline;
            while((thisline=in.readLine())!=null) {
                StringTokenizer st = new StringTokenizer(thisline, "\t");
                
                String protein = st.nextToken();
                String goterm = st.nextToken();
                
                
                for(String category : collections.keySet()) {
	                
	                try { 
		                ArrayList<GOTerm> collection = collections.get(category);
		                int index = indices.get(category).get(goterm).intValue();

		                collection.get(index).addProtein(protein);
		            } catch(NullPointerException npe) { }

                }

            }
            
            in.close();
            
        } catch(IOException ioe) { System.err.println(ioe.getMessage()); }

		
	}



    public void execute(String filename, ArrayList<ArrayList<GOTerm>> collections, ArrayList<HashMap<String, Integer>> indices) {
        
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            
            String thisline;
            while((thisline=in.readLine())!=null) {
                StringTokenizer st = new StringTokenizer(thisline, "\t");
                
                String protein = st.nextToken();
                String goterm = st.nextToken();

                for(int i=0;i<collections.size();i++) {
                	try {
	                	collections.get(i).get(indices.get(i).get(goterm).intValue()).addProtein(protein);
	                } catch(Exception e) { }
                }

            }
            
            in.close();
        } catch(IOException ie) { System.err.println(ie.getMessage()); }

    }



	// print out functional paths to a file
	public static void printPaths2File(String filename, ArrayList<String> paths) {
	
		try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filename));
            
            for(int i=0;i<paths.size();i++) {
                out.write(paths.get(i) + "\n");
            }
            
            out.close();
            
        } catch(IOException ioe) { System.err.println(ioe.getMessage()); }
	
	}


	// print out number of covered proteins for each GO term
	public static void printNumProteins2File(String[] categories, HashMap<String, ArrayList<GOTerm>> collections) {
		
		for(int i=0;i<categories.length;i++) {
		
			ArrayList<GOTerm> collection = collections.get(categories[i]);
		
			try {
		
				BufferedWriter out = new BufferedWriter(new FileWriter("protein.counts." + categories[i]));
				
				for(int j=0;j<collection.size();j++) {
					
					out.write(collection.get(j).getID() + "\t" + collection.get(j).getNumProteins() + "\n");
					
				}
				
				out.close();	
			
			} catch(IOException ioe) { }
		
		}
	
	}



	// load GeneOntology file and save the data
	public static void loadGOFile(String filename, HashMap<String, HashMap<String, Integer>> indices, HashMap<String, ArrayList<GOTerm>> collections, HashMap<String, ArrayList<String>> alternativeID) {
	
		
		try {
		    
		    
		    BufferedReader in = new BufferedReader(new FileReader(filename));
  
  			// parsing through metadata part            
            while(!in.readLine().equals("[Term]")) {}

            ArrayList<String> readLines = new ArrayList<String>();
            String thisline = in.readLine();  
            
			boolean done = false;
            while(!done) {
				
                if(thisline.equals("[Term]")) {
  					doProcessing(readLines, indices, collections, alternativeID);
       
                    // clear all the lines
                    readLines.clear();
                }
                else if(thisline.equals("[Typedef]")) { 
	                done = true; 
	            }
                else {
                   readLines.add(thisline);
                }
				
				thisline = in.readLine().trim();
            }
            doProcessing(readLines, indices, collections, alternativeID);

            in.close();
		
		
		} catch(IOException ioe) { System.err.println(ioe.getMessage()); }	
	
	
	}


	// process an arraylist of strings that are for one term, and store information to data structure
	private static void doProcessing(ArrayList<String> readLines, HashMap<String, HashMap<String, Integer>> indices, HashMap<String, ArrayList<GOTerm>> collections, HashMap<String, ArrayList<String>> alternativeID) {
	
		
                // do processing
                String id = readLines.get(0).replace("id: ", "").trim();
                String name = readLines.get(1).replace("name: ", "").trim();
                String namespace = readLines.get(2).replace("namespace: ", "").trim();
                
                                
                // create a GOTerm node
                GOTerm node = new GOTerm(id, namespace);
                
                
                // find right collection to save the node
                ArrayList<GOTerm> collection = collections.get(namespace);
                HashMap<String, Integer> index = indices.get(namespace);



                // store alternative ids
                ArrayList<String> alts = new ArrayList<String>();

                boolean obsolete = false;
                ArrayList<Integer> parents = new ArrayList<Integer>();
                
                
                for(int x=3;x<readLines.size();x++) {
                        String line = readLines.get(x);
                        
                        if(line.contains("alt_id:")) {
                                alts.add(line.replace("alt_id: ", "").trim());
                        }
                        
                        if(line.contains("is_obsolete:")) {
                                obsolete = true;
                        }

                        
                        if(line.contains("is_a:") || (line.contains("relationship:") && !line.contains("regulates"))) {
                        
                                int subind = line.indexOf("GO:");
                                String parentid = line.substring(subind, subind+10);
                        
                                if(index.containsKey(parentid)) {
                                        int ind = index.get(parentid).intValue();
                                        parents.add(ind);
                                }
                                else {
                                        GOTerm parent = new GOTerm(parentid, namespace);
                                        collection.add(parent);
                                        index.put(parentid, collection.size()-1);
                                        parents.add(collection.size()-1);
                                        
                                }
                        }
                        
                }


                // if the term is NOT obsolete and not in the collection
                if(!obsolete && !index.containsKey(id)) {
                           		
            		node.addParents(parents);
            		collection.add(node);
                    index.put(id, collection.size()-1);
 
                    // add this to parent's children list
                    for(int i=0;i<parents.size();i++) {
                        int parentind = parents.get(i).intValue();
                        collection.get(parentind).addChild(index.get(id));
                    }

                    alternativeID.put(id, alts);
                }
                else if(index.containsKey(id)) {
                    
                    int ind = index.get(id).intValue();
                    collection.get(ind).addParents(parents);
                    
                    // add this to parent's children list
                    for(int i=0;i<parents.size();i++) {
	                    int parentind = parents.get(i).intValue();
                        collection.get(parentind).addChild(index.get(id));
                    }                    

                    
                }		
	
	}






	public static void logging(String log) {
		System.err.println(log);
	}
	
	
	public static void logging(String filename, String log) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(filename));
			out.write(log);
			out.close();
		} catch(IOException ioe) { }
	
	}


	public static void printUsage() {
	
		System.err.println("Wrong Usage: consult README file.");
		System.err.println("Scenario 1: build all possible functional paths for three categories: BP, CC and MF.");
		System.err.println("            ./run [path_to_gene_ontology_file]");
		System.err.println();
		System.err.println("Scenario 2: build all possible functional paths for three categories: BP, CC and MF & identify informative nodes (GO terms that cover more than or equal to (t) number of proteins.");
		System.err.println("            ./run [path_to_gene_ontology_file] [t: threshold_for_identifying_informative_nodes] [functional_annotation_file: BP] [functional_annotation_file: CC] [functional_annotation_file: MF]");
		System.err.println("            Input file format for functional annotation file (tab delimited): [protein_name]     [GO term]");
		System.err.println("            i.e., Q0045   GO:0006123");
		System.err.println("                  Q0045   GO:0009060");
		System.err.println("                  Q0050   GO:0006315");
	
		
		System.err.println();
		System.err.println("Scenario 3: verbose version of Scenario 2");
		System.err.println("            Add -v flag to the end of arguments");
		System.err.println("            /run [path_to_gene_ontology_file] [t: threshold_for_identifying_informative_nodes] [functional_annotation_file: BP] [functional_annotation_file: CC] [functional_annotation_file: MF] -v");
		
	
		System.exit(1);
	
	}






}