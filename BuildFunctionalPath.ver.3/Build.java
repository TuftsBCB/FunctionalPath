import java.io.*;
import java.util.*;
import java.text.*;

public class Build {

	// store alternavie ids
    public static HashMap<String, ArrayList<String>> alternativeID;
        
    // store indices of node
    // 0: biological process  1: cellular component  2: molecular function 
    public static  HashMap<String, HashMap<String, Integer>> indices;
        
    // collection of terms
    // 0: biological process  1: cellular component  2: molecular function 
    public static HashMap<String, ArrayList<GOTerm>> collections;
 

	public static String gofile;
	public static ArrayList<String> fafiles;
	public static int threshold;


    public static void main(String[] args) {

        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Date date = new Date();
        IO.logging("Process started: " + dateFormat.format(date));


   		// namespaces in GO file
    	String[] categories = {"biological_process", "cellular_component", "molecular_function"};


    	boolean verbose = false;
    	boolean identifyInformative = false;
    	
    	
    	// taking arguments
    	switch(args.length) {
	    	
	    	case 1:
	    		gofile = args[0];
	    		break;
			case 3:
			case 4:
			case 5:
			case 6:
				gofile = args[0];
				threshold = Integer.parseInt(args[1]);
				identifyInformative = true;
				
				fafiles = new ArrayList<String>();
				for(int i=2;i<args.length;i++) {
					if(args[i].equals("-v")) {
						verbose = true;
					}
					else {
						fafiles.add(args[i]);
					}
				}	    		
				break;
			default:
				IO.printUsage();
				break;
					
    	}
    	
    	
    	
    	alternativeID = new HashMap<String, ArrayList<String>>();
        
        
        // There are three ArrayList<GOTerm>, one for each category: BP, CC and MF
        indices = new HashMap<String, HashMap<String, Integer>>();
        collections = new HashMap<String, ArrayList<GOTerm>>();




    	// initialize collections
		for(int i=0;i<categories.length;i++) {
			indices.put(categories[i], new HashMap<String, Integer>());
			collections.put(categories[i], new ArrayList<GOTerm>());
		}    	
    	
    	
    	
		// read gene ontology file, build data structure and store information
        // each term is recognized when "[Term]" is observed in the file
		IO.loadGOFile(gofile, indices, collections, alternativeID);
       
                if(verbose) {
                        date = new Date();
                        IO.logging("Gene Ontology file is loaded: " + dateFormat.format(date));
                }


 
        
        // load all the functional annotations
        for(int i=0;i<fafiles.size();i++) {
	        IO.loadCoveredProteins(fafiles.get(i), indices, collections);
        }


               if(verbose) {
                        date = new Date();
                        IO.logging("Functional annotations are loaded: " + dateFormat.format(date));
                }
        


        // record the number of proteins that are covered by each GO term
	    if(verbose) IO.printNumProteins2File(categories, collections);
        
        
        // build functional paths
        for(int i=0;i<categories.length;i++) {
        
        	FunctionalPath temp = new FunctionalPath(categories[i], collections.get(categories[i]));
        	
        	if(verbose) temp.isVerbose(true);
        	
        	
      		if(identifyInformative) {
	      		temp.execute(threshold);
      		}
      		else {
	      		// if we are just building the paths without informative nodes
	      		temp.execute();
      		}     	
        
      		int divind = categories[i].indexOf("_");
      		String filename = (categories[i].substring(0,1) + categories[i].substring(divind+1, divind+2)).toUpperCase() + ".path";
      		
      		IO.printPaths2File(filename, temp.getFuncPaths());
        

                if(verbose) { 
                        date = new Date();
                        IO.logging("Building functional paths for " + categories[i] + "is done: " + dateFormat.format(date));
                }



        }
    
    
        // record the number of proteins that are covered by each GO term
	    if(verbose) IO.printNumProteins2File(categories, collections);

        

    }




}



