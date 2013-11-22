import java.io.*;
import java.util.*;

public class FileIO {

   
    // load functional paths 
    public static ArrayList<LinkedList<String>> loadFunctionalPath(String filename, HashMap<String, ArrayList<Integer>> index) {
        
        int lengthOfTerm = 10;
        ArrayList<LinkedList<String>> paths = new ArrayList<LinkedList<String>>();

        
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            
            String previous = "";
            String thisline;
            
            while((thisline=in.readLine())!=null) {
                
		// tokens is a collection of GO terms that did not appear in the functional path already built
                ArrayList<String> tokens = new ArrayList<String>();
                StringTokenizer st = new StringTokenizer(thisline, "\t");
                
                // skip the root: we don't count the root
                st.nextToken();
                

		// path for thisline
                LinkedList<String> path = new LinkedList<String>();
                
                
                // find the longest common substring to utilize already built paths
		// for the longest common substring, all the possible functional paths are already built.
		// therefore, to prevent building all the functional paths that already existed in the "path", we collect newly appeared GO terms in "tokens" 
		StringBuffer existing = new StringBuffer("");

                while(st.hasMoreTokens()) {
                    String next = st.nextToken();
                    
		    // remove suffix if existed
                    if(next.length() >= lengthOfTerm) {
                        next = next.substring(0,lengthOfTerm);
                    }

		    existing.append(next);

                    int len = existing.length();
		
                    if(previous.length() < len) {
                    	// if previous path is shorter than existing, it is a new term
			tokens.add(next); 
                    }
                    else if(previous.substring(0,len).equals(existing.toString())) {
                        // next is the new term after the longest functional path that are already built
			// at this moment the existing is the longest common substring between previous and thisline
			path.add(next);
		
                    }
                    else { tokens.add(next); }
                    
                }
		// after this, "path" is the longest common substring between the previous and thisline
                
                
                previous = existing.toString();
                

                
                // build all possible paths for various lengths
                for(int i=0;i<tokens.size();i++) {
                    String term = tokens.get(i);
                    
		    // build a new path by adding a new term to the longest common substring between previous line and thisline
                    LinkedList<String> newpath = new LinkedList<String>(path);
                    newpath.add(term);
                   
		    // for newly discovered GO term, add the linkage between the term and the new functional paths to hashmap
                    try {
                        index.get(term).add(paths.size());
		    }
                    catch(NullPointerException ex) {
                        ArrayList<Integer> temp = new ArrayList<Integer>();
                        temp.add(paths.size());
                        index.put(term, temp);
                    }
                    
                    // add the new path to the collection
		    paths.add(newpath);
                    
                    path = newpath;
                }
                
            }
            
            in.close();
            
        } catch(Exception e) { System.err.println(e.getMessage()); }
        
        return paths;
    }

    
    // load prediction
    public static HashMap<String, ArrayList<String>> loadPredictions(String filename) {
        
        HashMap<String, ArrayList<String>> predicted = new HashMap<String, ArrayList<String>>();
        
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            
            String thisline;
            while((thisline=in.readLine())!=null) {
	
                StringTokenizer st = new StringTokenizer(thisline, "\t");
               
		String protein = st.nextToken().trim();
              
 
		ArrayList<String> predictedAnnotations = new ArrayList<String>();
                

		while(st.hasMoreTokens()) {

                    String next = st.nextToken().trim();
 
                    if(next.length() >= 10) {
                    
                        int divind = next.indexOf('(');
                        
                        String annot = next.substring(0, divind);
        
 
			/*
                        // use this part when having a threshold for confidence score
                        String scoreStr = next.substring(divind, next.length());
                        scoreStr = scoreStr.replace("(","").replace("*","").replace(")","").trim();
                        

                        double score = 0;
                        try { score = Double.parseDouble(scoreStr); }
                        catch(Exception e) { System.out.println("conversion to double failed."); } 
			*/


			// assumption: preiction set is sorted
			// just consider top 3
			if(!predictedAnnotations.contains(annot)) predictedAnnotations.add(annot);
			if(predictedAnnotations.size() == 3) break;
                        
                        
                    }
                }
                
                predicted.put(protein, predictedAnnotations);
 
            }
            
            in.close();
            
        } catch(Exception e) { System.err.println(e.getMessage()); }
        
        return predicted;
    }
    
    
    // load known annotation file
    public static HashMap<String, ArrayList<String>> loadKnownAnnotation(String filename) {
        
        HashMap<String, ArrayList<String>> annotations = new HashMap<String, ArrayList<String>>();
        
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            
            String thisline;
            while((thisline=in.readLine())!=null) {
                
                StringTokenizer st = new StringTokenizer(thisline, "\t");
                
                String protein = st.nextToken();
                protein = protein.substring(0, protein.length()-1);
                
                ArrayList<String> annotation = new ArrayList<String>();
                
                while(st.hasMoreTokens()) {
                    annotation.add(st.nextToken());
                }

                
                annotations.put(protein, annotation);
            }
            
            
        } catch(Exception e) { System.err.println(e.getMessage()); }
        
        return annotations;
    }

}
