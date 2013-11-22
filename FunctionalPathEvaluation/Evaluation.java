import java.io.*;
import java.util.*;

public class Evaluation {

    public static ArrayList<LinkedList<String>> fp;
    public static HashMap<String, ArrayList<Integer>> index;
    
    public static HashMap<String, ArrayList<String>> knownAnnotation;
    
    public static void main(String[] args) {
        
        if(args.length != 3) {
            System.out.println("wrong usage: consult README file");
            System.exit(1);
        }

	        
	        
        // directory containing prediction results
        String predictionDir = args[0];
        String outputDir = (predictionDir.charAt(predictionDir.length()-1)=='/')? predictionDir.substring(0, predictionDir.length()-2) + "eval/" : predictionDir + "_eval/";

        
        // load known annotations
        knownAnnotation = FileIO.loadKnownAnnotation(args[1]);
        
        
        // load functional paths
        // index contains functional path id that a functional annotation belongs to
        index = new HashMap<String, ArrayList<Integer>>();
        fp = FileIO.loadFunctionalPath(args[2], index);
        
                
        
        // collect names of input files
        File iDir = new File(predictionDir);
        File[] listOfFiles = iDir.listFiles();
        
        
        
        // create the output directory if not exist
        File oDir = new File(outputDir);
        if(!oDir.exists()) oDir.mkdirs();
        
        
        // Do evalution for all files in a given directory
        for(int i=0;i<listOfFiles.length;i++) {
            
            // name of input file and output file
            String inputfile = predictionDir + listOfFiles[i].getName();
            String outputfile = outputDir + listOfFiles[i].getName();
            

	    // load all the predicted GO annotations
            HashMap<String, ArrayList<String>> predicted = FileIO.loadPredictions(inputfile);


 	    // PRcalculation instance will calculate recall and precision for count and depth methods	
            PRcalculation c = new PRcalculation(knownAnnotation, predicted, fp, index);
            c.execute();
            

            // collect all the results
            double avgprec_count = c.getAveragedPrecisionCount();
            double avgrec_count = c.getAveragedRecallCount();
            
            HashMap<String, Double> precision_count = c.getPrecisionCount();
            HashMap<String, Double> recall_count = c.getRecallCount();
            
            double avgprec_depth = c.getAveragedPrecisionDepth();
            double avgrec_depth = c.getAveragedRecallDepth();
            
            HashMap<String, Double> precision_depth = c.getPrecisionDepth();
            HashMap<String, Double> recall_depth = c.getRecallDepth();
            
            
            
            /*
            // print summary to terminal
            System.out.println(listOfFiles[i].getName() + "\t" + avgprec_count + "\t" + avgrec_count + "\t" + avgprec_depth + "\t" + avgrec_depth);
            */
            
            
                               
            // print the result to files
            try {
                BufferedWriter outcount = new BufferedWriter(new FileWriter(outputfile.replace("prelist", "count")));
                BufferedWriter outdepth = new BufferedWriter(new FileWriter(outputfile.replace("prelist", "depth")));
                
                outcount.write("average precision: " + avgprec_count + "\taverage recall: " + avgrec_count + "\n");
                outdepth.write("average precision: " + avgprec_depth + "\taverage recall: " + avgrec_depth + "\n");
                                                             
                Iterator it = precision_count.entrySet().iterator();
                while(it.hasNext()) {
                    Map.Entry pair = (Map.Entry)it.next();
     
                    String protein = pair.getKey().toString();
                    
                    double prec = Double.parseDouble(pair.getValue().toString());
                    
                    double rec = 0;
                    try { rec = recall_count.get(protein).doubleValue(); }
                    catch(NullPointerException ex) { System.err.println("Count recall not found for: " + protein); } 
                
                    outcount.write(protein + "\t" + prec + "\t" + rec + "\n");
                    
                    prec = 0; rec = 0;
                    
                    try { prec = precision_depth.get(protein).doubleValue(); }
                    catch(NullPointerException ex) { System.err.println("Depth precision not found for: " + protein); } 
                    
                    try { rec = recall_depth.get(protein).doubleValue(); }
                    catch(NullPointerException ex) { System.err.println("Depth recall not found for: " + protein); } 
                    
                    outdepth.write(protein + "\t" + prec + "\t" + rec + "\n");
                }
                
                outcount.close();
                outdepth.close();
                
            } catch(Exception e) { System.err.println(e.getMessage()); }
            
            
        }
   
        
    }
    

}
