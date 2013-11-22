import java.io.*;
import java.util.*;

public class PRcalculation {

    // functional annotation: protein -> array list of functional annotations
    public HashMap<String, ArrayList<String>> known;
    public HashMap<String, ArrayList<String>> prediction;
    
    // functional path in a namespace (BP, CC or MF)
    public ArrayList<LinkedList<String>> paths;

    // hashmap to index of functional path in "paths" that contains the key "protein".
    public HashMap<String, ArrayList<Integer>> index;
    

    // precision and recall for counting method for each protein
    public HashMap<String, Double> precision_count;
    public HashMap<String, Double> recall_count;
    

    // precision and recall for depth method for each protein
    public HashMap<String, Double> precision_depth;
    public HashMap<String, Double> recall_depth;
   


    // constructor 
    public PRcalculation(HashMap<String, ArrayList<String>> known, HashMap<String, ArrayList<String>> prediction, ArrayList<LinkedList<String>> paths, HashMap<String, ArrayList<Integer>> index) {
        
        this.known = known;
        this.prediction = prediction;
        this.paths = paths;
        this.index = index;
        
	// the number of protein in prediction set and known set should be equal.
        assert(prediction.size() == known.size());

	// where the results will be stored        
        precision_count = new HashMap<String, Double>();
        recall_count = new HashMap<String, Double>();
        
        precision_depth = new HashMap<String, Double>();
        recall_depth = new HashMap<String, Double>();
    }

   

    // actual calculation
    public void execute() {
        
        int counter = 0;
        
	// iterate through proteins in known annotation set
        for(String protein : known.keySet()) {
            
            counter++;

            
            //testing
            //System.out.println(counter + "\t" + protein);
            

	    // for each protein, collect all the known and predicted functional annotations 
            ArrayList<String> knownAnnotations = known.get(protein);
            ArrayList<String> predictedAnnotations = prediction.get(protein);
            

            // indices of funtionalPaths that contain functional annotations in known and prediction set
            HashMap<Integer, Integer> fpIndexKnown = new HashMap<Integer, Integer>();
            HashMap<Integer, Integer> fpIndexPredicted = new HashMap<Integer, Integer>();
            




	    // fpKnown: collect indices of functional paths that cover the GO terms in known set
            for(int i=0;i<knownAnnotations.size();i++) {

	        // all the functional paths for a GO term (=knownAnnotations(i))
		ArrayList<Integer> fpaths = index.get(knownAnnotations.get(i));

                for(int j=0;j<fpaths.size();j++) {
                    Integer temp = fpaths.get(j);
                    if(!fpIndexKnown.containsKey(temp)) fpIndexKnown.put(temp, 1);
                }
            }
            

	    // fpPredicted: collect indices of functional paths that cover the GO terms in prediction set
            for(int i=0;i<predictedAnnotations.size();i++) {

		if(!index.containsKey(predictedAnnotations.get(i))) System.out.println(predictedAnnotations.get(i));

		// all the functional paths for a GO term (=predictedAnnotations(i))
		ArrayList<Integer> ppaths = index.get(predictedAnnotations.get(i));

                for(int j=0;j<ppaths.size();j++) {
                    Integer temp = ppaths.get(j);
                    if(!fpIndexPredicted.containsKey(temp)) fpIndexPredicted.put(temp, 1);
                }

            }
            
            

            // Reduce the number of functional paths by taking the longest possible ones to cover all annotations
	    // Remove the shorter functional paths if other functional paths can cover the proteins
	    // convert HashMap key set to ArrayList<Integer> to apply the consolidatePaths function
	    ArrayList<Integer> fpKnown = new ArrayList<Integer>(fpIndexKnown.keySet());
	    Collections.sort(fpKnown);
	    consolidatePaths(fpKnown);

	    ArrayList<Integer> fpPredicted = new ArrayList<Integer>(fpIndexPredicted.keySet());
	    Collections.sort(fpPredicted);
            consolidatePaths(fpPredicted);
            
            

            // The number of functional paths that were consolidated
            int fpKnownSize = fpKnown.size();
            int fpPredictedSize = fpPredicted.size();
            
            
            // store maximum level of overlaps between two lists for all possible pairs of lists
	    // matrix: for each pair of functioanl paths one from known and one from prediction, record the maximum level of overlaps
            int[][] mlevel = new int[fpKnownSize][fpPredictedSize];
            


	    
	    // maximum level of overlaps
            // pre_depth[x] = max mlevel[x][i] for 0<=i<fpKnownSize
	    int[] pre_depth = new int[fpKnownSize];

	    // annot_depth[x] = max mleve[i][x] for 0<=i<fpPredictedSize
	    int[] annot_depth = new int[fpPredictedSize];




            // record the height of each functional paths that we used for PR calculation
            int[] k1 = new int[fpKnownSize];
            int[] k2 = new int[fpPredictedSize];

            
            
	    // when we are comparing two paths, we are collecting all the GO terms in the functional path,
	    // to calcuate the precision and recall for "counting" method
	    // Collection all the GO terms from the functional paths will ensure that all the parent terms of annotated function will be included
            HashMap<String, Integer> overlaps = new HashMap<String, Integer>();
            HashMap<String, Integer> knownlist = new HashMap<String, Integer>();
            HashMap<String, Integer> predictionlist = new HashMap<String, Integer>();                                    
            
            
            for(int i=0;i<fpKnownSize;i++) {
                LinkedList<String> path1 = paths.get(fpKnown.get(i).intValue());
                
                k1[i] = path1.size();
                
                for(int j=0;j<fpPredictedSize;j++) {
                    LinkedList<String> path2 = paths.get(fpPredicted.get(j).intValue());
                    
                    k2[j] = path2.size();
                   

		    // compareLinkedList is find the maximum level where two paths overlap 
                    int olevel = compareLinkedList(path1, path2, overlaps, knownlist, predictionlist);


		    // getting pre_depth for each functional path in known set (maximum column value for each row)
                    if(pre_depth[i] < olevel) pre_depth[i] = olevel;
                   
		    // overlap level
                    mlevel[i][j] = olevel;
                }
            }
            

            /*
            //testing
            for(int i=0;i<paths.size();i++) {
                System.out.println(paths.get(i));
            }
            
            for(int i=0;i<fpKnownSize;i++) {
                System.out.print(fpKnown.get(i).intValue() + " ");
            }
            System.out.println();
            
            for(int i=0;i<fpPredictedSize;i++) {
                System.out.print(fpPredicted.get(i).intValue() + " ");
            }
            System.out.println();
 
            
            for(int i=0;i<mlevel.length;i++) {
                for(int j=0;j<mlevel[i].length;j++) {
                    System.out.print(mlevel[i][j] + " ");
                }
                System.out.println();
            }
            */
            

	    // calculating annot_depth for each functional path in prediction set (maximum row value for each column)
            for(int j=0;j<fpPredictedSize;j++) {
                int max = 0;
                for(int i=0;i<fpKnownSize;i++) {
                    if(mlevel[i][j] > max) max = mlevel[i][j];
                }
                annot_depth[j] = max;
            }
            

            
            /*
            //testing 
            System.out.println(overlaps.size() + "\t" + predictionlist.size() + "\t" + knownlist.size());
            if(knownlist.size() == 0 || predictionlist.size() == 0) {
                System.out.println(protein);
            }
            */
            


            // calculate precision and recall: "count" method
            double precision = (predictionlist.size()==0)? 0 : (double)overlaps.size()/(double)predictionlist.size();
            double recall = (knownlist.size()==0)? 0 : (double)overlaps.size()/(double)knownlist.size();
           

	    // for each protein, store the result for counting method 
            precision_count.put(protein, precision);
            recall_count.put(protein, recall);
            


            //testing
            //System.out.println("precision_count = " + precision + "\trecall_count = " + recall);

            
            
            // calculate precision and recall: "depth" method
            precision = 0; recall = 0;
            

	    // recall
            for(int i=0;i<fpKnownSize;i++) {
                int exp = k1[i]-pre_depth[i];
                
                recall += (1/Math.pow(4, exp));
            }
            
            recall = (fpKnownSize==0)? 0 : recall/(double)fpKnownSize;
            

	    // precision
            for(int j=0;j<fpPredictedSize;j++) {
                int exp = k2[j]-annot_depth[j];
                
                precision += (1/Math.pow(4, exp));
            }
            
            precision = (fpPredictedSize==0)? 0 : precision/(double)fpPredictedSize;
            


	    // record the results
            precision_depth.put(protein, precision);
            recall_depth.put(protein, recall);
            

            //testing
            //System.out.println("precision_depth = " + precision + "\trecall_depth = " + recall);

            
        }
    
        
    }

    
    
    // return the maximum level of overlapped part of two lists
    private int compareLinkedList(LinkedList<String> l1, LinkedList<String> l2, HashMap<String, Integer> overlaps, HashMap<String, Integer> elements1, HashMap<String, Integer> elements2) {
        
        int len1 = l1.size();
        int len2 = l2.size();
        
        int minlen = (len1<=len2)? len1 : len2;
        
        int i = 0;
        while(i < minlen) {
            String temp1 = l1.get(i);
            String temp2 = l2.get(i);
            
            if(!elements1.containsKey(temp1)) elements1.put(temp1, 1);
            if(!elements2.containsKey(temp2)) elements2.put(temp2, 1);
            
            if(temp1.equals(temp2)) { 
                i++;
                if(!overlaps.containsKey(temp1)) overlaps.put(temp1, 1);
            }
            else break;
        }
        
        for(int j=i;j<len1;j++) {
            String temp = l1.get(j);
            if(!elements1.containsKey(temp)) elements1.put(temp, 1);
        }
        
        for(int j=i;j<len2;j++) {
            String temp = l2.get(j);
            if(!elements2.containsKey(temp)) elements2.put(temp, 1);
        }
        
        return i;
    }

    
   
    // Assumption: the paths are sorted by branches and from the shortest to the longest 
    // minimize the number of functional paths to cover all annotations
    // iterate through the paths and compare the paths next to each other, and keep only the longest ones
    private void consolidatePaths(ArrayList<Integer> indices) {
        
        int x=0;
        while(x < indices.size()-1) {
            LinkedList<String> list1 = paths.get(indices.get(x).intValue());
            LinkedList<String> list2 = paths.get(indices.get(x+1).intValue());
            
            if(list1.size() < list2.size()) {
                int count_match = 0;
                for(int i=0;i<list1.size();i++) {
                    if(list1.get(i).equals(list2.get(i))) count_match++;
                }
                if(count_match == list1.size()) indices.remove(x);
                else x++;
            }
            else x++;
        }
    }
    


    // calculate the average of the precision or recall over all proteins
    public double average(HashMap<String, Double> scores) {
        double sum = 0;
        
        Iterator it = scores.entrySet().iterator();
        
        while(it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            sum += Double.parseDouble(pair.getValue().toString());
        }
        
        return sum/(double)scores.size();
    }
    

    
    // return averaged scores
    public double getAveragedPrecisionCount() {
        return average(this.precision_count);
    }
    public double getAveragedRecallCount() {
        return average(this.recall_count);
    }

    public double getAveragedPrecisionDepth() {
        return average(this.precision_depth);
    }
    public double getAveragedRecallDepth() {
        return average(this.recall_depth);
    }


    
    // get functions
    public HashMap<String, Double> getPrecisionCount() { return precision_count; }
    public HashMap<String, Double> getRecallCount() { return recall_count; }
    
    public HashMap<String, Double> getPrecisionDepth() { return precision_depth; }
    public HashMap<String, Double> getRecallDepth() { return recall_depth; }
}
