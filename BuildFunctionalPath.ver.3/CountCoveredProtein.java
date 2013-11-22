import java.io.*;
import java.util.*;

public class CountCoveredProtein {

    public CountCoveredProtein() { }

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
    
}
