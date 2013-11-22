import java.util.*;

public class FunctionalPath {

	private String category;
	private ArrayList<String> funcPaths;
	private ArrayList<GOTerm> nodes;
	
	private int root;
	private HashMap<Integer, Integer> terminals;
	

	private boolean verbose = false;


	public FunctionalPath(String c, ArrayList<GOTerm> nodes) {
		this.category = c;
		this.funcPaths = new ArrayList<String>();
		
		this.nodes = nodes;
	
		findExtremes();
	}


	public void isVerbose(boolean b) {
		this.verbose = b;
	}


	public ArrayList<String> getFuncPaths() {
		return this.funcPaths;
	}


	// find root and terminals
	private void findExtremes() {
		this.terminals = new HashMap<Integer, Integer>();
		
		for(int i=0;i<this.nodes.size();i++) {
			// if a node has no parent, it is a root
			if(this.nodes.get(i).getNumParents() == 0) {
				this.root = i;
			}
			else if(this.nodes.get(i).getNumChildren() == 0) {
				this.terminals.put(i,1);
			}
		}
	}



	// build functional paths without informative nodes
	public void execute() {
		this.build(false);
	}


	public void execute(int threshold) {
		this.carryProteinsToTop();
		this.whichIsInformative(threshold);
		this.build(true);
	}



	// covered proteins //////////////////////////////////////////////////////////////////////////////////
	
	// pull out proteins that are covered by children
	public void carryProteinsToTop() {
		
		for(Integer key : this.terminals.keySet()) {
			
			int termind = key.intValue();
			carry(termind);
			
		}
		
	}

	
	// recursively gather proteins
	private void carry(int terminalindex) {
		
		if(terminalindex == this.root) {
			return;
		}
		else {
			HashMap<Integer, Integer> parents = this.nodes.get(terminalindex).getParents();	
			
			for(Integer key : parents.keySet()) {
				int parentind = key.intValue();
				
				this.nodes.get(parentind).addProteins(this.nodes.get(terminalindex).getProteins());
				carry(parentind);
			}
			
		}
	
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
    // resetInformativeness of each node ///////////////////////////////////////////////////////////////////////

	public void whichIsInformative(int threshold) {
		
		for(int i=0;i<this.nodes.size();i++) {
			this.nodes.get(i).resetInformative(threshold, this.nodes);
		}
		
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
	

	// functional paths //////////////////////////////////////////////////////////////////////////////////


	// build functional paths
	public void build(boolean identifyInformative) {
		this.funcPaths = findPaths(this.root, identifyInformative);
	}



	// recurrence
	private ArrayList<String> findPaths(int rootindex, boolean identifyInformative) {
        
        int termlen = 13;
        
        
        String isuffix = "";
        String nsuffix = "";
        
        
        if(identifyInformative) {
	        isuffix = "(i)"; 
	        nsuffix = "(n)";
        }
       
  
        GOTerm rootnode = this.nodes.get(rootindex);
        String rootterm = rootnode.getID();

	
		if(identifyInformative) {
		
			rootterm = rootterm + (rootnode.isInformative()?isuffix:nsuffix);
	
			//testing ////////////////////////////////////////
			if(this.verbose) rootterm = rootterm + ":" + rootnode.getNumProteins();
			/////////////////////////////////////////////////
		}
	

		ArrayList<String> paths = new ArrayList<String>();
	
        
        // base case: when rootindex is terminal
        if(this.terminals.containsKey(rootindex)) {
	        paths.add(rootterm);
	        return paths;
        }
        else {
	        
	         HashMap<Integer, Integer> children = rootnode.getChildren();
	         
	         for(Integer key : children.keySet()) {
		         
		         int childindex = key.intValue();
		         String childterm = this.nodes.get(childindex).getID();
		         
		         ArrayList<String> subpaths = findPaths(childindex, identifyInformative);
		         
		         for(int j=0;j<subpaths.size();j++) {
			         paths.add(rootterm + "\t" + subpaths.get(j));
		         }
		         
	         }
	        
	         return paths;
        }
        
        
    }
        
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////    



}