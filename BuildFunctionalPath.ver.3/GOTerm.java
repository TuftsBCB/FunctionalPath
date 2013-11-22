import java.util.*;


public class GOTerm implements Comparable<GOTerm> {

	private String id;
    private String namespace;
    private boolean informative;
    private HashMap<String, Integer> proteins;


	private HashMap<Integer, Integer> children;
	private HashMap<Integer, Integer> parents;    
    

	// Constructor
	public GOTerm(String id, String namespace) {
	
		this.id = id;
		this.namespace = namespace;
		
	    this.informative = false;
        
        this.children = new HashMap<Integer, Integer>();
        this.parents = new HashMap<Integer, Integer>();
        
        this.proteins = new HashMap<String, Integer>();
	
	}







	// informativeness ////////////////////////////////////////////////////////////////////////////////////////
	
	
	// there is no public function to set informative
	// informative should be manipulated within the object
	public boolean isInformative() {
        return this.informative;
    }


    // ver3: 2 coditions should be met
    // 1) a term should cover at least "threshold" proteins
    // 2) a term should not cover the same number of proteins as its children 
    public void resetInformative(int threshold, ArrayList<GOTerm> collection) {

	
		// if >= 1 of children has same number of proteins, then this node turns to be non-informative
		boolean sameNumProteinsWithChildren = false;
		for(Integer key : this.children.keySet()) {
			int ind = key.intValue();
			if(collection.get(ind).getNumProteins() == this.getNumProteins()) sameNumProteinsWithChildren = true;
		}
			

    	if(this.getNumProteins() > threshold && !sameNumProteinsWithChildren) { this.informative = true; }


/* testing
        /////////////////////////////////////////
         if(this.getID().contains("16043")) {
                 System.out.println("#protein = " + this.getNumProteins());
                 System.out.print("childrens = ");
                 for(int x=0;x<this.children.size();x++) {
                         System.out.print(this.children.get(x).getID() + "(" + this.children.get(x).getNumProteins() + "), ");
                 }
                 System.out.println("sameNumProteinsWithChildren: " + sameNumProteinsWithChildren);
		System.out.println("is this informative? " + this.informative);
		 System.out.println();
                 System.exit(1);
         }
         ////////////////////////////
*/
	
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////



    
    // protein ////////////////////////////////////////////////////////////////////////////////////////////////

	// get number of proteins that this goterm covers
	public int getNumProteins() {
		return this.proteins.size();
	}


	// get proteins that are annotated by this term
    public HashMap<String, Integer> getProteins() {
    	return this.proteins;
    }
    
    
    public void addProtein(String p) {
    	if(!this.proteins.containsKey(p)) this.proteins.put(p,1);
    }
    

    public void addProteins(HashMap<String, Integer> proteins) {
    	for(String protein : proteins.keySet()) {
	    	this.addProtein(protein);
    	}
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // relationship //////////////////////////////////////////////////////////////////////////////////////////////

    
    // add relationship information
    // add a child
    public void addChild(int ind) {
    	if(!this.children.containsKey(ind)) this.children.put(ind, 1);
    }
    
    public void addChildren(ArrayList<Integer> children) {
		for(int i=0;i<children.size();i++) {
			this.addChild(children.get(i).intValue());
		}
	}

	// add a parent
	public void addParent(int ind) {
		if(!this.parents.containsKey(ind)) this.parents.put(ind, 1);
	}
	
    public void addParents(ArrayList<Integer> parents) {
		for(int i=0;i<parents.size();i++) {
			this.addParent(parents.get(i).intValue());
		}
	}



	// check the relationship
	public boolean isChild(int index) {
		return (this.children.containsKey(index));
	}

	public boolean isParent(int index) {
		return (this.parents.containsKey(index));
	}

	
    // return children or parents of this node
    public HashMap<Integer, Integer> getChildren() {
        return this.children;
    }

	public HashMap<Integer, Integer> getParents() {
		return this.parents;
	}


	// get number of parents or children of this node
	public int getNumChildren() {
		return this.children.size();
	}

      public int getNumParents() {
        return this.parents.size();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


	

    // metadata ///////////////////////////////////////////////////////////////////////////////////////////////////
    
	// set namespace
    public void setNamespace(String namespace) {
        this.namespace = namespace;
    }

    // return name of this node
    public String getID() {
        return this.id;
    }


    // get namespace of this node
    public String getNamespace() {
        return this.namespace;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // compareTO function
    public int compareTo(GOTerm other) {
        int comparison = this.namespace.compareTo(other.getNamespace());
        if(comparison != 0) return comparison;
        else return this.id.compareTo(other.getID());
    }



}