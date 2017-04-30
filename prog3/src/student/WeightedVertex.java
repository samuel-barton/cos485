package student;

import graph.Vertex;
import graph.Edge;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;

public class WeightedVertex implements Comparable<WeightedVertex>{
	private Vertex vert;
	private WeightedVertex prior;
	private int weight;
	public WeightedVertex(Vertex vert, int weight, WeightedVertex priorVertex){
		this.setWeight(weight);
		this.setVert(vert);
		this.setPrior(priorVertex);
	}	
	public Vertex getVert() {
		return vert;
	}
	public void setVert(Vertex vert) {
		this.vert = vert;
	}
	public int getWeight() {
		return weight;
	}
	public void setWeight(int weight) {
		this.weight = weight;
	}
	@Override
	public int compareTo(WeightedVertex o) {
		return weight -o.weight;
	}
	public String toString(){
		return "(" + vert.getLabel() + ", " + weight + ". Came from " + prior + ")";
	}
	public WeightedVertex getPrior() {
		return prior;
	}
	public void setPrior(WeightedVertex prior) {
		this.prior = prior;
	}

    /*===========================================================================
     *
     * Method name: getPath
     * 
     * Parameters: none
     * 
     * Returns: ArrayList<Edge> - the path back from the current vertes
     *
     * Description: This method returns the edges in the path back to the parent
     *              black node.
     *
     *=========================================================================*/
    public ArrayList<Edge> getPath()
    {
        ArrayList<Edge> path_back = new ArrayList<Edge>();
		Vertex start = this.vert;
		WeightedVertex p = this.prior;
		for (; p != null; p = p.prior)
        {
			Iterator<Edge> itr = p.getVert().iterator();
			while (itr.hasNext()) 
            {
                Edge curr_edge = itr.next();
				if (curr_edge.getOppositeVertexOf(p.getVert()).equals(start)) 
                {
                    // edge is in path
                    path_back.add(curr_edge);
                    break;
                }
            }
            start = p.vert;
        }
    
        return path_back;    
    }

    /*===========================================================================
     *
     * Method name: zeroPath
     * 
     * Parameters: none
     * 
     * Returns: void
     *
     * Description: This method zeros the values of the edges in the path back to
     *              the parent black node.
     *
     *=========================================================================*/
    public void zeroPath()
    {
		Vertex start = this.vert;
		WeightedVertex p = this.prior;
		for (; p != null; p = p.prior)
        {
			Iterator<Edge> itr = p.getVert().iterator();
			while (itr.hasNext()) 
            {
                Edge curr_edge = itr.next();
				if (curr_edge.getOppositeVertexOf(p.getVert()).equals(start)) 
                {
                    // edge is in path
                    curr_edge.setValue(0);
                    break;
                }
            }
            start = p.vert;
        }
    }

    /*===========================================================================
     *
     * Method name: replacePath
     * 
     * Parameters: HasMap<int,int> old_weights - edges which have been modified
     * 
     * Returns: void
     *
     * Description: This method removes the changes made by the path represented
     *              by this WeightedVertex using old_weights as a guide.
     *
     *=========================================================================*/
    public void replacePath(HashMap<Integer, Integer> old_weights)
    {
		Vertex start = this.vert;
		WeightedVertex p = this.prior;
		for (; p != null; p = p.prior)
        {
			Iterator<Edge> itr = p.getVert().iterator();
			while (itr.hasNext()) 
            {
                Edge curr_edge = itr.next();
				if (curr_edge.getOppositeVertexOf(p.getVert()).equals(start)) 
                {
                    // edge is in path
                    if (!old_weights.containsKey(curr_edge.getId()))
                        curr_edge.setValue(curr_edge.getWeight());
                    break;
                }
            }
            start = p.vert;
        }
    }

    public ArrayList<Vertex> getEnds()
    {
        ArrayList<Vertex> ends = new ArrayList<Vertex>();

        ends.add(vert);

        WeightedVertex p = prior;

        while (p.prior != null)
            p = p.prior;

        ends.add(p.vert);
        
        return ends;
    }
	
	public static WeightedVertex getPathOrigin(WeightedVertex vert){
		if(vert != null){
			for(; vert.getPrior() != null; vert = vert.getPrior());
			return vert;
		}
		return null;
	}
}
