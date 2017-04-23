package student;

import graph.Vertex;
import graph.Edge;
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
	
	// Can have an increment of 0 - no increment.
	// 1 - found a new path, want to reflect that
	// -1 - that new path was bad, so go along and undo that offset.
	public static ArrayList<Edge> getPath(WeightedVertex vert, int increment) {
		ArrayList<Edge> pathEdges = new ArrayList<>();
		WeightedVertex prior = vert.getPrior();
		Vertex start = vert.getVert();
		
		vert.getVert().setMark(vert.getVert().getMark() + increment);
		for (; prior != null; prior = prior.getPrior()) {
			// Find the edge going between the current edge and the prior.
			Iterator<Edge> itr = prior.getVert().iterator();
			while (itr.hasNext()) {
				Edge currEdge = itr.next();
				if (currEdge.getOppositeVertexOf(prior.getVert()).equals(start)) 
                {
					// Found a match.
					prior.getVert().setMark(prior.getVert().getMark() + increment);
					pathEdges.add(currEdge);
					break;
				}
			}
			start = prior.getVert();
		}
		return pathEdges;
	}
	
	public static WeightedVertex getPathOrigin(WeightedVertex vert){
		if(vert != null){
			for(; vert.getPrior() != null; vert = vert.getPrior());
			return vert;
		}
		return null;
	}
}
