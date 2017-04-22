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
	public static ArrayList<Edge> getPath(WeightedVertex vert) {
		ArrayList<Edge> pathEdges = new ArrayList<>();
		WeightedVertex prior = vert.getPrior();
		Vertex start = vert.getVert();

		for (; prior != null; prior = prior.getPrior()) {
			// Find the edge going between the current edge and the prior.
			Iterator<Edge> itr = prior.getVert().iterator();
			while (itr.hasNext()) {
				Edge currEdge = itr.next();
				if (currEdge.getOppositeVertexOf(prior.getVert()).equals(start)) 
                {
					// Found a match.
					pathEdges.add(currEdge);
					break;
				}
			}
			start = prior.getVert();
		}
		return pathEdges;
	}
}
