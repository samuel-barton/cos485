package student;

import graph.Vertex;

public class WeightedVertex implements Comparable<WeightedVertex>{
	private Vertex vert, prior;
	private int weight;
	public WeightedVertex(Vertex vert, int weight, Vertex priorVertex){
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
	public Vertex getPrior() {
		return prior;
	}
	public void setPrior(Vertex prior) {
		this.prior = prior;
	}
}