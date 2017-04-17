package student;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Random;

import graph.*;
import steinerTree.SteinerTreeTester;

/* 
 * This Student class is meant to contain your algorithm.
 * You should implement the static method:
 * 
 *   steinerTree - which finds a good Steiner Tree on the graph 
 *                 
 *   You do not need to find the optimal solution, but shorter is 
 *   better!
 *   It should set the mark field on edges that are part of your tree.
 *   It should return the sum of the edge weights of the selected edges.
 *   
 * The inputs are:
 *   1. Graph object, which has:
 *      an ArrayList of all vertices - use graph.vertexIterator()
 *      an ArrayList of all edges - use graph.edgeIterator()
 *      each vertex has an ArrayList of its edges - use 
 *      vertex.edgeIterator()
 *      see the documentation for: Graph, Vertex, and Edge for more 
 *      details
 *   2. An ArrayList of vertices that are the targeted vertices for 
 *      inclusion
 *      in your Steiner tree. The mark fields are also already set in 
 *      the graph for these vertices.  
 */ 
public class SteinerTree
{
    /*=================================================================
     *
     * Method Name: bartMont
     *
     * Parameters:  Graph g - The graph we are working with
     *              ArrayList<Vertex> targets - the target verticies
     * 
     * Returns: int - the total cost of our Steiner Tree
     *
     * Description: The algorithm works as follows:
     *
     *      1.  Calculate the neighbor edge cost sums for each 
     *          non-target vertex.
     *      2.  Find all direct path's between K-nodes
     *      3.  Add the N-K node with the minimum sum to K
     *      4.  Check for connections to K nodes
     *      5.  Repeat 3-4 until we generate a cycle or solution
     *      6.  If solution better than previous, update solution
     *      7.  If we have a cycle, pick the least expensive path back
     *          to the last visited K-node
     *      8.  Repeat 3-7 until all nodes and edges have been visited
     *
     *===============================================================*/
    private static int bartMont(Graph g, ArrayList<Vertex> targets)
    {
        HashMap<Vertex, Integer> costs = new HashMap<Vertex,Integer>();

        // Store the sums in a minHeap by sum
        // FIX BELOW
        Iterator<Vertex> itr = g.vertexIterator();

        // 1.

        while (itr.hasNext())
        {
            Vertex v = itr.next();
            if (!isTarget(v, targets)) 
            {
                int sum = 0;
                for (Edge e : v)
                    sum += e.getWeight();

                costs.put(v,sum);
                System.out.println("vertex "+v+", cost "+sum);
            }    
        }

        return 0;
    }

    private static boolean isTarget(Vertex v, ArrayList<Vertex> targets)
    {
        for (Vertex t : targets)
            if (t.equals(v))
                return true;

        return false;
    }


	// Simple example routine that just does a depth first search until 
    // it reaches all of the target vertices.
	public static int steinerTree(Graph g, ArrayList<Vertex> targets)
	{	
		// sort each vertex's edges shortest first, this should help a 
        // little
		g.sortVertexEdgeLists(new Graph.CompareEdgesSmallestFirst());

		// start at first target vertex
		Vertex v = targets.get(0);
        // search based on the number of remaining targets
		dfs(g, v, targets.size()); 

		// go add up the weights of all the marked edges
		int length = 0;
        // iterate over all edges in the graph
		Iterator<Edge> itr = g.edgeIterator();  
		while (itr.hasNext()) {
			Edge e = itr.next();
			if (e.getMark() == 1)
				length += e.getWeight();
		}

		return length;
	}

	// Recursively depth first search the graph until all targets are 
    // reached. As it searches, it tags vertices as reached by setting 
    // their value field to 1. It sets the mark field for edges that 
    // were used to reach new target vertices. We can tell this 
    // occurred when we return from a recursive search of an edge and 
	// the number of targets remaining has decreased.
	//
	// We exit the dfs early once all targets have been found.
	//
	// returns: the number of targets still remaining 
	public static int dfs(Graph g, Vertex v, int targetsRemaining) {
        // set value to indicated this vertex has been reached
		v.setValue(1); 
		if (v.getMark() == 1)
			targetsRemaining--;    // we found a target vertex
		if (targetsRemaining == 0)
			return 0;				// all targets found, we are done
		
		// iterate over all edges out of this vertex
		for (Edge e : v) {
			Vertex newv = e.getOppositeVertexOf(v);
			if (newv.getValue() == 0) { // found an unreached vertex
				e.setMark(1);              // we are considering
                // color it green for animation)
				e.setColor(Color.GREEN);     
				SteinerTreeTester.show(g);
                // recursively search
				int newRemaining = dfs(g, newv, targetsRemaining); 
                // did this edge lead to any new targets
				if (newRemaining < targetsRemaining) { 
					targetsRemaining = newRemaining;
                    // mark this edge as part of solution
					e.setMark(1);                 
                    // for animation show the graph at this point
					SteinerTreeTester.show(g);    
				}
				else {
                    // unmark, this lead to nothing used
					e.setMark(0);                 
					e.setColor(Color.RED);        // draw in red
					SteinerTreeTester.show(g);
				}
				if (targetsRemaining == 0)
                    // all targets found, we are done
					return 0;				
			}
		}

		return targetsRemaining;
	}
}
