package student;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.PriorityQueue;

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

public class SteinerTree {

    private class Path
    {
        public WeightedVertex start;
        public WeightedVertex end;
        public ArrayList<WeightedVertex> intermediates;

        public Path(WeightedVertex s, 
                    WeightedVertex e, 
                    ArrayList<WeightedVertex> i)
        {
            start = s;
            end = e;
            intermediates = i;
        }
    }

	public static int steinerTree(Graph g, ArrayList<Vertex> targets) 
    {
		return pathFinder(g, targets);
	}
	
	private static ArrayList<Edge> getPath(WeightedVertex vert){
		return null;
	}

    /*=========================================================================
     *
     * Method name: pathFinder
     *
     * Parameters: Graph g - the graph we are analyzing
     *             ArrayList<WeightedVertex> targets - the target verticies
     *
     * Returns: int - the cost of the final Stiener tree
     *
     * Description: Find the optimal Stiener tree between the target nodes 
     *              within this grpah.
     *
     *=======================================================================*/
    private static int pathFinder(Graph g, ArrayList<Vertex> targets)
    {
        /* First, we iteratively run Dijkstra's algorithm starting at each one
         * of the target nodes, and find the shortest path between any two 
         * Target nodes. */
       
		HashSet<Vertex> targetSet = new HashSet<Vertex>();
		targetSet.addAll(targets);

        // set the value field for each edge to be the same as the edge weight
        Iterator<Edge> itr = g.edgeIterator();

        Edge e;
        while (itr.hasNext())
        {
            e = itr.next();
            e.setValue(e.getWeight());
        }

        // find the shortest path between black nodes    
        WeightedVertex shortest = getShortestPath(targetSet, targets,g);

        // set the edge weights to zero on shortest path
        
        return 0;
    }

    /*==========================================================================
     *
     * Method name: getShortestPath
     *
     * Parameters: HashSet<Vertex> black_nodes - the set of target nodes
     *             ArrayList<Vertex> targets - the list or target nodes
     *
     * Returns: WeightedVertex - the shortest path between two black nodes
     *
     * Description: This function determines the shortest path between two 
     *              black nodes.
     *
     *=======================================================================*/
    private static WeightedVertex getShortestPath(HashSet<Vertex> black_nodes, 
                                                  ArrayList<Vertex> targets,
                                                  Graph g)
    {
        ArrayList<WeightedVertex> short_paths = new ArrayList<>();    

		for(Vertex target : targets)
        {
            ArrayList<WeightedVertex> paths = shortestPaths(g, target);

            ArrayList<WeightedVertex> options = 
                new ArrayList<WeightedVertex>();

            for (WeightedVertex path : paths)
            {
                if (targets.contains(path.getVert()) && 
                    !path.getVert().equals(target))
                {
                    options.add(path);
                }
            }

            Collections.sort(options);
            short_paths.add(options.get(0));
		}

        Collections.sort(short_paths);

        return short_paths.get(0);
    }

	/**
	 * This method takes in a graph and a starting vertex and runs an O(V lg E)
	 * implementation of Dijkstra's Shortest Paths algorithm.
	 * 
	 * @param g
	 *            An inputed Graph
	 * @param start
	 *            An arbitrary vertex in the graph to calculate the shortest
	 *            distances from.
	 * @return An ArrayList of WeightedVertex objects, each of which contain the
	 *         vertex to reach, the lowest cost to reach it, and it's parent in
	 *         the path to reach it.
	 */
	private static ArrayList<WeightedVertex> shortestPaths(Graph g, 
                                                           Vertex start) 
    {
		// Initialize storage that we'll need.
		ArrayList<WeightedVertex> outputs = new ArrayList<>();
		HashSet<Vertex> addedToQueue = new HashSet<Vertex>();
		PriorityQueue<WeightedVertex> heap = new PriorityQueue<>();

		// Use the heap object to keep track of what we've dealt with.
		// addedToQueue is a quick way of keeping track of what we've added over
		// the course of execution, even if it isn't in there right now.
		heap.add(new WeightedVertex(start, 0, null));
		addedToQueue.add(start);

		while (!heap.isEmpty()) 
        {
			// Pull a solution from the top of the heap.
			// This only works as the cheapest solutions rise to the top.
			WeightedVertex prior = heap.poll();
			outputs.add(prior);

			// Use the current solution to make the solution(s).
			Iterator<Edge> itr = prior.getVert().edgeIterator();
			while (itr.hasNext()) 
            {
				Edge currEdge = itr.next();
				Vertex vertToAdd = 
                    currEdge.getOppositeVertexOf(prior.getVert());

				// If this is a neighbor we haven't ever dealt with before
				// (We check in order to prevent re-adding a vertex we already
				// claimed was done.)
				if (!addedToQueue.contains(vertToAdd)) 
                {
					addedToQueue.add(vertToAdd);
					// The cost of this vertex must be equal to the cost of the
					// current path so far, plus the cost of the edge to get
					// from that path to here.
					heap.add(new WeightedVertex(vertToAdd, 
                                prior.getWeight() + currEdge.getValue(), 
                                prior));
				}
			}
		}
		return outputs;
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
			targetsRemaining--; // we found a target vertex
		if (targetsRemaining == 0)
			return 0; // all targets found, we are done

		// iterate over all edges out of this vertex
		for (Edge e : v) {
			Vertex newv = e.getOppositeVertexOf(v);
			if (newv.getValue() == 0) { // found an unreached vertex
				e.setMark(1); // we are considering
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
				} else {
					// unmark, this lead to nothing used
					e.setMark(0);
					e.setColor(Color.RED); // draw in red
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
