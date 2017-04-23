package student;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.HashMap;
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

public class SteinerTree 
{
    private static HashMap<Integer, Integer> old_weights;

	public static int steinerTree(Graph g, ArrayList<Vertex> targets) 
    {
        old_weights = new HashMap<Integer, Integer>();
		return pathFinder(g, targets);
	}

	/*
	 * =========================================================================
	 *
	 * Method name: pathFinder
	 *
	 * Parameters: Graph g - the graph we are analyzing
	 * ArrayList<WeightedVertex> targets - the target verticies
	 *
	 * Returns: int - the cost of the final Stiener tree
	 *
	 * Description: Find the optimal Stiener tree between the target nodes
	 * within this grpah.
	 *
	 * =======================================================================
	 */
	private static int pathFinder(Graph g, ArrayList<Vertex> targets) {
		/*
		 * First, we iteratively run Dijkstra's algorithm starting at each one
		 * of the target nodes, and find the shortest path between any two
		 * Target nodes.
		 */
		int sum = 0;
		HashSet<Vertex> targetSet = new HashSet<Vertex>();
		targetSet.addAll(targets);

		// set the value for each edge to be the same as the edge weight
		Iterator<Edge> itr = g.edgeIterator();

		Edge e;
		while (itr.hasNext()) {
			e = itr.next();
			e.setValue(e.getWeight());
		}

        // Find a base solution
		for (int i = 0; i < targets.size() - 1; i++) {
			// find the shortest path between black nodes
			WeightedVertex shortest = getShortestPath(targetSet, targets, g);
			// set the edge weights to zero on shortest path
			ArrayList<Edge> path_back = shortest.getPath();

			for (Edge edge : path_back) {
                old_weights.put(edge.getId(), 0);
				sum += edge.getValue();
				edge.setValue(0);
				edge.setMark(1);
			}
			SteinerTreeTester.show(g);
		}

        HashMap<Integer,Integer> clearing_search_space = 
            new HashMap<>(old_weights);

        // improve the solution by removing the path to one black node and
        // finding the minimum path to that node from the current solution
        for (Vertex target : targets)
            target.setValue(1);

        for (Vertex target : targets)
        {
            // remove the path to this vertex
            clearPath(target, clearing_search_space);
            WeightedVertex target_path = shortestPathTo(g, target, targetSet);
            System.out.println(target_path);
			ArrayList<Edge> path_back = target_path.getPath();
            for (Edge edge : path_back)
            {
                edge.setValue(0);
				edge.setMark(1);
            }
			SteinerTreeTester.show(g);
        }
/*
        System.out.println("---------------------------------------------");
        for (Vertex target : targets)
        {
            // remove the path to this vertex
            clearPath(target);
        }
*/
		return sum;
	}

    public static void clearPath(Vertex v, HashMap<Integer, Integer> space)
    {
        System.out.println("starting " + v);
        removeP(v, null, space);
    }

    public static boolean removeP(Vertex v, Vertex p, 
            HashMap<Integer, Integer> space)
    {
        Vertex next_v = v;
        for (Edge e : v)
        {
            if (space.containsKey(e.getId()) && 
                !e.getOppositeVertexOf(v).equals(p))
            {
                next_v = e.getOppositeVertexOf(v);
                System.out.println("visiting "+next_v);
                if (next_v.getMark() == 1)
                {
                    e.setValue(e.getWeight());
                    e.setMark(0);
                    space.remove(e.getId());
                    System.out.println("ending " + next_v);
                    return true;
                }

                if (removeP(next_v, v, space))
                {
                    e.setValue(e.getWeight());
                    return true;
                }
            }
        }

        return false;
    }

    public static boolean removePath(Vertex v, Vertex p)
    {
        // find the path to the nearest unvisited black node
        Vertex next_v;
        for (Edge e : v)
        {
            // find the next edge in the path
            if (e.getValue() == 0 && 
                !e.getOppositeVertexOf(v).equals(p))
            {
                next_v = e.getOppositeVertexOf(v);
                if (next_v.getMark() == 1 && 
                    next_v.getValue() < usedEdges(next_v))
                {
                    next_v.setValue(next_v.getValue()+1);
                    System.out.println("ending " + next_v);
                    return true;
                }
                else if (next_v.getMark() == 2 && 
                         next_v.getValue() < usedEdges(next_v))
                {
                    next_v.setValue(next_v.getValue()+1);
                    System.out.println("ending " + next_v);
                    return true;
                }
                else if (next_v.getMark() == 2)
                {
                    System.out.println("already visited " + next_v);
                    continue;
                }

                if (removePath(next_v, v))
                    return true;
            }
        }

        return false;
    }

	/*=========================================================================
	 * 
	 *
	 * Method name: getShortestPath
	 *
	 * Parameters: HashSet<Vertex> black_nodes - the set of target nodes
	 * ArrayList<Vertex> targets - the list or target nodes
	 *
	 * Returns: WeightedVertex - the shortest path between two black nodes
	 *
	 * Description: This function determines the shortest path between two black
	 * nodes.
	 *
	 *======================================================================*/
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
            WeightedVertex path = null;
            for (int i = 0; i < options.size(); i++)
            {
                if ((path = options.get(i)).getWeight() != 0){
                	// We're in the door.

                    // mark the edges in this path as zero
                    path.zeroPath();

                	if(!hasCycle(g, path.getVert(), null))
                    {
                        path.replacePath(old_weights);
                		break;
                	}
                    // set the path back to the way it was
                    path.replacePath(old_weights);
                }
            }

            short_paths.add(path);
		}

        Collections.sort(short_paths);

        WeightedVertex selected = short_paths.get(0);
        System.out.println("Selecting " + selected);
        return selected;
    }

    private static WeightedVertex shortestPathTo(Graph g, Vertex start,
                                                 HashSet<Vertex> targets)
    {
        ArrayList<WeightedVertex> paths = shortestPaths(g, start);

        Collections.sort(paths);

        for (WeightedVertex path : paths)
        {
            if (targets.contains(path.getVert()) && 
                !path.getVert().equals(start) && path.getWeight() > 0)
            {
                return path;
            }
        }

        return null;
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
	private static boolean hasCycle(Graph g, Vertex v, Vertex parent) {
		boolean cycle = innerHasCycle(g, v, parent);
		Iterator<Vertex> itr = g.vertexIterator();
		while (itr.hasNext()) {
			itr.next().setValue(0);
		}
		return cycle;
	}

    public static int usedEdges(Vertex v)
    {
        int count = 0;
        for (Edge e : v)
        {
            if (e.getValue() == 0)
                count++;
        }

        return count;
    }

    public static boolean vertexInPath(Vertex v)
    {
        for (Edge e : v)
        {
            if (e.getValue() == 0)
                return true;
        }

        return false;
    }

	private static boolean innerHasCycle(Graph g, Vertex v, Vertex parent) {
		// set value to indicated this vertex has been reached
		v.setValue(1);

		// iterate over all edges out of this vertex
		for (Edge e : v) {
			Vertex newv = e.getOppositeVertexOf(v);
			// If part of the selected path
			if (vertexInPath(newv) &&
                    !newv.equals(parent) && e.getValue() == 0 &&
                    (newv.getValue() == 1 || innerHasCycle(g, newv, v))) {
				return true;
			}
		}
		return false;
	}
}
