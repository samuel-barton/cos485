package student;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Set;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Comparator;
import java.util.PriorityQueue;

import graph.*;

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
    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    public static int steinerTree(Graph g, ArrayList<Vertex> targets) 
    {
        old_weights = new HashMap<Integer, Integer>();
        int sum = pathFinder(g, targets);
        return sum;
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
        // we could just use array accessing
        HashSet<Vertex> targetSet = new HashSet<Vertex>();
        targetSet.addAll(targets);


        // set the value for each edge to be the same as 
        // the edge weight
        Iterator<Edge> itr = g.edgeIterator();

        Edge e;
        while (itr.hasNext()) {
            e = itr.next();
            e.setValue(e.getWeight());
        }

        // Find a base solution
        boolean[] canVisit = new boolean[g.numVertices()];
        for(int i = 0; i < g.numVertices(); i++)
            canVisit[i] = true;
        for (int i = 0; i < targets.size() - 1; i++) {
            // find the shortest path between black nodes
            WeightedVertex shortest = getShortestPath(targetSet, 
                                                      targets, g, 
                                                      canVisit);

            // set the edge weights to zero on shortest path
            ArrayList<Edge> path_back = shortest.getPath();
            for (Edge edge : path_back) 
            {
                old_weights.put(edge.getId(), 0);
                sum += edge.getValue();
                edge.setValue(0);
            }
        }
        Iterator<Vertex> itrv = g.vertexIterator();
        Vertex[] vert = new Vertex[g.numVertices()];
        while(itrv.hasNext())
        {
            Vertex v = itrv.next();
            vert[v.getId()] = v;
        }
       
        Edge[] edgeArray = new Edge[g.numEdges()];
        Iterator<Edge> ii = g.edgeIterator();
        while(ii.hasNext())
        {
            Edge f = ii.next();
            edgeArray[f.getId()] = f;
        } 

        // set up for the iterative improvement
        int[] numberOfUsedEdgesTester = new int[g.numVertices()];
        int[] numberOfUsedEdges = new int[g.numVertices()];
        int bestSolution = sum;
        int tempSolution = 0;
        CmpVertex cmp = new CmpVertex();
        HashSet<Vertex> usedVertices = new HashSet<Vertex>();
        updateNumberOfConnections(numberOfUsedEdges, g, usedVertices);
        
        HashSet<Edge> usedEdges = new HashSet<Edge>();
        Iterator<Edge> itre = g.edgeIterator();
        while(itre.hasNext())
        {
            Edge e2 = itre.next();
            if(e2.getValue() == 0)
                usedEdges.add(e2);
        }

        boolean[] edgeMask = new boolean[g.numEdges()];
        boolean[] edgeMaskTemp = new boolean[g.numEdges()];
        Iterator<Edge> eI = g.edgeIterator();
        while(eI.hasNext())
        {
            Edge f = eI.next();
            if(f.getValue() == 0)
                edgeMask[f.getId()] = true;
        }

        for(int i = 0; i < g.numEdges(); i++)
        {
            if(edgeArray[i].getValue() == 0)
                edgeMask[i] = true;
            else
                edgeMask[i] = false;
        }
 
        // here is the iterative improvement bit prepare for loops
        while(tempSolution < bestSolution)
        {
            // tempSolution is the sum of the edges in the current graph, 
            // this will change as we add and remove things
            tempSolution = bestSolution;
            Iterator<Vertex> usedVItr = usedVertices.iterator();

            while(usedVItr.hasNext())
            {
                // Iterate through all the intermediates we used last time 
                // to get a solution. Throw them all out.
                Vertex v0 = usedVItr.next();
                if(targetSet.contains(v0))//don't remove target nodes
                    continue;

                // Else ...
                // this is the removal process, by throwing this flag we say 
                // that dijkstra cannot visit here
                canVisit[v0.getId()] = false;
                
                // set all of the removed vertices edges to Integer.MAX_VALUE 
                // so that we can recognize which ones need to be reset later
                for(Edge eg : v0)
                {
                    if(eg.getValue() == 0)
                        eg.setValue(Integer.MAX_VALUE);
                }
                
                int numLinks = numberOfUsedEdges[v0.getId()];
                
                for (int i = 0; i < targets.size() - 1; i++)
                {
                    // find the shortest path between black nodes
                    WeightedVertex shortest = getShortestPath(targetSet, 
                                                              targets, g, 
                                                              canVisit);
                    
                    // set the edge weights to zero on shortest path
                    ArrayList<Edge> path_back = shortest.getPath();
                    for (Edge edge : path_back) 
                    {
                        old_weights.put(edge.getId(), 0);
                        edge.setValue(0);
                    }
                }
                // we can figure out that tempSolution another 
                // way, lets do so
                tempSolution = 0;
                Iterator<Edge> tempSol = g.edgeIterator();
                while(tempSol.hasNext())
                {
                    Edge eeeee = tempSol.next();
                    if(eeeee.getValue() == 0)
                        tempSolution+=eeeee.getWeight();
                }
                
                // adjust the values of the numberOfUsedEdges to
                // reflect the number of links to that node
                // (determined by how many edges are zero)
                Iterator<Vertex> connectionCounter = g.vertexIterator();
                while(connectionCounter.hasNext())
                {
                    int count = 0;
                    Vertex vName = connectionCounter.next();
                    for(Edge currEdge : vName)
                        if(currEdge.getValue() == 0)
                            count++;
                    numberOfUsedEdgesTester[vName.getId()] = count;
                }
                
                // Clear cycles
                for (Vertex v_0 : targets)
                    while (resolveCycle(g,v_0));

                // envoke trimmer to modify our solution by removing 
                // white leaves as these will never be part of an 
                // optimal solution
                for(int i = 0; i < g.numVertices(); i++)
                    tempSolution = betterTrimmer(tempSolution, 
                                                 numberOfUsedEdgesTester, 
                                                 i, vert, targetSet);
                
                // allow our removed node to be considered in the future
                canVisit[v0.getId()] = true;
                
                // if we've bettered ourself then update best solution
                System.out.println(tempSolution + " <? " + bestSolution);
                if(tempSolution < bestSolution)
                {
                    bestSolution = tempSolution;
                    tempSolution = 0;
                    System.arraycopy(numberOfUsedEdgesTester, 0, 
                                     numberOfUsedEdges, 0, 
                                     g.numVertices());
                    for(Edge eee : v0)
                    // set all of the edges that were formally part 
                    // of a solution with the removed node to their weight
                        if(eee.getValue() == Integer.MAX_VALUE)
                            //since they are no longer part of a solution
                            eee.setValue(eee.getWeight());
                    // stop this pass of the vertex iteration to 
                    // try our new set(also the iterator may have 
                    // been invalidated by the removal of
                    break;
                    // successive white nodes that formed a chain 
                    // that was not needed in our solution
                }
                // Else, if we have reached this code then the answer we 
                // came up with by removing the vertex v0 did not 
                // improve our solution, so reset everything.
                System.arraycopy(numberOfUsedEdges, 0, 
                                 numberOfUsedEdgesTester, 0, 
                                 g.numVertices());
                
                int testingNumEdgesReset = 0;
                for(Edge currEdge : v0)
                {
                    if(currEdge.getValue() == Integer.MAX_VALUE)
                    {
                        currEdge.setValue(0);
                        testingNumEdgesReset++;
                    }
                }
            }     
        }
        
        for(int i = 0; i < g.numEdges(); i++)
            if(edgeMask[i])
                edgeArray[i].setMark(1); 

        return sum;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    public static int betterTrimmer(int sum, int[] numberOfConnections, 
                                    int index, Vertex[] v, 
                                    HashSet<Vertex> targetNodes)
    {
        if(numberOfConnections[index] == 1 && !targetNodes.contains(v[index]))
        {
            numberOfConnections[index] = 0;
            for(Edge e : v[index])
            {
                if(e.getValue() == 0)
                {
                    e.setValue(e.getWeight());
                    sum = betterTrimmer(sum-e.getWeight(),numberOfConnections, 
                                        e.getOppositeVertexOf(v[index]).getId(), 
                                        v, targetNodes);
                    break;
                }
            }
        }
        return sum;
    }
    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    public static void modifiedDFS(HashSet<Vertex> visited, Vertex child)
    {
        for(Edge e : child)
        {
            if(e.getValue() == 0)
            {
                Vertex canidate = e.getOppositeVertexOf(child);
                if(!visited.contains(canidate))
                {
                    visited.add(canidate);
                    modifiedDFS(visited, canidate);
                }
            }
        }
    }
    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    public static void updateNumberOfConnections(int[] numberOfUsedEdges, 
                                                 Graph g, 
                                                 Set<Vertex> UsedVertices)
    {//I was suppose to use this but I kept just writing this method instead
        Iterator<Vertex> vItr = g.vertexIterator();
        while(vItr.hasNext())
        {
            Vertex v = vItr.next();
            int count = 0;
            for(Edge e : v)
                if(e.getValue() == 0)
                    count++;
            if(count!=0)
            {
                numberOfUsedEdges[v.getId()] = count;
                UsedVertices.add(v);
            }
        }
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
            Graph g, boolean[] canVisit)
    {
        ArrayList<WeightedVertex> short_paths = new ArrayList<>();    
        for(Vertex target : targets)
        {
            ArrayList<WeightedVertex> paths = shortestPaths(g, target, canVisit);
            ArrayList<WeightedVertex> options = new ArrayList<WeightedVertex>();
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
        return selected;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static WeightedVertex shortestPathTo(Graph g, 
                                                 Vertex start, 
                                                 HashSet<Vertex> targets, 
                                                 boolean[] canVisit)
    {
        ArrayList<WeightedVertex> paths = shortestPaths(g, start, canVisit);
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
    /*make modified version of this tomorrow*/
    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static ArrayList<WeightedVertex> shortestPaths(Graph g, 
                                                           Vertex start, 
                                                           boolean[] canVisit) 
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
                if(!canVisit[vertToAdd.getId()])//the update to dijkstra
                    addedToQueue.add(vertToAdd);
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

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static boolean resolveCycle(Graph g, Vertex v)
    {
        ArrayList<Vertex> verticies = findCycle(g,v);

        if (verticies == null)
            return false;

        ArrayList<Edge> edges = getEdges(verticies);

        Collections.sort(edges, new CmpEdge());

        Edge e = edges.get(edges.size()-1);
        e.setValue(e.getWeight());
        e.setMark(0);
        System.out.println("We removed "+e);

        return true;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/

    private static ArrayList<Edge> getEdges(ArrayList<Vertex> verticies)
    {
        HashSet<Vertex> magical_verticies = new HashSet<>(verticies);
        ArrayList<Edge> edges = new ArrayList<Edge>();
        for (Vertex v : verticies)
        {
            for (Edge e : v)
                if (e.getValue() == 0 && magical_verticies.contains(e.getOppositeVertexOf(v)))
                    edges.add(e);
        }

        return edges;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static ArrayList<Vertex> findCycle(Graph g, Vertex v)
    {
        ArrayList<Vertex> cycle = new ArrayList<Vertex>(); 
        ArrayList<Vertex> actual_cycle = null;

        if (findInnerCycle(g,v,null,cycle))
        {
            Vertex cycle_bound = null;
            int cycle_index = 0;
            HashSet<Vertex> keys = new HashSet<>();
            for (Vertex v0 : cycle)
            {
                if (!keys.contains(v0))
                    keys.add(v0);
                else
                {
                    cycle_bound = v0; 
                    break;
                }
                cycle_index++;
            }

            if (cycle_bound == null)
                return null;   

            actual_cycle = new ArrayList<Vertex>();
            for (int i = cycle_index-1; i > 0 && !cycle.get(i).equals(cycle_bound); i--)
                actual_cycle.add(cycle.get(i));
            
            actual_cycle.add(cycle_bound); 
        }

        Iterator<Vertex> itr = g.vertexIterator();
        while (itr.hasNext()) 
        {
            itr.next().setValue(0);
        }

        return actual_cycle;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static boolean findInnerCycle(Graph g, Vertex v, Vertex parent, ArrayList<Vertex> cycle)
    {
        // set value to indicated this vertex has been reached
        v.setValue(1);

        // iterate over all edges out of this vertex
        for (Edge e : v) {
            Vertex newv = e.getOppositeVertexOf(v);
            // If part of the selected path
            if (vertexInPath(newv) &&
                    !newv.equals(parent) && e.getValue() == 0 &&
                    (newv.getValue() == 1 || findInnerCycle(g, newv, v, cycle))) 
            {
                cycle.add(newv);
                return true;
            }
        }
        return false;
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
    private static boolean hasCycle(Graph g, Vertex v, Vertex parent) 
    {
        boolean cycle = innerHasCycle(g, v, parent);
        Iterator<Vertex> itr = g.vertexIterator();
        while (itr.hasNext()) {
            itr.next().setValue(0);
        }
        return cycle;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
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

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    public static boolean vertexInPath(Vertex v)
    {
        for (Edge e : v)
        {
            if (e.getValue() == 0)
                return true;
        }

        return false;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static boolean innerHasCycle(Graph g, Vertex v, Vertex parent) 
    {
        // set value to indicated this vertex has been reached
        v.setValue(1);

        // iterate over all edges out of this vertex
        for (Edge e : v) {
            Vertex newv = e.getOppositeVertexOf(v);
            // If part of the selected path
            if (vertexInPath(newv) &&
                    !newv.equals(parent) && e.getValue() == 0 &&
                    (newv.getValue() == 1 || innerHasCycle(g, newv, v))) 
            {
                return true;
            }
        }
        return false;
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static class CmpEdge implements Comparator<Edge>
    {
        public int compare(Edge e1, Edge e2)
        {
            return e1.getWeight() - e2.getWeight();
        }
    }

    /*===========================================================================
     *
     * Method name: 
     * 
     * Parameters: 
     * 
     * Returns: 
     *
     * Description: 
     *
     *=========================================================================*/
    private static class CmpVertex implements Comparator<Vertex>
    {
        public int compare(Vertex v1, Vertex v2)
        {
                    return v1.getId() - v2.getId();
        }
    }
}

