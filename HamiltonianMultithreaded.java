import java.util.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.lang.*;

class HamiltonianThread extends Thread {
    private int numVertices, threadNum;
    private boolean ans = true;

    public HamiltonianThread(int numVertices, int i) {
        this.numVertices = numVertices;
        this.threadNum = i;
    }

    @Override
    public void run() {
        HamiltonianCycle hamiltonian = new HamiltonianCycle();
        int[][] graph = HamiltonianMultithreaded.createMatrix(numVertices, threadNum);
        if (!hamiltonian.hamCycle(graph)) {
            List<Edge> edges;
            for (int vertexNum; vertexNum < numVertices; vertexNum++) {
                for (int i; i < graph.length; i++) {
                    for (int j; j < graph.length; j++) {
                        if (graph[i][j] == 1) {
                            edges.add(new Edge(i,j));
                        }
                        // Set number of vertices in the graph
                
                        // create a graph from edges
                        Graph g = new Graph(edges, numVertices);
                
                        // starting node
                        int start = vertexNum;
                
                        // add starting node to the path
                        List<Integer> path = new ArrayList<>();
                        path.add(start);
                
                        // mark start node as visited
                        boolean[] visited = new boolean[numVertices];
                        visited[start] = true;
                
                        if (!HamiltonianPaths.checkHamiltonianPaths(g, start, visited, path, numVertices)); {
                            ans = false;
                        }                    
                    }
                }
            }
        } else {
            ans = false;
        }
    }

    public boolean getAns() {
        return ans;
    }
}


public class HamiltonianMultithreaded {
    /**
     * Find the max of the elements of an array.
     *
     * @param numVertices integer to max
     * @return true/false
     * @throws InterruptedException shouldn't happen
     */
    public static boolean hamiltonian(int numVertices) throws InterruptedException {
        boolean ans = false;
        
        int numThreads = (int) Math.pow((double) 2, (double) (numVertices*numVertices-3*numVertices)/2);
        // Create and start threads.
        HamiltonianThread[] ts = new HamiltonianThread[numThreads];
        for (int i = 0; i < numThreads; i++) {
            ts[i] = new HamiltonianThread(numVertices, i);
            ts[i].start();
        }

        // Wait for the threads to finish and sum their results.
        for (int i = 0; i < numThreads; i++) {
            ts[i].join();
            if (ts[i].getAns() == true) {
                ans = true;
                System.out.println(createMatrix(numVertices, i));
            }
        }
        return ans;
    }
    public static int[][] createMatrix(int numVertices, int id) {
        int[][] graph = new boolean[numVertices][numVertices];
        for (int i = 0; i < numVertices; i++) {
            graph[i] = new boolean[numVertices];
        }
        for (int i = 0; i < numVertices - 1; i++) {
            graph[i][i+1] = true;
            graph[i+1][i] = true;
            //graph[0][numVertices-1] = false, graph[numVertices-1][0] = false
        }
    }
}


// https://www.geeksforgeeks.org/hamiltonian-cycle-backtracking-6/
/* Java program for solution of Hamiltonian Cycle problem 
using backtracking */
class HamiltonianCycle 
{ 
	final int V = 5; 
	int path[]; 

	/* A utility function to check if the vertex v can be 
	added at index 'pos'in the Hamiltonian Cycle 
	constructed so far (stored in 'path[]') */
	boolean isSafe(int v, int graph[][], int path[], int pos) 
	{ 
		/* Check if this vertex is an adjacent vertex of 
		the previously added vertex. */
		if (graph[path[pos - 1]][v] == 0) 
			return false; 

		/* Check if the vertex has already been included. 
		This step can be optimized by creating an array 
		of size V */
		for (int i = 0; i < pos; i++) 
			if (path[i] == v) 
				return false; 

		return true; 
	} 

	/* A recursive utility function to solve hamiltonian 
	cycle problem */
	boolean hamCycleUtil(int graph[][], int path[], int pos) 
	{ 
		/* base case: If all vertices are included in 
		Hamiltonian Cycle */
		if (pos == V) 
		{ 
			// And if there is an edge from the last included 
			// vertex to the first vertex 
			if (graph[path[pos - 1]][path[0]] == 1) 
				return true; 
			else
				return false; 
		} 

		// Try different vertices as a next candidate in 
		// Hamiltonian Cycle. We don't try for 0 as we 
		// included 0 as starting point in hamCycle() 
		for (int v = 1; v < V; v++) 
		{ 
			/* Check if this vertex can be added to Hamiltonian 
			Cycle */
			if (isSafe(v, graph, path, pos)) 
			{ 
				path[pos] = v; 

				/* recur to construct rest of the path */
				if (hamCycleUtil(graph, path, pos + 1) == true) 
					return true; 

				/* If adding vertex v doesn't lead to a solution, 
				then remove it */
				path[pos] = -1; 
			} 
		} 

		/* If no vertex can be added to Hamiltonian Cycle 
		constructed so far, then return false */
		return false; 
	} 

	/* This function solves the Hamiltonian Cycle problem using 
	Backtracking. It mainly uses hamCycleUtil() to solve the 
	problem. It returns false if there is no Hamiltonian Cycle 
	possible, otherwise return true and prints the path. 
	Please note that there may be more than one solutions, 
	this function prints one of the feasible solutions. */
	boolean hamCycle(int graph[][]) 
	{ 
		path = new int[V]; 
		for (int i = 0; i < V; i++) 
			path[i] = -1; 

		/* Let us put vertex 0 as the first vertex in the path. 
		If there is a Hamiltonian Cycle, then the path can be 
		started from any point of the cycle as the graph is 
		undirected */
		path[0] = 0; 
		if (hamCycleUtil(graph, path, 1) == false) 
		{ 
			System.out.println("\nSolution does not exist"); 
			return false; 
		} 

		printSolution(path); 
		return true; 
	} 

	/* A utility function to print solution */
	void printSolution(int path[]) 
	{ 
		System.out.println("Solution Exists: Following" + 
						" is one Hamiltonian Cycle"); 
		for (int i = 0; i < V; i++) 
			System.out.print(" " + path[i] + " "); 

		// Let us print the first vertex again to show the 
		// complete cycle 
		System.out.println(" " + path[0] + " "); 
	} 
/*
	// driver program to test above function 
	public static void main(String args[]) 
	{ 
		HamiltonianCycle hamiltonian = 
								new HamiltonianCycle(); 
		/* Let us create the following graph 
		(0)--(1)--(2) 
			| / \ | 
			| / \ | 
			| /	 \ | 
		(3)-------(4) 
		int graph1[][] = {{0, 1, 0, 1, 0}, 
			{1, 0, 1, 1, 1}, 
			{0, 1, 0, 0, 1}, 
			{1, 1, 0, 0, 1}, 
			{0, 1, 1, 1, 0}, 
		}; 

		// Print the solution 
		hamiltonian.hamCycle(graph1); 

		/* Let us create the following graph 
		(0)--(1)--(2) 
			| / \ | 
			| / \ | 
			| /	 \ | 
		(3)	 (4) 
		int graph2[][] = {{0, 1, 0, 1, 0}, 
			{1, 0, 1, 1, 1}, 
			{0, 1, 0, 0, 1}, 
			{1, 1, 0, 0, 0}, 
			{0, 1, 1, 0, 0}, 
		}; 

		// Print the solution 
		hamiltonian.hamCycle(graph2); 
	} */
} 
// This code is contributed by Abhishek Shankhadhar 
 
// data structure to store graph edges
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest) {
        this.source = source;
        this.dest = dest;
    }
}
 
// class to represent a graph object
class Graph
{
    // A List of Lists to represent an adjacency list
    List<List<Integer>> adjList = null;
 
    // Constructor
    Graph(List<Edge> edges, int N)
    {
        adjList = new ArrayList<>();
 
        for (int i = 0; i < N; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the undirected graph
        for (Edge edge: edges)
        {
            int src = edge.source;
            int dest = edge.dest;
 
            adjList.get(src).add(dest);
            adjList.get(dest).add(src);
        }
    }
}
 
class HamiltonianPaths
{
    public static boolean checkHamiltonianPaths(Graph g, int v, boolean[] visited, List<Integer> path, int N) {
        // if all the vertices are visited, then
        // hamiltonian path exists
        if (path.size() == N) {
            // print hamiltonian path
            return true;
        }
 
        // Check if every edge starting from vertex v leads
        // to a solution or not
        for (int w : g.adjList.get(v))
        {
            // process only unvisited vertices as hamiltonian
            // path visits each vertex exactly once
            if (!visited[w])
            {
                visited[w] = true;
                path.add(w);
 
                // check if adding vertex w to the path leads
                // to solution or not
                if (checkHamiltonianPaths(g, w, visited, path, N)) {
                    return true;
                }
 
                // Backtrack
                visited[w] = false;
                path.remove(path.size()-1);
            }
        }
    return false;
    }
/*
    public static void main(String[] args)
    {
        // List of graph edges as per above diagram
        List<Edge> edges = Arrays.asList(
                new Edge(0, 1), new Edge(0, 2), new Edge(0, 3),
                new Edge(1, 2), new Edge(1, 3), new Edge(2, 3)
        );
 
        // Set number of vertices in the graph
        final int N = 4;
 
        // create a graph from edges
        Graph g = new Graph(edges, N);
 
        // starting node
        int start = 0;
 
        // add starting node to the path
        List<Integer> path = new ArrayList<>();
        path.add(start);
 
        // mark start node as visited
        boolean[] visited = new boolean[N];
        visited[start] = true;
 
        printAllHamiltonianPaths(g, start, visited, path, N);
    } */

    public static void main(String[] args) {
        //System.out.println(HamiltonianMultithreaded.hamiltonian(4));
        HamiltonianCycle hamiltonian = new HamiltonianCycle();
        int[][] graph = new int[4][4];
        int[] list = new int[4];
        list[0]=0;
        list[1]=0;
        list[2]=0;
        list[3]=0;
        graph[0] = list;
        graph[1] = list;
        graph[2] = list;
        graph[3] = list;
        int numVertices = 4;
        if (!hamiltonian.hamCycle(graph)) {
            List<Edge> edges;
            for (int vertex; vertex < numVertices; vertex++) {
                for (int i; i < graph.length; i++) {
                    for (int j; j < graph.length; j++) {
                        if (graph[i][j] == 1) {
                            edges.add(new Edge(i,j));
                        }
                        // Set number of vertices in the graph
                
                        // create a graph from edges
                        Graph g = new Graph(edges, numVertices);
                
                        // starting node
                        int start = vertex;
                
                        // add starting node to the path
                        List<Integer> path = new ArrayList<>();
                        path.add(start);
                
                        // mark start node as visited
                        boolean[] visited = new boolean[numVertices];
                        visited[start] = true;
                
                        if (!checkHamiltonianPaths(g, start, visited, path, numVertices)); {
                            System.out.println("false");
                        }                    
                    }
                }
            }
        } else {
            System.out.println("false");
        }
        System.out.println("true");
    }
}