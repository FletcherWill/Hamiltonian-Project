import java.util.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.lang.*;
import java.io.*;

/**
 * This class represents an object that checks a single graph to see if it is n-strung with no Hamiltonian cycle.
 * Graphs are passed in via their unique threadID. Each instance of this class is ran in a unique thread.
 */
class HamiltonianThread extends Thread {
    private int numVertices, threadNum;
    private boolean ans = true;

    /**
     * Initializer for HamiltonianThread class.
     * @param numVertices
     * @param i
     */
    public HamiltonianThread(int numVertices, int i) {
        this.numVertices = numVertices;
        this.threadNum = i;
    }

    /**
     * This is the function called when the thread is started. This function takes the thread ID, 
     * uses createMatrix function to generate the corresponding graph and then checks it for 
     * numVertices-strungness and having no Hamiltonian cycle.
     */
    @Override
    public void run() {
        HamiltonianCycle hamiltonian = new HamiltonianCycle(this.numVertices);
        int[][] graph = HamiltonianMultithreaded.createMatrix(this.numVertices, threadNum);
        if (!hamiltonian.hamCycle(graph)) {
            List<Edge> edges = new ArrayList<Edge>();
            for (int i = 0; i < graph.length; i++) {
                for (int j = 0; j < graph.length; j++) {
                    if (graph[i][j] == 1) {
                        edges.add(new Edge(i,j));
                    }
                }
            }
            for (int vertexNum = 0; vertexNum < numVertices; vertexNum++) {
                Graph g = new Graph(edges, numVertices);
                int start = vertexNum;
                List<Integer> path = new ArrayList<>();
                path.add(start);
                boolean[] visited = new boolean[numVertices];
                visited[start] = true;
                if (!HamiltonianPaths.checkHamiltonianPaths(g, start, visited, path, numVertices)) {
                    this.ans = false;
                }               
            }
        } else {
            this.ans = false;
        }
    }

    /**
     * Returns the answer set by the run method.
     * @return ans
     */
    public boolean getAns() {
        return this.ans;
    }
}

/**
 * This class tests all graphs of a given number of vertices, n, for the property of n-strungness but no Hamiltonian cycle.
 */
public class HamiltonianMultithreaded {
    /**
     * Tests all graphs of a certain number of vertices. Returns true if it finds an n-strung graph with no Hamiltonian cycle,
     * and prints the graph. Otherwise, returns false. Each graph is tested in its own thread for a faster runtime.
     *
     * @param numVertices integer
     * @return true/false
     * @throws InterruptedException shouldn't happen
     */
    public static boolean hamiltonian(int numVertices) throws InterruptedException {
        boolean ans = false;
        
        int numThreads = (int) Math.pow((double) 2, (double) (numVertices*numVertices-3*numVertices)/2);
        HamiltonianThread[] ts = new HamiltonianThread[numThreads];
        for (int i = 0; i < numThreads; i++) {
            ts[i] = new HamiltonianThread(numVertices, i);
            ts[i].start();
        }

        int counterExample = -1;
        for (int i = 0; i < numThreads; i++) {
            ts[i].join();
            if (ts[i].getAns() == true) {
                ans = true;
                System.out.println(Arrays.deepToString(createMatrix(numVertices, i)));
                counterExample = i;
            }
        }
        if (ans) {
            System.out.println(Arrays.deepToString(createMatrix(numVertices, counterExample)));
            System.out.println(counterExample);
        }
        return ans;
    }
    
    /**
     * Helper function for the createMatrix function. This cleans up the binary numbers to have the correct amount
     * of digits by adding 0's to the left of the string.
     * @param inputString
     * @param length
     * @return
     */
    public static String padLeftZeros(String inputString, int length) {
		if (inputString.length() >= length) {
			return inputString;
		}
		StringBuilder sb = new StringBuilder();
		while (sb.length() < length - inputString.length()) {
			sb.append('0');
		}
		sb.append(inputString);

		return sb.toString();
	  }
    
    /**
     * Takes in an ID and number of vertices to have in the graph. Generates a graph by converting the ID to a
     * binary number, assigning each entry of the adjacency matrix to a digit of the binary number. For
     * optimization, entries along the diagonal are always set to 0, entries to the immediate left and right of
     * the diagonal are set to 1 creating a single Hamiltonian path, since we know there must be at least one path,
     * so we assume it is this one without loss of generality. Finally the top right corner and bottom left corner
     * are set to 0, so that this Hamiltonian path does not form a cycle. Finally, we only change entries to the
     * right of the diagonal, since this is an undirected graph, so we reflect edges over the diagonal to finish
     * our adjacency matrix.
     * @param numVertices integer
     * @param id integer
     * @return graph integer matrix
     */
    public static int[][] createMatrix(int numVertices, int id) {
        int[][] graph = new int[numVertices][numVertices];
		for (int i = 0; i < numVertices; i++) {
			graph[i] = new int[numVertices];
		}
	
		for (int i = 0; i < numVertices - 1; i++) {
			graph[i][i+1] = 1;
		  graph[i+1][i] = 1;
        }
	
		int stringLength = (numVertices*numVertices-3*numVertices)/2;
		String binary = Integer.toBinaryString(id);
		binary = padLeftZeros(binary, stringLength);
		int insert = stringLength - numVertices + 3;
		binary = binary.substring(0, insert) + "0" + binary.substring(insert);
		stringLength++;

		int k = 0;
		for (int i = 0; i < numVertices - 2; i++ ) {
		  for (int j = i + 2; j < numVertices; j++) {
			graph[i][j] = (int)(binary.charAt(stringLength - k - 1)) - 48;
			graph[j][i] = (int)(binary.charAt(stringLength - k - 1)) - 48;
			k++;
		  }
		}
		return graph;
    }

    public static void main(String[] args) {
        int numVertices = 8;
        try {
            System.out.println(HamiltonianMultithreaded.hamiltonian(numVertices));
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    } 
}

/** 
 * Class for checking a graph for Hamiltonian cycles using the backtracking algorithm.
 * This is a modified version of the code found at:
 * https://www.geeksforgeeks.org/hamiltonian-cycle-backtracking-6/
 */
class HamiltonianCycle 
{ 
	int V; 
    int path[];
    
    public HamiltonianCycle(int V) {
        this.V = V;
    }

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
			//System.out.println("\nSolution does not exist"); 
			return false; 
		} 
		return true; 
	} 
} 

/**
 * data structure to store graph edges
 */ 
class Edge
{
    int source, dest;
 
    /**
     * Initializing function
     * @param source
     * @param dest
     */
    public Edge(int source, int dest) {
        this.source = source;
        this.dest = dest;
    }
}
 
/**
 * Class to represent a graph object.
 */
class Graph
{
    List<List<Integer>> adjList = null;
 
    /**
     * Initializer
     * @param edges
     * @param N
     */
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

/**
 * Class that checks for a Hamiltonian path starting at a certain vertex using the backtracking algorithm.
 */
class HamiltonianPaths
{
    /**
     * Recursive algorithm that tries to add unvisited vertices to the path and then recursively call itself
     * on each possible addition to the current path.
     * @param g
     * @param v
     * @param visited
     * @param path
     * @param N
     * @return true/false
     */
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
}
