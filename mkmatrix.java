import java.util.*; 
import java.lang.*; 
import java.io.*;

class Graph{ 
    static class Edge implements Comparable<Edge>{ 
        int src, dest, weight; 
        public int compareTo(Edge compareEdge) 
        { 
            return this.weight-compareEdge.weight; 
        } 
    }; 
    static class subset{ 
        int parent, rank; 
    }; 
    static int V, E;    // V-> no. of vertices & E->no.of edges 
    static Edge edge[]; // collection of all edges 
    Graph(int v, int e) 
    { 
        V = v; 
        E = e; 
        edge = new Edge[E]; 
        for (int i=0; i<e; ++i) 
            edge[i] = new Edge(); 
    } 
    static int find(subset subsets[], int i){ 
        // find root and make root as parent of i (path compression) 
        if (subsets[i].parent != i) 
            subsets[i].parent = find(subsets, subsets[i].parent); 
  
        return subsets[i].parent; 
    } 
    static void Union(subset subsets[], int x, int y){ 
        int xroot = find(subsets, x); 
        int yroot = find(subsets, y); 
  
        // Attach smaller rank tree under root of high rank tree 
        // (Union by Rank) 
        if (subsets[xroot].rank < subsets[yroot].rank) 
            subsets[xroot].parent = yroot; 
        else if (subsets[xroot].rank > subsets[yroot].rank) 
            subsets[yroot].parent = xroot; 
  
        // If ranks are same, then make one as root and increment 
        // its rank by one 
        else{ 
            subsets[yroot].parent = xroot; 
            subsets[xroot].rank++; 
        } 
    } 
  
    // The main function to construct MST using Kruskal's algorithm 
    static void KruskalMST(){ 
        Edge result[] = new Edge[V];  // Tnis will store the resultant MST 
        int e = 0;  // An index variable, used for result[] 
        int i = 0;  // An index variable, used for sorted edges 
        for (i=0; i<V; ++i) result[i] = new Edge();
        // Step 1:  Sort all the edges in non-decreasing order of their 
        // weight.  If we are not allowed to change the given graph, we 
        // can create a copy of array of edges 
        Arrays.sort(edge); 
  
        // Allocate memory for creating V ssubsets 
        subset subsets[] = new subset[V]; 
        for(i=0; i<V; ++i) 
            subsets[i]=new subset(); 
  
        // Create V subsets with single elements 
        for (int v = 0; v < V; ++v) 
        { 
            subsets[v].parent = v; 
            subsets[v].rank = 0; 
        } 
  
        i = 0;  // Index used to pick next edge 
  
        // Number of edges to be taken is equal to V-1 
        while (e < V - 1) 
        { 
            // Step 2: Pick the smallest edge. And increment  
            // the index for next iteration 
            Edge next_edge = new Edge(); 
            next_edge = edge[i++]; 
  
            int x = find(subsets, next_edge.src); 
            int y = find(subsets, next_edge.dest); 
  
            // If including this edge does't cause cycle, 
            // include it in result and increment the index  
            // of result for next edge 
            if (x != y)
            { 
                result[e++] = next_edge; 
                Union(subsets, x, y); 
            } 
            // Else discard the next_edge 
        } 
  
        // print the contents of result[] to display 
        // the built MST 
        System.out.println("Following are the edges in " +  
                                     "the constructed MST"); 
        for (i = 0; i < e; ++i) 
            System.out.println(result[i].src+" -- " +  
                   result[i].dest+" == " + result[i].weight); 
    }
}
public class mkmatrix{
	public static int distance(String a, String b) {
        a = a.toLowerCase();
        b = b.toLowerCase();
        int [] costs = new int [b.length() + 1];
        for (int j = 0; j < costs.length; j++)
            costs[j] = j;
        for (int i = 1; i <= a.length(); i++) {
            costs[0] = i;
            int nw = i - 1;
            for (int j = 1; j <= b.length(); j++) {
                int cj = Math.min(1 + Math.min(costs[j], costs[j - 1]), a.charAt(i - 1) == b.charAt(j - 1) ? nw : nw + 1);
                nw = costs[j];
                costs[j] = cj;
            }
        }
        return costs[b.length()];
    }
    public static int br(String a, String b){
    	int[][] FKP;
    	if(a.length() > b.length()){
    		String c = a;
    		a = b;
    		b = c;
    	}
    	int m = a.length();
    	int n = b.length();
    	int ZERO_K = n; //initialized with length of larger string
    	FKP = new int[2*n+1][n+2];
    	for(int k = -ZERO_K;k < 0;k++){
    		int p = -k-1;
    		FKP[k + ZERO_K][p + 1] = Math.abs(k) - 1;
    		FKP[k + ZERO_K][p] = -Integer.MAX_VALUE;
    	}
    	FKP[ZERO_K][0] = -1;
    	for(int k = 1; k <= ZERO_K; k++){
    		int p = k - 1;
            FKP[k + ZERO_K][p + 1] = -1;
            FKP[k + ZERO_K][p] = -Integer.MAX_VALUE;
    	}
    	int p = n - m - 1;
    	do {
            p++;

            for (int i = (p - (n-m))/2; i >= 1; i--) {
                f(a, b, FKP, ZERO_K, n-m+i, p-i);
            }

            for (int i = (n-m+p)/2; i >= 1; i--) {
                f(a, b, FKP, ZERO_K, n-m-i, p-i);
            }

            f(a, b, FKP, ZERO_K, n - m, p);
        } while (FKP[(n - m) + ZERO_K][p] != m);
        return p - 1;
    }

    public static void f(String x, String y, int[][] FKP, int ZERO_K, int k, int p){
    	int t = Math.max(FKP[k + ZERO_K][p] + 1, FKP[k - 1 + ZERO_K][p]);
    	t=Math.max(t,FKP[k + 1 + ZERO_K][p] + 1);

        while (t < Math.min(x.length(), y.length() - k) && x.charAt(t) == y.charAt(t + k)) {
            t++;
        }

        FKP[k + ZERO_K][p + 1] = t;
    }

    public static void main(String[] args){
    	String [] data = { "atcgatatggaattaaagccactgtagaaaagcg", "atcgtgggttggtaaactgaaccgacgtaactga", "atcgtgggttggtaaactgaaccgacgtagcagc", "gtctggctatcagagagaactttttgatatcatgggatg", "aggttctagattcaactgactggcatgct", "gtctggctatcagagagaactttttgatatcga" };
    	int[][] matrix = new int[data.length][data.length];
    	int [] arr_d = new int[data.length]; //store min distance for each sequence
    	int [] arr_v = new int[data.length]; //store corresponding sequance
    	HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();//distances->corresponding index
        HashMap<Integer,Integer> bmap = new HashMap<Integer,Integer>();//min_d->index of arr_d
    	for(int i = 0; i<data.length;i++){//Create a matrix
    		int min_d = Integer.MAX_VALUE;
    		for(int j = 0;j<data.length;j++){
    	 		matrix[i][j] = br(data[i],data[j]);
    	 		if(matrix[i][j]!=0){
    	 			min_d = Math.min(min_d,matrix[i][j]);
    	 			map.put(matrix[i][j],j);
    	 		}
    	 	}
    	 	arr_d[i]= min_d;
    	 	arr_v[i]=map.get(min_d);
    	}
        int V = data.length;
        int E = V*(V-1); //upper bound for number of edit distances to be calculated
        Graph graph = new Graph(V,E);
        for(int i = 0; i < data.length;i++){
            for(int j = 0; j < data.length;j++){
                if(i > j){
                    graph.edge[i*j].src = i;
                    graph.edge[i*j].dest = j;
                    graph.edge[i*j].weight = matrix[i][j];
                }
            }
        }
        
        graph.KruskalMST();
        
    	for(int i = 0; i<data.length;i++){
    	 	for(int j = 0;j<data.length;j++){
    	 		System.out.printf("%6d",matrix[i][j]);
    	 	}
    	 	System.out.println();
    	}
    	for(int i = 0; i < data.length;i++){
    		System.out.printf("%6d",arr_d[i]);
    	}
    	System.out.println();
    	for(int i = 0; i < data.length;i++){
    		System.out.printf("%6d",arr_v[i]);
    	}
    	System.out.println();
    }
}