
#pragma once


#ifndef __DGGRAPH_H__
#define __DGGRAPH_H__


#include <stack>
#include <vector>

using namespace std;


class DirGraph_Node{
public:
	int node_index;                   //Node index (unique)
	int *edges;                       //An array of the edges incident to the node
	int nedges;
	int levels;                       //the level during DFS
	int parent;                       //store the parent of current nodes for backward tracking to build the SCC
	int sscomp_index;
	unsigned char visited;
	int global_index;                 //Node index in the global 12/03/2009 by Qingqing
	bool rev_visited;

	DirGraph_Node()
	{
		edges = NULL;
		nedges = 0;
	}
};

class DirGraph_NodeList{
public:
	DirGraph_Node **dirnodes;
	int ndirnodes;
	int curMaxNumENodes;

	// The similar list operations
	DirGraph_NodeList(int initsize = 1000) //construction
	{
		//dirnodes = (DirGraph_Node **)malloc(sizeof(DirGraph_Node *)*initsize);
		dirnodes = new DirGraph_Node *[initsize];
		curMaxNumENodes = initsize;
		ndirnodes = 0;

		if(dirnodes == NULL)
		{
			char rout[256], var[256];
			sprintf(rout, "%s", "DirGraph_NodeList Constructor");
			sprintf(var, "%s", "dirnodes");

			//write_mem_error(rout, var, 0);
			curMaxNumENodes = 0;
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			dirnodes[i] = NULL;
	} 

	~DirGraph_NodeList()
	{
		if(dirnodes != NULL)
			//free(dirnodes);
		{
		//  we may need to get inside to delete each element
			for (int i=0; i<curMaxNumENodes; i++)
				delete dirnodes[i];
			delete [] dirnodes;
		}

	}

	void finalize()
	{
		if(dirnodes != NULL)
			//free(dirnodes);
		{
		//  we may need to get inside to delete each element
			printf("releasing graph node list...\n");
			for (int i=0; i<curMaxNumENodes; i++)
			{
				delete dirnodes[i];
			}
			delete [] dirnodes;

			dirnodes = NULL;
		}
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(DirGraph_Node *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		dirnodes[ndirnodes] = s;
		//copyElem(s, polist[nporbits]);
		ndirnodes++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		ndirnodes --;
		return true;
	} 

	inline void copy_Elem(DirGraph_Node *s, DirGraph_Node *d)
	{
		d->node_index = s->node_index;
		d->levels = s->levels;
		d->parent = s->parent;
		d->sscomp_index = s->sscomp_index;
		d->visited = s->visited;
		d->nedges = s->nedges;

		d->edges = new int[d->nedges];
		int i;
		for(i = 0; i < d->nedges; i++)
			d->edges[i] = s->edges[i];
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(DirGraph_Node *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < ndirnodes; i++)
		{
			if(dirnodes[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < ndirnodes-1; i++)
		{
			//we need a copy function
			copy_Elem(dirnodes[i], dirnodes[i+1]);
		}

		ndirnodes--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(ndirnodes == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(ndirnodes == curMaxNumENodes) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		DirGraph_Node **temp = dirnodes;
		//dirnodes = (DirGraph_Node **) malloc(sizeof(DirGraph_Node *) * (curMaxNumENodes + step));
		dirnodes = new DirGraph_Node *[curMaxNumENodes + step];  // modified on 07/05/2010
		if( dirnodes == NULL)
		{
			//fail
			char rout[256], var[256];
			sprintf(rout, "%s", "DirGraph_NodeList::extend");
			sprintf(var, "%s", "dirnodes");

			//write_mem_error(rout, var, 0);
			exit(-1);

			dirnodes = temp;
			return false;
		}

		int i;

		for(i = 0; i < curMaxNumENodes; i++)
			dirnodes[i] = temp[i];

		for(i=curMaxNumENodes; i<curMaxNumENodes+step; i++)
			dirnodes[i]=NULL;

		curMaxNumENodes += step;

		//free(temp);
		delete [] temp;

		return true;
	}

	inline void reset()
	{
		ndirnodes = 0;
	}

}; //end of DirGraph_NodeList class



/* Graph edge */
class Graph_Edge{
public:
	bool cancel;
	bool visited;
	int edge_index;                   //Edge index (unique)
	int node_index1, node_index2;     //Indices of the two nodes this edge connects
	//unsigned int traj;                          //the trajectory index of the corresponding separatrix 

	/*the variables to represent the region corresponds to this edge 08/29/2007*/
	int *triangles;
	int ntris;

	Graph_Edge()
	{
		triangles = NULL;
		ntris = 0;
	}

	~Graph_Edge()
	{
		if (triangles != NULL)
			delete [] triangles;
		triangles = NULL;
	}
};







class Graph_EdgeList{
public:
	Graph_Edge **edges;
	int nedges;
	int curMaxNumGedges;

	// The similar list operations as VertexList
	// The similar list operations
	Graph_EdgeList(int initsize = 1000) //construction
	{
		//edges = (Graph_Edge **)malloc(sizeof(Graph_Edge *)*initsize);
		edges = new Graph_Edge *[initsize]; // modified on 07/09/2010
 		curMaxNumGedges = initsize;
		nedges = 0;

		if(edges == NULL)
		{
			char rout[256], var[256];
			sprintf(rout, "%s", "Graph_EdgeList constructor");
			sprintf(var, "%s", "edges");

			//write_mem_error(rout, var, 0);
			curMaxNumGedges = 0;
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			edges[i] = NULL;
	
	} 

	~Graph_EdgeList()
	{
		if(edges != NULL)
		{
			int i;
			for(i = 0; i < curMaxNumGedges/*curMaxNumGedges*/; i++)
			{
				if(edges[i]!=NULL)
					//free(edges[i]);
					delete edges[i];
			}
			free(edges);
			delete [] edges;
		}
	}

	void finalize()
	{
		if(edges != NULL)
		{
			int i;
			printf("releasing graph edge list...\n");
			for(i = 0; i < curMaxNumGedges/*curMaxNumGedges*/; i++)
			{
				if(edges[i]!=NULL)
					//free(edges[i]);
					delete edges[i];
				edges[i] = NULL;
			}
			//free(edges);
			delete [] edges;
			edges = NULL;
		}
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(Graph_Edge *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		edges[nedges] = s;
		//copyElem(s, polist[nporbits]);
		nedges++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nedges --;
		return true;
	} 

	inline void copy_Elem(Graph_Edge *s, Graph_Edge *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(Graph_Edge *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nedges; i++)
		{
			if(edges[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nedges-1; i++)
		{
			//we need a copy function
			copy_Elem(edges[i], edges[i+1]);
		}

		nedges--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nedges == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nedges == curMaxNumGedges) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		Graph_Edge **temp = edges;
		//edges = (Graph_Edge **) malloc(sizeof(Graph_Edge *) * (curMaxNumGedges + step));
		edges = new Graph_Edge *[curMaxNumGedges + step];
		if( edges == NULL)
		{
			//fail
			char rout[256], var[256];
			sprintf(rout, "%s", "Graph_Edge::extend");
			sprintf(var, "%s", "edges");

			//write_mem_error(rout, var, 1);
			//exit(-1);

			edges = temp;
			return false;
		}

		int i;
		for(i = 0; i < curMaxNumGedges; i++)
			edges[i] = temp[i];
		
		for(i = curMaxNumGedges; i < curMaxNumGedges+step; i++)
			edges[i] = NULL;

		curMaxNumGedges += step;

		//free(temp);
		delete [] temp;

		return true;
	}

	inline void reset()
	{
		nedges = 0;
	}

}; //end of Graph_EdgeList class


/**Directed graph for encoding the dynamics of the flow **/
class DirGraph{
public:
	DirGraph_NodeList *nlist;
	Graph_EdgeList *elist;

	int *cur_nodes_order; /*use for finding strongly connected components*/
	int cur_sccnode_index;
	int num_sccomps;
	
	/*int *dfs_stack;
	int curMaxDFSStack;
	int top_dfsstack;*/
	stack<int> dfs_stack;



	/**Member functions**/
	DirGraph(int nnodes = 0, int nedges = 0)
	{
		//FILE *fp;
		//fp = fopen("detect_porbit_cooling.txt", "a");
		//fp = fopen("detect_porbit_swirl.txt", "a");
		//fp = fopen("detect_porbit.txt", "a");
		//fprintf(fp, "Start allocating nlist...\n");
		//fclose(fp);

		if(nnodes == 0)
			nlist = NULL;
		else
		{
			nlist = new DirGraph_NodeList(nnodes);
			nlist->curMaxNumENodes = nnodes;
			nlist->ndirnodes = nnodes;
		}

		//fp = fopen("detect_porbit_cooling.txt", "a");
		//fp = fopen("detect_porbit_swirl.txt", "a");
		//fp = fopen("detect_porbit.txt", "a");
		//fprintf(fp, "Start allocating elist...\n");
		//fclose(fp);

		if(nedges == 0)
			elist = NULL;
		else
		{
			elist = new Graph_EdgeList(nedges);
			elist->curMaxNumGedges = nedges;
			elist->nedges = 0;
		}

		//dfs_stack = NULL;
		while(!dfs_stack.empty())
			dfs_stack.pop();
	}

	~DirGraph()
	{
		//if(nlist != NULL)
		//	delete nlist;
		//if(elist != NULL)
		//	delete elist;
		nlist->finalize();
		elist->finalize();
	}

	void init_graph()
	{
		int i;
		for(i = 0; i < nlist->curMaxNumENodes; i++)
		{
			if(nlist->dirnodes[i] == NULL)
				nlist->dirnodes[i] = new DirGraph_Node();
		}

		for(i = 0; i < elist->curMaxNumGedges; i++)
		{
			if(elist->edges[i] == NULL)
				elist->edges[i] = new Graph_Edge();
		}
	}

	/*For the directed graph construction*/
	bool is_repeated_edge(int from, int to);
	bool is_repeated_edge_2(int from, int to);
	void add_to_edgelist(int node1, int node2, int &cur_index);
	void add_edge_to_node(int node, int edgeindex);

	//Graph related algorithms
	//Dijkstra
	//BFS
	//DFS
	//MST

	/*Find the strongly connected components*/
	void find_SCCS(); /*finding all SCCs using double DFS visiting*/
	void DFS_visit_nonrecursive(int u, int scc_index, int inverse); /*non-recursive DFS*/
	void DFS_visit(int u, int &time, int inverse); /*recursive DFS*/
	void add_to_dfsstack(int elem);

	void DFS(int inverse); /*DFS traverse all the nodes of given directed graph*/
	void reverse_edges();
	//void build_SCCElemList();

};




/** The SCC component class **/
class SCComponent{
public:
	//int *nodes;          //the list of pixels belonging to this component
	//int cur_nnodes;
	vector<pair<int, int>> pixels;
	int node_index;      //the index of Morse sets
	bool valid;          //true - valid SCC; false - invalid 
	
	//int nseppts;
	//int nattpts;
	//int nboundaries;
	//int *singular_tri;    //the list of triangles containing singularity 
	
	vector<pair<int, int>> singular_pixels;  // the list of pixels containing singularities
	//int nfixedpoints;
	//int *periodicorbits;
	//int nperiodicorbits;

	//conley index
	//added 12/03/2009 by Qingqing

	int XM;//Euler Characteristics of Mesh
	int XL;//Euler Characteristics of Boundary
	int Conley0;//Conley index 0
	int Conley1;//Conley index 1
	int Conley2;//Conley index 2
	int classification;//0=trival;1=source;2=sink;3=saddle
	double priority;//for auto refine
	double variance_vector;//vector variance

	//corresponding
	int global_SCC;//the SCC# of global scclist
	int local_SCC;//the SCC# of local scclist

	//function for finding a SCC component

	/*calculate all the separation and attachment points of a strongly connected component*/
	//void cal_sep_attp_pts();

	//void reset_edge_intersections();  /*reset the intersection information for all the edges
									  //of the strongly connected component*/
	bool is_connected();

};  // end of SCComponent class



#endif /* __DGGRAPH_H__ */