#pragma warning(disable : 4996)
#include "StdAfx.h"
#include "EGTree.h"
#include "utility.h"
#include "KeywordsGenerator.h"
#include<metis.h>
using namespace std;

// offset 
#define _FILE_OFFSET_BITS 64
// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
//#define WEIGHT_INFLATE_FACTOR 100000
#define WEIGHT_INFLATE_FACTOR 1
// egtree fanout
#define PARTITION_PART 4
// egtree leaf node capacity = tau(in paper)
#define LEAF_CAP 32
#define ATTRIBUTE_DIMENSION 6
#define SKY_PARTITION 5
#define edDis 0.5
#define covThre 100
#define splitBlock 5

/********************************DataStructure************************************/
//-----basic infromation of gtree----------------------


int noe; // number of edges
int nLeafNode;
vector<Node> Nodes;
int coversize = 0;

// init status struct
typedef struct {
	int tnid; // tree node id
	set<int> nset; // node set
}Status;

// use for metis
// idx_t = int64_t / real_t = double
idx_t nvtxs; // |vertices|
idx_t ncon; // number of weight per vertex 
idx_t* xadj; // array of adjacency of indices
idx_t* adjncy; // array of adjacency nodes
idx_t* vwgt; // array of weight of nodes
idx_t* adjwgt; // array of weight of edges in adjncy
idx_t nparts; // number of parts to partition
idx_t objval; // edge cut for partitioning solution
idx_t* part; // array of partition vector
idx_t options[METIS_NOPTIONS]; // option array

vector<TreeNode> EGTree;

// METIS setting options
void options_setting() {
	printf("EGtree-option_setting start.\n");
	METIS_SetDefaultOptions(options);
	
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // _RB
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // _VOL
	options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // _RM
	options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; // _GROW _EDGE _NODE
	options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // _GREEDY _SEP2SIDED _SEP1SIDED
	options[METIS_OPTION_NCUTS] = 1;
	options[METIS_OPTION_NITER] = 10;
	/* balance factor, used to be 500 */
	//options[METIS_OPTION_UFACTOR] = 500;
	// options[METIS_OPTION_MINCONN];
	options[METIS_OPTION_CONTIG] = 1;
	// options[METIS_OPTION_SEED];
	options[METIS_OPTION_NUMBERING] = 0;
	// options[METIS_OPTION_DBGLVL] = 0; 
}

// read the data from EdgeMap instead of file
void init_input(int nOfNode, EdgeMapType EdgeMap) {
	nLeafNode = 0;
	// process node information
	printf("EGtree-init_input start.\n");
	printf("PROCESSING NODE...");
	// note that vertex id starts from 0
	for (int i = 0; i<nOfNode; i++) {
		Node node;
		node.adjnodes.clear();
		node.adjweight.clear();
		node.egtreepath.clear();
		node.isborder = false;
		Nodes.push_back(node);
	}
	printf("COMPLETE. NODE_COUNT=%d\n", (int)Nodes.size());

	// load edge
	printf("PROCESSING EDGE...");
	int eid;
	int snid, enid;
	float weight;
	int iweight;
	noe = 0;
	EdgeMapType::iterator iter = EdgeMap.begin();
	int n = EdgeMap.size();
	printf("EdgeMap size is :%d", n);
	for (; iter != EdgeMap.end(); iter++) {
		edge* e = iter->second;
		noe++;
		snid = e->Ni;
		enid = e->Nj;
		weight = e->dist;

		iweight = (int)(weight * WEIGHT_INFLATE_FACTOR);
		Nodes[snid].adjnodes.push_back(enid);
		Nodes[snid].adjweight.push_back(iweight);
		Nodes[enid].adjnodes.push_back(snid);
		Nodes[enid].adjweight.push_back(iweight);
	}
	printf("COMPLETE.\n");
}

// transform original data format to that suitable for METIS
void data_transform_init(set<int> &nset) {
	// nvtxs, ncon
	nvtxs = nset.size();
	ncon = 1;

	xadj = new idx_t[nset.size() + 1];
	adjncy = new idx_t[noe * 2];
	adjwgt = new idx_t[noe * 2];


	int xadj_pos = 1;
	int xadj_accum = 0;
	int adjncy_pos = 0;

	// xadj, adjncy, adjwgt
	unordered_map<int, int> nodemap;
	nodemap.clear();

	xadj[0] = 0;
	int i = 0;
	for (set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++) {
		// init node map
		nodemap[*it] = i;

		int nid = *it;
		int fanout = Nodes[nid].adjnodes.size();
		for (int j = 0; j < fanout; j++) {
			int enid = Nodes[nid].adjnodes[j];
			// ensure edges within
			if (nset.find(enid) != nset.end()) {
				// xadj_accum used to record the start and end position of adjacent nodes for current visited nodes
				xadj_accum++;

				adjncy[adjncy_pos] = enid;
				adjwgt[adjncy_pos] = Nodes[nid].adjweight[j];
				adjncy_pos++;
			}
		}
		xadj[xadj_pos++] = xadj_accum;
	}

	// adjust nodes number started by 0  ###########这部分以及下面部分要修改，为何要做映射，为什么权重改为1？(可以继续使用)
	for (int i = 0; i < adjncy_pos; i++) {
		adjncy[i] = nodemap[adjncy[i]];
	}

	// adjwgt -> 1
	if (ADJWEIGHT_SET_TO_ALL_ONE) {
		for (int i = 0; i < adjncy_pos; i++) {
			adjwgt[i] = 1;
		}
	}

	// nparts
	nparts = PARTITION_PART;

	// part
	part = new idx_t[nset.size()];
}

void init(int nOfNode, const EdgeMapType EdgeMap) {
	printf("EGtree-init start.");
	init_input(nOfNode, EdgeMap);
	options_setting();
}

void finalize() {
	delete xadj;
	delete adjncy;
	delete adjwgt;
	delete part;
}

void freetreenode(TreeNode& tn) {
	//free borders
	vector<int>().swap(tn.borders);
	vector<int>().swap(tn.children);
	vector<int>().swap(tn.leafnodes);
	vector<int>().swap(tn.union_borders);
	vector<float>().swap(tn.mind);
	set<unsigned long long>().swap(tn.coverkwds);
	//vector<int>(tn.children).swap(tn.children);
}
// graph partition
// input: nset = a set of node id
// output: <node, node belong to partition id>
unordered_map<int, int> graph_partition(set<int> &nset) {
	unordered_map<int, int> result;

	// transform data to metis
	data_transform_init(nset);
	// partition, result -> part
	// k way partition
	int returnstate;
	for (int i = 0; i < 100; i++) {
		idx_t n = adjncy[i];
		n++;
	}
	returnstate = METIS_PartGraphKway(
		&nvtxs,
		&ncon,
		xadj,
		adjncy,
		NULL,
		NULL,
		adjwgt,
		&nparts,
		NULL,
		NULL,
		options,
		&objval,
		part
		);
	
	// return result state
	if (METIS_ERROR_INPUT == returnstate) {
		printf("Input Error!");
		return result;
	}
	else if (METIS_ERROR_MEMORY == returnstate) {
		printf("Memory Error!");
		return result;
	}
	else if (METIS_ERROR == returnstate) {
		printf("Other Type Error!");
		return result;
	}
	else {

	}

	// push to result
	result.clear();
	int i = 0;
	// 又将结果变回来了，在上面采用nodemap可能是为了满足一定的条件，part中的结果是什么呢？子部分用什么来表示呢？
	for (set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++) {
		result[*it] = part[i];
		//printf("The part is:%d\n", part[i]);
	}

	// finalize
	finalize();

	return result;
}

// egtree construction
void build(EdgeMapType EdgeMap) {
	printf("EGtree-build start.\n");
	// init root
	TreeNode root;
	root.isleaf = false;
	root.father = -1;
	root.unionkwds = 0;
	EGTree.push_back(root);

	// init stack
	stack<Status> buildstack;
	Status rootstatus;
	rootstatus.tnid = 0;
	rootstatus.nset.clear();
	// careful the index of node start from 0
	for (int i = 0; i < Nodes.size(); i++) {
		rootstatus.nset.insert(i);
	}
	buildstack.push(rootstatus);

	// start to build
	unordered_map<int, int> presult;
	set<int> childset[PARTITION_PART];

	while (buildstack.size() > 0) {
		// pop top
		Status current = buildstack.top();
		buildstack.pop();

		// update egtreepath
		for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
			Nodes[*it].egtreepath.push_back(current.tnid);
		}

		// check cardinality
		if (current.nset.size() <= LEAF_CAP) {
			// build leaf node
			nLeafNode++;
			EGTree[current.tnid].isleaf = true;
			EGTree[current.tnid].leafnodes.clear();
			for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
				EGTree[current.tnid].leafnodes.push_back(*it);				
			}
			continue;
		}

		// partition
		// printf("PARTITIONING...NID=%d...SIZE=%d...", current.tnid, (int)current.nset.size() );
		int n = 0;
		n++;
		presult = graph_partition(current.nset);
		//printf("graph_part.\n");
// test---
		

		// construct child node set
		for (int i = 0; i < PARTITION_PART; i++) {
			childset[i].clear();
		}
		int slot;
		// put the nodes into corresponding sub-partition(slot)
		for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
			slot = presult[*it];
			//printf("EGtree-slot %d\n",slot);
			childset[slot].insert(*it);
		}

		// generate child tree nodes
		int childpos;
		for (int i = 0; i < PARTITION_PART; i++) {
			TreeNode tnode;
			tnode.isleaf = false;
			tnode.father = current.tnid;
			tnode.unionkwds = 0;

			// insert to EGTree first
			EGTree.push_back(tnode);
			childpos = EGTree.size() - 1;
			EGTree[current.tnid].children.push_back(childpos);

			// calculate border nodes
			EGTree[childpos].borders.clear();
			for (set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++) {

				bool isborder = false;
				for (int j = 0; j < Nodes[*it].adjnodes.size(); j++) {
					if (childset[i].find(Nodes[*it].adjnodes[j]) == childset[i].end()) {
						isborder = true;
						break;
					}
				}
				if (isborder) {
					EGTree[childpos].borders.push_back(*it);
					// update globally
					Nodes[*it].isborder = true;
				}
			}

			// add to stack
			Status ongoingstatus;
			ongoingstatus.tnid = childpos;
			ongoingstatus.nset = childset[i];
			buildstack.push(ongoingstatus);

		}

	}
	printf("EGtree-build end.\n");

}

// dump EGTree index to file
void egtree_save(string filename) {
	// FILE_GTREE
	printf("making egtreeFile\n");
	FILE *fout = fopen(filename.c_str(), "wb");
	int *buf = new int[Nodes.size()];
	float *buff = new float[Nodes.size()];
	int loop = -1;
	for (int i = 0; i < EGTree.size(); i++) {
		// borders
		printf("The id is %d\n", i);
		loop++;
		int count_borders = EGTree[i].borders.size();
		fwrite(&count_borders, sizeof(int), 1, fout);
		copy(EGTree[i].borders.begin(), EGTree[i].borders.end(), buf);
		fwrite(buf, sizeof(int), count_borders, fout);
		// children
		int count_children = EGTree[i].children.size();
		fwrite(&count_children, sizeof(int), 1, fout);
		copy(EGTree[i].children.begin(), EGTree[i].children.end(), buf);
		fwrite(buf, sizeof(int), count_children, fout);
		// isleaf
		fwrite(&EGTree[i].isleaf, sizeof(bool), 1, fout);
		// leafnodes
		int count_leafnodes = EGTree[i].leafnodes.size();
		fwrite(&count_leafnodes, sizeof(int), 1, fout);
		copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		fwrite(buf, sizeof(int), count_leafnodes, fout);
		// father
		fwrite(&EGTree[i].father, sizeof(int), 1, fout);

		//**************************************************
		//保存EGBU和EGTD信息

		// union_border
		int count_unionBorder = EGTree[i].union_borders.size();
		fwrite(&count_unionBorder, sizeof(int), 1, fout);
		copy(EGTree[i].union_borders.begin(), EGTree[i].union_borders.end(), buf);
		fwrite(buf, sizeof(int), count_unionBorder, fout);
		// mind
		int count_mind = EGTree[i].mind.size();
		fwrite(&count_mind, sizeof(float), 1, fout);
		
		// divide the factor

		copy(EGTree[i].mind.begin(), EGTree[i].mind.end(), buff);
		fwrite(buff, sizeof(float), count_mind, fout);
		// unionkwds
		fwrite(&EGTree[i].unionkwds, sizeof(unsigned long long), 1, fout);
		
		if (!iscovertest) {
			// coverkwds;
			int count_coverkwds = EGTree[i].coverkwds.size();
			fwrite(&count_coverkwds, sizeof(int), 1, fout);
			unsigned long long tempKwd;
			set<unsigned long long> ::iterator itEdit = EGTree[i].coverkwds.begin();
			for (; itEdit != EGTree[i].coverkwds.end(); itEdit++) {
				unsigned long long temp = *itEdit;
				fwrite(&temp, sizeof(unsigned long long), 1, fout);
			}
		}
		else {
			// write the 50
			int count_coverkwds = EGTree[i].coverkwds_50.size();
			fwrite(&count_coverkwds, sizeof(int), 1, fout);
			unsigned long long tempKwd;
			set<unsigned long long> ::iterator itEdit = EGTree[i].coverkwds_50.begin();
			for (; itEdit != EGTree[i].coverkwds_50.end(); itEdit++) {
				unsigned long long templ = *itEdit;
				fwrite(&templ, sizeof(unsigned long long), 1, fout);
			}
			// write the 100
			count_coverkwds = EGTree[i].coverkwds_100.size();
			fwrite(&count_coverkwds, sizeof(int), 1, fout);
			//unsigned long long tempKwd;
			itEdit = EGTree[i].coverkwds_100.begin();
			for (; itEdit != EGTree[i].coverkwds_100.end(); itEdit++) {
				unsigned long long templ = *itEdit;
				fwrite(&templ, sizeof(unsigned long long), 1, fout);
			}
			// write the 150
			count_coverkwds = EGTree[i].coverkwds_150.size();
			fwrite(&count_coverkwds, sizeof(int), 1, fout);
			//unsigned long long tempKwd;
			itEdit = EGTree[i].coverkwds_150.begin();
			for (; itEdit != EGTree[i].coverkwds_150.end(); itEdit++) {
				unsigned long long templ = *itEdit;
				fwrite(&templ, sizeof(unsigned long long), 1, fout);
			}
			// write the 200
			count_coverkwds = EGTree[i].coverkwds_200.size();
			fwrite(&count_coverkwds, sizeof(int), 1, fout);
			//unsigned long long tempKwd;
			itEdit = EGTree[i].coverkwds_200.begin();
			for (; itEdit != EGTree[i].coverkwds_200.end(); itEdit++) {
				unsigned long long templ = *itEdit;
				fwrite(&templ, sizeof(unsigned long long), 1, fout);
			}
			// write the 250
			count_coverkwds = EGTree[i].coverkwds_250.size();
			fwrite(&count_coverkwds, sizeof(int), 1, fout);
			//unsigned long long tempKwd;
			itEdit = EGTree[i].coverkwds_250.begin();
			for (; itEdit != EGTree[i].coverkwds_250.end(); itEdit++) {
				unsigned long long templ = *itEdit;
				fwrite(&templ, sizeof(unsigned long long), 1, fout);
			}
			/*for (int size = 50; size < 300; ) {
				set<unsigned long long>& temp = EGTree[i].coverkwds_50;
				if (size == 50) temp = EGTree[i].coverkwds_50;
				if (size == 100) temp = EGTree[i].coverkwds_100;
				if (size == 150) temp = EGTree[i].coverkwds_150;
				if (size == 200) temp = EGTree[i].coverkwds_200;
				if (size == 250) temp = EGTree[i].coverkwds_250;
				
				int count_coverkwds = temp.size();
				fwrite(&count_coverkwds, sizeof(int), 1, fout);
				unsigned long long tempKwd;
				set<unsigned long long> ::iterator itEdit = temp.begin();
				for (; itEdit != temp.end(); itEdit++) {
					unsigned long long templ = *itEdit;
					fwrite(&templ, sizeof(unsigned long long), 1, fout);
				}
				size += 50;
			}*/
		}
			
	}
	fclose(fout);
	delete[] buf;
}

// load EGTree index from file
void egtree_load(const char* filename, vector<TreeNode>& EGTree) {
	// FILE_GTREE
	printf("loading egtreeFile\n");
	FILE *fin = fopen(filename, "rb");
	int *buf = new int[NodeNum];
	float *buff = new float[NodeNum];
	int count_borders, count_children, count_leafnodes;
	bool isleaf;
	int father;

	// clear EGTree
	EGTree.clear();
	while (fread(&count_borders, sizeof(int), 1, fin)) {
		TreeNode tn;
		// borders
		tn.borders.clear();
		fread(buf, sizeof(int), count_borders, fin);
		for (int i = 0; i < count_borders; i++) {
			tn.borders.push_back(buf[i]);
		}
		// children
		fread(&count_children, sizeof(int), 1, fin);
		fread(buf, sizeof(int), count_children, fin);
		for (int i = 0; i < count_children; i++) {
			tn.children.push_back(buf[i]);
		}
		// isleaf
		fread(&isleaf, sizeof(bool), 1, fin);
		tn.isleaf = isleaf;
		// leafnodes
		fread(&count_leafnodes, sizeof(int), 1, fin);
		fread(buf, sizeof(int), count_leafnodes, fin);
		for (int i = 0; i < count_leafnodes; i++) {
			tn.leafnodes.push_back(buf[i]);
		}
		// father
		fread(&father, sizeof(int), 1, fin);
		tn.father = father;

		// union_border
		int count_unionBorder;
		fread(&count_unionBorder, sizeof(int), 1, fin);
		fread(buf, sizeof(int), count_unionBorder, fin);
		for (int i = 0; i < count_unionBorder; i++) {
			tn.union_borders.push_back(buf[i]);
		}
		// mind
		int count_mind;
		fread(&count_mind, sizeof(int), 1, fin);
		fread(buff, sizeof(float), count_mind, fin);
		for (int i = 0; i < count_mind; i++) {
			tn.mind.push_back(buff[i]);
		}
		// unionkwds
		fread(&tn.unionkwds, sizeof(unsigned long long), 1, fin);
		

		if (!iscovertest) {
			// coverkwds;
			int count_coverkwds;
			fread(&count_coverkwds, sizeof(int), 1, fin);
			unsigned long long tempKwd;
			for (int i = 0; i < count_coverkwds; i++) {
				fread(&tempKwd, sizeof(unsigned long long), 1, fin);
				tn.coverkwds.insert(tempKwd);
			}
		}
		else {
			// read 50
			int count_coverkwds;
			fread(&count_coverkwds, sizeof(int), 1, fin);
			unsigned long long tempKwd;
			tn.coverkwds_50.clear();
			for (int i = 0; i < count_coverkwds; i++) {
				fread(&tempKwd, sizeof(unsigned long long), 1, fin);
				tn.coverkwds_50.insert(tempKwd);
			}
			// read 100
			fread(&count_coverkwds, sizeof(int), 1, fin);
			//unsigned long long tempKwd;
			tn.coverkwds_100.clear();
			for (int i = 0; i < count_coverkwds; i++) {
				fread(&tempKwd, sizeof(unsigned long long), 1, fin);
				tn.coverkwds_100.insert(tempKwd);
			}
			// read 150
			fread(&count_coverkwds, sizeof(int), 1, fin);
			//unsigned long long tempKwd;
			tn.coverkwds_150.clear();
			for (int i = 0; i < count_coverkwds; i++) {
				fread(&tempKwd, sizeof(unsigned long long), 1, fin);
				tn.coverkwds_150.insert(tempKwd);
			}
			// read 200
			fread(&count_coverkwds, sizeof(int), 1, fin);
			//unsigned long long tempKwd;
			tn.coverkwds_200.clear();
			for (int i = 0; i < count_coverkwds; i++) {
				fread(&tempKwd, sizeof(unsigned long long), 1, fin);
				tn.coverkwds_200.insert(tempKwd);
			}
			// read 250
			fread(&count_coverkwds, sizeof(int), 1, fin);
			//unsigned long long tempKwd;
			tn.coverkwds_250.clear();
			for (int i = 0; i < count_coverkwds; i++) {
				fread(&tempKwd, sizeof(unsigned long long), 1, fin);
				tn.coverkwds_250.insert(tempKwd);
			}

			/*for (int size = 50; size < 300; ) {
				set<unsigned long long>& temp = tn.coverkwds_50;
				if (size == 50) temp = tn.coverkwds_50;
				if (size == 100) temp = tn.coverkwds_100;
				if (size == 150) temp = tn.coverkwds_150;
				if (size == 200) temp = tn.coverkwds_200;
				if (size == 250) temp = tn.coverkwds_250;

				int count_coverkwds;
				fread(&count_coverkwds, sizeof(int), 1, fin);
				unsigned long long tempKwd;
				for (int i = 0; i < count_coverkwds; i++) {
					fread(&tempKwd, sizeof(unsigned long long), 1, fin);
					temp.insert(tempKwd);
				}
				size += 50;
			}*/
		}
		EGTree.push_back(tn);
	}
	fclose(fin);
	delete[] buf;
	delete[] buff;
}

// dijkstra search, used for single-source shortest path search WITHIN one EGTree leaf node!
// input: s = source node
//        cands = candidate node list
//        graph = search graph(this can be set to subgraph)
vector<int> dijkstra_candidate(int s, vector<int> &cands, vector<Node> &graph) {
	// init
	set<int> todo;
	todo.clear();
	todo.insert(cands.begin(), cands.end());

	unordered_map<int, int> result;
	result.clear();
	set<int> visited;
	visited.clear();
	unordered_map<int, int> q;
	q.clear();
	q[s] = 0;

	// start
	int minpos, adjnode;
	int min;
	int weight;
	while (!todo.empty() && !q.empty()) {
		min = -1;
		weight = 0;
		for (unordered_map<int, int>::iterator it = q.begin(); it != q.end(); it++) {
			if (min == -1) {
				minpos = it->first;
				min = it->second;
			}
			else {
				if (it->second < min) {
					min = it->second;
					minpos = it->first;
				}
			}
			/*if (it->second < min) {
			min = it->second;
			minpos = it->first;
			}*/
		}

		// put min to result, add to visited
		result[minpos] = min;
		visited.insert(minpos);
		q.erase(minpos);

		if (todo.find(minpos) != todo.end()) {
			todo.erase(minpos);
		}

		// expand
		for (int i = 0; i < graph[minpos].adjnodes.size(); i++) {
			adjnode = graph[minpos].adjnodes[i];
			if (visited.find(adjnode) != visited.end()) {
				continue;
			}
			weight = graph[minpos].adjweight[i];

			if (q.find(adjnode) != q.end()) {
				if ((min + weight) < q[adjnode]) {
					q[adjnode] = min + weight;
				}
			}
			else {
				q[adjnode] = min + weight;
			}

		}
	}

	// output
	vector<int> output;
	output.clear();
	for (int i = 0; i < cands.size(); i++) {
		output.push_back(result[cands[i]]);
	}
	return output;
}

bool sortBySize(const unsigned long long& left, const unsigned long long& right) {
	bitset<MAX_KEYWORDS> l(left);
	bitset<MAX_KEYWORDS> r(right);
	return l.count() > r.count();
}

void handleCovKwd(int tn, vector<unsigned long long> nodecoverkwds) {
	
	// first handle kwd
	//vector<unsigned long long> temp(nodecoverkwds);
	vector<unsigned long long> temp;
	temp.clear();
	temp.assign(nodecoverkwds.begin(), nodecoverkwds.end());
	sort(temp.begin(), temp.end(), sortBySize);

	vector<unsigned long long>::iterator iti = temp.begin();
	vector<unsigned long long>::iterator itj;
	// 首先删除被包含的关键字
	for (; iti != temp.end();) {
		unsigned long long ki = *iti;
		vector<unsigned long long>::iterator it = iti;
		it++;
		for (itj = it ; itj != temp.end();) {
			unsigned long long kj = *itj;
			if ((ki&kj)==kj) {
				itj = temp.erase(itj);
			}
			else { 
				itj++;
			}
		}
		iti++;
	}
	
	// 迭代合并关键字	
	while (temp.size() > coversize) {
		// 第一步合并重叠最多的关键字
		sort(temp.begin(), temp.end(), sortBySize);
		vector<unsigned long long>::iterator unioni;
		vector<unsigned long long>::iterator unionj;
		int max = 0;
		unsigned long long intesect;
		unsigned long long maxunion;
		for (iti = temp.begin(); iti != temp.end(); ) {
			unsigned long long ki = *iti;
			vector<unsigned long long>::iterator it = iti;
			bitset<MAX_KEYWORDS> l(ki);
			if (l.count() <= max) break;
			it++;
			for (itj = it; itj != temp.end();itj++) {
				unsigned long long kj = *itj;
				bitset<MAX_KEYWORDS> r(kj);
				if (r.count() <= max) break;

				intesect = ki & kj;

				bitset<MAX_KEYWORDS> intersc(intesect);
				if (intersc.count() > max) {
					max = intersc.count();
					unioni = iti;
					unionj = itj;
					maxunion = ki | kj;
				}				
			}
			iti++;
		}
		// 删除被合并后关键字包含的元素
		//temp.erase(unioni);
		//temp.erase(unionj);
		for (iti = temp.begin(); iti != temp.end();) {
			unsigned long long ki = *iti;
			if ((maxunion&ki) == ki) {
				iti = temp.erase(iti);
			}
			else {
				iti++;
			}
		}
		temp.push_back(maxunion);
	}
	
	//------------------------注意，这里面很可能会有问题，访问同时删除---------------------	
	// 保存kwd信息
	vector<unsigned long long>::iterator itKwd = temp.begin();
	if (!iscovertest) {
		for (; itKwd != temp.end(); itKwd++) {
			EGTree[tn].coverkwds.insert(*itKwd);
		}
	}
	else {
		if (coversize == 50) {
			if (!EGTree[tn].coverkwds_50.empty()) {
				EGTree[tn].coverkwds_50.clear();
			}
			for (; itKwd != temp.end(); itKwd++) {
				EGTree[tn].coverkwds_50.insert(*itKwd);
			}
			int tsize = EGTree[tn].coverkwds_50.size();
		} 
		else if(coversize == 100) {
			if (!EGTree[tn].coverkwds_100.empty()) {
				EGTree[tn].coverkwds_100.clear();
			}
			for (; itKwd != temp.end(); itKwd++) {
				EGTree[tn].coverkwds_100.insert(*itKwd);
			}
			int tsize = EGTree[tn].coverkwds_100.size();
		}
		else if (coversize == 150) {
			if (!EGTree[tn].coverkwds_150.empty()) {
				EGTree[tn].coverkwds_150.clear();
			}
			for (; itKwd != temp.end(); itKwd++) {
				EGTree[tn].coverkwds_150.insert(*itKwd);
			}
			int tsize = EGTree[tn].coverkwds_150.size();
		}
		else if (coversize == 200) {
			if (!EGTree[tn].coverkwds_200.empty()) {
				EGTree[tn].coverkwds_200.clear();
			}
			for (; itKwd != temp.end(); itKwd++) {
				EGTree[tn].coverkwds_200.insert(*itKwd);
			}
			int tsize = EGTree[tn].coverkwds_200.size();
		} 
		else {
			if (!EGTree[tn].coverkwds_250.empty()) {
				EGTree[tn].coverkwds_250.clear();
			}
			for (; itKwd != temp.end(); itKwd++) {
				EGTree[tn].coverkwds_250.insert(*itKwd);
			}
			int tsize = EGTree[tn].coverkwds_250.size();
		}
	}	
	temp.clear();
}


// calculate the distance matrix
void hierarchy_shortest_path_calculation() {
	// level traversal
	printf("EGtree-hier start.\n");
	vector< vector<int> > treenodelevel;

	vector<int> current; // record all the treenodes in current level 
	current.clear();
	current.push_back(0);
	treenodelevel.push_back(current);
	// put all the nodes into treenodelevel according to their levels
	vector<int> mid;
	while (current.size() != 0) {
		mid = current;
		current.clear();
		for (int i = 0; i < mid.size(); i++) {
			for (int j = 0; j < EGTree[mid[i]].children.size(); j++) {
				current.push_back(EGTree[mid[i]].children[j]);
			}
		}
		if (current.size() == 0) break;
		treenodelevel.push_back(current);
	}
	printf("EGtree-hier-bottom-up calculation start.\n");
	// bottom up calculation
	// temp graph
	vector<Node> graph;
	graph = Nodes;
	vector<int> cands;
	vector<int> result;
	unordered_map<int, unordered_map<int, float> > vertex_pairs;

	// do dijkstra
	int s, t, tn, nid, cid, weight;
	vector<int> tnodes, tweight;
	set<int> nset;
	vector<unsigned long long> nodecoverkwds; // 关键系信息等该节点访问完后统一处理

	for (int i = treenodelevel.size() - 1; i >= 0; i--) {
		for (int j = 0; j < treenodelevel[i].size(); j++) {

			tn = treenodelevel[i][j];

			nodecoverkwds.clear();
			cands.clear();
			vertex_pairs.clear();
			printf("EGtree-node %d.\n",tn);
			
			if (EGTree[tn].isleaf) {
				//sort lefenodes and union_borders
				sort(EGTree[tn].leafnodes.begin(), EGTree[tn].leafnodes.end(), less<int>());
				sort(EGTree[tn].borders.begin(), EGTree[tn].borders.end(), less<int>());
				// cands = leafnodes
				cands = EGTree[tn].leafnodes;
				EGTree[tn].union_borders = EGTree[tn].borders;

				//--------------record the mind information
				int sizel = cands.size();
				int sizeb = EGTree[tn].union_borders.size();
				int size = sizel*sizeb;
				vector<float> vd(size, 0);
				EGTree[tn].mind = vd;
				for (int k = 0; k < EGTree[tn].union_borders.size(); k++) {
					//printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, EGTree[tn].union_borders[k] );
					result = dijkstra_candidate(EGTree[tn].union_borders[k], cands, graph);
					//printf("DIJKSTRA...END\n");

					// save to map
					for (int p = 0; p < result.size(); p++) {
						//EGTree[tn].mind.push_back(result[p]);
						int pos = k*sizel + p;
						//EGTree[tn].mind.push_back(result[p]);
						EGTree[tn].mind[pos] = result[p];
						vertex_pairs[EGTree[tn].union_borders[k]][cands[p]] = result[p];
					}
				}

				//---------------extend for EGTD Algorithm-----------------
				__int64 keyID;
				for (int k = 0; k < EGTree[tn].leafnodes.size(); k++) {
					// for each leafnode we handle each adjacent edge of it
					int nodein = EGTree[tn].leafnodes[k];
					__int64 nodei = nodein;
					int adjSize = Nodes[nodei].adjnodes.size();
					for (int p = 0; p < adjSize; p++) {
						int nodejn = Nodes[nodei].adjnodes[p];
						__int64 nodej = nodejn;
						keyID = getKey(nodei, nodej);

						edge *e = EdgeMap[keyID];
						int n = e->pts.size();
						int m = n;
						vector<InerNode>::iterator inode = e->pts.begin();
						// handle unionkwds
						EGTree[tn].unionkwds = EGTree[tn].unionkwds | e->sumkwds;
// ***********modified this *****handle coverkwds
						for (; inode != e->pts.end(); inode++) {
							InerNode tempN = *inode;
							if (find(nodecoverkwds.begin(), nodecoverkwds.end(), tempN.vct) == nodecoverkwds.end()) {
								nodecoverkwds.push_back(tempN.vct);
							}
						}
											
					}					
				}
				// 统一处理所有的InterNode信息
				if (!iscovertest) {
					coversize = covThre;
					handleCovKwd(tn, nodecoverkwds);
				}
				else {
					for (int size = 50; size < 300; ) {
						coversize = size;
						//printf("size is %d\n", nodecoverkwds.size());
						int sizes = nodecoverkwds.size();
						if (sizes == 0) {
							printf("*****************************size is %d\n", nodecoverkwds.size());
						}
						handleCovKwd(tn, nodecoverkwds);
						size += 50;
					}
				}
			}
			else {
				nset.clear();
				for (int k = 0; k < EGTree[tn].children.size(); k++) {
					cid = EGTree[tn].children[k];
					nset.insert(EGTree[cid].borders.begin(), EGTree[cid].borders.end());
					//------------------------处理中间节点的信息,将多个孩子节点的值赋给它
					if (k == 0) { //直接初始化
						// handle union keywords
						EGTree[tn].unionkwds = EGTree[cid].unionkwds;
						set<unsigned long long> ::iterator unionk = EGTree[cid].coverkwds.begin();
						for (; unionk != EGTree[cid].coverkwds.end(); unionk++) {
							unsigned long long kwd = *unionk;
							if (find(nodecoverkwds.begin(), nodecoverkwds.end(), kwd) == nodecoverkwds.end()) {
								nodecoverkwds.push_back(kwd);
							}
						}
						// handle cover keywords 
						if (!iscovertest) {
							set<unsigned long long> ::iterator unionk = EGTree[cid].coverkwds.begin();
							for (; unionk != EGTree[cid].coverkwds.end(); unionk++) {
								unsigned long long kwd = *unionk;
								if (find(nodecoverkwds.begin(), nodecoverkwds.end(), kwd) == nodecoverkwds.end()) {
									nodecoverkwds.push_back(kwd);
								}
							}
						}
						else {
							EGTree[tn].coverkwds_50.clear();
							EGTree[tn].coverkwds_50 = EGTree[cid].coverkwds_50;
							EGTree[tn].coverkwds_100.clear();
							EGTree[tn].coverkwds_100 = EGTree[cid].coverkwds_100;
							EGTree[tn].coverkwds_150.clear();
							EGTree[tn].coverkwds_150 = EGTree[cid].coverkwds_150;
							EGTree[tn].coverkwds_200.clear();
							EGTree[tn].coverkwds_200 = EGTree[cid].coverkwds_200;
							EGTree[tn].coverkwds_250.clear();
							EGTree[tn].coverkwds_250 = EGTree[cid].coverkwds_250;

							//for (int size = 50; size < 300; ) {
							//	set<unsigned long long> ::iterator unionk;
							//	set<unsigned long long> cset;
							//	set<unsigned long long>& tempnode = EGTree[tn].coverkwds_50;
							//	//nodecoverkwds.clear();
							//	if (size == 50) { 
							//		//unionk = EGTree[cid].coverkwds_50.begin();  
							//		cset = EGTree[cid].coverkwds_50;									
							//		tempnode = EGTree[tn].coverkwds_50;
							//	}
							//	if (size == 100) {
							//		//unionk = EGTree[cid].coverkwds_100.begin();
							//		cset = EGTree[cid].coverkwds_100;
							//		tempnode = EGTree[tn].coverkwds_100;
							//	}
							//	if (size == 150) {
							//		//unionk = EGTree[cid].coverkwds_150.begin();
							//		cset = EGTree[cid].coverkwds_150;
							//		tempnode = EGTree[tn].coverkwds_150;
							//	}
							//	if (size == 200) {
							//		//unionk = EGTree[cid].coverkwds_200.begin();
							//		cset = EGTree[cid].coverkwds_200;
							//		tempnode = EGTree[tn].coverkwds_200;
							//	}
							//	if (size == 250) {
							//		//unionk = EGTree[cid].coverkwds_250.begin();
							//		cset = EGTree[cid].coverkwds_250;
							//		tempnode = EGTree[tn].coverkwds_250;
							//	}
							//	unionk = cset.begin();
							//	for (; unionk != cset.end(); unionk++) {
							//		unsigned long long kwd = *unionk;
							//		if (find(tempnode.begin(), tempnode.end(), kwd) == tempnode.end()) {
							//			tempnode.insert(kwd);
							//		}
							//	}
							//	size += 50;
							//}
						}
					
						
						
					}
					else {		
						// handle union keywords
						EGTree[tn].unionkwds |= EGTree[cid].unionkwds;
						// handle cover keywords						
						if (!iscovertest) {
							set<unsigned long long> ::iterator unionk = EGTree[cid].coverkwds.begin();
							for (; unionk != EGTree[cid].coverkwds.end(); unionk++) {
								unsigned long long kwd = *unionk;
								if (find(nodecoverkwds.begin(), nodecoverkwds.end(), kwd) == nodecoverkwds.end()) {
									nodecoverkwds.push_back(kwd);
								}
							}
						}
						//sort(EGTree[tn].coverkwds_50.begin(), EGTree[tn].coverkwds_50.end(), sortBySize);
						//sort(EGTree[cid].coverkwds_50.begin(), EGTree[cid].coverkwds_50.end(), sortBySize);
						//sort(EGTree[tn].coverkwds_100.begin(), EGTree[tn].coverkwds_100.end(), sortBySize);
						//sort(EGTree[cid].coverkwds_100.begin(), EGTree[cid].coverkwds_100.end(), sortBySize);
						//sort(EGTree[tn].coverkwds_150.begin(), EGTree[tn].coverkwds_150.end(), sortBySize);
						//sort(EGTree[cid].coverkwds_150.begin(), EGTree[cid].coverkwds_150.end(), sortBySize);
						//sort(EGTree[tn].coverkwds_200.begin(), EGTree[tn].coverkwds_200.end(), sortBySize);
						//sort(EGTree[cid].coverkwds_200.begin(), EGTree[cid].coverkwds_200.end(), sortBySize);
						//sort(EGTree[tn].coverkwds_250.begin(), EGTree[tn].coverkwds_250.end(), sortBySize);
						//sort(EGTree[cid].coverkwds_250.begin(), EGTree[cid].coverkwds_250.end(), sortBySize);
						
						//set_union(EGTree[tn].coverkwds_50.begin(), EGTree[tn].coverkwds_50.end(), EGTree[cid].coverkwds_50.begin(), EGTree[cid].coverkwds_50.end(), EGTree[tn].coverkwds_50.begin());
						//set_union(EGTree[tn].coverkwds_100.begin(), EGTree[tn].coverkwds_100.end(), EGTree[cid].coverkwds_100.begin(), EGTree[cid].coverkwds_100.end(), EGTree[tn].coverkwds_100.begin());
						//set_union(EGTree[tn].coverkwds_150.begin(), EGTree[tn].coverkwds_150.end(), EGTree[cid].coverkwds_150.begin(), EGTree[cid].coverkwds_150.end(), EGTree[tn].coverkwds_150.begin());
						//set_union(EGTree[tn].coverkwds_200.begin(), EGTree[tn].coverkwds_200.end(), EGTree[cid].coverkwds_200.begin(), EGTree[cid].coverkwds_100.end(), EGTree[tn].coverkwds_200.begin());
						//set_union(EGTree[tn].coverkwds_250.begin(), EGTree[tn].coverkwds_250.end(), EGTree[cid].coverkwds_250.begin(), EGTree[cid].coverkwds_250.end(), EGTree[tn].coverkwds_250.begin());
						set<unsigned long long> ::iterator unionk = EGTree[cid].coverkwds_50.begin();
						for (; unionk != EGTree[cid].coverkwds_50.end(); unionk++) {
							unsigned long long kwd = *unionk;
							if (find(EGTree[tn].coverkwds_50.begin(), EGTree[tn].coverkwds_50.end(), kwd) == EGTree[tn].coverkwds_50.end()) {
								EGTree[tn].coverkwds_50.insert(kwd);
							}
						}

						unionk = EGTree[cid].coverkwds_100.begin();
						for (; unionk != EGTree[cid].coverkwds_100.end(); unionk++) {
							unsigned long long kwd = *unionk;
							if (find(EGTree[tn].coverkwds_100.begin(), EGTree[tn].coverkwds_100.end(), kwd) == EGTree[tn].coverkwds_100.end()) {
								EGTree[tn].coverkwds_100.insert(kwd);
							}
						}

						unionk = EGTree[cid].coverkwds_150.begin();
						for (; unionk != EGTree[cid].coverkwds_150.end(); unionk++) {
							unsigned long long kwd = *unionk;
							if (find(EGTree[tn].coverkwds_150.begin(), EGTree[tn].coverkwds_150.end(), kwd) == EGTree[tn].coverkwds_150.end()) {
								EGTree[tn].coverkwds_150.insert(kwd);
							}
						}

						unionk = EGTree[cid].coverkwds_200.begin();
						for (; unionk != EGTree[cid].coverkwds_200.end(); unionk++) {
							unsigned long long kwd = *unionk;
							if (find(EGTree[tn].coverkwds_200.begin(), EGTree[tn].coverkwds_200.end(), kwd) == EGTree[tn].coverkwds_200.end()) {
								EGTree[tn].coverkwds_200.insert(kwd);
							}
						}

						unionk = EGTree[cid].coverkwds_250.begin();
						for (; unionk != EGTree[cid].coverkwds_250.end(); unionk++) {
							unsigned long long kwd = *unionk;
							if (find(EGTree[tn].coverkwds_250.begin(), EGTree[tn].coverkwds_250.end(), kwd) == EGTree[tn].coverkwds_250.end()) {
								EGTree[tn].coverkwds_250.insert(kwd);
							}
						}
						
						//for (int size = 50; size < 300; ) {
						//		set<unsigned long long> ::iterator unionk;
						//		set<unsigned long long> cset;
						//		set<unsigned long long>& tempnode = EGTree[tn].coverkwds_50;
						//		//nodecoverkwds.clear();
						//		if (size == 50) {
						//			//unionk = EGTree[cid].coverkwds_50.begin();
						//			cset = EGTree[cid].coverkwds_50;
						//			tempnode = EGTree[tn].coverkwds_50;
						//		}
						//		if (size == 100) {
						//			//unionk = EGTree[cid].coverkwds_100.begin();
						//			cset = EGTree[cid].coverkwds_100;
						//			tempnode = EGTree[tn].coverkwds_100;
						//		}
						//		if (size == 150) {
						//			//unionk = EGTree[cid].coverkwds_150.begin();
						//			cset = EGTree[cid].coverkwds_150;
						//			tempnode = EGTree[tn].coverkwds_150;
						//		}
						//		if (size == 200) {
						//			//unionk = EGTree[cid].coverkwds_200.begin();
						//			cset = EGTree[cid].coverkwds_200;
						//			tempnode = EGTree[tn].coverkwds_200;
						//		}
						//		if (size == 250) {
						//			//unionk = EGTree[cid].coverkwds_250.begin();
						//			cset = EGTree[cid].coverkwds_250;
						//			tempnode = EGTree[tn].coverkwds_250;
						//		}
						//		unionk = cset.begin();
						//		for (; unionk != cset.end(); unionk++) {
						//			unsigned long long kwd = *unionk;
						//			if (find(tempnode.begin(), tempnode.end(), kwd) == tempnode.end()) {
						//				tempnode.insert(kwd);
						//			}
						//		}
						//		size += 50;
						//}
					}
				}
			
				//---------------extend for EGTD---------------
				for (int size = 50; size < 300; ) {
					coversize = size;
					nodecoverkwds.clear();
					if (size == 50) nodecoverkwds.assign(EGTree[tn].coverkwds_50.begin(), EGTree[tn].coverkwds_50.end());
					if (size == 100) nodecoverkwds.assign(EGTree[tn].coverkwds_100.begin(), EGTree[tn].coverkwds_100.end());
					if (size == 150) nodecoverkwds.assign(EGTree[tn].coverkwds_150.begin(), EGTree[tn].coverkwds_150.end());
					if (size == 200) nodecoverkwds.assign(EGTree[tn].coverkwds_200.begin(), EGTree[tn].coverkwds_200.end());
					if (size == 250) nodecoverkwds.assign(EGTree[tn].coverkwds_250.begin(), EGTree[tn].coverkwds_250.end());
					//printf("size is %d\n", nodecoverkwds.size());
					int sizes = nodecoverkwds.size();
					if (sizes == 0) {
						printf("*****************************size is %d\n", nodecoverkwds.size());
					}
					handleCovKwd(tn, nodecoverkwds);
					size += 50;
				}
				//handleCovKwd(tn, nodecoverkwds);

				cands.clear();
				for (set<int>::iterator it = nset.begin(); it != nset.end(); it++) {
					cands.push_back(*it);
				}
				EGTree[tn].union_borders = cands;
				//sort lefenodes and union_borders
				sort(EGTree[tn].union_borders.begin(), EGTree[tn].union_borders.end(), less<int>());
				sort(cands.begin(), cands.end(), less<int>());

				//------------------------record the mind information			
				// for each border, do min dis
				//int cc = 0;
				int csize = cands.size();
				int size = ((csize + 1)*csize) / 2;
				vector<float> vd(size, 0);
				EGTree[tn].mind = vd;
				for (int k = 0; k < EGTree[tn].union_borders.size(); k++) {
					//printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, EGTree[tn].union_borders[k] );
					result = dijkstra_candidate(EGTree[tn].union_borders[k], cands, graph);
					//printf("DIJKSTRA...END\n");

					// save to map
					for (int p = 0; p < result.size(); p++) {
						int pos = 0;
						if (k <= p) {
							pos = (p*(p + 1)) / 2 + k;
						}
						else {
							pos = (k*(k + 1)) / 2 + p;
						}
						/*if (k <= p) {
							EGTree[tn].mind.push_back(result[p]);
						}*/			
						float minds = EGTree[tn].mind[pos];
						float res = result[p];
						if (EGTree[tn].mind[pos]>0 && EGTree[tn].mind[pos] != result[p]) {
							if (result[p] < EGTree[tn].mind[pos]) {
								EGTree[tn].mind[pos] = result[p];
							}
 							//printf("The distance compute Error!\n");
						}
						EGTree[tn].mind[pos] = result[p];
						vertex_pairs[EGTree[tn].union_borders[k]][cands[p]] = result[p];
						//vertex_pairs[EGTree[tn].union_borders[p]][cands[k]] = result[p];
					}
				}
			}			

			// IMPORTANT! after all border finished, degenerate graph,###用于简化Dijkstra算法距离计算
			// first, remove inward edges
			for (int k = 0; k < EGTree[tn].borders.size(); k++) {
				s = EGTree[tn].borders[k];
				tnodes.clear();
				tweight.clear();
				for (int p = 0; p < graph[s].adjnodes.size(); p++) {
					nid = graph[s].adjnodes[p];
					weight = graph[s].adjweight[p];
					// if adj node in same tree node
					// ????
					if (graph[nid].egtreepath.size() <= i || graph[nid].egtreepath[i] != tn) {
						// only leave those useful
						tnodes.push_back(nid);
						tweight.push_back(weight);
					}
				}
				// cut it
				graph[s].adjnodes = tnodes;
				graph[s].adjweight = tweight;
			}
			// second, add inter connected edges
			for (int k = 0; k < EGTree[tn].borders.size(); k++) {
				for (int p = 0; p < EGTree[tn].borders.size(); p++) {
					if (k == p) continue;
					s = EGTree[tn].borders[k];
					t = EGTree[tn].borders[p];
					graph[s].adjnodes.push_back(t);
					graph[s].adjweight.push_back(vertex_pairs[s][t]);
				}
			}
		}
	}
	printf("EGtree-hier end.\n");
}

int mainEgtree(int nOfNode, EdgeMapType EdgeMap) { 
	// initiate the input and gtree parameters	
	printf("init.\n");
	init(nOfNode, EdgeMap);
	printf("build.\n");
	// egtree_build
	build(EdgeMap);
	printf("hier.\n");
	// hierarchy_build the distance matrix
	hierarchy_shortest_path_calculation();
	//finalize();
	return 0;
}