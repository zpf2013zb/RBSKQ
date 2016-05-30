#ifndef __NET_SHARED_
#define __NET_SHARED_

#include <queue>
#include <bitset>
#include <vector>
#include <set>
#include "utility.h"

using namespace std;
// handle
#define MAX_DIST 9999999999.0
#define MAX_KEYWORDS 64

extern int NodeNum;
extern int EdgeNum;
extern int NNnum;
// define the query point
struct roadParameter {
	int avgNKwd;
	int avgNPt;
};

struct QueryPoint
{
	int Ni, Nj;
	float dist_Ni;
	float distEdge;
};

struct Query {
	int k;
	unsigned long long keywords;
	vector<QueryPoint> querypts;
};

/*
ostream& operator<<(ostream& os, Query& Q)
{
	os << bitset<64>(Q.keywords).to_string() << " "<< Q.k;
	vector<QueryPoint> ::iterator iter = Q.querypts.begin();
	for (; iter != Q.querypts.end(); iter++) {
		os << "," << iter->Ni << " " << iter->Nj << " " << iter->dist_Ni;
	}
	return os;
}

struct POI
{
	int Ni, Nj;
	int pre;//Refine use only
	unsigned long long keywords;
	float dist_Ni;
	float dist_Nj;
	float dist_toquery;
};


struct POIComparison
{
	bool operator () (const POI& left, const POI& right) const
	{
		return left.dist_toquery> right.dist_toquery;
	}
};


ostream& operator<<(ostream& os, const POI& poi)
{
	os << "(" << poi.Ni << "," << poi.Nj << ") ";
	os << bitset<64>(poi.keywords).to_string();
	os << " distNi:" << poi.dist_Ni;
	os << " distQ:" << poi.dist_toquery;
	os << " preNode:" << poi.pre << endl;

	return os;
}
*/
struct DStepEvent
{
    double dist;
    int node;
};

struct DStepComparison
{
    bool operator () (const DStepEvent& left, const DStepEvent& right) const
    {
        return left.dist > right.dist;
    }
};

typedef	priority_queue<DStepEvent,vector<DStepEvent>,DStepComparison> DStepQueue;


struct point
{
    int Ni,Nj,pos;
};

//Modified by Qin Xu
//Denote POI on RoadEdge

// record the information of POIcc
struct InerNode
{
	int pid;
	float dis;
	unsigned long long vct;//Vector of keywords denoted by 64-bit
};

struct edge
{
	int FirstRow;
	int Ni, Nj;
	float dist;
	unsigned long long sumkwds;
	FastArray<InerNode> pts;
	//float dLB,dUB;		// for gendata use only
};

typedef map<__int64,edge*> EdgeMapType; // map node to edges

// build AdjList on fly
extern FastArray<int>* AdjList;
extern FastArray<point> PtList;
extern EdgeMapType EdgeMap;	// key: i0*NodeNum+j0

// get the key of edge
inline __int64 getKey(__int64 i, __int64 j)
{
	__int64 i0=i<j?i:j,j0=j<i?i:j;	// be careful of the order in other places
	//unsigned int key = (i0*NodeNum + j0);
	__int64 key = (i0*NodeNum + j0);
    return key; // map this edge to unique key
}
// break the index of key of edge
inline void breakKey(__int64 key,int& Ni,int& Nj)
{
    Ni=key/NodeNum;
    Nj=key%NodeNum;
}

//-------
// Step
//-------
struct StepEvent {
	float dist;
	int node, ClusID;	// ClusID for multiple expansion
	float accDist;		// posDist for gendata.cc only

	float gdist, hdist;
	float xc, yc;
	bool isSortAcc;
};

struct StepComparison
{
	bool operator () (const StepEvent& left, const StepEvent& right) const
	{
		return left.dist > right.dist;
	}
};

typedef	priority_queue<StepEvent, vector<StepEvent>, StepComparison> StepQueue;

//void printStepEvent(StepEvent& event);

/*
void printEdgeKeys()
{
    int Ni,Nj;
    printf("EdgeKeys:\n");
    EdgeMapType::iterator p=EdgeMap.begin();
    while (p!=EdgeMap.end())
    {
        edge* e=p->second;
        breakKey(p->first,Ni,Nj);
        printf("%d %d %d\n",Ni,Nj,(e==NULL));
        p++;
    }
}
*/
#endif // __NET_SHARED


