#ifndef _DISK_H_
#define _DISK_H_

#include "btree.h"
#include <sstream>
using namespace std;

#define PtTid 	1000
#define PtFid 	2000
#define AdjTid 	3000
#define AdjFid 	4000
#define Ref(obj) ((void*)&obj)

extern FILE *PtFile;
extern FILE *AdjFile;
extern BTree *PtTree;

enum FixedF { SIZE_A, NI_P, NJ_P, DIST_P, SIZE_P };
enum VarE { ADJNODE_A, DIST_A, SUMKWD_A, PTKEY_A, PT_P, PT_DIST, PT_KWD };

char* getFlatBlock(FILE* myfile, int BlockId);
void getVarE(VarE type, void* buf, int BaseAddr, int pos);
void getFixedF(FixedF type, void* buf, int BaseAddr);
int getFileSize(FILE* f);
void OpenDiskComm(string fileprefix, int _cachesize);
void CloseDiskComm();

#endif

