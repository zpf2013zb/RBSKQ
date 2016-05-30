#include  "diskbased.h"
#include "utility.h"
#include "netshare.h"
#include <memory.h>

float alpha;
int PtNum;
int PtFileSize, AdjFileSize;

extern int i_capacity;
char *gBuffer;	// temp, global buffer for atomic use
FreqCache FC_A, FC_P; // FC for ensureCache exclusive use !
extern int BlkLen;
extern int  algorithmId;


// read the blockID from FC-cache-myfile
char* getFlatBlock(FILE* myfile, int BlockId)
{
	int UserId, FileSize;
	FreqCache* FC;
	int b = BlockId;
	if (myfile == PtFile)
	{
		FC = &FC_P;
		UserId = PtFid;
		FileSize = PtFileSize;
	}
	else
	{
		FC = &FC_A;
		UserId = AdjFid;
		FileSize = AdjFileSize;
	}
	if (inFreqCache(*FC, UserId, BlockId)) return (FC->buffer);

	if (!getCacheBlock(FC->buffer, UserId, BlockId))
	{
		int readsize = BlkLen<(FileSize - BlockId*BlkLen) ? BlkLen : (FileSize - BlockId*BlkLen);//Modified By Qin Xu
		fseek(myfile, BlockId*BlkLen, SEEK_SET);
		fread(FC->buffer, readsize, sizeof(char), myfile);
		storeCacheBlock(FC->buffer, UserId, BlockId, isData);
	}
	FC->UserId = UserId;
	FC->BlockId = BlockId;
	return (FC->buffer);
}

//Modified by Qin Xu
int PTGRP_HEADSIZE = 3 * sizeof(int) + sizeof(float);
int PTGRP_ITEMSIZE = sizeof(InerNode);
int ADJGRP_HEADSIZE = sizeof(int);
int ADJGRP_ITEMSIZE = 2 * sizeof(int) + sizeof(float) + sizeof(unsigned long long);	// changed

//enum FixedF { SIZE_A, NI_P, NJ_P, DIST_P, SIZE_P };
//enum VarE { ADJNODE_A, DIST_A, SUMKWD_A, PTKEY_A, PT_P, PT_DIST, PT_KWD };

// get the variant from the adjFile or PtFile, put the pos item from addr BaseAddr to buffer
void getVarE(VarE type, void* buf, int BaseAddr, int pos)
{
	// note default values !
	FILE* f = AdjFile;
	int size = sizeof(int), addr = -1;
	// -------------------M--modify the varBase part with flexible keyword number----
	int VarBaseA = BaseAddr + ADJGRP_HEADSIZE + pos*ADJGRP_ITEMSIZE;
	bool flag = 0; // use to exit the computation earlier

	// for VarE in AdjFile
	if (type == ADJNODE_A) addr = VarBaseA;
	if (type == DIST_A) {
		addr = VarBaseA + sizeof(int);
		size = sizeof(float);
	}
	if (type == SUMKWD_A) {
		addr = VarBaseA + sizeof(int) + sizeof(float);
		size = sizeof(unsigned long long);
	}
	if (type == PTKEY_A) {
		addr = VarBaseA + sizeof(int) + sizeof(float) + sizeof(unsigned long long);
	}

	// for VarE in PtFile
	int VarBaseP = BaseAddr + PTGRP_HEADSIZE + pos*PTGRP_ITEMSIZE;
	
	if (type == PT_P)
	{
		addr = VarBaseP;
		f = PtFile;
	}
	if (type == PT_DIST)
	{
		addr = VarBaseP + sizeof(int);
		size = sizeof(float);
		f = PtFile;
	}
	if (type == PT_KWD)
	{
		addr = VarBaseP + sizeof(int) + sizeof(float);
		size = sizeof(unsigned long long);
		f = PtFile;
	}

	char* BlockAddr = getFlatBlock(f, addr / BlkLen);
	memcpy(buf, BlockAddr + (addr%BlkLen), size);

}
// put the head information from addr BaseAddr to buf
void getFixedF(FixedF type, void* buf, int BaseAddr)
{
	// note default values !
	FILE* f = PtFile;
	int size = sizeof(int), addr = -1;

	// for FixedF in PtFile
	if (type == NI_P) addr = BaseAddr;
	if (type == NJ_P) addr = BaseAddr + sizeof(int);
	if (type == DIST_P)
	{
		addr = BaseAddr + 2 * sizeof(int);
		size = sizeof(float);
	}
	if (type == SIZE_P)
	{
		addr = BaseAddr + 2 * sizeof(int) + sizeof(float);
	}

	// for FixedF in AdjFile
	if (type == SIZE_A)
	{
		addr = BaseAddr;
		f = AdjFile;
	}

	char* BlockAddr = getFlatBlock(f, addr / BlkLen);
	memcpy(buf, BlockAddr + (addr%BlkLen), size);
}
// return the NodeID address in AdjList, note the NodeID from 1
/*
int getAdjListGrpAddr(int NodeID)  	// using AdjFile
{
return paID[NodeID].addr;
}
*/
int getFileSize(FILE* f)  	// side effect, setting f to begin
{
	fseek(f, 0, SEEK_END);
	int filesize = ftell(f);
	fseek(f, 0, SEEK_SET);
	return filesize;
}
// open the file and index
void OpenDiskComm(string fileprefix, int _cachesize)
{
	char tmpFileName[255];

	BlkLen = getBlockLength();
	int header_size = sizeof(char) + sizeof(int);
	gBuffer = new char[BlkLen];
	i_capacity = (BlkLen - header_size) / (sizeof(int) + sizeof(int));
	//printf("Blocklength=%d, i_cap=%d\n", BlkLen, i_capacity);

	InitFreqCache(FC_A);
	InitFreqCache(FC_P);
	InitCache(_cachesize);

	string name;
	name.clear();
	if (algorithmId > 2) {
		name = fileprefix + "\\pfile_tkde";
	}
	else {
		name = fileprefix + "\\pfile";
	}
	//name = fileprefix + "\\pfile";
	//sprintf(tmpFileName, "%s\\pfile", fileprefix);
	PtFile = fopen(name.c_str(), "rb");
	CheckFile(PtFile, name.c_str());
	PtFileSize = getFileSize(PtFile);

	name.clear();
	if (algorithmId > 2) {
		name = fileprefix + "\\pbtree_tkde";
	}
	else {
		name = fileprefix + "\\pbtree";
	}
	//name = fileprefix + "\\pbtree";
	char* c;
	int len = name.length();
	c = new char[len + 1];
	strcpy(c, name.c_str());
	//sprintf(tmpFileName, "%s\\pbtree", fileprefix);
	// load btree information from file and construct 128 cache block
	PtTree = new BTree(c, 128); // number of cached file entries: 128
	PtNum = PtTree->UserField;

	name.clear();
	if (algorithmId > 2) {
		name = fileprefix + "\\adjlist_tkde";
	}
	else {
		name = fileprefix + "\\adjlist";
	}
	//sprintf(tmpFileName, "%s\\adjlist", fileprefix);
	AdjFile = fopen(name.c_str(), "rb");
	CheckFile(AdjFile, tmpFileName);
	fread(Ref(NodeNum), 1, sizeof(int), AdjFile);	//	NodeNum=AdjTree->UserField;
	AdjFileSize = getFileSize(AdjFile);

	//printf("PtFileSize: %d, AdjFileSize: %d, NodeNum: %d, PtNum: %d\n", PtFileSize, AdjFileSize, NodeNum, PtNum);
}

void CloseDiskComm()
{
	fclose(PtFile);
	fclose(AdjFile);
	delete PtTree;
	delete[] gBuffer;

	//printPageAccess();
	DestroyCache();
	DestroyFreqCache(FC_A);
	DestroyFreqCache(FC_P);
}
