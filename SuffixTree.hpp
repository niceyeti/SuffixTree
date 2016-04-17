#pragma once

#ifndef _IOSTREAM_
#include <iostream>
#endif

#ifndef _STRING_
#include <string>
#endif

#ifndef _VECTOR_
#include <vector>
#endif

#ifndef _LIST_
#include <list>
#endif

#ifndef _IOMANIP_
#include <iomanip>
#endif

#ifndef _CHRONO_
#include <chrono>
#endif

#ifndef _CTIME_
#include <ctime>
#endif


#include "Util.hpp"

#define DBG

#define PERFTEST

using namespace std;

class TreeNode;

typedef struct edge {
    int i, j;
    TreeNode* Node;
}Edge;

/*
TreeNode is a primitive data structure in this algorithm, all public. This TreeNode is specialized
fr the SuffixTree, so don't remove it; treat it as a private class.
*/
class TreeNode {
    public:
        TreeNode();
        ~TreeNode();
        inline int NumEdges();
        Edge* GetEdge(char c, const string& input);
        Edge* GetAssociatedEdge();
        inline bool IsInternal();
        //has any child?
        inline bool HasChild();
        //has given child?
        bool HasChild(char c, const string& input);
        void AddEdge(int edge_i, int edge_j, TreeNode* edgeChild, const string& input);
        //For leaves, this represents the starting index of this suffix in input; less than zero for internal nodes
        int NodeID;
        //the character-length of all edges up to this node. 
        int StringDepth;

        //prg3 additions; should keep these. They are used for range querying leaves below some internal node, per prg3 (see spec for explanation)
        int StartLeafIndex; 
        int EndLeafIndex;

        //Edges MUST be stored in lexigraphic order
        vector<Edge> Edges;
        TreeNode* SuffixLink;
        TreeNode* Parent;
};

/*
The SuffixTree. This implementation is neither very space or time efficient in construction, but
is robust in terms of reliability/testing. Better performance could be gained. The primary bottleneck
is the use of vector<Edge> in the TreeNode class, in which a vector is used for list-api insertions (yuck).
Remember, the vector mem allocation strategy is likely to alloc more mem than used, up to some power of two,
so raw tree size stats are not correct wrt internal buffers.
For better performance a much easier, faster, smaller class cold be devised, it just depends what assumptions
you want to make. Rather than support flexible alphabets, personally I think a compile-time defined tree (possibly
using templating) with a fixed alphabet would be far better than a flexible alphabet treenode; the flexible alphabet
thing is a usability requirement not a dynamic one (eg, the alphabet would never be modified at runtime). Compiling
with fixed alphabets and a faster edge-lookup structure would make this class way faster, if needed. But its still
plenty fast (~15-20 second tree build time for a 1MB string, on an i5, consuming about 250 MB).

If code is ever used again:
1) better software patterns, mainly as described above
2) there are some mixed C/C++ idioms like pointers and references, need to get rid of these
*/
class SuffixTree {
public:
    SuffixTree();
    SuffixTree(string& input, const string& alphabet);
    ~SuffixTree();
    void Build(string* s, const string& alphabet);
    void Clear();
    void PrintBfs();
    void PrintDfs();
    int FindLoc(const string& read);
    bool IsEmpty();
    void PrintChildren(TreeNode* node);
    void PrintNativeSpaceStats();
    void PrintBWT();
    void Size(int* edges, int* internalNodes, int* leaves);
    void PrintSize();
    void WriteBWT(const string& outputFile);
    void SetAlphabet(const string& alphabet);
    void PrintLongestRepeatSubstring();

    //TODO: There is a direct relationship between 'A' and the BWT output of a given tree. In the future, this object
    //could be improved by refactoring the logic for each, since it was written for separate assignments.
    void PrepareST();
    vector<int> A;
private:
    int _numLeaves, _numInternalNodes, _numEdges, _nextIndex;
    string _alphabet;
    TreeNode* _root;
    TreeNode*  _deepestInternalNode;
    string* _input;
    void _prepareST_DFS(TreeNode* node, vector<int>& A);
    bool _isValidInput(const string* s);
    TreeNode* _nodeHops(TreeNode* u, const int suffixIndex);
    void _insertLeaf(TreeNode* parent, const int suffixIndex, const int edge_i);
    TreeNode* _splitEdge(TreeNode* parent, Edge* oldEdge, const int edgeSplitIndex);
    void _clear(TreeNode* node);
    void _printDfs(TreeNode* node);
    //void _printBWT(TreeNode* node);
    void _printBWT(TreeNode* node, ostream& outputStream);
    TreeNode* _findPath(TreeNode* v, const int startOffset, const int suffixIndex);
    TreeNode* _insertSuffix(TreeNode* lastInserted, const int suffixIndex);
};
