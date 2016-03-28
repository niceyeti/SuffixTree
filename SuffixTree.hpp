#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "Util.hpp"

#define DOLLAR 0
#define BASE_A 1
#define BASE_C 2
#define BASE_G 3
#define BASE_T 4
#define ALPHABET_SIZE 5

#define DBG

using namespace std;

class TreeNode;

typedef struct edge {
    int i, j;
    TreeNode* Node;
}Edge;

/*
TreeNode is a primitive data structure in this algorithm, all public.
Note this version has a fixed alphabet size and alphabet definition: 'ATCG' only.
*/
class TreeNode {
    public:
        TreeNode();
        ~TreeNode();
        bool InsertEdge(Edge* edge, char c);
        inline int NumEdges();
        Edge* GetEdge(char c);
        Edge* GetAssociatedEdge();
        bool HasChild();
        bool HasChild(char c);
        int AlphaToEdgeIndex(char c);
        char EdgeIndexToAlpha(int index);
        //leafid represents the starting index of this suffix in input; less than zero for internal nodes
        //this may need to become NodeID, for future assignment
        int LeafID;
        //the character-length of all edges up to this node. 
        int StringDepth;
        //Need |Sigma| edge pointers; here five, for the alphabet {a,t,c,g,$}
        Edge* Edges[ALPHABET_SIZE];
        TreeNode* SuffixLink;
        TreeNode* Parent;
};

class SuffixTree {
public:
    SuffixTree();
    SuffixTree(const string& alphabet);
    ~SuffixTree();
    void Build(const string* s);
    void Clear();
    void PrintBfs();
    void PrintDfs();
    bool IsEmpty();
    //string longestRepeatedSubstring();
    void PrintChildren(TreeNode* node);
    void PrintBWT();
    void Size(int* edges, int* internalNodes, int* leaves);
    void PrintSize();
    void SetAlphabet(const string& alphabet);
private:
    int numLeaves, numInternalNodes, numEdges;
    TreeNode* _root;
    string _alphabet;
    string const* _input;
    bool _isValidInput(const string* s);
    TreeNode* _nodeHops(TreeNode* u, const int suffixIndex);
    TreeNode* _nodeHopsOLD(TreeNode* u, const int suffixIndex);
    void _insertLeaf(TreeNode* parent, const int suffixIndex, const int edge_i);
    TreeNode* _splitEdge(TreeNode* parent, Edge* oldEdge, const int edgeSplitIndex);
    void _clear(TreeNode* node);
    void _printDfs(TreeNode* node);
    void _printBWT(TreeNode* node);
    TreeNode* _findPath(TreeNode* v, const int startOffset, const int suffixIndex);
    //TreeNode* _findPath1(TreeNode* v, const int startIndex);
    TreeNode* _insertSuffix(TreeNode* lastInserted, const int suffixIndex);
};
