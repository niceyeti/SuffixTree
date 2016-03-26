#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <list>

#define BASE_A 0
#define BASE_T 1
#define BASE_C 2
#define BASE_G 3
#define DOLLAR 4
#define ALPHABET_SIZE 5

#define DBG

using namespace std;

class TreeNode;

typedef struct edge {
    int i, j;
    TreeNode* Node;
}Edge;

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
        //only meaningful for leaves; -1 for internal nodes
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
    void Size(int* edges, int* internalNodes, int* leaves);
    void PrintSize();
    void SetAlphabet(const string& alphabet);
private:
    int numLeaves, numInternalNodes, numEdges;
    TreeNode* _root;
    string _alphabet;
    string const* _input;
    bool _isValidInput(const string* s);
    TreeNode* _getLinkedSubtree(TreeNode* u, const int suffixIndex);
    TreeNode* _nodeHops(TreeNode* u, const int suffixIndex);
    TreeNode* _nodeHopsOLD(TreeNode* u, const int suffixIndex);
    void _insertLeaf(TreeNode* parent, const int suffixIndex, const int edge_i);
    TreeNode* _splitEdge(TreeNode* parent, Edge* oldEdge, const int edgeSplitIndex);
    void _clear(TreeNode* node);
    void _printDfs(TreeNode* node);
    TreeNode* _findPath(TreeNode* v, const int startOffset, const int suffixIndex);
    //TreeNode* _findPath1(TreeNode* v, const int startIndex);
    TreeNode* _insertSuffix(TreeNode* lastInserted, const int suffixIndex);
};
