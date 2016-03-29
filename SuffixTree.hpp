#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include "Util.hpp"

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
        inline int NumEdges();
        Edge* GetEdge(char c, const string& input);
        Edge* GetAssociatedEdge();
        bool HasChild(char c, const string& input);
        //int AlphaToEdgeIndex(char c);
        //char EdgeIndexToAlpha(int index);
        //bool InsertEdge(Edge* edge, char c, const string& input);
        void AddEdge(int edge_i, int edge_j, TreeNode* edgeChild, const string& input);
        //For leaves, this represents the starting index of this suffix in input; less than zero for internal nodes
        int NodeID;
        //the character-length of all edges up to this node. 
        int StringDepth;
        //Need |Sigma| edge pointers; here five, for the alphabet {a,t,c,g,$}
        //Edge* Edges[ALPHABET_SIZE];
        //Edges MUST be stored in lexigraphic order
        vector<Edge> Edges;
        TreeNode* SuffixLink;
        TreeNode* Parent;
};

class SuffixTree {
public:
    SuffixTree();
    SuffixTree(string& input, const string& alphabet);
    ~SuffixTree();
    void Build(string* s, const string& alphabet);
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
    int _numLeaves, _numInternalNodes, _numEdges;
    string _alphabet;
    //stores translation from some character to its index in
    //vector<int> _alphaToIndexTable;
    //vector<char> _indexToAlphaTable;
    TreeNode* _root;
    string* _input;
    bool _isValidInput(const string* s);
    //char _edgeIndexToAlpha(int index);
    //int _alphaToEdgeIndex(char c);
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
