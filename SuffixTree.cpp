#include "SuffixTree.hpp"

/*
Shouldn't rely on a ctor this way, but this ctor guarantees a node is constructed
with all NULL ptrs and negative numeric vals.
*/
TreeNode::TreeNode()
{
    SuffixLink = NULL;
    Parent = NULL;
    StringDepth = -1;
    NodeID = -1;
}

TreeNode::~TreeNode()
{
}

int TreeNode::NumEdges()
{
    return Edges.size();
}

/*
Given a symbol, return its corresponding out-edge. Returns null if there is none.
Parameter *input* is required since edge labels correspond with indices in input.
*/
Edge* TreeNode::GetEdge(char c, const string& input)
{
    for (int i = 0; i < Edges.size(); i++) {
        if (input[Edges[i].i] == c) {
            return &Edges[i];
        }
    }

    return NULL;
}

bool TreeNode::IsInternal()
{
    return this->NodeID < 0;
}

/*
This is an unusual utility, but a child node needs a way to retrieve the edge with which
it is associated in its parent's edges. This searches the edge's of this node's parents
for the edge containing a pointer to this node.

An alternative is to just have each node also store a pointer to the edge that points to it.
Preconditon: Calling this on root would clearly be a logic error, by mccreight's algorithm this should
never occur.

TODO: This is a slow, structural requirement of not having each node store a pointer to its own associated edge. A child
must look this up in its parent, in time linear over the number of siblings it has.
*/
Edge* TreeNode::GetAssociatedEdge()
{
    for (int i = 0; i < this->Parent->NumEdges(); i++){
        if (this == this->Parent->Edges[i].Node) {
            return &this->Parent->Edges[i];
        }
    }

    #if defined DBG
        cout << "ERROR no edge found for child in GetChildsEdge! Fatal!" << endl;
    #endif

    return NULL;
}

bool TreeNode::HasChild()
{
    return this->NumEdges() == 0;
}

//Insert an edge into the node, at index in Edges given by c. Returns false if the edge already exists.
void TreeNode::AddEdge(int edge_i, int edge_j, TreeNode* edgeChild, const string& input)
{
    int i;

    #if defined DBG
    if (HasChild(input[edge_i], input)) {
        //cannot happen by structure of a suffix tree; if it does, we're in trouble
        cout << "ERROR edge already exists! FATAL. Attempted to insert edge with existing out-edge for " << input[edge_i] << endl;
    }
    #endif

    //this is terribly inefficient and abuses vector api to emulate list api, but works for now. On each insert, resize the vector and shift remaining items to maintain lexicographic order of edges
    //find index of insertion point
    i = 0;
    while(i < Edges.size() && input[Edges[i].i] < input[edge_i]){
        i++;
    }
    //i now points to position at which to insert; which may be one past end

    //reached end, so just append the new item
    if (i == Edges.size()) {
        Edges.resize(Edges.size() + 1); // resize is supposed to guarantee vector values do not change
        Edges[Edges.size() - 1].i = edge_i;
        Edges[Edges.size() - 1].j = edge_j;
        Edges[Edges.size() - 1].Node = edgeChild;
    }
    // else not at end, so insert and copy previous items right (awful)
    else {
        Edges.resize(Edges.size() + 1); // resize is supposed to guarantee vector values do not change
        //copy-shift the old elements one to the right
        for (int j = Edges.size() - 1; j > i; j--) {
            Edges[j].i = Edges[j - 1].i;
            Edges[j].j = Edges[j - 1].j;
            Edges[j].Node = Edges[j - 1].Node;
        }
        //insert the new edge at position i
        Edges[i].i = edge_i;
        Edges[i].j = edge_j;
        Edges[i].Node = edgeChild;
    }

    //lastly, set the child to point at the parent
    edgeChild->Parent = this;
}

/*
Given some symbol (char) in the overall alphabet (Sigma), this checks if the corresponding
child edge/node exists in Edges.

const string& input is required to correlate an edge-label (an index in input) with the given char c.
*/
bool TreeNode::HasChild(char c, const string& input)
{
    return GetEdge(c, input) != NULL;
}

SuffixTree::SuffixTree()
{
    _root = NULL;
    _input = NULL;
    _numEdges = 0;
    _numInternalNodes = 0;
    _numLeaves = 0;
}

/*
Input is left non-const since we may need to modify it (appending $, etc).
*/
SuffixTree::SuffixTree(string& input, const string& alphabet)
{
    _root = NULL;
    _input = NULL;
    _numEdges = 0;
    _numInternalNodes = 0;
    _numLeaves = 0;

    Build(&input,alphabet);
}

SuffixTree::~SuffixTree()
{
    Clear();
}

void SuffixTree::SetAlphabet(const string& alphabet)
{
    _alphabet = alphabet;
    if (_alphabet.find('$') == string::npos) {
        _alphabet += "$";
    }
}

bool SuffixTree::IsEmpty()
{
    return _root == NULL;
}

void SuffixTree::PrintSize()
{
    if (!IsEmpty()) {
        cout << "From input string of len " << _input->length() << ", suffix-tree has " << _numLeaves << " leaves, " << _numInternalNodes << " internal nodes, and " << _numEdges << " edges." << endl;
        cout << "SuffixTree=" << sizeof(SuffixTree) << " bytes  TreeNode=" << sizeof(TreeNode) << " bytes  Edge=" << sizeof(Edge) << " bytes" << endl;
        int totalBytes = sizeof(SuffixTree) + (_numLeaves + _numInternalNodes) * sizeof(TreeNode) + _numEdges * sizeof(Edge);
        cout << "Given the static class size overhead, the current tree consumes a minimum of " << totalBytes << " bytes." << endl;
        cout << "Actual size may differ as vector allocator allocates larger data capacity than usage per node." << endl;
        cout << "Space constant (tree size / input size): " << (totalBytes / _input->length()) << endl;
    }
    else {
        cout << "Tree empty" << endl;
    }
}

void SuffixTree::Size(int* edges, int* internalNodes, int* leaves)
{
    *edges = _numEdges;
    *internalNodes = _numInternalNodes;
    *leaves = _numLeaves;
}

void SuffixTree::Clear()
{
    if (_root != NULL) {
        cout << "Clearing tree... " << endl;
        _clear(_root);
        _root = NULL;
        _input = NULL;
        _numEdges = 0;
        _numInternalNodes = 0;
        _numLeaves = 0;
        _alphabet.clear();
    }
}

//Recursively deletes subtree under some root, and then root itself
void SuffixTree::_clear(TreeNode* root)
{
    if (root != NULL) {
        for (int i = 0; i < root->NumEdges(); i++) {
            #if defined DBG
            if (root->Edges[i].Node == NULL) {
                cout << "WARN _clear() encountered non-null Edge with null Node ptr??" << endl;
            }
            #endif
            //recursively delete the subtree beneath this child
            _clear(root->Edges[i].Node);
        }
        //delete the edges, after recursively deleting all subtrees
        if (!root->Edges.empty()) {
            root->Edges.clear();
        }
        delete root;
    }
}

void SuffixTree::PrintBfs()
{
    int i, maxDepth, averageDepth;
    TreeNode* node;
    list<TreeNode*> nodeQ;

    if (_root == NULL) {
        cout << "Tree empty." << endl;
        return;
    }

    PrintSize();

    nodeQ.push_back(_root);
    maxDepth = 0;
    averageDepth = 0;
    while (!nodeQ.empty()) {
        //deque the next node
        node = nodeQ.front();
        nodeQ.pop_front();
        //track depth statistics
        if (node->StringDepth > maxDepth)
            maxDepth = node->StringDepth;
        averageDepth += node->StringDepth;

        //print node info
        if (node->NodeID < 0)
            cout << "Internal node child edges: ";
        else
            cout << "Leaf " << node->NodeID << ": (id=" << _input->substr(node->NodeID,string::npos) << ", depth=" << node->StringDepth << ")";

        //print edge info (internal nodes only)
        for (i = 0; i < node->NumEdges(); i++) {
            string edgeLabel = _input->substr(node->Edges[i].i, node->Edges[i].j - node->Edges[i].i + 1);
            cout << "[" << node->Edges[i].i << "," << node->Edges[i].j << "](" << edgeLabel << ")   ";
            nodeQ.push_back(node->Edges[i].Node);
        }
        cout << endl;
    }

    cout << "String depth of deepest internal node: " << maxDepth << endl;
    cout << "Average string depth of internal nodes: " << (averageDepth / this->_numInternalNodes) << endl;

    cout << "\n" << endl;
}

void SuffixTree::PrintDfs() 
{
    PrintSize();
    _printDfs(_root);
    cout << "\n" << endl;
}

void SuffixTree::_printDfs(TreeNode* node)
{
    if (node != NULL) {
        for (int i = 0; i < node->NumEdges(); i++) {
            cout << node->NodeID << ":" << node->StringDepth << endl;
            _printDfs(node->Edges[i].Node);
        }
    }
}

/*
Just an assignment requirement. Uses DFS to print the BWT of the suffix tree.
The BWT is defined as traversing dfs, visiting children in lexicographic order,
and printing out the letter corresponding to the NodeID minus one. In this
way, the letter preceding any prefix given by some NodeID is is the bwt for that suffix.

BWT[i] = s[ NodeID - 1 ]  <--The previous letter! Not, say, a zero-based index adjustment.
*/
void SuffixTree::PrintBWT()
{
    cout << "SuffixTree bwt index: " << endl;
    _printBWT(_root,cout);
    cout << endl;
}

//Just an assignment requirement. Given some node, print its immediate children.
void SuffixTree::PrintChildren(TreeNode* node)
{
    for (int i = 0; i < node->NumEdges(); i++) {
            TreeNode* child = node->Edges[i].Node;
            cout << "child char: " << _input[node->Edges[i].i] << endl;
            cout << "child->ID: " << child->NodeID << endl;
            cout << "child->StringDepth: " << child->StringDepth << "\r\n" << endl;
    }
}

/*
This is implemented to correspond with the pdf notes in this repo. 

Entrant cases: 2A and 2B (from notes)
    

Cases:
    1) |Beta| chars exhausted, and we land on an internal node: return internal node as v, and call findPath on it
    2) |Beta| chars exhausted, and we land in the middle of an edge; so create an internal node v, and place leaf there (no findPath call)

Returns v, the lowest node under v_prime before |Beta| is exhausted (which may be v_prime itself).

In case 2B (see notes), u->Parent is root. If these nodes are connected by a single character edge (eg,
their difference in string depth is one), then this function should immediately return _root as v_prime,
and caller should then use v_prime (root) as its suffix link for u.

Precondition: This function assume's u's suffix link is null, and uses its parent's instead.

Params:
@u: The node u (in pdf notes), the parent internal node of the last leaf inserted
*/
TreeNode* SuffixTree::_nodeHops(TreeNode* u, const int suffixIndex)
{
    char c;
    int beta, betaIt, edgeLen;
    Edge* hoppedEdge, *uPrimeEdge;
    TreeNode* v_prime = NULL;

    //set up the hopping parameters
    uPrimeEdge = u->GetAssociatedEdge();
    betaIt = uPrimeEdge->i;
    // the difference in string depth of u and its parent, u_prime
    beta = uPrimeEdge->j - uPrimeEdge->i + 1; //plus one, since single char edges are stored as [i,i]
    v_prime = u->Parent->SuffixLink;
    //case 2B: u->parent is root, so use betaPrime instead of beta by subtracting one
    if (u->Parent == _root) {
        betaIt = suffixIndex;
        beta--;
        //If beta is empty, just return u (the root). This occurs for case 2B if |beta| = 1, and
        //|beta'| = 0, meaning the function has vacuously done its work of skipping nodes.
        if (beta <= 0) {
            return _root;
        }
    }

    //get the first outgoing edge by char
    c = _input->at(betaIt);
    hoppedEdge = v_prime->GetEdge(c,*_input);
    #if defined DBG
        if (hoppedEdge == NULL) {
            cout << "ERROR hoppedEdge NULL! Fatal." << endl;
        }
    #endif

    //get the len of this edge, to hop
    edgeLen = hoppedEdge->j - hoppedEdge->i + 1;
    //find lowest internal node below v_prime at or *before* the point |Beta| is exhausted (which may be v_prime itself)
    //while next hop would not make beta negative
    while ((beta - edgeLen) >= 0) {
        //decrement beta
        beta -= edgeLen;
        betaIt += edgeLen;
        //hop the edge
        v_prime = hoppedEdge->Node;
        //get the next edge
        c = _input->at(betaIt);
        hoppedEdge = v_prime->GetEdge(c,*_input);
        if (hoppedEdge == NULL) {
            //this is valid;
            //cout << "ERROR hopped edge NULL in nodeHops" << endl;
        }
        else {
            //get the len of this edge, to hop
            edgeLen = hoppedEdge->j - hoppedEdge->i + 1;
        }
    }
    /*Post-loop: (1) landed directly on an internal node (Beta was exausted), or (2) we stopped above an edge too
    long for the remainder of Beta. In case (1), we return v_prime. In case (2), an internal node needs to be
    created on this edge, and findPath need not be called.
    */

    if (hoppedEdge == NULL) {
        return v_prime;
    }

    #if defined DBG
    if (beta < 0) {
        cout << "ERROR beta less than zero in nodeHops(). Fatal." << endl;
    }
    #endif

    //if betaLen is > 0, a mismatch occurs on this edge. So split the edge there, and create a new internal node (v)
    if (beta > 0) {
        //beta chars remain, but edge is longer than beta; so split the edge at edge->i+beta
        v_prime = _splitEdge(v_prime, hoppedEdge, hoppedEdge->i + beta);
    }
    //elsif beta == 0, we landed on existing internal node 'v' directly, so just return it.

    return v_prime;
}

/*
  Given lastNode, go up to its parent
    1) If its suffix link is not null, follow it to starting subtree in which to insert
    2) If its suffix link is null, traverse to its parent and follow its suffix link. It is guaranteed to have a suffix link.

    On exit, this function does NOT set the suffix link of u (parent of lastLeaf) if it was null on ingress; it
    is expected that this will be done outside this function by the caller.

    Returns: Pointer to last internal node under which this suffix's leaf was inserted.
*/
TreeNode* SuffixTree::_insertSuffix(TreeNode* u, const int suffixIndex)
{
    TreeNode* v = NULL;

    #ifdef DBG
    if (u == NULL || u->Parent == NULL) {
        if (u == NULL)
            cout << "ERROR u == NULL in insertSuffix! Fatal" << endl;
        else if (u->Parent == NULL)
            cout << "ERROR u->parent == NULL in insertSuffix! Fatal" << endl;
        return NULL;
    }
    #endif

    //implements the four cases in the notes; these may be reducible, but this is clearer
    //case 1A: SL(u) known and u != root, follow link and findPath from it
    if (u->SuffixLink != NULL && u != _root) {
         //run findPath from v, with string index starting at suffixIndex plus the skipped chars of starting from v
        //cout << "dbg" << endl;
        v = _findPath(u->SuffixLink, suffixIndex + u->StringDepth - 1, suffixIndex);
    }
    //case 1B: SL(u) is known and u == root, just call findPath from root
    else if (u->SuffixLink != NULL && u == _root) {
        //cout << "dbg" << endl;
        v = _findPath(_root, suffixIndex, suffixIndex);
    }
    //case 2A: SL(u) unknown and u' != root, get v' from u', and nodeHop, then findPath
    else if (u->SuffixLink == NULL && u->Parent != _root) {
        //cout << "dbg" << endl;
        //see notes; u contains all info to get to u_prime then v_prime; nodeHops returns lowest internal node whose prefix comprises Beta
        v = _nodeHops(u, suffixIndex);
        //set link of u; this must be done after nodeHops
        if (u != _root) {
            u->SuffixLink = v;
        }
        //if v is new, expect findPath to immediately insert suffix as a leaf and return the same v
        v = _findPath(v, v->StringDepth+suffixIndex, suffixIndex);
    }
    //case 2B: SL(u) unknown and u' == root
    else if(u->SuffixLink == NULL && u->Parent == _root) {
        //cout << "dbg" << endl;
        v = _nodeHops(u, suffixIndex);
        if (u != _root) {
            u->SuffixLink = v;
        }
        v = _findPath(v, suffixIndex + v->StringDepth, suffixIndex);
    }

    #if defined DBG
        if (v == NULL) {
            cout << "ERROR inserted==NULL in insert()" << endl;
        }
    #endif

    return v;
}

/*
Utility code extracted from findPath; this function just makes the findPath() function more readable.

Given a parent node, the suffix' start index, and the edge-start index

parent: The parent under which to insert a leaf and its edge.
suffixIndex: The ABSOLUTE suffixIndex (not relative) of the new leaf.
edge_i: The start position at which the prefix of this suffix suffered a character mismatch.

*/
void SuffixTree::_insertLeaf(TreeNode* parent, const int suffixIndex, const int edge_i)
{
    //make the new leaf
    TreeNode* newLeaf = new TreeNode();
    newLeaf->NodeID = suffixIndex;
    newLeaf->Parent = parent;
    newLeaf->StringDepth = _input->length() - suffixIndex;

    //set the outlink of the parent node
    parent->AddEdge(edge_i, _input->length() - 1, newLeaf, *_input);
    _numLeaves++;
    _numEdges++;
}

/*
As with _insertLeaf, this is just code extracted from findPath so findPath is more readable.

Given the require parameters, split an edge into two edges, thread in a new internal node, and
return the pointer to new internal node. edgeSplitIndex is here defined to be the start
index of the NEW continuation edge; equivalently, it is where the letter mismatch occurred. So
if an edge with label [1-9] matches up to index 4, then edgeSplitIndex is 4, and the new edges
should be labelled [1-5] and [4-9]. Also note that the edgeSplitIndex corresponds with the edge label,
and is not an absolute value; its a relative value in range (i,j].
*/
TreeNode* SuffixTree::_splitEdge(TreeNode* parent, Edge* oldEdge, const int edgeSplitIndex)
{
    //alloc the new internal node (the new v)
    TreeNode* newInternalNode = new TreeNode();
    newInternalNode->Parent = parent;
    _numInternalNodes++;
    //set internal node id to the negation of the number of internal nodes
    newInternalNode->NodeID = -1 * _numInternalNodes;
    newInternalNode->StringDepth = parent->StringDepth + (edgeSplitIndex - oldEdge->i);

    newInternalNode->AddEdge(edgeSplitIndex, oldEdge->j, oldEdge->Node, *_input);
    //fix the original edge
    oldEdge->j = edgeSplitIndex - 1;
    oldEdge->Node = newInternalNode;
    _numEdges++;

    return newInternalNode;
}

/*
    Given a suffix index in _input, and a starting node, walk down the edges of the subtree below startNode
    comparing characters in _input[i:] to those on the edge label until a mismatch is found, and insert a new node there.
    There are two cases:
        1) New node occurs along an edge (if mismatch, or if chars are exhausted)
        2) New node occurs at an existing node

    This function assumes that the input string is terminated with '$'; as such, a mismatch will ALWAYS occur
    before the input string/suffix length is exhausted, since every preceding suffix will have a deeper string
    depth than the current one.

    Returns internal node (parent) of new leaf 'v_prime', which may be existing or newly created.

    Params:
        v_prime: some internal node under which to insert the suffix given by suffixIndex
        suffixOffset: The string index in _input of the suffix to be inserted; this value is NOT necessarily
            the original suffixIndex, since the caller may jump nodes, adding some amount to the suffixIndex accordingly.
            So suffixOffset = suffixIndex + (some number of skipped chars, per suffix links).
        suffixIndex: The ABSOLUTE start index of the suffix currently being inserted; this could be derived, but its
            clearer to pass the original value as const than to pollute the code with math expressions to derive it.
*/
TreeNode* SuffixTree::_findPath(TreeNode* v_prime, const int suffixOffset, const int suffixIndex)
{
    TreeNode* v = NULL;
    TreeNode* parent = v_prime;
    int i, edgeIt;
    Edge* edge = NULL;
    bool found;

    //walk edges, and insert new leaf when insertion point is found
    found = false;
    i = suffixOffset;
    edge = parent->GetEdge(_input->at(i),*_input);
    while (!found) {
        //if current edge (of parent) is null, create a new edge/leaf node and insert into the current parent node
        if (edge == NULL) {
            //make and insert the new edge and leaf
            found = true;
            _insertLeaf(parent, suffixIndex, i);
            v = parent;
        }
        //else walk the current edge comparing chars
        else {
            //walk the edge looking for a mismatch
            edgeIt = edge->i;
            while (i < _input->length() && edgeIt <= edge->j && _input->at(i) == _input->at(edgeIt)) {
                edgeIt++;  i++;
            }

            //see function header; this is should be unreachable; a mismatch should always occur before input/suffix exhaustion
            if (i >= _input->length()) {
                cout << "ERROR input exhausted in findPath(). Fatal." << endl;
            }

            //if edgeIt did not reach beyond end of edge there was a mismatch, so break edge and create leaf node below break
            if (edgeIt <= edge->j) {
                found = true;
                //break the existing edge and insert a new internal node there, returning it as v
                v = _splitEdge(parent, edge, edgeIt);
                //make a new leaf for the new internal node
                _insertLeaf(v, suffixIndex, i);
            }
            //else, all chars matched on this edge, so advance to the next edge
            else {
                parent = edge->Node;
                edge = parent->GetEdge(_input->at(i),*_input);
            }
        }
    }

    return v;
}

//Verify the input string ends with '$', and contains '$' only once, and no characters other than in
//nucleotide alphabet 'atcg' (plus dollar).
bool SuffixTree::_isValidInput(const string* s)
{
    if (s->length() == 0) {
        cout << "ERROR input string cannot be empty" << endl;
        return false;
    }

    for (int i = 0; i < s->length(); i++) {
        if (_alphabet.find(s->at(i)) < 0) {
            cout << "ERROR input string may contain only >" << _alphabet << "< but found char >" << s->at(i) << "<" << endl;
            return false;
        }
    }

    return true;
}

/*
Constructs a suffix-tree using suffix-links, per McCreight's algorithm.
See the pdf in this repo for an explanation of the variables and private functions;
they were written to correspond with conventions in the notes.
*/
void SuffixTree::Build(string* s, const string& alphabet)
{
    int i;
    TreeNode* v;

    //clear existing tree, if any (this nulls _input, so do it first)
    if (!this->IsEmpty()) {
        Clear();
    }

    
    #ifdef PERFTEST
    //track build time stats
    clock_t c_start = clock();
    #endif

    _input = s;
    SetAlphabet(alphabet);
    if (_input == NULL) {
        cout << "ERROR input string ptr null in Build()" << endl;
        return;
    }
    if (_input->length() == 0) {
        cout << "ERROR input string empty. Build() aborted" << endl;
        return;
    }
    if (!_isValidInput(_input)) {
        cout << "ERROR suffix tree input improperly formatted. Build() aborted" << endl;
        return;
    }

    //format requirement: make sure strings ends with '$'
    if (_input->at(_input->length() - 1) != '$') {
        (*_input) += "$";
    }

    //initialize the root node, with its suffix link and parent ptr pointing to itself
    _root = new TreeNode();
    _root->Parent = _root;
    _root->SuffixLink = _root;
    _root->NodeID = -1;
    _root->StringDepth = 0;
    _numInternalNodes = 1;

    for (i = 0, v = _root; i < s->length(); i++) {
        v = _insertSuffix(v,i);

        //report progress, for very large trees
        if (i % 10000 == 9999) {
            cout << "\rInserted " << i << " of " << s->length() << " suffixes, " << (int)(((float)i / (float)s->length()) * 100) << "% complete          " << flush;
        }
    }
    cout << "\r100% complete.                                       "<< endl;

    #ifdef PERFTEST
    clock_t c_end = clock();
    cout << fixed << setprecision(2) << "CPU time used: " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms" << endl;;
    #endif
}

/*
Simply finds the lowest internal node, in terms of string depth. This path label from root to this node necssarily represents
the longest repeated substring in the tree/input. The deepest internal node could also be tracked in the
build procedure, this just the work explicitly after building the tree, for the hw api.
*/
void SuffixTree::PrintLongestRepeatSubstring()
{
    int i, maxDepth;
    TreeNode* node, *deepest;
    list<TreeNode*> nodeQ;

    if (_root == NULL) {
        cout << "Tree empty." << endl;
        return;
    }

    cout << "Searching for longest repeated substring..." << endl;
    //search BFS for the deepest internal node
    nodeQ.push_back(_root);
    deepest = _root;
    while (!nodeQ.empty()) {
        //deque the next node
        node = nodeQ.front();
        nodeQ.pop_front();
        
        //track deepest internal node
        if (node->IsInternal() && node->StringDepth > deepest->StringDepth) {
            deepest = node;
        }

        //queue up child nodes
        for (i = 0; i < node->NumEdges(); i++) {
            nodeQ.push_back(node->Edges[i].Node);
        }
    }

    //print the deepest node
    int len = deepest->StringDepth;
    //deepest internal node necessarily has exactly two leaf children;
    int startIndex = deepest->Edges[0].Node->NodeID; //start index is not unique
    cout << "Longest repeating substring, of length " << deepest->StringDepth << " starting at string index " << startIndex << ": ";
    //print the path label, the longest repeating substring
    cout << _input->substr(startIndex, len) << endl;
}

//merge with overloaded version
void SuffixTree::_printBWT(TreeNode* node, ostream& outputStream)
{
    if (node != NULL) {
        //print leaf info
        if (node->NodeID >= 0) {
#if defined DBG
            if (node->NodeID >= _input->length()) {
                cout << "ERROR NodeID > input length in _printBWT FATAL" << endl;
            }
#endif
            //print the leaf id character by its bwt index char (leadId - 1)
            int bwtIndex = node->NodeID == 0 ? (_input->length() - 1) : (node->NodeID - 1);
            //for printing, output to screen without line breaks
            if (&outputStream == &cout) {
                outputStream << _input->at(bwtIndex);
            }
            //else, we're writing to some file ostream, so append linebreaks like Ananths test files
            else {
                outputStream << _input->at(bwtIndex) << "\n";
            }
        }
        else {
            //else this is an internal node, so print the rest of nodes using dfs, in lexicographic order
            for (int i = 0; i < node->NumEdges(); i++) {
                _printBWT(node->Edges[i].Node,outputStream);
            }
            outputStream << flush;
        }
    }
    else {
        cout << "ERROR node NULL in _printBWT FATAL" << endl;
    }
}

void SuffixTree::WriteBWT(const string& ofname)
{
    if (fileExists(ofname)) {
        ofstream outputFile(ofname, ofstream::trunc);
        if (outputFile.good()) {
            _printBWT(_root, outputFile);
            outputFile.close();
        }
        else {
            cout << "Failed to open BWT output file: " << ofname;
        }
    }
    else {
        cout << "ERROR BWT outputFile not found: " << ofname << endl;
    }
}

void SuffixTree::PrintNativeSpaceStats()
{
    cout << "Native space stats by class (bytes):\r\n\tsizeof(SuffixTree): " << sizeof(SuffixTree) << "\r\n\tsizeof(TreeNode): " << sizeof(TreeNode);
    cout << "\r\n\tsizeof(Edge): " << sizeof(Edge) << endl;
}

