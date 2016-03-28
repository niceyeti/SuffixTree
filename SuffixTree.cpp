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
    LeafID = -1;

    for (int i = 0; i < 5; i++) {
        Edges[i] = NULL;
    }
}

TreeNode::~TreeNode()
{
    //Verifies that all children of a node are deleted before the node; this requires deletion traversal sets a nodes child pointers to null before deleting the node
    for (int i = 0; i < NumEdges(); i++) {
        if (Edges[i] != NULL) {
            cout << "ERROR TreeNode child pointer not null in ~TreeNode (bottom up Node deletion assumed)" << endl;
        }
    }
}

/*
CRITICAL: This is fixed, and does not count the number of non-null child edges.
IF nodes are made to support flexible alphabets, this needs to be changed.
*/
int TreeNode::NumEdges()
{
    return sizeof(Edges) / sizeof(Edge*);
}

/*
Given a symbol, return its corresponding out-edge. Returns null if there is none.
*/
Edge* TreeNode::GetEdge(char c)
{
    return Edges[AlphaToEdgeIndex(c)];
}

/*
This is an unusual utility, but a (child) node needs a way to retrieve the edge with which
it is associated in its parent's edges. This searches the edge's of this node's parents
for the edge containing a pointer to this node.
*/
Edge* TreeNode::GetAssociatedEdge()
{
    Edge* childsEdge = NULL;

    for (int i = 0; i < this->Parent->NumEdges(); i++){
        if (this->Parent->Edges[i] != NULL && this == this->Parent->Edges[i]->Node) {
            childsEdge = this->Parent->Edges[i];
            break;
        }
    }

    #if defined DBG
    if (childsEdge == NULL) {
        cout << "ERROR no edge found for child in GetChildsEdge! Fatal!" << endl;
    }
    #endif

    return childsEdge;
}

//Insert an edge into the node, at index in Edges given by c. Returns false if the edge already exists.
bool TreeNode::InsertEdge(Edge* edge, char c)
{
    bool success = false;

    if (!HasChild(c)) {
        Edges[AlphaToEdgeIndex(c)] = edge;
        success = true;
    }
    else {
        cout << "ERROR edge could not be inserted in InsertEdge()! Child already exists. Fatal" << endl;
    }
//    bones seboomboom! skeletor's    groom boon room womb tomb doom loom zoom-zoom vroom boom-boom moon rheum assume 
    // bones seboomboom! can fit through the 

    return success;
}

/*
Given some character in the input string, map it to its index in each node's Edges array.
This is fine as-is for small alphabets (here four), but a lookup table would be better/faster
for larger alphabets.
*/
char TreeNode::EdgeIndexToAlpha(int index)
{
     char c = '\0';

    switch (index) {
        case BASE_A:
            c = 'A';
        break;

        case BASE_T:
            c = 'T';
        break;

        case BASE_C:
            c = 'C';
        break;

        case BASE_G:
            c = 'G';
        break;

        case DOLLAR:
            c = '$';
        break;

        default:
            cout << "ERROR unmapped char in alphaToEdgeIndex: >" << c << "< FATAL" << endl;
        break;
    }

    return c;
}

/*
Given some character in the input string, map it to its index in each node's Edges array.
This is fine as-is for small alphabets (here four), but a lookup table would be better/faster
for larger alphabets.

*/
int TreeNode::AlphaToEdgeIndex(char c)
{
    int index = -1;

    switch ((int)c) {
    case 'a':
    case 'A':
        index = BASE_A;
        break;

    case 'c':
    case 'C':
        index = BASE_C;
        break;

    case 't':
    case 'T':
        index = BASE_T;
        break;

    case 'g':
    case 'G':
        index = BASE_G;
        break;

    case '$':
        index = DOLLAR;
        break;

    default:
        cout << "ERROR unmapped char in alphaToEdgeIndex: >" << c << "< FATAL" << endl;
        break;
    }

    return index;
}

/*
Given some symbol (char) in the overall alphabet (Sigma), this checks if the corresponding
child edge/node exists in Edges.
*/
bool TreeNode::HasChild(char c)
{
    return Edges[AlphaToEdgeIndex(c)] != NULL;
}

//TreeNode has a child if any of its child edge pointers are non-null
bool TreeNode::HasChild()
{
    for (int i = 0; i < NumEdges(); i++) {
        if (Edges[i] != NULL) {
            return true;
        }
    }
    return false;
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
    _alphabet = alphabet;
    _root = NULL;
    _input = NULL;
    _numEdges = 0;
    _numInternalNodes = 0;
    _numLeaves = 0;

    Build(&input);
}

SuffixTree::~SuffixTree()
{
    Clear();
}

bool SuffixTree::IsEmpty()
{
    return _root == NULL;
}

void SuffixTree::SetAlphabet(const string& alphabet)
{
    _alphabet = alphabet;
}

void SuffixTree::PrintSize()
{
    if (!IsEmpty()) {
        cout << "From input string of len " << _input->length() << ", suffix-tree has " << _numLeaves << " leaves, " << _numInternalNodes << " internal nodes, and " << _numEdges << " edges." << endl;
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
        delete _root;
        _root = NULL;
        _input = NULL;
        _numEdges = 0;
        _numInternalNodes = 0;
        _numLeaves = 0;
    }
}

//Recursively deletes subtree under some root
//It is assumed the caller will delete the root itself, and set to null as required.
void SuffixTree::_clear(TreeNode* root)
{
    if (root != NULL) {
        for (int i = 0; i < root->NumEdges(); i++) {
            if (root->Edges[i] != NULL) {
                if (root->Edges[i]->Node == NULL) {
                    cout << "WARN _clear() encountered non-null Edge with null Node ptr??" << endl;
                }
                //recursively delete the substree beneath this child
                _clear(root->Edges[i]->Node);
                //at this point, we know the subtree has been recursively deleted; so delete this edge and the node it points to
                delete root->Edges[i]->Node;
                delete root->Edges[i];
                root->Edges[i] = NULL;
            }
        }
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
        if (node->LeafID < 0)
            cout << "Internal node child edges: ";
        else
            cout << "Leaf " << node->LeafID << ": (id=" << _input->substr(node->LeafID,string::npos) << ", depth=" << node->StringDepth << ")";

        //print edge info (internal nodes only)
        for (i = 0; i < node->NumEdges(); i++) {
            if (node->Edges[i] != NULL) {
                string edgeLabel = _input->substr(node->Edges[i]->i, node->Edges[i]->j - node->Edges[i]->i + 1);
                cout << "[" << node->Edges[i]->i << "," << node->Edges[i]->j << "](" << edgeLabel << ")   ";
                nodeQ.push_back(node->Edges[i]->Node);
            }
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
            if (node->Edges[i] != NULL && node->Edges[i]->Node != NULL) {
                cout << node->LeafID << ":" << node->StringDepth << endl;
                _printDfs(node->Edges[i]->Node);
            }
        }
    }
}

/*
Just an assignment requirement. Uses DFS to print the BWT of the suffix tree.
The BWT is defined as traversing dfs, visiting children in lexicographic order,
and printing out the letter corresponding to the leafId minus one. In this
way, the letter preceding any prefix given by some leafId is is the bwt for that suffix.

BWT[i] = s[ leafId - 1 ]  <--The previous letter! Not, say, a zero-based index adjustment.
*/
void SuffixTree::PrintBWT()
{
    cout << "SuffixTree bwt index: " << endl;
    _printBWT(_root);
    cout << endl;
}

//Inner utility for printing the bwt of a suffix tree. BWT is defined as post-order outputting of the leaf ids.
void SuffixTree::_printBWT(TreeNode* node)
{
    if (node != NULL) {
        //print leaf node info
        if (node->LeafID >= 0) {
            #if defined DBG
                if (node->LeafID > _input->length()) {
                    cout << "ERROR LeafID > input length in _printBWT FATAL" << endl;
                }
            #endif

            //print the leaf id character by its bwt index char (leadId - 1)
            int bwtIndex = node->LeafID == 0 ? (_input->length() - 1) : (node->LeafID - 1);
            //cout << node->LeafID << ":" << _input->at(bwtIndex) << endl;
            cout << _input->at(bwtIndex);
        }
        else {
            //else this is an internal node, so print the rest of nodes using dfs, in lexicographic order
            for (int i = 0; i < node->NumEdges(); i++) {
                if (node->Edges[i] != NULL && node->Edges[i]->Node != NULL) {
                    _printBWT(node->Edges[i]->Node);
                }
            }
            cout << flush;
        }
    }
}

//Just an assignment requirement. Given some node, print its immediate children.
void SuffixTree::PrintChildren(TreeNode* node)
{
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (node->Edges[i] != NULL) {
            TreeNode* child = node->Edges[i]->Node;
            cout << "child char: " << child->EdgeIndexToAlpha(i) << endl;
            cout << "child->ID: " << child->LeafID << endl;
            cout << "child->StringDepth: " << child->StringDepth << "\r\n" << endl;
        }
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
    hoppedEdge = v_prime->GetEdge(c);
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
        hoppedEdge = v_prime->GetEdge(c);
        if (hoppedEdge == NULL) {
            //this is valid;
            cout << "ERROR hopped edge NULL in nodeHops" << endl;
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

TreeNode* SuffixTree::_nodeHopsOLD(TreeNode* u, const int suffixIndex)
{
    int beta = u->StringDepth - u->Parent->StringDepth;
    TreeNode* v_prime = u->Parent->SuffixLink;
    char c = _input->at(v_prime->StringDepth + suffixIndex);
    Edge* hoppedEdge = v_prime->GetEdge(c);

    //find lowest internal node below v_prime at or *before* the point |Beta| is exhausted (which may be v_prime itself)
    //while next hop would not make beta negative
    while ((beta - (hoppedEdge->j - hoppedEdge->i + 1)) >= 0) {
        //decrement beta
        beta -= (hoppedEdge->j - hoppedEdge->i + 1);
        //hop the edge
        v_prime = hoppedEdge->Node;
        //get the next edge symbol and its corresponding edge
        c = _input->at(v_prime->StringDepth + suffixIndex);
        hoppedEdge = v_prime->GetEdge(c);
        if (hoppedEdge == NULL) {
            cout << "ERROR hopped edge NULL in nodeHops" << endl;
        }
    }
    //post-loop: landed directly on an internal node (Beta was exausted), or we stopped above an edge too long for the remainder of Beta

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

    //error check
    if (u == NULL || u->Parent == NULL) {
        if (u == NULL)
            cout << "ERROR u == NULL in insertSuffix! Fatal" << endl;
        else if (u->Parent == NULL)
            cout << "ERROR u->parent == NULL in insertSuffix! Fatal" << endl;
        return NULL;
    }

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
    newLeaf->LeafID = suffixIndex;
    newLeaf->Parent = parent;
    //cout << "string depth???" << endl;
    newLeaf->StringDepth = _input->length() - suffixIndex;
    //make the new edge
    Edge* newEdge = (Edge*)malloc(sizeof(Edge));
    newEdge->i = edge_i;
    newEdge->j = _input->length() - 1;
    newEdge->Node = newLeaf;
    //set the outlink of the parent node
    parent->InsertEdge(newEdge, _input->at(newEdge->i));
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
    newInternalNode->LeafID = -1;
    newInternalNode->StringDepth = parent->StringDepth + (edgeSplitIndex - oldEdge->i);
    _numInternalNodes++;

    //split this edge, allocing a new one (continuation) that the new internal node will contain
    Edge* continuation = (Edge*)malloc(sizeof(Edge));
    continuation->i = edgeSplitIndex;
    continuation->j = oldEdge->j;
    continuation->Node = oldEdge->Node;
    continuation->Node->Parent = newInternalNode;
    newInternalNode->InsertEdge(continuation, _input->at(continuation->i));
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
    edge = parent->GetEdge(_input->at(i));
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
                edge = parent->GetEdge(_input->at(i));
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
void SuffixTree::Build(string* s)
{
    int i;
    TreeNode* v;

    //clear existing tree, if any (this nulls _input, so do it first)
    if (!this->IsEmpty()) {
        Clear();
    }

    _input = s;
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
    _root->LeafID = -1;
    _root->StringDepth = 0;
    _numInternalNodes = 1;

    for (i = 0, v = _root; i < s->length(); i++) {
        v = _insertSuffix(v,i);

        //report progress, for very large trees
        if (i % 1000 == 999) {
            cout << "\rInserted " << i << " of " << s->length() << " suffixes, " << (int)(((float)i / (float)s->length()) * 100) << "% complete          " << flush;
        }
    }
    cout << "\r100% complete.                              "<< endl;
}
