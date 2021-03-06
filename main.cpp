#include "SuffixTree.hpp"

int max(int l, int r)
{
    if (l > r)
        return l;
    return r;
}

/*
Main driver for suffix tree testing.
*/
int main(int argc, char* argv[])
{
    string inputFile, alphaFile;
    Sequence inputSequence;
    string alphabet;
    string input;
    SuffixTree st;
    Sequence seq;
    string sigma = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ$";

    
    /*For api for the hw
    if (argc != 3) {
        cout << "Incorrect number of parameters. Usage: './st [local sequence file] [local alphabet file]'" << endl;
        return -1;
    }
    */

    //inputFile = argv[1];
    //alphaFile = argv[2];

    inputFile = "Test2.txt";
    alphaFile = "english.txt";
    if (fileExists(inputFile)){
        if (fileExists(alphaFile)) {
            if (parseAlphabetFile(alphaFile, alphabet)) {
                if (parseFastaFile(inputFile, inputSequence, alphabet)) {
                    SuffixTree suffixTree(inputSequence.seq, alphabet);
                    
                    //do other suffix tree stuff...
                    //suffixTree.PrintBfs();
                    //suffixTree.PrintBWT();
                    string outputFile = inputFile.substr(0, inputFile.find_first_of('.'));
                    outputFile += "_BWT.txt";
                    cout << "outputting to file: " << outputFile << endl;
                    suffixTree.WriteBWT(outputFile);
                    suffixTree.PrintBWT();
                    //Uncomment the next line to print the BWT to the console, which is awful for large strings.
                    //suffixTree.PrintBWT();
                    suffixTree.PrintNativeSpaceStats();
                    suffixTree.PrintSize();
                    suffixTree.PrintLongestRepeatSubstring();

                    //prg3 testing
                    /*
                    suffixTree.PrepareST();
                    string read = "TTTTTTTT";
                    TreeNode* result = suffixTree.FindLoc(read,1);
                    cout << "found candidates:" << endl;
                    for (int i = result->StartLeafIndex; i < result->EndLeafIndex; i++) {
                        //extract subtrings +/-read.length() to each side of leafIndex
                        cout << inputSequence.seq.substr(max(suffixTree.A[i] - read.length(), 0), read.length() * 2) << endl;
                    }
                    */
                    suffixTree.Clear();
                }
                else {
                    cout << "ERROR could not parse FASTA file" << endl;
                }
            }
            else {
                cout << "ERROR could not parse alphabet file" << endl;
            }
        }
        else {
            cout << "ERROR alphabet file not found: " << alphaFile << endl;
        }
    }
    else {
        cout << "ERROR input file not found: " << inputFile << endl;
    }
    

    /*
    //build input from test files
    input = "a$";
    st.Build(&input);
    //st.PrintBfs();
    //st.PrintDfs();

    input = "aa$";
    st.Build(&input);
    st.PrintBfs();
    //st.PrintDfs();
    */
    /*
    input = "aaaaa$";
    st.Build(&input);
    st.PrintBfs();
    

    input = "tattt$";
    st.Build(&input);
    st.PrintBfs();
    */
    
    /*
    input = "aaaa$";
    st.Build(&input);
    st.PrintBfs();
    */

    /*
    input = "aaaa$";
    st.Build(&input,sigma);
    st.PrintBfs();
    st.Clear();

    
    input = "tattt$";
    st.Build(&input,sigma);
    st.PrintBfs();
    st.Clear();

    
    input = "tatttcgtagtcgaaaaatatagctagctcgctgtatagctctgaagcccgtagctaaccggtgaagcgcgt$";
    st.Build(&input,sigma);
    st.PrintBfs();
    st.PrintBWT();

    
    input = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa$";
    st.Build(&input,sigma);
    st.PrintBfs();

    input = "aaaaaaaaaaaaaatttttttttttttaaaaaaaaaaaaacccccccccccccccccgggggggggatcg$";
    st.Build(&input,sigma);
    st.PrintBfs();
    */  

    



    /*
    //input = "GTGGCGCG$";
    input = "GGCGCG$";
    st.Build(&input,sigma);
    st.PrintBfs();
    */
    /*
    parseFastaFile("Test3.txt", seq, sigma);
    seq.seq += "$";
    st.Build(&seq.seq);
    //st.PrintBfs();
    st.PrintBWT();
    cout << "pause" << endl;
    cin >> seq.seq;

    parseFastaFile("Tomato.txt", seq, sigma);
    seq.seq += "$";
    st.Build(&seq.seq);
    //st.PrintBfs();
    //st.PrintBWT();
    
    parseFastaFile("Yeast.txt", seq, sigma);
    seq.seq += "$";
    st.Build(&seq.seq);
    //st.PrintBfs();
    //st.PrintBWT();
    */
    return 0;
}

