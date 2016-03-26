#include "SuffixTree.hpp"
//#include "FastaParser.hpp"

/*
Main driver for suffix tree testing.
*/
int main(void)
{
    string input;
    SuffixTree st;

    /*
    //build input from test files
    input = "a$";
    st.Build(&input);
    st.PrintBfs();
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


    input = "tattt$";
    st.Build(&input);
    st.PrintBfs();
    

    
    input = "tatttcgtagtcgaaaaatatagctagctcgctgtatagctctgaagcccgtagctaaccggtgaagcgcgt$";
    st.Build(&input);
    st.PrintBfs();

    input = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa$";
    st.Build(&input);
    st.PrintBfs();

    input = "aaaaaaaaaaaaaatttttttttttttaaaaaaaaaaaaacccccccccccccccccgggggggggatcg$";
    st.Build(&input);
    st.PrintBfs();

    

    return 0;
}





