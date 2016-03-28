#include "Util.hpp"

//returns a string filtered of all but characters in validChars
string filter(const string& s, const string& validChars)
{
    string filtered;

    for (int i = 0; i < s.length(); i++) {
        if (validChars.find(tolower(s[i])) != string::npos) {
            filtered += s[i];
        }
    }

    return filtered;
}

/*
Parses a simple, single-sequence fasta file, filling in a simple Sequence data structure.
Look up fasta format; its a simple file format with a description line prepended '>' followed by lines of a raw sequence.

sigma: the valid chars; anything not in these will be omitted
*/
void parseFastaFile(const string fname, Sequence& output, const string& sigma)
{
    int seqnum;
    ifstream inputFile;
    string line;

    if (fileExists(fname)) {
        inputFile.open(fname);
        if (inputFile.is_open()) {
            //init
            output.seq.clear();
            seqnum = 0;

            while (getline(inputFile, line)) {
                //check if line begins with '>'; if so, advance state, and grab the description line
                if (seqnum == 0 && line.length() > 0 && line[0] == '>') {
                    output.desc = line.substr(1, line.length());
                    seqnum++;
                }
                else {
                    output.seq += filter(line, sigma);
                }
            }
        }
        inputFile.close();
    }
    else {
        cout << "ERROR, file not found: " << fname << endl;
    }
}

bool fileExists(const string& path)
{
    ifstream myStream(path);

    bool fileAccessible = !myStream.fail();
    if (fileAccessible) {
        myStream.close();
    }

    return fileAccessible;
}
