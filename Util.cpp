#include "Util.hpp"

//returns a string filtered of all but characters in validChars
string filter(const string& s, const string& validChars)
{
    string filtered;

    for (int i = 0; i < s.length(); i++) {
        if (validChars.find(s[i]) != string::npos) {
            filtered += s[i];
        }
    }

    return filtered;
}

/*
Parses an alphabet (sigma) from a file containing a single line of space-delimited chars: "A B C D E..."
*/
bool parseAlphabetFile(const string& alphaFile, string& alphabet)
{
    bool success = false;
    ifstream inputFile;
    string line;

    if (fileExists(alphaFile)) {
        inputFile.open(alphaFile);
        if (inputFile.is_open()) {
            //get the single line in the alphabet file
            if (getline(inputFile, line) && line.length() > 0) {
                alphabet.clear();
                for (int i = 0; i < line.length(); i++) {
                    if (line[i] != ' ') {
                        alphabet += line[i];
                    }
                }
                success = alphabet.length() > 0;

                if (!success) {
                    cout << "ERROR alphabet file contained only blanks" << endl;
                }
            }
            else {
                cout << "ERROR alphabet file empty or corrupted." << endl;
            }
            inputFile.close();
        }
    }
    else {
        cout << "ERROR, file not found: " << alphaFile << endl;
    }

    return success;
}

/*
Parses a simple, single-sequence fasta file, filling in a simple Sequence data structure.
Look up fasta format; its a simple file format with a description line prepended '>' followed by lines of a raw sequence.

sigma: the valid chars; anything not in these will be omitted
*/
bool parseFastaFile(const string fname, Sequence& output, const string& sigma)
{
    bool success = false;
    int seqnum;
    ifstream inputFile;
    string line;

    if (fileExists(fname)) {
        inputFile.open(fname);
        if (inputFile.is_open()) {
            //init
            output.seq.clear();
            seqnum = 0;
            cout << "Parsing sequence from file: " << fname << ". This may require a few seconds, for input files > 200kb." << endl;
            while (getline(inputFile, line)) {
                //check if line begins with '>'; if so, advance state, and grab the description line
                if (seqnum == 0 && line.length() > 0 && line[0] == '>') {
                    output.desc = line.substr(1, line.length());
                    seqnum++;
                }
                else {
                    output.seq += filter(line, sigma);
                }
                success = true;
            }
            inputFile.close();
        }
        else {
            cout << "ERROR file could not be opened! " << fname << endl;
        }
    }
    else {
        cout << "ERROR, file not found: " << fname << endl;
    }

    return success;
}

size_t getFilesize(const string& filename)
{
    struct stat st;

    if (stat(filename.c_str(), &st) != 0) {
        return 0;
    }

    return st.st_size;
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
