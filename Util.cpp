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
    int seqnum, thousands = 0;
    ifstream inputFile;
    string line;
    int fsize;

    if (fileExists(fname)) {
        fsize = getFilesize(fname);
        output.seq.reserve(fsize);

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

                if ((output.seq.size() & 0x00001FFF) == 0) {
                    thousands = output.seq.size() / 5000;
                    cout << "\r" << (int)(100 * (float)output.seq.size() / (float)fsize) << "% complete                " << flush;
                }
            }
            cout << "\r100% complete" << endl;
            inputFile.close();
            success = true;
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

//Dumb utility for parsing a file that containing 2 fasta files concatenated, as required in prg1.
void parse2FastaFile(const string fname, Sequence& seq1, Sequence& seq2)
{
    int seqnum;
    ifstream myReadFile;
    myReadFile.open(fname);
    string line;
    string validBases = "atcgATCG";

    string& s1 = seq1.seq;
    string& s2 = seq2.seq;

    if (myReadFile.is_open()) {
        //init
        s1.clear();
        s2.clear();
        seqnum = 0;

        while (getline(myReadFile, line)) {
            //check if line begins with '>'; if so, advance state, and grab the description line
            if (line.length() > 0 && line[0] == '>') {
                if (seqnum == 0)
                    seq1.desc = line.substr(1, line.length());
                else
                    seq2.desc = line.substr(1, line.length());
                seqnum++;
            }
            else {
                switch (seqnum) {
                    //parsing sequence 1
                case 1:
                    s1 += filter(line, validBases);
                    break;
                    //parsing sequence 2
                case 2:
                    s2 += filter(line, validBases);
                    break;
                default:
                    cout << "ERROR unknown seqnum: " << seqnum << endl;
                    break;
                }
            }
        }
    }
    myReadFile.close();

    //cout << "Parsed sequences:" << endl;;
    //cout << "s1 >>" << s1 << endl;
    //cout << "s2 >>" << s2 << endl;
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
