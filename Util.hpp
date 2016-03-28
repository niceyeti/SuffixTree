#pragma once
#ifndef _STRING_
#include <string>
#endif

#ifndef _IOSTREAM_
#include <iostream>
#endif

#ifndef _FSTREAM_
#include <fstream>
#endif

using namespace std;

typedef struct sequence {
    string seq;
    string desc;
}Sequence;

bool parseAlphabetFile(const string& alphaFile, string& alphabet);
bool parseFastaFile(const string fname, Sequence& output, const string& sigma);
bool fileExists(const string& path);
string filter(const string& s, const string& validChars);