#pragma once
#ifndef _STRING_
#include <string>
#endif

//needed for stat'ing file sizes of very large files
#ifndef __stat64
#include <sys/stat.h>
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

size_t getFilesize(const string& filename);
bool parseAlphabetFile(const string& alphaFile, string& alphabet);
bool parseFastaFile(const string fname, Sequence& output, const string& sigma);
bool fileExists(const string& path);
string filter(const string& s, const string& validChars);