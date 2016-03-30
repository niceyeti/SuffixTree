This folder contains all my source code and outputs. I've manually moved sample test outputs into their own folder.

To compile, just use 'sh ./compile'. The compile folder will compile for your machine, as long as you have g++ with c++ 11 and O2 options.

Usage: './st [local fasta input file] [local alphabet file]

The program is currently configured to run and output statistics for build time, space, and the longest repeated substring.
It will also output the BWT of the input string to a file named '[input file name]_BWT.txt'. THIS FILE MUST EXIST BEFORE
RUNNING THE PROGRAM. So for example, if the input file is "someString.txt", then the BWT output will be written to
"someString_BWT.txt", a file which must exist beforehand.

You can also view and comment/uncomment certain tasks in main. The current procedure does the following:
	-build the tree
	-write the BWT to file as described above
	-output space stats
	-output maximal repeating substring
	-clear/destroy tree
	
To execute all tests:
sh ./compile
sh ./runTests