/*****************************************************************/
  RAIphy
  Version 2.0, 2/22/2010

  http://bioinfo.unl.edu/
  sway@unl.edu

  This program is freely available for academic use, without any
  warranty.  Commercial distribution of this program, in whole or
  in part, requires prior agreement with the authors.
/*****************************************************************/

1. INSTALL

RAIphy is written in C++, and so should build without error on any
platform with a C++ compiler.  Use the following commands to build
the program.

cd <path_to_source_code_directory>
make clean
make

At this point, the executable "raiphy" (on linux/unix/macosx)
will reside in the src directory as well.
Copy the executable to anywhere in your path.

/*****************************************************************/

2. USAGE

-i : Specify input file. 
Input files should be in FASTA format. 

-I : Specify input folder.
Specifies as input all file in the supplied directory with the file extension specified by "-e". 

-o : Specify output file. (Default: "./outfile")
-d : Specify database file. (Default: "./defaultDb")
-e : Specify file extension. (Default: ".fasta")
-c : Specify classification threshold. (Default: "100")
-h : Print help screen.
-m : Specify run mode. (Default: "0")
	0 : Bin DNA Fragments.
	1 : Bin DNA Fragments with Iterative Refinement.
	2 : Create Database.
	3 : Add Items to Database.