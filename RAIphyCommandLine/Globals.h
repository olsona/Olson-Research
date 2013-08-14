/*****************************************************************/
/*  University of Nebraska-Lincoln                               */
/*  Department of Electrical Engineering                         */
/*  Bioinformatics Group                                         */
/*  Sam Way                                                      */
/*  2/22/10                                                      */
/*****************************************************************/

#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

/******************************************************************/
/* Public Module Constants                                        */
/******************************************************************/
#define VERSION_MAJOR            2
#define VERSION_MINOR            0
#define VERSION_DATE             "02/22/11"

#define FALSE                    false
#define TRUE                     true 

#define COIN_FLIP                ((rand() & 0x01) == 0x01)
#define MIN(a,b)                 ((a)<(b)?(a):(b))
#define MAX(a,b)                 ((a)>(b)?(a):(b))

/******************************************************************/
/* Public Type Definitions                                        */
/******************************************************************/
typedef char            BYTE;
typedef unsigned char   UBYTE;

typedef short           SHORT;
typedef unsigned short  USHORT;

typedef long            LONG;
typedef unsigned long   ULONG;

typedef float           FLOAT;
typedef double          DOUBLE;

typedef int             INT;
typedef unsigned int    UINT;

typedef bool            BOOLEAN;

struct Sequence
{
    string sequenceName;
    vector<DOUBLE> sequenceVector;
    ULONG sequenceLength;
};

struct SequenceCounts
{
    vector<ULONG> nonzeroIndex;
    vector<ULONG> nonzeroCount;
};

/******************************************************************/
/* Public Module Variables */
/******************************************************************/
/******************************************************************/
/* Public Module Function Prototypes */
/******************************************************************/

#endif // GLOBALS_H

/********************************* END OF FILE ***********************************/
