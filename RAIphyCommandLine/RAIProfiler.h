/*****************************************************************/
/*  University of Nebraska-Lincoln                               */
/*  Department of Electrical Engineering                         */
/*  Bioinformatics Group                                         */
/*  Sam Way                                                      */
/*  2/22/10                                                      */
/*****************************************************************/

#ifndef RAIPROFILER_H
#define RAIPROFILER_H

/*********************************************************************************/
/* Included Header Files                                                         */
/*********************************************************************************/
#include "Globals.h"
#include "Parms.h"
#include <map>

/*********************************************************************************/
/* Class Definition                                                              */
/*********************************************************************************/
class RAIProfiler
{
public:
    RAIProfiler();
    string ProfileFastaFile(string filename, INT wordLength, vector<Sequence> *list, BOOLEAN normalize);
    void    PrepareVector(DOUBLE* vector, INT vectorSize, INT wordLength);
    BOOLEAN ProfileNextSequence(fstream *inStream, Sequence *tempSequence, INT wordLength, SequenceCounts *tempCounts);

    map<BYTE, INT> m_hash;
};

#endif // PROFILER_H

/********************************* END OF FILE ***********************************/
