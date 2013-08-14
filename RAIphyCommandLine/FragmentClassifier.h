/*****************************************************************/
/*  University of Nebraska-Lincoln                               */
/*  Department of Electrical Engineering                         */
/*  Bioinformatics Group                                         */
/*  Sam Way                                                      */
/*  2/22/10                                                      */
/*****************************************************************/

#ifndef FRAGMENTCLASSIFIER_H
#define FRAGMENTCLASSIFIER_H

/*********************************************************************************/
/* Included Header Files                                                         */
/*********************************************************************************/
#include "Globals.h"
#include "RAIphyDatabase.h"
#include "RAIProfiler.h"

/*********************************************************************************/
/* Class Definition                                                              */
/*********************************************************************************/

class FragmentClassifier
{

public:
    FragmentClassifier();
    INT ClassifySequence ( RAIphyDatabase *database, SequenceCounts* tempCounts, DOUBLE *matchScore);
    BOOLEAN ClassifyFiles ( RAIphyDatabase *database,
                            vector<string> *inFilenames,
                            string outFilename,
                            INT classificationThreshold);
    BOOLEAN ClassifyFilesAndUpdateDatabase ( RAIphyDatabase *database,
                                             vector<string> *inFilenames,
                                             string outFilename,
                                             INT* updatedCounts,
                                             DOUBLE** updatedDatabase );

private:
    RAIProfiler m_profiler;

    void CleanOutput(string outputFilename, DOUBLE* bestScores, RAIphyDatabase* database, INT classificationThreshold);

};

#endif // FRAGMENTCLASSIFIER_H

/********************************* END OF FILE ***********************************/

