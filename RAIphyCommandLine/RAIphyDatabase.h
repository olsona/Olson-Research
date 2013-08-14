/*****************************************************************/
/*  University of Nebraska-Lincoln                               */
/*  Department of Electrical Engineering                         */
/*  Bioinformatics Group                                         */
/*  Sam Way                                                      */
/*  2/22/10                                                      */
/*****************************************************************/

#ifndef RAIPHYDATABASE_H
#define RAIPHYDATABASE_H

/*********************************************************************************/
/* Included Header Files                                                         */
/*********************************************************************************/
#include "Globals.h"
#include "Parms.h"
#include "RAIProfiler.h"

/*********************************************************************************/
/* Class Defintion                                                               */
/*********************************************************************************/

class RAIphyDatabase 
{
    public:
        RAIphyDatabase(); 
        ~RAIphyDatabase(); 

        vector<Sequence> RAIProfiles;

        // Public Functions 
        BOOLEAN IsValid(void); 
        ULONG WordLength(void); 
        INT GetIndex(string modelName); 
        void ReadDatabase(string fileName); 
        BOOLEAN SaveDatabase(void); 
        void AddItems(vector<string> *inFileNames); 
        void CreateDatabase(vector<string> *inFilenames, string outFilename, INT wordLength);
        BOOLEAN UpdateDatabase(string outputFilename);
        void UpdateVector(INT vectorIndex, DOUBLE *tempVector);
        BOOLEAN UpdateVector(INT vectorIndex, DOUBLE *tempVector, INT vectorLength, string outputFilename); 
    
    private:
        string m_fileName; 
        BOOLEAN m_isValid;
        ULONG m_updateCount; 
        ULONG m_wordLength; 
        
        RAIProfiler m_profiler;
    
}; 

#endif // CLASSIFIERDATABASE_H

/********************************* END OF FILE ***********************************/
