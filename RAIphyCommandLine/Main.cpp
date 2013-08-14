/*****************************************************************/
/*  University of Nebraska-Lincoln                               */
/*  Department of Electrical Engineering                         */
/*  Bioinformatics Group                                         */
/*  Sam Way                                                      */
/*  2/22/10                                                      */
/*****************************************************************/

/*********************************************************************************/
/* Included Header Files */
/*********************************************************************************/
#include "Globals.h"
#include "Parms.h"
#include "RAIphyDatabase.h"
#include "FragmentClassifier.h"

/*********************************************************************************/
/* Private Module Constants */
/*********************************************************************************/
/*********************************************************************************/
/* Private Type Definitions */
/*********************************************************************************/
/*********************************************************************************/
/* Private Module Function Prototypes */
/*********************************************************************************/
void ClassifyFiles(vector<string> *inFilenames, 
                   string outFilename, 
                   INT classificationThreshold, 
                   BOOLEAN updateDatabase); 

/*********************************************************************************/
/* Private Module Variables */
/*********************************************************************************/
static RAIphyDatabase m_Database; 
static FragmentClassifier m_Classifier;

/*********************************************************************************/
/* Public and Private Module Functions */
/*********************************************************************************/
INT main(INT argc, BYTE **argv)
{
    BOOLEAN OkToRun; 

    PARMS_Initialize();
    OkToRun = PARMS_ProcessCommandLine(argc, argv); 
    
    if (OkToRun)
    {
        switch (PARMS_GetRunMode())
        {
            case MODE_BIN:
            case MODE_BIN_WITH_REFINEMENT:
                m_Database.ReadDatabase( PARMS_GetDatabaseFileName() ); 
                if (m_Database.IsValid())
                {
                    ClassifyFiles( PARMS_GetInputFileNames(), 
                                   PARMS_GetOutputFileName(),
                                   PARMS_GetClassificationThreshold(),
                                   (BOOLEAN)(PARMS_GetRunMode() == 
                                   MODE_BIN_WITH_REFINEMENT)); 
                }                   
                break;
                
            case MODE_CREATE_DATABASE:
                std::cout << "Creating database..." << endl; 
                m_Database.CreateDatabase( PARMS_GetInputFileNames(),
                                           PARMS_GetDatabaseFileName(),
                                           DEFAULT_OLIGMER_SIZE ); 
                break;
                
            case MODE_ADD_TO_DATABASE:
                m_Database.ReadDatabase( PARMS_GetDatabaseFileName() ); 
                if (m_Database.IsValid())
                {
                    m_Database.AddItems( PARMS_GetInputFileNames() ); 
                }
                break;
            default:
                std::cout << "Invalid run mode. Exiting program." << endl; 
                break;
        }
    }

    return 0; 
}

/*********************************************************************************/

void ClassifyFiles(vector<string> *inFilenames, 
                   string outFilename, 
                   INT classificationThreshold, 
                   BOOLEAN updateDatabase)
{
    BOOLEAN successful = TRUE;
    INT updates = 0;
    ULONG i, j;
    INT* updatedCounts;
    DOUBLE** updatedDatabase;
    
    if (updateDatabase)
    {
        updatedCounts = new INT[m_Database.RAIProfiles.size()];
        updatedDatabase = new DOUBLE*[m_Database.RAIProfiles.size()];
        for (i = 0; i < m_Database.RAIProfiles.size(); i++) 
        {
            updatedDatabase[i] = new DOUBLE[(INT)pow(NUM_BASES,m_Database.WordLength())];
        }

        for (updates=0; updates<DB_MAX_UPDATES; updates++)
        {
            std::cout << "Pass #" << updates+1 << ":" << endl; 
            for (i = 0; i < m_Database.RAIProfiles.size(); i++)
            {
                updatedCounts[i] = 0;
                for (j = 0; j < (INT)pow(NUM_BASES,m_Database.WordLength()); j++) updatedDatabase [i][j] = 0;
            }
            
            successful = m_Classifier.ClassifyFilesAndUpdateDatabase( &m_Database,
                                                                     inFilenames,
                                                                     outFilename,
                                                                     updatedCounts,
                                                                     updatedDatabase );
            if (!successful)
            {
                cout << "ERROR: Could not classify file(s)." << endl;
                return;
            }
        }
        
        for (i = 0; i < m_Database.RAIProfiles.size(); i++)
        { 
            delete [] updatedDatabase[i];
        }
        delete [] updatedDatabase;
        delete [] updatedCounts;
    }
    
    std::cout << "Pass #" << updates+1 << ":" << endl; 
    m_Classifier.ClassifyFiles( &m_Database,
                               inFilenames,
                               outFilename,
                               classificationThreshold );
    if (updateDatabase)
    {
        m_Database.SaveDatabase();
    }
}
