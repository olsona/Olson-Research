/*****************************************************************/
/*  University of Nebraska-Lincoln                               */
/*  Department of Electrical Engineering                         */
/*  Bioinformatics Group                                         */
/*  Sam Way                                                      */
/*  2/22/10                                                      */
/*****************************************************************/

/*********************************************************************************/
/* Included Header Files                                                         */
/*********************************************************************************/
#include <fstream>
#include "RAIphyDatabase.h"

using namespace std;

/*********************************************************************************/
/* Private Module Constants                                                      */
/*********************************************************************************/

/*********************************************************************************/
/* Constructors / Destructors                                                    */
/*********************************************************************************/
RAIphyDatabase::RAIphyDatabase()
{
    RAIProfiles.clear();
    m_fileName = ""; 
    m_isValid = FALSE;
    m_updateCount = 0; 
    m_wordLength = -1;
}

/*********************************************************************************/

RAIphyDatabase::~RAIphyDatabase()
{
}

/*********************************************************************************/
/* Private / Public Functions                                                    */
/*********************************************************************************/
BOOLEAN RAIphyDatabase::IsValid(void)
{
    return m_isValid;
}

/*********************************************************************************/

ULONG RAIphyDatabase::WordLength(void)
{
    return m_wordLength;
}

/*********************************************************************************/

INT RAIphyDatabase::GetIndex(string modelName)
{

    ULONG i;

    for (i=0; i<RAIProfiles.size(); i++)
    {
        if (RAIProfiles.at(i).sequenceName.find(modelName) != string::npos)
        {
            return i;
        }
    }
    
    return -1;
}

/*********************************************************************************/

void RAIphyDatabase::ReadDatabase(string fileName)
{
    ifstream inFile(fileName.c_str());
    vector<string> tempStringParts;
    string tempString, statusString;
    Sequence tempSequence;
    ULONG i, j;

    if (inFile)
    {
        try
        {
            RAIProfiles.clear();
            m_isValid = FALSE;
            j = 0;

            getline(inFile, tempString);
            
            PARMS_GetStringParts(tempString, &tempStringParts, SEPARATOR); 
            if (tempStringParts.size() != 2 || tempStringParts.at(0) != HEADER_WORD_LENGTH)
            {
                m_isValid = FALSE;
                std::cout << "ERROR: Problem reading database. Database missing word length field." << endl;
                return;
            }
            
            m_wordLength = atoi(tempStringParts.at(1).c_str());
            m_fileName = fileName;
            getline(inFile, tempString); 

            while(!inFile.eof())
            {
                std::cout << "Loading item " << j << " of database..." << "...           \r";

                PARMS_GetStringParts(tempString, &tempStringParts, SEPARATOR); 
                tempSequence.sequenceName = tempStringParts.at(0);
                tempSequence.sequenceVector.clear();
                for (i=1; i < (ULONG)tempStringParts.size(); i++)
                {
                    tempSequence.sequenceVector.push_back(atof(tempStringParts.at(i).c_str()));
                        
                }
                RAIProfiles.push_back(tempSequence);
                getline(inFile, tempString); 
                j++;
            }
            
            m_isValid = TRUE;
            std::cout << "Database loaded successfully!  (" << RAIProfiles.size()
                        << " items)." << endl;
        }
        catch (INT e)
        {
            m_isValid = FALSE;
            std::cout << "ERROR: Problem reading database." << endl;
        }
    }
    else
    {
        m_isValid = FALSE;
        std::cout << "ERROR: Could not open database file for reading." << endl;
    }
    
    inFile.close();
}

/*********************************************************************************/

BOOLEAN RAIphyDatabase::SaveDatabase(void)
{
    BOOLEAN updateStatus = TRUE;
    ULONG i, j;
    ofstream outFile; 
    
    std::cout << "Saving updated database..." << endl; 
    
    outFile.open(m_fileName.c_str()); 
    
    if (outFile.is_open())
    {
        try
        {
            outFile << HEADER_WORD_LENGTH << SEPARATOR << m_wordLength << endl;
            
            for (i=0; i < (ULONG)RAIProfiles.size(); i++)
            {
                outFile << RAIProfiles.at(i).sequenceName;
                for (j=0; j < (ULONG)RAIProfiles.at(i).sequenceVector.size(); j++)
                {
                    outFile << SEPARATOR << RAIProfiles.at(i).sequenceVector.at(j);
                }
                outFile << endl;
            }
            
            m_updateCount = 0;
            cout << "Database successfully updated." << endl;
        }
        catch (INT e)
        {
            cout << "ERROR: Problem updating database file." << endl;
        }
        outFile.close();
    }
    else
    {
        cout << "ERROR: Could not open database file for updating." << endl;
        updateStatus = FALSE;
    }
    return updateStatus;
}

/*********************************************************************************/

void RAIphyDatabase::AddItems(vector<string> *inFileNames)
{
    ULONG i;
    vector<Sequence> tempList;
    
    std::cout << "Adding new items to the database..." << endl;
    
    for (i=0; i<inFileNames->size(); i++)
    {
        cout << "Adding sequences from file " << i+1 << " of " << inFileNames->size() 
             <<  "..." << endl; 
        m_profiler.ProfileFastaFile(inFileNames->at(i), m_wordLength, &tempList, TRUE);
    }
    
    for (i=0; i<tempList.size(); i++)
    {
        RAIProfiles.push_back(tempList[i]);
    }
    
    std::cout << "Added " << tempList.size() << " items to the database." << endl;     
    SaveDatabase();
}


/*********************************************************************************/

void RAIphyDatabase::CreateDatabase(vector<string> *inFilenames, string outFilename, INT wordLength)
{
    ULONG i, j;
    ofstream outFile;
    
    outFile.open(outFilename.c_str()); 
    
    if (outFile.is_open())
    {
        m_isValid = FALSE; 
        RAIProfiles.clear();
        
        try
        {
            outFile << HEADER_WORD_LENGTH << SEPARATOR << wordLength << endl;
            
            for (i=0; i < (ULONG)inFilenames->size(); i++)
            {
                std::cout << "Profiling sequences from file " << i+1 << " of " << inFilenames->size()
                        << "..." << endl; 
                m_profiler.ProfileFastaFile(inFilenames->at(i), wordLength, &RAIProfiles, TRUE);
            }
            
            for (i=0; i < (ULONG)RAIProfiles.size(); i++)
            {
                std::cout << "Writing database vector " << i+1 << " of " << RAIProfiles.size() 
                    << "..." << endl; 
                outFile << RAIProfiles.at(i).sequenceName;
                for (j=0; j < (ULONG)RAIProfiles.at(i).sequenceVector.size(); j++)
                {
                    outFile << SEPARATOR << RAIProfiles.at(i).sequenceVector.at(j);
                }
                outFile << endl;
            }
            
            m_isValid = TRUE;
            m_wordLength = wordLength;
            m_fileName = outFilename;
            m_updateCount = 0;
            std::cout << "Database successfully created (" << RAIProfiles.size() << " items)." << endl;
        }
        catch (INT e)
        {
            m_isValid = FALSE;
            std::cout << "ERROR: Problem creating database file." << endl;
        }

        outFile.close();
    }
    else
    {
        std::cout << "ERROR: Could not open file to create database." << endl;
    }
}

/*********************************************************************************/

BOOLEAN RAIphyDatabase::UpdateDatabase(string outputFilename)
{
    ifstream inFile(outputFilename.c_str());
    BYTE tempChar = '\0'; 
    string tempString;
    vector<string> tempStringParts;
    INT counts[RAIProfiles.size()];
    ULONG runningSums[RAIProfiles.size()];
    DOUBLE tempVector[(INT)pow(NUM_BASES,m_wordLength)];
    ULONG i, tempIndex, updates;
    BOOLEAN updateSuccess;
    
    for (i=0; i<RAIProfiles.size(); i++)
    {
        counts[i] = 0;
        runningSums[i] = 0;
    }
    
    if (inFile)
    {
        try
        {
            updates = 0;
            getline(inFile, tempString);
            PARMS_GetStringParts(tempString, &tempStringParts, SEPARATOR); 
            if (tempStringParts.size() != 2 || tempStringParts.at(0) != HEADER_WORD_LENGTH)
            {
                std::cout << "ERROR: Problem updating database. Output file missing word length field." << endl;
                return FALSE;
            }
            if (atoi(tempStringParts.at(1).c_str()) != (INT)this->m_wordLength)
            {
                std::cout << "ERROR: Could not update database.  Output file differs in word length." << endl;
                return FALSE;
            }

            tempChar = tolower(inFile.get());
            while(inFile.good())
            {
                while (inFile.good() && tempChar != '>') tempChar = tolower(inFile.get());
                if (inFile.good())
                {
                    getline(inFile, tempString);            // Get the header line (">Example_Label")
                    getline(inFile, tempString);            // Get the matching db vector ("Vector_Label")
                    tempIndex = GetIndex(tempString);       // Get index corresponding to the db vector
                    getline(inFile, tempString);            // Get the fragment length
                    
                    if (tempIndex >= 0)
                    {
                        counts[tempIndex] += 1;
                        runningSums[tempIndex] += strtoul(tempString.c_str(), NULL, 0);
                    }
                    tempChar = tolower(inFile.get());
                }
            }
        }
        catch (INT e)
        {
            std::cout << "ERROR: Could not update database.  An unexpected error has occurred." << endl;
            return FALSE;
        }
        inFile.close();
        
        for (i=0; i<RAIProfiles.size(); i++)
        {
            if (runningSums[i] >= UPDATE_THRESHOLD)
            {
                updateSuccess = UpdateVector(i, tempVector, (INT)pow(NUM_BASES,m_wordLength), outputFilename);
                updates++;
                if (!updateSuccess)
                {
                    std::cout << "ERROR: Could not update database.  Failed to update vector." << endl;
                    return FALSE;
                }
            }
        }
        
        std::cout << "Updated " << updates << " item(s) in the database." << endl; 
        return TRUE;
    }
    return FALSE;
}

/*********************************************************************************/

void RAIphyDatabase::UpdateVector(INT vectorIndex, DOUBLE *tempVector)
{
    INT i;
    
    m_profiler.PrepareVector(tempVector, (INT)pow(NUM_BASES,m_wordLength), m_wordLength);
    
    for ( i = 0; i < (INT)pow(NUM_BASES,m_wordLength); i++)
    {
        RAIProfiles[vectorIndex].sequenceVector[i] = tempVector[i];
    }
}

/*********************************************************************************/

BOOLEAN RAIphyDatabase::UpdateVector(INT vectorIndex, DOUBLE *tempVector, INT vectorLength, string outputFilename)
{
    ifstream inFile(outputFilename.c_str());
    BYTE tempChar = 0; 
    string tempString;
    vector<string> tempStringParts;
    INT i, j, tempInt;
    DOUBLE tempDouble;
    
    for (i=0; i<vectorLength; i++) tempVector[i] = 0;
    
    if (inFile)
    {
        try
        {
            while(inFile.good())
            {
                while (inFile.good() && tempChar != '>') tempChar = tolower(inFile.get());
                if (inFile.good())
                {
                    getline(inFile, tempString);                    // Get rest of the header line
                    getline(inFile, tempString);                    // Get the db vector name
                    tempString = PARMS_Trim(tempString); 
                    
                    if (RAIProfiles[vectorIndex].sequenceName.find(tempString)!=string::npos)
                    {
                        getline(inFile, tempString);           // Get fragment length
                        getline(inFile, tempString);           // Get index:count fields
                        PARMS_GetStringParts(tempString, &tempStringParts, SEPARATOR);
                        
                        for (j=0; j<(INT)(tempStringParts.size()/2); j++)
                        {
                            tempInt = atoi(tempStringParts[j*2].c_str());
                            tempDouble = atof(tempStringParts[(j*2)+1].c_str());
                            tempVector[tempInt] += tempDouble;
                        }
                    }
                }
            }
            
            m_profiler.PrepareVector(tempVector, vectorLength, m_wordLength);
            for (i=0; i<vectorLength; i++) 
            {
                RAIProfiles[vectorIndex].sequenceVector[i] = tempVector[i];
            }
            return TRUE;
        }
        catch (INT e) {  }
        inFile.close();
    }
    return FALSE;
}

/********************************* END OF FILE ***********************************/

