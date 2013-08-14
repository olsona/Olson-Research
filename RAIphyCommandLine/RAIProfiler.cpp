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
#include "RAIProfiler.h"

/*********************************************************************************/
/* Private Module Constants                                                      */
/*********************************************************************************/

/*********************************************************************************/
/* Constructors / Destructors                                                    */
/*********************************************************************************/
RAIProfiler::RAIProfiler()
{
    m_hash.clear();
    m_hash['a'] = BASE_VALUE_A;
    m_hash['g'] = BASE_VALUE_G;
    m_hash['c'] = BASE_VALUE_C;
    m_hash['t'] = BASE_VALUE_T;
}

/*********************************************************************************/
/* Private / Public Functions                                                    */
/*********************************************************************************/
string RAIProfiler::ProfileFastaFile(string filename, INT wordLength, vector<Sequence> *list, BOOLEAN normalize)
{
    string errorMessage = "";
    BYTE tempChar = 0; 
    Sequence tempSequence;
    ifstream inFile(filename.c_str()); 
    DOUBLE databaseVector[(INT)pow(NUM_BASES,wordLength)];
    ULONG readBuffer[wordLength];
    ULONG head, tail, indexForward, indexReverse;

    if (inFile)
    {
        while(inFile.good())
        {
            if (inFile.good())
            {   
                while (inFile.good() && tempChar != '>') tempChar = tolower(inFile.get());
                
                if (!inFile.good() || tempChar != '>')                                          // Read until header line
                {
                    inFile.close();
                    return "ERROR: Invalid Fasta file.";
                }
                 
                getline(inFile, tempSequence.sequenceName);                                     // Extract sequence name
                tempSequence.sequenceVector.clear();
                tempSequence.sequenceLength = 0;

                tail = 0;
                indexForward = 0;
                indexReverse = 0;
                
                for (head = 0; head < (INT)pow(NUM_BASES,wordLength); head++)                   // Initialize DB vector
                {
                    databaseVector[head] = 0;
                }
                
                tempChar = tolower(inFile.get());
                for (head = 0; head < (ULONG)wordLength; head++)                                // Grab the first n-mer
                {
                    readBuffer[head] = m_hash[tempChar];
                    indexForward += (INT)(readBuffer[head]*pow(NUM_BASES,wordLength-1-head));
                    indexReverse += (INT)((NUM_BASES-1-readBuffer[head])*pow(NUM_BASES,head));  // Reverse gets complimented
                    tempChar = tolower(inFile.get());
                }
                
                tempSequence.sequenceLength = wordLength;
                databaseVector[indexForward] += 1;                                              // Increment DB accordingly
                databaseVector[indexReverse] += 1;
                
                while(inFile.good() && tempChar != '>')
                {
                    indexForward -= (INT)(readBuffer[tail]*pow(NUM_BASES,wordLength-1));        // Remove last base
                    indexReverse -= (NUM_BASES-1-readBuffer[tail]);
                    
                    head++;  if (head >= (ULONG)wordLength) head = 0;
                    tail++;  if (tail >= (ULONG)wordLength) tail = 0;
                    
                    readBuffer[head] = m_hash[tempChar];
                    indexForward = (INT)(indexForward*NUM_BASES) + readBuffer[head];            // Shift and add new base
                    indexReverse = ((INT)(indexReverse/NUM_BASES) +
                                    (NUM_BASES-1-readBuffer[head])*pow(NUM_BASES, wordLength-1));
                    
                    databaseVector[indexForward] += 1;                                          // Increment vectors
                    databaseVector[indexReverse] += 1;
                    tempSequence.sequenceLength += 1;
                    
                    do
                    {
                        tempChar = tolower(inFile.get());                                       // Grab next character
                    } while(inFile.good() && m_hash.find(tempChar) == m_hash.end() && tempChar != '>');
                }
                
                if (normalize) PrepareVector(databaseVector, (INT)pow(NUM_BASES,wordLength), wordLength);
                
                tail = 0;
                for (head = 0; head < (INT)pow(NUM_BASES,wordLength); head++)                   // Copy results
                {
                    tempSequence.sequenceVector.push_back((databaseVector[head]));
                }
                list->push_back(tempSequence);                                                  // Add results to list
            }
        }
    }
    else 
    {
        std::cout << "ERROR: Could not open file." << endl;
    }


    inFile.close();
    return errorMessage;
}

/*********************************************************************************/

void RAIProfiler::PrepareVector(DOUBLE* vector, INT vectorSize, INT wordLength)
{
    DOUBLE normalizeVector[(INT)pow(NUM_BASES,wordLength-1)];
    ULONG tempValue;
    INT i, j;

    for (i=0; i<vectorSize; i+=NUM_BASES)
    {
        tempValue = 0;
        for (j=0; j<NUM_BASES; j++)
        {
            tempValue += vector[i+j];
        }

        normalizeVector[i/NUM_BASES] = tempValue;                       // Build second vector

        if (tempValue != 0)
        {
            for (j=0; j<NUM_BASES; j++)
            {
                vector[i+j] = vector[i+j] / tempValue;                  // Normalize first vector
            }
        }
    }

    for (i=0; i < (INT)pow(NUM_BASES,wordLength-1); i+=NUM_BASES)
    {
        tempValue = 0;
        for (j=0; j<NUM_BASES; j++)
        {
            tempValue += normalizeVector[i+j];
        }

        if (tempValue != 0)                                             // Normalize second vector
        {
            for (j=0; j<NUM_BASES; j++)
            {
                normalizeVector[i+j] = normalizeVector[i+j] / tempValue;
            }
        }
    }

    for (i=0; i<vectorSize; i+=NUM_BASES)
    {
        for (j=0; j<NUM_BASES; j++)
        {
            if (vector[i+j] != 0 && normalizeVector[i/NUM_BASES] != 0)
            {
                vector[i+j] = 11*log(vector[i+j]) + log(normalizeVector[i/NUM_BASES]);
            }
            else
            {
                vector[i+j] = ZERO_ENTRY_VALUE;
            }
        }
    }
}

/*********************************************************************************/

BOOLEAN RAIProfiler::ProfileNextSequence(fstream *inStream, Sequence *tempSequence, INT wordLength, SequenceCounts *tempCounts)
{
    ULONG readBuffer[wordLength];
    BYTE tempChar;
    DOUBLE tempVector[(INT)pow(NUM_BASES,wordLength)];
    ULONG head, tail, indexForward, indexReverse;

    try
    {
        tempSequence->sequenceName = "";
        tempSequence->sequenceLength = 0;
        tempSequence->sequenceVector.clear();

        tempCounts->nonzeroCount.clear();
        tempCounts->nonzeroIndex.clear();

        tail = 0;
        indexForward = 0;
        indexReverse = 0;
        for (head = 0; head < (ULONG)wordLength; head++) readBuffer[head] = 0;
        for (head = 0; head < (ULONG)pow(NUM_BASES,wordLength); head++) tempVector[head] = 0;

        getline((*inStream), tempSequence->sequenceName); 

        tempChar = tolower(inStream->get());
        for (head = 0; head < (ULONG)wordLength; head++)                                // Grab the first n-mer
        {
            readBuffer[head] = m_hash[tempChar];
            indexForward += (INT)(readBuffer[head]*pow(NUM_BASES,wordLength-1-head));
            indexReverse += (INT)((NUM_BASES-1-readBuffer[head])*pow(NUM_BASES,head));  // Reverse gets complimented
            tempChar = tolower(inStream->get());
        }
        tempSequence->sequenceLength = wordLength;

        tempVector[indexForward] += 1;                                              // Increment DB accordingly
        tempVector[indexReverse] += 1;

        while(inStream->good() && tempChar != '>')
        {
            indexForward -= (INT)(readBuffer[tail]*pow(NUM_BASES,wordLength-1));        // Remove last base
            indexReverse -= (NUM_BASES-1-readBuffer[tail]);

            head++;  if (head >= (ULONG)wordLength) head = 0;
            tail++;  if (tail >= (ULONG)wordLength) tail = 0;

            readBuffer[head] = m_hash[tempChar];
            indexForward = (INT)(indexForward*NUM_BASES) + readBuffer[head];            // Shift and add new base
            indexReverse = ((INT)(indexReverse/NUM_BASES) +
                (NUM_BASES-1-readBuffer[head])*pow(NUM_BASES, wordLength-1));

            tempVector[indexForward] += 1;                                          // Increment vectors
            tempVector[indexReverse] += 1;
            tempSequence->sequenceLength += 1;

            do
            {
                tempChar = tolower(inStream->get());                                   // Grab next character
            } while(inStream->good() && m_hash.find(tempChar) == m_hash.end() && tempChar != '>');
        }

        for (head = 0; head < (INT)pow(NUM_BASES,wordLength); head++)                   // Copy results
        {
            if (tempVector[head] > 0)
            {
                tempCounts->nonzeroIndex.push_back(head);
                tempCounts->nonzeroCount.push_back(tempVector[head]);
            }
        }

        return TRUE;
    }
    catch (int e)
    {
        return FALSE;
    }
}

/********************************* END OF FILE ***********************************/
