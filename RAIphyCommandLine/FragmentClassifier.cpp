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
#include "FragmentClassifier.h"
#include <fstream>

/*********************************************************************************/
/* Private Module Constants                                                      */
/*********************************************************************************/
#define NOT_CLASSIFIED "(unknown)"

/*********************************************************************************/
/* Constructors / Destructors                                                    */
/*********************************************************************************/
FragmentClassifier::FragmentClassifier() {}

/*********************************************************************************/
/* Private / Public Functions                                                    */
/*********************************************************************************/

INT FragmentClassifier::ClassifySequence(RAIphyDatabase *database, SequenceCounts *tempCounts, DOUBLE *matchScore)
{
    ULONG i,j;
    INT bestIndex = 0;
    DOUBLE bestScore = 0;
    DOUBLE secondBestScore = 0;
    DOUBLE tempScore = 0;

    for (i=0; i < tempCounts->nonzeroCount.size(); i++)
    {
        bestScore += (database->RAIProfiles.at(0).sequenceVector.at(tempCounts->nonzeroIndex.at(i)) * tempCounts->nonzeroCount.at(i));
    }

    for (i=1; i < database->RAIProfiles.size(); i++)
    {
        tempScore = 0;

        for (j=0; j < tempCounts->nonzeroCount.size(); j++)
        {
            tempScore += (database->RAIProfiles.at(i).sequenceVector.at(tempCounts->nonzeroIndex.at(j)) * tempCounts->nonzeroCount.at(j));
        }

        if (tempScore > bestScore)
        {
            secondBestScore = bestScore;
            bestScore = tempScore;
            bestIndex = i;
        }
    }

    *matchScore = fabs(bestScore - secondBestScore);

    return bestIndex;
}

/*********************************************************************************/

BOOLEAN FragmentClassifier::ClassifyFiles( RAIphyDatabase *database,
                                           vector<string> *inFilenames,
                                           string outFilename,
                                           INT classificationThreshold )
{
    Sequence tempSequence;
    SequenceCounts tempCounts;
    ULONG i, j, matchID;
    DOUBLE matchScore;
    DOUBLE *bestScores = NULL;
    string tempString;
    stringstream displayString;
    BYTE tempChar = 0; 
    fstream inFile;
    ofstream outFile;

    std::cout << "Binning DNA fragments..." << endl;

    if (classificationThreshold < 100)
    {
        bestScores = new DOUBLE[database->RAIProfiles.size()];
        for (i=0; i < database->RAIProfiles.size(); i++)
        { 
            bestScores[i] = 0;
        }
        outFile.open(string(outFilename+".temp").c_str());
    }
    else
    {
        outFile.open(outFilename.c_str()); 
    }

    if ((BOOLEAN)outFile.is_open())
    {
        outFile << HEADER_WORD_LENGTH << SEPARATOR << database->WordLength() << endl;
        for (i=0; i < inFilenames->size(); i++)
        {
            inFile.open(inFilenames->at(i).c_str(), ios::in); 
            j = 0;

            if (inFile.is_open())
            {
                do { tempChar = tolower(inFile.get()); } while (inFile.good() && tempChar != '>');
                if (!inFile.good())
                {
                    std::cout << "ERROR: No input sequences found in file! Program requires valid FASTA format." << endl;
                    outFile.close();
                    inFile.close();
                    return FALSE;
                }
                
                displayString.str("\r"); 
                
                while (inFile.good())
                {
                    displayString.str().replace(1, displayString.str().length(), " "); 
                    std::cout << displayString.str(); 
                    displayString.str(""); 
                    displayString << "\rProcessing sequence " << j+1 << " of file " << i+1 << "...";
                    std::cout << displayString.str(); 

                    m_profiler.ProfileNextSequence(&inFile, &tempSequence, database->WordLength(), &tempCounts);
                    matchID = ClassifySequence(database, &tempCounts, &matchScore);

                    if (classificationThreshold < 100)
                    {
                        outFile << ">" << tempSequence.sequenceName << endl
                                  << matchID << endl
                                  << matchScore << endl
                                  << database->RAIProfiles.at(matchID).sequenceName << endl;
                        if (matchScore > bestScores[matchID]) bestScores[matchID] = matchScore;
                    }
                    else
                        outFile << ">" << tempSequence.sequenceName << endl
                                  << database->RAIProfiles.at(matchID).sequenceName << endl;
                    j++;
                }
                inFile.close();
                std::cout << "Done!" << endl;
            }
        }

        std::cout << "Processing complete!" << endl; 
        outFile.close();

        if (classificationThreshold < 100)
        {
            CleanOutput(outFilename, bestScores, database, classificationThreshold);
            delete [] bestScores;
        }

        return TRUE;
    }
    else
    {
        std::cout << "ERROR: Could not create results/output file!" << endl;
        return FALSE;
    }
}

/*********************************************************************************/

void FragmentClassifier::CleanOutput(string outFilename, DOUBLE* bestScores, RAIphyDatabase* database, INT classificationThreshold)
{
    ULONG i, tempID;
    DOUBLE tempDouble;
    string tempString, headerLine;
    ofstream outFile; 
    ifstream inFile(string(outFilename+".temp").c_str()); 
    
    outFile.open(outFilename.c_str()); 

    if (outFile.is_open())
    {
        if (inFile)
        {
            tempDouble = (100 - (DOUBLE)classificationThreshold) / 100;
            for (i=0; i < database->RAIProfiles.size(); i++) 
            {
                bestScores[i] *= tempDouble;
            }

            getline(inFile, tempString); 
            if (!inFile.good())
            {
                std::cout << "ERROR: Could not reformat output." << endl;
                outFile.close();
                inFile.close();
                return;
            }

            outFile << tempString << endl;

            getline(inFile, headerLine);            // Sequence name
            while (inFile.good())
            {
                getline(inFile, tempString);        // Match Index
                tempID = atoi(tempString.c_str()); 
                getline(inFile, tempString);        // Match Score
                tempDouble = atof(tempString.c_str()); 
                getline(inFile, tempString);        // Match Name

                outFile << headerLine << endl;

                if (tempDouble >= bestScores[tempID])
                {
                    outFile << tempString << endl;
                }
                else
                {
                    outFile << NOT_CLASSIFIED << endl;
                }
                getline(inFile, headerLine);        // Sequence name
            }

            outFile.close();
            inFile.close();
            remove (string(outFilename+".temp").c_str());
        }
    }
}

/*********************************************************************************/

BOOLEAN FragmentClassifier::ClassifyFilesAndUpdateDatabase ( RAIphyDatabase *database,
                                                             vector<string> *inFilenames,
                                                             string outFilename,
                                                             INT* updatedCounts,
                                                             DOUBLE** updatedDatabase )
{
    ULONG i, j, k;
    INT matchID;
    DOUBLE matchScore;
    Sequence tempSequence;
    SequenceCounts tempCounts;
    fstream  inFile;
    ofstream outFile;
    string tempString;
    stringstream displayString; 
    BYTE tempChar; 
    
    outFile.open(outFilename.c_str()); 
    
    std::cout << "Binning DNA fragments..." << endl;

    if (outFile.is_open())
    {
        outFile << HEADER_WORD_LENGTH << SEPARATOR << database->WordLength() << endl;
        for (i=0; i < inFilenames->size(); i++)
        {
            j = 0;
            inFile.open(inFilenames->at(i).c_str(), ios::in); 

            if (inFile.is_open())
            {   
                do {  tempChar = tolower(inFile.get()); } while (inFile.good() && tempChar != '>');
                if (!inFile.good())
                {
                    std::cout << "ERROR: No input sequences found in " << inFilenames->at(i) 
                    << "!\nProgram requires valid FASTA format." << endl;
                    outFile.close();
                    inFile.close();
                    return FALSE;
                }

                displayString.str("\r"); 
                
                while (inFile.good())
                {
                    displayString.str().replace(1, displayString.str().length(), " "); 
                    std::cout << displayString.str(); 
                    displayString.str(""); 
                    displayString << "\rProcessing sequence " << j+1 << " of file " << i+1 << "...";
                    std::cout << displayString.str(); 

                    m_profiler.ProfileNextSequence(&inFile, &tempSequence, database->WordLength(), &tempCounts);
                    matchID = ClassifySequence(database, &tempCounts, &matchScore);
                    outFile << ">" << tempSequence.sequenceName << endl << database->RAIProfiles.at(matchID).sequenceName << endl;

                    updatedCounts[matchID] += tempSequence.sequenceLength;
                    for (k=0; k < tempCounts.nonzeroCount.size(); k++)
                    {
                        updatedDatabase[matchID][tempCounts.nonzeroIndex.at(k)] += tempCounts.nonzeroCount.at(k);
                    }
                    j++;
                }
                std::cout << "Done!" << endl;
                inFile.close();
            }

            std::cout << "Updating the database..." << endl;
            j = 0;
            displayString.str("\r"); 

            for ( i=0; i < database->RAIProfiles.size(); i++ )
            {
                if (updatedCounts[i] >= UPDATE_THRESHOLD)
                {
                    j++;
                    displayString.str().replace(1, displayString.str().length(), " "); 
                    std::cout << displayString.str(); 
                    displayString.str(""); 
                    displayString << "\rUpdating database item " << i+1 << "...";
                    std::cout << displayString.str(); 
                    database->UpdateVector(i, updatedDatabase[i]);
                }
            }
            std::cout << "Done!\nUpdated " << j << " item(s) in the database." << endl;
        }

        std::cout << "Processing complete!" << endl;
        outFile.close();
        return TRUE;
    }
    else
    {
        std::cout << "ERROR: Could not create results/output file!" << endl;
        return FALSE;
    }
}

/********************************* END OF FILE ***********************************/
