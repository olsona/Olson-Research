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
#include <sys/types.h>
#include <dirent.h>

/*********************************************************************************/
/* Private Module Constants */
/*********************************************************************************/
#define PARMS_DEFAULT_OUTFILE       "./outfile"
#define PARMS_DEFAULT_DATABASE      "./defaultDb"
#define PARMS_DEFAULT_EXTENSION     ".fasta"
#define PARMS_DEFAULT_THRESHOLD     100
#define PARMS_DEFAULT_MODE          MODE_BIN

/*********************************************************************************/
/* Private Type Definitions */
/*********************************************************************************/
/*********************************************************************************/
/* Private Module Function Prototypes */
/*********************************************************************************/
void PARMS_PrintHelpScreen(void); 
void PARMS_AddFilesFromDir(string dir);
string& trimleft(string& s);
string& trimright(string& s);
string& trim(string& s);

/*********************************************************************************/
/* Private Module Variables */
/*********************************************************************************/
static vector<string> m_InputFileNames;  
static string  m_OutputFileName;
static string  m_DatabaseFileName; 
static string  m_Extension; 
static INT     m_ClassificationThreshold; 
static INT     m_RunMode; 

/*********************************************************************************/
/* Public and Private Module Functions */
/*********************************************************************************/

void PARMS_Initialize(void)
{
    m_InputFileNames.clear();
    m_OutputFileName = PARMS_DEFAULT_OUTFILE;
    m_DatabaseFileName = PARMS_DEFAULT_DATABASE;
    m_Extension = PARMS_DEFAULT_EXTENSION; 
    m_ClassificationThreshold = PARMS_DEFAULT_THRESHOLD; 
    m_RunMode = PARMS_DEFAULT_MODE;
}

/*********************************************************************************/

BOOLEAN PARMS_ProcessCommandLine(INT argc, BYTE **argv)
{
    INT i; 
    BOOLEAN okToRun = TRUE;
    BOOLEAN useDirectory = FALSE; 
    string directory = "default";  

    for(i=1; i<argc; i++)
    {
        /* Input File ****************/ 
        if (strcmp(argv[i], "-i") == 0)
        {
            if (i+1 < argc)
            {
                m_InputFileNames.push_back(argv[i+1]); 
                i++; 
            }
            else
            {   
                std::cout << "No input file supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* Input Folder ****************/ 
        else if (strcmp(argv[i], "-I") == 0)
        {
            if (i+1 < argc)
            {
                directory = string(argv[i+1]); 
                useDirectory = TRUE; 
                i++; 
            }
            else
            {   
                std::cout << "No input folder supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* Output File ***************/ 
        else if (strcmp(argv[i], "-o") == 0)
        {
            if (i+1 < argc)
            {
                m_OutputFileName = string(argv[i+1]); 
                i++; 
            }
            else
            {
                std::cout << "No output file supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* Database File ***************/ 
        else if (strcmp(argv[i], "-d") == 0)
        {
            if (i+1 < argc)
            {
                m_DatabaseFileName = string(argv[i+1]); 
                i++; 
            }
            else
            {
                std::cout << "No database file supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* File Extension *************/ 
        else if (strcmp(argv[i], "-e") == 0)
        {
            if (i+1 < argc)
            {
                m_Extension = string(argv[i+1]); 
                i++; 
            }
            else
            {
                std::cout << "No file extension supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* Classification Threshold ***/ 
        else if (strcmp(argv[i], "-c") == 0)
        {
            if (i+1 < argc)
            {
                m_ClassificationThreshold = atoi(argv[i+1]); 
                i++; 
            }
            else
            {
                std::cout << "No threshold supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* Update Mode ***************/ 
        else if (strcmp(argv[i], "-m") == 0)
        {
            if (i+1 < argc)
            {
                m_RunMode = atoi(argv[i+1]); 
                i++; 
            }
            else
            {
                std::cout << "No threshold supplied" << endl;
                okToRun = FALSE; 
            }
        }
        /* Print Help ****************/ 
        else if (strcmp(argv[i], "-h") == 0)
        {
            PARMS_PrintHelpScreen();
            okToRun = FALSE; 
        }
        /* Unknown Parameter *********/ 
        else
        {
            std::cout << "Unknown Parameter: " << argv[i] << endl;
            okToRun = FALSE; 
        }
    }
    
    if (useDirectory)
    {
        PARMS_AddFilesFromDir(directory); 
    }

    if (okToRun && m_InputFileNames.size() < 1)
    {
        std::cout << "No input files provided." << endl;
        okToRun = FALSE; 
    }
    
    return okToRun; 
}

/*********************************************************************************/

void PARMS_AddFilesFromDir(string dir)
{
    struct dirent *dirElement=NULL;
    DIR *directory=NULL;
    
    directory=opendir(dir.c_str());
    if(directory == NULL)
    {
        std::cout << "Couldn't open input file directory" << endl; 
        exit(1); 
    }
        
    m_InputFileNames.clear(); 
        
    while((dirElement = readdir(directory)))
    {
        if (string(dirElement->d_name).find(m_Extension) != string::npos)
        {
            m_InputFileNames.push_back(dir + string(dirElement->d_name));
        }
    }
        
    closedir(directory);
}

/*********************************************************************************/

vector<string> *PARMS_GetInputFileNames(void)
{
    return &m_InputFileNames; 
}

/*********************************************************************************/

string PARMS_GetOutputFileName(void)
{
    return (m_OutputFileName); 
}

/*********************************************************************************/

string PARMS_GetDatabaseFileName(void)
{
    return (m_DatabaseFileName); 
}

/*********************************************************************************/

INT PARMS_GetClassificationThreshold(void)
{
    return m_ClassificationThreshold;
}

/*********************************************************************************/

INT PARMS_GetRunMode(void)
{
    return m_RunMode; 
} 

/*********************************************************************************/

void PARMS_PrintHelpScreen(void)
{
    std::cout << endl; 
    std::cout << "------------------------------------------------" << endl; 
    std::cout << "                RAIphy (v." << VERSION_MAJOR << "." << VERSION_MINOR << ")" << endl; 
    std::cout << "------------------------------------------------" << endl;
    std::cout << "-i : Specify input file." << endl; 
    std::cout << "-I : Specify input folder." << endl; 
    std::cout << "-o : Specify output file. (Default: \"" << PARMS_DEFAULT_OUTFILE << "\")" << endl; 
    std::cout << "-d : Specify database file. (Default: \"" << PARMS_DEFAULT_DATABASE << "\")" << endl;
    std::cout << "-e : Specify file extension. (Default: \"" << PARMS_DEFAULT_EXTENSION << "\")" << endl;
    std::cout << "-c : Specify classification threshold. (Default: \"" << PARMS_DEFAULT_THRESHOLD << "\")" << endl;
    std::cout << "-h : Print help screen." << endl;  
    std::cout << "-m : Specify run mode. (Default: \"" << PARMS_DEFAULT_MODE << "\")" << endl;
    std::cout << "\t" << MODE_BIN << " : Bin DNA Fragments." << endl; 
    std::cout << "\t" << MODE_BIN_WITH_REFINEMENT << " : Bin DNA Fragments with Iterative Refinement." << endl; 
    std::cout << "\t" << MODE_CREATE_DATABASE << " : Create Database." << endl; 
    std::cout << "\t" << MODE_ADD_TO_DATABASE << " : Add Items to Database." << endl; 
    std::cout << endl; 
}

/*********************************************************************************/

string& trimleft( string& s )
{
    string::iterator it;
    
    for( it = s.begin(); it != s.end(); it++ )
        if( !isspace( *it ) )
            break;
    
    s.erase( s.begin(), it );
    return s;
}

/*********************************************************************************/

string& trimright( string& s )
{
    string::difference_type dt;
    string::reverse_iterator it;
    
    for( it = s.rbegin(); it != s.rend(); it++ )
        if( !isspace( *it ) )
            break;
    
    dt = s.rend() - it;
    
    s.erase( s.begin() + dt, s.end() );
    return s;
}

/*********************************************************************************/

string& trim( string& s )
{
    trimleft( s );
    trimright( s );
    return s;
}

/*********************************************************************************/

string PARMS_Trim( const string& s )
{
    string t = s;
    return trim( t );
}

/*********************************************************************************/

void PARMS_GetStringParts(string originalString, vector<string> *stringParts, string delimiter)
{
    BYTE *cString, *token; 
    (*stringParts).clear();
    
    cString = new BYTE[originalString.size()+1]; 
    strcpy(cString, originalString.c_str()); 
    
    token = strtok(cString, SEPARATOR); 
    while (token != NULL)
    {
        (*stringParts).push_back(token); 
        token = strtok (NULL, delimiter.c_str());
    }
    
    delete[] cString; 
}

