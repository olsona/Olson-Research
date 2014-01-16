myMethods = ['raiphy']

usageString = 'bootstrap.py -i <input file> -o <output file> -p <RAIphy path> '\
            '[-r <reference database> -c <cooling schedule>]'
                
bigUsageString = '-i, --ifile\t\tInput metagenomics FASTA file\n'+\
                '-o, --ofile\t\tOutput file\n'+\
                '-c, --schedule\tCooling schedule (optional)\n'+\
                '\t\t\t  format: [5,2,1]\n'+\
                '\t\t\t  [1,2,5] implies the cooling schedule >= 5kbp, >= 2kbp,\n'+\
                '\t\t\t    >= 1kbp\n'+\
                '-p, --path\t\tPath to RAIphy'
                
defaultSchedule = [100,50,20,10,5,2,1]