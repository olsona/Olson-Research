myMethods = ['raiphy']

usageString = 'bootstrap.py -i <input file> -o <output file> -p <RAIphy path> '\
            '[-r <reference database> -c <cooling schedule>]'
                
bigUsageString = '-i, --ifile\t\tInput metagenomics FASTA file\n'+\
                '-o, --ofile\t\tOutput file\n'+\
                '-r, --db, --reference\tReference database (optional)\n'+\
                '-c, --schedule\tCooling schedule (optional)\n'+\
                '\t\t\t  format: [1,2,5]\n'+\
                '\t\t\t  [1,2,5] implies the cooling schedule <= 1kbp, <= 2kbp,\n'+\
                '\t\t\t    <= 5kbp, > 5kbp\n'+\
                '-p, --path\t\tPath to RAIphy'
                
defaultSchedule = [1,2,5,10,20,50,100]