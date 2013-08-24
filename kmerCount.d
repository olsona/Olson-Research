import std.stdio;
import std.file;
import std.string;
import std.algorithm;
import std.conv;
import std.math;
import std.string;

void main(string[] args) {

    int i;
    char[] h;
    char[] b;
    char[] s;
    uint[string] freqs;
    uint[int] histo;
    int klen;
    int arg;
    int lineCount;
	int numEvents;

	if (args.length > 1) {
		if (isNumeric(args[1])) arg = to!int(args[1]);
  //      writefln("args[1] =%s", args[1]);
  //      writefln("arg = %d", arg);
		if (arg > 0 && arg < 40) {
			klen = arg;
      	}
      	else {
			writefln("bad kmer length -- must be between 1 and 40");
			return;
		} 
	}
    else {
    	klen = 29;
    }
    
    char[char] rules;
    rules['a'] = 't';
    rules['c'] = 'g';
    rules['g'] = 'c';
    rules['t'] = 'a';
    rules.rehash;
    
    char[] invalid = ['N','n'];
    
    foreach (line; stdin.byLine()) {
    	if (line.length > 0 && line[0] == '>') {
			if (b.length != 0 ) {
	  			lineCount++;
	  			if (lineCount % 100000) {
	    			writefln("Read %s", lineCount);
	  			}
	  			for (i = 0; b.length > i + klen; i++) {
	    			s = b[i..i+klen];
	    			if (isCorrect(s, invalid) != 0) {
      					++freqs[s.idup];
      					char[] rc = reverseComplement(s.idup, rules);
	    				++freqs[rc.idup];
	    				numEvents += 2;
	    			}
	    			else
	    				i += klen;
	  			}
	  			b.length = 0;
			}
			h = line;
      	}
      	else {
			b = b ~ line;
      	}
    }
    
    for (i = 0; b.length > i + klen; i++) {
		s = b[i..i+klen];
		if (isCorrect(s, invalid) != 0) {
      		++freqs[s.idup];
      		char[] rc = reverseComplement(s.idup, rules);
	    	++freqs[rc.idup];
	    	numEvents += 2;
	    }
    }
    
    writefln("Total:\t%u", numEvents);
    
    string[] kmers = freqs.keys;		    
    long nkmer = 0;
    int k;
//	sort!((a,b) {return freqs[a] > freqs[b]; })(kmers);
	sort!((a,b) {return a < b;})(kmers);
    
    double D_kl = 0.0;
    double q = 1.0/(cast(real)pow(4,klen));
    
    foreach(kmer; kmers) {
		if (freqs[kmer] > 0) {
			double pi = cast(float)freqs[kmer]/cast(float)numEvents;
			D_kl += log(pi/q)*pi;
		}
//		k = freqs[kmer];
//		++histo[k];
		nkmer++;
    } 
    writefln("Kullback-Leibler divergence:\t%.12f",D_kl);
    
	writefln("Fraction of kmers that are represented:\t%.16f\n", cast(double)nkmer/cast(double)pow(4,klen));
    
    foreach(kmer; kmers) {
    	writefln("%s\t%u", kmer, freqs[kmer]);
    }
    /+
    int[] histos = histo.keys;
    
    sort!((a,b) { return a < b; })(histos);

    writefln("Count\tNumber of Occurances\n");
    foreach(his; histos) {
    	writefln("%u\t%u", his, histo[his]);
    }
    +/
}

int isCorrect(char[] input, char[] invalid) {
	int corr = 1;
	foreach (char c; input) {
		foreach (char i; invalid) {
			if (i == c) {
				corr = 0;
				break;
			}
		}
	}
	return corr;
}

int isCharCorrect(char input, char[] invalid) {
	int corr = 1;
	foreach (char i; invalid) {
		if (i == input) {
			corr = 0;
			break;
		}
	}
	return corr;
}

char[] reverseComplement(string str, char[char] rules) {
	ulong l = str.length;
	char[] res = new char[l];
	int i;
	char cres;
	for (i = 0; i < l; i++) {
		cres = rules.get(str[i], 'N');
		res[l-i-1] = cres;
	}
	return res;
} 