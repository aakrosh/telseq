/*
 * telseq.h
 *
 *  Created on: 6 Sep 2013
 *      Author: zd1
 */

#ifndef TELSEQ_H_
#define TELSEQ_H_


#include <getopt.h>
#include <vector>
#include <map>
#include <sstream>
#include <cstdint>
#include "config.h"
#include "Util.h"


//
// typedef
//
typedef std::vector<std::string> StringVector;

namespace ScanParameters{

	static std::string FIELD_SEP="\t";
	// this is in basepair
	const unsigned int GENOME_LENGTH_AT_TEL_GC =  332720800;
	// output unit, it kb by default (1000), change it to 1 if in bp
	const unsigned int LENGTH_UNIT = 1000;
	static unsigned int READ_LENGTH = 100;
	static std::string PATTERN="TTAGGG";
	static std::string PATTERN_REVCOMP="CCCTAA";
	const uint64_t TELOMERE_ENDS = 46;
	const float GC_LOWERBOUND = 0.4;
	const float GC_UPPERBOUND = 0.6;
	const float GC_BINSIZE = 0.02;

	const float GC_TELOMERIC_LOWERBOUND = 0.48;
	const float GC_TELOMERIC_UPPERBOUND = 0.52;

	// maximum motif counts. add 1 to include 0 count.
	static uint64_t TEL_MOTIF_N = READ_LENGTH/PATTERN.size() +1;
	const uint64_t TEL_MOTIF_CUTOFF = 7;
	const int GC_BIN_N = (int) ((GC_UPPERBOUND-GC_LOWERBOUND)/GC_BINSIZE+0.5);

	const std::string LABEL_RG="ReadGroup";
	const std::string LABEL_LB="Library";
	const std::string LABEL_SAMPLE="Sample";
	const std::string LABEL_BAM="BAM";
	const std::string LABEL_TOTAL="Total";
	const std::string LABEL_MAPPED="Mapped";
	const std::string LABEL_DUP="Duplicates";
	const std::string LABEL_TEL="TEL";
	const std::string LABEL_GC="GC";
	const std::string LABEL_LEN="LENGTH_ESTIMATE";

	const std::string SCAN_FILE_SUFFIX = "bamscan";

};

struct ScanResults
{
	std::string sample;
	std::string lib;
	std::string bam;
	std::vector<uint64_t> telcounts;
	std::vector<uint64_t> gccounts;
	uint64_t numTotal;
	uint64_t numMapped;
	uint64_t numDuplicates;
	double telLenEstimate;

	// logging data
	uint64_t n_totalunfiltered = 0;
	uint64_t n_exreadsExcluded = 0;
	uint64_t n_exreadsChrUnmatched=0;

	ScanResults() { setDefaults(); }

    // Set reasonable default values for the qc filters
    void setDefaults()
    {
		telcounts = std::vector<uint64_t>(ScanParameters::TEL_MOTIF_N,0);
		gccounts = std::vector<uint64_t>(ScanParameters::GC_BIN_N,0);
        numTotal = 0;
        numMapped = 0;
        numDuplicates = 0;
        telLenEstimate = 0;
        n_exreadsExcluded = 0;
        n_exreadsChrUnmatched=0;
        n_totalunfiltered = 0;
    }


};

// headers for the output

struct Headers{

	Headers(){setHeaders();}
	void setHeaders(){
		headers.push_back(ScanParameters::LABEL_RG);
		headers.push_back(ScanParameters::LABEL_LB);
		headers.push_back(ScanParameters::LABEL_SAMPLE);
		headers.push_back(ScanParameters::LABEL_TOTAL);
		headers.push_back(ScanParameters::LABEL_MAPPED);
		headers.push_back(ScanParameters::LABEL_DUP);
		headers.push_back(ScanParameters::LABEL_LEN);

		for(unsigned int i=0;i < ScanParameters::TEL_MOTIF_N;i++){
			std::string h = ScanParameters::LABEL_TEL + NumberToString(i);
			headers.push_back(h);
		}
		for(int i=0;i<ScanParameters::GC_BIN_N;i++){
			std::string h = ScanParameters::LABEL_GC + NumberToString(i);
			headers.push_back(h);
		}
	}
	std::vector<std::string> headers;

};


void parseScanOptions(int argc, char** argv);
double calcGC(const std::string& seq);
int countMotif(std::string &read, std::string pattern, std::string pattern_revcomp);
double calcTelLength(ScanResults results);
int scanBam();
int outputresults(std::vector< std::map<std::string, ScanResults> > );
void printlog(std::vector< std::map<std::string, ScanResults> > resultlist);
void printout(std::string, ScanResults, std::ostream*);


//std::ifstream::pos_type getFilesize(const std::string& filename);
//std::istream* createReader(const std::string& filename);
//std::ostream* createWriter(const std::string& filename);

#endif /* TELSEQ_H_ */
