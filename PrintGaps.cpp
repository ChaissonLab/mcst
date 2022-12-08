#include "htslib/sam.h"
#include "htslib/hts.h"
#include <iostream>
#include <string>
#include <stdlib.h>

#include <fstream>
#include "htslib/kseq.h"


#include "htslib/faidx.h"
using namespace std;


int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "Usage: pg input.bam ref.fasta [--region region] [ --out outFile ] [--minsig n(40)]" << endl;
		exit(1);
  }
	string filename = argv[1];
	string refFileName=argv[2];
	string outFileName="";
	bool printSeq=false;
	string iterChrom="";
	int minSig=40;
	for (int argi=3; argi < argc; argi++) {
		if (argv[argi] == "--seq") {
			printSeq = true;
		}
		if (string(argv[argi]) == "--region") {
			++argi;
			iterChrom=argv[argi];
			//			cerr << "iterating " << iterChrom << endl;
		}			
		if (string(argv[argi]) == "--out") {
			++argi;
			outFileName = argv[argi];
		}
		if (string(argv[argi]) == "--minsig") {
		  ++argi;
		  minSig=atoi(argv[argi]);
		}
		++argi;
	}
	ofstream outFile;
	ostream* outPtr;
	if (outFileName == "") {
		outPtr = &cout;
	}
	else {
		outFile.open(outFileName.c_str());
		outPtr = &outFile;
	}

	htsFile *htsfp;
	bam_hdr_t *samHeader;			
  char *refSeq=NULL;
  faidx_t *ref;
  ref=fai_load(refFileName.c_str());
  htsfp = hts_open(filename.c_str(), "r");
	const htsFormat *fmt = hts_get_format(htsfp);

  hts_idx_t *bamIndex=sam_index_load(htsfp, filename.c_str());

	samHeader = sam_hdr_read(htsfp);
	hts_itr_t *itr;
	if (iterChrom != "") {
		itr=sam_itr_querys(bamIndex, samHeader, iterChrom.c_str());
	}

	bam1_t *b = bam_init1();
	int res=1;
	if (iterChrom != "") {
		res= sam_itr_next(htsfp, itr, b);
	}
	else {
		res = sam_read1(htsfp, samHeader, b);
	}
	long readIndex=0;
	long nCounted=0;
	string chrom="";
	const char* ops[] = {"insertion", "deletion" };

	while (res > 0) {
		readIndex+=1;
		long refPos = b->core.pos;
		int tid=b->core.tid;
		if (tid < 0) {
			break;
		}

		string alignedChrom(sam_hdr_tid2name(samHeader, tid));
		if (refSeq == NULL || chrom != alignedChrom) {
			if (refSeq != NULL) {
				free(refSeq);
			}
			int len;
			refSeq = fai_fetch(ref, alignedChrom.c_str(), &len);
			cerr << "pg:\t" << alignedChrom << "\t" << strlen(refSeq)  << endl;
			chrom=alignedChrom;
		}
		int cigarLen = b->core.n_cigar;
		long qPos=0;
		
		unsigned int *cigar=bam_get_cigar(b);
		string query;
		query.resize(b->core.l_qseq);
		if (query.size() == 0) {
			res = sam_read1(htsfp, samHeader, b);
			continue;
		}
		uint8_t *q = bam_get_seq(b);
		for (int qi=0; qi < query.size(); qi++) {query[qi]=seq_nt16_str[bam_seqi(q,qi)];	}
		for (int i = 0; i < cigarLen; i++) {
			char op = bam_cigar_opchr(cigar[i]);
			int  opLen = bam_cigar_oplen(cigar[i]);
			int opi=0;
			if (op == 'I' or op == 'D') {
				int rEnd;
				int qEnd;
				if (op == 'I') {
					rEnd=refPos+1;
					qEnd=qPos + opLen;
					opi=0;
				}
				else {
					rEnd = refPos + opLen;
					qEnd = qPos+1;
					opi=1;
				}
				//				cerr << i << "\t" << op << "\t" << opLen << endl;
				if (opLen>= minSig && (op == 'I' || op == 'D')) {					
					
					string qSeq(&query[qPos], qEnd-qPos);
					string rSeq(&refSeq[refPos], rEnd-refPos);
					(*outPtr) << chrom << "\t" << refPos << "\t" << rEnd << "\t" << ops[opi] << "\t" 
							 << opLen << "\t" << bam_get_qname(b) << "\t" << qPos << "\t" << qEnd ;
					if ( printSeq ) {
						(*outPtr) << "\t" << rSeq << "\t" << qSeq;
					}
					(*outPtr) << endl;
				}
			}
			if (op == 'S') {
				qPos += opLen;
			}
			else if (op == 'M' || op == 'X' || op == '=') {
				qPos += opLen;
				refPos  += opLen;
			}
			else if (op == 'D') {
				refPos += opLen;
			}
			else if (op == 'I') {
				qPos += opLen;
			}
		}
		if (iterChrom != "") {
			res= sam_itr_next(htsfp, itr, b);
		}
		else {
			res = sam_read1(htsfp, samHeader, b);
		}

		if (readIndex %10000 == 0) {
			cerr << "pg: " << readIndex / 1000 << "k" << endl;
		}
	}
}
