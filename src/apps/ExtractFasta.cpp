//#define SA_USE_PARASAIL
//#define SA_USE_STATIC_DI

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <chrono>

#include "Sequence.h"
#include "Misc.h"
#include "Query.h"
#include "Reader.h"


int main(int argc, char* argv[]) {
	using namespace std;
	if (argc != 2) {
		cout << "Usage:<max nof sequences>" << endl;
		std::this_thread::sleep_for(std::chrono::seconds(4));
		exit(-1);
	}
	unsigned int maxSeqs = 500000;
	//stringstream sstream(argv[1]);
	//sstream >> maxSeqs;

	string seqRecsFileName{ "C:/data/UniProt/uniprot_sprot.fasta" };
	string seqRecsFileName2{ "C:/data/UniProt/TREMBL/uniprot_trembl.fasta" };
	string seqRecsFileName3{ "C:/data/UniProt/TREMBL/mammals/uniprot_trembl_mammals.fasta" };

	std::clock_t start = std::clock();
	FastaReader reader;
	reader.sampleSubsetAndCopy(seqRecsFileName3, maxSeqs, 2000000);
	//reader.sampleSubsetAndCopy(seqRecsFileName2, 5, 20);
	auto readTime = dTimeSeconds(start);
	cout << "sampleTime=" << readTime <<endl ;

	return 0;
}


