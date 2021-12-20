//#define SA_USE_STATIC_DI

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <chrono>
#include <regex>

#include "EditDistanceTests.h"

#ifdef __WIN32
unsigned int DistanceIntervalSM<float>::used = 0;
std::vector< DistanceIntervalSM<float>> DistanceIntervalSM<float>::data;
#endif

int main(int argc, char* argv[]) {
	using namespace std;

	string bfn {"/home/mzuniga/data"};
	string fileNamePrefix = { "mtreeRunTimes_EDM" };
	string cmFileName{ "checkMetrics.txt" };
	string sweFileName{ "SWExamples.txt" };
	string pairsFileName{ "pairsFileName.txt" };
	string sprotFileName{ bfn + "/UniProt/uniprot_sprot.fasta"};
	string tremblFileName{ bfn + "/UniProt/TREMBL/uniprot_trembl_500000.fasta" };
	string uproMammalsFileName{ bfn + "/UniProt/TREMBL/mammals/uniprot_trembl_mammals_500000.fasta" };
	string influenzaFileName{ bfn + "/HA/influenza.fna" };
    string haEncodingFileName{ bfn + "/HA/enc_09/ha_enc09.fasta"};
	string hivDirName{ bfn + "/HIV1" };
	string pridePrejFileName(bfn + "/text/words.txt");
	string dictionFileName(bfn + "/dictionary/dic_en_us_gb_au_lcunique.txt");
			
	/*   EDIT DISTANCE TESTS */

	//RADIUS (RANGE) SEARCH TESTS
	//radiusSearchTestEDM(fileNamePrefix , dictionFileName);
	//radiusSearchTestEDM_S(fileNamePrefix +"_sp_mo1_", sprotFileName);
	//radiusSearchTestEDM_S(fileNamePrefix + "_sprot_radius", sprotFileName);


	// KNN SEARCH TESTS

	//nkSearchTestEDM_EXQ(fileNamePrefix + "_spro_nk" , sprotFileName);
	//KNN+RANGE EQUIV SEARCH TESTS
	//nkRsSearchTestEDM(fileNamePrefix + "_SwissProt_" , sprotFileName);
	//nkSearchTestEDM_EXQ_NP(fileNamePrefix + "_spro_np_nk" , sprotFileName);

	//COMPARE AWNSERS FOR CORRECTNESS TESTS:
	//kNNCompareEDM(fileNamePrefix + "_knncmp", dictionFileName, 100000, 100, 10 , PivotType::RAN, PartType::BOM, true);

	//KNN search with initial radius
	nkSearchTestEDM_EXQ(fileNamePrefix + "_spro_rknn" , sprotFileName);

	//nkIncreasingDensityTestEM(100000, 10000, PivotType::RAN, PartType::PIV, 5, fileNamePrefix + "_id");
	//nkIncreasingDensityTestEM(fileNamePrefix + "_id");


	//collectSearchTestEDM_V2(fileNamePrefix + "_sprot" , sprotFileName, 0, 24, 50);
	//collectSearchTestEDM_V2(fileNamePrefix + "_sprot_zoom" , sprotFileName, 0, 10, 10);
	

	//collectCountSearchTestEDM(fileNamePrefix + "_Sprot_count_10K" , sprotFileName);


	//treeStatsTests("treeStrucStats_kimmo_rs" , dictionFileName);
 
	return 0;

}





