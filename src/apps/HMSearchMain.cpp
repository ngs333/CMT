#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <chrono>
#include "HMDistanceTest.h"

int main(int argc, char* argv[]) {
	using namespace std;
	int start, end; // for parrel
	string filename;
	if(argc!=3)
	{
		cout<<"error!"<<endl;
	}
	else{
		start = stoi(argv[1]);
		end = stoi(argv[2]);
	}
	string fileNamePrefix =  "HM_distance_k="+string(argv[1]) + "_n=" + string(argv[2]);
	radiusSearchCompareEM(fileNamePrefix start, end);//zhidaoyylenianbizimeixiestrigkunbizi

	/*  EUCLIDIAN METRIC TESTS */
	
	//radiusSearchTestEM(fileNamePrefix);
	//radiusSearchTestEM_S(fileNamePrefix);

	//the count-only version of radius search
	//radiusSearchTestEMCount(fileNamePrefix);


	//radiusSearchTestEM_Brin95  (fileNamePrefix + "_Brin95");

	//nkSearchTestEM(fileNamePrefix + "_DIM10");
	//nkSearchTestEM(fileNamePrefix + "_DIM10_nk");
	//nkSearchTestEM(fileNamePrefix + "_LCMTs_DIM10");
	//nkSearchTestEM_S(fileNamePrefix + "_DIMSS_IQI_V2");
	
	//---------------------------------------------------------------------------------
	//radiusSearchCompareEM(fileNamePrefix);
	//radiusSearchCompareEM(10000, 4, PivotType::RAN, PartType::DMR, fileNamePrefix);
	//----------------------------------------------------------------------------------
	
	//---------------------------------------------
	//nkRsSearchTestEM(fileNamePrefix + "_DIM10_D");

	//--------------------------------------------------------
	//collectSearchTestEM(fileNamePrefix + "_pf");
	//collectSearchTestEMAutoRad(fileNamePrefix + "_DIM10b_zoom", 30, 0.3);
	// collectSearchTestEMAutoRad(fileNamePrefix + "_DIM10b", 10, 0.9);
	
	//collectCountSearchTestEM(fileNamePrefix);
	//collectCountSearchPlusTestEM(fileNamePrefix + "_2parts");

	
	//------------------------------------------------------
	// kNNSearchCompareEM(100000, 1000, PivotType::RAN, PartType::DMR, 1, fileNamePrefix, true);
	// kNNSearchCompare(fileNamePrefix,k,n);

	//nkIncreasingDensityTestEM(100000, 10000, PivotType::RAN, PartType::PIV, 5, fileNamePrefix + "_id");
	//nkIncreasingDensityTestEM(fileNamePrefix + "_id");
	return 0;
}






