#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "EuclidianTests.h"
#include "Metric.h"
#include "SearchCommon.h"

#ifdef __WIN32
unsigned int DistanceIntervalSM<float>::used = 0;
std::vector<DistanceIntervalSM<float>> DistanceIntervalSM<float>::data;
#endif

int main(int argc, char* argv[]) {
    using namespace std;

    string fileNamePrefix = {"mtreeRunTimes_EM"};

    /*  EUCLIDIAN METRIC TESTS */

    // radiusSearchTestEM(fileNamePrefix);
    // radiusSearchTestEM_S(fileNamePrefix);

    // the count-only version of radius search
    // radiusSearchTestEMCount(fileNamePrefix);

    // radiusSearchTestEM_Brin95  (fileNamePrefix + "_Brin95");

    // nkSearchTestEM(fileNamePrefix + "_DIM10");
    // nkSearchTestEM(fileNamePrefix + "_DIM10_nk");
    // nkSearchTestEM(fileNamePrefix + "_LCMTs_DIM10");
    // nkSearchTestEM(fileNamePrefix + "_DIM10");

    //---------------------------------------------------------------------------------
    // radiusSearchCompareEM(fileNamePrefix);
    // radiusSearchCompareEM(10000, 4, PivotType::RAN, PartType::DMR, fileNamePrefix);
    //----------------------------------------------------------------------------------

    //---------------------------------------------
    // nkRsSearchTestEM(fileNamePrefix + "_DIM10_D");

    //--------------------------------------------------------
    // collectSearchTestEM(fileNamePrefix + "_pf");
    // collectSearchTestEMAutoRad(fileNamePrefix + "_DIM10b_zoom", 30, 0.3);
    // collectSearchTestEMAutoRad(fileNamePrefix + "_DIM10b", 10, 0.9);

    // collectCountSearchTestEM(fileNamePrefix);
    // collectCountSearchPlusTestEM(fileNamePrefix + "_2parts");

    //------------------------------------------------------
    kNNSearchCompareEM(100000, 10000, PivotType::RAN, PartType::BOM, 1, fileNamePrefix, true);

    // nkIncreasingDensityTestEM(100000, 10000, PivotType::RAN, PartType::PIV, 5, fileNamePrefix + "_id");
    // nkIncreasingDensityTestEM(fileNamePrefix + "_id");
    return 0;
}
