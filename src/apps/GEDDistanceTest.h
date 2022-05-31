#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>
#include <chrono>
#include <numeric>
#include <unordered_map>
#include <array>
#include <map>

#include "Misc.h"
#include "Query.h"
#include "Metric.h"
#include "SearchCommon.h"
#include "SPMTree.h"
#include "CMTree.h"
#include "BruteForceSearch.h"

#include "TestCommon.h"
#define DIM 100
/*
	Various search algorithm tests with variations on the Hanming metric
*/

std::tuple<std::vector<GED>, std::vector<GED>> //return value
generatePointsQD(unsigned int nPoints, unsigned int nQueries) {
using namespace std;

	cout<<nQueries<<endl;
		srand(456789);
		using namespace std;
		
		std::vector<GED> graphs;
		graphs.reserve(nPoints); //For vector

		std::vector<GED> qGraphs;
		qGraphs.reserve(nPoints); //For vector
		ifstream myfile;
		myfile.open("../Graph_Edit_Distance/filelist.txt");
		string line;
		for(int i=0; i<nPoints; i++)
		{
			getline(myfile,line);
			// cout<<line<<endl;
			GED* word = new GED(line);
			graphs.push_back(*word);
			free(word);
		}

		for(int i=0; i<nQueries; i++)
		{
			qGraphs.push_back(graphs[rand()%nPoints]);
		}

		// cout<<qWords[0]<<endl;
		cout<<qGraphs.size()<<endl;
		return { graphs, qGraphs };
}

void kNNSearchCompare(unsigned int nPoints, unsigned int nQueries,
	PivotType pivT, PartType partT, unsigned int nofResults,
	const std::string& fileNamePrefix, bool csvh = true) {
	using namespace std;

	unsigned int bTime, sTime;
	using MetricType = GEDMetric;
	MetricType met;

	constexpr unsigned int dim = EuclidianPointDim;
	
	//auto [points, qPoints] = generatePointsQD(nPoints, nQueries);
	auto [points, qPoints] = generatePointsQD(nPoints, nQueries);
	//duzizhelicaishengchengxuyaozhongdi

	std::clock_t start = std::clock();
	SPMTree<GED, MetricType> stree(points, met, pivT, partT, kxBalancedTreeHeight(1, points.size()));
	CMTree<GED, MetricType> stree2(points, met, pivT, partT, kxBalancedTreeHeight(1, points.size()));
	// BruteForceSearch<GED, MetricType> stree2(points, met, pivT, partT);//, kxBalancedTreeHeight(1, points.size()));

	bTime = dTimeSeconds(start);
	cout << "firstkSearchTest btime=" << bTime << endl;

	start = std::clock();
	unsigned int diffCount = 0;
	unsigned int nFound = 0;
	unsigned int nqActual = 0;
	auto rad{ std::numeric_limits<float>::max() };
	for ( auto& qp : qPoints) {
		NearestKQuery<GED> rq(qp,  nofResults, rad);
		stree.search(rq);
		nqActual++;
		nFound += rq.getNeighbors().size();

		NearestKQuery<GED> rq2(qp,  nofResults, rad);
		stree2.search(rq2);

		if (!rq.hasSameNeighbors(rq2)) {
			diffCount++;
		}

		if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
	}
	sTime = dTimeSeconds(start);

	//AND the CSV file:
	ofstream csvFile(fileNamePrefix + ".csv", ios::app);
	// if (csvh == true) {
	// 	csvFile << ";;TEST NEARK_COMPARE_EM: " << currentDateTime()
	// 		<< ";;[CLasses]:" << endl
	// 		<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
	// 		<< "DiffCnt,name,MADI,dim,dbSize,nQueries,maxResults,nfound,aveNfound,Pivot,Partition,nodesVisited,numDistCalls," << endl;
	// }
	csvFile << diffCount << "," << stree.shortName() << "," << stree.getMADIorK() <<"," << dim << "," << points.size() << "," << nqActual << "," << nofResults << ","
		<< nFound << "," << ((nFound * 1.0f) / nqActual) << ","
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
		<< (1.0f * stree2.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
		<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
		<< static_cast<float>(stree2.getPerfStats().getDistanceCalls()) / nqActual <<endl;
		// << (1.0f * stree.getPerfStats().getDistanceCalls()) / (1.0f * nFound)  << ","
		// << (1.0f * stree2.getPerfStats().getDistanceCalls()) / (1.0f * nFound) <<  endl;
	csvFile.close();
}



void kNNSearchCompare(const std::string& fileNamePrefix) {
	//std::map<unsigned int, unsigned int> nofPoints{ {100,1},{1000,10}, {10000,10}, {1000000,100} };
	std::map<unsigned int, unsigned int> nofPoints{ {24440,100} };
	for (const auto& [np, nQueries] : nofPoints) {
		for (const auto& [pivType, pivVal] : pivotTypeMap) {
			for (const auto& [parType, parVal] : partTypeMap) {
				for (int mr=1;mr<1000;mr+=10) {
				// for (int mr=k;mr<k+n;mr=+5) {
					kNNSearchCompare(np, nQueries, pivType, parType, mr, fileNamePrefix );
				}
			}
		}
	}
}

