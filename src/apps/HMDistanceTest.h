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

std::tuple<std::vector<HMPoint>, std::vector<HMPoint>> //return value
generatePointsQD(unsigned int nPoints, unsigned int nQueries) {
using namespace std;

	cout<<nQueries<<endl;
		srand(456789);
		using namespace std;
		constexpr unsigned int dim = 100; //string length
		
		std::vector<HMPoint> words;
		words.reserve(nPoints); //For vector

		std::vector<HMPoint> qWords;
		qWords.reserve(nPoints); //For vector
		int nQCount = 0;
		// 65-90

		char *myStr = (char*)malloc(dim * sizeof(char)); // For somereason, use array will cause additional char at the end
		for(int j=0; j<dim; j++)
		{
			myStr[j] = rand()%(2)+65;
			// cout<<myStr[j]<<endl;
		}
		string s = myStr;
		HMPoint* word = new HMPoint(s, s);
		words.push_back(*word);
		free(word);
		for(int i=1; i<nPoints; i++)
		{
			// else
			// {
				for(int j=0; j<dim; j++)
				{
					myStr[j] = rand()%2+65;
				}
				string s = myStr;
			// cout<<s<<endl;

				// string s = "aaaaaa";
				// cout<<s<<endl;
				HMPoint* word = new HMPoint(s, s);
				words.push_back(*word);
				free(word);
			// }

		}

		// int notInSet = nPoints; //rand()%(nPoints/50); // not 10 ()  not 10* xkoukun
		int notInSet = 0;rand()%(nQueries/50); // not 10 ()  not 10* xkoukun
		cout<<"-----------"<<endl;
		// qWords.push_back(words[rand()%nPoints]);
		for(int i=0; i<nQueries-notInSet; i++)
		{
			qWords.push_back(words[rand()%nPoints]);
		}
		cout<<"notInSet: "<<notInSet<<endl;
		for(int i=0; i<notInSet; i++)
		{
			for(int j=0; j<dim; j++)
			{
				myStr[j] = rand()%(2)+65;
			}
			string s = myStr;
			// string s = "aaaaaa";
			// cout<<s<<endl;
			HMPoint* word = new HMPoint(s, s);
			qWords.push_back(*word);
			free(word);
		}
		free(myStr);
		// cout<<qWords[0]<<endl;
		cout<<qWords.size()<<endl;
		return { words, qWords };
}

void kNNSearchCompare(unsigned int nPoints, unsigned int nQueries,
	PivotType pivT, PartType partT, unsigned int nofResults,
	const std::string& fileNamePrefix, bool csvh = true) {
	using namespace std;

	unsigned int bTime, sTime;
	using MetricType = HMMetric;
	MetricType met;

	constexpr unsigned int dim = EuclidianPointDim;
	
	//auto [points, qPoints] = generatePointsQD(nPoints, nQueries);
	auto [points, qPoints] = generatePointsQD(nPoints, nQueries);
	//duzizhelicaishengchengxuyaozhongdi

	std::clock_t start = std::clock();
	CMTree<HMPoint, MetricType> stree(points, met, pivT, partT, kxBalancedTreeHeight(1, points.size()));
	SPMTree<HMPoint, MetricType> stree2(points, met, pivT, partT, kxBalancedTreeHeight(1, points.size()));

	bTime = dTimeSeconds(start);
	cout << "firstkSearchTest btime=" << bTime << endl;

	start = std::clock();
	unsigned int diffCount = 0;
	unsigned int nFound = 0;
	unsigned int nqActual = 0;
	auto rad{ std::numeric_limits<float>::max() };
	for ( auto& qp : qPoints) {
		NearestKQuery<HMPoint> rq(qp,  nofResults, rad);
		stree.search(rq);
		nqActual++;
		nFound += rq.getNeighbors().size();

		NearestKQuery<HMPoint> rq2(qp,  nofResults, rad);
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
		<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
		<< static_cast<float>(stree2.getPerfStats().getDistanceCalls()) / nqActual <<endl;
		// << (1.0f * stree.getPerfStats().getDistanceCalls()) / (1.0f * nFound)  << ","
		// << (1.0f * stree2.getPerfStats().getDistanceCalls()) / (1.0f * nFound) <<  endl;
	csvFile.close();
}



void kNNSearchCompare(const std::string& fileNamePrefix, int k, int n) {
	//std::map<unsigned int, unsigned int> nofPoints{ {100,1},{1000,10}, {10000,10}, {1000000,100} };
	std::map<unsigned int, unsigned int> nofPoints{ {100000,1000} };
	// std::vector<unsigned int> maxResults{ 1,2,3,4,5,6,7,8,9,10  };
	//
	for (const auto& [np, nQueries] : nofPoints) {
		for (const auto& [pivType, pivVal] : pivotTypeMap) {
			for (const auto& [parType, parVal] : partTypeMap) {
				for (int mr=k;mr<k+n;mr*=2) {
				// for (int mr=k;mr<k+n;mr=+5) {
					kNNSearchCompare(np, nQueries, pivType, parType, mr, fileNamePrefix );
				}
			}
		}
	}
}
