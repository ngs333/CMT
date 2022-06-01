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
		int notInSet = 0;//rand()%(nQueries/50); // not 10 ()  not 10* xkoukun
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
		<< (1.0f * stree2.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
		<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
		<< static_cast<float>(stree2.getPerfStats().getDistanceCalls()) / nqActual <<endl;
		// << (1.0f * stree.getPerfStats().getDistanceCalls()) / (1.0f * nFound)  << ","
		// << (1.0f * stree2.getPerfStats().getDistanceCalls()) / (1.0f * nFound) <<  endl;
	csvFile.close();
}



void kNNSearchCompare(const std::string& fileNamePrefix, int k, int n) {
	//std::map<unsigned int, unsigned int> nofPoints{ {100,1},{1000,10}, {10000,10}, {1000000,100} };
	std::map<unsigned int, unsigned int> nofPoints{ {1000,100} };
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



void radiusSearchCompareEM(unsigned int nPoints, const unsigned int nQueries, PivotType pivT, PartType partT, 
	const std::string& fileNamePrefix, float rad,  std::vector<HMPoint> points, std::vector<HMPoint> qPoints) {
	using namespace std;
	unsigned int bTime, sTime;
	using MetricType = HMMetric;
	MetricType met;
	unsigned int nqActual = 0;

	// auto [points, qPoints] = generatePointsQD(nPoints, nQueries);

	auto start = std::clock();
	CMTree<HMPoint, MetricType> stree(points, met,pivT, partT, kxBalancedTreeHeight(1,points.size()));
	// BruteForceSearch<HMPoint, MetricType> stree2(points, met); //This is not baseline weism huilaile xiamian de zhixian shi stree not 2
	SPMTree<HMPoint, MetricType> stree2(points, met, pivT, partT, kxBalancedTreeHeight(1,points.size()));  //Is that because I didn't put in enough arguments?
	bTime = dTimeSeconds(start);
	cout << "radiusSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	start = std::clock();
	unsigned int nFound = 0;
	unsigned int diffCount = 0;
	const unsigned int maxResults = 4294967295;
	// auto rad = 1;
	for (const auto& qp : qPoints) {
		RadiusQuery<HMPoint> rq(qp, rad, maxResults);
		stree.search(rq);

		// NearestKQuery<HMPoint> rq2(qp, rad, maxResults);
		RadiusQuery<HMPoint> rq2(qp, rad, maxResults);
		stree2.search(rq2); //guaibude0zxkousuankunhuangkun

		nFound += rq.getNeighbors().size();


	    nqActual++;
		if ((nqActual % 10000) == 0) { cout << "finished search i= " << nQueries << endl; }
	}
	sTime = dTimeSeconds(start);

	string runTimeFileName{ fileNamePrefix + "_rs_compare.csv"};
	ofstream theFileA(runTimeFileName, ios::app);
	theFileA 
		<< rad << "," << points.size() << "," << diffCount << "," << nqActual << "," << nFound << "," << radiusSum << ","
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (stree.getPerfStats().getNodesVisited() / nQueries) << ","
		<< (stree.getPerfStats().getDistanceCalls() / nQueries) << ","
		<< (stree2.getPerfStats().getDistanceCalls() / nQueries)<<endl;
		// << bTime << "," << sTime << endl;

	theFileA.close();	
}


void radiusSearchCompareEM(const std::string& fileNamePrefix, float start, float end) {
	std::map<unsigned int, unsigned int> nofPoints{  {10000000, 20} };
	// std::vector<float> rads{2,4,8,16,32,64,128,256,512,1024,2048,4096,8192}; 
	for (const auto& [np, nSkip] : nofPoints) {
		auto [points, qPoints] = generatePointsQD(np, nSkip);

		for (const auto& [pivType, pivVal] : pivotTypeMap) {
			for (const auto& [parType, parVal] : partTypeMap) {
				for (float rad = start; rad<end; rad++) {
					radiusSearchCompareEM(np, nSkip, pivType, parType, fileNamePrefix, rad, points, qPoints);
				}
			}
		}
	}
}
