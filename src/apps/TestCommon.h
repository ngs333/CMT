#ifndef TEST_COMMOM
#define TEST_COMMON

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
#include "TreeNodes.h"
#include "CMTree.h"
#include "BruteForceSearch.h"


template <class ObjType, class MetricType, class TreeType>
void radiusSearchTest(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects, 
    unsigned int avgDim, const std::vector<float>& radii, 
	const PivotType pivT,  const PartType partT, float madi, 
	const std::string& fileNamePrefix, bool csvh) {
	using namespace std;
	unsigned int bTime, sTime;

    MetricType met;

    std::clock_t start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "radiusSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	unsigned int maxResults = 1000000;
	for (auto& radius : radii) {
		cout << "starting radiusSearchTest radius=" << radius << endl;
		stree.getPerfStats().reset();
		start = std::clock();
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		for (auto& qp : qObjects) {
			RadiusQuery<ObjType> rq(qp, radius, maxResults);
			stree.search(rq);
			nqActual++;
			nFound += rq.getNeighbors().size();
			
			if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
			rq.getNeighbors().clear();
		}
		sTime = dTimeSeconds(start);

		ofstream csvFile(fileNamePrefix + "_rs.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST RADIUS: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim,dbSize,collection, nQueries,radius,nfound,aveNfound,radiusSum,Pivot,Partition,NNodesVistd,NDistCalls,btime,stime,stime_n1000" << endl;
			csvh = false;
		}
		
		csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << "," 
			<<"F" << "," << nqActual << "," << radius << ","
			<< nFound << "," << ((nFound * 1.0f) / nqActual) << "," << radiusSum << ","
			<< stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();
	}
}

//In this version the results are not stored, only a count of the objects meeting the search radius.
//Used for some large radius complexity testing.
template <class ObjType, class MetricType, class TreeType>
void radiusSearchTestCount(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects, 
    unsigned int avgDim, const std::vector<float>& radii, 
	const PivotType pivT,  const PartType partT, float madi, 
	const std::string& fileNamePrefix, bool csvh) {
	using namespace std;
	unsigned int bTime, sTime;

    MetricType met;

    std::clock_t start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "radiusSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	unsigned int maxResults = 1000000;
	for (auto& radius : radii) {
		cout << "starting radiusSearchTest radius=" << radius << endl;
		stree.getPerfStats().reset();
		start = std::clock();
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		for (auto& qp : qObjects) {
			RadiusCountQuery<ObjType> rq(qp, radius);
			stree.search(rq);
			nqActual++;
			nFound += rq.getCount();
			
			if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
		}
		sTime = dTimeSeconds(start);

		ofstream csvFile(fileNamePrefix + "_rs_count.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST RADIUS: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim,dbSize,nQueries,radius,nfound,aveNfound,radiusSum,Pivot,Partition,NNodesVistd,NDistCalls,btime,stime,stime_n1000" << endl;
			csvh = false;
		}
		
		csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << ","
			<< nqActual << "," << radius << ","
			<< nFound << "," << ((nFound * 1.0f) / nqActual) << "," << radiusSum << ","
			<< stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();
	}
}



/* Collection count search - the query does not store the IDs or references to the 
individial results - it just stores a count of the number of results.
*/
template <class ObjType, class MetricType, class TreeType>
void collectCountSearchTest(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects, 
    unsigned int avgDim, const std::vector<float>& radii, bool radiusMult,
	const PivotType pivT,  const PartType partT, float madi, 
	const std::string& fileNamePrefix, bool csvh ) {
	using namespace std;
	unsigned int bTime, sTime;

	MetricType met;
    
    std::clock_t start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "radiusSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	//unsigned int maxResults = 100000000;
	for (auto& radiusM : radii) {
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		cout << "starting collectSearchTest collect component radius=" << radiusM << endl;
		stree.getPerfStats().reset();
		start = std::clock();
		nFound = 0;
		nqActual = 0;
		for (auto& qp : qObjects) {
			auto searchRadius = radiusM;
			if(radiusMult == true){
				searchRadius = qp.getValue().size() * radiusM;
			}
			RadiusCountQuery<ObjType> rq(qp, searchRadius);
			stree.searchCollect(rq);
			nqActual++;
			nFound += rq.getCount();
			if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
		
		}
		sTime = dTimeSeconds(start);
		ofstream csvFile(fileNamePrefix + "_count_cs.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST COLLECTION COUNT: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim,dbSize,collection, nQueries,radius,nfound,aveNfound,radiusSum,Pivot,Partition,NNodesVistd,NDistCalls,btime,stime,stime_n1000" << endl;
			csvh = false;
		}
		
		csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << "," <<"T" << ","
			<< nqActual << "," << radiusM << ","
			<< nFound << "," << ((nFound * 1.0f) / nqActual) << "," << radiusSum << ","
			<< stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();
	}
}

/**
  Test of the standard collection tests where the results are stored by the query.
  If input flag collectPlusRange is true, then both collection and the regualr range/radius
   search are performed.
**/
template <class ObjType, class MetricType, class TreeType>
void collectSearchTest(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects, 
    unsigned int avgDim, const std::vector<float>& radii, bool radiusMult,
	const PivotType pivT,  const PartType partT, float madi, 
	const std::string& fileNamePrefix, bool csvh, bool collectPlusRange = false ) {
	using namespace std;
	unsigned int bTime, sTime;

	MetricType met;
    
    std::clock_t start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "radiusSearchTest btime=" << bTime << endl;

    if (csvh == true) {
      ofstream csvFile(fileNamePrefix + "_crs.csv", ios::app);
      csvFile << ";;TEST COLLECTION: " << currentDateTime() << ";;[CLasses]:" << endl
              << ";;[" << typeid(stree).name() << "],"
              << "[" << met << "]" << endl
              << "sname,MADI,dim,dbSize,collection,nQueries,radius,nfound,aveNfound,radiusSum,Pivot,Partition,"
                 "NNodesVistd,NDistCalls,btime,stime,stime_n1000"
              << endl;
		csvFile.close();
      csvh = false;
    }

    std::set<unsigned int> depths;
    double radiusSum = stree.radiusSumAndDepths(depths);

    unsigned int maxResults = 100000000;
    unsigned int nFound;
    unsigned int nqActual;
    for (auto& radiusM : radii) {
      if (collectPlusRange == true) {  // DO NORMAL RADIUS SEARCH IF SO SELECTED.
        cout << "starting collectSearchTest; regualr range search. radiusM=" << radiusM << endl;


        stree.getPerfStats().reset();
        start = std::clock();
        nFound = 0;
        nqActual = 0;
        for (auto& qp : qObjects) {
          auto searchRadius = radiusM;
          if (radiusMult == true) {
            searchRadius = qp.getValue().size() * radiusM;
          }
          RadiusQuery<ObjType> rq(qp, searchRadius, maxResults);
          stree.search(rq);
          nqActual++;
          nFound += rq.getNeighbors().size();
          rq.getNeighbors().clear();
          if ((nqActual % 1000) == 0) {
            cout << "Finished search i= " << nqActual << endl;
          }
        }
        sTime = dTimeSeconds(start);

		ofstream csvFile(fileNamePrefix + "_crs.csv", ios::app);
        csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << ","
                << "F"
                << "," << nqActual << "," << radiusM << "," << nFound << "," << ((nFound * 1.0f) / nqActual) << ","
                << radiusSum << "," << stree.getPivotType() << "," << stree.getPartType() << ","
                << (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
                << static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << "," << bTime << ","
                << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();
      }

      // DO COLLECTION SEARCH TEST
      cout << "starting collectSearchTest commect compone t radius=" << radiusM << endl;
      stree.getPerfStats().reset();
      start = std::clock();
      nFound = 0;
      nqActual = 0;
      for (auto& qp : qObjects) {
        auto searchRadius = radiusM;
        if (radiusMult == true) {
          searchRadius = qp.getValue().size() * radiusM;
        }
        RadiusQuery<ObjType> rq(qp, searchRadius, maxResults);
        stree.searchCollect(rq);
        nqActual++;
        nFound += rq.getNeighbors().size();
        rq.getNeighbors().clear();
        if ((nqActual % 1000) == 0) {
          cout << "Finished search i= " << nqActual << endl;
        }
      }
      sTime = dTimeSeconds(start);

      if (csvh == true) {
		ofstream csvFile(fileNamePrefix + "_crs.csv", ios::app);
        csvFile << ";;TEST COLLECTION: " << currentDateTime() << ";;[CLasses]:" << endl
                << ";;[" << typeid(stree).name() << "],"
                << "[" << met << "]" << endl
                << "sname,MADI,dim,dbSize,collection, "
                   "nQueries,radius,nfound,aveNfound,radiusSum,Pivot,Partition,NNodesVistd,NDistCalls,btime,stime,"
                   "stime_n1000"
                << endl;
        csvh = false;
      }

	  ofstream csvFile(fileNamePrefix + "_crs.csv", ios::app);
      csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << ","
              << "T"
              << "," << nqActual << "," << radiusM << "," << nFound << "," << ((nFound * 1.0f) / nqActual) << ","
              << radiusSum << "," << stree.getPivotType() << "," << stree.getPartType() << ","
              << (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
              << static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << "," << bTime << "," << sTime
              << "," << (sTime * 1000.0 / nqActual) << endl;
      csvFile.close();
    }
}


/**
 * 	Perform kNN search with the query object qObjects on dataset objects.
 *  The radii vector should be tthe size of qObjects.
 *  
 **/
template <class ObjType, class MetricType, class TreeType>
void nkSearchTest(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects,  std::vector<float>& radii,
	 unsigned int avgDim, const std::vector<unsigned int> & maxResultsVec, PivotType pivT, PartType partT, 
	 float madi,  const std::string& fileNamePrefix, bool csvh, float radMultOut = std::numeric_limits<float>::max()) {
	using namespace std;
	unsigned int bTime;

	MetricType met;
    
    auto start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "nkSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	for (auto& maxResults : maxResultsVec) {
		cout << "starting nkSearchTest maxResults=" << maxResults << endl;
		stree.getPerfStats().reset();
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		auto sTime = std::clock();
		for (int i = 0; i< qObjects.size() ; i++ ) {
			auto rad = radii[i];
			NearestKQuery<ObjType> nkQ(qObjects[i], maxResults, rad );
			//nkQ.setSoftMaxResults(true);
			stree.search(nkQ);
			nqActual++;
			nFound += nkQ.getNeighbors().size();

			if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
		     std::string theId = qObjects[i].getId();
			//txtFile << theId <<":";
			//for( const auto& nb : nkQ.getNeighbors()){
				//txtFile << " " << nb;
			//}
		    // txtFile << std::endl;
			nkQ.getNeighbors().clear();
		}
		//txtFile.close();
		sTime = dTimeSeconds(sTime);

		ofstream csvFile(fileNamePrefix + "_knn.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST KNN: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim,dbSize,nQueries,radMult,maxResults, nFound,aveNfound, radiusSum,Pivot,Partition,NNodesVstd,NDistCalls,btime,stime,stime_N1000"
				<< endl;
			csvh = false;
		}

		

		csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << "," 
			<< nqActual << "," << radMultOut << "," << maxResults << "," << nFound << "," << ((nFound * 1.0f) / nqActual) << "," 
			<< radiusSum << "," << stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();

	}
}

template <class ObjType, class MetricType, class TreeType>
void rnkSearchTest(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects,  std::vector<float>& radiM,
	 unsigned int avgDim, const std::vector<unsigned int> & maxResultsVec, PivotType pivT, PartType partT, 
	 float madi,  const std::string& fileNamePrefix, bool csvh) {
	using namespace std;
	unsigned int bTime;

	MetricType met;
    
    auto start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "nkSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	for (const auto& radm : radiM){

	for (auto& maxResults : maxResultsVec) {
		cout << "starting nkSearchTest maxResults=" << maxResults << endl;
		stree.getPerfStats().reset();
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		auto sTime = std::clock();
		for (int i = 0; i< qObjects.size() ; i++ ) {
			auto rad = qObjects[i].getValue().size() * radm;
			NearestKQuery<ObjType> nkQ(qObjects[i], maxResults, rad );
			//nkQ.setSoftMaxResults(true);
			stree.search(nkQ);
			nqActual++;
			nFound += nkQ.getNeighbors().size();

			if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
		     std::string theId = qObjects[i].getId();
			//txtFile << theId <<":";
			//for( const auto& nb : nkQ.getNeighbors()){
				//txtFile << " " << nb;
			//}
		    // txtFile << std::endl;
			nkQ.getNeighbors().clear();
		}
		//txtFile.close();
		sTime = dTimeSeconds(sTime);

		ofstream csvFile(fileNamePrefix + "_rnk.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST KNN: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim|avgLen,dbSize,nQueries,radMult,maxResults, nFound,aveNfound, radiusSum,Pivot,Partition,NNodesVstd,NDistCalls,btime,stime,stime_N1000"
				<< endl;
			csvh = false;
		}

		csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << "," 
			<< nqActual << "," << radm << "," << maxResults << "," << nFound << "," << ((nFound * 1.0f) / nqActual) << "," 
			<< radiusSum << "," << stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();

	}
	}
}

//////////// NK+RS Search tests ///////

template <class ObjType, class MetricType, class TreeType>
void nkRsSearchTest(std::vector<ObjType> & objects, std::vector<ObjType>& qObjects, 
	 unsigned int avgDim, const std::vector<unsigned int> & maxResultsVec, PivotType pivT, PartType partT, 
	 float madi,  const std::string& fileNamePrefix, bool csvh) {
	using namespace std;

	unsigned int bTime, sTime;
	MetricType  met;

	 auto start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "nkSearchTest btime=" << bTime << endl;

	std::map<std::string, double> radMap;
	std::map<std::string, unsigned int > cntMap;
	
	for (auto& maxResults : maxResultsVec) {
		start=std::clock();
		cout << "starting nkSearchTest  maxResults=" << maxResults << endl;
		stree.getPerfStats().reset();
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		radMap.clear();
		cntMap.clear();
		auto rad{ std::numeric_limits<float>::max() };
		sTime = std::clock();
		for (auto& qp : qObjects) {
			NearestKQuery<ObjType> nkQ(qp,  maxResults, rad);
			nkQ.setSoftMaxResults(true);
			stree.search(nkQ);
			nqActual++;
			nFound += nkQ.getNeighbors().size();
			radMap[qp.getId()] = nkQ.searchRadius();
			cntMap[qp.getId()] = nkQ.getNeighbors().size();
			if ((nqActual % 1000) == 0) { cout << "Finished NK search i= " << nqActual << endl; }
			nkQ.getNeighbors().clear();
		}
		sTime = dTimeSeconds(sTime);

		ofstream csvFile(fileNamePrefix + "_nk_rs.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST nK_rs: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MDI,AvgDIM,dbSize,nQueries,radius, maxResults,nfound,aveNfound,Pivot,Partition,NNodesVstd,NDistCalls,btime,stime,stime_1000" << endl;
			csvh = false;
		}
		csvFile << stree.shortName() << "," <<stree.getMADIorK() <<"," << avgDim << "," << objects.size() << "," << nqActual << "," 
		<< rad <<"," << maxResults << "," << nFound << "," << ((nFound * 1.0f) / nqActual) << "," 
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
		<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
		<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;

		//Do the equivalent range query.
		unsigned int nFoundRS = 0;
		double rad_div_len = 0;
		stree.getPerfStats().reset();
		sTime = std::clock();
		for (auto& qp : qObjects) {
			start=std::clock();
			auto radius = radMap[qp.getId()];
			rad_div_len += radius / avgDim;
			//TODO: Correct radius for acccuracy radius += 0.0000001;//TODO: 
			//radius += 0.0000003;
			RadiusQuery<ObjType> rq(qp, radius, 1000000);
			rq.setSoftMaxResults(false);
			stree.search(rq);
			nFoundRS += rq.getNeighbors().size();
			if (cntMap[qp.getId()] != rq.getNeighbors().size()) {
				std::cout << "kKNN= "<<cntMap[qp.getId()]<<" kRS= "<<rq.getNeighbors().size()
						<<	" KrAD= " << radius << std::endl;
				//error("nkQ != nkQ2");
			}
		//	if ((nqActual % 1000) == 0) { cout << "Finished nk search RS i= " << nqActual << endl; }
			//rq.getNeigbors().clear();
		}
		sTime = dTimeSeconds(sTime);
		csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << "," << nqActual << ","
			<< (rad_div_len/nqActual) << "," << maxResults << "," << nFoundRS << "," << ((nFoundRS * 1.0f) / nqActual) << "," 
			<< stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << "," << (sTime * 1000.0 / nqActual) << endl;
		csvFile.close();
	}
}


template <class ObjType, class MetricType, class TreeType>
void treeStatsTest(std::vector<ObjType> & objects, unsigned int avgDim,
	const PivotType pivT,  const PartType partT, float madi, 
	const std::string& fileNamePrefix, bool csvh) {
	using namespace std;
	unsigned int bTime, sTime;

	MetricType met;
    
    std::clock_t start = std::clock();
	TreeType stree(objects, met, pivT, partT, madi);
	bTime = dTimeSeconds(start);
	cout << "treeStatsTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	
	CMTStructStat<ObjType> st;
	
	st = stree.getStructureStats() ;

	ofstream csvFile(fileNamePrefix + "_rs.csv", ios::app);
	if (csvh == true) {
			csvFile << ";;STUCTURE STAT: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim,dbSize, log2dbSize,Pivot,Partition,nNodes, minDepth,maxDepth,nPivots,AvgNPivots,nAdis,AvgNAdis,radSum,bTime" << endl;
			csvh = false;
	}

	csvFile << stree.shortName() << "," << stree.getMADIorK() << "," << avgDim << "," << objects.size() << "," <<std::log2(objects.size()) <<","
		
			<< stree.getPivotType() << "," << stree.getPartType() << "," << st << "," << bTime  <<","  << endl;
		csvFile.close();
		
	}



#endif
