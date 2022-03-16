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

/*
	Various search algorithm tests with variations on the Euclidian metric, partition strategies and
	dataset sizes - mostly in support of the original 2019 paper.
*/


void nkIncreasingDensityTestEM(unsigned int nofPoints, unsigned int nQueries, PivotType pivT, PartType partT, unsigned int nofResults,
	const std::string& fileNamePrefix = { "mtreeRunTimes_id" }, bool csvh = false);
void nkIncreasingDensityTestEM(const std::string& fileNamePrefix = { "mtreeRunTimes_id" });


/**
	Generate the Eucclidian Metric DB points and query points.
	nPoints - the number of DB points to generate.
	nQueries - the number of query points desired. 
	
	The query points are selected from the BD points and also
	will be interior (away from the walls) of the box containing the DB points.
	Because of this the number of query points generated may be less than nQueries.
**/
std::tuple<std::vector<EuclidianPoint>, std::vector<EuclidianPoint>>
generatePointsQD(unsigned int nPoints, unsigned int nQueries, const double density = 1.0f) {
		using namespace std;
		constexpr unsigned int dim = EuclidianPointDim;
		auto width = ceil(pow((nPoints / density), 1.0 / dim));

		std::array< EuclidianPointPointType, dim> boxO;
		std::array< EuclidianPointPointType, dim> boxP;

		//Make a box that sorrounds the interior points - if any.
		//In high dimensions (e.g. DIM=9) almost all points are at or near 
		//the exterior shell enclosing the data. Query points near the interior make 
		// for slightly harder queries.
		/* auto frac = 0.15;
		if(dim > 10){
			frac = 0.01;
		}else if (dim > 5) {
			if(nPoints < 1000){
				frac = 0.1;
			}
			if(nPoints > 5000){
				frac = 0.20;
			}
		}
		*/
	
		auto frac = 0.01;
		
		for (int i = 0; i < dim; i++) {
			boxO[i] = 0.0f + width * frac;
			boxP[i] = (1.0 - frac) * width ;
		}

		std::clock_t start = std::clock();
		std::vector<EuclidianPoint> dbPoints =
			genEuclidianPointsUniform<dim, EuclidianPointPointType>(nPoints, 0.0f, width);
		std::cout << "dataGenTime=" << dTimeSeconds(start) << std::endl;

		//And the query points ...
		auto qPoints = getInteriorPoints<float, dim>(dbPoints, boxO, boxP);
		qPoints = getSubsetRandom(qPoints, nQueries);

		//TODO: verify we dont need std::move usage.
		return { dbPoints, qPoints };
}

std::tuple<std::vector<EuclidianPoint>, std::vector<EuclidianPoint>>
generatePointsBrin95(unsigned int nPoints, unsigned int nQueries) {
		using namespace std;
	
		std::clock_t start = std::clock();
		std::vector<EuclidianPoint> dbPoints =
			genEuclidianPointsUniform<EuclidianPointDim, EuclidianPointPointType>(nPoints, 0.0f, 1.0f);
		std::cout << "dataGenTime=" << dTimeSeconds(start) << std::endl;

		//And the query points ...
		//auto qPoints = getInteriorPoints<float, dim>(dbPoints, boxO, boxP);
		//qPoints = getSubsetRandom(qPoints, nQueries);
		auto qPoints = getSubsetRandom(dbPoints, nQueries);

		//TODO: verify we dont need std::move usage.
		return { dbPoints, qPoints };
}

std::tuple<std::vector<EuclidianPoint>, std::vector<EuclidianPoint>>
generatePointsQDClustered(unsigned int nClusters, unsigned int nPointsPer, unsigned int nQueries, 
		const double density = 1.0f) {
	using namespace std;
	unsigned int nPoints = nClusters * nPointsPer;
	constexpr unsigned int dim = EuclidianPointDim;
	auto width = ceil(pow((nPoints / density), 1.0 / dim));
	auto clusterWidth = ceil(pow((nPointsPer / density), 1.0 / dim));
	float sigma = width / (4.0 * clusterWidth);

	std::array< EuclidianPointPointType, dim> boxO;
	std::array< EuclidianPointPointType, dim> boxP;

	//Make a box that sorrounds the interior points.
	//TODO: get better delta
	auto delta = 1;
	for (int i = 0; i < dim; i++) {
		boxO[i] = 0.0f + delta;
		boxP[i] = width - delta;
	}

	using dataType = EuclidianPointPointType;
	
	std::clock_t start = std::clock();
	std::vector<EuclidianPoint> dbPoints =
		generateEuclidianPointsClustered<EuclidianPointDim, dataType>(nClusters, nPointsPer, 0.0f, width, sigma);
	//std::vector<EuclidianPoint> dbPoints = genEuclidianPointsUniform<dim, EuclidianPointPointType>(nPoints, 0.0f, width);
	std::cout << "dataGenTime=" << dTimeSeconds(start) << endl;

	//And the query points ...
	auto qPoints = getInteriorPoints<float, dim>(dbPoints, boxO, boxP);
	qPoints = getSubsetRandom(qPoints, nQueries);

	//TODO: verify we dont need std::move usage.
	return { dbPoints, qPoints };
}


void getDataBounds(std::vector<EuclidianPoint>& points, std::array<EuclidianPointPointType, EuclidianPointDim>& boxO,
              std::array<EuclidianPointPointType, EuclidianPointDim>& boxP) {
  constexpr unsigned int dim = EuclidianPointDim;
  for (int i = 0; i < dim; i++) {
    boxO[i] = boxP[i] = points[0][i];
  }

  for (int j = 1; j < points.size(); j++) {
    for (int i = 0; i < dim; i++) {
      if (points[j][i] < boxO[i]) boxO[i] = points[j][i];
      if (points[j][i] > boxP[i]) boxP[i] = points[j][i];
    }
  }
}

/**
	Use the search tree(s) to perform range search over various combinations of the
	search options (pivot types, radii, db sizes, etc).
	This function is used primarily to gather performance comparison statistics
**/
void radiusSearchTestEM(const std::string& fileNamePrefix) {
	//Set the DB sizes and number of queries that we want to test
	//std::map<unsigned int, unsigned int> nofPoints = getTestSizes(2, 10,18,1, 10000) ;
	//std::map<unsigned int, unsigned int> nofPoints{ {1000,100}, {10000,1000 } ,{100000,1000 } ,{1000000,1000 } };
	std::map<unsigned int, unsigned int> nofPoints{  {10000,1000 }  };

	//std::vector<float> radii {0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f};
	std::vector<float> radii {1.0f};

	//std::vector<float> lcmt_npivs {0.25, 0.5, 1.0, 2.0, 3.0 };
	//std::vector<int> cmt_madis {1000, 1};

	std::vector<float> lcmt_npivs ;
	std::vector<int> cmt_madis {1000};

	std::set<PivotType> includePivotTypes{ PivotType::RAN};
	std::set<PartType>  includePartTypes{ PartType::BOM};
	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto& [np, nQueries] : nofPoints) {
		auto [points, qPoints] = generatePointsQD(np, nQueries);
		//auto [points, qPoints] = generatePointsQDClustered(nClusters ,nPointsPer, nQueries);
		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					radiusSearchTest<EuclidianPoint, EuclidianMetric,CMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, EuclidianPointDim, radii, pivType, partType, madi, fileNamePrefix, hflag);
					hflag = false; //Dont print header after 1st time.
				}
			}
		}
	}
}

/**
	Use the search tree(s) to perform range search over various combinations of the
	search options (pivot types, radii, db sizes, etc). Does not store resuts - only count them.
	This function is used primarily to gather performance comparison statistics
**/
void radiusSearchTestEMCount(const std::string& fileNamePrefix) {
	//Set the DB sizes and number of queries that we want to test
	//std::map<unsigned int, unsigned int> nofPoints = getTestSizes(2, 10,18,1, 10000) ;
	std::map<unsigned int, unsigned int> nofPoints{ {1000,100}, {10000,1000 } ,{100000,1000 } ,{1000000,1000 } };
	
	std::vector<float> radii {0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 4.0f, 5.0f };

	std::vector<int> cmt_madis {1000, 10, 5, 1, 0};

	std::set<PivotType> includePivotTypes{ PivotType::RAN};
	std::set<PartType>  includePartTypes{ PartType::BOM };
	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto& [np, nQueries] : nofPoints) {
		auto [points, qPoints] = generatePointsQD(np, nQueries);
		//auto [points, qPoints] = generatePointsQDClustered(nClusters ,nPointsPer, nQueries);
		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					radiusSearchTestCount<EuclidianPoint, EuclidianMetric,CMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, EuclidianPointDim, radii, pivType, partType, madi, fileNamePrefix, hflag);
					hflag = false; //Dont print header after 1st time.
				}
			
			}
		}
	}
}




////// NK Searches ////////////////////////////

void nkSearchTestEM(const std::string &fileNamePrefix){
	//std::map<unsigned int, unsigned int> nofPoints = getTestSizes(2, 10,18,1, 10000) ;
	std::map<unsigned int, unsigned int> nofPoints{{1000, 100}, {10000, 100}, {100000, 100}, {1000000, 100},
		{10000000, 100}};
	//std::map<unsigned int, unsigned int> nofPoints{ {10000000, 100}};
	
	std::vector<unsigned int> maxResults{1, 2 , 4, 6, 8,  10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000};
	
	std::vector<float> lcmt_npivs; //{0.5, 1, 2, 4}
	std::vector<int> cmt_madis {1,1000};
	std::vector<int> spmt_madis{0};


	std::set<PivotType> includePivotTypes{PivotType::RAN};
	std::set<PartType> includePartTypes{PartType::BOM};

	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto &[np, nQueries] : nofPoints){
		auto [points, qPoints] = generatePointsQD(np, nQueries);
		std::vector<float> radii;	
		for (int i = 0; i < nQueries; i++) {
       		radii.push_back(std::numeric_limits<float>::max());
     	}
		for (const auto &pivType : includePivotTypes){
			for (const auto &partType : includePartTypes){
				for (const auto& madi : cmt_madis){
					nkSearchTest<EuclidianPoint, EuclidianMetric,CMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, radii, EuclidianPointDim, maxResults, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) {	hflag = false; 	}
				}
				
				for (const auto& madi : spmt_madis){
					nkSearchTest<EuclidianPoint, EuclidianMetric,SPMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, radii, EuclidianPointDim, maxResults, pivType, partType,madi, fileNamePrefix, hflag);
					if (hflag == true) {hflag = false;}
				}
			}
		}
	}
}




void nkRsSearchTestEM(const std::string& fileNamePrefix) {
	unsigned int bTime, sTime;
	//std::map<unsigned int, unsigned int> nofPoints{ {1000,100}, {10000,1000 } ,{100000,1000 } ,{1000000,1000 } };
	std::map<unsigned int, unsigned int> nofPoints{ {1000000,100}  };

	std::vector<unsigned int> maxResults{1, 2, 5, 10, 20};
	//std::vector<unsigned int> maxResults{1};

	std::vector<float> lcmt_npivs ;// {1, 2};
	std::vector<int> cmt_madis {1000, 1};
	std::vector<int> spmt_madis {0};

	std::set<PivotType> includePivotTypes{ PivotType::RAN};
	std::set<PartType>  includePartTypes{ PartType::BOM};
	bool hflag = true;
	for (const auto& [np, nQueries] : nofPoints) {
		auto [allPoints, qPoints2] = generatePointsQD(np + nQueries, nQueries);

		std::vector<EuclidianPoint> points (allPoints.begin(), allPoints.begin() + (allPoints.size() - nQueries));
		std::vector<EuclidianPoint> qPoints(allPoints.begin() + (allPoints.size() - nQueries), allPoints.end());

		for (const auto& pivType : includePivotTypes) {
            for (const auto& partType : includePartTypes) {
              for (const auto& madi : cmt_madis) {
                nkRsSearchTest<EuclidianPoint, EuclidianMetric, CMTree<EuclidianPoint, EuclidianMetric>>(
                    points, qPoints, EuclidianPointDim, maxResults, pivType, partType, madi, fileNamePrefix, hflag);
                hflag = false;
              }
              for (const auto& madi : spmt_madis) {
                nkRsSearchTest<EuclidianPoint, EuclidianMetric, SPMTree<EuclidianPoint, EuclidianMetric>>(
                    points, qPoints, EuclidianPointDim, maxResults, pivType, partType, madi, fileNamePrefix, hflag);
                hflag = false;
              }
            }
        }
     }
}


///////// NK search tests ///////////////////
void nkSearchTestEMC(unsigned int nClusters, unsigned int nPointsPer, unsigned int nQueries,
	PivotType pivT, PartType partT, const std::vector<unsigned int> & maxResultsVec ,
	const std::string& fileNamePrefix, bool csvh) {
	using namespace std;

	unsigned int bTime, sTime;
	using MetricType = EuclidianMetric;
	MetricType  met;

	constexpr unsigned int dim = EuclidianPointDim;
	unsigned int nPoints = nClusters * nPointsPer;
	
	auto [points, qPoints] = generatePointsQD(nPoints, nQueries);
	//auto [points, qPoints] = generatePointsQDClustered(nClusters, nPointsPer, nQueries);

	bTime = std::clock();
	CMTree<EuclidianPoint, MetricType> stree(points, met, pivT, partT, 3);
	bTime = dTimeSeconds(bTime);
	cout << "firstkSearchTest btime=" << bTime << endl;


	for (auto& maxResults : maxResultsVec) {
		cout << "starting nkSearchTestEM maxResults=" << maxResults << endl;
		stree.getPerfStats().reset();
		unsigned int nFound = 0;
		unsigned int nqActual = 0;
		auto rad{ std::numeric_limits<float>::max() };
		sTime = std::clock();
		for (auto& qp : qPoints) {
			NearestKQuery<EuclidianPoint> nkQ(qp,  maxResults, rad);
			stree.search(nkQ);
			nqActual++;
			nFound += nkQ.getNeighbors().size();
			if ((nqActual % 10000) == 0) { cout << "Finished search i= " << nqActual << endl; }
		}
		sTime = dTimeSeconds(sTime);

		ofstream csvFile(fileNamePrefix + "_EM__nk.csv", ios::app);
		if (csvh == true) {
			csvFile << ";;TEST NEARK_EM: " << currentDateTime()
				<< ";;[CLasses]:" << endl
				<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
				<< "sname,MADI,dim,dbSize,nClusters,nQueries,maxResults,nfound,aveNfound,Pivot,Partition,nodesVisited,numDistCalls,btime,stime" << endl;
		}
		csvh = false;

		csvFile << stree.shortName() << "," <<stree.getMADIorK() <<"," << dim << "," << points.size() << "," << nClusters << "," 
			<< nqActual << "," << maxResults << ","<< nFound << "," << ((nFound * 1.0f) / nqActual) << "," 
			<< stree.getPivotType() << "," << stree.getPartType() << ","
			<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
			<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
			<< bTime << "," << sTime << endl;
	csvFile.close();
	}
}



void nkSearchTestEMC(const std::string& fileNamePrefix) {
//Set the DB sizes and number of queries that we want to test
	//std::vector<unsigned int> maxResults{ 1,5 };
	std::vector<unsigned int> maxResults{ 1,3,5,10};
	unsigned int nPoints = 1000000;	
	unsigned int nQueries = 1000;
	unsigned const int divisor = 10;
	unsigned int nPointsPer = nPoints;
	unsigned int nClusters = 1;

	std::set<PivotType> includePivotTypes{ PivotType::RAN };
	std::set<PartType>  includePartTypes{PartType::BOM };
	bool hflag = true;
	while (nPointsPer >= 100){
		for (const auto &pivType : includePivotTypes){
			for (const auto &partType : includePartTypes){
				nkSearchTestEMC(nClusters, nPointsPer, nQueries, pivType, partType, maxResults, fileNamePrefix, hflag);
				if (hflag == true){
					hflag = false; //only print cvs header info for first data set
				}
			}
		}
		nClusters *= divisor;
		nPointsPer /= divisor;
	}
}


/**

**/

void collectCountSearchTestEM(const std::string& fileNamePrefix) {
	std::map<unsigned int, unsigned int> nofPoints{ {100000,100 } };
	
	//std::vector<float> radii { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 8.5f, 9.0f, 9.5f ,10.0f};
	std::vector<float> radii { 1.0f, 2.0f, 4.0f,  6.0f,  8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 20.0f, 22.0f, 24.0f, 26.0f, 28.0f, 30.0f};

//std::vector<float> radii {5.0f};

	std::set<PivotType> includePivotTypes{ PivotType::RAN}; 
	std::set<PartType>  includePartTypes{ PartType::BOM};

	std::vector<int> cmt_madis {1,1000};

	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto& [np, nQueries] : nofPoints) {
		auto [points, qPoints] = generatePointsQD(np, nQueries);
		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					collectCountSearchTest<EuclidianPoint, EuclidianMetric,CMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, EuclidianPointDim, radii, false, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) {hflag = false;}
				}
			}
		}
	}
}

// Perform collection count with CMTs, and also radius counts
void collectCountSearchPlusTestEM(const std::string& fileNamePrefix) {
  std::map<unsigned int, unsigned int> nofPoints{{100000, 100}};

   std::vector<float> radii { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 8.5f, 9.0f, 9.5f ,10.0f};
  // std::vector<float> radii
  // { 1.0f, 2.0f, 4.0f,  6.0f,  8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 20.0f, 22.0f, 24.0f, 26.0f, 28.0f, 30.0f};


  std::set<PivotType> includePivotTypes{PivotType::RAN};
  std::set<PartType> includePartTypes{PartType::BOM};

  std::vector<float> lcmt_npivs{1, 2, 4};
  std::vector<int> cmt_madis{1, 1000};

  bool hflag = true;
  std::vector<EuclidianPoint> points, qPoints;

  // The counting queries
  for (const auto& [np, nQueries] : nofPoints) {
    auto [points, qPoints] = generatePointsQD(np, nQueries);
    for (const auto& pivType : includePivotTypes) {
      for (const auto& partType : includePartTypes) {
        for (const auto& madi : cmt_madis) {
          collectCountSearchTest<EuclidianPoint, EuclidianMetric, CMTree<EuclidianPoint, EuclidianMetric>>(
              points, qPoints, EuclidianPointDim, radii, false, pivType, partType, madi, fileNamePrefix, hflag);
          if (hflag == true) {
            hflag = false;
          }
        }
      }
    }
  }

  // The range queries
  hflag = true;
  for (const auto& [np, nQueries] : nofPoints) {
    auto [points, qPoints] = generatePointsQD(np, nQueries);
    for (const auto& pivType : includePivotTypes) {
      for (const auto& partType : includePartTypes) {
        for (const auto& madi : cmt_madis) {
          radiusSearchTest<EuclidianPoint, EuclidianMetric, CMTree<EuclidianPoint, EuclidianMetric>>(
              points, qPoints, EuclidianPointDim, radii, pivType, partType, madi, fileNamePrefix, hflag);
          hflag = false;
        }
      }
    }
  }
}

void collectSearchTestEM(const std::string& fileNamePrefix) {
	std::map<unsigned int, unsigned int> nofPoints{ {1000000,100 } };
	//std::vector<float> radii {0.25f, 0.5f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f}; //For DIM=10
    std::vector<float> radii {1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f, 120.0f, 130.0f, 140.0f}; //For DIM=3 1M



	std::set<PivotType> includePivotTypes{ PivotType::RAN}; 
	std::set<PartType>  includePartTypes{ PartType::BOM};

	std::vector<int> cmt_madis {1, 1000};

	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto& [np, nQueries] : nofPoints) {
		auto [points, qPoints] = generatePointsQD(np, nQueries);
		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					collectSearchTest<EuclidianPoint, EuclidianMetric,CMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, EuclidianPointDim, radii, false, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) {hflag = false;}
				}
			}
		}
	}
}


/*
	Collectionand normal rad search. Radii are auto determined from the
	boundaries of the dataset. There are nrad (e.g. 10) evently psaced small radii and
	nrad evenaly spaced large radii. The largest radii is equal to the 
	largest distance possible in the hyper-cubed space.

*/
void collectSearchTestEMAutoRad(const std::string& fileNamePrefix, int nrad, float maxRadPct) {
	std::map<unsigned int, unsigned int> nofPoints{{10000,100 } ,{100000,100 } , {1000000,100 } , {10000000,100 }};
	//std::map<unsigned int, unsigned int> nofPoints{{10000,100 } ,{100000,100 } , {1000000,100 }};
	//std::map<unsigned int, unsigned int> nofPoints{{10000000,100 }};
    std::vector<float> radii;
	bool rangeSearchToo = true;
	

	std::set<PivotType> includePivotTypes{ PivotType::RAN}; 
	std::set<PartType>  includePartTypes{ PartType::BOM };

	std::vector<int> cmt_madis{1,1000}; // {1};
	std::vector<int> spmt_madis{0};// {0};


	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto& [np, nQueries] : nofPoints) {
		auto [points, qPoints] = generatePointsQD(np, nQueries);
		std::array<EuclidianPointPointType, EuclidianPointDim> boxO,boxP;
		getDataBounds(points, boxO, boxP);
		EuclidianPoint pO("0", boxO);
		EuclidianPoint pP("P", boxP);
		auto maxDist = pO.distance(pP);
		maxDist *= maxRadPct;

		float delta = 3.0 / nrad;

		radii.clear();
		for (int i = 0; i <= nrad; i++){
			radii.push_back(i * delta * maxDist);
		}


		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					collectSearchTest<EuclidianPoint, EuclidianMetric,CMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, EuclidianPointDim, radii, false, pivType, partType, madi, fileNamePrefix, hflag, rangeSearchToo);
					if (hflag == true) {hflag = false;}
				}
				for (const auto& madi : spmt_madis){
					collectSearchTest<EuclidianPoint, EuclidianMetric,SPMTree<EuclidianPoint, EuclidianMetric>>
					(points, qPoints, EuclidianPointDim, radii, false, pivType, partType, madi, fileNamePrefix, hflag, rangeSearchToo);
					if (hflag == true) {hflag = false;}
				}
			}
		}
	}
}



void radiusSearchCompareEM(unsigned int nPoints, const unsigned int nQueries, PivotType pivT, PartType partT, 
		const std::string& fileNamePrefix) {
	using namespace std;
	unsigned int bTime, sTime;
	using MetricType = EuclidianMetric;
	MetricType  met;
	unsigned int nqActual = 0;

	auto [points, qPoints] = generatePointsBrin95(nPoints, nQueries);

	auto start = std::clock();
	//SPMTree<EuclidianPoint, MetricType> stree(points, met);
	CMTree<EuclidianPoint, MetricType> stree(points, met,pivT, partT, kxBalancedTreeHeight(1,points.size()));
	BruteForceSearch<EuclidianPoint, MetricType> stree2(points, met);
	bTime = dTimeSeconds(start);
	cout << "radiusSearchTest btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	start = std::clock();
	unsigned int nFound = 0;
	unsigned int diffCount = 0;
	const unsigned int maxResults = 1000000;
	auto rad = 0.0;
	for (const auto& qp : qPoints) {
		/*if(nqActual == 2) {
			cout << "pq id="<<qp.getId() << endl;
			RadiusQuery<EuclidianPoint> rqt(qp, rad, maxResults);
			stree.search(rqt);
		}*/
		RadiusQuery<EuclidianPoint> rq(qp, rad, maxResults);
		stree.search(rq);
	
		nFound += rq.getNeighbors().size();

		RadiusQuery<EuclidianPoint> rq2(qp, rad, maxResults);
		stree2.search(rq2);

		if (!rq.hasSameNeighbors(rq2)) {
			diffCount++;
			/*
			auto missingIn1 = rq.missingNeighbors(rq2);
			auto missingIn2 = rq2.missingNeighbors(rq);
			std::cout << "id=" << rq.getTarget() << std::endl;
			std::cout << "missingIn1" << std::endl;
			for (const auto& n : missingIn1) {
				std::cout << "d =" << n.distance << " id=" << n.id << std::endl;
			}
			std::cout << "missingIn2" << std::endl;
			for (const auto& n : missingIn2) {
				std::cout << "d =" << n.distance << " id=" << n.id << std::endl;
			}
			
			std::cout << "---listing nbrs---:" << std::endl;
			auto nba = rq.getNeigbors();
			auto nbb = rq2.getNeigbors();
			for (const auto& n : nba) {
				std::cout << "d =" << n.distance << " id=" << n.id << std::endl;
			}
			for (const auto& n : nbb) {
				std::cout << "bd=" << n.distance << " id=" << n.id << std::endl;
			}
			
			
			std::cout << "-----" << std::endl;
			int i;
			std::cin >> i;
			*/
		}

	    nqActual++;
		if ((nqActual % 10000) == 0) { cout << "finished search i= " << nQueries << endl; }
	}
	sTime = dTimeSeconds(start);

	string runTimeFileName{ fileNamePrefix + "_rs_compare.csv"};
	ofstream theFileA(runTimeFileName, ios::app);
	theFileA << ";;TEST EMSC: " << currentDateTime()
		<< ";;[CLasses]:" << endl
		<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
		<< ";;[CLasses # 2]:" << endl
		<< ";;[" << typeid(stree2).name() << "]" << endl
		<< ";; Tree depths, min depth, max depth :" << depths.size() << ","
		<< *(std::min_element(depths.begin(), depths.end())) << "," << *(std::max_element(depths.begin(), depths.end())) << endl
		<< "dbSize,DiffCount,nQueries,nfound,radiusSu,Pivo,Partitio,<nodesVisited>,<numDistanceCalls>,btime,stime" << endl
		<< points.size() << "," << diffCount << "," << nqActual << "," << nFound << "," << radiusSum << ","
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (stree.getPerfStats().getNodesVisited() / nQueries) << ","
		<< (stree.getPerfStats().getDistanceCalls() / nQueries) << ","
		<< bTime << "," << sTime << endl;

	theFileA.close();
}
//////////////////////


void kNNSearchCompareEM(unsigned int nPoints, unsigned int nQueries,
	PivotType pivT, PartType partT, unsigned int nofResults,
	const std::string& fileNamePrefix, bool csvh = true) {
	using namespace std;

	unsigned int bTime, sTime;
	using MetricType = EuclidianMetric;
	MetricType  met;

	constexpr unsigned int dim = EuclidianPointDim;
	
	//auto [points, qPoints] = generatePointsQD(nPoints, nQueries);
	auto [points, qPoints] = generatePointsBrin95(nPoints, nQueries);


	std::clock_t start = std::clock();
	CMTree<EuclidianPoint, MetricType> stree(points, met, pivT, partT, kxBalancedTreeHeight(1, points.size()));
	SPMTree<EuclidianPoint, MetricType> stree2(points, met, pivT, partT, kxBalancedTreeHeight(1, points.size()));
	// BruteForceSearch<EuclidianPoint, MetricType> stree2(points, met);

	bTime = dTimeSeconds(start);
	cout << "firstkSearchTest btime=" << bTime << endl;

	start = std::clock();
	unsigned int diffCount = 0;
	unsigned int nFound = 0;
	unsigned int nqActual = 0;
	auto rad{ std::numeric_limits<float>::max() };
	for ( auto& qp : qPoints) {
		NearestKQuery<EuclidianPoint> rq(qp,  nofResults, rad);
		stree.search(rq);
		nqActual++;
		nFound += rq.getNeighbors().size();

		NearestKQuery<EuclidianPoint> rq2(qp,  nofResults, rad);
		stree2.search(rq2);

		if (!rq.hasSameNeighbors(rq2)) {
			diffCount++;
		}

		if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
	}
	sTime = dTimeSeconds(start);

	//AND the CSV file:
	ofstream csvFile(fileNamePrefix + ".csv", ios::app);
	if (csvh == true) {
		csvFile << ";;TEST NEARK_COMPARE_EM: " << currentDateTime()
			<< ";;[CLasses]:" << endl
			<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
			<< "DiffCnt,name,MADI,dim,dbSize,nQueries,maxResults,nfound,aveNfound,Pivot,Partition,nodesVisited,numDistCalls," << endl;
	}
	csvFile << diffCount << "," << stree.shortName() << "," << stree.getMADIorK() <<"," << dim << "," << points.size() << "," << nqActual << "," << nofResults << ","
		<< nFound << "," << ((nFound * 1.0f) / nqActual) << ","
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
		<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
		<< static_cast<float>(stree2.getPerfStats().getDistanceCalls()) / nqActual << ","<<endl;
		// << (1.0f * stree.getPerfStats().getDistanceCalls()) / (1.0f * nFound)  << ","
		// << (1.0f * stree2.getPerfStats().getDistanceCalls()) / (1.0f * nFound) <<  endl;
	csvFile.close();
}


void radiusSearchCompareEM(const std::string& fileNamePrefix) {
	//std::map<unsigned int, unsigned int> nofPoints{ {100,1},{1000,10}, {10000,10}, {1000000,100} };
	std::map<unsigned int, unsigned int> nofPoints{ {100000,10} };
	for (const auto& [np, nSkip] : nofPoints) {
		for (const auto& [pivType, pivVal] : pivotTypeMap) {
			for (const auto& [parType, parVal] : partTypeMap) {
				///cout << np << " " << pivType << " " << pivVal << endl; break;
				radiusSearchCompareEM(np, nSkip, pivType, parType, fileNamePrefix);
			}
		}
	}
}

void kNNSearchCompareEM(const std::string& fileNamePrefix) {
	//std::map<unsigned int, unsigned int> nofPoints{ {100,1},{1000,10}, {10000,10}, {1000000,100} };
	std::map<unsigned int, unsigned int> nofPoints{ {100000,1000} };
	std::vector<unsigned int> maxResults{ 1,5, 10 };
	//
	for (const auto& [np, nQueries] : nofPoints) {
		for (const auto& [pivType, pivVal] : pivotTypeMap) {
			for (const auto& [parType, parVal] : partTypeMap) {
				for (const auto& mr : maxResults) {
					kNNSearchCompareEM(np, nQueries, pivType, parType, mr, fileNamePrefix );
				}
			}
		}
	}
}

void PCPETest() {
	using namespace std;


	struct Node {
		EuclidianPoint object;
		Node(EuclidianPoint object) : object(object) {}
	};
	using NodePtr = std::shared_ptr<Node>;
	using NodeItr = std::vector<NodePtr>::iterator;

	using MetricType = EuclidianMetric;
	MetricType  met;

	//Pivot will be the point in the "center":

	constexpr unsigned nPoints = 100;
	constexpr unsigned int dim = EuclidianPointDim;
	constexpr double density = 1.0;
	auto width = ceil(pow((nPoints / density), 1.0 / dim));
	std::vector<EuclidianPoint> points =
		genEuclidianPointsUniform<dim, EuclidianPointPointType>(nPoints, 0.0f, width);
	std::vector<NodePtr> nodes;
	for (auto& point : points) {
		nodes.push_back(std::make_shared <Node>(point));
	}

	auto nodeItr = findCenter<EuclidianPoint, NodeItr>(nodes.begin(), nodes.end(), met, 10);

	cout << (*nodeItr)->object << endl;

	NodeItr extrema = findExtrema<EuclidianPoint, NodeItr>(nodes.begin(), nodes.end(), met);

	cout << (*extrema)->object << endl;

	std::sort(points.begin(), points.end(), [](EuclidianPoint& lp, EuclidianPoint& rp) {
		auto& l = lp.getValue();
		auto& r = rp.getValue();
		if (l[0] < r[0]) return true;
		else if ((l[0] == r[0]) && (l[1] < r[1])) return true;
		else return false; });
	for (auto& point : points) {
		cout << point << endl;
	}

	string runTimeFileName{ "PCPETest.txt" };
	ofstream theFileA(runTimeFileName + "ARSS", ios::app);
	theFileA << "EMS Test: " << currentDateTime() << "classes:" << endl;
	theFileA.close();
}

/*
	Generate a data file for plotting and visual inspection.
*/
void clusteredDataGenTest() {
	using namespace std;
	// distType = std::normal_distribution<dataType>;
	//std::vector<EuclidianPoint> points =
	//	generateEuclidianPointsClustered<EuclidianPointDim, EuclidianPointPointType>(2, 10, 0, 10, 0.25);
	std::vector<EuclidianPoint> points =
		generateEuclidianPointsClustered<EuclidianPointDim, EuclidianPointPointType>(4, 25, 0, 10, 0.25);
	ofstream pFile("points_2x10B.csv", ios::app);
	pFile << "x,y" << endl; //panda plotter needs the labes
	for (const auto& p : points) {
		pFile << p[0] << "," << p[1] << endl;
	}
}



