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
#include "CMTree.h"
#include "SPMTree.h"
#include "APMTree.h"
#include "BruteForceSearch.h"
#include "GPCharacters.hpp"
#include "Reader.h"
#include "TestCommon.h"

auto generateWords(unsigned int nWords, unsigned int aveWordLength, unsigned int sigma) {
	std::vector<Sequence> words;
	words.reserve(nWords);
	std::ostringstream idStream;
	//std::normal_distribution<float> dist(aveWordLength, sigma);
	std::uniform_int_distribution<> dist(1, 2 * aveWordLength);
	std::default_random_engine generator;
	//std::mt19937 generator(3001);
	 int wSize = 0;
	for (int id = 0; id < nWords; id++) {
		wSize = dist(generator);
		while (wSize < 1) {
			wSize += aveWordLength;
		}
		idStream << id;
		words.push_back(Sequence(idStream.str(), CharacterSets::randomEnglish(wSize)));
		idStream.str(std::string());
		idStream.clear();
	}
	return words;
}

/**
	Use the search tree(s) to perform range search over the various combinations of the
	search options (pivot types, radii, db sizes, etc).
	This function is used primarily to gather performance comparison statistics
**/
void radiusSearchTestEDM(const std::string& fileNamePrefix, const std::string& inputDataFName) {
	unsigned int bTime;
	std::vector<float> radii{ 1.0f, 2.0f};

	//std::map<unsigned int, unsigned int> nofPoints = getTestSizes(2, 10,18,1, 10000) ;
	//std::map<unsigned int, unsigned int> nofPoints = getTestSizes(10, 2, 7, 1, 10000);
	//std::map<unsigned int, unsigned int> nofPoints{ {1000,1000 }, {10000,1000 }, {100000,1000 },  {1000000,1000 } };
	std::map<unsigned int, unsigned int> nofPoints{ {100000,1000 } };
	
	
	//std::set<PivotType> includePivotTypes{ PivotType::RAN, PivotType::EXT};//FOR CMT
	//std::set<PartType>  includePartTypes{ PartType::DMR,PartType::OM, PartType::BOM };

	std::set<PivotType> includePivotTypes{ PivotType::RAN, PivotType::EXT };
	std::set<PartType>  includePartTypes{  PartType::BOM, PartType::DMR };

	
	std::vector<int> cmt_madis {1000, 1};
	std::vector<int> amt_madis {0};

	auto start = std::clock();
	Reader reader;
	auto allWords = reader.read(inputDataFName, 100000);
	bTime = dTimeSeconds(start);
	std::cout << "radiusSearchTestEDM datagen time=" << bTime << std::endl;
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();

	bool hflag = true;
	for (const auto& [np, nQueries] : nofPoints) {
		auto words = getSubsetRandom(allWords, np);
		auto qSeqs= getSubsetRandom(words, nQueries);
		for (const auto& pivType : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					radiusSearchTest<Sequence, EditDistMetric,CMTree<Sequence, EditDistMetric>>
					(words, qSeqs, aveWordLength, radii, pivType, partType, madi, fileNamePrefix, hflag);
					hflag = false; //Dont print header after 1st time.
				}
				for (const auto& madi : amt_madis){
					radiusSearchTest<Sequence, EditDistMetric,APMTree<Sequence, EditDistMetric>>
					(words, qSeqs, aveWordLength, radii, pivType, partType, madi, fileNamePrefix, hflag);
					hflag = false; //Dont print header after 1st time.
				}
			}
		}
	}
}



void collectSearchTestEDM(const std::string& fileNamePrefix, const std::string& inputDataFName) {
	std::vector<unsigned int> nofPoints{ 1000000};
	const int nQueries = 100;
	//These radii with radiusMult == false are actuall radii
	//std::vector<float> radii {  1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
	//std::vector<float> radii {  1.0f,  40.0f, 80.0f, 120.0f, 160.0f, 200.0f, 240.0f, 280.0f, 320.0f, 360.0f};
	//bool radiusMult = false;


	//Yhese radii with radiusMult==true are actually factors to the radius: rad = lengt * (1.0 - fac)
	//std::vector<float> 	radii {0.0, 0.001, 0.002, 0.003, 0.004, 0.005 };		
	bool radiusMult = true;
	std::vector<float> 	radii {0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};	


	std::set<PivotType> includePivotTypes{ PivotType::RAN}; 
	std::set<PartType>  includePartTypes{ PartType::BOM};

	std::vector<int> cmt_madis {1,1000};

	auto start = std::clock();
	Reader reader;
	auto allWords = reader.read(inputDataFName, 10000000);
	allWords = randomize(allWords);
	auto bTime = dTimeSeconds(start);
	std::cout << "radiusSearchTestEDM datagen time=" << bTime << std::endl;
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();
	allWords = getSubsetRandom(allWords, 1000000);
	allWords = randomize(allWords);	

	std::vector<Sequence> words(allWords.begin(), allWords.begin() + (allWords.size() - nQueries));
	std::vector<Sequence> qSeqs(allWords.begin() + (allWords.size() - nQueries), allWords.end());
    //auto qSeqs= getSubsetRandom(words, nQueries);

	bool hflag = true;
	for (const auto& np: nofPoints) {
		auto words = getSubsetRandom(allWords, np);
		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					collectSearchTest<Sequence, EditDistMetric,CMTree<Sequence, EditDistMetric>>
					(words, qSeqs, aveWordLength, radii, radiusMult, pivType, partType, madi, fileNamePrefix + "_EXT", hflag);
					if (hflag == true) {hflag = false;}
				}
			}
		}
	}
}

void collectSearchTestEDM_V2(const std::string& fileNamePrefix, const std::string& inputDataFName, int fRad, int lRad,
                             float multRad) {
  unsigned int maxWords = 10000000;
  const int nQueries = 100;
  bool rangeFlag = true;  // false for collection only; true for both collection and normal range search.

  // Radii are claculated here.
  bool radiusMult = false;  // collectSearchTest is to not multiply rad
  std::vector<float> radii;
  for (int i = fRad; i < lRad; i++) {
    radii.push_back(i * multRad);
  }

  std::set<PivotType> includePivotTypes{PivotType::RAN};
  std::set<PartType> includePartTypes{PartType::BOM};

  std::vector<int> cmt_madis{1, 1000};
  std::vector<int> spmt_madis{0};

  auto start = std::clock();
  Reader reader;
  auto allWords = reader.read(inputDataFName, maxWords);
  allWords = randomize(allWords);
  auto bTime = dTimeSeconds(start);
  std::cout << "radiusSearchTestEDM datagen time=" << bTime << std::endl;
  auto lenSum = 0;
  for (const auto& s : allWords) {
    lenSum += s.getValue().size();
  }
  auto aveWordLength = lenSum / allWords.size();
  allWords = getSubsetRandom(allWords, maxWords);
  allWords = randomize(allWords);

  std::vector<Sequence> words(allWords.begin(), allWords.begin() + (allWords.size() - nQueries));
  std::vector<Sequence> qSeqs(allWords.begin() + (allWords.size() - nQueries), allWords.end());

  bool hflag = true;
  for (const auto& pivType : includePivotTypes) {
    for (const auto& partType : includePartTypes) {
      for (const auto& madi : cmt_madis) {
        collectSearchTest<Sequence, EditDistMetric, CMTree<Sequence, EditDistMetric>>(
            words, qSeqs, aveWordLength, radii, radiusMult, pivType, partType, madi, fileNamePrefix + "_ext", hflag,
            rangeFlag);
        if (hflag == true) {
          hflag = false;
        }
      }
      for (const auto& madi : spmt_madis) {
        collectSearchTest<Sequence, EditDistMetric, SPMTree<Sequence, EditDistMetric>>(
            words, qSeqs, aveWordLength, radii, radiusMult, pivType, partType, madi, fileNamePrefix + "_ext", hflag,
            rangeFlag);
        if (hflag == true) {
          hflag = false;
        }
      }
    }
  }
}

void collectCountSearchTestEDM(const std::string& fileNamePrefix, const std::string& inputDataFName) {
	std::map<unsigned int, unsigned int> nofPoints{ {100000,100 } };
	//These radii with radiusMult == false are actuall radii
	//std::vector<float> radii {  1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f. 7.0f};
	//For the SwissProt DB:
	std::vector<float> radii {1.0f,100.0f, 200.0f, 300.0f,400.0f,500.f,600.0f,700.0f,800.0f,900.0f, 1000.0f, 1100.0f,1200.0f, 1300.0f};
	bool radiusMult = false;

	//These radii with radiusMult==true are actually factors to the radius: rad = lengt * (1.0 - fac)
	//bool radiusMult = true;
	//std::vector<float> 	radii {0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};	

	std::set<PivotType> includePivotTypes{ PivotType::RAN}; 
	std::set<PartType>  includePartTypes{ PartType::BOM};

	std::vector<int> cmt_madis {1, 100};

	auto start = std::clock();
	Reader reader;
	auto allWords = reader.read(inputDataFName, 10000000);
	allWords = randomize(allWords);
	auto bTime = dTimeSeconds(start);
	std::cout << "rcollectCountSeatchEDM datagen time=" << bTime << std::endl;
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();

	bool hflag = true;
	std::vector<EuclidianPoint> points, qPoints;
	for (const auto& [np, nQueries] : nofPoints) {
		
		auto words = getSubsetRandom(allWords, np);
		auto qSeqs= getSubsetRandom(words, nQueries);

		for (const auto& pivType  : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					collectCountSearchTest<Sequence, EditDistMetric,CMTree<Sequence, EditDistMetric>>
					(words, qSeqs, aveWordLength, radii, radiusMult, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) {hflag = false;}
				}
			}
		}
	}
}


/*
	Use the search tree(s) to perform KNN search over the various combinations of the
	 search options (pivot types, radii, db sizes, etc). This function is used primarily to 
	 gather performance comparison statistics
*/

void nkSearchTestEDM(const std::string& fileNamePrefix, const std::string& inputDataFName) {
	unsigned int bTime, sTime;
	std::map<unsigned int, unsigned int> nofPoints{ {1000,100 }, {10000,100 }, {100000,100 } , {500000,100 } };//for dictionary
	//std::map<unsigned int, unsigned int> nofPoints{ {100000,100 } };
	//std::vector<unsigned int> maxResults{ 1,2,3, 4, 5,};//for dictionary
	std::vector<unsigned int> maxResults{ 1, 2,3,4, 5, 10};//for bio seqs

	//std::vector<float> lcmt_npivs {0.25, 0.5, 1.0, 2.0, 3.0};
	//std::vector<float> cmt_madis {1000, 10, 5, 1};

	std::vector<float> lcmt_npivs {1.0};
	std::vector<float> cmt_madis {1000, 1};

	auto start = std::clock();
	Reader reader;
	auto allWords = reader.read(inputDataFName, 900000);
	bTime = dTimeSeconds(start);
	std::cout << "KNNSearchTestEDM datagen time=" << bTime << std::endl;
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();

	//std::set<PivotType> includePivotTypes{ PivotType::RAN, PivotType::EXT};
	std::set<PivotType> includePivotTypes{ PivotType::RAN};
	std::set<PartType>  includePartTypes{ PartType::BOM };

	std::vector<float> radii;

	bool hflag = true;
	for (const auto& [np, nQueries] : nofPoints) {
		auto words = getSubsetRandom(allWords, np);
		auto qSeqs= getSubsetRandom(words, nQueries);
		for (const auto& pivType : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					nkSearchTest<Sequence, EditDistMetric,CMTree<Sequence, EditDistMetric>>
					(words, qSeqs, radii, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) { hflag = false; }
				}
			}
		}
	}}



/*
	Use the search tree(s) to perform KNN search. Query objects are not in
	the dataset being queried.
	The Small combinations version - used for exploratory purpose.
	Use an ifinite radius if radMultiple <0, otherwise the radius is
		seqeunce.size() * radMultiple
*/
void nkSearchTestEDM_EXQ(const std::string& fileNamePrefix, const std::string& inputDataFName) {
	const int maxWordsRead = 1000000;
	const int nQueries = 100;
	unsigned int bTime;
	std::vector<int> nofPoints{ 1000000};
	std::vector<unsigned int> maxResults{ 1, 2, 4, 6, 8, 10};
	//std::vector<float> radMults{0.01, 0.02, 0.03, 0.04, .05, 0.06, .08, 0.10};
	//std::vector<float> radMults{-1.0};// One negative for infinite radii
	std::vector<float> radMults;
	for (int i = 1; i<= 10; i++){
		radMults.push_back ( 0.01 * i);
	}
	
	//radMults.push_back (.015 );

	//for (int i = 2; i<=8 ; i++){
		//radMults.push_back ( 0.1 * i);
	//}
	
	auto start = std::clock();
	Reader reader;
	auto allWords = reader.read(inputDataFName, maxWordsRead);
	allWords = getSubsetRandom(allWords, maxWordsRead);
	allWords = randomize(allWords);

	bTime = dTimeSeconds(start);
	std::cout << "KNNSearchTestEDM datagen time=" << bTime << std::endl;
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();

	//ofstream fFile(fileNamePrefix + "_wordfreq_gcc.txt", ios::app);
	//printFrequency(allWords, fFile);

	std::set<PivotType> includePivotTypes {PivotType::RAN};
	std::set<PartType>  includePartTypes{PartType::BOM};
	
	std::vector<float> lcmt_npivs ;
	std::vector<int> cmt_madis {1,1000};
	std::vector<int> spmt_madis;
	std::vector<int> apmt_madis {0};
	
	//Below assumes words are already in random order
	std::vector<Sequence> words(allWords.begin(), allWords.begin() + (allWords.size() - nQueries));
    std::vector<Sequence> qSeqs(allWords.begin() + (allWords.size() - nQueries), allWords.end());

    for (const auto& radMult : radMults) {
      std::vector<float> radii;
	  float radMultOut = radMult;  //used for output purpose
      if (radMult > 0) {
        for (int i = 0; i < qSeqs.size(); i++) {
          radii.push_back(qSeqs[i].getValue().size() * radMult);
        }
      } else {
		radMultOut = std::numeric_limits<float>::max();
        radii.push_back(std::numeric_limits<float>::max());
      }

      bool hflag = true;
      for (const auto np : nofPoints) {
        auto dWords = getSubsetRandom(words, np);
        for (const auto& pivType : includePivotTypes) {
          for (const auto& partType : includePartTypes) {
            for (const auto& madi : cmt_madis) {
              nkSearchTest<Sequence, EditDistMetric, CMTree<Sequence, EditDistMetric>>(
                  dWords, qSeqs, radii, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix, hflag,
                  radMultOut);

              hflag = false;
            }
             for (const auto& madi : spmt_madis) {
              nkSearchTest<Sequence, EditDistMetric, SPMTree<Sequence, EditDistMetric>>(
                  dWords, qSeqs, radii, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix, hflag,
                  radMultOut);
              hflag = false;
            }
            for (const auto& madi : apmt_madis) {
              nkSearchTest<Sequence, EditDistMetric, APMTree<Sequence, EditDistMetric>>(
                  dWords, qSeqs, radii, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix, hflag,
                  radMultOut);
              hflag = false;
            }
          }
        }
	  }
    }
}

/*
	Use the search tree(s) to perform KNN search. Query objects are not in
	the dataset being queriedd
	The Small combinations version - used for exploratory purpose.
	Use an ifinite radius if radMultiple <0, otherwise the radius is
		seqeunce.size() * radMultiple
*/
void nkSearchTestEDM_EXQ_NP(const std::string& fileNamePrefix, std::vector<Sequence>& allWords, const int nQueries) {
	//std::vector<int> nofPoints{ 1000, 10000, 100000, 1000000};
	//std::vector<unsigned int> maxResults{1 , 2, 4, 6, 8, 10};
	std::vector<int> nofPoints{ 1000000 };//for dictionaries 
	std::vector<unsigned int> maxResults{1, 2 , 4, 6, 8,  10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000};
	
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();

	//ofstream fFile(fileNamePrefix + "_wordfreq_gcc.txt", ios::app);
	//printFrequency(allWords, fFile);

	std::set<PivotType> includePivotTypes {PivotType::RAN};
	std::set<PartType>  includePartTypes{PartType::BOM};

	//std::vector<float> lcmt_npivs{1,2, 4} ;
	std::vector<float> lcmt_npivs;
	std::vector<int> cmt_madis {1,1000};
	std::vector<int> apmt_madis {0};

			
	bool hflag = true; 
	
	for ( auto np : nofPoints) {
		if (np + nQueries > allWords.size()){
			np = allWords.size() - nQueries;
		}
		//Below assumes words are already in random order
		std::vector<Sequence> npWords = getSubsetRandom(allWords, np + nQueries);
		std::vector<Sequence> dWords(npWords.begin(), npWords.begin() + np);//the first np words from the ransom subset
		std::vector<Sequence> qSeqs(npWords.begin() + np, npWords.begin() + np + nQueries);// the remaining in npWords, which should

		std::vector<float> radii;
		for (int i = 0; i < qSeqs.size(); i++) {
       		radii.push_back(std::numeric_limits<float>::max());
     	}

		for (const auto& pivType : includePivotTypes) {
			for (const auto& partType : includePartTypes) {
				for (const auto& madi : cmt_madis){
					nkSearchTest<Sequence, EditDistMetric,CMTree<Sequence, EditDistMetric>>
					(dWords, qSeqs, radii, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) { hflag = false; }
				}
				
				for (const auto& madi : apmt_madis){
					nkSearchTest<Sequence, EditDistMetric,APMTree<Sequence, EditDistMetric>>
					(dWords, qSeqs, radii, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix, hflag);
					if (hflag == true) { hflag = false; }
				}
			}
		}
	}
}
void nkSearchTestEDM_EXQ_NP(const std::string& fileNamePrefix, const int nSeqs, const int nQueries = 100) {
  std::vector<Sequence> result;
  auto start = std::clock();
  for (int i = 1; i <= nSeqs + nQueries; i++) {
    std::string sId = std::to_string(i);
	unsigned int sz = randomUnsignedInteger(2, 50);
   // std::string sVal = CharacterSets::randomEnglish(sz);
	 std::string sVal = CharacterSets::randomGenomic(sz);
    result.push_back(Sequence(sId, sVal));
  }
  auto bTime = dTimeSeconds(start);
  std::cout << "KNNSearchTestEDM datagen time=" << bTime << std::endl;
  nkSearchTestEDM_EXQ_NP(fileNamePrefix, result, nQueries);
}

void nkSearchTestEDM_EXQ_NP(const std::string& fileNamePrefix, const std::string& inputDataFName, const int nQueries = 100) {
   int MAX_WORDS = 1000000;
  auto start = std::clock();
  Reader reader;
  auto allWords = reader.read(inputDataFName, MAX_WORDS);
  allWords = getSubsetRandom(allWords, MAX_WORDS);
  allWords = randomize(allWords);
  auto bTime = dTimeSeconds(start);
  std::cout << "KNNSearchTestEDM datagen time=" << bTime << std::endl;
  nkSearchTestEDM_EXQ_NP(fileNamePrefix, allWords, nQueries);
}



void nkRsSearchTestEDM(const std::string fileNamePrefix, const std::string inputDataFName) {
  unsigned int bTime, sTime;
  std::vector<unsigned int> nofPoints{1000000};
  int nQueries = 100;
  std::vector<unsigned int> maxResults{1, 5, 10, 20};

  auto start = std::clock();
  Reader reader;
  auto allWords = reader.read(inputDataFName, 1000000);
  bTime = dTimeSeconds(start);
  std::cout << "radiusSearchTestEDM datagen time=" << bTime << std::endl;
  auto lenSum = 0;
  for (const auto& s : allWords) {
    lenSum += s.getValue().size();
  }
  auto aveWordLength = lenSum / allWords.size();

  allWords = getSubsetRandom(allWords, 1000000);
  allWords = randomize(allWords);

  std::vector<Sequence> words(allWords.begin(), allWords.begin() + (allWords.size() - nQueries));
  std::vector<Sequence> qSeqs(allWords.begin() + (allWords.size() - nQueries), allWords.end());

  std::set<PivotType> includePivotTypes{PivotType::RAN};
  std::set<PartType> includePartTypes{PartType::BOM};
  std::vector<int> cmt_madis{1,1000}; //{1000};
  std::vector<int> spmt_madis{0};

  bool hflag = true;
  for (const auto& np : nofPoints) {
    auto dWords = getSubsetRandom(words, np);
    for (const auto& pivType : includePivotTypes) {
      for (const auto& partType : includePartTypes) {
        for (const auto& madi : cmt_madis) {
          nkRsSearchTest<Sequence, EditDistMetric, CMTree<Sequence, EditDistMetric>>(
              words, qSeqs, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix + "_EXT", hflag);
          hflag = false;
        }
        for (const auto& madi : spmt_madis) {
          nkRsSearchTest<Sequence, EditDistMetric, SPMTree<Sequence, EditDistMetric>>(
              words, qSeqs, aveWordLength, maxResults, pivType, partType, madi, fileNamePrefix + "_EXT", hflag);
          hflag = false;
        }
      }
    }
  }
}

void radiusSearchCompareEDM(unsigned int nPoints, unsigned int nQueries,
	PivotType pivT, PartType partT, float radius,
	const std::string& fileNamePrefix, bool csvh) {
	using namespace std;
	unsigned int bTime, sTime;

	unsigned int aveWordLength = 20;
	unsigned int sigma = aveWordLength / 4;

	std::vector<Sequence> words = generateWords(nPoints, aveWordLength, sigma);
	auto qSeqs = getSubsetRandom(words, nQueries);
	EditDistMetric met;
	auto start = std::clock();
	CMTree<Sequence, EditDistMetric> stree(words, met, pivT, partT, 1);
	BruteForceSearch<Sequence, EditDistMetric> stree2(words, met);


	bTime = dTimeSeconds(start);
	cout << "radiusSearchTestEDM btime=" << bTime << endl;

	std::set<unsigned int> depths;
	double radiusSum = stree.radiusSumAndDepths(depths);

	start = std::clock();
	unsigned int nFound = 0;
	unsigned int nqActual = 0;
	const unsigned int maxResults = 100000000;
	auto rad{ radius };
	unsigned int diffCount = 0;
	for ( auto& qp : qSeqs) {
		RadiusQuery<Sequence> rq(qp, radius, maxResults);
		stree.search(rq);
		nqActual++;
		nFound += rq.getNeighbors().size();

		RadiusQuery<Sequence> rq2(qp, rad, maxResults);
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
			*/
			/*
			auto nba = rq.getNeigbors();
			auto nbb = rq2.getNeigbors();
			for (const auto& n : nba) {
				std::cout << "d =" << n.distance << " id=" << n.id << std::endl;
			}
			for (const auto& n : nbb) {
				std::cout << "bd=" << n.distance << " id=" << n.id << std::endl;
			}
			*/
			/*
			std::cout << "-----" << std::endl;
			int i;
			std::cin >> i;
			*/
		}
		if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
		nQueries++;
		nFound += rq.getNeighbors().size();
	}
	sTime = dTimeSeconds(start);

	string runTimeFileName{ "treeResultsCompare.txt" };
	ofstream theFileA(runTimeFileName, ios::app);
	theFileA << ";;TEST EMSC: " << currentDateTime()
		<< ";;[CLasses]:" << endl
		<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
		<< ";;[CLasses # 2]:" << endl
		<< ";;[" << typeid(stree2).name() << "]" << endl
		<< "dbSize nQueries nfound  radiusSum Pivot Partition <nodesVisited> <numDistanceCalls> " << endl
		<< words.size() << "," << nqActual << "," << nFound << "," << radiusSum << ","
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (stree.getPerfStats().getNodesVisited() / nqActual) << ","
		<< (stree.getPerfStats().getDistanceCalls() / nqActual) << endl
		<< ";;btime, stime:" << endl
		<< bTime << "," << sTime << endl
		<< "Tree depths, min depth, max depth :" << depths.size() << ","
		<< *(std::min_element(depths.begin(), depths.end())) << "," << *(std::max_element(depths.begin(), depths.end())) << endl
		<< "DiscrepanciesTotals= " << diffCount << endl;

	theFileA.close();

	words.clear();
	qSeqs.clear();
}



void kNNCompareEDM(const std::string fileNamePrefix, const std::string inputDataFName,
	unsigned int nWords, unsigned int nQueries, unsigned int nofResults,
	PivotType pivT, PartType partT, bool csvh) {
	using namespace std;

	unsigned int bTime, sTime;

	auto start = std::clock();
	Reader reader;
	auto allWords = reader.read(inputDataFName, 200000);
	bTime = dTimeSeconds(start);
	std::cout << "KNNSearchTestEDM datagen time=" << bTime << std::endl;
	auto lenSum = 0;
	for (const auto& s : allWords) {
		lenSum += s.getValue().size();
	}
	auto aveWordLength = lenSum / allWords.size();

	auto words = getSubsetRandom(allWords, nWords);
	auto qSeqs = getSubsetRandom(words, nQueries);

	EditDistMetric met;
	start = std::clock();
	CMTree<Sequence, EditDistMetric> stree(words, met, pivT, partT);
	//APMTree<Sequence, EditDistMetric> stree(words, met, pivT, partT);
	BruteForceSearch<Sequence, EditDistMetric> stree2(words, met);
	bTime = dTimeSeconds(start);
	cout << "firstkSearchTest btime=" << bTime << endl;

	start = std::clock();
	unsigned int diffCount = 0;
	unsigned int nFound = 0;
	unsigned int nqActual = 0;
	auto rad{ std::numeric_limits<float>::max() };
	for ( auto& qp : qSeqs) {
		NearestKQuery<Sequence> nkQ(qp,  nofResults, rad);
		stree.search(nkQ);
		nqActual++;
		nFound += nkQ.getNeighbors().size();

		NearestKQuery<Sequence> nkQ2(qp,  nofResults, rad);
		stree2.search(nkQ2);

		if (!nkQ.hasSameNeighbors(nkQ2)) {
			diffCount++;
			auto nAs = nkQ.getNeighbors();
			for (const auto& nn : nAs) {
				std::cout << "na : " << nn << std::endl;
			}
			auto nBs = nkQ2.getNeighbors();
			for (const auto& nn : nBs) {
				std::cout << "nb : " << nn << std::endl;
			}
		}

		if ((nqActual % 1000) == 0) { cout << "Finished search i= " << nqActual << endl; }
	}
	sTime = dTimeSeconds(start);

	//AND the CSV file:
	ofstream csvFile(fileNamePrefix + ".csv", ios::app);
	if (csvh == true) {
		csvFile << ";;TEST KNN_COMPARE_EDM: " << currentDateTime()
			<< ";;[CLasses]:" << endl
			<< ";;[" << typeid(stree).name() << "]," << "[" << met << "]" << endl
			<< "DiffCnt,name,MDI,dim,dbSize,nQueries,maxResults,nfound,aveNfound,Pivot,Partition,nodesVisited,numDistCalls," << endl;
	}
	csvFile << diffCount << "," << stree.shortName() << "," << stree.getMADIorK() <<"," << aveWordLength<< "," 
		<<words.size() << "," << nqActual << "," << nofResults << ","
		<< nFound << "," << ((nFound * 1.0f) / nqActual) << ","
		<< stree.getPivotType() << "," << stree.getPartType() << ","
		<< (1.0f * stree.getPerfStats().getNodesVisited()) / (1.0f * nqActual) << ","
		<< static_cast<float>(stree.getPerfStats().getDistanceCalls()) / nqActual << ","
		<< (1.0f * stree.getPerfStats().getDistanceCalls()) / (1.0f * nFound) <<  endl;
	csvFile.close();
}

