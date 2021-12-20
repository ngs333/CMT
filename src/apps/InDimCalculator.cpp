#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <chrono>

#include "Misc.h"
#include "Metric.h"
#include "SearchCommon.h"
#include "Reader.h"
#include "GPCharacters.hpp"

#ifdef __WIN32
unsigned int DistanceIntervalSM<float>::used = 0;
std::vector<DistanceIntervalSM<float>> DistanceIntervalSM<float>::data;
#endif

string bfn{"/home/mzuniga/data"};
string sprotFileName{bfn + "/UniProt/uniprot_sprot.fasta"};
string tremblFileName{bfn + "/UniProt/TREMBL/uniprot_trembl_500000.fasta"};
string uproMammalsFileName{bfn + "/UniProt/TREMBL/mammals/uniprot_trembl_mammals_500000.fasta"};
string influenzaFileName{bfn + "/HA/influenza.fna"};
string haEncodingFileName{bfn + "/HA/enc_09/ha_enc09.fasta"};
string hivDirName{bfn + "/HIV1"};
string dictionFileName(bfn + "/text/words.txt");


template <class T>
int  findCenter(std::vector<T> & objects,  Metric<T>& met, unsigned int NSamples) {
	int center{ objects.size()}; //Iterator (pointer) to the object at the center
	double minDist{ std::numeric_limits<double>::max() };

	std::vector<double> dv;
	dv.reserve(objects.size());
	auto rintSet = randomIntegerSet(0, objects.size() - 1, NSamples);
	for (auto rit = rintSet.begin(); rit != rintSet.end(); ++rit) {
		for(int i = 0; i< objects.size(); i++){
			dv [i] = met.distance(objects[*rit], objects[i]);
		}
		auto farDist = * std::max_element(dv.begin(), dv.end());
		if (farDist <= minDist) { //a new min farthest neighbor
			minDist = farDist;
			center= *rit;
		}
	}
	return  center;
}
template<class T>
double intrinsicDimAlt(std::vector<T>&  objects, Metric<T>& metric, unsigned int NSamples= 100){
	if (NSamples > objects.size())
		NSamples = objects.size();

	auto center = findCenter<T>(objects, metric, NSamples);

	std::vector<double> dv;
	for(const auto& obj : objects){
		dv.push_back( metric.distance(obj, objects[center]) );
	}
	unsigned int size = dv.size();
    double mean = std::accumulate(dv.begin(), dv.end(), 0.0) / size;

	double variance = 0;
	for(const auto & val : dv){
		variance += (val - mean) * (val - mean);
	}
	variance /= (size -1);

	return mean * mean / (2 * variance);
}

template<class T>
double intrinsicDim(std::vector<T>&  objects, Metric<T>& metric, unsigned int NSamples= 100){
	if (NSamples > objects.size())
		NSamples = objects.size();

	auto intsA= randomIntegerSet(0, objects.size() -1, NSamples);
	auto intsB= randomIntegerSet(0, objects.size() -1, NSamples);

	std::vector<double> dv;
	for (const auto & ia: intsA){
		for (const auto & ib : intsB){	
			dv.push_back( metric.distance(objects[ia], objects[ib]) );
		}
	}
	unsigned int size = dv.size();
    double mean = std::accumulate(dv.begin(), dv.end(), 0.0) / size;

	double variance = 0;
	for(const auto & val : dv){
		variance += (val - mean) * (val - mean);
	}
	variance = variance / (size -1);

	return mean * mean / (2 * variance);

}
double FileDataInDim(const std::string fname, const unsigned int size){
	EditDistMetric metric;
	Reader reader;
	auto words = reader.read(fname, size);
	double inDim = intrinsicDim<Sequence>(words, metric, 400);
	return inDim;
}
/**
 * Intrinsic calculator application.
 **/
int main(int argc, char *argv[])
{
	using namespace std;

	string str = CharacterSets::randomEnglish(26);

	cout << "string= " << str << endl;
	


	std::array<EuclidianPointPointType, EuclidianPointDim> pa;
	EuclidianMetric met;
	float density = 1.0;
	auto width = ceil(pow((1000 / density), 1.0 / EuclidianPointDim));
	auto  dbPoints = genEuclidianPointsUniform<EuclidianPointDim, EuclidianPointPointType>(1000, 0.0f, width);

	double inDim = intrinsicDim<EuclidianPoint>(dbPoints, met);
	std::cout << "For DIM="<<EuclidianPointDim<<" Euclidian indDim=" << inDim << std::endl;

	inDim = FileDataInDim(hivDirName, 10000);
	std::cout << "For "<< hivDirName<< ", indDim=" << inDim << std::endl;

	inDim = FileDataInDim(haEncodingFileName, 10000);
	std::cout << "For "<< haEncodingFileName<< ", indDim=" << inDim << std::endl;

	inDim = FileDataInDim(sprotFileName, 10000);
	std::cout << "For "<< sprotFileName<< ", indDim=" << inDim << std::endl;
}


