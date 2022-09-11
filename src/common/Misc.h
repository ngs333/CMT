#ifndef MISC_H
#define MISC_H

#include <string>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <ratio>
#include <sstream>
#include <random>
#include <cmath>
#include <numeric>
#include <ostream>
#include <set>
#include <map>
#include <functional>
#include "Sequence.h"
#include "EuclidianPoint.h"

void error(std::string s)
{
	throw std::runtime_error(s);
}

void error(std::string s1, std::string s2)
{
	throw std::runtime_error(s1 + s2);
}


std::string currentDateTime() {
	using namespace std::chrono;
	//duration<int, std::ratio<60 * 60 * 24> > one_day(1);
	system_clock::time_point today = system_clock::now();
	time_t tt = system_clock::to_time_t(today);
	char str[26];
#ifdef _WIN32
	ctime_s(str, sizeof str, &tt);	
#elif __linux__
	ctime_r(&tt, str);
#else
	ctime_r(&tt, str);//TODO: Not really sure about this
#endif

	std::string string(str);
	return string;
}

/*
Time leapsed in milliseconds
*/
std::clock_t  dTimeSeconds(std::clock_t start) {
	std::clock_t rTime = 1000 * (std::clock() - start) / CLOCKS_PER_SEC;
	return rTime;
}

template <typename T>
std::string NumberToString(T Number)
{
	std::ostringstream ss;
	ss << Number;
	return ss.str();
}

bool endsWith(const std::string& mainStr, const std::string& toMatch)
{
	if (mainStr.size() >= toMatch.size() &&
		mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
		return true;
	else
		return false;
}


/*
	Generate and return vector to a std::array<PointType, DIM> from the passed in distribution object.
*/
template <unsigned int DIM, class PointType, class DistType>
std::vector<std::array<PointType,DIM>>
generatePoints(unsigned int nPoints, DistType& dist) {
	std::default_random_engine generator;
	std::vector<std::array<PointType, DIM>> points; 
	points.reserve(nPoints);
	std::array<PointType, DIM>  point;
	for (auto i = 0; i < nPoints; i++) {
		for (auto j = 0; j < DIM; j++) {
			point[j] = dist(generator);
		}
		points.push_back(point);
	}
	return points;
}


/*
     Generate nPoints euclidian points in the box of 1D ranges [start,width)
*/
template <unsigned int DIM, class PointType>
std::vector<EuclidianPoint>
genEuclidianPointsUniform(std::vector<int>::size_type nPoints, PointType start, PointType width) {
	//TODO: Below uses or _real_ or _int_ . Make compile time switch
	using DistribTypeUnif = std::uniform_real_distribution<EuclidianPointPointType>;
	DistribTypeUnif distribution(start, width);
	std::ostringstream idStream;
	std::vector<EuclidianPoint> EuPoints;
	auto points = generatePoints<DIM,PointType, DistribTypeUnif>(nPoints, distribution);
	unsigned int id = 0;
	EuPoints.reserve(nPoints);
	for (const auto& pt : points) {

		constexpr unsigned int nCopies = 1;
		for(int i = 0; i < nCopies; i++){  //carbon copies with different ids
			idStream << id;
			if(i != 0) idStream << "_" << i;
			EuPoints.push_back(EuclidianPoint(idStream.str(), pt));
			idStream.str(std::string());
			idStream.clear();
		}
		++id;
	}
	return EuPoints;
}
/*
	 Generate nPoints euclidian of equal spacing on a line. Usefull for
	 1D debugging
*/
template <unsigned int DIM, class PointType>
std::vector<EuclidianPoint>
genEuclidianPointsOnLine(const unsigned int nPoints){
	std::ostringstream idStream;
	std::vector<EuclidianPoint> EuPoints;
	EuPoints.reserve(nPoints);
	std::array<PointType, DIM> point;
	for (auto i = 0; i < nPoints; i++){
		for (auto j = 0; j < DIM; j++){
			point[j] = i;
		}
		idStream.str("");
		idStream.clear();
		idStream << "_" << i;
		EuPoints.push_back(EuclidianPoint(idStream.str(), point));
	}
	return EuPoints;
}

/**
	Get a vector of DIM distribution objects.
	The distribution function constructor is expected to take two values : mean and sigma.
**/
template<unsigned int DIM, class PointType, class DistType>
std::vector<DistType>
getDistributionVector(const EuclidianPoint& p, float sigma) {
	std::vector<DistType> distributions;
	for (int i = 0; i < DIM; i++) {
		auto mean = p[i];
		distributions.push_back(DistType(mean, sigma));
	}
	return distributions;
}


/***
	generate clustered euclidian points
	generateEuclidianPointsClustered<3,float,std::normal_distribution<float>
	
	N = p * Vol;  w = pow(N/p, 1/dim)
***/

template <unsigned int DIM, class PointType>
std::vector<EuclidianPoint>
generateEuclidianPointsClustered(unsigned int nClusters, unsigned int nPointsPer,
	PointType start, PointType width, float sigma) {
	std::vector<EuclidianPoint> points;  //resultant points.
	auto centerPoints = 
		genEuclidianPointsUniform<DIM, EuclidianPointPointType>(nClusters, start, width);

	using DistType = std::normal_distribution<float>;  //must be FLOAT or DOUBLE
	std::array<PointType, DIM> tempPoint;
	unsigned int id = 1;
	std::mt19937 gen(3001);
	std::ostringstream idStream;
	//for each center point, generate "nPointsPer" points:
	for (const auto& point : centerPoints) {
		auto distribVec = getDistributionVector<DIM, PointType, DistType>(point, sigma);
		for (auto i = 0; i < nPointsPer; i++) {
			idStream << id;
			for (auto j = 0; j < DIM; j++) {
				tempPoint[j] = (distribVec[j])(gen);
			}
			points.push_back(EuclidianPoint(idStream.str(), tempPoint));
			idStream.str(std::string());
			idStream.clear();
			++id;
		}
	}
	return points;
}

//Make a vector of sequences of equal lengths from a vector of seqeunces of 
//varying lengths.
std::vector<Sequence> resize(std::vector<Sequence> &isv,
	const double mean, const double sdev, const unsigned int bound = 3) {
	std::vector<Sequence> result;
	std::string rStr; //a result string

	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<> d{ mean,sdev };

	unsigned int namec = 0;

	//The lengh of the first one
	unsigned int wLength = std::round(d(gen));
	while ((wLength < (mean - bound * sdev)) ||
		(wLength > (mean + bound * sdev))) wLength = std::round(mean);

	//std::cout << "resize input" << std::endl;
	//for (auto elem : isv) { std::cout << elem << std::endl; }

	for (auto elem : isv) {
		std::string iStr = elem.getValue(); //an input string
		namec = 0;
		while (iStr.length() > 0) { //wile iStr has characters remaining
			auto needed = wLength - rStr.length();
			if (needed > 0) {
				if (iStr.length() <= needed) {
					rStr += iStr;
					iStr.clear();
				}
				else {
					rStr += iStr.substr(0, needed);
					iStr = iStr.substr(needed, iStr.length() - 1);
					//iStr.erase(iStr.begin(), needed);
				}
			}
			else {
				// make a unique new name for the new sequence.
				//assuming the input element names are also unique
				std::stringstream ss;
				ss << namec;
				std::string newName = elem.getId() + "_" + ss.str();
				result.push_back(Sequence(newName, rStr));
				namec++;
				rStr.clear();
				wLength = std::round(d(gen));
				while ((wLength < (mean - bound * sdev)) ||
					(wLength > (mean + bound * sdev))) wLength = std::round(mean);
			}
		}

	}

	if (rStr.length() > 0) {//last wasnt saved.
		std::stringstream ss;
		ss << namec;
		std::string newName = isv.back().getId() + "_" + ss.str();
		result.push_back(Sequence(newName, rStr));
		rStr.clear();
	}

	//std::cout << "resize ouput" << std::endl;
	//for (auto elem : result) {std::cout << elem << std::endl;}

	return result;
}

/** Split one sequence into a vector of subseqences of roughly equal lengths.
The two adjacent subsequences will have overlapps (tail to head). SUbsequences
will also be given IDs derived from the parents ID.

mL is the length of the middle part of the resulting sequences.
oL is the overlap length or length of the tail and the head  for
most resulting seqeunces.

TODO: can the splitting be done by score?
**/
std::vector<Sequence> splitSequenceWO(Sequence& s, unsigned int mLen, unsigned int oLen) {
	auto iStr = s.getValue(); //A copy of the input sequence value;
	std::string rStr;       //a result string
	std::vector<Sequence> result;
	int mL(mLen);
	int oL(oLen);
	int namec = 0;
	int cStart = 0;
	int cEnd = 0;
	int remain = iStr.length();
	int mStart = 0;
	int mEnd = oL;
	while (remain > 0) {
		//Advance the middle part by mL.
		mStart = mEnd;
		mEnd += mL;
		//overlap to the left - with check:
		cStart = std::max(0, mStart - oL);
		//overlap to the right - with check
		cEnd = std::min(mEnd + oL, (int)iStr.length());
		remain = iStr.length() - cEnd;
		if ((cEnd == iStr.length()) && ((cEnd - cStart) < (mL + oL))) {
			//if the last one is shorter than desired, start even further back
			cStart = std::max(0, cEnd - (mL + oL));
		}

		rStr = iStr.substr(cStart, cEnd - cStart);

		//Calc the name of the new substring
		std::stringstream ss;
		ss << namec;
		std::string newName = s.getId() + "_" + ss.str();
		Sequence ns(newName, rStr);
		result.push_back(ns);
		//rStr.clear();
		namec++;
	}
	return  result;
}

std::vector<Sequence> splitSequenceWO(std::vector<Sequence>& iSeqs, unsigned int mL, unsigned int oL) {
	std::vector<Sequence> result;
	for (auto s : iSeqs) {
		auto sv = splitSequenceWO(s, mL, oL);
		//Add the the new vector of subsequences to the result vector:
		//std::move(pts.begin(), pts.end(), std::inserter(result, result.end()));
		//result.insert(result.end(), std::make_move_iterator(pts.begin()),
			//std::make_move_iterator(pts.end()));
		for (auto seq : sv) {
			result.push_back(seq);
		}
		sv.clear();
	}
	return result;
}


void testSplitWOverlap(Sequence& s) {
	using std::endl;
	auto result = splitSequenceWO(s, 4, 2);
	std::cout << "for seq:" << s << " " << s.getValue() << std::endl;
	for (auto e : result) {
		std::cout << e << " " << e.getValue() << std::endl;
	}
}
void testSplitWOverlap() {
	Sequence ts2{ "SA", "a" };
	testSplitWOverlap(ts2);
	Sequence ts3{ "SA", "abcdef" };
	testSplitWOverlap(ts3);
	Sequence ts4{ "SA", "abcdefgh" };
	testSplitWOverlap(ts4);
	
	Sequence ts5{ "SA", "abcdefghi" };
	testSplitWOverlap(ts5);
	
	Sequence ts{ "SA", "01234567890123456789" };
	testSplitWOverlap(ts);
	Sequence ts6{ "SA", "abcdefghijklmnopqr" };
	testSplitWOverlap(ts6);
	
}
void checkSequenceVector(std::vector<Sequence> & seqs) {
	for (auto s : seqs) {
		if (s.getId() == "abcd") {
			std::cout << s << std::endl;
		}
		if (s.getValue() == "abcd") {
			std::cout << s << std::endl;
		}
	}
}

//Seq Stats
class SeqVecStats {
public:
	unsigned int size;
	unsigned int minL;
	unsigned int maxL;
	double   meanL;
	double   sdevL;
	double   minS;
	double   maxS;
	double   meanS;
	double   sdevS;
	friend std::ostream& operator<<(std::ostream& os, const SeqVecStats& s) {
		os << "SVS size=" << s.size << std::endl
			<< " lengths:[" << s.minL << "," << s.maxL << "];( " << s.meanL << "," << s.sdevL << ")"
			<< std::endl
			<< " scores:[" << s.minS << "," << s.maxS << "];( " << s.meanS << "," << s.sdevS << ")";
		return os;
	}
};

unsigned int addLength(unsigned int v, Sequence& s) {
	return v + s.getValue().length();
}
//calculate the mean and standard deviation of a sequence vector.
SeqVecStats 
calcSeqVecStats(std::vector<Sequence> & v, std::vector<double> vss) {

	SeqVecStats stat{ 0 };
	if (v.empty()) {
		return stat;
	}
	stat.size = v.size();
	stat.minL = v.begin()->getValue().length();
	stat.maxL = stat.minL;
	double sumL = std::accumulate(v.begin(), v.end(), 0, addLength);
	stat.meanL = sumL / stat.size;
	double product = 0;
	for (auto & e : v) {
		auto l = e.getValue().length();
		if (l < stat.minL) stat.minL = l;
		if (l > stat.maxL) stat.maxL = l;
		product += (l - stat.meanL) * (l - stat.meanL);
	}
	stat.sdevL = sqrt(product / (stat.size - 1));

	if (vss.empty()) {
		return stat;
	}
	stat.minS = *(vss.begin());
	stat.maxS = stat.minL;
	double sumS = std::accumulate(vss.begin(), vss.end(), 0.0);
	stat.meanS = sumS / stat.size;
	product = 0.0;
	for (auto & ss : vss) {
		if (ss < stat.minS) stat.minS = ss;
		if (ss > stat.maxS) stat.maxS = ss;
		product += (ss - stat.meanS) * (ss - stat.meanS);
	}
	stat.sdevS = sqrt(product / (stat.size - 1));

	return stat;
}

double pctDifference(double a, double b) {
	double p = (a - b) * 100.0 / (a + b);
	if (p < 0.0) p = -p;
	return p;
}

unsigned int
randomUnsignedInteger(const unsigned int min, const unsigned int max) {
	//std::random_device rd;  //Will be used to obtain a seed for the random number engine
	//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen(3001); //Standard mersenne_twister_engine seeded with 3001
	std::uniform_int_distribution<> dis(min, max);
	return dis(gen);
}

std::vector<int>
randomIntegerSet(const unsigned int min, const unsigned int max, const unsigned int nSamples, bool ordered=true) {
	//std::set<int> intSet;
	std::vector<int> intSet;

	//std::random_device rd;  //Will be used to obtain a seed for the random number engine
	//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen(3001); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(min, max);

	const unsigned int MaxMMin = max - min + 1;

	if (nSamples > (MaxMMin)) {
		error("randomIntegerSet", "nSamples > (max - min + 1)");
	}

	std::vector<unsigned int> vints;
	vints.reserve(MaxMMin);
	for (unsigned int i = min; i <= max; i++) {
		vints.push_back(i);
	}

	for (unsigned int i = min; i <= max; i++) {
		auto j = dis(gen);
		std::swap(vints[i], vints[j]);
	}

	for (unsigned int i = 0; i < nSamples; i++) {
		//intSet.insert(vints[i]);
		intSet.push_back(vints[i]);
	}

	if(ordered == true){
		std::sort(intSet.begin(), intSet.end());
	}

	return intSet;
}

/*
	Select and return the subset of points that are interior to
	the axis-aling bounding bod defined by [boxO,boxP] (the boxes origin and far points)
*/
template<class PT, unsigned int DIM>
std::vector<EuclidianPoint>
getInteriorPoints(std::vector<EuclidianPoint>& points, std::array<PT, DIM> const& boxO, std::array<PT, DIM>  const& boxP) {
	std::vector<EuclidianPoint> iPoints;
	for (auto& pt : points) {
		if (pt.isInsideBox(boxO, boxP)) {
			iPoints.push_back(pt);
		}
	}
	return iPoints;
}

	
/*
	Return an order randomized copy of the input vrector
	*/
template<class T>
std::vector<T>
randomize(std::vector<T>& pts) {
	std::vector<T> result;
	auto rInts = randomIntegerSet(0, pts.size() - 1, pts.size(), false);
	for (const auto& ri : rInts) {
		result.push_back(pts[ri]);
	}
	return result;
}

/*
	Return a vector that contains a random subset of size nSubset of the pts object vector.
*/
template <class T>
std::vector<T>
getSubsetRandom(std::vector<T>& pts, const unsigned int nSubset, bool randomizeSmall= false) {
	std::vector<T> qPoints;
	if (pts.size() <= nSubset && randomizeSmall == false) {
			qPoints = pts;
			return qPoints;
	}
	auto nInts = nSubset;
	if (pts.size() <= nSubset) nInts = pts.size();
	auto rInts = randomIntegerSet(0, pts.size() - 1, nInts);
	for (const auto& ri : rInts) {
		qPoints.push_back(pts[ri]);
	}
	return qPoints;
}

/*
	Select and return nPoints from the interior of the point cube sorrounding
	the points vector "points".

	For NSeqs < pts.size(), the returned subset is chosen by selecting entries from
	pts that are roughly equaly spaced (in order within pts) and with the interior.

*/
template<class PT, unsigned int DIM>
std::vector<EuclidianPoint>
getInteriorQueryPoints(std::vector<EuclidianPoint>& pts, unsigned int nQueries, std::array<PT, DIM> const& o, std::array<PT, DIM>  const& p) {
	std::vector<EuclidianPoint> qPoints;
	qPoints.reserve(nQueries);
	if (pts.size() <= nQueries) {
		qPoints = pts;
		return qPoints;
	}
	unsigned int skip = std::floor(pts.size() / (1.0 * nQueries));
	if (skip < 1) skip = 1;
	for (int k = 0; k < skip;  ++k) {
		if (qPoints.size() >= nQueries) break;
		for (auto j = k; j < pts.size(); j += skip) {
			if (qPoints.size() >= nQueries) break;
			if (pts[j].isInsideBox(o, p)) {
				qPoints.push_back(pts[j]);
			}
		}
	}
	//TODO: if (qPoints.size() < nQueries) ???
	return qPoints;
}

/**
	Generate a map of DB of pairs of [DB size, Number of  DB queries].
	The purpose of this is to help automate test results.
	Normally base should be 2 or 10.
	All pairs with DB sizes less than maxNQ will get the number of queries to be
	half of the DB size (since interested in finding non-surface points). All others
	will get for the number of queries the first number that exceeded maxNQ
**/
std::map<unsigned int, unsigned int>
getTestSizes(unsigned int base,  unsigned int startP, unsigned int endP, unsigned int skipP, unsigned int maxNQ = 10000) {
	std::map<unsigned int, unsigned int> sizes; //<target size of DB; not of queries desired
	unsigned int delta = 0;
	if (base == 2) delta = 1;
	auto lastMaxNQ = 0;
	for (auto p = startP; p <= endP; p+= skipP) {
		auto np = std::pow(base, p) - delta;
		auto nq = np/2;
		if (nq > maxNQ) {
			if (lastMaxNQ == 0) {
				lastMaxNQ =  maxNQ;
			}
			if (lastMaxNQ != 0) {
				nq = lastMaxNQ;
			}
		}
		sizes[np] = nq;
	}
	return sizes;
}


template <class T, size_t ROW, size_t COL>
using NativeMatrix = T[ROW][COL];


#endif
