#ifndef SEARCH_COMMON_H
#define SEARCH_COMMON_H

#include <vector>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_set>

#include "Metric.h"


/***
	File of miscellaneous strunctures and algorithms common to various search trees
***/

// #define USE_BOOST_ARCHIVE   //If you define USE_BOOST_ARCHIVE, then you will need the boost archive libs


namespace PerfStatsNS {
	long long int distanceCalls;
	long long int nodesVisited;
	long long int maxPDCtr[10000];
};

/**
	OM => Object-Median
	BOM => Balanced (tree) OM
	DMR => Distance mid-range
	EXT => extremum
	CENT => center
*/

enum class PivotType { RAN, EXT, RXT, MIN, MED, CENT, N4, TSP};//Random, Extreme, MINIMUM, median, center, N/4
enum class PartType {OM, BOM, DMR, EXT, MIN  };

static std::map<PivotType, std::string> pivotTypeMap{ {PivotType::RAN, "RAN"},{PivotType::EXT, "EXT"},{PivotType::RXT, "RXT"},
	{PivotType::MIN,"MIN"},{PivotType::MED, "MED"}, {PivotType::CENT, "CENT"},{PivotType::N4, "N4"} };

static std::map<PartType, std::string> partTypeMap{ {PartType::OM, "OM"}, {PartType::BOM, "BOM"}, {PartType::DMR, "DMR"},
	{PartType::EXT, "EXT"},{PartType::MIN,"MIN"} };

static std::map<std::string, PivotType> pivotTypeRevMap{ {"RAN", PivotType::RAN},{"EXT", PivotType::EXT},{"RXT", PivotType::RXT},
	{"MIN", PivotType::MIN},{"MED", PivotType::MED}, {"CENT", PivotType::CENT},{"N4", PivotType::N4} };

static std::map<std::string, PartType> partTypeRevMap{ {"OM", PartType::OM}, {"BOM", PartType::BOM}, {"DMR", PartType::DMR,},
	{"EXT", PartType::EXT,}, {"MIN", PartType::MIN} };


std::ostream& operator<<(std::ostream& os, PivotType p) {
	os << pivotTypeMap.at(p);
	return os;
}
std::ostream& operator<<(std::ostream& os, PartType p) {
	os << partTypeMap.at(p);
	return os;
}

//template <class T, size_t ROW, size_t COL >
//using Matrix = std::array<std::array<T, COL>, ROW>;


/*
	constant k times the balanced tree heigh form tree with numNodes nodes.
	This function is sometimes used for the LCMT number of pivots decision, 
	and for which k is USUALLY element of { 0.25, 0.5, 1, 1.5,2, 2.5. 3 }
*/
unsigned int kxBalancedTreeHeight( double k, unsigned int numNodes){ 
	return std::ceil(k * std::log2(numNodes)) ;
	}

/*
	return a random iterator to one in the set [begin,end) of iterators.
*/
template <class TPtr>
inline TPtr findRandomIter(const TPtr begin, const TPtr end)  {
	auto randIter = begin + randomUnsignedInteger(0, end - begin - 1);
	return  randIter;
}

/*
	From the set of points pointed to by the iterator set [begin,end), find and return the 
	largest distance from the input iterator (sIter) object.
*/
template <class T, class TPtr>
double findFurthestDist(const TPtr begin, const TPtr end, const TPtr sIter,  Metric<T>& mp) {
	auto farDist = 0.0;
	auto farthest = sIter;
	for (auto iter = begin; iter != end; ++iter) {
		if (&(*sIter) == &(*iter)) {
			continue;
		}
		double dist = mp.distance((*sIter)->object, (*iter)->object);
		if (dist > farDist) {
			farDist = dist;
			farthest = iter;
		}
	}
	return farDist;
}

/*
	Find the center object in the set [begin,end).
	essintially by brue force search by picking the object out of NSamples that has the
	minimum furthest distance to the NPoints.
	Note the computational complexity is NSamples * NPoints, or quadratic for NSamples ~ Npoints.
	NSamples should be <= NPoints.

*/
template <class T, class TPtr>
TPtr findCenter(const TPtr begin, const TPtr end, Metric<T>& mp, unsigned int NSamples) {
	auto NPoints = end - begin;
	if (NPoints == 1 || NPoints == 2) {
		return begin;
	}
	if (NSamples > NPoints) {
		NSamples = NPoints;
	}

	auto centerIter{ end }; //Iterator (pointer) to the object at the center
	double minDist{ std::numeric_limits<double>::max() };

	auto rintSet = randomIntegerSet(0, NPoints - 1, NSamples);
	for (auto rit = rintSet.begin(); rit != rintSet.end(); ++rit) {
		auto sIter = begin + *rit; //Iterator to the sample object
		auto farDist = findFurthestDist<T,TPtr>(begin, end, sIter, mp);
		if (farDist <= minDist) { //a new min farthest neighbor
			minDist = farDist;
			centerIter = sIter;
		}
	}
	return  centerIter;
}


/*
	Find (and return a pointer to) an extrema for the set of object in the container
	of nodes in range [begin, end).
	TODO: 1. Does the benefits of calling findExtrema() a seecond time outweeigh the cost 
			(in the context of tree build and tree search)
		  2. Choose between returning pIter and qIter dependeing of "left-of" and/or "smaller" concepts.

*/
template <typename T, typename TPtr>
TPtr  findExtrema(const TPtr begin, const TPtr end, Metric<T>& mp) {
	TPtr  extB;
	return findExtrema(begin, end, mp , extB);
}

template <typename T, typename TIter>
TIter findExtrema(const TIter begin, const TIter end, Metric<T>& mp,  TIter& extB) {
	auto NPoints = end - begin;
	if (NPoints == 1 || NPoints == 2) {
		return begin;
	}
	TIter rIter = findRandomIter(begin, end);
	auto pIter = findExtrema(begin, end, rIter, mp);
	auto qIter = findExtrema(begin, end, pIter, mp);
	extB = pIter;//the other extrema
	return qIter;
}

template <typename T, typename TIter>
TIter findRandomExtrema(const TIter begin, const TIter end, Metric<T>& mp) {
	auto NPoints = end - begin;
	if (NPoints == 1 || NPoints == 2) {
		return begin;
	}
	TIter rIter = findRandomIter(begin, end);
	auto pIter = findExtrema(begin, end, rIter, mp);
	return pIter;
}

/*
Find and return an iterator the the furthest object from the one pointed to by peItr.
The set of points are the points in the container of nodes in range [begin, end).
It is assumded distances are greater than or equal to zero.]
TODO: Refactor the various find extrema to have one commen inner function.
*/
template <typename T, typename TPtr>
TPtr  findExtrema(const TPtr begin, const TPtr end, const TPtr peItr, Metric<T>& mp) {
	if (end - begin == 1) {
		return begin;
	}
	auto farDist = 0.0;
	auto qIter = begin;
	for (auto sIter = begin; sIter != end; sIter++) {
		if (&(*sIter) == &(*peItr))  continue;
		auto dist = mp.distance((*sIter)->object, (*peItr)->object);
		if (dist > farDist) {
			qIter = sIter;
			farDist = dist;
		}
	}
	return qIter;
}



/*
template <typename T>
std::vector<T>  
findExtrema(std::vector<T> & points,  unsigned int NExt, const Metric<T>& mp) {
	std::vector<T> extrema;
	std::set<std::string> usedIdSet;
	auto intSet = randomIntegerSet(0, points.size() - 1, NExt);
	for(const auto & ri : intSet ){
		auto eIter = points.begin() + ri;
		//std::cout << (*eIter).getId() <<std::endl;
		eIter = findExtrema(points.begin(), points.end(), eIter, usedIdSet,  [&](T& a, T& b) {return mp.distance(a, b); });
		extrema.push_back(*eIter);
		usedIdSet.insert((*eIter).getId());
	}
	return extrema;
}
*/
/*
template <class T, class RIter, class M>
unsigned int findFurthest(std::vector<T>& objs,  const RIter refI  ,  M & met){
	auto farDist = 0.0;
	unsigned int  farthest = objs.size();
	for (unsigned int i = 0; i < objs.size(); ++i){
		double dist = met.distance(*refI, objs[i]);
		if (dist > farDist) {
			farDist = dist;
			farthest = i;
		}
	}
	return farthest;
}
*/
/*
	Return a vector that contains a set NS of extrema objects.
	The extrema objects are first determined by choosing NS
	random objects, then finding the NSE objects furthest from them.
	if NSE < NS, then NS - NSE objects from the original random set are 
	the the final result set.
*/
/*
template <class T, class M>
std::vector<T>
getSubsetRandomExtrema(std::vector<T> &pts, const unsigned int NS,  M& mp){
	//Simpler with std::unordered_set<T, MyHasher, CompareById> sSet; ?
	std::unordered_set<std::string> idSet;
	//unsigned int ns2 = std::min(4 * NS, (unsigned int) pts.size());
	auto randoms = getSubsetRandom(pts,NS);
	std::vector<T> subset;
		
	for (auto sIter = randoms.begin(); sIter != randoms.end(); ++sIter){
		unsigned int furthest = findFurthest(pts, sIter, mp);
		if(idSet.find( pts[furthest].getId() ) == idSet.end()){
			idSet.insert(pts[furthest].getId());
			subset.push_back(pts[furthest]);
		}
		if(subset.size() == NS)break;
	}
	for (auto sIter = randoms.begin(); sIter != randoms.end(); ++sIter){
		if(subset.size() == NS)break;
		if(idSet.find((*sIter).getId() ) == idSet.end()){
			idSet.insert((*sIter).getId());
			subset.push_back(*sIter);
		}
		
	}
	return subset;
}
*/
/*
	Return the vector that contains the min-max set of pivots as per mico 1994 paper.
	TODO: This is a work in progress. 
*/
template <class T, class M>
std::vector<T>
getSubsetMaxMin(std::vector<T> &pts, const unsigned int NS,  M& mp){
	std::unordered_set<std::string> pivIds;
	std::vector<T> subset;

	unsigned int fi =  randomUnsignedInteger(0, pts.size() -1);
	subset.push_back( pts [ fi ] );
	pivIds.insert(pts[fi].getId());
	T* optr = & (pts[fi]);

	while(subset.size() < NS){
		double maxAcc = 0;
			//Find the object that is furthest away from all the existing pivots.
			//Furthest means having the largest sum of distances to the existing pivots.
			for( auto & obj : pts ){
				double acc = 0.0;	
				for (const auto & piv : subset){
					if(pivIds.find(obj.getId()) == pivIds.end()){
						//I.e. if object is not yet a pivot:
						acc += mp.distance(obj, piv);
					}
				}
				if(acc >= maxAcc){
					maxAcc = acc;
					optr  = &obj;
				}
			}
			subset.push_back( *optr );
	}
	return subset;
}

/*
	Return the vector that contains a min-max set of pivots as per S. Brinn 1995 paper
	TODO: This is a work in progress. 
*/
template <class T, class M>
std::vector<T>
getSubsetRanMaxMin(std::vector<T> &pts, const unsigned int NS,  M& mp){
	std::unordered_set<std::string> pivIds;
	std::vector<T> pivots;

	//Get a set of candidate pivots about three times the size
	// of the target pivot sizes
	auto candidates = getSubsetRandom(pts, 3 * NS);
	
	pivots = getSubsetRandom(candidates, std::ceil(NS/2));
	for (const auto & piv : pivots){
		pivIds.insert(piv.getId());	
	}

	while (pivots.size() < NS){
		//Find the next candidate obj that is furthest away from all the existing pivots
		auto maxDist = 0.0;
		T* optr = nullptr; //ptr to candidate object furthest from the pivots
		for (auto & cObj : candidates){
			if (pivIds.find(cObj.getId()) == pivIds.end()){
				//std::cout << "L:" << cObj.getValue()<<std::endl;
				//For the current candidate object, find its distance to the nearest pivot
				auto minDist = std::numeric_limits<float>::max();
				for (const auto & piv : pivots){
					auto dist = mp.distance(cObj, piv);
					if (dist <= minDist){
						minDist = dist;
					}
				}
				if (minDist >= maxDist){
					maxDist = minDist;
					optr = &cObj;
				}
			}
		}
		//std::cout<<"md="<<maxDist<<std::endl;
		pivots.push_back(*optr);
		pivIds.insert((*optr).getId());
	}
	return pivots;
}

/*
Find and return an iterator the 
the object the is "most different than ancestors",
"most different than" possibly defined by:
a) The sum of distances to each ancestors is greatest
b) The distance to nearest ancestor is the greatest.
The set of points are the points in the container of nodes in range [begin, end).
It is assumded distances are greater than or equal to zero.
*/
template <typename T, typename TPtr>
TPtr  furthestFromAncestors(const TPtr begin, const TPtr end, std::vector<T*>& ancestors, Metric<T>& mp) {
	if (end - begin == 1) {
		return begin;
	}
	auto farDist = 0.0;
	auto qIter = begin;
	for (auto sIter = begin; sIter != end; sIter++) {
		auto dist = 1.0e10; //sum of distances;
		for (const auto& an : ancestors) {
			 auto temp = mp.distance((*sIter)->object, *an);
			 if (temp < dist) dist = temp;
		}
		if (dist > farDist) {
			qIter = sIter;
			farDist = dist;
		}
	}
	return qIter;
}


/*
Select a pivot based on the strategy choosen via argument pivotType.
*/
template <class T, class M, class NodeItr, class Comparator>
inline NodeItr selectPivot(const NodeItr begin, const NodeItr median,
		const NodeItr end, PivotType pivotType, Comparator& compare, M& mp) {
	auto pivotItr = end;
	auto mmIter = begin + (median + 1 - begin) / 2;
	switch (pivotType) {
	case PivotType::RAN:
		pivotItr = findRandomIter<NodeItr>(begin, end);
		break;
	case PivotType::EXT:
		pivotItr = findExtrema<T, NodeItr>(begin, end, mp);
		break;
	case PivotType::RXT:
		pivotItr = findRandomExtrema<T, NodeItr>(begin, end, mp);
		break;
	case PivotType::MIN:
		pivotItr = std::min_element(begin, end, compare);
		break;
	case PivotType::MED:
		std::nth_element(begin, median, end, compare);
		pivotItr = median;
		break;
	case PivotType::CENT:
		pivotItr = findCenter<T, NodeItr>(begin, end, mp, 10);
		break;
	case PivotType::N4:   //Pivot = median of the median (by SS) on the lef
		std::nth_element(begin, mmIter, median + 1, compare);
		pivotItr = mmIter;
		break;
	default:
		error("selectPivot function default");
	}
	return pivotItr;
}



// Collect the Tree the nodes per level. Nodes on the LHS are visited 1st.
template <typename T>
void collectTreeNodes(T ndPtr, std::map<int, std::vector<T> >& map, int level){
	if (ndPtr){
		++level;
		collectTreeNodes(ndPtr->left, map, level);
		collectTreeNodes(ndPtr->right, map, level);
		map[level].push_back(ndPtr);
	}
}

/***
	Print in level order a binary tree.
	ndPtr is the root node pointer.
	Each node type is assumed to have the ostream<< operator defined.

	*** WARNING: Trees of more than 10 levels are not printed!
***/

template <typename T>
void printBinaryTree(T ndPtr, std::ostream& os)
{
	//In map levels, each level of the tree will have a vector of nodes.
	//The map is an ordered map and should be traversable by level value.
    std::map<int, std::vector<T>> levels;
	collectTreeNodes<T>(ndPtr, levels, 0);

	auto maxLevel = levels.size();
	if (maxLevel > 10) { exit(-1); }  // Dont even try to print huge trees
		
	for (auto& ln : levels) {
		auto levelNo  = ln.first;
		std::string str(10 * (maxLevel - levelNo), ' ');// blank chars printed before first node of level!
		os << "\nLevel " << levelNo << ": " << std::endl << str;
		for (auto& nd : ln.second) {  //for each node in the vector of the level
			os << *nd << "  ";
		}
	}
	os << std::endl;
}



template < typename T>
void collectNodesInMap(T ndPtr, std::map<std::string, T>& nmap) {
	if (ndPtr == nullptr)
		return;
	nmap[ndPtr->object.getId()] = ndPtr;
	collectNodesInMap(ndPtr->left, nmap);
	collectNodesInMap(ndPtr->right, nmap);
}


/**
	Given a query object q, an object p, and a set of objects C that are within interval 
	DI=[nearDistance, farDistance] from p, then the pruning distance is the distance that q is outside the DI.
	NearDisatnce and farDistance are also the distances of the nearest and furthest objects in C from object p.
	pqDistance is the distnace btween objects p and q.
	If if the query radius is less than this distance, then the search tree can prune that node.
**/

template < typename TD>
inline
TD pruningDistance(TD pqDistance, TD nearDistance, TD farDistance) {
	TD pruningDist = 0.0; //prDist of zero means object is inside the interval.
	if (pqDistance < nearDistance) {
		pruningDist = nearDistance - pqDistance;
	}
	else if (pqDistance > farDistance) {
		pruningDist = pqDistance - farDistance;
	}
	return pruningDist;
}



/**
	interiorDistance:
	This is a heuristic for the likelyhood of the object very similar to target being 
	inside an interval (of a node) 
	The function was used as part of an alternative priority in kNN searches.
	**/
template < typename TD>
inline
TD interiorDistance(TD pqDistance, TD nearDistance, TD farDistance) {
	TD pruningDist = 0.0; //prDist of zero means object is inside the annulus.
	if ((pqDistance > farDistance) || (pqDistance < nearDistance)) {
		return 0.0;
	}
	else {
		if (farDistance == nearDistance) {
			return 1;
		}
		else {
			auto num = 1.0 * (std::min((pqDistance - nearDistance), (farDistance - pqDistance)));
			auto den = std::pow(farDistance - nearDistance, 2);
			auto result = num / den;
			return result;
		}
	}
}
/**
distanceLowerBound:
This was for calculating a lower bound on how far a taarget object is outside an interval
(of a node).
The function was used as part of an alternative priority in kNN searches, and PrunningDistance was a better
alternative.
*/
template < typename TD>
inline
TD distanceLowerBound(TD pqDistance, TD nearDistance, TD farDistance) {
	return 	std::max(pqDistance - farDistance, 0.0);
}



#endif // !SEARCH_COMMON_H
