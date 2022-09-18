#ifndef SPMTREE_H
#define SPMTREE_H

/**
 *  @file    SPMTree.h
 *  @author  Miguel R. Zuniga, Jeff K Uhlmann
 *  @date    April 2018 - December 2021
 *
 *  @brief A C++ implementation of the priority metric tree.
 *
 *  @section DESCRIPTION
 *
 *  This is a C++ implementation of the priority metric tree. As a metric tree it is a binary search
 *  tree for ball (k-d range) queries in mertric spaces.  This file also provides a search function
 *	that prioritizes the order of node visitation for kNN searches. Other metric searches (e.g. range, 
 *  collection, range-bound kNN) are supported.
 *
**/

#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <queue>

#include "DistanceInterval.h"
#include "TreeNodes.h"
#include "Query.h"
#include "PerfStats.h"
#include "SearchCommon.h"
#include "DistanceInterval.h"
#include "SearchPQ.h"

template < typename N, typename T, typename M>
class SPMTree_Base {
protected:
	const std::string myShortName{ "SPMT" };
	using Node = N;
	using NodePtr = Node*;
	using NodeItr =typename std::vector <SPMTree_Base <N, T, M>::Node*>::iterator;

	using PQNode = SPQNode<Node>;
	using PQNodePtr = PQNode*;
	//Comparator for the priority queue
	class CompareResultPQ {
	public:
		bool operator()(const PQNode* lhs, const PQNode* rhs) const {
			return lhs->pruningDist > rhs->pruningDist;
		}
	};
	using PQType = SearchPQ<Node, PQNode, CompareResultPQ >;

	void searchK(PQType& pq, const T& target, NearestKQuery<T>& sq, M& met);
	void searchR(NodePtr& nd, const T & target, RadiusQuery<T> & sq, M & met);
	void searchCollect(NodePtr& nd, const T& target, SimilarityQuery<T>& sq, M& met);
	void collect(NodePtr& nd, SimilarityQuery<T>& sq, M& met, const double dist);
	Node* buildTreePCPE(const NodeItr begin, const NodeItr end);
	Node* buildTreeAlt(const NodeItr begin, const NodeItr end);
	auto moveOverlap(const NodeItr start, const NodeItr median, const NodeItr end);
	void partition(const NodeItr firstC, NodeItr& median, const NodeItr end);
	void radiusSumAndDepths(NodePtr& nd, unsigned int depth, double& dsum, std::set<unsigned int>& depths);
	void calculateDI(NodeItr ndPtr, const NodeItr begin,  const NodeItr end);


    Node* root;
	PerfStats perfStats;
	M  metric;

	bool doRangeQueryWithPQ{ true };

	PivotType pivotType{ PivotType::RAN };
	PartType  partitionType{ PartType::BOM };

public:
	//The constructor below also builds the tree.
	SPMTree_Base(std::vector<T>& objects, const M& met, const PivotType pivT, const PartType partT);
	~ SPMTree_Base(){
		std::cout << "SPMTree_Base destructor" << std::endl;
		for(const auto& nd : nodes){
			delete nd;
		}
		nodes.clear();
	}

	void searchCollect(SimilarityQuery<T> & q);	
	void search(RadiusQuery<T> & q);
	void search(NearestKQuery<T>& q);

	double radiusSum();
	PerfStats  & getPerfStats() { return perfStats; }
	std::vector < Node*> nodes;//TODO:move

	double radiusSumAndDepths(std::set<unsigned int>& depths);

	void setPivotType(PivotType pt) { pivotType == pt; }
	PivotType getPivotType() { return pivotType; }
	void setPartType(PartType pt) { partitionType = pt; }
	PartType getPartType() { return partitionType; }

	std::string shortName() const { return myShortName; }

	int getMADIorK() { return 0; } //

};

/*
  Tree building constructor.
  Calls one of the versions of the buildTree() algorithm.
  Calls buildTreeAlt() untill we can settle on the more efficient 
  of pivot selection+ partitioning.
*/
template < typename N, typename T, typename M>
SPMTree_Base <N, T, M>::SPMTree_Base(std::vector < T > & objects, const M& met, const PivotType pivT, const PartType partT) :
	metric{ met }, pivotType{ pivT }, partitionType{ partT } {

	unsigned int nofDIExpected = objects.size();
	
	nodes.reserve(objects.size());
	for (auto& object : objects) {
		nodes.push_back(new Node(&object));
	}

#ifdef SA_USE_STATIC_DI
	DI::initialize(nofDIExpected);
#endif

}

/* Repartition to remove overlaps is any.  */
template < typename N, typename T, typename M>
auto
SPMTree_Base < N, T, M >::moveOverlap(const NodeItr start, const NodeItr median, const NodeItr end) {
	auto newMedian = median;
	if (end - start > 3) {
		auto countL = std::count_if(start, median, EqualByVal<N>((*median)->sstemp));
		//If there are any on the LHS with value equal to the median, repartition.
		if(countL > 0){
			LessThanVal<N> func((*median)->sstemp);
			newMedian = std::partition(start, end, func);
		}
	}
	return newMedian;
}

/*
Build the tree by selecting a center object as the pivot and 
partitioning based on distance to extrema
*/
template < typename N, typename T, typename M>
typename SPMTree_Base < N, T, M >::Node*
SPMTree_Base <N, T, M >::buildTreePCPE(const NodeItr begin, const NodeItr end) {
	if (begin == end) { return nullptr; }
	if ((end - begin) == 1) {
		(*begin)->setLeaf();
		return *begin;
	}

	const auto firstC = begin + 1; //iterator pointing to first child
	const auto median = firstC + (end - firstC) / 2; //median of the children

	const auto nofSamples = 10;//=std::max(10, (int)std::ceil(std::log2(end - begin)));
	NodeItr pivot = findCenter<T, NodeItr>(begin, end, metric, nofSamples);
	std::iter_swap(begin, pivot); //The begin object is now the pivot object.

	NodeItr extrema = findExtrema<T, NodeItr>(firstC, end, metric);

	//First store (temporarily) the distances to the extrema somewhere:
	for (auto nPtr = firstC; nPtr != end; nPtr++) {
		(*nPtr)->sstemp = metric.distance((*nPtr)->object, (*extrema)->object);
	}
	//Then partiton by distance to the extrema
	LessThanTemp<N> ltt;
	std::nth_element(firstC, median, end, ltt);

	//Reset the innerRadiius data to be distance to the pivot object
	for (auto nPtr = firstC; nPtr != end; nPtr++) {
		(*nPtr)->sstemp = metric.distance((*begin)->object, (*nPtr)->object);
	}

	//Set the nodes data (innerRad, outerRad, ,and maxLeft)
	calculateDI(begin, firstC, end);

	(*begin)->left = buildTreePCPE(firstC, median);
	(*begin)->right = buildTreePCPE(median, end);

	return *begin;
}

/**
buildTreeAlt() is the an alternaltive tree building algorithm allows for the selection of different
pivot selection and partitioning algorithms  - mostly for experimental purposes. The other buildTree
algorithms are mostly "fixed" in selection and partitioning algorithms.
**/
template <typename N, typename T, typename M>
typename SPMTree_Base <N, T, M>::Node*
SPMTree_Base <N, T, M >::buildTreeAlt(const NodeItr begin, const NodeItr end) {
	if (begin == end) { return nullptr; }
	if ((end - begin) == 1) {
		(*begin)->setLeaf();
		return *begin;
	}

	const auto firstC = begin + 1; //iterator pointing to first child
	 auto median = firstC + (end - firstC) / 2; //median of the children

	LessThanLen<N> ltl;
	auto pivotItr = selectPivot<T, M, NodeItr, LessThanLen<N>>(begin, median, end, pivotType,ltl, metric);
	std::iter_swap(begin, pivotItr);

	for (auto itr = firstC; itr != end; itr++) {
		(*itr)->sstemp =metric.distance((*begin)->object, (*itr)->object);
	}

	partition(firstC, median, end);
	
	calculateDI(begin, firstC, end);

	(*begin)->left = buildTreeAlt(firstC, median);
	(*begin)->right = buildTreeAlt(median, end);

	return *begin;
}


/*
	Partiton the set in interval [fistC,end) based upon the user selected PartitioningType.
	Note: far is set hold a temp value for patitioning but it is not set back.
*/
template < typename N, typename T, typename M>
void SPMTree_Base < N, T, M >::partition(const NodeItr firstC, NodeItr& median, const NodeItr end) {
if (partitionType == PartType::BOM) {
		//Balanced (|LHS|=|RHS|) Object Median Partition.
		std::nth_element(firstC, median, end, LessThanTemp<N>());
	}
	else {
		error("partition function default");
	}
}

/*
	Find the near and far member variable values of a node. We asume the distances to parents are stored
	in the value far.
*/
template <typename N, typename T, typename M>
void SPMTree_Base<N, T, M>::calculateDI(const NodeItr nd, const NodeItr begin, const NodeItr end) {
  if (begin == end) {
    // for leaf nodes
    (*nd)->di.setNear(0);
    (*nd)->di.setFar(0);

  } else {
    auto near = std::numeric_limits<float>::max();
    auto far = 0.0;
    for (auto p = begin; p != end; p++) {
      if ((*p)->sstemp < near) near = (*p)->sstemp;
      if ((*p)->sstemp > far) far = (*p)->sstemp;
    }
    (*nd)->di.setNear(near);
    (*nd)->di.setFar(far);
  }
}

/*
Users public interface search method which will call the basic search methods.
*/
template <typename N, typename T, typename M>
void SPMTree_Base<N, T, M>::search(RadiusQuery<T>& q) {
  PerfStatsNS::nodesVisited = 0;
  PerfStatsNS::distanceCalls = 0;

  if (root != nullptr) {
    searchR(root, q.getTarget(), q, metric);
  }

  perfStats.incNodesVisited(PerfStatsNS::nodesVisited);
  perfStats.incDistanceCalls(PerfStatsNS::distanceCalls);
}

template <typename N, typename T, typename M>
void SPMTree_Base<N, T, M>::search(NearestKQuery<T>& q) {
  PerfStatsNS::nodesVisited = 0;
  PerfStatsNS::distanceCalls = 0;

  if (root != nullptr) {
    PQType queue;
    auto dist = metric.distance(&(q.getTarget()), root->object);
	auto pd = pruningDistance<double>(dist, root->di.getNear(), root->di.getFar());
    PerfStatsNS::distanceCalls++;
	PerfStatsNS::nodesVisited++;
    queue.push(queue.newNode(root, nullptr, dist, pd));
    searchK(queue, q.getTarget(), q, metric);
  }

  perfStats.incNodesVisited(PerfStatsNS::nodesVisited);
  perfStats.incDistanceCalls(PerfStatsNS::distanceCalls);
}

//SearchCollect should only work with RadiusQuery or RadiusCountQuery
template < typename N,  typename T, typename M>
void SPMTree_Base < N, T, M >::searchCollect(SimilarityQuery<T>& q) {
	if (this->root != nullptr) {
		std::vector <double> distStack;
		searchCollect(root, q.getTarget(), q, metric);
	}
}
//KNN search implementation.
template < typename N, typename T, typename M>
void SPMTree_Base <N,  T, M >::searchK(PQType& pq, const T& target, NearestKQuery<T>& sq, M& met) {
	while (!pq.empty()) {
		PQNodePtr qn = pq.top();
		pq.pop();
		sq.addResult(qn->node->object, qn->distance);	
		if(qn->pruningDist <= sq.searchRadius()) {
			auto lnd = qn->node->left;
			auto rnd = qn->node->right;
			if (lnd != nullptr) {
				PerfStatsNS::nodesVisited++; 
				auto dist = met.distance(&target, lnd->object);
				PerfStatsNS::distanceCalls++;
				auto pd = pruningDistance<double>(dist, lnd->di.getNear(), lnd->di.getFar());
				pq.push(pq.newNode(lnd, qn, dist, pd));
			}
			if (rnd != nullptr) {
				PerfStatsNS::nodesVisited++; 
				auto dist = met.distance(&target, rnd->object);
				PerfStatsNS::distanceCalls++;
				auto pd = pruningDistance<double>(dist, rnd->di.getNear(), rnd->di.getFar());
				pq.push(pq.newNode(rnd, qn, dist, pd));
			}
		}
	}
}
/*
   A simple version of the search function - traditional without the use of the priority queue.
   Mostly for legacy and test purposes.
*/

template <typename N, typename T, typename M>
void SPMTree_Base<N, T, M>::searchR(NodePtr& nd, const T& target, RadiusQuery<T>& sq, M& met) {
  PerfStatsNS::nodesVisited++;

  if (sq.searchRadius() < 0) {
    return;
  }

  const auto dist = met.distance(&target, nd->object);
  
  PerfStatsNS::distanceCalls++;
  sq.addResult(nd->object, dist);

  auto pd = pruningDistance<double>(dist, nd->di.getNear(), nd->di.getFar());
  if (pd <= sq.searchRadius()) {
    if (nd->left != nullptr) {
      searchR(nd->left, target, sq, met);
    }
    if (nd->right != nullptr) {
      searchR(nd->right, target, sq, met);
    }
  }
}

template < typename N, typename T, typename M>
void SPMTree_Base < N, T, M >::searchCollect(NodePtr& nd, const T& target, SimilarityQuery<T>& sq, M& met) {
	if (nd == nullptr) return;
		perfStats.incNodesVisited();
		auto dist = met.distance(&target, nd->object);
		perfStats.incDistanceCalls();

		if(dist + nd->di.getFar() <=  sq.searchRadius()){
			collect(nd, sq, met, dist + nd->di.getFar() );
		}else{
			sq.addResult(nd->object, dist);
			auto pd = pruningDistance<double>(dist, nd->di.getNear(), nd->di.getFar());
			if (pd <= sq.searchRadius()) {
				searchCollect(nd->left, target, sq, met);
				searchCollect(nd->right, target, sq, met);
			}	
		}
	return;
	}

/**
	Collect all objects of the passed subtree, which are assumed to be 
	within distance dist from the query object. The same dist argument is reported for each 
	iindividual result(i.e. ) its not refined here at all) and therefore for an individual result it 
	may merely be an upper bound of the pq distance and "dist <= sq.searchRadius()" should hold.
	TODO: Make a version that reports (in addResult) when dist is from collect and when it is not,
	so user can know id the pq dist is actual or an upper bound.
**/
template < typename N, typename T, typename M>
void SPMTree_Base <N, T, M >::collect(NodePtr& nd, SimilarityQuery<T>& sq, M& met, const double dist) {
	if(nd == nullptr){
		return;
	}else{
		sq.addResult(nd->object, dist);
		collect(nd->left,sq, met, dist);
		collect(nd->right, sq, met, dist);
	}
}

/*
	
*/
template < typename N, typename T, typename M>
double SPMTree_Base < N, T, M>::radiusSumAndDepths(std::set<unsigned int>& depths) {
	double sum = 0.0;
	radiusSumAndDepths(root, 0, sum, depths);
	return sum;
}

/*
	Traverse the tree to collect the sum of the radius and the depths of the tree.
	(These are used as measures of tree building correctness)
*/

template < typename N, typename T, typename M>
void SPMTree_Base < N, T, M >::radiusSumAndDepths(NodePtr& nd,  unsigned int depth, double& dsum, std::set<unsigned int>& depths) {
	if (nd->isLeaf()) {
		depths.insert(depth);
		return;
	}
	//if (di.size() > 0)
	dsum += nd->di.getFar();
	if (nd->left != nullptr)
		radiusSumAndDepths(nd->left, depth + 1, dsum, depths);
	if (nd->right != nullptr)
		radiusSumAndDepths(nd->right, depth + 1, dsum, depths);
}


/*
	SPMTree Does not support priority score queries, so just  accept all of the functionality of SPMTree_BAase.
*/
#include "TreeNodes.h"
template < typename T, typename M>
class SPMTree : public  SPMTree_Base<MNode<T>, T, M> {
	const std::string myShortName{ "SPMT" };

	//using SPMTree_Base<MNode<T>, T, M>::SPMTree_Base;
public:
	SPMTree(std::vector < T >& objects, const M& met,
		PivotType pivT = PivotType::RAN, PartType partT = PartType::BOM, const unsigned int madi = 0) :
	      SPMTree_Base<MNode<T>, T, M>::SPMTree_Base(objects, met, pivT, partT) {
		std::cout << "SPMTree buildTree starting" << std::endl;
		this->root = this->buildTreeAlt(this->nodes.begin(), this->nodes.end());
		std::cout << "SPMTree buildTree finished" << std::endl;

	}

	std::string shortName() const { return myShortName; }

};

#endif
