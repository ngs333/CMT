#ifndef APMTREE_H
#define APMTREE_H

/**
 *  @file    APMTree.h
 *  @date    June 2022
 *
 *  @brief A C++ implementation of the priority metric tree.
 *
 *  @section DESCRIPTION
 *
 *  This is a C++ implementation of the priority ball metric tree.  Tree nodes store the bounding interval
 *  of pivot distances of the left child and right child. This is roughly equivalent to the original 
 *  Ball Metric Tree of Uhlmann (1991) and the VP-1 tree of Yanilos (1993) since the internal distances
 *  of the two intervals are roughly equivalent to the partitioning distance and its the most usefull
 *  info for pruning (The external distances of the two intervals are significantly less usefull). 
 *  Added is a priority queue used for determining node order visitaation in kNN search. This version 
 *  also supports collection search and range bounded kNN search.
 *
**/

#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <queue>

#include "TreeNodes.h"
#include "Query.h"
#include "PerfStats.h"
#include "SearchCommon.h"
#include "DistanceInterval.h"
#include "SearchPQ.h"
#include "SPMTree.h"
#include "PartitionFunction.h"

template < typename T, typename M>
class APMTree : public  SPMTree_Base<ANode<T>, T, M> {
protected:
	const std::string myShortName{ "APMT" };
	using Node = ANode<T>;  
	using NodePtr = Node*;
	using NodeItr = typename std::vector<Node*>::iterator;
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
	//void collect(NodePtr& nd, SimilarityQuery<T>& sq, M& met, const double dist);
	Node* buildTreePiv(const NodeItr begin, const NodeItr end);
	Node* buildTreeAlt(const NodeItr begin, const NodeItr end);
	void partition(const NodeItr firstC, NodeItr& median, const NodeItr end);
	void radiusSumAndDepths(NodePtr& nd, unsigned int depth, double& dsum, std::set<unsigned int>& depths);
	void calculateDI(NodeItr ndPtr, const NodeItr begin,  const NodeItr median, const NodeItr end);
	PartitionFunction<Node,NodeItr,T,M>* pFunction;

public:
	APMTree(std::vector < T >& objects, const M& met,
		PivotType pivT = PivotType::RAN, PartType partT = PartType::BOM, const unsigned int madi = 0) :
	      SPMTree_Base<ANode<T>, T, M>::SPMTree_Base(objects, met, pivT, partT) {
		std::cout << "APMTree buildTree starting" << std::endl;
		if( (pivT == PivotType::CENT) && (partT==PartType::EXT)){
			this->pFunction = getPartitionFunctor<Node,NodeItr,T,M>( PartType::DMR );
			this->root = this->buildTreeAlt(this->nodes.begin(), this->nodes.end());
		} else{
			this->pFunction = getPartitionFunctor<Node,NodeItr,T,M>( this->partitionType );
			this->root = this->buildTreePiv(this->nodes.begin(), this->nodes.end());
		}
		std::cout << "APMTree buildTree finished" << std::endl;
	}
~APMTree(){
		std::cout << "APMTree destructor" << std::endl;
		delete pFunction;
}
	void searchCollect(SimilarityQuery<T> & q);	
	void search(RadiusQuery<T> & q);
	void search(NearestKQuery<T>& q);
	double radiusSumAndDepths(std::set<unsigned int>& depths);

	std::string shortName() const { return myShortName; }



};

/**
buildTreePiv() tree building algorithm where the distances for partitioning are calculated relative
to the pivot. The combinations of PivotType/PartType that is suppoerted here includes these six:
[RAN,CEN,EXT]/[BOM,DMR].
**/

template <typename T, typename M>
typename APMTree<T, M>::Node *
APMTree<T, M>::buildTreePiv(const NodeItr begin, const NodeItr end)
{
	if (begin == end){
		return nullptr;
	}
	int size =end - begin;

	if ( size == 1){
		(*begin)->setLeaf();
		return *begin;
	}

	const auto firstC = begin + 1;			   // iterator pointing to first child
	auto median = firstC + (end - firstC) / 2; // median of the children

	LessThanLen<Node> ltl;
	auto pivotItr = selectPivot<T, M, NodeItr, LessThanLen<Node>>(begin, median, end,
												this->pivotType, ltl, this->metric);
	std::iter_swap(begin, pivotItr);
	(*begin)->size= size;

	for (auto itr = firstC; itr != end; itr++){
		(*itr)->sstemp = this->metric.distance((*begin)->object, (*itr)->object);
	}

	(*pFunction)(firstC, median, end);

	this->calculateDI(begin, firstC, median, end);

	(*begin)->left = this->buildTreePiv(firstC, median);
	(*begin)->right = this->buildTreePiv(median, end);

	return *begin;
}


/**
buildTreeAlt() is the tree building algorithm where the distances for partitioning are 
relative to an object (an Alternative) that is not pivot. The PivotType/PartType 
combination(s) suppported here is: [CENT/EXT]. Not that in the partitioing step,
DMR is used once distances are calculated relative to the extrema.
**/
template <typename T, typename M>
typename APMTree<T, M>::Node *
APMTree<T, M>::buildTreeAlt(const NodeItr begin, const NodeItr end)
{
	if (begin == end){
		return nullptr;
	}
	int size =end - begin;

	if ( size == 1){
		(*begin)->setLeaf();
		return *begin;
	}

	const auto firstC = begin + 1;			   // iterator pointing to first child
	auto median = firstC + (end - firstC) / 2; // median of the children

	LessThanLen<Node> ltl;
	auto pivotItr = selectPivot<T, M, NodeItr, LessThanLen<Node>>(begin, median, end,
												this->pivotType, ltl, this->metric);
	std::iter_swap(begin, pivotItr);
	(*begin)->size= size;

	//Find an extrema among [fristC, end)]
	auto extItr = findExtrema<T, NodeItr>(firstC, end, this->metric);
	//Set distances to the extrema
	for (auto itr = firstC; itr != end; itr++) {
		(*itr)->sstemp  = this->metric.distance((*extItr)->object, (*itr)->object);
	}
	//partition by those distances
	//For EXT, the partition function was set tp BOM or DMR
	(*pFunction)(firstC, median, end);

	//Now calc and set distances to the pivot.
	//Note that the distance to tree root is stored at dpivots[0] for any node.
	for (auto itr = firstC; itr != end; itr++){
		(*itr)->sstemp = this->metric.distance((*begin)->object, (*itr)->object);
	}

	this->calculateDI(begin, firstC, median, end);

	(*begin)->left = buildTreeAlt(firstC, median);
	(*begin)->right = buildTreeAlt(median, end);

	return *begin;
}


/*
	Find and store the bounding intervals of te left and right child. 
	We asume the distances to parents are stored
	in the value far.
*/
template <typename T, typename M>
void APMTree<T, M>::calculateDI(const NodeItr nd, const NodeItr begin, 
	const NodeItr median, const NodeItr end){
	(*nd)->diL.reset();
	(*nd)->diR.reset();

	if (begin == end)
		return;
	LessThanTemp<Node> cf;
	calculateDBI<Node,NodeItr> ( (*nd)->diL, begin, median, cf);
	calculateDBI<Node,NodeItr> ( (*nd)->diR, median, end, cf);
}

/*
Users public interface search method which will call the basic search methods.
*/
template <typename T, typename M>
void APMTree<T, M>::search(RadiusQuery<T>& q) {
  if (this->root != nullptr) {
    searchR(this->root, q.getTarget(), q, this->metric);
  }
}

template <typename T, typename M>
void APMTree<T, M>::search(NearestKQuery<T>& q) {
  if (this->root != nullptr) {
    PQType queue;
    auto dist = this->metric.distance(&(q.getTarget()), this->root->object);
	 this->perfStats.incNodesVisited();
  	this->perfStats.incDistanceCalls();
	auto pd = pruningDistance<double>(dist, 0, 1.0e30);
    queue.push(queue.newNode(this->root, nullptr, dist, pd));
	searchK(queue, q.getTarget(), q, this->metric);
  }
}

//SearchCollect should only work with RadiusQuery or RadiusCountQuery
template <typename T, typename M>
void APMTree<T, M >::searchCollect(SimilarityQuery<T>& q) {
	if (this->root != nullptr) {
		std::vector <double> distStack;
		searchCollect(this->root, q.getTarget(), q, this->metric);
	}
}

//KNN search implementation.
template <typename T, typename M>
void APMTree<T, M >::searchK(PQType& pq, const T& target, NearestKQuery<T>& sq, M& met) {
	while (!pq.empty()) {
		PQNodePtr qn = pq.top();
		pq.pop();
		auto pnd = qn->node; // tree node
		sq.addResult(pnd->object, qn->distance);	
		if(qn->pruningDist <= sq.searchRadius()) {
			auto pqDist = qn->distance;
			auto lnd = pnd->left;
			auto rnd = pnd->right;
			if (lnd != nullptr) {
				if (pnd->diL.rangeOverlaps(pqDist, sq.searchRadius()) == true) {
					this->perfStats.incNodesVisited();
					auto dist = met.distance(&target, lnd->object);
					this->perfStats.incDistanceCalls();
					auto pd = pruningDistance<double>(dist, lnd->getNear(), lnd->getFar());
					pq.push(pq.newNode(lnd, qn, dist, pd));
				}
			}
			if (rnd != nullptr) {
				if (pnd->diR.rangeOverlaps(pqDist, sq.searchRadius()) == true){
					this->perfStats.incNodesVisited();
					auto dist = met.distance(&target, rnd->object);
					this->perfStats.incDistanceCalls();
					auto pd = pruningDistance<double>(dist, rnd->getNear(), rnd->getFar());
					pq.push(pq.newNode(rnd, qn, dist, pd));
				}
			}
		}
	}
}

/*
Range search function
**/
template <typename T, typename M>
void APMTree<T, M>::searchR(NodePtr& nd, const T& target, RadiusQuery<T>& sq, M& met) {
    if (nd == nullptr)
        return;
    this->perfStats.incNodesVisited();
    auto dist = met.distance(&target, nd->object);
    this->perfStats.incDistanceCalls();
    sq.addResult(nd->object, dist);

    if (nd->left != nullptr) {
		if (nd->diL.rangeOverlaps(dist, sq.searchRadius())) {	
            searchR(nd->left, target, sq, met);
        }
    }

    if (nd->right != nullptr) {
        if (nd->diR.rangeOverlaps(dist, sq.searchRadius())) {
            searchR(nd->right, target, sq, met);
        }
    }
}

template <typename T, typename M>
void APMTree<T, M>::searchCollect(NodePtr& nd, const T& target, SimilarityQuery<T>& sq, M& met) {
    if (nd == nullptr)
        return;

    this->perfStats.incNodesVisited();
    auto dist = met.distance(&target, nd->object);
    this->perfStats.incDistanceCalls();

    auto fd = nd->getFar();
    if (dist + fd <= sq.searchRadius()) {
        this->collect(nd, sq, met, dist + fd);
    } else {
        sq.addResult(nd->object, dist);
        if (nd->left != nullptr) {
            if (nd->diL.rangeOverlaps(dist, sq.searchRadius())) {
                searchCollect(nd->left, target, sq, met);
            }
        }
        if (nd->right != nullptr) {
            if (nd->diR.rangeOverlaps(dist, sq.searchRadius())) {
                searchCollect(nd->right, target, sq, met);
            }
        }
    }
    return;
}

/*
	
*/
template <typename T, typename M>
double APMTree<T, M>::radiusSumAndDepths(std::set<unsigned int>& depths) {
	double sum = 0.0;
	radiusSumAndDepths(this->root, 0, sum, depths);
	return sum;
}
/*
	Traverse the tree to collect the sum of the radius and the depths of the tree.
	(These are used as measures of tree building correctness)
*/

template <typename T, typename M>
void APMTree<T, M >::radiusSumAndDepths(NodePtr& nd,  unsigned int depth, double& dsum, std::set<unsigned int>& depths) {
	if (nd->isLeaf()) {
		depths.insert(depth);
		return;
	}
	//if (di.size() > 0)
	dsum += nd->getFar();
	if (nd->left != nullptr)
		radiusSumAndDepths(nd->left, depth + 1, dsum, depths);
	if (nd->right != nullptr)
		radiusSumAndDepths(nd->right, depth + 1, dsum, depths);
}


#endif
