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

template < typename T, typename M>
class APMTree : public  SPMTree_Base<ANode<T>, T, M> {
	protected:
	const std::string myShortName{ "APMT" };
	using Node = ANode<T>;  
	using NodePtr = Node*;
	using NodeItr = typename std::vector<APMTree<T, M>::Node*>::iterator;
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
	Node* buildTreeAlt(const NodeItr begin, const NodeItr end);
	void partition(const NodeItr firstC, NodeItr& median, const NodeItr end);
	void radiusSumAndDepths(NodePtr& nd, unsigned int depth, double& dsum, std::set<unsigned int>& depths);
	void calculateDI(NodeItr ndPtr, const NodeItr begin,  const NodeItr median, const NodeItr end);
public:
	APMTree(std::vector < T >& objects, const M& met,
		PivotType pivT = PivotType::RAN, PartType partT = PartType::BOM, const unsigned int madi = 0) :
	      SPMTree_Base<ANode<T>, T, M>::SPMTree_Base(objects, met, pivT, partT) {
	std::cout << "APMTree buildTree starting" << std::endl;
	this->root = this->buildTreeAlt(this->nodes.begin(), this->nodes.end());
	std::cout << "APMTree buildTree finished" << std::endl;

	}

	void searchCollect(SimilarityQuery<T> & q);	
	void search(RadiusQuery<T> & q);
	void search(NearestKQuery<T>& q);
	double radiusSumAndDepths(std::set<unsigned int>& depths);

	std::string shortName() const { return myShortName; }

};

/**
buildTreeAlt() is the an alternaltive tree building algorithm allows for the selection of different
pivot selection and partitioning algorithms  - mostly for experimental purposes. The other buildTree
algorithms are mostly "fixed" in selection and partitioning algorithms.
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

	for (auto itr = firstC; itr != end; itr++){
		(*itr)->sstemp = this->metric.distance((*begin)->object, (*itr)->object);
	}

	this->partition(firstC, median, end);

	this->calculateDI(begin, firstC, median, end);

	(*begin)->left = buildTreeAlt(firstC, median);
	(*begin)->right = buildTreeAlt(median, end);

	return *begin;
}

/*
	Partiton the set in interval [fistC,end) based upon the user selected PartitioningType.
	Note: far is set hold a temp value for patitioning but it is not set back.
*/
template <typename T, typename M>
void APMTree<T, M >::partition(const NodeItr firstC, NodeItr& median, const NodeItr end) {
if (this->partitionType == PartType::BOM) {
		//Balanced (|LHS|=|RHS|) Object Median Partition.
		std::nth_element(firstC, median, end, LessThanTemp<Node>());
	}
	else {
		error("partition function default");
	}
}

/*

/*
	Find and store the bounding intervals of te left and right child. 
	We asume the distances to parents are stored
	in the value far.
*/
template <typename T, typename M>
void APMTree<T, M>::calculateDI(const NodeItr nd, const NodeItr begin, const NodeItr median, const NodeItr end)
{
	(*nd)->diL.setNear(0);
	(*nd)->diL.setFar(0);
	(*nd)->diR.setNear(0);
	(*nd)->diR.setFar(0);

	if (begin == end)
		return;

	auto near = std::numeric_limits<float>::max();
	auto far = 0.0;
	for (auto p = begin; p != median; p++){
		if ((*p)->sstemp < near)
			near = (*p)->sstemp;
		if ((*p)->sstemp > far)
			far = (*p)->sstemp;
	}
	(*nd)->diL.setNear(near);
	(*nd)->diL.setFar(far);
	// And same for RHS child node:
	near = std::numeric_limits<float>::max();
	far = 0.0;
	for (auto p = median; p != end; p++){
		if ((*p)->sstemp < near)
			near = (*p)->sstemp;
		if ((*p)->sstemp > far)
			far = (*p)->sstemp;
	}
	(*nd)->diR.setNear(near);
	(*nd)->diR.setFar(far);
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
	auto pd = pruningDistance<double>(dist, this->root->diL.getNear(), this->root->diR.getFar());
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
		auto tnd = qn->node; // tree node
		sq.addResult(tnd->object, qn->distance);	
		if(qn->pruningDist <= sq.searchRadius()) {
			auto pqDist = qn->distance;
			auto lnd = tnd->left;
			auto rnd = tnd->right;
			if (lnd != nullptr) {
				if ( pqDist <= tnd->diL.getFar() + sq.searchRadius())  {
					this->perfStats.incNodesVisited();
					auto dist = met.distance(&target, lnd->object);
					this->perfStats.incDistanceCalls();
					auto pd = pruningDistance<double>(dist, 
								 std::min(lnd->diL.getNear(),  lnd->diR.getNear()), 
								 std::max(lnd->diL.getFar(),   lnd->diR.getFar()));
					pq.push(pq.newNode(lnd, qn, dist, pd));
			}
			}
			if (rnd != nullptr) {
				if (pqDist >= tnd->diR.getNear() - sq.searchRadius()){
					this->perfStats.incNodesVisited();
					auto dist = met.distance(&target, rnd->object);
					this->perfStats.incDistanceCalls();
					auto pd = pruningDistance<double>(dist,
													  std::min(rnd->diL.getNear(), rnd->diR.getNear()),
													  std::max(rnd->diL.getFar(), rnd->diR.getFar()));
					pq.push(pq.newNode(rnd, qn, dist, pd));
				}
		}
	}
	}
}

	// auto pd = pruningDistance<double>(dist, std::min(nd->diL.getNear(), nd->diR.getNear()), nd->diR.getFar());
	/// Worked:
	/*
	if (nd->left != nullptr){
		auto pd = pruningDistance<double>(dist, nd->diL.getNear(), nd->diL.getFar());
		if (pd <= sq.searchRadius()){
			searchR(nd->left, target, sq, met);
		}
	}
	if (nd->right != nullptr){
		auto pd = pruningDistance<double>(dist, nd->diR.getNear(), nd->diR.getFar());
		if (pd <= sq.searchRadius()) {
			searchR(nd->right, target, sq, met);
		}
	}
	****/
template <typename T, typename M>
void APMTree<T, M>::searchR(NodePtr& nd, const T& target, RadiusQuery<T>& sq, M& met) {
 	if (nd == nullptr)
		return;
  	
	this->perfStats.incNodesVisited();

 	auto dist = met.distance(&target, nd->object);
 	this->perfStats.incDistanceCalls();

  	sq.addResult(nd->object, dist);

	if (nd->left != nullptr){
		if (dist <= nd->diL.getFar() + sq.searchRadius()){
			searchR(nd->left, target, sq, met);
		}
  	}
  	if (nd->right != nullptr){
	  	if (dist >= nd->diR.getNear() - sq.searchRadius()){
			  searchR(nd->right, target, sq, met);
	  	}
  	}
}

template <typename T, typename M>
void APMTree<T, M >::searchCollect(NodePtr& nd, const T& target, SimilarityQuery<T>& sq, M& met) {
	std::cout <<"This function is not finished and  tested" << std::endl; 
	std::exit(-1);
	if (nd == nullptr) return;
		
	this->perfStats.incNodesVisited();
	auto dist = met.distance(&target, nd->object);
	this->perfStats.incDistanceCalls();

	if(dist + nd->diR.getFar() <=  sq.searchRadius()){
		collect(nd, sq, met, dist + nd->diR.getFar() );
	}else{
		sq.addResult(nd->object, dist);
		auto pd = pruningDistance<double>(dist, nd->diR.getNear(), nd->diR.getFar());
		if (pd <= sq.searchRadius()) {
			if (nd->left != nullptr){
				if (dist <= nd->diL.getFar() + sq.searchRadius()){
						searchCollect(nd->left, target, sq, met);
				}
			}
			if (nd->right != nullptr){
	  			if (dist >= nd->diR.getNear() - sq.searchRadius()){
						searchCollect(nd->right, target, sq, met);
				}	
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
	dsum += nd->diR.getFar();
	if (nd->left != nullptr)
		radiusSumAndDepths(nd->left, depth + 1, dsum, depths);
	if (nd->right != nullptr)
		radiusSumAndDepths(nd->right, depth + 1, dsum, depths);
}


#endif
