#ifndef CMTREE_H
#define CMTREE_H

/**
 *  @file    CMTree.h
 *  @author  Miguel R. Zuniga, Dr. Jeffrey K. Uhlmann
 *  @date    Sep 2018 - December 2021
 *
 *  @brief This is the reference implementation of the Cascaded Metric Tree (CMT)
 *
 *  @section DESCRIPTION
 *
 *  This is the reference implementation of The Cascaded Metric Tree (CMT) for performing metric search queries.
 *  The data structures are presented in papers on the tree by Jeffrey K. Uhlmann and Miguel R. Zuniga, e.g. 
 *  "The Cascading Metric Tree", submitted 2021.
 * 
 *  The original idea of augmenting the Metric Tree with the added storage of the NlogN tree construction 
 *  metric evaluations  is from Dr. Jeff Uhlmann. Some C++ implementation ideas is borrowed from the work of
 *  Seth Weisman and Yeshwanthi Pachalla.
 *  
 * 
**/

#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ostream>

//#include <boost/circular_buffer.hpp>

#include "DistanceInterval.h"
#include "TreeNodes.h"
#include "SearchCommon.h"
#include "Query.h"
#include "PerfStats.h"
#include "SearchPQ.h"

template < typename N, typename T, typename M>
class CMTree_Base {
protected:
	using Node = N;
	using NodePtr = Node*;
	using NodeItr = typename std::vector<CMTree_Base<N, T, M>::Node*>::iterator;

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

	
	NodePtr root;
	std::vector < NodePtr> nodes;
	//std::vector < Node> nodeMemory;
	NodePtr buildTreeBOM(const NodeItr begin, const NodeItr end);
	NodePtr buildTreeAlt(const NodeItr begin, const NodeItr end);
	NodePtr buildTreeDMR(const NodeItr begin, const NodeItr end);
	NodePtr buildTreePCPE(const NodeItr begin, const NodeItr end);
	void partition(const NodeItr begin,  NodeItr& median, const NodeItr end);
	auto moveOverlap(const NodeItr begin, const NodeItr median, const NodeItr end);
	void searchK(PQType& nodeQueue, T& target, SimilarityQuery<T>& sq, M& met);
	void searchR(NodePtr& nd, T& target, SimilarityQuery<T>& q, M& met, std::vector<double>& dpivots);
	void searchCollect(NodePtr& nd, T& target, SimilarityQuery<T>& q, M& met, std::vector<double>& dpivots);
	void collect(NodePtr& nd, SimilarityQuery<T>& sq, M& met, const double dist);
	void computeADIVec(const NodeItr begin, const NodeItr end);
	void calculateDI(const NodeItr nd, const NodeItr fc, const NodeItr end);
	void radiusSumAndDepths(NodePtr& nd, unsigned int depth, double& dsum, std::set<unsigned int>& depths);
	double maxPruningDistance(PQNodePtr& pqNodePtr, NodePtr& nd);
	double maxPruningDistance(NodePtr& nd, std::vector <double>& distStack);
	double minCollectionDistance(NodePtr& nd, std::vector <double>& distStack);
	void collectSpacePivots(NodePtr& nd,const unsigned int nSpace, unsigned int& count,  std::vector<T> &newPivots, std::vector<T> &newNodes);
	static bool IntervalInsideQueryVol(double pqDist, double qRadius, DI& adi);
	void freePivDistances(NodePtr& nd);

	PerfStats perfStats;
	M  metric;

	//Some configuration info with defaults:
	bool doRangeQueryWithPQ{ false };
	PivotType pivotType{ PivotType::RAN };
	PartType  partitionType{ PartType::BOM};
	//const unsigned int PQ_MEM_SIZE{ 40000 };
	std::string myShortName{ "CMT_B" };

	//Maximum number of distance intervals per node in unlimited case:
	inline static unsigned int MAX_NDI{ 1000000 }; 
	//Maximum number of ancestral distance intervals per node (excluding the the intervals of the
	//child nodes of the node. This should be 0 for the metric tree equivalent, and
	//1  if only once ancestor is used for the CMT: Note: MDI = MADI + 1
	unsigned int  MADI; 
	//TODO: INTERVAL_DX has to be thought thorough 
	//const float INTERVAL_DX = 0.0000001;
	const float INTERVAL_DX = 0;

	//NodeItr globalExt ;

public:
	CMTree_Base(std::vector < T >& objects, const M& met, 
		PivotType pivT = PivotType::RAN, PartType partT = PartType::BOM, unsigned int maxdi = MAX_NDI);
	virtual ~CMTree_Base() {
		for (auto& nd : nodes) {
			delete nd;
		}
		//	std::vector<Node*>().swap(nodes);
		nodes.clear();
	}

	void searchCollect(SimilarityQuery<T> & q);	
	void search(RadiusQuery<T> & q);
	void search(NearestKQuery<T>& q);
	void search(RadiusCountQuery<T> & q);
	
	PerfStats& getPerfStats() { return perfStats; }
	double radiusSumAndDepths(std::set<unsigned int>& depths);
    void collectSpacePivots(const unsigned int nPivots, std::vector<T> &newPivots, std::vector<T> &newNodes);
	void printTree(std::ostream& os){ printBinaryTree<NodePtr>(root, os); }
	void printDI(std::ostream& os);
	CMTStructStat<T> getStructureStats();

	void setPivotType(PivotType pt) { pivotType == pt; }
	PivotType getPivotType() { return pivotType; }
	void setPartType(PartType pt) { partitionType = pt; }
	PartType getPartType() { return partitionType; }
	int getMADIorK() { return MADI; }

	std::string shortName() const { return myShortName; }

};


template < typename N, typename T, typename M>
CMTree_Base <N, T, M>::CMTree_Base(std::vector < T > & objects, const M& met, PivotType pivT, PartType partT, unsigned int maxadi) :
	metric{ met }, pivotType{ pivT }, partitionType{ partT }, MADI{ maxadi }{

	//Tree height (max depth) if perfectly balanced
	unsigned int treeHeight = std::ceil(std::log2(objects.size()));
	unsigned int nofDIExpected = (treeHeight + 1) * objects.size();

	std::cout << "nofDIExpected=" << nofDIExpected
		<< ", sizeof(Node*)=" << sizeof(Node*)
		<< ", sizeof(Node)=" << sizeof(Node)
		<< ", sizeof(VectorFM<DI>)=" << sizeof(FSVector<DI>)
		<< ", sizeof(DI)=" << sizeof(DI)
		<< ", sizeof(std::vector<float>)=" << sizeof(std::vector<float>)
		<< std::endl;

	nodes.reserve(objects.size());
	for (auto& object : objects) {
		nodes.push_back(new Node(&object));
	}

	for (auto& node : nodes) {
		//node->dpivots.reserve(treeHeight + 1);
		node->sstemp = 0;
	}


#ifdef SA_USE_STATIC_DI
	DI::initialize(nofDIExpected);
#endif

	return;

	
}

template < typename N, typename T, typename M>
typename CMTree_Base <N, T, M >::Node*
CMTree_Base < N, T, M >::buildTreeBOM(const NodeItr begin, const NodeItr end) {
	if (begin == end) { return nullptr; }
	if ((end - begin) == 1) {
		(*begin)->setLeaf();
		computeADIVec(begin, end);
		return *begin;
	}

	//Pivot will be a random object
	auto rndIter = findRandomIter< NodeItr>(begin, end);
	std::iter_swap(begin, rndIter);

	const NodeItr firstC = begin + 1; //iterator pointing to first child
	NodeItr median = firstC + (end - firstC) / 2; //median of the children

	/**
	Store the child node distances to the pivot! If a nodeis at level L (where root is at level 0),
	the distance to its immediate parent is stored at dpivots[l-1], and to its grandparent at dpivots[l-2], etc.
	**/
	for (auto itr = firstC; itr != end; itr++) {
		auto dist = metric.distance((*begin)->object, (*itr)->object);
		(*itr)->dpivots.push_back(dist);
		(*itr)->sstemp = dist;
	}
	calculateDI(begin, firstC, end);

	//Partition into two sets by distance to the pivot object just calculated.
	LessThanTemp<N> ltt;
	std::nth_element(firstC, median, end, ltt);

	(*begin)->left = buildTreeBOM(firstC, median);
	(*begin)->right = buildTreeBOM(median, end);

	//Compute and set the adi for the begin node;
	computeADIVec(begin, end);

	return *begin;
}

/* Repartition to remove overlaps is any.  */
template < typename N, typename T, typename M>
auto
CMTree_Base < N, T, M >::moveOverlap(const NodeItr start, const NodeItr median, const NodeItr end) {
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


/**
Compute the ancestral distance interval vector (ADI[]) for the node pointed to by the begin 
iterator. There are L distance intervals for node at level L. The root node is at level 0.
**/
template < typename N, typename T, typename M>
void CMTree_Base < N, T, M >::computeADIVec(const NodeItr begin, const NodeItr end) {
	struct ComparePD {
		unsigned int ld;
		ComparePD(unsigned int l) : ld{ l } {};
		bool operator()(NodePtr& lhs, NodePtr& rhs) const
		{ return lhs->dpivots[ld] < rhs->dpivots[ld]; }
	};

	const unsigned int depth = (*begin)->dpivots.size();
	const unsigned int intervals = std::min(depth, MADI);
	if (intervals == 0) return;//tree root has no ancestral intervals

	(*begin)->adi.reserve(intervals);

	for (auto ii = 0; ii < intervals; ii++) {//Index (and counter) to distance interval
		unsigned int pdi = depth - ii - 1;
		//unsigned int pdi = ii;
		ComparePD cmp(pdi);
		const auto nearest  = std::min_element(begin, end, cmp);
		const auto furthest = std::max_element(begin, end, cmp);
		(*begin)->adi.push_back(DI((*nearest)->dpivots[pdi], ((*furthest)->dpivots[pdi]) + INTERVAL_DX));
	}
	return;
}


template <typename N,  typename T, typename M>
typename CMTree_Base <N, T, M>::Node*
CMTree_Base < N, T, M >::buildTreePCPE(const NodeItr begin, const NodeItr end) {
	if (begin == end) { return nullptr; }
	if ((end - begin) == 1) {
		(*begin)->setLeaf();
		computeADIVec(begin, end);
		return *begin;
	}
	const auto firstC = begin + 1; //iterator pointing to first child
	auto median = firstC + (end - firstC) / 2;

	//Pivot will be the object in the "center":
	const auto nofSamples = 10;//=std::max(10, (int)std::ceil(std::log2(end - begin)));
	NodeItr pivot = findCenter<T, NodeItr>(begin, end, metric, nofSamples);
	std::iter_swap(begin, pivot); //The begin object is now the pivot object.

	NodeItr extrema = findExtrema<T, NodeItr>(firstC, end, metric);

	//if(globalExt == end){globalExt = findExtrema<T, NodeItr>(begin, end, metric);}

	//With respect to the extrema, find the radius that partitions into equal sized subsets:
	//First store (temporarily) the distances to the extrema somewhere:
	for (auto itr = firstC; itr != end; itr++) {
		//better to use Node temp variable thatn dpivots->push_back(d) followed by dpivots.pop_back();
		(*itr)->sstemp = metric.distance((*itr)->object, (*extrema)->object);
		//(*itr)->sstemp = metric.distance((*itr)->object, (*globalExt)->object);
	}

	//Partiton by distance to the extrema
	//std::nth_element(firstC, median, end,
	//	[](auto left, auto right) {
	//		return left->sstemp < right->sstemp; });
	// median = moveOverlap(firstC, median, end);

	//Partition into two sets by distance to 
	auto minElem = std::min_element(firstC, end, LessThanTemp<N>());
	auto maxElem = std::max_element(firstC, end, LessThanTemp<N>());
	double medDist = 0.5 * ((*maxElem)->sstemp + (*minElem)->sstemp);
	LessThanVal<N> ltm(medDist);
	median = std::partition(firstC, end, ltm);
	//One way to take care of degenerate case:
	if	(median == firstC || median == end) {
		median = firstC + (end - firstC) / 2;
	}


	//And reset node->innerRadius to hold distance to pivot
	for (auto itr = firstC; itr != end; itr++) {
		auto dist = metric.distance((*begin)->object, (*itr)->object);
		(*itr)->dpivots.push_back(dist);
		(*itr)->sstemp = dist;
	}

	calculateDI(begin, firstC, end);

	(*begin)->left = buildTreePCPE(firstC, median);
	(*begin)->right = buildTreePCPE(median, end);

	//Compute and set the DI vector for the begin node;
	computeADIVec(begin, end);

	return *begin;
}

/**
buildTreeAlt()allows for the selection of different pivot selection and partitioning algorithms in 
building the tree; the methods are to be specified inthe constructor.
At time of writing there is no universal "best" set of partitioning and pivot selection
methods, and this routine allows for experimentation with so far..
**/
template < typename N, typename T, typename M>
typename CMTree_Base <N, T, M>::Node*
CMTree_Base < N, T, M >::buildTreeAlt(const NodeItr begin, const NodeItr end) {
	if (begin == end) { return nullptr; }
	if ((end - begin) == 1) {
		(*begin)->setLeaf();
		computeADIVec(begin, end);
		return *begin;
	}

	const NodeItr firstC = begin + 1; //iterator pointing to first child
	NodeItr median = firstC + (end - firstC) / 2; //median of the children

	LessThanLen<N> ltlf;
	auto pivotItr = selectPivot<T, M, NodeItr, LessThanLen<N>>(begin, median, end, pivotType, ltlf, metric);
	std::iter_swap(begin, pivotItr);

	//Note that the distance to tree root is stored at dpivots[0] for any node. 
	for (auto itr = firstC; itr != end; itr++) {
		auto dist = metric.distance((*begin)->object, (*itr)->object);
		(*itr)->dpivots.push_back(dist);
		(*itr)->sstemp = dist;
	}
	calculateDI(begin, firstC, end);

	partition(firstC, median, end);

	(*begin)->left = buildTreeAlt(firstC, median);
	(*begin)->right = buildTreeAlt(median, end);

	//Compute and set the DI vector for the begin node;
	computeADIVec(begin, end);
	return *begin;
}



template < typename N, typename T, typename M>
typename CMTree_Base <N, T, M>::Node*
CMTree_Base < N, T, M >::buildTreeDMR(const NodeItr begin, const NodeItr end) {
	if (begin == end) { return nullptr; }
	if ((end - begin) == 1) {
		(*begin)->setLeaf();
		computeADIVec(begin, end);
		return *begin;
	}

	//Pivot will be a random object
	auto rndIter = findRandomIter< NodeItr>(begin, end);
	std::iter_swap(begin, rndIter);

	const NodeItr firstC = begin + 1; //iterator pointing to first child
	NodeItr median = firstC + (end - firstC) / 2; //median of the children

	/**
	Store the child node distances to the pivot! If a nodeis at level L (where root is at level 0),
	the distance to its immediate parent is stored at dpivots[l-1], and to its grandparent at dpivots[l-2], etc.
	**/
	for (auto itr = firstC; itr != end; itr++) {
		auto dist = metric.distance((*begin)->object, (*itr)->object);
		(*itr)->dpivots.push_back(dist);
		//(*itr)->dpivots.push_front(dist); 
		(*itr)->sstemp = dist;
	}
	calculateDI(begin, firstC, end);

	LessThanTemp<N> ltt;
	auto minElem = std::min_element(firstC, end, ltt );
	auto maxElem = std::max_element(firstC, end, ltt); 
	double medDist = 0.5 * ((*maxElem)->sstemp + (*minElem)->sstemp);
	LessThanVal<N> ltm(medDist);
	median = std::partition(firstC, end, ltm);

	(*begin)->left = buildTreeDMR(firstC, median);
	(*begin)->right = buildTreeDMR(median, end);

	//Compute and set the DI vector for the begin node;
	computeADIVec(begin, end);
	return *begin;
}


//PartType { PIV, EXT, MIN };
template < typename N, typename T, typename M>
void CMTree_Base < N, T, M >::partition(const NodeItr firstC, NodeItr& median, const NodeItr end) {
if (partitionType == PartType::BOM) {
		//Balanced (|LHS|=|RHS|) Object Median Partition.
		std::nth_element(firstC, median, end, LessThanTemp<N>());
	}
	else {
		error("partition function default");
	}
}

template < typename N,  typename T, typename M>
void CMTree_Base < N, T, M >::search(RadiusQuery<T>& q) {
	if (this->root != nullptr) {
		if (doRangeQueryWithPQ == true) {
			//This method is for experimental purposes only
			PQType pq;
			pq.push(pq.newNode(root, nullptr, 0.0, 0.0));
			searchK(pq, q.getTarget(), q, metric);
		}
		else { //This is the standard method
			std::vector <double> distStack;
			searchR(root, q.getTarget(), q, metric, distStack);
		}
	}
}

template <typename N, typename T, typename M>
void CMTree_Base<N, T, M>::search(RadiusCountQuery<T>& q) {
  std::vector<double> distStack;
  searchR(root, q.getTarget(), q, metric, distStack);
}

//SearchCollect should only work with RadiusQuery or RadiusCountQuery
template < typename N,  typename T, typename M>
void CMTree_Base < N, T, M >::searchCollect(SimilarityQuery<T>& q) {
	if (this->root != nullptr) {
		std::vector <double> distStack;
		searchCollect(root, q.getTarget(), q, metric, distStack);
	}
}

template < typename N, typename T, typename M>
void CMTree_Base < N, T, M >::search(NearestKQuery<T>& q) {
	if (this->root != nullptr) {
		PQType pq;
		perfStats.incNodesVisited();
		pq.push(pq.newNode(root, nullptr, 0.0, 0.0));
		searchK(pq, q.getTarget(), q, metric);
	}
}

	/*
		The nearest-k priority que based search.
		Can also be used as radius search unders certain conditions.
	*/
	template < typename N, typename T, typename M>
	void CMTree_Base <N, T, M >::searchK(PQType & queue, T & target, SimilarityQuery<T> & sq, M & met) {
		while (!queue.empty()) {
			PQNodePtr qn = queue.top();
			//perfStats.incNodesVisited();
			queue.pop();
			auto nd = qn->node;
			if (qn->pruningDist <= sq.searchRadius()) {//this is the MAX pruning distance.
				auto dist = met.distance(&target, nd->object);
				perfStats.incDistanceCalls();
				sq.addResult(nd->object, dist);
				qn->distance = dist;
				auto pd = pruningDistance<double>(dist, nd->di.getNear(), nd->di.getFar());									
				if (pd <= sq.searchRadius()) {
					if (nd->left != nullptr) {
						perfStats.incNodesVisited();
						auto maxPruningDist = maxPruningDistance(qn, nd->left);
						if (maxPruningDist <= sq.searchRadius()) {
							queue.push(queue.newNode(nd->left, qn, 0, maxPruningDist));
						}
					}
					if (nd->right != nullptr) {
						perfStats.incNodesVisited();
						auto maxPruningDist = maxPruningDistance(qn, nd->right);
						if (maxPruningDist <= sq.searchRadius()) {
							queue.push(queue.newNode(nd->right, qn, 0, maxPruningDist));
						}
					}
				}
			}
		}
	}
	
	
	

/*
  The maximum prunning distance is the furthest distance that the query object is outside
  the set of node adi. A query radius less than the max prunning distance means that at least one
  of the inequality test will fail (the query object with given radius is too far) and
  the node can be prunned.
  Node: the prunning distance to a leaf node is zero as the leaf node does not have an adi
        other than [0,0]
  cnd - pointer to the current node. The distance between the query object and each of the
	parents is cnd is assumed known.
  pqNodePtr - PQ node pointer to current nodes parent.
*/
	template < typename N, typename T, typename M>
	double CMTree_Base < N, T, M >::maxPruningDistance(PQNodePtr& pqNodePtr, NodePtr& nd) {
		auto ndPtr = pqNodePtr;
		auto ld = 0;  //level difference between nd.adi[0] and node pointed to be pqNodePtr
		auto maxLD = 0;
		double maxPD{ 0 };
		auto nIntervals = nd->adi.size();
		while (ndPtr != nullptr && ld < nIntervals) {
			auto pd = pruningDistance<float>(ndPtr->distance,
				nd->adi[ld ].getNear(),
				nd->adi[ld ].getFar());
			if (pd > maxPD) {
				maxPD = pd;
				maxLD = ld;
			}
			ndPtr = ndPtr->parent;
			++ld;
		}
		++PerfStatsNS::maxPDCtr[maxLD];
		return maxPD;
	}


template < typename N, typename T, typename M>
double CMTree_Base < N, T, M >::maxPruningDistance(NodePtr& nd, std::vector <double>& distStack) {
	auto nLevels = nd->adi.size();
	auto stackSize = distStack.size();
	double maxPD{ 0 };
	for (int lv = 0; lv < nLevels; lv++){
		auto pd = pruningDistance<float>(distStack[stackSize - lv - 1 ],
				nd->adi[lv].getNear(), nd->adi[lv].getFar());
		if (pd > maxPD) { maxPD = pd; }
	}
	return maxPD;
}


template < typename N, typename T, typename M>
double CMTree_Base < N, T, M >::minCollectionDistance(NodePtr& nd, std::vector <double>& distStack) {
	auto nLevels = nd->adi.size();
	auto stackSize = distStack.size();
	auto minCD{  std::numeric_limits<float>::max() };
	for (int lv = 0; lv < nLevels; lv++){
		auto cd = distStack[stackSize - lv - 1 ]  +  nd->adi[lv].getFar();
		if (cd < minCD) { minCD = cd; }
	}
	return minCD;
}

/**
Return true if the node object plus the distance interval adi is completely inside the query volume.
pqDist = distance from the node with adi to the Query center.
radius = the query radius
**/
template < typename N, typename T, typename M>
bool CMTree_Base <N,  T, M >::IntervalInsideQueryVol(double pqDist, double radius, DI& adi) {
	if (pqDist + adi.getFar() <= radius) {
		return true;
	}
	else {
		return false;
	}
}

//Radius (range) search for CMT, without the priority queue
template < typename N, typename T, typename M>
void CMTree_Base < N, T, M >::searchR(NodePtr& nd, T& target, SimilarityQuery<T>& sq, M& met,
	std::vector <double>& distStack) {
	if (nd == nullptr) return;
	perfStats.incNodesVisited();
	auto maxPD = maxPruningDistance(nd, distStack);
	if (maxPD <= sq.searchRadius()) {
		auto dist = met.distance(&target, nd->object);
		perfStats.incDistanceCalls();
		sq.addResult(nd->object, dist);
		auto pd = pruningDistance<double>(dist, nd->di.getNear(), nd->di.getFar());
		if (pd <= sq.searchRadius()) {
			distStack.push_back(dist);
			searchR(nd->left, target, sq, met, distStack);
			searchR(nd->right, target, sq, met, distStack);
			distStack.pop_back();
		}
	}
	return;
}


/*
	Find the near and far member variable values of a node. We asume the distances to parents are stored
	in the value far.
*/
template < typename N, typename T, typename M>
void CMTree_Base < N, T, M >::calculateDI(const NodeItr nd, const NodeItr firstC, const NodeItr end) {
	auto near = std::numeric_limits<float>::max();
	auto far = 0.0;
	for (auto p = firstC; p != end; p++) {
		if ((*p)->sstemp < near) near = (*p)->sstemp;
		if ((*p)->sstemp > far) far = (*p)->sstemp;
	}
	(*nd)->di.setNear(near);
	(*nd)->di.setFar(far + INTERVAL_DX);
}



template < typename N, typename T, typename M>
void CMTree_Base < N, T, M >::searchCollect(NodePtr& nd, T& target, SimilarityQuery<T>& sq, M& met,
	std::vector <double>& distStack) {
	if (nd == nullptr) return;
	perfStats.incNodesVisited();
	auto maxPD = maxPruningDistance(nd, distStack); //maxPD=0 to test CMT-0 with CMT-1 build.
	if (maxPD <= sq.searchRadius()) {
		auto minCD = minCollectionDistance(nd, distStack);
		if(minCD <= sq.searchRadius()){
			collect(nd, sq, met, minCD );
			return;
		}

		auto dist = met.distance(&target, nd->object);
		perfStats.incDistanceCalls();
		if(dist + nd->di.getFar() <=  sq.searchRadius()){
			collect(nd, sq, met, dist + nd->di.getFar() );
		}else{
			sq.addResult(nd->object, dist);
			auto pd = pruningDistance<double>(dist, nd->di.getNear(), nd->di.getFar());
			if (pd <= sq.searchRadius()) {
				distStack.push_back(dist);
				searchCollect(nd->left, target, sq, met, distStack);
				searchCollect(nd->right, target, sq, met, distStack);
				distStack.pop_back();
			}	
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
void CMTree_Base <N, T, M >::collect(NodePtr& nd, SimilarityQuery<T>& sq, M& met, const double dist) {
	if(nd == nullptr){
		return;
	}else{
		sq.addResult(nd->object, dist);
		collect(nd->left,sq, met, dist);
		collect(nd->right, sq, met, dist);
	}
}





template < typename N, typename T, typename M>
double CMTree_Base < N, T, M>::radiusSumAndDepths(std::set<unsigned int>& depths) {
	double sum = 0.0;
	radiusSumAndDepths(root,0, sum, depths);
	return sum;
}

template < typename N, typename T, typename M>
void CMTree_Base <N, T, M >::radiusSumAndDepths(NodePtr& nd, unsigned int depth, double& dsum, std::set<unsigned int>& depths ) {
	if (nd->isLeaf()) {
		depths.insert(depth);
		return;
	}
	if(nd->adi.size() > 0)
		dsum += nd->adi.back().getFar();
	if(nd->left != nullptr)
		radiusSumAndDepths(nd->left,depth+1,dsum, depths);
	if(nd->right!=nullptr)
		radiusSumAndDepths(nd->right,depth+1,dsum,depths);
}


template < typename N, typename T, typename M>
void CMTree_Base <N, T, M >::freePivDistances(NodePtr& nd ) {
	if (nd->isLeaf()) {
		return;
	}
	nd->dpivots.clear();
	nd->dpivots.shrink_to_fit();
	if(nd->left != nullptr)
		freePivDistances(nd->left);
	if(nd->right!=nullptr)
		freePivDistances(nd->right);
}

template < typename N, typename T, typename M>
CMTStructStat <T>  CMTree_Base < N, T, M>::getStructureStats() {
	CMTStructStat<T> s;
	s.processNode(root, 0);
	return s;
}

// Print the adi of each node. Print one node per line, line ordered by node is
template < typename N, typename T, typename M>
void CMTree_Base<N, T, M >::printDI(std::ostream& os) {
	std::map<std::string, NodePtr> nodes;
	collectNodesInMap<NodePtr>(root, nodes);
	for (const auto& pair : nodes) {
		auto nd = pair.second;
		os << nd->object.getId();
		for (const auto& an : nd->adi) {
			os << " " << an;
		}
		os << std::endl;
	}
	os << std::endl;
}


// Collect the objects of evenly spaced leaf nodes
//nPivots is the amount of pivots the user want to find from the leaf nodes.
template < typename N, typename T, typename M>
void CMTree_Base < N, T, M>::collectSpacePivots(const unsigned int nPivots,  std::vector<T>& newPivots, std::vector<T>& newObjects) {
	unsigned int nNodes = nodes.end() - nodes.begin();
	unsigned int nSpace = std::ceil( nNodes / nPivots ); //space to skip between leaf nodes
	unsigned int count = 0;
	collectSpacePivots(root, nSpace, count, newPivots, newObjects);
}

template < typename N, typename T, typename M>
void CMTree_Base <N, T, M >::collectSpacePivots(NodePtr& nd,  const unsigned int nSpace, unsigned int& count, 
	std::vector<T>& newPivots, std::vector<T>& newNodes ) {
	if (nd->isLeaf() ) {
		if (count % nSpace == 0 ){
			newPivots.insert(nd->object);
		}else{
			newNodes.insert(nd->object);
		}
		count++;
		return;
	}else{
		newNodes.insert(nd->object);
	}

	if(nd->left != nullptr)
		collectSpacePivots(nd->left,nSpace, count, newPivots, newNodes);
	if(nd->right!=nullptr)
		collectSpacePivots(nd->right,nSpace, count, newPivots, newNodes);
}

template <typename T, typename M>
class CMTree : public  CMTree_Base<CNode<T>, T, M> {
	using CMTree_Base<CNode<T>, T, M>::CMTree_Base;
	using CMTree_Base<CNode<T>, T, M>::MAX_NDI;
	using CMTree_Base<CNode<T>, T, M>::NodeItr;
public:
	CMTree(std::vector < T >& objects, const M& met, PivotType pivT, PartType partT, unsigned int maxdi = MAX_NDI):
		CMTree_Base<CNode<T>, T, M>::CMTree_Base( objects, met, pivT, partT, maxdi){
	    this->myShortName = "CMT";
	
		std::cout << "Starting CMTree::buildTree() piV parT Nobjects:" << 
		this->pivotType << " " << this->partitionType << " " << objects.size() << std::endl;
		// 
		this->root = this->buildTreeAlt(this->nodes.begin(), this->nodes.end());

		this->freePivDistances (this->root);
	
		std::cout << "Finished CMTree::buildTree()" << std::endl;

		return;
	}

	std::string shortName() const { return this->myShortName; }

};
#endif

