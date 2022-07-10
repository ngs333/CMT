#ifndef CMT_PARTITION_FUNCTION
#define CMT_PARTITION_FUNCTION

#include <algorithm>
#include <vector>

#include "SearchCommon.h"
#include "TreeNodes.h"

// Abstract base class
template < typename N, typename  NodeItr, typename T, typename M>                                                                                                                                                                                                 
class PartitionFunction {
public:
  PartitionFunction() {};
  virtual void operator() (const NodeItr firstC, NodeItr& median, const NodeItr end) = 0;
};

// Add two doubles 
template < typename N, typename  NodeItr, typename T, typename M>                                                                                                                                                                                                     
class BOM : public PartitionFunction<N,NodeItr,T,M> {
public:
  BOM() {};
  virtual void operator() (const NodeItr firstC, NodeItr& median, const NodeItr end) { 
    std::nth_element(firstC, median, end, LessThanTemp<N>());
    return;
     }
};

template < typename N, typename  NodeItr, typename T, typename M>
class DMR : public PartitionFunction<N,NodeItr,T,M>{
public:
    DMR(){};
    virtual void operator()(const NodeItr firstC, NodeItr &median, const NodeItr end){
        auto minElem = std::min_element(firstC, end, LessThanTemp<N>());
        auto maxElem = std::max_element(firstC, end, LessThanTemp<N>());
        double medDist = 0.5 * ((*maxElem)->sstemp + (*minElem)->sstemp);

        LessThanVal<N> ltm(medDist);
        median = std::stable_partition(firstC, end, ltm);
        // One way to take care of degenerate case:
        if (median == firstC || median == end)
        {
            median = firstC + (end - firstC) / 2;
        }
        return;
    }
};


template < typename N, typename  NodeItr, typename T, typename M>
PartitionFunction<N,NodeItr,T,M> * getPartitionFunctor(PartType pType){
  PartitionFunction<N,NodeItr,T,M>  * ftor = nullptr;
    if(pType == PartType::BOM){
      ftor = new BOM<N,NodeItr,T,M>();
    }else if(pType == PartType::DMR){
      ftor = new DMR<N,NodeItr,T,M>();
    }else if(pType == PartType::EXT){
      ftor = new BOM<N,NodeItr,T,M>();
    }else{
      std::cout <<"Error in getPartitionFunctor"<<std::endl;
      exit(-1);
    }
    return ftor;
}


template<class T, class Node, class  NodeItr, class Comparator>
inline std::tuple<T,T> getNearFar( const NodeItr begin, const NodeItr end, Comparator& compare){
  auto nearFar = std::minmax_element(begin, end, compare);
         //[] (Node const* lhs, Node const* rhs) {return lhs->sstemp < rhs->sstemp;});
	T near = (*nearFar.first)->sstemp;
  T far = (*nearFar.second)->sstemp;
  return {near, far};
}


/*
  Determine and set the bouding distance interval di from the
  node->sstemp variable in the container of nodes, using nodes in [begin,end)
*/
template <class Node, class NodeItr, class Comparator>
inline void calculateDBI(DI & di, const NodeItr begin, const NodeItr end, Comparator& compare){
  auto [near, far] = getNearFar<float, Node, NodeItr>(begin, end, compare);
  di.setNear(near);
  di.setFar(far);
}

#endif
