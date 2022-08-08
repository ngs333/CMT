#ifndef TREE_NODES
#define TREE_NODES


#include <vector>
#include "DistanceInterval.h"
#include "MemoryAux.h"

#ifdef SA_USE_STATIC_DI
using DI = DistanceIntervalSM<float>; //Static memory version
#else
using DI = DistanceInterval<float>;
#endif // SA_USE_STATIC_DI


//Cascading Metric tree nodes:
template <class T>
class CNode {
public:
	T* object; 
	CNode* left;
	CNode* right;
	int size;
	FSVector<DI> adi;//the ancestral distance interrval
	//std::vector<DI> adi;
	DI   diL, diR; //Normal child distance intervals

	//std::vector<float> dpivots;//Note that sizeof(std::vector<float>) =24 bytes on a 64 bit Windows 
	std::deque<float> dpivots; //TODO: fix above
	float sstemp;
	CNode(T* object) : object(object) {}
	virtual ~CNode() {} //DI should be cleared if not static; pivods cleared in Tree::buildTree().
	bool isLeaf() { return ((left == nullptr) && (right == nullptr)); }
	void setLeaf() {
		left = nullptr;  right = nullptr;
		diL.setNear(0); diL.setFar(0); diR.setNear(0); diR.setFar(0);
	}
	float getNear(){
		if (left != nullptr){
			return diL.getNear();
		}else{
			return diR.getNear();
		}
	}
	// get the furthest from among both left and right.	
	float getFar(){
		if (right != nullptr){
			return diR.getFar();
		}else{
			return diL.getFar();
		}
	}

	/*
		Return true if the query overlaps the left
		distace interval (i.e if the distance(p,q) is not too
		far to the right of the left interval)).
	*/
	bool rangeOverlapsLeft( const double pqDistance, const double radius){
		if(left == nullptr){
			return false;
		}else if (pqDistance <= diL.getFar() + radius){
			return true;
		}else{
			return false;
		}
	}
	/*
		Return true if the query overlaps the right
		distace interval (i.e if the distance(p,q) is not too
		far to the left of the right interval)).
	*/
	bool rangeOverlapsRight( const double pqDistance, const double radius){
		if(right == nullptr){
			return false;
		}else if (pqDistance + radius >= diR.getNear()){
			return true;
		}else{
			return false;
		}
	}

	friend std::ostream& operator<<(std::ostream& os, CNode& nd) {
		os << "(" << nd.object << ")";
		return os;
	}
};

//Metric tree nodes.
template <class T>
class MNode {
public:
	T* object;           // Point or object selected as pivot - one per node.
	int size;
	DistanceInterval<float> di; 
	float sstemp;//Used a temp variable and aslo for sef score after tree is done.
	MNode* left;  //I.e. left child, etc.
	MNode* right;
	MNode(T* object) : object{ object } {}
	bool isLeaf() { return ((left == nullptr) && (right == nullptr)); }
	void setLeaf() { left = nullptr;  right = nullptr; }
};

//Metric tree nodes.
template <class T>
class ANode {
public:
	T* object;           // Point or object selected as pivot - one per node.
	int size;
#ifdef USE_HALF_INTERVALS
	LeftHDI<float> diL;
	RightHDI<float> diR;
#else
	DI diL, diR;   //Left and right distance intervals, distances relative to pivot.
#endif

	float sstemp;//Used a temp variable and aslo for sef score after tree is done.
	ANode* left;  //I.e. left child, etc.
	ANode* right;
	ANode(T* object) : object{ object } {}
	bool isLeaf() { return ((left == nullptr) && (right == nullptr)); }
	void setLeaf() { left = nullptr;  right = nullptr; }
	// get the nearest from among both left and right.
	float getNear(){
		if (left != nullptr){
			return diL.getNear();
		}else{
			return diR.getNear();
		}
	}
	// get the furthest from among both left and right.	
	float getFar(){
		if (right != nullptr){
			return diR.getFar();
		}else{
			return diL.getFar();
		}
	}
	 bool rangeOverlaps(const double pqDistance, const double radius) {
		if (((left != nullptr) && (diL.rangeOverlaps( pqDistance, radius) == true)) ||
			((right != nullptr) && (diR.rangeOverlaps( pqDistance, radius) == true)))
			return true;
		else
			return false;
	 }

	/*
		Return true if the query overlaps the left
		distace interval (i.e if the distance(p,q) is not too
		far to the right of the left interval)). Note this
		coniders the left distance interval di = [near, far] == [0, far]
	*/
	/*bool rangeOverlapsLeft( const double pqDistance, const double radius){
		if(left == nullptr){
			return false;
		}else if (pqDistance <= diL.getFar() + radius){
			return true;
		}else{
			return false;
		}
	}
	*/
	/*
		Return true if the query overlaps the right
		distace interval (i.e if the distance(p,q) is not too
		far to the left of the right interval)).  Note this
		coniders the right distance interval di = [near, far] == [near, infinity]
	*/
	/*
	bool rangeOverlapsRight( const double pqDistance, const double radius){
		if(right == nullptr){
			return false;
		}else if (pqDistance + radius >= diR.getNear()){
			return true;
		}else{
			return false;
		}
    }
	*/

    bool rangeOverlaps(DistanceInterval<float> & di, const double pqDistance, const double radius) {
        if ((pqDistance >= di.getNear() - radius) &&
            (pqDistance <= di.getFar() + radius)) {
            return true;
        } else {
            return false;
        }
    }
};

/**
 * Class for storing and getting the CMT strucure statistics.
*  Can be called once the tree is build.
*  Warning, for LCMT, nPivots calculated as below does not make sense and it is too
*  small by a factor of two.
*  For unbalanced trees, depths will have more than one or two elements, and usually we are
*  interested in the min depth and the max depth.
*/
template <class T>
class CMTStructStat {
protected:
	unsigned int nNodes;  //number of tree nodes
	std::set<unsigned int> depths ;  
	double dSum;
	unsigned int nAdis;
	unsigned int nPivots;
public:
	friend std::ostream& operator<<(std::ostream& os, CMTStructStat& s) {
			os <<s.nNodes<<","<<*(s.depths.cbegin()) <<","<< *(s.depths.crbegin()) <<","
			<<s.nPivots<<","<< (1.0 * s.nPivots)/s.nNodes << ","<< s.nAdis<<","<< (1.0 * s.nAdis)/s.nNodes << ","
			<<s.dSum;
		return os;
	}
	CMTStructStat(){
		dSum = 0;
		nAdis = 0;
		nPivots = 0;
		depths.clear();
		nNodes = 0;
	}

	void processNode(CNode<T>* nd, unsigned int depth){
		if(nd == nullptr){
			return;
		}
		nNodes++;
		nPivots += nd->dpivots.size();
		nAdis += nd->adi.size();
		if (nd->isLeaf()) {
			depths.insert(depth);
			return;
		}
		if(nd->adi.size() > 0)
			dSum += nd->adi.back().getFar();
		processNode(nd->left, depth + 1);
		processNode(nd->right, depth + 1);
	}
};

/***
	Following are several comparators that use the
	"sstemp" node memeber variable.
***/
template <class Nd>
class LessThanTemp {
	public:
		LessThanTemp(){ }
		bool operator()(const Nd* left, const Nd* right) const {
			return left->sstemp < right->sstemp;
		}
	};

template <class Nd>
class LessThanLen {
public:
	LessThanLen() { }
	bool operator()(const Nd* left, const Nd* right) const {
		return left->object->getValue().size() < right->object->getValue().size();
	}
};

template <class Nd>
class LessThanVal {
		double val;
	public:
		LessThanVal(double m) : val { m }{ }
		bool operator()(const Nd* p) const {
			return (p->sstemp < val); //TODO < or <=
		}
	};

template <class Nd>
class LessThanEqVal {
	double val;
public:
	LessThanEqVal(double m) : val{ m } { }
	bool operator()(const Nd* p) const {
		return (p->sstemp <= val); //TODO < or <=
	}
};

template <class Nd>
class GreaterThanEqTo {
	double val;
public:
	GreaterThanEqTo(double m) : val{ m } { }
	bool operator()(const Nd* p) const {
		return (p->sstemp >= val);
	}
};

template <class Nd>
class EqualByVal {
	double val;
public:
	EqualByVal(double m) : val{ m } { }
	bool operator()(const Nd* left) const {
		return (left->sstemp == val); //TOD: By pctError
	}
};

template <class Nd>
class PivDistLessThanAtD {
	unsigned int depth;
public:
	PivDistLessThanAtD(unsigned int d) : depth{ d } { }
	bool operator()(const Nd* left, const Nd* right) const {
		return left->dpivots[depth] < right->dpivots[depth];
	}
};

template <class Nd>
class PivDistLessThanEqValAtD {
	double val;
	unsigned int depth;
public:
	PivDistLessThanEqValAtD(double v, unsigned int d) : val{ v }, depth{ d } { }
	bool operator()(const Nd* p) const {
		return (p->dpivots[depth] <= val); //TODO < or <=
	}
};

template <class Nd>
class PivDistLessThanValAtD {
	double val;
	unsigned int depth;
public:
	PivDistLessThanValAtD(double v, unsigned int d) : val{ v }, depth{ d } { }
	bool operator()(const Nd* p) const {
		return (p->dpivots[depth] < val); //TODO < or <=
	}
};



template <class Nd>
class PivDistEqValAtD {
	double val;
	unsigned int depth;
public:
	PivDistEqValAtD(double v, unsigned int d) : val{ v }, depth{ d } { }
	bool operator()(const Nd* p) const {
		return (p->dpivots[depth] == val); //TODO < or <=
	}
};


template <class Nd>
class PivDistLessThanAtDCyc {
	unsigned int depth;
	unsigned int maxd;
public:
	PivDistLessThanAtDCyc(unsigned int d, unsigned int md) : depth{ d }, maxd{ md }{ }
	bool operator()(const Nd* left, const Nd* right) const {
		return left->dpivots[depth % maxd] < right->dpivots[depth % maxd];
	}
};

template <class Nd>
class PivDistLessThanValAtDCyc {
	double val;
	unsigned int depth;
	unsigned int maxd;
public:
	PivDistLessThanValAtDCyc(double v, unsigned int d, unsigned int md ): val{ v }, depth{ d }, maxd{ md } { }
	bool operator()(const Nd* p) const {
		return (p->dpivots[depth % maxd] < val); //TODO < or <=
	}
};

template <class Nd>
class PivDistLessThanEqValAtDCyc {
	double val;
	unsigned int depth;
	unsigned int maxd;
public:
	PivDistLessThanEqValAtDCyc(double v, unsigned int d, unsigned int md ): val{ v }, depth{ d }, maxd{ md } { }
	bool operator()(const Nd* p) const {
		return (p->dpivots[depth % maxd] <= val); //TODO < or <=
	}
};

/*
class LessThanLen {
public:
	LessThanLen() { }
	bool operator()(const NodePtr& left, const NodePtr& right) const {
		return (left->object.getValue().size() < right->object.getValue().size());
	}
};
*/

#endif
