#ifndef SEARCH_PQ
#define SEARCH_PQ

#include <queue>
#include <vector>
#include <algorithm>

/*
   The priority queue used by the metric trees in the search function.
   Also self-manages the PQNode memory. 
   This was origianlly desinged ans an alternative to using std::shared_pointer<PQNode>
   as shared pointers were too slow.
*/

template <class Node>
class SPQNode {
public:
	Node* node;
	SPQNode* parent;
	float distance;
	float pruningDist;
	SPQNode() {}
	explicit SPQNode(Node* nd, SPQNode* pp, float ds, float pd) :
		node{ nd }, parent{ pp }, distance{ ds }, pruningDist{ pd }{
		error("PQNode ctor");
	}
	void reset(Node* nd, SPQNode* pp, float ds, float pd) {
		node = nd;
		parent = pp;
		distance = ds;
		pruningDist = pd;
	}
};


template <class Node, class PQNode, class Comp>
class SearchPQ {
protected:
	inline static unsigned int NODES_PER_ARRAY{ 512 };
	//using PQNode = SPQNode<Node>;
	//Comparator for the priority queue
	/*struct CompareResultPQ {
		bool operator()(const PQNode* lhs, const PQNode* rhs) const {
			return lhs->pruningDist > rhs->pruningDist;
		}
	};
	*/

	unsigned int CURRENT_ARRAY;
	unsigned int used;
	std::priority_queue<PQNode*, std::vector<PQNode*>, Comp> pq;
	std::vector<PQNode*> dataVec;
public:
	SearchPQ() {
		PQNode* nd = new PQNode[NODES_PER_ARRAY];
		dataVec.push_back(nd);
		CURRENT_ARRAY = 0; //0 is the first array
		used = 0;
	}
	~SearchPQ() {
		for (const auto& na : dataVec) {
			delete[] na;
		}
		dataVec.clear();
	}

	//TODO: Modify design so that memory used can grow.
	PQNode* newNode(Node* nd, PQNode* pp, float ds, float pd) {
		PQNode* ed = nullptr;
		if (used >= NODES_PER_ARRAY) {
			//Get the next array of nodes.
			if ((CURRENT_ARRAY + 1) < dataVec.size()) {
				++CURRENT_ARRAY;
				used = 0;

			}
			else {
				PQNode* nd = new PQNode[NODES_PER_ARRAY];
				dataVec.push_back(nd);
				++CURRENT_ARRAY;
				used = 0;
			}
		}
		ed = dataVec[CURRENT_ARRAY];
		ed += used;
		++used;
		ed->reset(nd, pp, ds, pd);
		return ed;
	}

	bool empty() {
		return pq.empty();
	}
	PQNode* top() {
		return pq.top();
	}

	void pop() {
		pq.pop();
	}

	void push(PQNode* ptr) {
		pq.push(ptr);
	}

};

#endif
