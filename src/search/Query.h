#ifndef CMT_QUERY
#define CMT_QUERY

#include <string>
#include <vector>
#include <queue>
#include <set>

#include "Sequence.h"
#include "Neighbor.h"
#include "Metric.h"


//TODO: Untill the class PQWith Access can be debugged:
//#define USE_ACCESS_DERIVED_PQ = 1

using SequenceNeighbor = Neighbor<Sequence>;


/*
	This priority queue class is used because multiset is too slow.
	It was derived from the std::priority_queue in order to have access
	the the "c" memeber variable - which is a vector. 
*/
#ifdef USE_ACCESS_DERIVED_PQ
 template <class CT, class Container = std::vector<CT>,
 	 class Compare = std::less<typename Container::value_type> > 
  class  PQWithAcess : public std::priority_queue<CT,Container,Compare> {
	protected:
		Container & getDataVector() {
			 return c;
		}
};
#endif

// SimilarityQuery is the interface used by the metric search data structures.
// TODO: Better name it metric query ?
template <class T>
class SimilarityQuery {
public:
	virtual double searchRadius() = 0;
	virtual void   addResult(T * r, double distance) = 0;
	virtual T &   getTarget() = 0; //The query sequence record
};

// A class where the query search radius is simply the metric radius
// and THE RESULTS ARE NOT STORED -  just the count of the results are maintained.
// Used purely for certain computaional complexity testing at large radii.
template <class T>
class RadiusCountQuery : public SimilarityQuery<T> {
protected:
	T   target;
	double     radius;
	int        count; // the count of results
public:
	RadiusCountQuery(const T& si, double rad):
		target{ si }, radius{ rad }, count{0} {}

	double searchRadius() { 
		return radius;
	}
	//double searchRadius(double maxSS) { return radius; }

	T & getTarget() { return target; }

	void setSearchRadius(double r) { radius = r; }

	void   addResult(T * obj, double dist) {
		if (dist > searchRadius()) {
			return;
		}else{
			count++;
		}
	}	

	int getCount(){ return count;}	
};

// A class where the query search radius is simply the metric radius
template <class T>
class RadiusQuery :  public SimilarityQuery<T> {
protected:
	T   target;
	double     radius;
	int        maxResults;
	bool    softMaxResults;//If true, allow for more than maxResults if last ones have same value as the maxResult one..
	struct distanceLessThan
	{
		bool operator()(const Neighbor<T>& lhs, const Neighbor<T>& rhs) const{
			return lhs.distance < rhs.distance;
		}
	};
	

#ifdef  USE_ACCESS_DERIVED_PQ
	//std::priority_queue<Neighbor<T>, std::vector<Neighbor<T>>, distanceLessThan> neighbors2;
	PQWithAccess <Neighbor<T>, std::vector<Neighbor<T>>, distanceLessThan> neighbors;
#else
	std::multiset<Neighbor<T>, distanceLessThan> neighbors;
#endif

public:
	RadiusQuery(const T& si, double rad, const int mr = 10000000) :
		target{ si }, radius{ rad }, maxResults{ mr }, softMaxResults{ true } {}

	double searchRadius() { 
		return radius;
	}
	double searchRadius(double maxSS) { return radius; }

	T & getTarget() { return target; }

	void setSoftMaxResults(bool soft) {
		softMaxResults = soft;
	}

	void setSearchRadius(double r) { radius = r; }

#ifdef  USE_ACCESS_DERIVED_PQ
	std::vector<Neighbor<T>> & getNeighbors() {
	return neighbors.getDataVector();
	}
#else
	std::multiset<Neighbor<T>, distanceLessThan> &  getNeighbors() {
		return neighbors;
	}
#endif

	/**
		Add to the result set. Note if maxResults is exceeded, the worst on the
		 result set is deleted - even if it has the same distance as the next worst.
		 NOTE: We are  verifying neigbor.distance <= radius

	**/
	void   addResult(T * seq, double dist) {
		const double score = 0.0; //Score is not used.
		if (dist > searchRadius()) {
			return;
		}
		
		neighbors.insert(Neighbor<T>(seq, dist, score));
			
#ifdef USE_ACCESS_DERIVED_PQ
		if (neighbors.size() > maxResults) {
			 if (softMaxResults == false) {
				while(neighbors.size() > maxResults){
					Neighbor n = neighbors.pop();
				}
			}else {
				// delete only those larger than the nth;
				while(neighbors.size() > maxResults){
					Neighbor n = neighbors.pop();
					Neighbor np = neightbors.top();
					if (n.distance == np.distance){
						neighbors.insert(n);
					break;
					}
				}
			}
		}
				
#else
		if (neighbors.size() > maxResults) {
			auto ptr = neighbors.begin();
			if (softMaxResults == false) {
				//delete all after the nth - even if they have same value as nth.
				std::advance(ptr, maxResults);
				neighbors.erase(ptr , neighbors.end());
			}else {
				// delete only those larger than the nth;
				std::advance(ptr, maxResults - 1);
				neighbors.erase(neighbors.upper_bound(*ptr), neighbors.end());
			}
		}
#endif
	}

	bool hasSameNeighbors(RadiusQuery & qB){
		if (neighbors.size() != qB.getNeighbors().size()) {
			return false;
		}
		std::multiset<std::string> idsB;
		for (const auto & n : qB.getNeighbors()) {
			idsB.insert(n.id);
		}
		for (const auto & n : neighbors) {
			if (idsB.find(n.id) == idsB.end()) {
				return false;
			}
		}
		return true;
	}

	/*
		Return the neighbors of qB that are not present in the this object.
	*/
	auto missingNeighbors(RadiusQuery& qB) {
		std::multiset<Neighbor<T>, distanceLessThan> missing;
		std::multiset<std::string> ids;
		for (const auto& n : neighbors) {
			ids.insert(n.id);
		}
		for (const auto& n : qB.getNeigbors()) {
			if (ids.find(n.id) == ids.end()) {
				missing.insert(n);
			}
		}
		return missing;
	}
};




template <class T>
class NearestKQuery :  public RadiusQuery<T> {
public:
	//Constructor: Note the order of max results and radius is reverssed vs. Radius query
	NearestKQuery(const T& si, const int mr, double rad = std::numeric_limits<double>::max()) 
		: RadiusQuery<T> ::RadiusQuery(si,  rad, mr ){}

	
	double searchRadius() override {
		if (this->neighbors.size() < this->maxResults) {
			return this->radius; //Assuming all results added are less then radius.
		}
		else {
			auto ptr = this->neighbors.crbegin();
			double mrad = ptr->distance;
			//double mrad = neighbors.top().distance;
			return mrad;
		}
	}
	//double searchRadius(double maxSS) override {return searchRadius();}
};

//A class where the search radius depends upon the score function in addition to the metric.
//Stores results based on both high scores and high disctance.
//Note that the same object may be reused on a different search tree(and thereby adding to the previously 
//stored results) by resetting the maxSelfScore.
class SequenceScoreQuery : public SimilarityQuery<Sequence> {
private:
	Sequence & target;
	double     selfScore; //the self score of the member or "query sequence"
	double     minScore;// Minimum acceptable score, or p times selfScore allowed in the results found.
	int        maxResults;//The maximum number of results to save. 
	double     worstScore;//The worst score among the results saved. Zero if maxResults
							//have not yet been found.

	ScoreMetric<Sequence> * metric;
	unsigned int  minSequenceLength;//the min self score of the sequences in the SET being searched
	unsigned int  maxSequenceLength;//the maximum self score of the sequences in the set being searched
	double        minSelfScore;//the min self score of the sequences in the SET being searched.
	double        maxSelfScore;//the maximum self score of the sequences in the SET being searched.
	
	struct scoreGreaterThan {
		bool operator()(const SequenceNeighbor& lhs, const SequenceNeighbor& rhs) const
		{
			return lhs.score > rhs.score;
		}
	};
	struct distanceLessThan {
		bool operator()(const SequenceNeighbor& lhs, const SequenceNeighbor& rhs) const
		{
			return lhs.distance < rhs.distance;
		}
	};
	// The sets used to store the neighbors (or results) found by the query. At most
	//  maxResults will be stored in them.
	std::multiset<SequenceNeighbor, scoreGreaterThan> scoreNeighbors;
	std::multiset<SequenceNeighbor, distanceLessThan> distanceNeighbors;
public:
	SequenceScoreQuery(Sequence & si, double pss, double ss, const int mr = 10000000) :
		target{ si }, minScore(pss * ss), selfScore{ ss }, maxResults{ mr } {
		metric = nullptr;
		maxSelfScore = 0.0; //this must be reset before use be new partition search
		worstScore = 0.0;
	}

	//The search radius is in part based on an inequality derived from the function:
	//     D_ij = S_ii + S_jj - 2 * S_ij. 
	//The search radius in this case depends on:
	//		1) The self score of the current target
	//      2) The maximum self score in the database or bin being searched.
	//      3) The gap extension costs (possibly configered in the DB during DB creation)
	//      4) Other parameters set by the user for this query
	double searchRadius(double maxSS) {
		double radius;
		if (worstScore > minScore) {
			radius = selfScore + maxSS - 2 * worstScore;
		}
		else {
			radius = selfScore + maxSS - 2 * minScore;
		}
		double offset = metric->minGapCost(target.getValue().length(),
			minSequenceLength, maxSequenceLength);
		radius -= offset;
		return radius;
	}
	double searchRadius() {
		return searchRadius(maxSelfScore);
	}
	Sequence & getTarget() { return target; }
	std::multiset<SequenceNeighbor, distanceLessThan> & getDistanceNeighbors() {
		return distanceNeighbors;
	}

	std::multiset<SequenceNeighbor, scoreGreaterThan> &  getScoreNeighbors() {
		return scoreNeighbors;
	}

	double getMinScore() {return minScore;}

	void setMinSelfScore(double minSS) { minSelfScore = minSS; }
	void setMaxSelfScore(double maxSS) { maxSelfScore = maxSS; }
	void setMaxSequenceLength(double maxSL) { maxSequenceLength = maxSL; }
	void setMinSequenceLength(double minSL) { minSequenceLength = minSL; }
	void setMetricParameters(ScoreMetric<Sequence> & mp) { metric = &mp; }

	bool hasSameNeighbors(SequenceScoreQuery & qB) {
		if (scoreNeighbors.size() != qB.getScoreNeighbors().size()) {
			return false;
		}
		std::multiset<std::string> idsB;
		for (const auto & n : qB.getScoreNeighbors()) {
			idsB.insert(n.id);
		}
		for (const auto & n : scoreNeighbors) {
			if (idsB.find(n.id) == idsB.end()) {
				//printNeighbors(getScoreNeighbors(), qB.getScoreNeighbors());
				// int ii; std::cin >> ii;
				return false;
			}
		}
		return true;
	}


	void printNeighbors(std::multiset<SequenceNeighbor, scoreGreaterThan>& na, std::multiset<SequenceNeighbor, scoreGreaterThan>& nb) {
		std::cout << std::endl << "Neighbot set A:" << std::endl;
		for (auto n : na) {
			std::cout << n << std::endl;
		}
		std::cout << std::endl << "Neighbot set B:" << std::endl;
		for (auto n : nb) {
			std::cout << n << std::endl;
		}
		std::cout << "END NS AB " << std::endl;
	}



/**
	Add the sequence to the two neighbors list
	The score is not know by the search trees, so calculate it here.
**/
	void   addResult(Sequence * seq, double dist) {
		if (dist > searchRadius()) {
			return;
		}
		//score from pre-calculated distance :
		auto pairScore = metric->score(target, *seq, dist);
		if (pairScore  < minScore) {
			return;
		}

		Neighbor<Sequence> r(seq, dist, pairScore); //Note neighbor mmust be a C++ element

		distanceNeighbors.insert(r);
		if (distanceNeighbors.size() > maxResults) {
			auto ptr = distanceNeighbors.begin();
			// Soft_max_results case: delete only those larger than the nth;
			std::advance(ptr, maxResults - 1);
			distanceNeighbors.erase(distanceNeighbors.upper_bound(*ptr), distanceNeighbors.end());
		}

		scoreNeighbors.insert(r);
		if (scoreNeighbors.size() > maxResults) {
			auto ptr = scoreNeighbors.begin();
			// Soft_max_results case: delete only those larger than the nth;
			std::advance(ptr, maxResults - 1);
			scoreNeighbors.erase(scoreNeighbors.upper_bound(*ptr), scoreNeighbors.end());
		
		}

		//Save the value worst score already found into tempScore;
		double tempScore = (scoreNeighbors.crbegin())->score;
		if ((distanceNeighbors.size() >= maxResults) && 
			(tempScore > worstScore)) {
			worstScore = tempScore;
		}
	}

};
#endif // !CMT_QUERY
