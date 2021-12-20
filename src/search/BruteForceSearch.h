#ifndef BRUTE_FORCE_S_H
#define BRUTE_FORCE_S_H


#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>

#include "Sequence.h"
#include "Query.h"
#include "PerfStats.h"

template < typename T, typename M>
class BruteForceSearch {
private:
	PerfStats perfStats;
	M  metric;
	std::vector < T > & points;
public:
	BruteForceSearch(std::vector <T>& pv, const M& met) :
		points{ pv }, metric{ met }{}
	~BruteForceSearch() {};
	PerfStats& getPerfStats() { return perfStats; }

	void search(SimilarityQuery<T> & sq) {
		auto target = sq.getTarget();
		for (auto& point : points) {
			//TODO: If CMT is converting to floats before recording result, so must BF to agree.
			//double dist = metric.distance(target, point);
			float dist = metric.distance(target, point);
			sq.addResult(&point, dist);
		}
	}
};

#endif // !BRUTE_FORCE_S_H
