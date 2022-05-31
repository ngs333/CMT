#ifndef METRIC_PARAMS
#define METRIC_PARAMS

#include <ostream>
#include <cstdlib>
#include <cmath>

#include "Sequence.h"
#include "ParasailInterface.h"
#include "SmithWatermanUM.h"
#include "Misc.h"
#include "EditDistance.h"
#include "edlib.h"
#include "HM.h"
#include "ged.h"

template <class T>
class Metric {
public:
	virtual double distance(const T& a, const T& b) = 0;
	virtual double distance(const T* a, const T* b) = 0;
	virtual ~Metric() {};
};
template <class T>
class ScoreMetric : public Metric<T>  {
public:
	virtual double score(const T& a, const T& b) = 0;
	//score from distance :
	virtual double score(const T& a, const T& b, const double dist) = 0;
	virtual unsigned int minGapCost(unsigned int ql, unsigned int minL, unsigned int maxL) = 0;
	virtual ~ScoreMetric() {};
};

class EuclidianMetric : public  Metric<EuclidianPoint> {
public:
	double distance(const EuclidianPoint& p1, const EuclidianPoint& p2) {
		return p1.distance(p2);
	}
	double distance(const EuclidianPoint* p1, const EuclidianPoint* p2) {
		return distance(*p1, *p2);
	}
	friend std::ostream& operator <<(std::ostream& os, EuclidianMetric& mp) {
		os << typeid(mp).name() << ", ";
		return os;
	}
};

class EditDistMetric : public  Metric<Sequence>{
public:
	double distance(const Sequence& p1, const Sequence& p2) {

		auto result = edlibAlign(p1.getValue().c_str(), p1.getValue().length(),
			p2.getValue().c_str(), p2.getValue().length(), edlibDefaultAlignConfig());
		if (result.status != EDLIB_STATUS_OK) {
			error("EDLIB error");
		}
		auto dist =  result.editDistance;
		edlibFreeAlignResult(result);
		return dist;
		
		//return GeneralizedLevensteinDistance(p1.getValue(), p2.getValue());
	}
	double distance(const Sequence* p1, const Sequence* p2) {
		return distance(*p1, *p2);
	}

	friend std::ostream & operator <<(std::ostream & os, EditDistMetric  mp) {
		os << typeid(mp).name() << ", ";
		return os;
	}
};
// using HMPoint = HMPointG<float, 100>;

class HMMetric : public  Metric<HMPoint> {
public:
	double distance(const HMPoint& p1, const HMPoint& p2) {
		return p1.distance(p2);
	}
	double distance(const HMPoint* p1, const HMPoint* p2) {
		return distance(*p1, *p2);
	}
};

class GEDMetric : public  Metric<GED> {
public:
	double distance(const GED& p1, const GED& p2) {
		return p1.distance(p2);
	}
	double distance(const GED* p1, const GED* p2) {
		return distance(*p1, *p2);
	}
};

#endif
