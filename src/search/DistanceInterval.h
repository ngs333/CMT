#ifndef DISTANCE_INTERVAL
#define DISTANCE_INTERVAL

#include <stack>
#include <vector>
#include <memory>
#include <iostream>
#include "Misc.h"
#include "MemoryAux.h"

template <typename T>
class DistanceInterval {
private:
	T d[2]; //Also possible: T near, far;
public:
	DistanceInterval(T nd, T fd) {
		d[0] = nd; d[1] = fd;
	}
	DistanceInterval() {
		d[0] = 0; d[1] = 0;
	}
	inline void setNear(T nd) { d[0] = nd; }
	inline void setFar(T fd) { d[1] = fd; }
	inline T getNear() { return d[0]; }  //TODO: return T or T&
	inline T getFar() { return d[1]; }
	DistanceInterval& operator = (const DistanceInterval& that) {
		d[0] = that.d[0];
		d[1] = that.d[1];
		return *this;
	}
	
};

//A DistanceInterval with data of instances in one large static rarray:
template <typename T>
class DistanceIntervalSM: public MemoryFSA< DistanceIntervalSM<T>> {
private:
	T d[2]; //Also possible: T near, far;
public:
	DistanceIntervalSM(T nd, T fd) {
		d[0] = nd; d[1] = fd;
	}
	DistanceIntervalSM() {
		d[0] = 0; d[1] = 0;
	}
	inline void setNear(T nd) { d[0] = nd; }
	inline void setFar(T fd) { d[1] = fd; }
	inline T getNear() { return d[0]; }  //TODO: return T or T&
	inline T getFar() { return d[1]; }
	DistanceIntervalSM& operator = (const DistanceIntervalSM& that) {
		d[0] = that.d[0];
		d[1] = that.d[1];
		return *this;
	}
};

#endif

