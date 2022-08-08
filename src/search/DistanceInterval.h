#ifndef DISTANCE_INTERVAL
#define DISTANCE_INTERVAL

#include <stack>
#include <vector>
#include <memory>
#include <iostream>
#include <cmath>
#include "Misc.h"
#include "MemoryAux.h"

#define USE_NEXTAFTER 1

template <typename T>
class DistanceInterval {
private:
	T d[2]; //Also possible: T near, far;
public:
	DistanceInterval(const T nd, const T fd) {//Note std::nextafter not used in contructor.
		d[0] = nd; d[1] = fd;
	}
	DistanceInterval() {
		d[0] = 0; d[1] = 0;
	}
	inline void setNear(const T nd) { 
	#ifndef USE_NEXTAFTER
		d[0] = nd; 
	#else
		d[0] = std::nextafter(nd,  0.0f);
	#endif
	}
	inline void setFar(const T nd) { 
	#ifndef USE_NEXTAFTER
		d[1] = nd; 
	#else
		d[1] = std::nextafter(nd,  std::numeric_limits<float>::max());
	#endif
	}
	inline const T getNear() { return d[0]; } 
	inline const T getFar() { return d[1]; }
	inline void reset() { d[0] = 0; d[1] = 0; }
	DistanceInterval& operator = (const DistanceInterval& that) {
		d[0] = that.d[0];
		d[1] = that.d[1];
		return *this;
	}
	inline bool rangeOverlaps(const double pqDistance, const double radius) {
        if ((pqDistance >= getNear() - radius) &&
            (pqDistance <= getFar() + radius)) {
            return true;
        } else {
            return false;
        }
	}
};
/*
	Half open distance intervals for rhs data
*/
template <typename T>
class LeftHDI {
private:
	T d; 
public:
	LeftHDI(const T dist) {
		d = dist;
	}
	LeftHDI() {
		d = 0;
	}
	inline void setFar(const T nd) { 
	#ifndef USE_NEXTAFTER
		d = nd; 
	#else
		d = std::nextafter(nd,  std::numeric_limits<float>::max());
	#endif
	}
	inline const T getFar() { return d; }  
	inline const T getNear() { return 0; }  //used in pruing priority
	inline void reset() { d = 0; }
	LeftHDI& operator = (const LeftHDI& that) {
		d = that.d;
		return *this;
	}
	inline bool rangeOverlaps( const double pqDistance, const double radius){
		if (pqDistance <= getFar() + radius){
			return true;
		}else{
			return false;
		}
	}
};

template <typename T>
class RightHDI {
private:
	T d; 
public:
	RightHDI(const T dist) {
		d = dist;
	}
	RightHDI() {
		d = 0;
	}
	inline void setNear(const T nd) {
	#ifndef USE_NEXTAFTER
		d = nd; 
	#else
		d = std::nextafter(nd, 0.0f);
	#endif
	}
	inline const T getNear() { return d; }
	inline const T getFar() { return  std::numeric_limits<float>::max(); }//used in pruning priority

	inline void reset() { d = 0; }
	RightHDI& operator = (const RightHDI& that) {
		d = that.d;
		return *this;
	}
	inline bool rangeOverlaps( const double pqDistance, const double radius){
			 if (pqDistance + radius >= getNear()){
			return true;
		}else{
			return false;
		}
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

