#ifndef EUCLIDIAN_POINT_H
#define EUCLIDIAN_POINT_H

#include <array>

/*
	HM DIstance class.
	Plase note the "using" statement right after the end of the class definition,
	which essentially fixes the dimention and data type for the project.
*/


template <class PT, unsigned int length>
class HMPointG {
private:
	std::string str;
	std::string id;
public:
	HMPointG(const std::string ident, const std::string data) : {
			str= data;
			id - ident;
			//defined in std to overrite every member
	}
	const std::string& getId() const { return id; }
	const std::string& getValue() { return str; }
	const unsigned int dim() { return length; }
	PT operator[](const int i) const { return str[i]; }

	double distance(const EuclidianPointG& p) const {
		auto sum = 0.0;
		for (int i = 0; i < DIM; i++) {
			if(str[i] != p->str[i])
			{
				sum++;
			}
		}
		return sum;
	}

	// friend std::ostream& operator<<(std::ostream& os, EuclidianPointG& n) {
	// 	os << "(" << n.id << ":<";
	// 	os << n.point[0];
	// 	for (int i = 1; i < DIM; i++) {
	// 		os << "," << n.point[i];
	// 	}
	// 	os << ">)";
	// 	return os;
	// }
};
static const unsigned int EuclidianPointDim = 10;
using EuclidianPointPointType = float;
using EuclidianPoint = EuclidianPointG<EuclidianPointPointType, EuclidianPointDim>;
using HMPoint = HMPointG<EuclidianPointPointType, 100>;


#endif // !EUCLIDIAN_POINT_H