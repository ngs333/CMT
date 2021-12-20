#ifndef EUCLIDIAN_POINT_H
#define EUCLIDIAN_POINT_H

#include <array>

/*
	Euclidian point class.
	Plase note the "using" statement right after the end of the class definition,
	which essentially fixes the dimention and data type for the project.
*/


template <class PT, unsigned int DIM>
class EuclidianPointG {
private:
	std::array<PT, DIM> point;
	std::string id;
public:
	EuclidianPointG(const std::string ident, const std::array<PT, DIM>&  pData) : id{ ident } {
			point= pData;//defined in std to overrite every member
	}
	const std::string& getId() const { return id; }
	const std::array<PT, DIM>& getValue() { return point; }
	const unsigned int dim() { return DIM; }
	PT operator[](const int i) const { return point[i]; }

	double distance(const EuclidianPointG& p) const {
		auto sum = 0.0;
		for (int i = 0; i < DIM; i++) {
			sum += (point[i] - p.point[i]) * (point[i] - p.point[i]);
		}
		return std::sqrt(sum);
	}
	bool isInsideBox(std::array<PT, DIM> const& o, std::array<PT, DIM>  const& p) {
		for (int i = 0; i < DIM; i++) {
			if ((point[i] < o[i]) || (point[i] > p[i])) {
				return false;
			}
		}
		return true;
	}
	friend std::ostream& operator<<(std::ostream& os, EuclidianPointG& n) {
		os << "(" << n.id << ":<";
		os << n.point[0];
		for (int i = 1; i < DIM; i++) {
			os << "," << n.point[i];
		}
		os << ">)";
		return os;
	}
};
static const unsigned int EuclidianPointDim = 10;
using EuclidianPointPointType = float;
using EuclidianPoint = EuclidianPointG<EuclidianPointPointType, EuclidianPointDim>;


#endif // !EUCLIDIAN_POINT_H