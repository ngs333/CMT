#ifndef NEIGHBOR_H
#define NEIGHBOR_H
// TODO:
// a) Add functions to qualify it as C++ Element (Stroustrup)
template <class T>
class Neighbor {
public:
	std::string id;
	double distance;
	double score;
	//std::string alignment;
	friend std::ostream& operator<<(std::ostream& os, const Neighbor& n) {
		os << "[NS:" << n.id << "|" << n.distance << "|" << n.score << "]";
		return os;
	}
	explicit Neighbor(const T * seq, double ds, double sc) :
		id{ seq->getId() }, distance{ ds }, score{ sc } {}
	explicit Neighbor(const T& seq, double ds, double sc) :
		id{ seq.id }, distance{ ds }, score{ sc } {}
	//copy constructor 
	Neighbor(const Neighbor & n) :
		id{ n.id }, distance{ n.distance }, score{ n.score } {}
	//copy assignment
	Neighbor & operator=(const Neighbor & n) {
		id = n.id;
		distance = n.distance;
		score = n.score;
	}
};
#endif
