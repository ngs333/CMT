
#ifndef CMT_SEQUENCE
#define CMT_SEQUENCE

#include <string>
#include <vector>
#include "Neighbor.h"

//MZ include headers that implement a archive in simple text format
//#include <boost/archive/tmpdir.hpp>


// Class Sequence is (originally) a representation of Swiss Prot .FA record, which
// includes the record unique identifier, the record value (i.e. the "sequence"), and
// the record comments. Also see class Query
// The final output of one search on a sequence database is the report all of these :
// a) The query information (such as the queryy seqeunce id or possibly a query id)
// b) An ordered list of "neighbors" which includes:
//      i) Neighbor id
//      ii) Distance to Neighbor
//      iii) Alignment Score to Neighbor
//      IV) Alignment string to Neighbor *
// * Alignment string storage is future work, though its may be calculated anyways
//     in the search process.
//TODO:
// a) Add functions to qualify it as C++ Element (Stroustrup)
// b) Place in a Namespace (Sequence, Neighbor and Query are very common words in CS)


class Sequence {
	//friend class boost::serialization::access;
	//friend class Neighbor<Sequence>;
	std::string id;
	std::string value; 
	//double selfScore; //SelfScore is substitution matrix dependent.
	// std::string comments; //Not yet used
	//The output operator is really just printing the metadata
	friend std::ostream& operator<<(std::ostream& os, Sequence& s) {
		os << "[SS:" << s.getId() << "|" << s.getValue().length() << "]";
		return os;
	}
	/*
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & id;
		ar & value;
	}
	*/
public:
	 Sequence(const std::string & sId, const std::string & sVal) : id{ sId }, value{ sVal } {}
	 const std::string & getId() const { return id; }
	 const std::string & getValue() const { return value; }
	//double getSelfScore() { return selfScore; }

	Sequence() {} //used by ar

};


#endif // !CMT_SEQUENCE

