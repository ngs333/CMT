
#include <array>
#include "gedLib/ged.h"


class GED {
private:
	std::string ChEBIid;
public:
	GED(const std::string ident) {
			ChEBIid = ident;
	}
	const std::string& getId() const { return ChEBIid; }

	long long distance(GED p) const {
		return getDistance(p.ChEBIid,ChEBIid);
	}

};

