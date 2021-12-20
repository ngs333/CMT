#ifndef FASTA_TEST_CHARS
#define FASTA_TEST_CHARS

#include <random>
#include <string>
#include <algorithm>

class CharacterSets {
protected:
	inline static const char genomicCharSet[] = "ACTGU";
	inline static const char proteonomicCharSet[] = "ACDEFGHIKLMNOPQRSTUVWY";
	inline static const char englishCharSet[] = "abcdefghijklmnopqrstuvwxyz";

	class RandChar{
	public:
		unsigned int maxIdx;
		const char* & charSet;
		RandChar(unsigned int n, const char*& cs) : maxIdx{ n }, charSet{ cs }{}
			 char operator()(){ 
					return charSet[rand() % maxIdx]; 
			 }
	};

	static std::string randomString(unsigned int length, const char* charSet, unsigned int maxId) {
		std::string str(length, 0);
		std::generate_n(str.begin(), length, RandChar(maxId, charSet));
		return str;
	}
public:
	static std::string randomEnglish(unsigned int length)
	{
		auto maxId = sizeof(englishCharSet) - 1;
		return randomString(length, englishCharSet, maxId);
	}


	static std::string randomGenomic(size_t length)
	{
		auto maxId = sizeof(genomicCharSet) - 1;
		return randomString(length, genomicCharSet, maxId);
	}

	static std::string randomProtein(size_t length)
	{
		auto maxId = sizeof(proteonomicCharSet) - 1;
		return randomString(length, proteonomicCharSet, maxId);
	}
};

#endif