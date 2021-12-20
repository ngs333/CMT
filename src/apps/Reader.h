#ifndef READER_FASTA
#define READER_FASTA

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <regex>

#include "Sequence.h"
#include "Misc.h"

using namespace std;

/*
Some functions for reading and partitioning FASTA file
*/
class FastaReader {
public:
	/*
	Get and return the sequence identifier from a line assumed to start
	with the new record deliminator '>'
	*/
	std::string extractIdentifier(const std::string& line)
	{
		std::string id = line.substr(0, line.find_first_of(' '));
		if (id.length() < 1) {
			throw std::invalid_argument("extractIdentifier cant parse: " + line);
		}
		return id;
	}

	unsigned int  countFastaSeqs(const std::string& filename)
	{
		static constexpr char fasta_delim{ '>' };
		unsigned int counter = 0;
		std::string line;
		ifstream ifs(filename.c_str());
		if (ifs.is_open())
		{
			while (getline(ifs, line))
			{
				if (line.size() == 0)continue;
				if (line.front() == fasta_delim) {//line.find( ... )
					counter++;
				}
			}
		}
		ifs.close();
		return counter;
	}

	std::vector<Sequence> 
		read(const std::string& fname, const unsigned int maxSeqs) {
		std::vector<Sequence> nullEV;
		filesystem::path pathObj(fname);
		if (filesystem::exists(pathObj)) {
			if (filesystem::is_regular_file(pathObj)) {
				return readFile(fname, maxSeqs);
			}
			else if (filesystem::is_directory(pathObj)) {
				return readDirectory(fname, maxSeqs);
			}
			else {
				error(fname, "is not fasta file or direcotry");
				return nullEV;
			}
		}
		else {
			error(fname, "does not exists");
			return nullEV;
		}
	}

	std::vector<Sequence> readDirectory(const std::string& fname,  const unsigned int maxSeqs) {
		std::vector<Sequence> seqs;
		for (const auto& de : std::filesystem::recursive_directory_iterator(fname)) {
			std::regex regex_text("(.*)(\.fasta|\.fna| \.faa)$");
			if (regex_match(de.path().string(), regex_text)) {
				auto seqsFromOne = read(de.path().string(), maxSeqs - seqs.size());
				seqs.insert(seqs.end(), seqsFromOne.begin(), seqsFromOne.end());
			}
		}
		return seqs;
	}

	std::vector<Sequence> readFile(const std::string& fname, const unsigned int maxSeqs) {
		static constexpr char fasta_delim{ '>' };
		//static constexpr char fastq_delim{ '@' };
		std::ifstream file;
		std::vector<Sequence> result;

		//auto nLines = countFastaSeqs(fname);

		file.open(fname);
		if (!file.is_open()) {
			error("Unable to open file  ", fname);
		}

		std::string seqId;
		std::string line;
		std::string seqVal;
		auto resultsStart = result.size();;
		while (std::getline(file, line) &&
			((result.size() - resultsStart) < maxSeqs)) {
			if (line.size() == 0) continue;
			if (line.front() == fasta_delim) {
				seqId = extractIdentifier(line);
				if (seqVal.empty() != true) {
					result.push_back(Sequence(seqId, seqVal));
					seqVal.clear();
				}
			}
			else {
				seqVal.append(line);
			}
		}
		if (seqVal.empty() != true) {
			result.push_back(Sequence(seqId, seqVal));
		}

		file.close();
		return result;
	}

	std::vector<std::vector<Sequence>*> partitionStrings(std::vector<Sequence> v) {
		double binSize = 0;
		for (auto s : v) {
			binSize += s.getValue().length();
		}

		binSize /= 10;
		std::cout << "Bin Size = " << binSize << std::endl;

		std::sort(v.begin(), v.end(), [](Sequence& r, Sequence& l) {return r.getValue().size() < l.getValue().size(); });

		std::vector<std::vector<Sequence>*> partitions;
		std::vector<Sequence>* bin = new std::vector<Sequence>();
		partitions.push_back(bin);
		auto binBytes = 0;
		auto minLen = v[0].getValue().length();
		auto priorLength = minLen;
		for (auto s : v) {
			auto sLength = s.getValue().length();
			if ((((binBytes + sLength) >= binSize) && (sLength != priorLength)) ||
				((minLen > 100) && (sLength > 2 * minLen))) {
				//done with the current bin; start a new one
				bin = new std::vector<Sequence>();
				partitions.push_back(bin);
				bin->push_back(s);
				binBytes = s.getValue().length();
				minLen = s.getValue().length();
			}
			else {
				binBytes += s.getValue().length();
				bin->push_back(s);
				//maxs = sLength;
				priorLength = s.getValue().length();
			}
		}

		//print out some binnin stats:
		std::cout << "\nminLen maxLen  maxLen/minLen bin.size()  " << std::endl;
		for (auto bin : partitions) {
			auto minLen = 10000000000;
			auto maxLen = 0;
			for (auto s : *bin) {
				auto sLength = s.getValue().length();
				if (sLength <= minLen) { minLen = sLength; }
				if (sLength >= maxLen) { maxLen = sLength; }
			}
			std::cout << minLen << " " << maxLen << " " << (1.0 * maxLen / minLen) << " " << bin->size() << std::endl;
		}

		return partitions;
	}

	/*
		Read nofReads sequences from a file and extract sampleSize seqences from it.
		nofReads neds to be k * sampleSize, wiith k>= 1
		*/
	void sampleSubsetAndCopy(const std::string& fname, const unsigned int sampleSize, const unsigned int nofReads) {
		static constexpr char fasta_delim{ '>' };  // '@'
		std::ifstream file;
		std::ofstream ofile;
		std::set<std::string> uSeqIds;//unique string seqeunce ids
		std::mt19937 gen(3001); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_int_distribution<> dis(1, nofReads);

		auto sizeStr = to_string(sampleSize);
		string ofname = fname.substr(0, fname.find(".fasta"));
		ofname = ofname + "_" + sizeStr + ".fasta";

		file.open(fname);
		ofile.open(ofname);
		if (!file.is_open()) {
			error("Unable to open file  ", fname);
		}
		if (!ofile.is_open()) {
			error("Unable to open file  ", ofname);
		}

		if (nofReads < sampleSize) {
			error("nofReads < sampleSilze");
		}


		std::string seqId;
		std::string line;
		std::string seqVal;
		unsigned int sCount = 0;
		bool writeSeq = false;
		while (std::getline(file, line) && (sCount <= sampleSize)) {
			if (line.front() == fasta_delim) {
				auto ranNo = dis(gen);
				auto id = extractIdentifier(line);
				if ((ranNo <= sampleSize) && (uSeqIds.find(id) == uSeqIds.end())) {
					++sCount;
					uSeqIds.insert(line);
					writeSeq = true;
					if (sCount > sampleSize) {//but not the line beyond the last one
						writeSeq = false;
					}
				}
				else {
					writeSeq = false;
				}
			}
			if (writeSeq == true) {
				ofile << line << endl;
			}
		}

		file.close();
		ofile.close();
	}
};

/*
	A reader of simple text file assuming one word per line.
	Read up to max words. Ignores words with numbers and blank chars.
*/
class SimpleTextReader {
public:
	std::vector<Sequence> read(const std::string& fname, const unsigned int maxWords) {
		std::ifstream file;
		std::vector<Sequence> result;

		file.open(fname);
		if (!file.is_open()) {
			error("Unable to open file  ", fname);
		}

		std::ostringstream idStream;
		unsigned int id = 1;
		std::string line;
		auto resultsStart = result.size();;
		while (std::getline(file, line) && ((result.size() - resultsStart) < maxWords)) {
			if (line.find_first_of("1234567890 ") == std::string::npos) {
				idStream << id;
				result.push_back(Sequence(idStream.str(), line));
				idStream.str(std::string());
				idStream.clear();
				++id;
			}
		}

		file.close();
		return result;
	}
};

/*
	Creates a and uses a SimpleTextReader if the file name (path) argument ends in .txt,
	otherwise it creates and uses a Fasta reader (for files or directories)

*/
class Reader {
public:
	std::vector<Sequence> read(const std::string& fname, const unsigned int maxWords) {
		if (fname.find(".txt") != std::string::npos) {
			SimpleTextReader reader;
			return reader.read(fname, 1000000);
		}
		else {
			FastaReader reader;
			return  reader.read(fname, maxWords);
		}
	}
};

void updateCharFrequency(const std::string& str, unordered_map<char, unsigned int> & cmap){
	for (std::string::size_type i = 0; i< str.size(); ++i) {
		if (cmap.find(str[i]) == cmap.end()) {
			cmap.insert(make_pair(str[i], 1));
		}else {
			cmap[str[i]]++;
		}
	}
}

void printFrequency(const std::vector<Sequence>& seqs, ostream& stream) {
	unordered_map<char, unsigned int> cmap;
	for (const auto& s : seqs) {
		updateCharFrequency(s.getValue(), cmap);
	}

	for (auto& it : cmap) {
		stream << it.first << ' ' << it.second << '\n';
	}
}


#endif
