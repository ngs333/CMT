5#ifndef SW_METRIC_TESTS
#define SW_METRIC_TESTS

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>

#include "Misc.h"
#include "Sequence.h"
#include "Metric.h"

extern void checkMetrics(const std::vector<std::string> & gsv, const std::string & fileName, Metric<Sequence>& mps);
extern void testTrigInequality(const std::vector<std::string>  & gsv, std::ofstream & theFile, Metric<Sequence>& mps);
extern void testMetricSymmetry(const std::vector<std::string>  & gsv, std::ofstream & theFile, Metric<Sequence>& mps);
extern void checkTwoMetricDiff(const std::vector<std::string>  & gsv, std::ofstream & theFile, Metric<Sequence>& mps);

void testTrigInequality( const std::vector<Sequence>  & gsv, std::ofstream & theFile,Metric<Sequence>& met ) {
	std::clock_t start = std::clock();
	unsigned int npassed = 0;
	unsigned int nfailed = 0;
	double maxdd = 0;
	double dd;
	double maxss = 0.0;

	for (auto i = 0; i < gsv.size(); i++) {
		//cout << "I=" << i << endl;
		for (auto j = 0; j < gsv.size(); j++) {
			unsigned int sk = rand() % (gsv.size() - 1);
			double sdik = met.distance (gsv[i], gsv[sk]);
			double sdkj = met.distance(gsv[sk], gsv[j]);
			double sdij = met.distance(gsv[i], gsv[j]);
			if (sdik > maxss) { maxss = sdik; }
			if (sdkj > maxss) { maxss = sdkj; }
			if (sdij > maxss) { maxss = sdij; }

			if (sdij > sdik + sdkj) {
				//cout << " Metric Error " << endl;
				nfailed++;
				dd = sdij - (sdik + sdkj);
				if (dd > maxdd) {
					maxdd = dd;
				}
			}
			else {
				//cout << "sidk= " << sdik << endl;
				npassed++;
			}
		}
	}
	unsigned int swtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	theFile << "MetricTest1:\n"
		<< "npassed=" << npassed << " nfailed=" << nfailed << " maxdd= " << maxdd
		<< " maxss= " << maxss << " time=" << swtime << std::endl;

}

void testMetricSymmetry(const std::vector<std::string>  & strs, 
	std::ofstream & theFile, Metric<Sequence>& met) {

	std::vector<Sequence> gsv;
	int counter = 1;
	for (auto elem : strs) {
		gsv.push_back(Sequence(NumberToString(counter), elem));
			counter++;
	}
	std::clock_t start = std::clock();
	unsigned int npassed = 0;
	unsigned int nfailed = 0;

	for (auto i = 0; i < gsv.size(); i++) {
		for (auto j = 0; j < gsv.size(); j++) {
			double d1 = met.distance(gsv[i], gsv[j]);
			double d2 = met.distance(gsv[j], gsv[i]);

			if (std::abs(d1 - d2) > 0.0000001) {
				std::cout << " Metric Error " << std::endl;
				nfailed++;
			}
			else {
				npassed++;
			}
		}
	}
	unsigned int swtime =   dTimeSeconds(start);
	theFile << "MetricTest2:\n"
		<< "npassed=" << npassed << " " << " nfailed=" << nfailed << " time=" << std::endl;

	//ofstream theFile("swruns.txt", ios::app);
	//theFile << "\nSWtest: maxNofStrings maxStringSize nfailed npassed time" << endl;
	//theFile << maxNofStrings << ", " << maxStringSize << ", " << nfailed
	//<< ", " << npassed << ", " << swtime << endl;
	///theFile.close();
}

void checkTwoMetricDiff(const std::vector<std::string>  & strs, 
	std::ofstream & theFile, Metric<Sequence>& met) {

	std::clock_t start = std::clock();
	unsigned int npassed = 0;
	unsigned int nfailed = 0;

	std::vector<Sequence> gsv;
	int counter = 1;
	for (auto elem : strs) {
		gsv.push_back(Sequence(NumberToString(counter), elem));
		counter++;
	}

	for (int i = 0; i < gsv.size(); i++) {
		for (int j = 0; j < gsv.size(); j++) {
			double d1 = met.distance(gsv[i], gsv[j]);
			double d2 = met.distance(gsv[i], gsv[j]);
			if (std::abs(d1 - d2) > 2) {
				//std::cout << " Metric Error " << std::endl;
				nfailed++;
			}
			else {
				npassed++;
			}
		}
	}
	unsigned int swtime = dTimeSeconds(start);
	theFile << "MetricTest3:\n"
		 << "npassed=" << npassed << " " << " nfailed=" << nfailed << " time=" << std::endl;


}
void checkMetrics(const std::vector<Sequence>& gsv, const std::string & fileName,
	Metric<Sequence>& mps) {
	std::ofstream theFile(fileName, std::ios::app);
	theFile << "\ncheckMetrics" << std::endl;
	testTrigInequality( gsv, theFile, mps);
	//metricTest2(gsv, theFile);
	//metricTest3(gsv, theFile);
	theFile.close();


}

#endif // !SW_METRIC_TESTS