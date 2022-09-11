#ifndef PARASAIL_INTERFACE
#define PARASAIL_INTERFACE

#include <string>

#ifdef USE_PARASAIL_LIB

#include "parasail.h"
#include "parasail/matrices/blosum62.h"
//#include "parasail/matrices/pam500.h"
#include "parasail/matrix_lookup.h"

//Interfaces to the Parasail Library Scoring Functions.

double parasailSWS(const std::string & s1, const std::string& s2, const int gapOpen, const int gapExtend) {
	const char *s1c = s1.c_str();
	const char *s2c = s2.c_str();
	int s1Len = (int)s1.length();
	int s2Len = (int)s2.length();
	parasail_result_t *result = NULL;
	const parasail_matrix_t *matrix = NULL;

	result = parasail_sw_scan_32(s1c, s1Len, s2c, s2Len, gapOpen, gapExtend, &parasail_blosum62);
																																														 
	double score = result->score;
	parasail_result_free(result);
	return score;
}

double parasailNWS(const std::string & s1, const std::string& s2, const int gapOpen, const int gapExtend) {
	const char *s1c = s1.c_str();
	const char *s2c = s2.c_str();
	int s1Len = (int)s1.length();
	int s2Len = (int)s2.length();
	parasail_result_t *result = NULL;
	const parasail_matrix_t *matrix = NULL;

	//result = parasail_nw(s1c, s1Len, s2c, s2Len, gapOpen, gapExtend, &parasail_blosum62);//RA
	result = parasail_nw_scan_32(s1c, s1Len, s2c, s2Len, gapOpen, gapExtend, &parasail_blosum62);
	double score = result->score;
	parasail_result_free(result);
	return score;
}
#else

double parasailSWS(const std::string& s1, const std::string& s2, const int gapOpen, const int gapExtend) {
	return 0.0f;
}

double parasailNWS(const std::string& s1, const std::string& s2, const int gapOpen, const int gapExtend) {
	return 0.0f;
}
#endif

/* parasail_result_ssw_t* parasail_ssw(
const char * const restrict s1, const int s1Len,
const char * const restrict s2, const int s2Len,
const int open, const int gap,
const parasail_matrix_t* matrix)
*/


/* typedef struct parasail_result_ssw {
uint16_t score1;
int32_t ref_begin1;
int32_t ref_end1;
int32_t read_begin1;
int32_t read_end1;
uint32_t *cigar;
int32_t cigarLen;
} parasail_result_ssw_t;
*/



//result = parasail_sw_scan_32(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//RA. 7X faster
//result = parasail_sw(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);
//result = parasail_sw_striped_sse41_128_32(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//Slight WA. 7X faster
//result = parasail_sw_scan_32(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//RA. 7X faster
//result = parasail_sw_trace_striped_sat(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//Slight WA. 30X slower
//result = parasail_sw_striped_sse41_128_64(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//Slight WA.  30X SLOWER
//result = parasail_sw_striped_avx2_256_32(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);////WA 10x faster
//result = parasail_sw_trace_striped_avx2_256_32(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//RA 50% SLOWee
//result = parasail_nw(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62); //RA.this nw was about 30% faster than sw
//result = parasail_nw_striped_sse41_128_16(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//WA
//result = parasail_nw_striped_sse41_128_32(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//WA
//result = parasail_nw_striped_avx2_256_64(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);//WA
//result = parasail_nw(s1c, s1Len, s2c, s2Len, 1, 1, &parasail_blosum62);
//double score = result->score;

#endif

