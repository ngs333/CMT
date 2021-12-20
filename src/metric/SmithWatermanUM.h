#ifndef JU_SMITH_WATERMAN
#define JU_SMITH_WATERMAN

//Appendix A
//Univeristy of Misssouri CS Smith-Waterman Distance Function
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>


double smithWatermanUM(const std::string & _s1, const std::string & _s2,
		unsigned int gap, unsigned int mismatch) {
	auto s1 = _s1;
	auto s2 = _s2;
	int l1 = s1.length();
	int l2 = s2.length();
	int i, j, compareValue;
	int returnValue = 0;
	//Allocating the space for matrix
	double ** r = (double **)malloc(sizeof(double *) * l1);
	for (i = 0; i < l1; i++) {
		r[i] = (double *)malloc(sizeof(double) * l2);
	}
	
	
	//initialization
	for (i = 0; i < l1; i++) {
		r[i][0] = gap * i;
	}
	for (i = 1; i < l2; i++) {
		r[0][i] = gap * i;
	}
	//filling out the matrix
	for (i = 1; i < l1; i++) {
		for (j = 1; j < l2; j++) {
			;  //BUG?
				//if Match
				if (s1[i - 1] == s2[j - 1]) {
					compareValue = r[i - 1][j - 1];
				}
				else { //Its a Mis-match
					compareValue = r[i - 1][j - 1] + mismatch;
				}
				//check if the cell on the left to the position is greater than the cell value
				if (compareValue > (r[i][j - 1] + gap)) {
					compareValue = r[i][j - 1] + gap;
				}
				//check if the cell on the top to the position is greater than the cell value
				if (compareValue > (r[i - 1][j] + gap)) {
					compareValue = r[i - 1][j] + gap;
				}
				//set current cell to highest possible score and move ahead
				r[i][j] = compareValue;
		}
	}

	returnValue = r[l1 - 1][l2 - 1];
	for (i = 0; i < l1; i++) {
		// std::cout << i << std::endl;
		free(r[i]);
	}
	free(r);
	return returnValue;
}

#endif 