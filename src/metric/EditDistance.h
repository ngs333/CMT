#include <algorithm>
#include <vector>

/**
Edit distance function implementation obtained from wikibooks.org at:
https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
**/
template<typename T>
typename T::size_type  GeneralizedLevensteinDistance(const T& source, const T& target,
	typename T::size_type insert_cost = 1, typename T::size_type delete_cost = 1, 
	typename T::size_type replace_cost = 1) {
	if (source.size() > target.size()) {
		return GeneralizedLevensteinDistance(target, source, delete_cost, insert_cost, replace_cost);
	}

	using TSizeType = typename T::size_type;
	const TSizeType min_size = source.size(), max_size = target.size();
	std::vector<TSizeType> lev_dist(min_size + 1);

	lev_dist[0] = 0;
	for (TSizeType i = 1; i <= min_size; ++i) {
		lev_dist[i] = lev_dist[i - 1] + delete_cost;
	}

	for (TSizeType j = 1; j <= max_size; ++j) {
		TSizeType previous_diagonal = lev_dist[0], previous_diagonal_save;
		lev_dist[0] += insert_cost;

		for (TSizeType i = 1; i <= min_size; ++i) {
			previous_diagonal_save = lev_dist[i];
			if (source[i - 1] == target[j - 1]) {
				lev_dist[i] = previous_diagonal;
			}
			else {
				lev_dist[i] = std::min(std::min(lev_dist[i - 1] + delete_cost, lev_dist[i] + insert_cost), previous_diagonal + replace_cost);
			}
			previous_diagonal = previous_diagonal_save;
		}
	}

	return lev_dist[min_size];
}
/** 
	from wikibooks.org
	https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
**/
unsigned int edit_distance(const std::string& s1, const std::string& s2)
{
	const std::size_t len1 = s1.size(), len2 = s2.size();
	std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

	d[0][0] = 0;
	for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
	for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

	for (unsigned int i = 1; i <= len1; ++i)
		for (unsigned int j = 1; j <= len2; ++j)
			// note that std::min({arg1, arg2, arg3}) works only in C++11,
			// for C++98 use std::min(std::min(arg1, arg2), arg3)
			d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
	return d[len1][len2];
}

/*
#include <stdio.h>
#include <math.h>
#include <string.h>

//from sanfoundry.com
int d[100][100];
#define MIN(x,y) ((x) < (y) ? (x) : (y))
int main()
{
	int i, j, m, n, temp, tracker;
	char s[] = "Sanfoundry";
	char t[] = "Education";
	m = strlen(s);
	n = strlen(t);

	for (i = 0; i <= m; i++)
		d[0][i] = i;
	for (j = 0; j <= n; j++)
		d[j][0] = j;

	for (j = 1; j <= m; j++)
	{
		for (i = 1; i <= n; i++)
		{
			if (s[i - 1] == t[j - 1])
			{
				tracker = 0;
			}
			else
			{
				tracker = 1;
			}
			temp = MIN((d[i - 1][j] + 1), (d[i][j - 1] + 1));
			d[i][j] = MIN(temp, (d[i - 1][j - 1] + tracker));
		}
	}
	printf("the Levinstein distance is %d\n", d[n][m]);
	return 0;
}
*/
