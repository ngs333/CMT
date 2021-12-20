#ifndef PERF_STATS
#define PERF_STATS

class PerfStats {
private:
	long long int nodesVisited;
	long long int distanceCalls;
public:
	PerfStats() : nodesVisited{ 0 }, distanceCalls{ 0 } {}
	long long int getNodesVisited() { return nodesVisited; }
	long long int getDistanceCalls() { return distanceCalls; }
	void incNodesVisited() { nodesVisited++; }
	void incDistanceCalls() { distanceCalls++; }
	void incNodesVisited(long long int v) { nodesVisited += v;  }
	void incDistanceCalls(long long int v) { distanceCalls += v; }

	void addStats( PerfStats & s) {
		nodesVisited += s.getNodesVisited();
		distanceCalls += s.getDistanceCalls();
	}
	void reset() {
		nodesVisited = 0;
		distanceCalls = 0;
	}
};
#endif // !PERF_STATS

