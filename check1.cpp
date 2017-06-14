#include <cstdlib>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

using namespace std;

#ifndef NUM_MASKS
#define NUM_MASKS 5
#endif

typedef unsigned int Item;

#if NUM_MASKS==5
Item masks[NUM_MASKS] = {
	0xFFFFFFFFu,
	0xFFFFFF00u,
	0xFFFF0000u,
	0xFF000000u,
	0x00000000u
};
#endif
#if NUM_MASKS==33
Item masks[NUM_MASKS] = {
	0xFFFFFFFFu << 0,
	0xFFFFFFFFu << 1,
	0xFFFFFFFFu << 2,
	0xFFFFFFFFu << 3,
	0xFFFFFFFFu << 4,
	0xFFFFFFFFu << 5,
	0xFFFFFFFFu << 6,
	0xFFFFFFFFu << 7,
	0xFFFFFFFFu << 8,
	0xFFFFFFFFu << 9,
	0xFFFFFFFFu << 10,
	0xFFFFFFFFu << 11,
	0xFFFFFFFFu << 12,
	0xFFFFFFFFu << 13,
	0xFFFFFFFFu << 14,
	0xFFFFFFFFu << 15,
	0xFFFFFFFFu << 16,
	0xFFFFFFFFu << 17,
	0xFFFFFFFFu << 18,
	0xFFFFFFFFu << 19,
	0xFFFFFFFFu << 20,
	0xFFFFFFFFu << 21,
	0xFFFFFFFFu << 22,
	0xFFFFFFFFu << 23,
	0xFFFFFFFFu << 24,
	0xFFFFFFFFu << 25,
	0xFFFFFFFFu << 26,
	0xFFFFFFFFu << 27,
	0xFFFFFFFFu << 28,
	0xFFFFFFFFu << 29,
	0xFFFFFFFFu << 30,
	0xFFFFFFFFu << 31,
	0x00000000u
};
#endif

void printip(Item item, FILE * fp = stdout) {
	fprintf(fp, "%d.%d.%d.%d", (int)(255&(item >> 24)), (int)(255&(item >> 16)), (int)(255&(item >> 8)), (int)(255&(item >> 0)));
}

class Counter {
	int x;
public:
	Counter() : x(0) {}
	void inc() {x++;}
	int get() const {return x;}
};

typedef pair<int, Item> Pair;

int conditionedcount(map<Item, Counter> * m, const vector<Pair> &p, Item item, int mask) {
	vector<Pair> descendants; //the direct descendants of item,mask
	//int countcandidates=0;
	for (vector<Pair>::const_iterator it = p.begin(); it != p.end(); ++it) {
		if ((masks[it->first] & masks[mask]) == masks[mask] && (it->second & masks[mask]) == item) {
			//it is a descendant, now eliminate non-direct descendants
			vector<Pair> newdescendants;
			for (vector<Pair>::const_iterator iter = descendants.begin(); iter != descendants.end(); ++iter) {
				if ((masks[iter->first] & masks[it->first]) != masks[it->first] || (iter->second & masks[it->first]) != it->second) {
					newdescendants.push_back(*iter);
				}
			}
			descendants = newdescendants;
			descendants.push_back(*it);
		}
	}
	
	int s = m[mask][item].get();
	for (vector<Pair>::iterator iter = descendants.begin(); iter != descendants.end(); ++iter) {
		s -= m[iter->first][iter->second].get();
	}

	return s;
}

char dummy[] = {'\0'};

int main(int argc, char * argv[]) {
	if (argc != 6 && argc != 7 && argc != 3 && argc != 4) {
		printf("Usage: %s <n> <c> <t> <data> <output> [<exact>]\n", argv[0]);
		printf("Alternative usage: %s <data> <output> [<exact>]\n", argv[0]);
		printf("\t<n>=number of IP pairs\n");
		printf("\t<c>=counters (1/eps)\n");
		printf("\t<data>=data file (IP pairs)\n");
		printf("\t<output>=output of HHH algorithm to check\n");
		printf("\t<exact>=file to output exact HHH\n");
		return 0;
	}
	int n, counters, threshold, num;
	float time=-1; char * algname = dummy;
	FILE * data, * output, * exactfile;
	if (argc > 4) {
		n = atoi(argv[1]);
		counters = atoi(argv[2]);
		threshold = atoi(argv[3]);
		data = fopen(argv[4], "r");
		output = fopen(argv[5], "r");
		exactfile = (argc == 7 ? fopen(argv[6], "w") : NULL);
		fscanf(output, "%d", &num);
	} else {
		data = fopen(argv[1], "r");
		output = fopen(argv[2], "r");
		exactfile = (argc == 4 ? fopen(argv[3], "w") : NULL);
		fscanf(output, "%as %d %d %d %d %f\n", &algname, &n, &counters, &threshold, &num, &time);
	}

	if (getenv("DONTCHECKABOVE") != NULL && atoi(getenv("DONTCHECKABOVE")) < n) {
		fclose(data);
		fclose(output);
		if (exactfile != NULL) fclose(exactfile);
		return 0;
	}
	FILE * fsummary = (getenv("SUMMARYFILE") != NULL ? fopen(getenv("SUMMARYFILE"), "a") : NULL);
	map<Item, Counter> m[NUM_MASKS];
	for (int i = 0; i < n; i++) {
		Item x;
		Item item = 0;
		for (int j = 0; j < 4; j++) {
			fscanf(data, "%u", &x);
			item = ((item << 8) | x);
		}
		for (int j = 0; j < NUM_MASKS; j++) {
			m[j][item & masks[j]].inc();
		}
	}	

	vector<Pair> exact;
	vector<int> exactcount;
	for (int i = 0; i < NUM_MASKS; i++) {
		for (map<Item, Counter>::const_iterator it = m[i].begin(); it != m[i].end(); ++it) {
			if (it->second.get() >= threshold && conditionedcount(m, exact, it->first, i) >= threshold) {
				exact.push_back(Pair(i, it->first));
				exactcount.push_back(it->second.get());
			}
		}
	}
	if (exactfile != NULL) {
		fprintf(exactfile, "%d\n", (int)exact.size()); 
		vector< Pair >::const_iterator it1 = exact.begin();
		vector< int >::const_iterator it2 = exactcount.begin();
		while (it1 != exact.end())  {
			fprintf(exactfile, "%d ", it1->first);
			printip(it1->second, exactfile);
			fprintf(exactfile, " %d %d\n", *it2, *it2);
			++it1;
			++it2;
		}
	}
	if (fsummary != NULL && argc > 4) {
		fprintf(fsummary, "algorithm=%s nitems=%d counters=%d threshold=%d outputsize=%d\n",
		                   argv[0],     n,        counters,   threshold,   (int)exact.size()); 
	}

	int accErrors=0;
	int covErrors=0;

	vector<Pair> p;
	set<Pair> s;
	double eps = -1; //the true accuracy
	double epsil = -1;

	for (int i = 0; i < num; i++) {
		int mask;
		fscanf(output, "%d", &mask);
		Item item = 0;
		Item w, x, y, z;
		fscanf(output, "%u.%u.%u.%u", &w, &x, &y, &z);
		item = ((item << 8)|w);
		item = ((item << 8)|x);
		item = ((item << 8)|y);
		item = ((item << 8)|z);
		if (getenv("NONVERBOSE") == NULL){
			printip(item);
			printf("---\n");
		}
		int fmin, fmax;
		fscanf(output, "%d%d", &fmin, &fmax);;

		int truecount = m[mask][item].get();
		if (truecount < fmin - n/counters || truecount > fmax + n/counters) {
			accErrors++;
			if (getenv("NONVERBOSE") == NULL){
				printf("Accuracy: [%d] ", mask);
				printip(item);
				printf(" [%d] %d %d\n", truecount, fmin, fmax);
			}
		}

		double epscand = ((double) (fmax-fmin)) / truecount;
		if (epscand > eps) eps = epscand;
		epsil = max(epsil, ((double)(fmax-fmin))/n);

		p.push_back(Pair(mask, item));
		s.insert(Pair(mask, item));
	}

	for (int i = 0; i < NUM_MASKS; i++) {
		for (map<Item, Counter>::const_iterator it = m[i].begin(); it != m[i].end(); ++it) {
			if (it->second.get() >= threshold && conditionedcount(m, p, it->first, i) >= threshold && 0==s.count(Pair(it->first, i))) {
				covErrors++;
				if (getenv("NONVERBOSE") == NULL){
					printf("Coverage: [%d] ", i);
					printip(it->first);
					printf(" [%d/%d]\n", it->second.get(), conditionedcount(m, p, it->first, i));
				}
			}
		}
	}

	if (accErrors+covErrors > 0 && fsummary != NULL) {
		fprintf(fsummary, "Error: algorithm=%s nitems=%d counters=%d threshold=%d accErrors=%d covErrors=%d\n",
			                  argv[0],     n,        counters,   threshold, accErrors, covErrors); 
	}

	if (fsummary != NULL && argc <= 4) {
		fprintf(fsummary, "check=true  algorithm=%-14s nitems=%-10d counters=%-5d threshold=%-9d outputsize=%-3d exactsize=%-3d time=%f accuracy=%lf epsilon=%lf\n",
		                   algname,        n,           counters,     threshold,     num,        (int)exact.size(), time,   eps, epsil); 
	}

	if (exactfile != NULL) fclose(exactfile);
	if (fsummary != NULL) fclose(fsummary);
	fclose(data);
	fclose(output);
	if (argc <= 4) free(algname);

	printf("%d accErrors %d covErrors %d exact hhhs\n", accErrors, covErrors, (int)exact.size());

	return 0;
}
