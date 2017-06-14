#include <cstdlib>
#include <cstdio>
#include <map>
#include <set>
#include <vector>
#include <cassert>

using namespace std;

#define NUM_MASKS 33*33

typedef unsigned long long Item;


Item masks[NUM_MASKS] = {
	//255.255.255.255
	0xFFFFFFFFFFFFFFFFull, //255.255.255.255
	0xFFFFFF00FFFFFFFFull, //255.255.255.0
	0xFFFF0000FFFFFFFFull, //255.255.0.0
	0xFF000000FFFFFFFFull, //255.0.0.0
	0x00000000FFFFFFFFull, //0.0.0.0

	//255.255.255.0
	0xFFFFFFFFFFFFFF00ull, //255.255.255.255
	0xFFFFFF00FFFFFF00ull, //255.255.255.0
	0xFFFF0000FFFFFF00ull, //255.255.0.0
	0xFF000000FFFFFF00ull, //255.0.0.0
	0x00000000FFFFFF00ull, //0.0.0.0

	//255.255.0.0
	0xFFFFFFFFFFFF0000ull, //255.255.255.255
	0xFFFFFF00FFFF0000ull, //255.255.255.0
	0xFFFF0000FFFF0000ull, //255.255.0.0
	0xFF000000FFFF0000ull, //255.0.0.0
	0x00000000FFFF0000ull, //0.0.0.0

	//255.0.0.0
	0xFFFFFFFFFF000000ull, //255.255.255.255
	0xFFFFFF00FF000000ull, //255.255.255.0
	0xFFFF0000FF000000ull, //255.255.0.0
	0xFF000000FF000000ull, //255.0.0.0
	0x00000000FF000000ull, //0.0.0.0

	//0.0.0.0
	0xFFFFFFFF00000000ull, //255.255.255.255
	0xFFFFFF0000000000ull, //255.255.255.0
	0xFFFF000000000000ull, //255.255.0.0
	0xFF00000000000000ull, //255.0.0.0
	0x0000000000000000ull  //0.0.0.0
};

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

int conditionedcount(const map<Item, Counter> * m, const vector<Pair> &p, Item item, int mask) {
	vector<Pair> descendants; //the direct descendants of item,mask
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
	
	map<Item, Counter>::const_iterator ci = m[mask].find(item);
	int s = (ci != m[mask].end() ? ci->second.get() : 0);
	for (vector<Pair>::const_iterator iter = descendants.begin(); iter != descendants.end(); ++iter) {
		map<Item, Counter>::const_iterator ci = m[iter->first].find(iter->second);
		if (ci != m[iter->first].end()) s -= ci->second.get();
	}
	for (vector<Pair>::const_iterator iter = descendants.begin(); iter != descendants.end(); ++iter) {
		for (vector<Pair>::const_iterator it = iter + 1; it != descendants.end(); ++it) {
			//compute glb of iter and it
			if ( (iter->second & masks[it->first]) != (it->second & masks[iter->first]) ) continue;
			int glbmask = 0; while (masks[glbmask] != (masks[it->first] | masks[iter->first])) glbmask++;
			//glb has been computed now check that there is no third parent
			int numparents = 0;
			for (vector<Pair>::const_iterator it3 = descendants.begin(); it3 != descendants.end(); ++it3) {
				if (((masks[it3->first] | masks[glbmask]) == masks[glbmask])
					&& (it3->second == ((it->second | iter->second) & masks[it3->first]))) numparents++;
			}
			if (numparents < 2) {
				printf("%d ", it->first); printip(it->second >> 32); printf(" "); printip(it->second); printf("\n");
				printf("%d ", iter->first); printip(iter->second >> 32); printf(" "); printip(iter->second); printf("\n");
				printf("%d\n", numparents);
				assert(false);
			}
			if (numparents == 2) {
				map<Item, Counter>::const_iterator ci = m[glbmask].find(it->second | iter->second);
				if (ci != m[glbmask].end()) s+= ci->second.get();
			}
		}
	}

	return s;
}

char dummy[] = {'\0'};

int main(int argc, char * argv[]) {
	if (argc != 6 && argc != 7 && argc != 3 && argc != 4) {
		printf("Usage: %s <n> <c> <t> <data> <output> [<exact>]\n", argv[0]);
		printf("Usage: %s <data> <output> [<exact>]\n", argv[0]);
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
		for (int j = 0; j < 8; j++) {
			fscanf(data, "%llu", &x);
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
			printip(it1->second >> 32, exactfile);
			fprintf(exactfile, " ");
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
	double eps=-1;
	double epsil=-1;
	
	vector<Pair> p;
	set<Pair> s;
	for (int i = 0; i < num; i++) {
		int mask;
		fscanf(output, "%d", &mask);
		Item item = 0;
		Item w, x, y, z;
		fscanf(output, "%llu.%llu.%llu.%llu", &w, &x, &y, &z);
		item = ((item << 8)|w);
		item = ((item << 8)|x);
		item = ((item << 8)|y);
		item = ((item << 8)|z);
		fscanf(output, "%llu.%llu.%llu.%llu", &w, &x, &y, &z);
		item = ((item << 8)|w);
		item = ((item << 8)|x);
		item = ((item << 8)|y);
		item = ((item << 8)|z);
		int fmin, fmax;
		fscanf(output, "%d%d", &fmin, &fmax);

		int truecount = m[mask][item].get();
		if (truecount < fmin - n/counters || truecount > fmax + n/counters) {
			accErrors++;
			if (getenv("NONVERBOSE") == NULL){
				printf("Accuracy: [%d] ", mask);
				printip(item >> 32);
				printf(" ");
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
