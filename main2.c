/*
The main function for the two-dimensional HHH program
- Ran Ben Basat, June 2017.
Based on an earlier version:
- Thomas Steinke (tsteinke@seas.harvard.edu) 2010-11
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>	
#include "randhhh2D.h"
#include <sys/timeb.h>

#ifndef NRRUNS
#define NRRUNS 5
#endif

#ifndef VMULT
#define VMULT 1
#endif

#if VMULT>1
#define PROB
#endif



#ifndef CLK_PER_SEC
#ifdef CLOCKS_PER_SEC
#define CLK_PER_SEC CLOCKS_PER_SEC
#endif
#endif

//the masks
LCLitem_t masks[NUM_COUNTERS] = {
	//255.255.255.255
	//0-4
	0xFFFFFFFFFFFFFFFFull, //255.255.255.255
	0xFFFFFF00FFFFFFFFull, //255.255.255.0
	0xFFFF0000FFFFFFFFull, //255.255.0.0
	0xFF000000FFFFFFFFull, //255.0.0.0
	0x00000000FFFFFFFFull, //0.0.0.0

	//255.255.255.0
	//5-9
	0xFFFFFFFFFFFFFF00ull, //255.255.255.255
	0xFFFFFF00FFFFFF00ull, //255.255.255.0
	0xFFFF0000FFFFFF00ull, //255.255.0.0
	0xFF000000FFFFFF00ull, //255.0.0.0
	0x00000000FFFFFF00ull, //0.0.0.0

	//255.255.0.0
	//10-14
	0xFFFFFFFFFFFF0000ull, //255.255.255.255
	0xFFFFFF00FFFF0000ull, //255.255.255.0
	0xFFFF0000FFFF0000ull, //255.255.0.0
	0xFF000000FFFF0000ull, //255.0.0.0
	0x00000000FFFF0000ull, //0.0.0.0

	//255.0.0.0
	//15-19
	0xFFFFFFFFFF000000ull, //255.255.255.255
	0xFFFFFF00FF000000ull, //255.255.255.0
	0xFFFF0000FF000000ull, //255.255.0.0
	0xFF000000FF000000ull, //255.0.0.0
	0x00000000FF000000ull, //0.0.0.0

	//0.0.0.0
	//20-24
	0xFFFFFFFF00000000ull, //255.255.255.255
	0xFFFFFF0000000000ull, //255.255.255.0
	0xFFFF000000000000ull, //255.255.0.0
	0xFF00000000000000ull, //255.0.0.0
	0x0000000000000000ull  //0.0.0.0
};

double dblmainmax(double a, double b) {return (a >= b ? a : b);}

int main(int argc, char * argv[]) {
		int m; //number of heavy hitters in output
		//scanf("%d%d%d", &counters, &threshold, &n);
		int counters = 100;
		int threshold = 1000;
		int n = 100000;
		double time, walltime;
		double epsil;
		HeavyHitter * ans;
		int i,j;
		int w, x, y, z;
		clock_t begint, endt;
		struct timeb begintb, endtb;
		unsigned long long ip, _ip;
		unsigned long long * data;
		double measurements[NRRUNS] = {0};
		double avg = 0, stdev = 0, upCI = 0, dnCI= 0, sos = 0;		
		double t_4 = 2.572; //Change if NRRUNS!=5
		FILE * fp = NULL;
		FILE * fsummary = NULL; //summary file	
		double insertProb = 1. / VMULT;	
		if (getenv("SUMMARYFILE") != NULL) fsummary = fopen(getenv("SUMMARYFILE"), "a");

		if (argc > 1) n = atoi(argv[1]);
		if (argc > 2) counters = atoi(argv[2]);
		if (argc > 3) threshold = atoi(argv[3]);
		if (argc > 4) fp = fopen(argv[4], "w");

		if(n/counters >= threshold) {
			printf("Unacceptable parameters: eps*n >= theshold\n");
			return 0;
		}

		

		data = (unsigned long long *) malloc(sizeof(unsigned long long) * n);

		for (i = 0; i < n; i++) {
			scanf("%d%d%d%d", &w, &x, &y, &z);
			ip = (unsigned long long)256*((unsigned long long)256*((unsigned long long)256*w + x) + y) + z;
			scanf("%d%d%d%d", &w, &x, &y, &z);
			_ip = (unsigned long long)256*((unsigned long long)256*((unsigned long long)256*w + x) + y) + z;
			data[i] = (ip << 32 | _ip);			
		}	
		for (j=0; j < NRRUNS; ++j){
			//printf("j = %d\n", j);
#ifdef RANDHHH
		init((double)1/(double)counters, insertProb);
#else
		init((double)1/(double)counters);
#endif
			begint = clock();
			ftime(&begintb);
			for (i = 0; i < n; i++) {
				update(data[i]);
			}
			endt = clock();
			ftime(&endtb);

			time = ((double)(endt-begint))/CLK_PER_SEC;
			walltime = ((double) (endtb.time-begintb.time))+((double)endtb.millitm-(double)begintb.millitm)/1000;
			measurements[j] = time;
			avg += time;
			if (j < NRRUNS - 1)
				deinit();
		}
		
		avg /= NRRUNS;
		
		for (j=0; j < NRRUNS; ++j){
			sos += pow(measurements[j] - avg, 2);
		}
		stdev = pow(sos / (NRRUNS - 1), 0.5);
		
		upCI = avg + stdev * t_4 / pow(NRRUNS, 0.5) ;
		dnCI = avg - stdev * t_4 / pow(NRRUNS, 0.5) ;

		free(data);

		printf("%d pairs took %fs (%fs-%fs) ", n, avg, dnCI, upCI);
		ans = output2(threshold, &m);

		printf("%d HHHs\n", m);

		deinit();

		if (fp != NULL) {
			//fprintf(fp, "%d\n", m);
			fprintf(fp, "%s %d %d %d %d %lf\n", argv[0], n, counters, threshold, m, time);
			for (i = 0; i < m; i++) {
				fprintf(fp, "%d %d.%d.%d.%d %d.%d.%d.%d %d %d\n",
				ans[i].mask,
				(int)((ans[i].item >> 56) & (LCLitem_t)255),
				(int)((ans[i].item >> 48) & (LCLitem_t)255),
				(int)((ans[i].item >> 40) & (LCLitem_t)255),
				(int)((ans[i].item >> 32) & (LCLitem_t)255),
				(int)((ans[i].item >> 24)& (LCLitem_t)255),
				(int)((ans[i].item >> 16) & (LCLitem_t)255),
				(int)((ans[i].item >> 8) & (LCLitem_t)255),
				(int)((ans[i].item >> 0) & (LCLitem_t)255),
#ifdef RANDHHH				
				ans[i].lower * ( (int) (NUM_COUNTERS / insertProb)), ans[i].upper * ( (int) (NUM_COUNTERS / insertProb)));
#else
				ans[i].lower, ans[i].upper);
#endif		
			}
			fclose(fp);
		}

		epsil = -1;
		for (i = 0; i < m; i++) {
			epsil = dblmainmax(epsil, ((double)(ans[i].upper-ans[i].lower))/n);
		}

		free(ans);

		if (fsummary != NULL) {
			fprintf(fsummary, "check=false algorithm=%-14s nitems=%-10d counters=%-5d threshold=%-9d outputsize=%-3d time=%lf walltime=%lf epsilon=%lf",
		                           argv[0],        n,           counters,     threshold,     m,              time, walltime,   epsil);
			fprintf(fsummary, "\n");
			fclose(fsummary);       
		}
	return 0;
}

