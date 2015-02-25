#ifndef __RNGEN__
#define __RNGEN__
#include <fstream>

class RNumgen
{
	public:
		RNumgen(int *seeds, int set1, int set2, bool which);
		void flushGen();
		double runGen();
		double gaussGen(double, double);
	private:
		double step;
		int numb;
		int m[4], n[4],l[4];
		double current;
		int keys[4];
		std::ifstream *inf;
};


#endif
