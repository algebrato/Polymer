#ifndef __MPMOVES__
#define __MPMOVES__
#include "PPolimerView.h"
#include "PTimeSectorView.h"
#include "RNumgen.h"
#include "DAverages.h"
#include "IVector.h"
#include "IKinds.h"
#include <cmath>
#include <fstream>
#include <iostream>

#define DIAGONAL 1
#define OFFDIAGONAL 0
#define ATTEMPTS 1000
/*
********* This class provides monte carlo moves to sample *********
********* an arbitrary probability distribution.          *********
********* It supports both zero and finite temperature 	  *********
********* calculations and contains identical particle	  *********
********* exchanges moves for bosons manybody systems.    ********* 
********* MPolymerMoves class is also provided with       *********
********* canonical worm functions.                       *********

*/


class MPolymerMoves
{
	public:
		MPolymerMoves(PPolimerView*,PTimeSectorView*, int=1);
		MPolymerMoves(PPolimerView*,PTimeSectorView*,RNumgen*,int=0);
		void setParameters(ifstream&);
		void polymerTranslation(unsigned);
		void polymerSingleMove(unsigned,unsigned, unsigned=0);
		void brownianBridge(unsigned,unsigned,unsigned);
		void IBrownianBridge(unsigned,unsigned,unsigned=0);
		virtual double polymerEvalDMatrix(unsigned, unsigned);
		int polymerTerminateMove(bool = false);
		double acceptedTranslations(int=0);
		double acceptedSingleMoves(int=0);
		double acceptedShadowMoves(int=0);
		double acceptedBB(int=0);
		double acceptedWormOpen(int=0);
		double acceptedWormClose(int=0);
		double acceptedWormSwap(int=0);
		void displayMovesAcceptation();
		double WORMsectorsRatio();
		void WORMInit(double);
		virtual void gotoZSector();
		void displayExchanges();
		void associateDAverages(DAverages *src){observables=src;src->setSSFTypes(typemove);}
		virtual void mcMove(unsigned) = 0;
		IVector* buildWORMPermutationCycle(int=0,int=1);
		IVector* buildPermutationCycle(int&,int=0,int =0, int=-1);
		
		int terminatePermutation(IVector*, int,int=0);
		void permutePolymers(IVector*,int,int=0);
		Point pbc(double,double,double);
		double pbc(double,int);
		Point pbc(Point&);
		Point pbc(double*);
		double abc(double*,double*);
		double abc(double,double,int);
		double abc(Point&,Point&);
		// WORM FUNCTIONS
		int WORMopen(int,int,int=0);
		int WORMclose();
		void WBB(int,int,int);  // Improved Brownian Bridge 
		int WORMterminate(bool=false);
		int WORMswap(int=0);
		int isDiagonal(){return ZoG;}
		double calculateSquareWindingNumber();
		double gaussianPSI(int,int,int);
		void gaussianPSIinit(double gc=8);
		void setJastrow(double b=2.81160,double m=5);
		
		
	private:
		double** dallmax;
		double** dsingmax;
		double** dshadmax;
		int *totTrasl,*totSM,*totSH,*totBB,*totPerm;
	
		double** temporaryCoords;
		PLink* ibbstart;
		int ibblength;
		
		void initializeCounters();
		bool metropolisAcceptance();
		int particle;
		void checkLinks();
		double perm_k_element(int,int,int,int,int*);
		double perm_get_C(int,int,int,int,int*);
		
		double sum_over_kelem(int, int,int,int,int*);
		double perm_k_elementWORM(int,int,int,int*,int=0);
		double sum_over_kelemWORM(int,int,int,int*,int=0);
		
		int faked_roulette(int, int,int,int,int*);
		int faked_WORMroulette(int,int,int,int*,int=0);
		double wclspring;
		bool can_accept_permute(int, int,int,int,int*);
		
		double swapkinetic(int,int,int);
		double** PSIs;
		double** PSILs;
		double gaussC;
		double jastrow_b, jastrow_m;
		bool jastrow;
		double pseudo_jastrow(double);
		double expo_jastrow(double r){return pow(jastrow_b/r,jastrow_m);}
		int fperm;
		double jashm;
		
	protected:
		
		bool shadow,optshadow;
		int *acpSM,*acpTrasl,*acpBB,*acpSH,particleNumber;
		int *acpPermute, **exchanges;
		int *bbMax, *excMax;
		int *bbTry, *excTry;
		int *wormM;
		bool* typemove;
		PTimeSectorView* timeSectorView; // hope's gonna work
		PPolimerView* polymerView;
		DAverages *observables;
		RNumgen* generator;
		IKinds* pkinds;
		int dimension;
		double box_dimensions[3];
		double pacc, kinmod;
		double varmod;
		int p1,p2;
		double kineticSpring(double*,double*,double);
		void undoBrownianBridges(IVector*,int,int);
		void releaseOldSegments(IVector*);
		void performExchanges(IVector*, int); // IF YOU FIRE A SEG.FAULT I WILL .... GRRRRRR !!!!
		void resetExchanges(IVector*, int); // working at first stroke? trascending murphy's law?????
		void performBrownianBridges(IVector*, int,int);
		void checkPBC(PLink&);
		void checkPBC(double*);
		void checkPBC(double&,int);
		int terminateIBB();
		int* pmoved;
		int npmoved;
		double shiftc;
		double windX, windY, windZ;
		double density;
		int rank;
		bool permuting;
		
		
		// WORM VARIABLES
		int ZoG;  // diagonal or off-diagonal?
		int nwormed;   // polymer opened
		int gnum, znum;
		int posworm;
		int* totWOpen;
		int* totWClose;
		int* totSwap;
		int *acpWOpen, *acpWClose, *acpSwap;
		bool worming;
		double wormacpconstant;
		void WORMunswap();
		bool swapping,swapping_reverse;
		IVector* wormswap;
		double wctemp[3], wotemp[3], scoord[3];
		double swmod;
		
		// repulsive factor function and variables
		double repulsive(double);
		double old_IMdistance;
		double new_IMdistance;
		double alfa_rep, beta_rep, gamma_rep, delta_rep;
		double rep_norm;
		
		double pi;
		// variational wave functions...
		double PSItry(int,int,int);
		double PSIJastrow(int,int,int);
		double PSIOptReal(int,int,int);
		void swapGaussian(int,int);
		
		double box_less;
		int shnumb;
		
		
		
		// kinetic variational estimators
		double KElocEstimator(int,int);
		double KElocJWF(int,int);
		double KElocSWF(int,int);
		double KElocORE(int,int);
		double KElocOSH(int,int);
		
		double KElocGWF(int,int);
		double KElocSH(int,int);
		int nshint;
		double shelocL, shelocR;
		double** frprimer;
		double*** pbcr;
		double** frprimel;
		double*** pbcl;
		double** abcr;
		double** abcl;
		
		
		
};

#endif
