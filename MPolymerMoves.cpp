#include "MPolymerMoves.h"
using namespace std;

// Constructor initializes random number generator.
MPolymerMoves::MPolymerMoves(PPolimerView* pview,PTimeSectorView* ptview, int rnk)
{
	polymerView=pview;
	timeSectorView=ptview;
	pkinds=pview->kinds;
	// read from file....!!
	int seeds[]={0,0,0,1};
	ifstream primes("PRIMES");
	
	int p1,p2;
	for(int i=0;i<rnk;i++)
	{
		primes>>p1;
		primes>>p2;
	}
	primes.close();
	generator = new RNumgen(seeds,p1,p2,true);
	particleNumber= polymerView->particleNumber();
	pbcr = new double**[particleNumber];
	pbcl = new double**[particleNumber];
	abcr = new double*[particleNumber];
	abcl = new double*[particleNumber];
	
	frprimer = new double*[particleNumber];
	frprimel = new double*[particleNumber];
	for(int i=0;i<particleNumber;i++)
	{
		pbcr[i]= new double*[particleNumber];
		pbcl[i]= new double*[particleNumber];
		abcr[i]= new double[particleNumber];
		abcl[i]= new double[particleNumber];
		
		frprimer[i]= new double[particleNumber];
		frprimel[i]= new double[particleNumber];
		for(int d=0;d<particleNumber;d++)
		{
			pbcr[i][d]=new double[polymerView->getDimension()];
			pbcl[i][d]=new double[polymerView->getDimension()];
		}
	}
	rank=rnk-1;
	initializeCounters();
	
	
}

MPolymerMoves::MPolymerMoves(PPolimerView* pview,PTimeSectorView* ptview,RNumgen* gen,int rnk)
{
	polymerView=pview;
	pkinds=pview->kinds;
	timeSectorView=ptview;
	// read from file....!!
	generator = gen;
	rank=rnk;
	initializeCounters();
	
}

// initialization of acceptance ratio counters.
void MPolymerMoves::initializeCounters()
{
	jastrow=false;
	permuting=false;
	shadow = polymerView->shadow;
	optshadow = polymerView->optshadow;
	
	shnumb=0;
	if(shadow)
		shnumb=1;
	varmod=1;
	gaussC=0;
	int knumb=pkinds->species_number;
	
	dsingmax=new double*[knumb];
	dshadmax=new double*[knumb];
	dallmax=new double*[knumb];
	for(int i=0;i<knumb;i++)
	{
		dsingmax[i]=new double[3];
		dshadmax[i]=new double[3];
		
		dallmax[i]=new double[3];
		for(int j=0;j<3;j++)
		{
			dsingmax[i][j]=0;
			dshadmax[i][j]=0;
			dallmax[i][j]=0;
		}
	}
	
	bbMax= new int[knumb];
	bbTry= new int[knumb];
	excMax= new int[knumb];
	excTry= new int[knumb];
	wormM = new int[knumb];
	
	acpTrasl=new int[knumb];
	acpSM=new int[knumb];
	acpSH=new int[knumb];
	
	acpBB=new int[knumb];
	acpWOpen=new int[knumb];
	acpWClose=new int[knumb];
	acpSwap=new int[knumb];
	totBB=new int[knumb];
	totTrasl=new int[knumb];
	totSM=new int[knumb];
	totSH=new int[knumb];
	acpPermute=new int[knumb];
	totPerm=new int[knumb];
	totWOpen=new int[knumb];
	totWClose=new int[knumb];
	totSwap=new int[knumb];
	typemove=new bool[knumb];
	exchanges = new int*[knumb];
	for(int i=0;i<knumb;i++){
		exchanges[i]=new int[particleNumber];
		for(int j=0;j<particleNumber;j++)
			exchanges[i][j]=0;
		
		acpTrasl[i]=0;
		totTrasl[i]=0;
		acpSM[i]=0;
		acpSH[i]=0;
		acpBB[i]=0;
		acpWClose[i]=0;
		acpWOpen[i]=0;
		acpSwap[i]=0;
		totBB[i]=0;
		totSM[i]=0;
		totSH[i]=0;
		acpPermute[i]=0;
		totPerm[i]=0;
		totWOpen[i]=0;
		totWClose[i]=0;
		totSwap[i]=0;
	}
	
	polymerView->singTot=totSM;
	polymerView->shTot=totSH;
	polymerView->acpSing=acpSM;
	polymerView->acpSh=acpSH;
	polymerView->tranTot=totTrasl;
	polymerView->acpTran=acpTrasl;
	polymerView->bbTot=totBB;
	polymerView->acpBB=acpBB;
	polymerView->wopenTot=totWOpen;
	polymerView->acpwopen=acpWOpen;
	polymerView->wcloseTot=totWClose;
	polymerView->acpwclose=acpWClose;
	polymerView->wswapTot=totSwap;
	polymerView->acpwswap=acpSwap;
	polymerView->acpPermute=acpPermute;
	polymerView->totPerm=totPerm;
	polymerView->exchanges=exchanges;
	polymerView->excTry = excTry;
		
	Point np = polymerView->getBoxDimensions();
	for(int i=0;i<3;i++)
		box_dimensions[i]=np.coords[i];
		
	double vol=1;
	int dim = polymerView->getDimension();
	for(int i=0;i<dim;i++)
		vol *= box_dimensions[i];
	
	box_less = box_dimensions[0];
	for(int i=0;i<dim;i++)
		if(box_less>box_dimensions[i])
			box_less=box_dimensions[i];
	
	density=particleNumber/vol;
	kinmod=1;
	ibbstart=NULL;
	ZoG=DIAGONAL;
	gnum=0;
	znum=0;
	windX=0;
	windY=0;
	windZ=0;
	worming=false;
	
	wormacpconstant=0;
	swapping=false;
	swapping_reverse=false;
	
	ifstream repf("input/repulsive_factor");
	repf>>alfa_rep;
	repf>>beta_rep;
	repf>>gamma_rep;
	repf>>delta_rep;
	repf.close();
	rep_norm = 1 + alfa_rep + gamma_rep;
	pi=acos(-1.);
	old_IMdistance=0;
	new_IMdistance=0;
	
	
	// jastrow L/2 value
	jashm = expo_jastrow(box_less/2);
	

	
}
// functions used to display acceptances values
double MPolymerMoves::acceptedTranslations(int i)
{
	return (double)acpTrasl[i]/totTrasl[i];
}

double MPolymerMoves::acceptedSingleMoves(int i)
{
	return (double)acpSM[i]/totSM[i];
}

double MPolymerMoves::acceptedShadowMoves(int i)
{
	return (double)acpSH[i]/totSH[i];
}


double MPolymerMoves::acceptedBB(int i)
{
	return (double)acpBB[i]/totBB[i];
}

double MPolymerMoves::acceptedWormOpen(int i)
{
	return (double)(acpWOpen[i])/totWOpen[i];
}

double MPolymerMoves::acceptedWormClose(int i)
{
	return (double)(acpWClose[i])/totWClose[i];
}

double MPolymerMoves::WORMsectorsRatio()
{
	return ((double)znum)/gnum;
}

double MPolymerMoves::acceptedWormSwap(int i)
{
	return ((double)acpSwap[i])/totSwap[i];
}
// This function reads moves parameters from "deltapar" file.
void MPolymerMoves::setParameters(ifstream& in)
{
	in>>dimension;
	timeSectorView->setDim(dimension);
	int snumb=pkinds->species_number;
	polymerView->typemove=typemove;
	for(int i=0;i<snumb;i++)
		typemove[i]=false;
	int mcount=0;
	bool first=true;
	char* tag = new char[10];
	in>>tag;
	while(!in.eof())
	{	
		
		int nkind=-1;
		
		for(int i=0;i<snumb;i++)
		{
			if(strcmp(tag,pkinds->species[i])==0)
			{
				nkind=i;
				mcount++;
				break;
			}
		}
		if(nkind==-1 && first && rank==0){
			first = false;
			cout<<"WARNING: some specie ain't gonna be moved."<<endl;
			
			continue;
		}
		else if(nkind==-1)
			continue;
		else
		{
			typemove[nkind]=true;
			
			for(int i=0;i<dimension;i++)
				in>>dallmax[nkind][i];
			for(int i=0;i<dimension;i++)
				in>>dsingmax[nkind][i];
			for(int i=0;i<dimension;i++)
				in>>dshadmax[nkind][i];
		
			in>>nshint;
			in>>bbMax[nkind];
			in>>bbTry[nkind];
			in>>excMax[nkind];
			in>>excTry[nkind];
			in>>wormM[nkind];
		}
		in>>tag;
		
	}
	delete tag;
	shiftc=0;
			
	if(rank==0)
		cout<<"I'm gonna move "<<mcount<<" species !"<<endl;
}

// Performs a single link move on a polymer. 'which' variable determines whether moving
// left or right bead.
void MPolymerMoves::polymerSingleMove(unsigned npart, unsigned tslice, unsigned which)
{
	
	double **move = dsingmax;
		
	particle=npart;
	npmoved=1;
	pmoved=new int;
	*pmoved=particle;
	
	p1=tslice;
	p2=tslice;
	PLink mylink = *((*polymerView)[npart][tslice]);
	
	int type=mylink.type;
	double *Imove1=NULL;
	double *Imove2=NULL;
	double *Icorr1,*Icorr2;
	double k1=1;
	double k2=1;
	
	int nlinks = (*polymerView)[npart].linkCount();
	bool isopt=false;
	if(shadow && (tslice == 0 || tslice==nlinks-1))
		move=dshadmax;
	
	if((tslice==0 || tslice==nlinks-1) && optshadow)
		k1=pkinds->itabs[type][type]->evaluateOPTShadow(sqrt(abc(mylink.firstBeadCoordinates,mylink.secondBeadCoordinates)),1);
	else
		k1=kineticSpring(mylink.firstBeadCoordinates,mylink.secondBeadCoordinates,mylink.kineticPrefactor);
	
	if(!shadow && ((tslice==0 && which==LEFT) || (tslice==nlinks-1 && which==RIGHT)) && polymerView->isPIGS())
		varmod=1./PSItry(npart,tslice,which);
	else if(shadow && ((tslice==0 && which==RIGHT) || (tslice==1 &&	which==LEFT) || (tslice==nlinks-1 && which==LEFT)|| (tslice==nlinks-2 && which==RIGHT)) && polymerView->isPIGS())
		varmod=1./PSItry(npart,tslice,which);
	
	double ad_kin;
	
	if(which==LEFT)
	{
		Imove1=mylink.firstBeadCoordinates;
		if(mylink.linked){
			ad_kin=mylink.pprev->kineticPrefactor;
			Imove2=mylink.pprev->secondBeadCoordinates;
			if(tslice==1 && optshadow){
				k2=pkinds->itabs[type][type]->evaluateOPTShadow(sqrt(abc(mylink.firstBeadCoordinates,mylink.pprev->firstBeadCoordinates)),1);
				isopt=true;
			}
			else
				k2=kineticSpring(mylink.firstBeadCoordinates,mylink.pprev->firstBeadCoordinates,ad_kin);
			Icorr1=mylink.firstBeadCoordinates;
			Icorr2=mylink.pprev->firstBeadCoordinates;
		}
	}
	else
	{
		Imove1=mylink.secondBeadCoordinates;
		if(mylink.pnext!=NULL && mylink.pnext->linked){
			ad_kin=mylink.pnext->kineticPrefactor;
			Imove2=mylink.pnext->firstBeadCoordinates;
			if(tslice==nlinks-2 && optshadow){
				k2=pkinds->itabs[type][type]->evaluateOPTShadow(sqrt(abc(mylink.secondBeadCoordinates,mylink.pnext->secondBeadCoordinates)),1);
				isopt=true;
			}
			else
				k2=kineticSpring(mylink.secondBeadCoordinates,mylink.pnext->secondBeadCoordinates,ad_kin);
			Icorr1=mylink.secondBeadCoordinates;
			Icorr2=mylink.pnext->secondBeadCoordinates;
		}
	
	}
	
	
	kinmod*=k1*k2;  

	double ncoords[dimension];
	
	for(int i=0;i<dimension;i++){
		ncoords[i]=(generator->runGen()-0.5)*move[type][i];
		Imove1[i]+=ncoords[i];
		Imove1[i]=pbc(Imove1[i],i);
	}
	if(Imove2!=NULL)
		for(int i=0;i<dimension;i++)
			Imove2[i]=Imove1[i];

	
	if(optshadow && (tslice==0 || tslice == nlinks-1))
		k1=pkinds->itabs[type][type]->evaluateOPTShadow(sqrt(abc(mylink.firstBeadCoordinates,mylink.secondBeadCoordinates)),1);
	else
		k1=kineticSpring(mylink.firstBeadCoordinates,mylink.secondBeadCoordinates,mylink.kineticPrefactor);
	if(Imove2!= NULL)
		if(isopt)
			k2=pkinds->itabs[type][type]->evaluateOPTShadow(sqrt(abc(Icorr1,Icorr2)),1);
		else
			k2=kineticSpring(Icorr1,Icorr2,ad_kin);
		
		
	
	(*polymerView)[npart].cutPolymer(tslice,tslice);
	(*polymerView)[npart].insertLink(mylink,tslice);

	
	kinmod /=k1*k2; 
	if(!shadow && ((tslice==0 && which==LEFT) || (tslice==nlinks-1 && which==RIGHT)) && polymerView->isPIGS())
		varmod*=PSItry(npart,tslice,which);
	else if(shadow && ((tslice==0 && which==RIGHT) || (tslice==1 &&	which==LEFT) || (tslice==nlinks-1 && which==LEFT)|| (tslice==nlinks-2 && which==RIGHT)) && polymerView->isPIGS())
		varmod*=PSItry(npart,tslice,which);

	if(shadow && (tslice==0 || tslice == nlinks-1))
		totSH[type]++;
	else
		totSM[type]++;
	

}

// Performs a rigid translation of a single polymer.   (note: need to fix spelling...forgot an n Oops)
// If, while scanning polymer, this function finds out a worm, new translation coordinates will be generated
// so that two polymer unlinked parts will translate indipendently
void MPolymerMoves::polymerTranslation(unsigned nparticle)
{
  
	particle=nparticle;
	npmoved=1;
	pmoved=new int;
	*pmoved=particle;
	int type=(*polymerView)[particle][0]->type;
	double ncoords[dimension];
	for(int i=0;i<dimension;i++)
		ncoords[i]=(generator->runGen()-0.5)*dallmax[type][i];
	totTrasl[type]++;
	int nlinks=(*polymerView)[particle].linkCount();
	p1=0;
	p2=nlinks-1;
	PLink oldlinks[nlinks];
	
	if(polymerView->isPIGS())
		if(!shadow)
			varmod = 1./(PSItry(particle,0,LEFT)*PSItry(particle,nlinks-1,RIGHT));
		else
			varmod = 1./(PSItry(particle,0,RIGHT)*PSItry(particle,nlinks-1,LEFT));
		
	
	
	for(int i=0;i<nlinks;i++)
	{
		oldlinks[i]=*((*polymerView)[particle][i]);
		if(!oldlinks[i].linked)
			for(int j=0;j<dimension;j++)
				ncoords[j]=(generator->runGen()-0.5)*dallmax[type][j];	
		
		for(int j=0;j<dimension;j++)
		{
			oldlinks[i].firstBeadCoordinates[j]+=ncoords[j];
			oldlinks[i].secondBeadCoordinates[j]+=ncoords[j];
			checkPBC(oldlinks[i].firstBeadCoordinates[j],j);
			checkPBC(oldlinks[i].secondBeadCoordinates[j],j);
			
		}
		//checkPBC(oldlinks[i].firstBeadCoordinates);
		//checkPBC(oldlinks[i].secondBeadCoordinates);

	}
	
	
	(*polymerView)[particle].cutPolymer(0,nlinks-1);
	
	for(int i=0;i<nlinks;i++)
		(*polymerView)[particle].insertLink(oldlinks[i],i);
	
	if(polymerView->isPIGS())
		if(!shadow)
			varmod *= (PSItry(particle,0,LEFT)*PSItry(particle,nlinks-1,RIGHT));
		else
			varmod *= (PSItry(particle,0,RIGHT)*PSItry(particle,nlinks-1,LEFT));
		
	
}

// This function sorts a random number and determines the success of a move
// with probability pacc.
bool MPolymerMoves::metropolisAcceptance()
{
	pacc=pacc/kinmod;
	pacc=pacc*varmod;
	varmod=1;
	kinmod=1;
	
	if(generator->runGen()<=pacc){
		old_IMdistance=new_IMdistance;
		return true;
	}
	else{
		new_IMdistance=old_IMdistance;
		return false;
	}
}
// This is a function which returns the correlation value between 
// pos1 and pos2 timeslices. You HAVE to overload it in a derived class
// to adapt it on your needs.
double MPolymerMoves::polymerEvalDMatrix(unsigned pos1, unsigned pos2)
{
	
	pacc=0;
	return 0;
		
}
// This method terminates a polymer move (single, translation or brownian
// bridge) handling acceptance and rejection with probability pacc
// determined by polymerEvalDMatrix method.
// The scheme of a move is therefore: <move>() -> polymerEvalDMatrix
// -> polymerTerminateMove();. This implementation is essential
// if you have to set a particular density matrix  for every system
// in inspection because you can derive this class and overload 
// polymerEvalDMatrix function.
int MPolymerMoves::polymerTerminateMove(bool opt)
{
	npmoved=0;
	delete pmoved;
  
        if(ibbstart!=NULL)
		return terminateIBB();
		
	int nlnks = (*polymerView)[particle].linkCount();
	int np1=p1;
	int np2=p2;
	if(!opt)
	{
		if(np1>0)
			np1--;
		if(np2<nlnks-1)
			np2++;
	}
	if(metropolisAcceptance())
	{
		
		(*polymerView)[particle].freeAnchor();
		
	
		for(int i=np1;i<=np2;i++){
			
			(*timeSectorView)[i].setChanged(true);  // I accepted move: I need to change time sector correlations
			(*timeSectorView)[i].refreshCorrelationValue();
			
			
		}
		
		return 1;
	}
	else{
		
		(*polymerView)[particle].deletePolymer(p1,p2);
		
		(*polymerView)[particle].restorePolymer();
		for(int i=np1;i<=np2;i++){
			//(*timeSectorView)[i].setChanged(true); 
			(*timeSectorView)[i].undoChanges();
						
		}
		
		

		return 0;
	}
	
}
// This is a function which upgrades links between two polymers.
// It is automatically called in MC moves which require it.
// This function was used for debug, now every move handles first and second bead coordinates and this function
// has come to be even malevolous. That's why I placed a return on first line (just in case in a drunken status I'm gonna
// call it...)
void MPolymerMoves::checkLinks()
{	
	return;
	for(int i=p1;i<=p2;i++)
	{
		PLink* mylink = (*polymerView)[particle][i];
		if(mylink->linked)
			for(int j=0;j<dimension;j++)
				mylink->pprev->secondBeadCoordinates[j]=mylink->firstBeadCoordinates[j];
	}
}


// This function crafts a brownian bridge between two slices of a polymer.
// It is meant to modify inner beads of a link adjacent selection.
// Explain: >-<>-<>-<>-<>-<>-<>-< 
//          1 22 33 44 55 66 77 8 
// if I set start=1, end = 7, B.B will be fashioned from bead 2 to bead 7.
// P.S.: another implementation has been done. See WBB.

void MPolymerMoves::brownianBridge(unsigned start, unsigned end, unsigned npart)
{
	particle=npart;
	int firstCut=start-1;
	int secondCut=end+1;
	int b_length = end-start+1;

	p1=start;
	p2=end;
	
	PLink oldlinks[b_length];
	
	int k=0;
	PLink* clink = ((*polymerView))[particle][p1]; // hardcopying old links. 
	totBB[clink->type]++;
	
	for(int i=p1;i<=p2;i++){
		oldlinks[k++]=*clink;
		clink=clink->pnext;
	}
	
	double* endCoords=oldlinks[b_length-1].secondBeadCoordinates;
	(*polymerView)[particle].cutPolymer(start,end);
	double kp = oldlinks[b_length-1].kineticPrefactor;
	
	// This B.B. implementation is not difficult to get.
	// Just keep in mind: link simmetry enables me to scan polymer interval
	// two by two (but if interval length's even...)
	// Also remind that: looping two by two means that I modify 
	// two links in a single iteration. This allows to "combine"
	// togather dimension loop  of link k and k+1.  
	k=1;
	int cn=2;
	for(int i=start+1;i<=end-1;i+=cn)  
	{
		double rappL = ((double)(end - i + 1))/(end-i+2);
		double rappR = ((double)(end - i))/(end-i+1);
		
		
		for(int j=0;j<dimension;j++)
		{
			double xprev=oldlinks[k-1].firstBeadCoordinates[j];
			double xvar= pbc(endCoords[j]-xprev,j);
			double xavg= xprev + xvar/(end-i+2);
			double newcoord =pbc(generator->gaussGen(xavg,rappL/kp),j);
			
			oldlinks[k-1].secondBeadCoordinates[j]= newcoord;
			oldlinks[k].firstBeadCoordinates[j]=newcoord;
			// chain to next
			if (cn==1)
				continue;
				
			xprev=newcoord;
			xvar= pbc(endCoords[j]-xprev,j);
			xavg= xprev + xvar/(end-i+1);
			newcoord =pbc(generator->gaussGen(xavg,rappR/kp),j);
			oldlinks[k].secondBeadCoordinates[j]=newcoord;
			oldlinks[k+1].firstBeadCoordinates[j]=newcoord;
			
		}
		if(i==end-2)
			cn=1;
		k+=cn;
			
	}
	k=0;
	for(int i=start;i<=end;i++)
		(*polymerView)[particle].insertLink(oldlinks[k++],i);
	
	

	
	
}
// from start, try to rebuild hmany beads. IBrownianBridge is the function you call
// to perform a Brownian Bridge. PIMC, PIGS, WORM: you don't have to worry about it.
// IBB's gonna pick up the right function for you :)
// (note to self: for now in PIMC case it's only gonna pick up a segmentation fault)
void MPolymerMoves::IBrownianBridge(unsigned start, unsigned npart, unsigned hm)
{
	int type=(*polymerView)[npart][0]->type;
	int hmany=hm;
	if(hm==0)
		hmany=bbMax[type];
		
	if(!permuting)
	{
		npmoved=1;
		pmoved=new int;
		*pmoved=npart;
	}
	
	if(start>=0&&((*polymerView)[npart][start+hmany-1]!=NULL)) // am I inside a polymer?
	{
		//brownianBridge(start,start+hmany-1,npart);
		WBB(start,start+hmany-1,npart);
		
		return;
	}
	
	
	cout<<"IBB error"<<" "<<npart<<" "<<start<<" "<<start+hmany-1<<" "<<hm<<endl;
	// not inside a polymer... then... go as long as you can !!!!!
	// NOTA: DA QUI IN POI E' TUTTO PIMC. NON ANCORA DEBUGGATO. E SICURAMENTE DA MODIFICARE (ottimizzazione...).
	int knt=0;
	PLink* clink=(*polymerView)[npart][start];
	ibbstart=clink;
	particle=npart;
	p1=start;
	while(clink->pnext!=NULL){
		knt++;
		if(knt>=hmany)
			break;
		clink=clink->pnext;
	}
	p2=start+knt;
	double *endpoint = clink->firstBeadCoordinates;
	clink=(*polymerView)[npart][start];
	temporaryCoords=new double*[knt];
	for(int i=0;i<knt;i++)
	{
		temporaryCoords[i]=new double[dimension];
		for(int j=0;j<dimension;j++)
			temporaryCoords[i][j]=clink->firstBeadCoordinates[j];
		
		clink=clink->pnext;
	}
	clink=(*polymerView)[npart][start];
	double *rbefore= clink->pprev->firstBeadCoordinates;
	totBB[clink->type]++;
		
	for(int i=0;i<knt;i++)
	{
		
		double kp= clink->kineticPrefactor;
		double xvar = (knt-i)/((knt-i+1)*kp);
		double xavg[dimension];
		//double *rbefore= clink->pprev->firstBeadCoordinates;
		for(int j=0;j<dimension;j++){
			xavg[j]= rbefore[j]+(pbc(endpoint[j]-rbefore[j],j))/(knt-i);
			clink->firstBeadCoordinates[j]=generator->gaussGen(xavg[j],xvar);
			checkPBC(clink->firstBeadCoordinates[j],j);
		}
		//checkPBC(clink->firstBeadCoordinates);
		rbefore= clink->firstBeadCoordinates;
		clink=clink->pnext;
	}
	ibblength=knt;
	
	
	
}
// restores 'interpolymeran' brownian bridge   ANCHE QUESTA E' DA SISTEMARE NEL PIMC.
// Con le permutazioni il BB e' completamente diverso: si estende su piu' polimeri, per cui
// la gestione e' necessariamente diversa (ma fortunatamente le correlazioni sono identiche...
// una volta capite le particelle coinvolte...)
int MPolymerMoves::terminateIBB()
{
	PLink* clink = ibbstart;
	if(metropolisAcceptance())
	{
		
		int nlinks=(*polymerView)[particle].linkCount();
		for(int i=0;i<nlinks;i++){
			(*timeSectorView)[i].setChanged(true);
			(*timeSectorView)[i].refreshCorrelationValue();
		}
		for(int i=0;i<ibblength;i++)
			delete temporaryCoords[i];
		
		
		delete temporaryCoords;
		ibbstart=NULL;
		return 1;
		
	}
	else{
		for(int i=0;i<ibblength;i++)
		{
			for(int j=0;j<dimension;j++)
				clink->firstBeadCoordinates[j]=temporaryCoords[i][j];
			
			clink=clink->pnext;
			delete temporaryCoords[i];
		}
		delete temporaryCoords;
		ibbstart=NULL;
		return 0;
		
	}
	
}

// In a permutation move you determine a permutation cycle with kinetic tests.
// This function returns a kinetic correlation between two point coordinates.
double MPolymerMoves::kineticSpring(double *x1, double* x2, double pref)
{
	double xs = abc(x1,x2);	
	return exp(-xs*pref);
}

// This method builds a permutation cycle and stores it in a IVector object.
// An n-particle closed exchanged will be composed of n+1 elements
// with last element equals to first. Only closed cycles are suitable
// for a permutation move. Only opened cycles will do for WORM swap.
IVector* MPolymerMoves::buildPermutationCycle(int &k, int type,int s, int pa)
{
	// don't use type.....for now
	// port to distinguishable particles case by changing 
	// random particle check/selection with a selective function
	// to exclude different real atoms.
	
	
	// cycle polymers if k+s > polymerlength or k < 0 !!
	if(s==0)
		s=excMax[type];
		
	int l=k+s-(*polymerView)[0].linkCount();
	if(l>=0)
	{
		k=k-l-1;
		polymerView->clockWork(-l);
		
	}
	if(k<0)
	{
		k=0;
		polymerView->clockWork(-k);
	}
	
	int partnumb=polymerView->particleNumber();
	int prevID;
	if(pa==-1)
		prevID= generator->runGen()*(partnumb-1);
	else
		prevID=pa;
		
	totPerm[(*polymerView)[prevID][k]->type]++;
	
	int firstID=prevID;
	IVector* cycle = new IVector(prevID);
	int bitmask[partnumb];
	int** neighCenter = timeSectorView->returnVerlet()->returnNeighCenter();
	int* neighNumber = timeSectorView->returnVerlet()->returnNeighNumber();
	
	for(int i=0;i<partnumb;i++)
		bitmask[i]=false;
	for(int i=0;i<neighNumber[firstID];i++)
		bitmask[neighCenter[firstID][i]]=true;
	bitmask[firstID]=true;
	
	/*if(ZoG==OFFDIAGONAL && (k==posworm || k == posworm-1)){
		bitmask[nwormed]=false;
		if(firstID==nwormed)
			return cycle;
	}*/
		
	//bitmask[firstID]=false; self-permutations already excluded in K
	bool loop=true;
	fperm=prevID;
	loop=can_accept_permute(prevID,k,s,partnumb,bitmask);
	int knt=1;
	while(loop)
	{
		int cID=faked_roulette(prevID,k,s,partnumb,bitmask);
		
		
			
		//if(loop){
			cycle->add(cID);
			bitmask[cID]=false;
			/*if(cID==firstID)
			{
				cycle->closed=true;
				loop=false;
				acpPermute[(*polymerView)[prevID][k]->type]++;
			}*/
			prevID=cID;
			knt++;
		//}
		loop = can_accept_permute(cID,k,s,partnumb,bitmask);
		if(cID==firstID)
		{
			cycle->closed=true;
			loop=false;
			acpPermute[(*polymerView)[prevID][k]->type]++;
		}
	}
		
	return cycle;	
	
}

IVector* MPolymerMoves::buildWORMPermutationCycle(int wise,int s)
{
	// don't use type.....for now
	// port to distinguishable particles case by changing 
	// random particle check/selection with a selective function
	// to exclude different real atoms.
	
	
	// cycle polymers if k+s > polymerlength or k < 0 !!
	int k= posworm-1;	
/*	int l=k+s-(*polymerView)[0].linkCount();
	if(l>=0)
	{
		k=k-l-1;
		polymerView->clockWork(-l);
		
	}
	if(k<0)
	{
		k=0;
		polymerView->clockWork(-k);
	}
	*/
	
	int partnumb=polymerView->particleNumber();
	int prevID;
	prevID=nwormed;
		
	int firstID=prevID;
	IVector* cycle = new IVector(prevID);
	int bitmask[partnumb];
	int** neighCenter = timeSectorView->returnVerlet()->returnNeighCenter();
	int* neighNumber = timeSectorView->returnVerlet()->returnNeighNumber();
	
	for(int i=0;i<partnumb;i++)
			bitmask[i]=true;
	for(int i=0;i<neighNumber[firstID];i++)
		bitmask[neighCenter[firstID][i]]=true;
	bitmask[firstID]=true;
		
	//bitmask[firstID]=false; self-permutations already excluded in K
	bool loop=true;
	
	int knt=1;
	while(loop)
	{
		int cID=faked_WORMroulette(prevID,s,partnumb,bitmask,wise);
		
			
		if(loop){
			cycle->add(cID);
			bitmask[cID]=false;
			if(cID==firstID)
			{
				cycle->closed=false;
				loop=false;
				break;
			}
			prevID=cID;
			knt++;
			if(knt==2){
				loop=false;
				cycle->closed=true;
				cycle->add(firstID); 
			// we force its closure because I'm lazy and wanna use
			// exchange functions even for worm swap :P
			}
		}
	}
	
	return cycle;	
	
}

// Theese are functions belonging to permutation determination algorytm.
// See extended documentation for description
double MPolymerMoves::perm_k_element(int p1, int p2,int k, int s, int* bmask)
{
	if(p1==p2)
		return 0;
	
	if(!bmask[p2])
		return 0;
	double kpre= (*polymerView)[p1][k]->kineticPrefactor;
	return kineticSpring((*polymerView)[p1][k]->firstBeadCoordinates,(*polymerView)[p2][k+s]->secondBeadCoordinates,kpre/(s+1));  // or first ?
	
	
}

double MPolymerMoves::perm_get_C(int pid,int k,int s, int ptot, int* bmask)
{
	double sum=0;
	double kpre= (*polymerView)[pid][k]->kineticPrefactor;
	
	for(int i=0;i<ptot;i++)
	{
		if(i==pid)
		{
			sum+=kineticSpring((*polymerView)[pid][k]->firstBeadCoordinates,(*polymerView)[i][k+s]->secondBeadCoordinates,kpre/(s+1)); // or first ?
			continue;
		}
		else if(i!=fperm)
			if(!bmask[i])
				sum+=kineticSpring((*polymerView)[pid][k]->firstBeadCoordinates,(*polymerView)[i][k+s]->secondBeadCoordinates,kpre/(s+1)); // or first?
		
	}
	
	return sum;

}

double MPolymerMoves::perm_k_elementWORM(int p1, int p2, int s, int* bmask, int wise)
{
	//if(p1==p2)
	//	return 0;
	
	if(!bmask[p2])
		return 0;
	double kpre= (*polymerView)[p1][posworm-1]->kineticPrefactor;
	/*if(wise==1)
		return	kineticSpring((*polymerView)[p1][posworm-1]->secondBeadCoordinates,(*polymerView)[p2][posworm+s-1]->firstBeadCoordinates,kpre/s);
	else
		return	kineticSpring((*polymerView)[p2][posworm-1]->secondBeadCoordinates,(*polymerView)[p1][posworm+s-1]->firstBeadCoordinates,kpre/s);
	*/
	if(wise==1)
		return kineticSpring((*polymerView)[p1][posworm]->firstBeadCoordinates,(*polymerView)[p2][posworm-s]->firstBeadCoordinates,kpre/s);
	else
		return kineticSpring((*polymerView)[p1][posworm-1]->secondBeadCoordinates,(*polymerView)[p2][posworm+s]->firstBeadCoordinates,kpre/s);
	
	
}

double MPolymerMoves::sum_over_kelem(int pid,int k,int s, int ptot, int* bmask)
{
	double sum=0;
	for(int i=0;i<ptot;i++)
		sum += perm_k_element(pid,i,k,s,bmask);

	return sum;
	
}

double MPolymerMoves::sum_over_kelemWORM(int pid,int s, int ptot, int* bmask,int wise)
{
	double sum=0;
	for(int i=0;i<ptot;i++)
		sum += perm_k_elementWORM(pid,i,s,bmask,wise);

	return sum;
	
}

int MPolymerMoves::faked_roulette(int pid, int k, int s, int ptot, int* bitmask)
{
	// partitioning [0;1] in probabilities.
	double prob_table[ptot];
	double Q=sum_over_kelem(pid,k,s,ptot,bitmask);
	double p= perm_k_element(pid,0,k,s,bitmask)/Q;
	prob_table[0]=p;
	prob_table[ptot-1]=1;
	for(int i=1;i<ptot-1;i++)
	{
		p= perm_k_element(pid,i,k,s,bitmask)/Q;
		prob_table[i]=prob_table[i-1]+p;
	}
	
	// now shuffle time... some say luck is blind... 
	// ... so we are gonna drive it ...
	
	double my_dice = generator->runGen();
	int winner=ptot-1;
	for(int i=0;i<ptot;i++){
		if(i==0 || prob_table[i]!=prob_table[i-1])
			winner=i;
			
		if(my_dice<prob_table[i])
			break;
		
	}
	
	return winner;
	
}

int MPolymerMoves::faked_WORMroulette(int pid, int s, int ptot, int* bitmask, int wise)
{
	// partitioning [0;1] in probabilities.
	double prob_table[ptot];
	double Q=sum_over_kelemWORM(pid,s,ptot,bitmask,wise);
	double p= perm_k_elementWORM(pid,0,s,bitmask,wise)/Q;
	
	prob_table[0]=p;
	prob_table[ptot-1]=1;
	for(int i=1;i<ptot-1;i++)
	{
		p= perm_k_elementWORM(pid,i,s,bitmask,wise)/Q;
		prob_table[i]=prob_table[i-1]+p;
	}
	
	// now shuffle time... some say luck is blind... 
	// ... so we are gonna drive it ...
	
	double my_dice = generator->runGen();
	int winner=ptot-1;
	for(int i=0;i<ptot;i++){
		if(i==0 || prob_table[i]!=prob_table[i-1])
			winner=i;
			
		if(my_dice<prob_table[i])
			break;
		
	}
	return winner;
	
}


bool MPolymerMoves::can_accept_permute(int pid, int k, int s, int ptot, int* bitmask)
{
	double Q=sum_over_kelem(pid,k,s,ptot,bitmask);
	//double kpre= (*polymerView)[pid][0]->kineticPrefactor;
		
	//double C=Q/(Q+kineticSpring((*polymerView)[pid][k]->firstBeadCoordinates,(*polymerView)[pid][k+s]->firstBeadCoordinates,kpre/s));
	double Cden = perm_get_C(pid,k,s,ptot,bitmask);
	double C = Q/(Q+Cden);
	if(generator->runGen()<C)
		return true;  // il polimero ci sta
	else
		return false; // il polimero gli da il due di picche
	
}

// Internal function of the permute method. Executes exchanges determined by
// buildPermutationCycle method.
void MPolymerMoves::performExchanges(IVector* permute, int k)
{
	// note that switchPolimerSlice performs a double exchange
	// therefore first element is automatically used to close loop
	// => no need to count it twice.
	int nexc = permute->elementsNumber(); 
	npmoved=nexc-1;
	pmoved=new int[npmoved];
	for(int i=0;i<npmoved;i++)
		pmoved[i]=permute->at(i);
	for(int i=0;i<nexc-2;i++)
	{
		(*timeSectorView)[k].switchPolimerSlice(permute->at(i),permute->at(i+1));
		timeSectorView->exchangeParticlesCorrelation(permute->at(i),permute->at(i+1),k);
		// swap gaussian WF !!!
		swapGaussian(permute->at(i),permute->at(i+1));
		
		if(permuting && ZoG==OFFDIAGONAL && k<posworm){
			if(permute->at(i+1)==nwormed)
				nwormed=permute->at(i);
			else if(permute->at(i)==nwormed)
				nwormed=permute->at(i+1);
		}
	}
}
// Undoes performed exchanged in rejection permutation case. Internal function.
// You don't have to invoke it unless you'd like to fire segmentation fault errors ;)
void MPolymerMoves::resetExchanges(IVector* permute, int k)
{
	int nexc= permute->elementsNumber();
	for(int i=nexc-3;i>=0;i--)
	{
		(*timeSectorView)[k].switchPolimerSlice(permute->at(i),permute->at(i+1));
		timeSectorView->exchangeParticlesCorrelation(permute->at(i),permute->at(i+1),k);
		swapGaussian(permute->at(i),permute->at(i+1));
	
	
		if(permuting && ZoG==OFFDIAGONAL && k<posworm){
			if(permute->at(i+1)==nwormed)
				nwormed=permute->at(i);
			else if(permute->at(i)==nwormed)
				nwormed=permute->at(i+1);
		}
	}
}

// This function builds brownian bridges of exchanged polymers.
void MPolymerMoves::performBrownianBridges(IVector* permute, int k, int s)
{
	int nint= permute->elementsNumber();
	pmoved = new int[nint-1];
	npmoved=nint-1;
	// remember: IVector is closed, last and first elements are the same !!!
	for(int i=0;i<nint-1;i++){
		pmoved[i]=permute->at(i);
		IBrownianBridge(k,permute->at(i),s);
		//WBB(k+1,permute->at(i),k+s-1);
		totBB[(*polymerView)[permute->at(i)][0]->type]--;  // what's that?? b.b. method increments totBB...but this is a permute move!
	}
}

// This function undoes the changes apported by the previous function.
void MPolymerMoves::undoBrownianBridges(IVector* permute,int k, int s)
{
	int nint= permute->elementsNumber();
	for(int i=0;i<nint-1;i++){
		//(*polymerView)[permute->at(i)].deletePolymer(k+1,k+s-1);
		(*polymerView)[permute->at(i)].deletePolymer(p1,p2);
		(*polymerView)[permute->at(i)].restorePolymer();
	}	
	
		
}
// In case of permute acceptance this internal method frees old polymers memory.
void MPolymerMoves::releaseOldSegments(IVector* permute)
{
	int nint= permute->elementsNumber();
	for(int i=0;i<nint-1;i++)
		(*polymerView)[permute->at(i)].freeAnchor();
		
		
}
// Call this function to terminate a permutation in MC step definition.
// This method will invoke needed functions to accept or reject changes.
int MPolymerMoves::terminatePermutation(IVector* permute, int k, int sm)
{
	
	if(!permute->closed)
		return 0;
		
	permuting = false;
	int s=sm;
	int type = (*polymerView)[permute->at(0)][k]->type;
	if(s==0)
		s = excMax[type];
	
	npmoved=0;
	delete pmoved;
	double val=generator->runGen();
	int nl=permute->elementsNumber();
	if(val<pacc){
		releaseOldSegments(permute);

		//for(int j=k;j<=k+s;j++)	{
		for(int j=p1;j<=p2;j++)	{
			(*timeSectorView)[j].setChanged(true);	
			(*timeSectorView)[j].refreshCorrelationValue();
		}	
		exchanges[type][permute->elementsNumber()-1]++;
		return 1;
		
	}	
	else{
		
		//for(int i=k;i<=k+s;i++)
		for(int i=p1;i<=p2;i++)
			(*timeSectorView)[i].undoChanges();
		
		undoBrownianBridges(permute,k,s);
		resetExchanges(permute,k);

		for(int i=0;i<nl-1;i++)
			(*polymerView)[permute->at(i)].updateEnd();

		return 0;
		
		
		
	}

	
}
// Call this method to begin polymers permutation. 
void MPolymerMoves::permutePolymers(IVector* permute, int k, int s)
{
	
	if(!permute->closed)
		return;
	
	permuting=true;
	//int s = excMax[(*polymerView)[permute->at(0)][0]->type];
	if(s==0)
		s=excMax[0];
	
	//totPerm[(*polymerView)[permute->at(0)][k]->type]++;	
	
	performExchanges(permute,k);
	
	performBrownianBridges(permute,k,s);
	
}
// A console output to display how many exchanges took place.
void MPolymerMoves::displayExchanges()
{
	int snumb=pkinds->species_number;
	for(int i=0;i<snumb;i++)
	{
		if(!typemove)
			continue;
		
		double ratio = ((double)acpPermute[i])/totPerm[i];
		std::cout<<"Exchanges acceptation ratio for "<<pkinds->species[i]<<": "<<ratio<<endl;
		std::cout<<endl<<"******Exchanges details******"<<endl<<endl;
		std::cout<<"\t Particles involved \t \t \t Successful exchanges"<<endl;
		for(int j=0;j<particleNumber;j++)
			if(exchanges[i][j]!=0)
				std::cout<<"\t \t"<<j<<"\t \t \t \t \t "<<exchanges[i][j]<<endl;
	}
}


// Periodic Boundary Conditions method. PBC are usefull to
// quench unwanted surface effects. This function is invoked
// during status-file loading to adapt data for new simbox dimensions.
Point MPolymerMoves::pbc(double x, double y, double z)
{
	Point pbced(dimension);
	double oldp[]={x, y, z};
	double use[dimension];
	for(int i=0;i<dimension;i++){
		double x_trasl = oldp[i] + box_dimensions[i]/2;
		double hmany = x_trasl/box_dimensions[i];
		use[i] = hmany - (int)hmany;
		if(use[i]<0)
			use[i]=1+use[i];
	pbced.coords[i]=box_dimensions[i]*use[i] - box_dimensions[i]/2;
	
	}
	
	
	return pbced;
	
}

// Overloaded pbc: gives the coordinate 'val' in pbc
// system on coordinate j (x=0,y=1,z=2)
double MPolymerMoves::pbc(double val,int j)
{
	double x_trasl= val + box_dimensions[j]/2;
	double hmany = x_trasl/box_dimensions[j];
	double use = hmany - (int)hmany;
	if(use<0)
		use+=1;
	
	return box_dimensions[j]*use-box_dimensions[j]/2;
}

// Overloading pbc function.
Point MPolymerMoves::pbc(Point& coords)
{
	return pbc(coords.coords[0],coords.coords[1],coords.coords[2]);
}

Point MPolymerMoves::pbc(double* coords)
{
	return pbc(coords[0],coords[1],coords[2]);
}



// absolute boundary condition, it is essentially a pbc squared.
double MPolymerMoves::abc(double* c1, double* c2)
{
	double dist=0;
	for(int i=0;i<dimension;i++)
	{

		double delta = c1[i]-c2[i];
		double sq = (box_dimensions[i]/2- fabs(box_dimensions[i]/2 - fabs(delta)));	
		dist+=sq*sq;
	}
	
	return dist;
	
}

double MPolymerMoves::abc(double c1, double c2, int i)
{
	double sq = (box_dimensions[i]/2- fabs(box_dimensions[i]/2 - fabs(c1-c2)));	
	return sq*sq;
	
}

double MPolymerMoves::abc(Point& coords1, Point& coords2)
{
	return abc(coords1.coords,coords2.coords);
}
// takes src first bead coordinates and rewrites them in pbc. Use it as less as you can. It's slow :(
void MPolymerMoves::checkPBC(PLink& src)
{
	
	double coords[dimension];
	for(int i=0;i<dimension;i++)
		coords[i]=src.firstBeadCoordinates[i];
	
	Point np = pbc(coords);
	for(int i=0;i<dimension;i++)
		src.firstBeadCoordinates[i]=np.coords[i];
}

// Straightforward overloading of previous method.
void MPolymerMoves::checkPBC(double* coords)
{
	Point np = pbc(coords);
	for(int i=0;i<dimension;i++)
		coords[i]=np.coords[i];
} 

void MPolymerMoves::checkPBC(double& coords, int i)
{
	coords=pbc(coords,i);
}
// legend: npart -> polymer to cut
//	   ncenter -> link to cut
//         jbb -> B.B length on both sides (if possible)
//	   nhole -> grand canonicle worm. links removed around ncenter. =0 by default (and only possibility for now)
int MPolymerMoves::WORMopen(int npart, int ncenter, int nhole)
{
	if(ZoG==OFFDIAGONAL){
		gnum++;
		return 0;
		
	}
	int type = (*polymerView)[npart][ncenter]->type;
	int jbb=wormM[type];
	jbb =1+ generator->runGen()*(jbb-1);
	
	ZoG=OFFDIAGONAL;
	worming=true;
	//int iendSX=ncenter-nhole/2;
	//int iendDX=ncenter+nhole/2;
	
	
	//int istartSX=iendSX-jbb;
	//int istartDX=iendDX+jbb;
	pmoved=new int;
	*pmoved=npart;
	npmoved=1;
	posworm=ncenter;
	
	(*polymerView)[npart][ncenter]->linked=false;
	nwormed=npart;
	int nprev=ncenter-1;
	if(ncenter==0)
		nprev=(*polymerView)[npart].linkCount()-1;
	
	
	
	double *x1 =(*polymerView)[nwormed][posworm-1-(jbb-1)]->firstBeadCoordinates;
	double k=(*polymerView)[nwormed][posworm-1]->kineticPrefactor;
	
	int nlinks = (*polymerView)[nwormed].linkCount();
		
	
	if(posworm==nlinks-shnumb)  // this works only if SWF is beeing used. Cannot apply worm at variational ends (can't open/close!!)
		varmod=1./PSItry(nwormed,posworm-1,RIGHT);
	
	
	for(int i=0;i<dimension;i++){
		wotemp[i]=(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i];
		(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i]=pbc(generator->gaussGen(x1[i],1./(k/jbb)),i);
	}
	
	if(posworm==nlinks-shnumb)
		varmod*=PSItry(nwormed,posworm-1,RIGHT);
	
	
	if(abc(x1,(*polymerView)[nwormed][posworm]->firstBeadCoordinates)*(k/jbb)>4)
	{
		for(int i=0;i<dimension;i++)
			(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i]=wotemp[i];
		
		(*polymerView)[npart][ncenter]->linked=true;
		worming=false;
		ZoG=DIAGONAL;
		delete pmoved;
		znum++;
		varmod=1;
		return 0;	
	}
	
	
	
/*	int* nopen = new int[particleNumber];
	for(int i=0;i<particleNumber;i++)
		nopen[i]=i;
	

	(*timeSectorView)[nprev].unbindSecMatrix((*timeSectorView)[ncenter].leftMLink(),(*timeSectorView)[ncenter].leftNMLink(),nopen,particleNumber);
	delete nopen;
*/	
	
	//WBB(istartSX,istartDX,npart);
	WBB(posworm-1-(jbb-1),posworm-1,npart);
	//brownianBridge(posworm-1-(jbb-1),posworm-1,npart);
	
	double *x2 =(*polymerView)[nwormed][posworm]->firstBeadCoordinates;
	wclspring=pow(k/(pi*jbb),((double)dimension)/2)*kineticSpring(x1,x2,k/jbb);
	totBB[type]--;
	totWOpen[type]++;
	
	
	
	return 1;
}

int MPolymerMoves::WORMclose()
{
	if(ZoG==DIAGONAL){
		znum++;
		return 0;
		
	}
	
	ZoG=DIAGONAL;
	worming=true;
	int type = (*polymerView)[nwormed][posworm]->type;
	int lnt=wormM[type];
	lnt =1+ generator->runGen()*(lnt-1);
	//int istartSX=posworm-bbMax/2;  // set as worm parameter !!
	//int istartDX=posworm+bbMax/2;
	pmoved=new int;
	*pmoved=nwormed;
	npmoved=1;
	
	double *x2 = (*polymerView)[nwormed][posworm]->firstBeadCoordinates;
	double *x1 = (*polymerView)[nwormed][posworm-1-(lnt-1)]->firstBeadCoordinates;
	double k = (*polymerView)[nwormed][posworm-1]->kineticPrefactor;
	
	
	if(abc(x1,x2)*(k/lnt)>4)
	{
		gnum++;
		ZoG=OFFDIAGONAL;
		worming=false;
		delete pmoved;
		return 0;	
	}
	int nlinks = (*polymerView)[nwormed].linkCount();
		
	
	if(posworm==nlinks-shnumb)  // this works only if SWF is beeing used. Cannot apply worm at variational ends (can't open/close!!)
		varmod=1./PSItry(nwormed,posworm-1,RIGHT);

	
	(*polymerView)[nwormed][posworm]->linked=true;
	
	
	for(int i=0;i<dimension;i++){
		wctemp[i]=(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i];
		(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i]=(*polymerView)[nwormed][posworm]->firstBeadCoordinates[i];
	}	
	
	int nprev=posworm-1;
	if(posworm==0)
		nprev=(*polymerView)[nwormed].linkCount()-1;
	
	//PLink* c1 = (*polymerView)[nwormed][posworm-1];
	//(*timeSectorView)[nprev].setLinked(true);	
	
	if(posworm==nlinks-shnumb)  // this works only if SWF is beeing used. Cannot apply worm at variational ends (can't open/close!!)
		varmod=1./PSItry(nwormed,posworm-1,RIGHT);
	
	
	WBB(posworm-1-(lnt-1),posworm-1,nwormed);
	
	wclspring=pow(k/(pi*lnt),((double)dimension)/2)*kineticSpring(x1,x2,k/lnt);
	totBB[type]--;
	totWClose[type]++;
	//new_IMdistance=0;  // HOLY CRAP !!!!!
	return 1;
}

int MPolymerMoves::WORMterminate(bool opt)
{
	int ans;
	int type = (*polymerView)[nwormed][posworm]->type;
		
	if(swapping || swapping_reverse)
	{
		int wise=0;
		if(swapping_reverse)
			wise=1;
			
		pacc*=swmod;
		ans=polymerTerminateMove(opt);
		if(ans)
		{
			acpSwap[type]++;
			if(swapping)
				nwormed=wormswap->at(1); // worm found a new apple (I'd like to find an Apple, too...)
				
		}
		else
		{
			WORMunswap();
						
		}
		delete wormswap;
		swapping=false;
		swapping_reverse=false;
		return 1;
	
	}
	
	
	worming=false;
	
		
	if(ZoG==OFFDIAGONAL)  // you are trying to open a worm !
	{
		
		pacc *=wormacpconstant*density;
		pacc /=wclspring;
		ans=polymerTerminateMove(opt);
		if(!ans)  // don't open ! 
		{
			// restore old config
			(*polymerView)[nwormed][posworm]->linked=true;
			
			for(int i=0;i<dimension;i++)
				(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i]=wotemp[i];
	
			
			int nprev=posworm-1;
			if(posworm==0)
				nprev=(*polymerView)[nwormed].linkCount()-1;
			
		/*	int* nopen = new int[particleNumber];
			for(int i=0;i<particleNumber;i++)
				nopen[i]=i;
				
			(*timeSectorView)[nprev].bindSecMatrix((*timeSectorView)[posworm].leftMLink(),(*timeSectorView)[posworm].leftNMLink(),nopen,particleNumber);
			//(*timeSectorView)[nprev].unbindSecMatrix((*timeSectorView)[posworm].leftMLink(),(*timeSectorView)[posworm].leftNMLink(),nopen,particleNumber);
			delete nopen;
			*/
			znum++;
			ZoG=DIAGONAL;
		}
		else
		{
			ZoG=OFFDIAGONAL;
			gnum++;
			acpWOpen[type]++;
		}
		
	}
	else  // trying to close....
	{
		pacc /=(wormacpconstant*density);
		pacc *= wclspring;
		ans=polymerTerminateMove();
		int nprev=posworm-1;
			if(posworm==0)
				nprev=(*polymerView)[nwormed].linkCount()-1;
				
		if(!ans)  // don't close !
		{
			// restore old config
			(*polymerView)[nwormed][posworm]->linked=false;
			//(*timeSectorView)[nprev].setLinked(false);
			for(int i=0;i<3;i++)
				(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates[i]=wctemp[i];
				
	
			gnum++;
			ZoG=OFFDIAGONAL;
		}
		else
		{
			ZoG=DIAGONAL;
			
			
			/*int* nopen = new int[particleNumber];
			for(int i=0;i<particleNumber;i++)
				nopen[i]=i;
			cout<<"CLOSING*************************"<<endl;	
			(*timeSectorView)[nprev].bindSecMatrix((*timeSectorView)[posworm].leftMLink(),(*timeSectorView)[posworm].leftNMLink(),nopen,particleNumber);
			//(*timeSectorView)[nprev].unbindSecMatrix((*timeSectorView)[posworm].leftMLink(),(*timeSectorView)[posworm].leftNMLink(),nopen,particleNumber);
			delete nopen;
			*/znum++;
			acpWClose[type]++;
			new_IMdistance=0;
		}
	}
	return ans;

}

// This Brownian Bridge is the same of (and reduces to) BrownianBridge function but is wise enough
// to get aware of "broken links" (worms) and build two indipendent BBs following WORM scheme 
// (determine a bb extremity sampling a gaussian distribution that will be the new tail or head coordinates)
// WBB has same execution time of BrownianBridge method, therefore the latter function won't be used anymore,
// btw, it won't be removed because it could be a starting point for future new implementations.

void MPolymerMoves::WBB(int start, int end, int npart)
{
	int hmany;
	int nlnk=(*polymerView)[npart].linkCount();
	int type = (*polymerView)[npart][start]->type;
	int reverse, icut=-1;
	particle=npart;
	totBB[type]++;
	p1=start;
	p2=end;
	// I can try also to move PIGS variational extremities in the same way.
	// In that case you have to move [start;nlinks-1] U [0;end] where 0 is unlinked.
	
	if(start<=end)   
	{
		hmany = end-start+1;
		reverse=0;
		p1=start;
		p2=end;
	}
	else
	{
		hmany = nlnk-start+end;
		reverse=1;
		p2=start;
		p1=end;
	}
	
	PLink clinks[hmany];
	clinks[0]=*(*polymerView)[npart][start];
	double kp = clinks[0].kineticPrefactor;

	int keyvector[hmany];
	int kn=1;
	keyvector[0]=start;
	// cut polymer in this block !!!
	if(!reverse){
		for(int i=start+1;i<=end;i++)
		{
			clinks[kn]=*clinks[kn-1].pnext;
			keyvector[kn]=i;
			if(clinks[kn++].linked==false)
				icut=kn-2;
		}
	}
	else
	{	
		kn=1;
		for(int i=start+1;i<nlnk;i++){
			clinks[kn]=*clinks[kn-1].pnext;
			keyvector[kn]=i;
			if(clinks[kn++].linked==false)
				icut=kn-2;
			
		}
		clinks[kn++]= *(*polymerView)[npart][0];
		for(int i=1;i<=end;i++){
			clinks[kn]=*clinks[kn-1].pnext;
			keyvector[kn]=i;
			if(clinks[kn++].linked==false)
				icut=kn-2;
			
		}
			
	}
	(*polymerView)[npart].cutPolymer(start,end);
	// now I've just ordered my links and I also know where a cut is (if any)
	// I don't care of effective polymer ordering as long as I map it into a 'keyvector'
	double firstendCoords[dimension];
	double secondstartCoords[dimension];
	double *secondend = clinks[hmany-1].secondBeadCoordinates;
	
	if(icut>=0){
		for(int i=0;i<dimension;i++)
		{
		
			//firstendCoords[i]=pbc(generator->gaussGen(clinks[0].secondBeadCoordinates[i],1./(kp/(icut+1))),i);
			firstendCoords[i]=pbc(generator->gaussGen(clinks[0].firstBeadCoordinates[i],1./(kp/(icut+1))),i);
			clinks[icut].secondBeadCoordinates[i]=firstendCoords[i];
			
			//firstendCoords[i]=clinks[icut].secondBeadCoordinates[i];
			
		
			//secondstartCoords[i]=pbc(generator->gaussGen(clinks[hmany-1].firstBeadCoordinates[i],1./(kp/(hmany-icut+1))),i);
			secondstartCoords[i]=pbc(generator->gaussGen(clinks[hmany-1].secondBeadCoordinates[i],1./(kp/(hmany-icut-1))),i);
			
			//secondstartCoords[i]=clinks[icut+1].firstBeadCoordinates[i];
			clinks[icut+1].firstBeadCoordinates[i]=secondstartCoords[i];
		}
	}
	else
	{
		for(int i=0;i<dimension;i++)
			secondstartCoords[i]=clinks[0].firstBeadCoordinates[i];
	}
	int cn=2;
	int k=1;
	for(int i=0;i<=icut-1;i+=cn) // BB between start and icut. If polymer isn't linked icut=-1... 
	{
		double rappL = ((double)(icut - i + 1))/(icut-i+2);
		double rappR = ((double)(icut - i))/(icut-i+1);
		for(int j=0;j<dimension;j++)
		{
			
			double xprev;
			double xvar;
			double xavg;
			double newcoord;
			
			if (cn==2)  
			{
				xprev=clinks[k-1].firstBeadCoordinates[j];
				xvar= pbc(firstendCoords[j]-xprev,j);
				xavg= xprev + xvar/(icut-i+2);
				newcoord =pbc(generator->gaussGen(xavg,rappL/kp),j);
			
				clinks[k-1].secondBeadCoordinates[j]= newcoord;
				clinks[k].firstBeadCoordinates[j]=newcoord;
				xprev=newcoord;
			}
			else
				xprev=clinks[k].firstBeadCoordinates[j];
			
			xvar= pbc(firstendCoords[j]-xprev,j);
			xavg= xprev + xvar/(icut-i+1);
			newcoord =pbc(generator->gaussGen(xavg,rappR/kp),j);
			clinks[k].secondBeadCoordinates[j]=newcoord;
			
			if(k!=icut)
				clinks[k+1].firstBeadCoordinates[j]=newcoord;
			
		
		
		
		}
		if(i==icut-2)
			cn=1;
		k+=cn;
		
			
	}
	cn=2;
	k=icut+2;
	//cout<<keyvector[0]<<" "<<keyvector[hmany-1]<<endl;
	// second polymer part (if worm doesn't exist this is the only loop performed...)
	for(int i=icut+2;i<=hmany-2;i+=cn)
	{
		//cout<<keyvector[i]<<" ";	
		double rappL = ((double)(hmany - i))/(hmany-i+1);
		double rappR = ((double)(hmany - i-1))/(hmany-i);
		for(int j=0;j<dimension;j++)
		{
			double xprev;
			double xvar;
			double xavg;
			double newcoord;
			
			if (cn==2)
			{
				xprev=clinks[k-1].firstBeadCoordinates[j];
				xvar= pbc(secondend[j]-xprev,j);
				xavg= xprev + xvar/(hmany-i+1);
			
				newcoord =pbc(generator->gaussGen(xavg,rappL/kp),j);
				clinks[k-1].secondBeadCoordinates[j]= newcoord;
				clinks[k].firstBeadCoordinates[j]=newcoord;
				xprev=newcoord;
			}
			else
				xprev=clinks[k].firstBeadCoordinates[j];
			// chain to next
			//if (cn==1 || clinks[k+1].linked==false)
			//	continue;
				
			//xprev=newcoord;
			xvar= pbc(secondend[j]-xprev,j);
			xavg= xprev + xvar/(hmany-i);
			newcoord =pbc(generator->gaussGen(xavg,rappR/kp),j);
			clinks[k].secondBeadCoordinates[j]=newcoord;
			if(clinks[k+1].linked)
				clinks[k+1].firstBeadCoordinates[j]=newcoord;
			
		}
		if(i==hmany-3)
			cn=1;
		k+=cn;	
	}
	// let's insert theese new links !
	for(int i=0;i<hmany;i++)
		(*polymerView)[npart].insertLink(clinks[i],keyvector[i]);
	

}

// PIMC winding number. REALLY an heavy computation. Sure to use it ?
double MPolymerMoves::calculateSquareWindingNumber()
{
	int npart= polymerView->particleNumber();
	PLink* clinks[npart];
	for(int i=0;i<npart;i++)
		clinks[i]=(*polymerView)[i][0];
		
	windX=0;
	windY=0;
	windZ=0;	
		
	for(int i=0;i<npart;i++)
	{
		int nlinks = (*polymerView)[i].linkCount();
		for(int j=0;j<nlinks;j++)
		{
			windX+=sqrt(abc(clinks[j]->firstBeadCoordinates[0],clinks[j]->secondBeadCoordinates[0],0));
			windY+=sqrt(abc(clinks[j]->firstBeadCoordinates[1],clinks[j]->secondBeadCoordinates[1],1));
			windZ+=sqrt(abc(clinks[j]->firstBeadCoordinates[2],clinks[j]->secondBeadCoordinates[2],2));
			clinks[j]=clinks[j]->pnext;
		}
	}
	
	return windX*windX+windY*windY+windZ*windZ;
}

int MPolymerMoves::WORMswap(int wise)
{
	if(ZoG==DIAGONAL)
		return 0;
		
	int nlnk=(*polymerView)[nwormed].linkCount();
	int type =(*polymerView)[nwormed][posworm]->type;
	int jmax= wormM[type];
	int nprev = posworm -1;
	int hmany =/*jmax;//*/1 + generator->runGen()*(jmax-1);
	if(nprev<0)
		nprev+=nlnk;
	
	wormswap = buildWORMPermutationCycle(wise,hmany);
	if(!wormswap->closed)
	{		
		delete wormswap;
		return 0;   // nobody wanna exchange with me :(((
	}
	
		
	int spart = wormswap->at(1);
	totBB[type]--;
	totSwap[type]++;
	pmoved = new int[2];
	npmoved=2;
	pmoved[0]=nwormed;
	pmoved[1]=spart;
	
	swmod=swapkinetic(hmany,spart,wise);
				
	if(wise==0)
	{	
		performExchanges(wormswap,posworm);
		WBB(posworm,posworm+hmany-1,nwormed);
		swapping=true;
	}
	else
	{
		swapping_reverse=true;
		
		(*polymerView)[nwormed][posworm]->linked=true; // switching "holes"
		(*polymerView)[spart][posworm]->linked=false;
		//double scoord[dimension];
		for(int i=0;i<dimension;i++)
			scoord[i]=(*polymerView)[nwormed][posworm]->firstBeadCoordinates[i];
		
		performExchanges(wormswap,posworm);
		//for(int i=0;i<dimension;i++)
		//	(*polymerView)[nwormed][posworm]->firstBeadCoordinates[i]=scoord[i];
		
		for(int i=0;i<dimension;i++){
			(*polymerView)[spart][posworm]->firstBeadCoordinates[i]=scoord[i];
			(*polymerView)[spart][posworm-1]->secondBeadCoordinates[i]=scoord[i];
		}
		WBB(posworm-hmany,posworm-1,spart);
		
	}
	

	
	return 1;
}

void MPolymerMoves::WORMunswap()
{
	
	int nlnk=(*polymerView)[nwormed].linkCount();
	int nprev = posworm -1;
	
	if(nprev<0)
		nprev=nlnk-1;
	
	if(swapping)
		resetExchanges(wormswap,posworm);
	else if(swapping_reverse)
	{
		int spart = wormswap->at(1);
		//double scoord[dimension];
		for(int i=0;i<dimension;i++)
			scoord[i]=(*polymerView)[nwormed][posworm]->firstBeadCoordinates[i];
			
		(*polymerView)[nwormed][posworm]->linked=true;
		(*polymerView)[spart][posworm]->linked=false;
		
		resetExchanges(wormswap,posworm);	
	//	for(int i=0;i<dimension;i++)
	//		(*polymerView)[nwormed][posworm]->firstBeadCoordinates[i]=scoord[i];
		
		for(int i=0;i<dimension;i++){
			(*polymerView)[spart][posworm]->firstBeadCoordinates[i]=scoord[i];
			(*polymerView)[spart][posworm-1]->secondBeadCoordinates[i]=scoord[i];
			
		}
		//resetExchanges(wormswap,0);
		//resetExchanges(wormswap,posworm-1);
		
		
		
	}
	
}

double MPolymerMoves::swapkinetic(int hmany,int spart,int wise)
{
	double s1=0;
	double s2=0;
	int** neighCenter = timeSectorView->returnVerlet()->returnNeighCenter();
	int* neighNumber = timeSectorView->returnVerlet()->returnNeighNumber();
	if(wise==0)
	{
		double *x1=(*polymerView)[nwormed][posworm-1]->secondBeadCoordinates;
		double *x1b=(*polymerView)[spart][posworm]->firstBeadCoordinates;
		double k = (*polymerView)[nwormed][posworm-1]->kineticPrefactor;
		for(int i=0;i<neighNumber[nwormed];i++)
		{
			int ip = neighCenter[nwormed][i];
			if(ip==nwormed)
				continue;
		
			double *x2 = (*polymerView)[ip][posworm+hmany]->firstBeadCoordinates;
			s1+=kineticSpring(x1,x2,k/hmany);
			s2+=kineticSpring(x1b,x2,k/hmany);
		
		}
	}
	else
	{
		double *x1b=(*polymerView)[spart][posworm-1]->secondBeadCoordinates;
		double *x1=(*polymerView)[nwormed][posworm]->firstBeadCoordinates;
		double k = (*polymerView)[spart][posworm-1]->kineticPrefactor;
		for(int i=0;i<neighNumber[nwormed];i++)
		{
			int ip = neighCenter[nwormed][i];
			if(ip==nwormed)
				continue;
		
			double *x2 = (*polymerView)[ip][posworm-hmany]->firstBeadCoordinates;
			s1+=kineticSpring(x1,x2,k/hmany);
			s2+=kineticSpring(x1b,x2,k/hmany);
		
		}
	}
	return s1/s2;
}



void MPolymerMoves::displayMovesAcceptation()
{
	for(int i=0;i<pkinds->species_number;i++)
	{
		if(typemove[i])
		{
			cout<<"Moves acceptance for "<<pkinds->species[i]<<endl<<endl;
			if(dsingmax[i][0]!=0 || dsingmax[i][1]!=0 || dsingmax[i][2]!=0)
				cout<<"Single moves acceptance ratio: "<<acceptedSingleMoves(i)<<endl;
			if(dallmax[i][0]!=0 || dallmax[i][1]!=0 || dallmax[i][2]!=0)
				cout<<"Translation moves acceptance ratio: "<<acceptedTranslations(i)<<endl;
			
			if(bbTry[i]!=0)
				cout<<"Brownian Bridges acceptance ratio: "<<acceptedBB(i)<<endl;
			cout<<"******************"<<endl;
		}
	}
}

double MPolymerMoves::gaussianPSI(int part,int link, int wise)
{
	
	if(gaussC==0)
		return 1;
	double val=0;
	double* r;
	if(wise==LEFT){
		r = (*polymerView)[part][link]->firstBeadCoordinates;
		val = abc(r,PSIs[part]);
	}
	else{
		r = (*polymerView)[part][link]->secondBeadCoordinates;
		val = abc(r,PSILs[part]);
	}
	
	
	return exp(-gaussC*val);

}


void MPolymerMoves::gaussianPSIinit(double gc)
{
	gaussC=gc;
	PSIs = new double*[particleNumber];
	PSILs = new double*[particleNumber];
	
	for(int i=0;i<particleNumber;i++)
	{	
		PSIs[i]=new double[dimension];
		PSILs[i]=new double[dimension];
		for(int j=0;j<dimension;j++){
			PSIs[i][j]=(*polymerView)[i][0]->firstBeadCoordinates[j];
			PSILs[i][j]=(*polymerView)[i][0]->firstBeadCoordinates[j];
		}

	}
}


void MPolymerMoves::setJastrow(double b, double m)
{
	jastrow=true;
	jastrow_b=b;
	jastrow_m=m;
}

double MPolymerMoves::PSItry(int part,int link, int wise)
{
	if(optshadow)
		return PSIOptReal(part,link,wise);
	if(jastrow)
		return PSIJastrow(part,link,wise);
	if(gaussC!=0)
		return gaussianPSI(part,link,wise); 
	
	return 1.;
}

double MPolymerMoves::PSIJastrow(int part, int link, int wise)
{
	double expo=0;
	double* rp;
	if(wise==LEFT)
		rp= (*polymerView)[part][link]->firstBeadCoordinates;
	else
		rp= (*polymerView)[part][link]->secondBeadCoordinates;
	
	
	for(int i=0;i<particleNumber;i++)
	{
		if(i==part)
			continue;
		double *ri;
		if(wise==LEFT)
			ri=(*polymerView)[i][link]->firstBeadCoordinates;
		else
			ri=(*polymerView)[i][link]->secondBeadCoordinates;
		
		expo += pseudo_jastrow(sqrt(abc(rp,ri)));
	}
	
	return exp(-0.5*expo);
}

double MPolymerMoves::PSIOptReal(int part, int link, int wise)
{
	double expo=0;
	double* rp;
	int t1 = (*polymerView)[part][link]->type;
	if(wise==LEFT)
		rp= (*polymerView)[part][link]->firstBeadCoordinates;
	else
		rp= (*polymerView)[part][link]->secondBeadCoordinates;
	
	
	for(int i=0;i<particleNumber;i++)
	{
		if(i==part)
			continue;
		double *ri;
		int t2 =(*polymerView)[i][link]->type;
		if(wise==LEFT)
			ri=(*polymerView)[i][link]->firstBeadCoordinates;
		else
			ri=(*polymerView)[i][link]->secondBeadCoordinates;
		
		
		double r = sqrt(abc(rp,ri));
		double vr = pkinds->itabs[t1][t2]->evaluateOPTShadow(r,0);
		//double vlmr =pkinds->itabs[t1][t2]->evaluateOPTShadow(box_less-r,0);  // cutof at L/2 already taken into account during table generation!
		//double vl2 = pkinds->itabs[t1][t2]->evaluateOPTShadow(box_less/2,0);
		expo +=  vr;// + vlmr - 2*vl2;
	}
	
	return exp(-expo);
}


double MPolymerMoves::pseudo_jastrow(double r)
{
	jashm = expo_jastrow(box_less/2);
	return expo_jastrow(r)+expo_jastrow(box_less-r)-2*jashm;
}


void MPolymerMoves::gotoZSector()
{
	if(rank==0)
		cout<<"Forcing switch to diagonal sector..."<<endl;
	int counter=0;
	wormacpconstant=0.1;
	while(ZoG==OFFDIAGONAL)
	{
		int ans=WORMclose();
		if(ans){
			pacc=1;
			WORMterminate();
		}
		
		polymerView->isDiagonal(ZoG);
		counter++;
		if(counter==ATTEMPTS)
		{
			cout<<"Node "<<rank<<": ARGH! Okay... I give up!"<<endl;
			break;
		}
	}
	if(counter<ATTEMPTS)
		cout<<"Node "<<rank<<": done."<<endl;
}



void MPolymerMoves::swapGaussian(int p1, int p2)
{
	if(gaussC!=0)
	{
		double* v1 = PSILs[p1];
		PSILs[p1]=PSILs[p2];
		PSILs[p2]=v1;
	}
}

double MPolymerMoves::repulsive(double r)
{
	return rep_norm/(1+alfa_rep*exp(-beta_rep*r*r)+gamma_rep*exp(-delta_rep*r));
	
}

void MPolymerMoves::WORMInit(double c)
{
	wormacpconstant=c;
	if(rank==0 && c != 0 && (alfa_rep!=0 || gamma_rep != 0))
		cout<<"A repulsive factor is beeing used !"<<endl;
}
		
// Local kinetic energy estimators

double MPolymerMoves::KElocEstimator(int link, int wise)
{

	if(optshadow)
		KElocSWF(link,wise);
	
	if(shadow)
		return KElocSWF(link,wise);
	if(jastrow)
		return KElocJWF(link,wise);
	if(gaussC!=0)
		return KElocGWF(link,wise); 
	
	return 0.;
}


double MPolymerMoves::KElocJWF(int link, int wise)
{
	if(optshadow)
		return KElocORE(link,wise);
	
	double ksum=0;
	double* masses = pkinds->masses;
	for(int i=0;i<particleNumber;i++)
	{
		double* ri;
		if(wise==LEFT)
			ri= (*polymerView)[i][link]->firstBeadCoordinates;
		else
			ri= (*polymerView)[i][link]->secondBeadCoordinates;
		
		double fac= -hbar*hbar/(2*masses[(*polymerView)[i][link]->type]);
		fac/=boltzk;
		fac*=metric_conversion;
		
		double keloc[dimension];
		for(int l=0;l<dimension;l++)
			keloc[l]=0;
		double keloc4=0;
	
		
		for(int j=0;j<particleNumber;j++)
		{
			
			if(i==j)
				continue;
				
			double* rj;
			if(wise==LEFT)
				rj= (*polymerView)[j][link]->firstBeadCoordinates;
			else
				rj= (*polymerView)[j][link]->secondBeadCoordinates;
		
			double r = sqrt(abc(ri,rj));
			
			if(r>box_less/2)
				continue;

			
			double vprime =	-jastrow_m*pow(jastrow_b/r,jastrow_m)/r;
			vprime+=jastrow_m*pow(jastrow_b/(box_less-r),jastrow_m)/(box_less-r);
			
			
			for(int l=0;l<dimension;l++)
				keloc[l]+=vprime*pbc(ri[l]-rj[l],l)/r;
			
			if(j<i)
			{
				double vsec=jastrow_m*(jastrow_m+1)*(pow(jastrow_b/r,jastrow_m)/(r*r));
				vsec+=jastrow_m*(jastrow_m+1)*(pow(jastrow_b/(box_less-r),jastrow_m)/((box_less-r)*(box_less-r)));
				keloc4+=vsec;
				keloc4+=(dimension-1)*vprime/r;
			}
		}
		for(int l=0;l<dimension;l++)
			ksum+=fac*keloc[l]*keloc[l]/4;
			
		ksum-= fac*keloc4;
	
	}
		
		
	
	return ksum/particleNumber;
	
}

double MPolymerMoves::KElocSWF(int link, int wise)
{
	int nlnk=(*polymerView)[0].linkCount();
	if(link==0)
		return shelocL;
	if(link==nlnk-1 && wise == RIGHT)
		return shelocR;
	
	return 0;
}

double MPolymerMoves::KElocGWF(int link, int wise)
{
	double el=0;
	for(int i=0;i<particleNumber;i++)
	{
		double fac= -hbar*hbar/(2*pkinds->masses[(*polymerView)[i][link]->type]);
		fac/=boltzk;
		fac*=metric_conversion;
		
		double* ri;
		double* si;
		if(wise==LEFT)
		{
			ri= (*polymerView)[i][link]->firstBeadCoordinates;
			si= PSIs[i];
		}
		else
		{
			ri= (*polymerView)[i][link]->secondBeadCoordinates;
			si=PSILs[i];
		}
		el += fac*2*gaussC*(2*gaussC*abc(ri,si)-dimension);
	}
	
	return el/particleNumber;	
}

double MPolymerMoves::KElocSH(int link, int wise)
{
	
	if(optshadow)
		return KElocOSH(link,wise);
	
	double ksum=0;
	double shC = polymerView->sh_C;
	for(int k=0;k<particleNumber;k++)
	{
		double* rk;
		double* sk;
		if(wise==LEFT){
			rk= (*polymerView)[k][link]->secondBeadCoordinates;
			sk= (*polymerView)[k][link]->firstBeadCoordinates;
		}
		else{
			rk= (*polymerView)[k][link]->firstBeadCoordinates;
			sk= (*polymerView)[k][link]->secondBeadCoordinates;
			
		}
		double fac= -hbar*hbar/(2*pkinds->masses[(*polymerView)[k][link]->type]);
		fac/=boltzk;
		fac*=metric_conversion;
		
		double cross[dimension];
		
		for(int l=0;l<dimension;l++)
			cross[l]=0;
		
	
		
		for(int j=0;j<particleNumber;j++)
		{
			
			if(k==j)
				continue;
				
			double* rj;
			if(wise==LEFT)
				rj= (*polymerView)[j][link]->secondBeadCoordinates;
			else
				rj= (*polymerView)[j][link]->firstBeadCoordinates;
		
			double r = sqrt(abc(rk,rj));
			
			if(r>box_less/2)
				continue;
			
			double vprime=-jastrow_m*pow(jastrow_b/r,jastrow_m)/r;
			vprime+=jastrow_m*pow(jastrow_b/(box_less-r),jastrow_m)/(box_less-r);
			
			
			for(int l=0;l<dimension;l++)
				cross[l]+=vprime*pbc(rk[l]-rj[l],l)/r;
		}
		
		
		for(int l=0;l<dimension;l++)
			ksum+=2*shC*fac*cross[l]*pbc(rk[l]-sk[l],l);
			
			
	
	ksum+= fac*2*shC*(2*shC*abc(rk,sk)-dimension);
	}
	
	return ksum/particleNumber;


}


double MPolymerMoves::KElocOSH(int link, int wise)
{
	
	double ksum=0;
	double accoM[particleNumber][dimension];
	double* masses = pkinds->masses;
	
	for(int k=0;k<particleNumber;k++)
	{
		double* rk;
		double* sk;
		if(wise==LEFT){
			rk= (*polymerView)[k][link]->secondBeadCoordinates;
			sk= (*polymerView)[k][link]->firstBeadCoordinates;
		}
		else{
			rk= (*polymerView)[k][link]->firstBeadCoordinates;
			sk= (*polymerView)[k][link]->secondBeadCoordinates;
			
		}
		int t1 =(*polymerView)[k][link]->type;
		
		
		double val[dimension];
			for(int i=0;i<dimension;i++)
				val[i]=0;
		
		double r = sqrt(abc(rk,sk));
				
		double fsp = pkinds->itabs[t1][t1]->evaluateOPTShadow(r,1);
		double fsp_prime = pkinds->itabs[t1][t1]->evaluateOPTShadow(r,5);
		double fsp_sec = pkinds->itabs[t1][t1]->evaluateOPTShadow(r,6);
		
		double acco = fsp_prime/(r*fsp);
		for(int d=0;d<dimension;d++)
			accoM[k][d]=acco*pbc(rk[d]-sk[d],d);
		
		for(int j=0;j<k;j++)
		{
			double* rj;
			double fr_prime;
			double pbcmy[dimension];
			double rik;
			if(wise==LEFT){
				rj= (*polymerView)[j][link]->secondBeadCoordinates;
				fr_prime = frprimel[k][j];
				for(int d=0;d<dimension;d++)
					pbcmy[d]=pbcl[k][j][d];
				rik = abcl[k][j];
			}
			else{
				rj= (*polymerView)[j][link]->firstBeadCoordinates;
				fr_prime = frprimer[k][j];
				for(int d=0;d<dimension;d++)
					pbcmy[d]=pbcr[k][j][d];
				rik = abcr[k][j];
			}
		
			
			
			if(rik>box_less/2)
				continue;
			
			int t2 = (*polymerView)[j][link]->type;
			
			
			for(int d=0;d<dimension;d++)
				val[d]+=(accoM[k][d]-accoM[j][d])*fr_prime*pbcmy[d]/rik;
		
		}
				
		
		double fac= -hbar*hbar/(2*masses[(*polymerView)[k][link]->type]);
		fac/=boltzk;
		fac*=metric_conversion;
		
		for(int d=0;d<dimension;d++)
			ksum+= -2*fac*val[d];
		
		
		ksum+= fac*(dimension-1)*fsp_prime/(fsp*r);
		
		ksum+=fac*fsp_sec/fsp;
				
		
		
	}
	
	return ksum/particleNumber;
	
	
	
	
	
	
	/*
	
	
	double ksum=0;
	double accoM[particleNumber][dimension];
	double* masses = pkinds->masses;
	
	for(int k=0;k<particleNumber;k++)
	{
		double* rk;
		double* sk;
		if(wise==LEFT){
			rk= (*polymerView)[k][link]->secondBeadCoordinates;
			sk= (*polymerView)[k][link]->firstBeadCoordinates;
		}
		else{
			rk= (*polymerView)[k][link]->firstBeadCoordinates;
			sk= (*polymerView)[k][link]->secondBeadCoordinates;
			
		}
		int t1 =(*polymerView)[k][link]->type;
		
		
		double val[dimension];
			for(int i=0;i<dimension;i++)
				val[i]=0;
		
		double r = sqrt(abc(rk,sk));
				
		double fsp = pkinds->itabs[t1][t1]->evaluateOPTShadow(r,1);
		double fsp_prime = pkinds->itabs[t1][t1]->evaluateOPTShadow(r,5);
		double fsp_sec = pkinds->itabs[t1][t1]->evaluateOPTShadow(r,6);
		
		double acco = fsp_prime/(r*fsp);
		for(int d=0;d<dimension;d++)
			accoM[k][d]=acco*pbc(rk[d]-sk[d],d);
		
		for(int j=0;j<k;j++)
		{
			double* rj;
			if(wise==LEFT)
				rj= (*polymerView)[j][link]->secondBeadCoordinates;
			else
				rj= (*polymerView)[j][link]->firstBeadCoordinates;
		
			double rik = sqrt(abc(rk,rj));
			
			if(rik>box_less/2)
				continue;
			
			int t2 = (*polymerView)[j][link]->type;
			
			double fr_prime = pkinds->itabs[t1][t2]->evaluateOPTShadow(rik,3);
			
			for(int d=0;d<dimension;d++)
				val[d]+=(accoM[k][d]-accoM[j][d])*fr_prime*pbc(rk[d]-rj[d],d)/rik;
		
		}
				
		
		double fac= -hbar*hbar/(2*masses[(*polymerView)[k][link]->type]);
		fac/=boltzk;
		fac*=metric_conversion;
		
		for(int d=0;d<dimension;d++)
			ksum+= -2*fac*val[d];
		
		
		ksum+= fac*(dimension-1)*fsp_prime/(fsp*r);
		
		ksum+=fac*fsp_sec/fsp;
				
		
		
	}
	
	return ksum/particleNumber;
*/

}


double MPolymerMoves::KElocORE(int link, int wise)
{
	
	double ksum=0;
	
	for(int i=0;i<particleNumber;i++)
	{
		double* ri;
		if(wise==LEFT)
			ri= (*polymerView)[i][link]->firstBeadCoordinates;
		else
			ri= (*polymerView)[i][link]->secondBeadCoordinates;
		
		double fac= -hbar*hbar/(2*pkinds->masses[(*polymerView)[i][link]->type]);
		int t1= (*polymerView)[i][link]->type;
		fac/=boltzk;
		fac*=metric_conversion;
		
		double keloc1=0;
		double squad[dimension];
		for(int d=0;d<dimension;d++)
			squad[d]=0;	
		
		for(int j=0;j<particleNumber;j++)
		{
			
			if(i==j)
				continue;
				
			double* rj;
			int t2= (*polymerView)[j][link]->type;
		
			if(wise==LEFT)
				rj= (*polymerView)[j][link]->firstBeadCoordinates;
			else
				rj= (*polymerView)[j][link]->secondBeadCoordinates;
		
			double r = sqrt(abc(ri,rj));
			
			if(r>box_less/2)
				continue;
			
			double fr_prime=pkinds->itabs[t1][t2]->evaluateOPTShadow(r,3);
			double fr_sec=pkinds->itabs[t1][t2]->evaluateOPTShadow(r,4);
			
			

			for(int d=0;d<dimension;d++)
				squad[d]+=fr_prime*pbc(ri[d]-rj[d],d)/r;
			
			keloc1+=(-fr_sec - (dimension-1)*fr_prime/r);
			
		
		}
			
		for(int d=0;d<dimension;d++)
			ksum+=fac*squad[d]*squad[d];
			
		ksum+= fac*keloc1;
	
	}
		
		
	
	return ksum/particleNumber;
	
}
