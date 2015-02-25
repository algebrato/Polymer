#include "RNumgen.h"
#include <cmath>
#include <iostream>
RNumgen::RNumgen(int *seeds, int set1, int set2, bool which)
{
	step=1./4096;
	numb=4096;
	m[0]=502;
	m[1]=1521;
	m[2]=4071;
	m[3]=2107;
	current=0;
	for(int i=0;i<4;i++)
	{
		l[i]=seeds[i];
		keys[i]=seeds[i];
	}
	n[0]=0;
	n[1]=0;
	if(which)
	{
		n[2]=set1;
		n[3]=set2;
	}
	else
	{
		n[2]=2896;
		n[3]=1263;
	}	
	
	//inf = new std::ifstream("rnumbers.rnd");
	
}
// returns uniform distributed random number
double RNumgen::runGen()
{
  int ik[4];
  ik[0]=l[0]*m[3]+l[1]*m[2]+l[2]*m[1]+l[3]*m[0]+n[0];
  ik[1]=l[1]*m[3]+l[2]*m[2]+l[3]*m[1] + n[1];
  ik[2]=l[2]*m[3]+l[3]*m[2]+n[2];
  ik[3]=l[3]*m[3] + n[3];
  l[3]=ik[3]%numb;
  ik[2] += ik[3]/numb;
  l[2]=ik[2]%numb;
  ik[1]+=ik[2]/numb;
  l[1]=ik[1]%numb;
  l[0]=(ik[0]+ik[1]/numb)%numb; 
  current=step*(l[0]+step*(l[1]+step*(l[2]+step*l[3])));
  return current; 
}

// This method saves seeds on a file.
void RNumgen::flushGen()
{
  std::fstream  outafile("random.out", std::ios::trunc);
  for(int i=0;i<4;i++)
	  {
		keys[i]=l[i];
		outafile << "Seed " << i+1 <<" = "<<keys[i]<<"\n";
	  }
  outafile.close();
}

// This method runs Box-Muller scheme to yeld a gaussian distributed
// random value. Be careful: variance is 2sigma !!!
double RNumgen::gaussGen(double avg, double variance)
{
	double x1 = runGen();
	double x2= runGen();
	double dpi = 2*acos(-1);
	return avg + sqrt(-variance*log(1-x1))*cos(dpi*x2);
}

