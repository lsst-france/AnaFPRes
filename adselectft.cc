/*  Main program to several fits files containing a table with 
    the same structure into a single table and writes it into 
    an output fits file
    R. Ansari - September 2013 */

//-------- standard C++ includes 
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>

//-------- SOPHYA include 
#include "array.h"    // for using arrays 
#include "histats.h"   // for using histograms and datatables 
#include "fitsioserver.h"   // FITS files ...

#include "ctimer.h"   // for using class Timer
#include "timing.h"   // for using PrtTim() and InitTim()

#include "sopnamsp.h"

#include "histinit.h"   // for insuring proper sophya initialization  

#include "analsstdc.h"   // useful functions

//------------------------------------------------------
//------------------  Main program ---------------------
//------------------------------------------------------
int main(int narg, char* arg[])
{

if (narg<5) {
  cout << " adselectft : ra,dec region select from FITS tables produced by fpfscan \n"
       << "   Usage: adselectft ra_min,max dec_min,max inputFits outlistfile \n" 
       << "      ra,dec values in degree \n" << endl;
  return 1;
 }
// We handle exception at the high level
try {
  // Some initialization 
  SophyaInit();
  InitTim();
  Timer tm("adselectft.cc");

  double ramin=0.,ramax=0.;
  sscanf(arg[1],"%lg,%lg",&ramin,&ramax);
  double decmin=0.,decmax=0.;
  sscanf(arg[2],"%lg,%lg",&decmin,&decmax);
  string infile=arg[3];
  string outfile=arg[4];
  cout << " ------------ adselectft: ra_min,max= " << ramin<<","<<ramax 
       << " dec_min,max= " << decmin<<","<<decmax << endl; 
  cout << " Input fits file= " << infile << " Output list file= " << outfile << endl;
  ofstream of(outfile.c_str());
  DataTable dta;
  FitsInOutFile fis(infile, FitsInOutFile::Fits_RO);
  fis >> dta; 
  cout << " read FITS file -> NRows=" << dta.NRows() << endl;
  DataTableRowPtr rowp=dta.EmptyRowPtr();
  double ra0,ra1,dec0,dec1;
  int_8 nst, count=0;
  int run, field, camcolfilt;
  char camcol, filter;
  for(size_t k=0; k<dta.NRows(); k++) {
    dta.GetCstRowPtr(k,rowp);
    nst=rowp(0);
    if (nst<1) continue;
    run=rowp(0);    field=rowp(1);
    camcolfilt=rowp(2);
    camcol=CamColId2CamCol(camcolfilt/10); 
    filter=FilterId2Filter(camcolfilt%10);
    ra0=rowp(4)[0];    ra1=rowp(4)[1]; 
    dec0=rowp(5)[0];   dec1=rowp(5)[1]; 
    bool fgs=false;
    if ((ra0>=ramin)&&(ra0<=ramax)&&(dec0>=decmin)&&(dec0<=decmax))  fgs=true;
    else if ((ra0>=ramin)&&(ra0<=ramax)&&(dec1>=decmin)&&(dec1<=decmax))  fgs=true;
    else if ((ra1>=ramin)&&(ra1<=ramax)&&(dec0>=decmin)&&(dec0<=decmax))  fgs=true;
    else if ((ra1>=ramin)&&(ra1<=ramax)&&(dec1>=decmin)&&(dec1<=decmax))  fgs=true;
    else if ((ramin>=ra0)&&(ramin<=ra1)&&(decmin>=dec0)&&(decmin<=dec1))  fgs=true;
    else if ((ramin>=ra0)&&(ramin<=ra1)&&(decmax>=dec0)&&(decmax<=dec1))  fgs=true;
    else if ((ramax>=ra0)&&(ramax<=ra1)&&(decmin>=dec0)&&(decmin<=dec1))  fgs=true;
    else if ((ramax>=ra0)&&(ramax<=ra1)&&(decmax>=dec0)&&(decmax<=dec1))  fgs=true;
    else fgs=false;

    if (fgs) {
      of << run << ' ' << field << ' ' << camcol << ' ' << filter << endl;
      count++;
    }
  }	
  cout << "===== adselect END of " << infile  << " scan ==== \n" 
       << "  ---> " << count << " [run field camccol filter] selected (outfile= " << outfile << ")" << endl;
}
 catch (std::exception & exc) {
   cerr << " adselect.cc: Catched Exception " << (string)typeid(exc).name() 
	<< " - Msg= " << exc.what() << endl;
   return 98;
 }
 catch (...)  {
   cerr << " adselect.cc: some other exception was caught ! " << endl;
   return 99;
 }
 return 0;
}

