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


//------------------------------------------------------
//------------------  Main program ---------------------
//------------------------------------------------------
int main(int narg, char* arg[])
{

if (narg<4) {
  cout << " mergeft : merging FITS tables into a single file/table \n"
       << "   Usage: mergeft file1 file2 [file3 file4 ...] outfile " << endl;
  return 1;
 }
// We handle exception at the high level
try {
  // Some initialization 
  SophyaInit();
  InitTim();
  Timer tm("mergeft.cc");

  vector<string> infiles;
  string outfile;
  for(int i=1; i<narg-1; i++) infiles.push_back(arg[i]);
  outfile=arg[narg-1];
  cout << " ------------ mergeft: reading " << infiles.size() << " input files ..." << endl; 
  DataTable dta;
  for(size_t i=0; i<infiles.size(); i++) {
    DataTable dt;
    FitsInOutFile fis(infiles[i], FitsInOutFile::Fits_RO);
    fis >> dt; 
    cout << i << "- read FITS file: " << infiles[i] << " -> NRows=" << dt.NRows() << endl;
    dta.CopyMerge(dt);   // merge rows into the dta table 
  }
  cout << " Saving merged table to file: " << outfile << endl;
  FitsInOutFile fos(outfile, FitsInOutFile::Fits_Create);
  fos << dta;
  dta.SetShowMinMaxFlag(true);
  cout << dta;
}
 catch (std::exception & exc) {
   cerr << " mergeft.cc: Catched Exception " << (string)typeid(exc).name() 
	<< " - Msg= " << exc.what() << endl;
   return 98;
 }
 catch (...)  {
   cerr << " mergeft.cc: some other exception was caught ! " << endl;
   return 99;
 }
 cout << " ----------------------- END of mergeft.cc ------------------------ " << endl;
 return 0;
}

