/*  Main program to scan forced photometry output files - 
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
#include "ctimer.h"   // for using class Timer
#include "timing.h"   // for using PrtTim() and InitTim()

#include "sopnamsp.h"

#include "histinit.h"   // for insuring proper sophya initialization  

//----- include files specific to this program 
#include "analsstdc.h"

//------------------------------------------------------
//------------------  Main program ---------------------
//------------------------------------------------------
int main(int narg, char* arg[])
{
  // We handle exception at the high level
  try {
  // This macro initialize the library
  // static objects handle this - However, not all loader call
  // the constructor for static objects
  SophyaInit();
  InitTim();
  Timer tm("fpfscan.cc");
  cout << " ------------ Program fpfscan.cc   ------------ " << endl;
  
  string action, outpath, basedir;
  vector<RunFieldMinMax> runlist;
  vector<RunFieldCamCol> runfcclist;
  vector<int> fieldlist;
  vector<char> camcols, filters;
  RADecLim radec;
  int prtlev;
  pair<double,double> fluxcut;
  int rc=Ana_LSST_DC_DecodeArgs(narg, arg, action, basedir, runlist, runfcclist, fieldlist, camcols, filters, 
				radec, outpath, prtlev, fluxcut);
  if (rc!=0) {
    cout<<"fpfscan non zero RC from DecodeArgs() -> exit"<<endl;
    return rc;
  }

  if ((action=="LIST")||(action=="SUMMFT")) {
    ForcedPhotFileScan fps(basedir, runlist, fieldlist, camcols, filters);
    ForcedPhotFileScan* fpsp=&fps;
    ForcedPhotFileScan fps2(basedir, runfcclist, filters);
    if (runfcclist.size()>0) {
      cout << " fpfscan/info: ForcedPhotFileScan using [run field camcol] triplet list ... " << endl;
      fpsp=&fps2;
    }
    else cout << " fpfscan/info: ForcedPhotFileScan using runn field camcol list ... " << endl;
    fpsp->SetPrintLevel(prtlev);
    fpsp->SetOutPath(outpath);
    if ((radec.nbra>0)&&(radec.nbdec>0)) fpsp->SetRADecLimits(radec);
    if (fluxcut.second>fluxcut.first)  fpsp->SetFluxLimits(fluxcut.first, fluxcut.second);
    if (action == "SUMMFT")  {
      cout << " fpfscan/info: Activating Summary info table fill -> reading forced-photometry output files ..." << endl;
      fpsp->ActivateFPFDTFill();
    }
    fpsp->ScanFiles();
  }
  else if (action == "SRCLIST") {
    FPObjectListBuilder fpolb(basedir, runlist, fieldlist, camcols, filters);
    FPObjectListBuilder* fplp=&fpolb;
    FPObjectListBuilder fpolb2(basedir, runfcclist, filters);
    if (runfcclist.size()>0) {
      cout << " fpfscan/info: FPObjectListBuilder using [run field camcol] triplet list ... " << endl;
      fplp=&fpolb2;
    }
    else cout << " fpfscan/info: FPObjectListBuilder using runn field camcol list ... " << endl;

    fplp->SetPrintLevel(prtlev);
    fplp->SetOutPath(outpath);
    if ((radec.nbra>0)&&(radec.nbdec>0)) fplp->SetRADecLimits(radec);
    if (fluxcut.second>fluxcut.first)  fplp->SetFluxLimits(fluxcut.first, fluxcut.second);

    fplp->ScanFiles();
  }
  else cout << "fpfscan/ERROR: unknown action requested ACT=" << action << endl;
  }
  catch (std::exception & exc) {
    cerr << " fpfscan.cc: Catched Exception " << (string)typeid(exc).name() 
	 << " - Msg= " << exc.what() << endl;
  }
  catch (...)  {
    cerr << " fpfscan.cc: some other exception was caught ! " << endl;
  }
  cout << " ----------------------- END of fpfscan.cc ------------------------ " << endl;
  return 0;
}

