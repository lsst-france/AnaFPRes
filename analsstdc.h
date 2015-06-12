/*  -----------------------------------------
  class and functions to scan and analyze forced photometry output files - 
  R. Ansari - September 2013 
---------------------------------------------------- */

#ifndef ANALSSTDC_H_SEEN
#define ANALSSTDC_H_SEEN

//-- C++ std include files 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

//-- this module include files (will include SOPHYA include files )
#include "fpfscanb.h"
#include "fpfscanmes.h"

using namespace std;
using namespace SOPHYA;
// namespace ANALSSTDC {

void Ana_LSST_DC_Usage();
int Ana_LSST_DC_DecodeArgs(int narg, char* arg[], string& action, string& basedir, 
			   vector<RunFieldMinMax>& runlist, vector<RunFieldCamCol>& runfcclist, 
			   vector<int>& fieldlist, vector<char>& camcols, vector<char>& filters,  
			   RADecLim& radec, string& outfile, int& prtlev, pair<double,double> & fluxcut);
int Ana_LSST_DC_read_listfile(const char* filename, vector<RunFieldMinMax>& list);
int Ana_LSST_DC_read_listfile2(const char* filename, vector<RunFieldCamCol>& list);


// }  Fin du namespace

#endif
