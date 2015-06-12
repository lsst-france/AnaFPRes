/*  -----------------------------------------
  class and functions to scan and analyze forced photometry output files - 
  R. Ansari - September 2013 
---------------------------------------------------- */

#ifndef FPFSCANB_H_SEEN
#define FPFSCANB_H_SEEN

//-- C++ std include files 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

//-- SOPHYA include files 
#include "sopnamsp.h"
#include "histats.h"
#include "fitsioserver.h"

using namespace std;
using namespace SOPHYA;
// namespace ANALSSTDC {
// -------------------------------------------------------
// ---- 
// -------------------------------------------------------

//! useful structure to represent a portion of sky
typedef struct radecst {double ramin,ramax; double decmin,decmax; int nbra,nbdec; double delra, deldec; } RADecLim; 
//! a structure to represent an SDSS run with a range of field numbers  
typedef struct runfminmax {int run, fieldmin, fieldmax; } RunFieldMinMax;
//! a structure to represent a given SDSS [run field camcol] 
typedef struct runfieldcc {int run, field, camcol; } RunFieldCamCol;


//-------- Quelques fonctions utiles ----------
//! Convert filter index (1...6) to the SDSS filter character 
// u->1 , g->2 , r->3 , i->4 , z->5 , y->6
inline char FilterId2Filter(int f)
{
  char filts[6]={'u','g','r','i','z','y'};
  if ((f>0)&&(f<7))  return filts[f-1];
  return '?';
}

//! convert camcol character to number  
// 1->'1' 2->'2' ... 6->'6'
inline char CamColId2CamCol(int camcol)
{  
  return ((char)camcol+'0');
}

/*! 
  return the filter number of index (1...6) for the SDSS filter character. 
  return u->1 g->2 r->3 i->4 z->5 y->6, 0 otherwise 
*/
inline int FilterId(char filter) 
{
  char filts[6]={'u','g','r','i','z','y'};
  for(int i=0; i<6; i++)
    if (filter==filts[i])  return (i+1);
  return 0;
}

//! return 1,2,3,4,5,6  for camcol='1','2'... 
inline int CamColId(char camcol) 
{
  int ccid=(int)camcol-(int)'0';
  if ((ccid<1)||(ccid>6))  return 0;
  return ccid;
  }

/*!
  base class to scan a set of ForcedPhotometry output fits files produced by the LSST stack.
  The pure virtual methods PrepareScan() and ProcessFile() should be implemented by the derived 
  classes.
*/ 
class ForcedPhotFileScanBase {
public:
  ForcedPhotFileScanBase(string const & basedir, vector<RunFieldMinMax> const& runlist, 
			 vector<int> const& fieldlist, vector<char> const& camcols, vector<char> const& filters);
  ForcedPhotFileScanBase(string const & basedir, vector<RunFieldCamCol> const& runfcclist, vector<char> const& filters);

  virtual ~ForcedPhotFileScanBase()  {  }

  inline void SetPrintLevel(int lev) { prtlev_=lev; }
  inline void SetOutPath(string const& outp) { outp_=outp; }

  inline void SetRADecLimits(RADecLim const& radec) 
  { 
    radec_=radec; 
    radec_.delra=(radec.ramax-radec.ramin)/(double)radec.nbra; 
    radec_.deldec=(radec.decmax-radec.decmin)/(double)radec.nbdec; 
  }

  inline void SetFluxLimits(double flxmin=0., double flxmax=5.e5)
    { minflxok_=flxmin; maxflxok_=flxmax; }

  inline size_t NbScannedFiles() { return totnfiles_; }
  inline size_t NbMissingFiles() { return totmissingfiles_; }
  inline size_t NbBadFiles() { return totbadfiles_; }
  inline size_t TotalFileSize() { return totfilesize_; }

  virtual int ScanFiles(); 

protected:
  //! Method called just before starting file scan 
  virtual void PrepareScan()=0;
  //! Method called for each existing FP file scanned, with non zero file size 
  virtual int ProcessFile(const char * filename, char filter, char camcol, int field, int run)=0;
  //! Method called at the end of the scan (to save results ...)
  virtual void FinalizeScan()=0;

  void ScanFilesA();   // boucle a partir de runlist_, fieldlist_, camcols_
  void ScanFilesB();   // boucle a partir de runfcclist_

  string basedir_;
  vector<RunFieldMinMax> runlist_;
  vector<int> fieldlist_;   // this fieldlist_ is used if fieldmin or fieldmax in RunFieldMinMax is negative
  vector<char> camcols_;

  vector<RunFieldCamCol> runfcclist_;
  bool fg_use_rfcc;  // true -> use runfcclist instead of runlist_

  vector<char> filters_;

  int prtlev_;

  RADecLim radec_;
  double minflxok_, maxflxok_;

  string outp_;     // output path 

  size_t totnfiles_;
  size_t totmissingfiles_;
  size_t totbadfiles_;
  size_t totfilesize_;
};

/*!
  Class to scan a set of ForcedPhotometry output fits files produced by the LSST stack.
  It can be used to identify missing or bad files, and make a FITS table with summary information.
*/ 

class ForcedPhotFileScan : public  ForcedPhotFileScanBase {
public:
  ForcedPhotFileScan(string const & basedir, vector<RunFieldMinMax> const& runlist, vector<int> const& fieldlist,
		     vector<char> const& camcols, vector<char> const& filters)
    : ForcedPhotFileScanBase(basedir, runlist, fieldlist, camcols, filters), dtrunfield_(256), fg_fill_dtrf_(false)
  {
  }
  ForcedPhotFileScan(string const & basedir, vector<RunFieldCamCol> const& runfcclist, vector<char> const& filters)
    : ForcedPhotFileScanBase(basedir, runfcclist, filters), dtrunfield_(256), fg_fill_dtrf_(false)
  {
  }

  virtual ~ForcedPhotFileScan()
  {
  }

  //! To activate filling of the table with per file summary information 
  inline void ActivateFPFDTFill() { fg_fill_dtrf_=true; }

  inline DataTable& GetSummaryTable() { return dtrunfield_; } 

protected:
  //! Method called just before starting file scan 
  virtual void PrepareScan();
  //! Method called for each existing FP file scanned, with non zero file size 
  virtual int ProcessFile(const char * filename, char filter, char camcol, int field, int run);
  //! Method called at the end of the scan (to save results ...)
  virtual void FinalizeScan();

  map<int, RunFieldMinMax> rfm_;   // list of run number with minimum and maximum valid field values 
  DataTable dtrunfield_; 
  bool fg_fill_dtrf_;  // true -> fill DataTable dtrunfield_ and make (run & field) count sky map
};



// }  Fin du namespace

#endif
