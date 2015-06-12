/*  -----------------------------------------
  class and functions to scan and analyze forced photometry output files - 
  R. Ansari - September 2013 
---------------------------------------------------- */

#ifndef FPFSCANMES_H_SEEN
#define FPFSCANMES_H_SEEN

//-- C++ std include files 
#include "fpfscanb.h"


using namespace std;
using namespace SOPHYA;
// namespace ANALSSTDC {

//! utility structure to hold individual (or sum) flux measurement
typedef struct fluxrec {double sum_mean, sqsum_sigma; int nmes; } FluxRec;

//! utility class to hold agregated multi-color information extracted from forced-photometry fits files for a given object
class FPMesures {
public:
  //! constructor
  FPMesures(int_8 id, int_8 oid, double ra, double dec, int nfilt) 
    : id_(id), oid_(oid), ra_(ra), dec_(dec), vsflx_(nfilt) 
  { 
    for(size_t i=0; i<vsflx_.size(); i++) {
      vsflx_[i].sum_mean=vsflx_[i].sqsum_sigma=0.; vsflx_[i].nmes=0; 
    }
  }
  //! default constructor
  FPMesures() 
    : id_(0), oid_(0), ra_(0.), dec_(0.) 
  { 
  }
  //! copy constructor
  FPMesures(FPMesures const& a)
  : id_(a.id_), oid_(a.oid_), ra_(a.ra_), dec_(a.dec_), vsflx_(a.vsflx_)
  {
  } 
  inline void AddMes(int jf, double flux)   // Attention - non protege pour bornes de jf 
  { 
    vsflx_[jf].sum_mean += flux; vsflx_[jf].sqsum_sigma += (flux*flux);  vsflx_[jf].nmes++; 
  }

  void ComputeMean()
  {
    for(size_t i=0; i<vsflx_.size(); i++) {
      if (vsflx_[i].nmes < 1) continue;
      vsflx_[i].sum_mean /= (double)vsflx_[i].nmes;  
      vsflx_[i].sqsum_sigma = sqrt(vsflx_[i].sqsum_sigma/(double)vsflx_[i].nmes - 
				      vsflx_[i].sum_mean*vsflx_[i].sum_mean);
    }
  }
  ostream& Print(ostream& os, const char* sep=" ") const {
    os << id_ << sep << oid_ << sep << ra_ << sep << dec_ << sep;
    for(size_t i=0; i<vsflx_.size(); i++) 
      os << vsflx_[i].sum_mean << sep << vsflx_[i].sqsum_sigma << sep << vsflx_[i].nmes << sep;
    return os;
  }
  int_8 id_,oid_;
  double ra_,dec_;
  vector<FluxRec> vsflx_;
};

inline ostream& operator << (ostream& os, FPMesures const& a) { return a.Print(os); }


/*!
  class which represents a list of objects build from forced-photometry results 
*/
typedef map<int_8, FPMesures> FPMesList;

/*!
  Class to scan a set of ForcedPhotometry output fits files produced by the LSST stack
  and create a list of sources: position, flux (mean,sigma) - single or multi colors
*/ 

class FPObjectListBuilder : public ForcedPhotFileScanBase
{
public:
  FPObjectListBuilder(string const & basedir, vector<RunFieldMinMax> const& runlist, vector<int> const& fieldlist,
		     vector<char> const& camcols, vector<char> const& filters)
    : ForcedPhotFileScanBase(basedir, runlist, fieldlist, camcols, filters), 
      totnobjs_(0), totnmes_(0), totbadmes_(0), totnmes_in_(0), nerr_oid_(0), nerr_radec_(0)
  {
  }

  FPObjectListBuilder(string const & basedir, vector<RunFieldCamCol> const& runfcclist, vector<char> const& filters)
    : ForcedPhotFileScanBase(basedir, runfcclist, filters), 
      totnobjs_(0), totnmes_(0), totbadmes_(0), totnmes_in_(0), nerr_oid_(0), nerr_radec_(0)
  {
  }

  inline size_t TotNbObjects() { return totnobjs_; }
  inline size_t TotNbMes() { return totnmes_; }
  inline size_t TotNbBadMes() { return totbadmes_; }
  inline size_t TotNbMesInRADecRegion() { return totnmes_in_; }

protected:
  //! Method called just before starting file scan 
  virtual void PrepareScan();
  //! Method called for each existing FP file scanned, with non zero file size 
  virtual int ProcessFile(const char * filename, char filter, char camcol, int field, int run);
  //! Method called at the end of the scan (to save results ...)
  virtual void FinalizeScan();

  //! method called for each object to update vmeslist_
  int UpdateList(int jf, int_8 id, int_8 oid, double ra, double dec, double flx, double refflx);

  vector< FPMesList > vmeslist_; // vector FPMesList , one for each ra,dec bin 
  map<int, int> filtindex_;
  size_t nfilt_;  // number of filters 
  size_t totnobjs_, totnmes_, totbadmes_, totnmes_in_;

  size_t nerr_oid_, nerr_radec_;
};

// }  Fin du namespace

#endif
