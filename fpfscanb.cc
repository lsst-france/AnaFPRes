/*  -----------------------------------------
  class and functions to scan and analyze forced photometry output files - 
  R. Ansari - September 2013 
---------------------------------------------------- */
#include "fpfscanb.h"
#include "ctimer.h"

// some specific c includes (for stat function)
#include <unistd.h>
#include <sys/stat.h>


using namespace std;
using namespace SOPHYA;
// namespace ANALSSTDC {

//--------------------------------------------------------------------------
// ---- Class ForcedPhotFileScanBase methods 
//--------------------------------------------------------------------------

/* --Methode-- */
ForcedPhotFileScanBase::ForcedPhotFileScanBase(string const & basedir, vector<RunFieldMinMax> const& runlist, 
					       vector<int> const& fieldlist, vector<char> const& camcols, 
					       vector<char> const& filters)
    : basedir_(basedir), runlist_(runlist), fieldlist_(fieldlist), camcols_(camcols), 
      fg_use_rfcc(false), filters_(filters), 
      prtlev_(0), totnfiles_(0), totmissingfiles_(0), totbadfiles_(0), totfilesize_(0)
{
  if (basedir_.length()>0) 
    if (basedir_[basedir_.length()-1]!='/')  basedir_ += '/';
  RADecLim radec;
  radec.ramin=0.;  radec.ramax=360.;  radec.nbra=36;
  radec.decmin=-90.;  radec.decmax=90.;  radec.nbra=18;
  SetRADecLimits(radec);
  SetFluxLimits();
}

/* --Methode-- */
ForcedPhotFileScanBase::ForcedPhotFileScanBase(string const & basedir, vector<RunFieldCamCol> const& runfcclist, 
					       vector<char> const& filters) 
  : basedir_(basedir), runfcclist_(runfcclist), fg_use_rfcc(true), filters_(filters), 
      prtlev_(0), totnfiles_(0), totmissingfiles_(0), totbadfiles_(0), totfilesize_(0)
{
  if (basedir_.length()>0) 
    if (basedir_[basedir_.length()-1]!='/')  basedir_ += '/';
  RADecLim radec;
  radec.ramin=0.;  radec.ramax=360.;  radec.nbra=36;
  radec.decmin=-90.;  radec.decmax=90.;  radec.nbra=18;
  SetRADecLimits(radec);
  SetFluxLimits();
}





/* --Methode-- */
int ForcedPhotFileScanBase::ScanFiles() 
{
  if (prtlev_>2) cout << "ForcedPhotFileScanBase::ScanFiles(): calling PrepareScan() ..." << endl;
  PrepareScan();
  if (fg_use_rfcc)  ScanFilesB();
  else ScanFilesA(); 
  cout << " ----- End of ForcedPhotFileScanBase::ScanFiles()  Number of files:" 
       << "\n Total= " << totnfiles_ << " Missing= " << totmissingfiles_ << " Bad= " << totbadfiles_ 
       << " NOK=" << totnfiles_-(totmissingfiles_+totbadfiles_)
       << "\n Total Size= " <<  totfilesize_/1024. << " kbytes (" << (double)totfilesize_/(1024.*1.e6) << " GB)" << endl;
  if (prtlev_>2) cout << "ForcedPhotFileScanBase::ScanFiles(): calling FinalizeScan() ..." << endl;
  FinalizeScan();
  return (totnfiles_-(totmissingfiles_+totbadfiles_));
}

/* --Methode-- */
void ForcedPhotFileScanBase::ScanFilesA() 
{
  Timer tm("ScanFiles",false);
  char filename[1024];  // complete filename
  char ccfd[8], ccfn[8];   // camcol and filter for the directory and in the file name 
  strcpy(ccfd,"1/i/");  strcpy(ccfn,"i1");
  struct stat lfst;    int rcfproc;    // to check the existence of the file
  vector<int> rfieldlist;
  for(size_t i=0; i<runlist_.size(); i++) {
    int run=runlist_[i].run;
    int fmin=runlist_[i].fieldmin; 
    int fmax=runlist_[i].fieldmax;
    vector<int> * fieldlistp=&fieldlist_;
    if ((fmin>=0)&&(fmax>=0))  {
      rfieldlist.clear();
      for(int jf=fmin; jf<=fmax; jf++)  rfieldlist.push_back(jf);
      fieldlistp=&rfieldlist;
    }
    for(size_t j=0; j<camcols_.size(); j++) {
      char camcol=camcols_[j];
      if (CamColId(camcol)==0)  {
	cout << "ForcedPhotFileScanBase::ScanFiles()/Warning - bad camcol value: " << camcol << " skipping ..." << endl;  
	continue;
      }
      ccfd[0]=ccfn[1]=camcol;
      for(size_t k=0; k<filters_.size(); k++) {
	char filter=filters_[k]; 
	if (FilterId(filter)==0)  {
	  cout << "ForcedPhotFileScanBase::ScanFiles()/Warning - bad filter value: " << filter << " skipping ..." << endl;  
	  continue;
	}
	ccfd[2]=ccfn[0]=filter; 
	for(size_t l=0; l<fieldlistp->size(); l++) {
	  int field=(*fieldlistp)[l];
	  sprintf(filename,"%s%d/%sforcedsources-%06d-%s-%04d.fits",basedir_.c_str(),run,ccfd,run,ccfn,field);
	  totnfiles_++;
	  if (prtlev_>2)  cout << totnfiles_ << "- ForcedPhotFileScan::ScanFiles() filename=" << filename << endl;
	  int rcfs=stat(filename,&lfst);
	  if (rcfs==0) {
	    totfilesize_ += lfst.st_size;  // computing total file size
	    rcfproc=ProcessFile(filename,filter,camcol,field,run);  // processing file 
	    if (rcfproc!=0)  totbadfiles_++;
	  }
	  else {
	    if (prtlev_>1)  cout << "ForcedPhotFileScanBase::ScanFiles()/Warning missing file: " << filename << endl;
	    totmissingfiles_++;
	  }
	}
      }
    }
    if (prtlev_>0) {
      tm.SplitQ();
      cout << i << "- done run=" << run << " files: Total= " << totnfiles_ << " Missing= " << totmissingfiles_ 
	   << " Bad= " << totbadfiles_ << " NOK=" << totnfiles_-(totmissingfiles_+totbadfiles_) 
	   << " CPU= " << tm.PartialCPUTime() << " Elapsed=" << tm.PartialElapsedTime() << " s" << endl;
    }
  }
  return;
} //-- fin de ScanFilesA

/* --Methode-- */
void ForcedPhotFileScanBase::ScanFilesB() 
{
  Timer tm("ScanFiles",false);
  char filename[1024];  // complete filename
  char ccfd[8], ccfn[8];   // camcol and filter for the directory and in the file name 
  strcpy(ccfd,"1/i/");  strcpy(ccfn,"i1");
  struct stat lfst;    int rcfproc;    // to check the existence of the file
  vector<int> rfieldlist;
  for(size_t i=0; i<runfcclist_.size(); i++) {
    int run=runfcclist_[i].run;
    int field=runfcclist_[i].field; 
    int ccol=runfcclist_[i].camcol;
    char camcol=CamColId2CamCol(ccol);
    ccfd[0]=ccfn[1]=camcol;
    for(size_t k=0; k<filters_.size(); k++) {
      char filter=filters_[k]; 
      if (FilterId(filter)==0)  {
	cout << "ForcedPhotFileScanBase::ScanFiles()/Warning - bad filter value: " << filter << " skipping ..." << endl;  
	continue;
      }
      ccfd[2]=ccfn[0]=filter; 
      sprintf(filename,"%s%d/%sforcedsources-%06d-%s-%04d.fits",basedir_.c_str(),run,ccfd,run,ccfn,field);
      totnfiles_++;
      if (prtlev_>2)  cout << totnfiles_ << "- ForcedPhotFileScan::ScanFiles() filename=" << filename << endl;
      int rcfs=stat(filename,&lfst);
      if (rcfs==0) {
	totfilesize_ += lfst.st_size;  // computing total file size
	rcfproc=ProcessFile(filename,filter,camcol,field,run);  // processing file 
	if (rcfproc!=0)  totbadfiles_++;
      }
      else {
	if (prtlev_>1)  cout << "ForcedPhotFileScanBase::ScanFiles()/Warning missing file: " << filename << endl;
	totmissingfiles_++;
      }
    }
    if (prtlev_>0) {
      tm.SplitQ();
      cout << i << "- done run=" << run << " files: Total= " << totnfiles_ << " Missing= " << totmissingfiles_ 
	   << " Bad= " << totbadfiles_ << " NOK=" << totnfiles_-(totmissingfiles_+totbadfiles_) 
	   << " CPU= " << tm.PartialCPUTime() << " Elapsed=" << tm.PartialElapsedTime() << " s" << endl;
    }
  }  
  return;
} //-- fin de ScanFilesB

//--------------------------------------------------------------------------
// ---- Class ForcedPhotFileScan methods 
//--------------------------------------------------------------------------

/* --Methode-- */
void ForcedPhotFileScan::PrepareScan()
{
  if (!fg_fill_dtrf_)  return;
  // summary info data table
  dtrunfield_.AddIntegerColumn("run");
  dtrunfield_.AddIntegerColumn("field");
  dtrunfield_.AddIntegerColumn("camcol_filter");
  dtrunfield_.AddIntegerColumn("nsrc");
  // mnx -> min and max (columns with a vector of 2 values, min and max) 
  dtrunfield_.AddDoubleColumn("ra_mnx",2);    
  dtrunfield_.AddDoubleColumn("dec_mnx",2);
  dtrunfield_.AddLongColumn("id_mnx",2);
  dtrunfield_.AddLongColumn("oid_mnx",2);
  dtrunfield_.AddDoubleColumn("flux_mnx",2);
  dtrunfield_.AddIntegerColumn("nfluxok");
  dtrunfield_.AddDoubleColumn("fluxmean");
  return;
}


/* --Methode-- */
int ForcedPhotFileScan::ProcessFile(const char * filename, char filter, char camcol, int field, int run) 
{
  if (prtlev_>3) 
    cout << " ForcedPhotFileScan::ProcessFile( " << filename << " )-> filterId=" << FilterId(filter) 
	 << " CamColId= " << CamColId(camcol) << endl; 
  map<int, RunFieldMinMax>::iterator it = rfm_.find(run);
  // On met a jour la liste des runs avec limite des numeros de fields
  if (it == rfm_.end())  {
    RunFieldMinMax rfv; 
    rfv.run=run; rfv.fieldmin=rfv.fieldmax=field;
    rfm_[run]=rfv;
  }
  else {
    if ((*it).second.fieldmin > field)  (*it).second.fieldmin=field;
    if ((*it).second.fieldmax < field)  (*it).second.fieldmax=field;
  }
  if (!fg_fill_dtrf_)  {
    string fitsflnm=filename;   fitsflnm += "[1]";
    FitsInOutFile fis(fitsflnm, FitsInOutFile::Fits_RO);
    int_8 nrows=fis.GetNbRows();
    if (prtlev_>2) 
      cout << filename << " -> NRows= " << nrows << endl;
    if (nrows<1) {
      if (prtlev_>1) cout << "ForcedPhotFileScan::ProcessFile()/Warning NRows=0 for " << filename << endl;
      return 1;
    }
    return 0;
  }
  // on ouvre le fichier fits et on remplit le data-table dtrunfield_
  {
    DataTable dts;
    string fitsflnm=filename;   fitsflnm += "[1]";
    fitsflnm += "[col id ; coord ; flux_psf ; refFlux ; objectId]";
    FitsInOutFile fis(fitsflnm, FitsInOutFile::Fits_RO);
    fis >> dts; 
    if (prtlev_>1)  cout << filename << " : NSrc= " << dts.NRows() << endl;
    DataTableRowPtr rowp=dtrunfield_.EmptyRowPtr();
    dtrunfield_.NextRowPtr(rowp);
    rowp(0)=run;
    rowp(1)=field;
    rowp(2)=10*CamColId(camcol)+FilterId(filter); 
    rowp(3)=(int_4)dts.NRows();  // nombre de source

    if (dts.NRows() > 0)  {  // on boucle sur les sources et on calcule quelques quantites 
      DataTableRowPtr rsp=dts.EmptyRowPtr();
      double flx, refflx, flxmin=9.e19, flxmax=-9.e19, flxmoy=0.; 
      int_4 nflxok=0;
      double ra, ramin=9.e9, ramax=-9.e9;
      double dec, decmin=9.e9, decmax=-9.e9;
      int_8 id, idmin=-1, idmax=-2;
      int_8 oid, oidmin=-1, oidmax=-2;
      for(size_t k=0; k<dts.NRows(); k++) {
	dts.GetCstRowPtr(k,rsp);
	id = rsp(0);   oid=rsp(4);
	ra=rsp(1)[0];   dec=rsp(1)[1];
	ra *= (180./M_PI);  dec *= (180./M_PI);  // convert radian to degree
	flx=rsp(2);  refflx=rsp(3);
	if (id<idmin) idmin=id;   if (id>idmax) idmax=id;
	if (oid<oidmin) oidmin=oid;   if (oid>oidmax) oidmax=oid;
	if (ra<ramin) ramin=ra;   if (ra>ramax) ramax=ra;
	if (dec<decmin) decmin=dec;   if (dec>decmax) decmax=dec;
	if (flx<flxmin) flxmin=flx;   if (flx>flxmax) flxmax=flx;
	if ((flx>minflxok_)&&(flx<maxflxok_))  {nflxok++; flxmoy+=flx; }
      }
      if (nflxok>0)  flxmoy /= (double)nflxok;
      rowp(4)[0]=ramin;    rowp(4)[1]=ramax; 
      rowp(5)[0]=decmin;   rowp(5)[1]=decmax; 
      rowp(6)[0]=idmin;    rowp(6)[1]=idmax; 
      rowp(7)[0]=oidmin;   rowp(7)[1]=oidmax; 
      rowp(8)[0]=flxmin;   rowp(8)[1]=flxmax; 
      rowp(9)=nflxok;      rowp(10)=flxmoy;
      return 0;
    }
    else {
      rowp(4)[0]=0.;    rowp(4)[1]=0.; 
      rowp(5)[0]=0.;    rowp(5)[1]=0.; 
      rowp(6)[0]=0;     rowp(6)[1]=0; 
      rowp(7)[0]=0;     rowp(7)[1]=0; 
      rowp(8)[0]=0;     rowp(8)[1]=0; 
      rowp(9)=0;        rowp(10)=0.;
      return 1;
    }
  }
  return 0;
}

/* --Methode-- */
void ForcedPhotFileScan::FinalizeScan() 
{
  bool fg_path_only=false;
  if (outp_.length()>0) 
    if (outp_[outp_.length()-1]=='/')  fg_path_only=true;
  {
    string flnm = outp_;
    if (fg_path_only)  flnm += "runfield_list.txt";
    else flnm += ".txt";
    cout << "ForcedPhotFileScan::FinalizeScan() : saving run  fieldmin fieldmax values to \n ->" << flnm << endl; 
    ofstream of(flnm.c_str());
    map<int, RunFieldMinMax>::iterator it;
    for(it=rfm_.begin(); it!=rfm_.end(); it++) 
      of << (*it).second.run << " " << (*it).second.fieldmin << " " << (*it).second.fieldmax << endl;
  }
  if (fg_fill_dtrf_)  {
    string flnm = "!" + outp_;    // the ! character ensures file overwrite
    if (fg_path_only)  flnm += "fpfsumm.fits"; 
    else flnm += ".fits";
    cout << "ForcedPhotFileScan::FinalizeScan(): saving forcedphotometry files summary information to \n ->" << flnm << endl; 
    FitsInOutFile fos(flnm, FitsInOutFile::Fits_Create);
    fos << dtrunfield_;
    dtrunfield_.SetShowMinMaxFlag(true);
    cout << dtrunfield_;
  }

}

// }  Fin du namespace
