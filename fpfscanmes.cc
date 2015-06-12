/*  -----------------------------------------
  class and functions to scan and analyze forced photometry output files - 
  R. Ansari - September 2013 
---------------------------------------------------- */
#include "fpfscanmes.h"

// namespace ANALSSTDC {

//--------------------------------------------------------------------------
// ---- Class FPObjectListBuilder methods 
//--------------------------------------------------------------------------

/* --Methode-- */
void FPObjectListBuilder::PrepareScan()
{
  int k=0;
  for(size_t i=0; i<filters_.size(); i++) {
    int filtid=FilterId(filters_[i]);
    if (filtid<1)  continue;
    filtindex_[filtid]=k;   k++; 
  }
  nfilt_=k;
  for(int i=0; i<radec_.nbra*radec_.nbdec; i++) {
    vmeslist_.push_back( FPMesList() );
  }
  if (prtlev_>1)  
    cout << "FPObjectListBuilder::PrepareScan() vmeslist_.size()=" << vmeslist_.size() << " NFilt=" << nfilt_ 
	 << "\n ra:" << radec_.ramin<<","<<radec_.ramax << " dec:" << radec_.decmin<<","<<radec_.decmax 
	 << " nb:" << radec_.nbra<<","<<radec_.nbdec << " -> del:" << radec_.delra<<","<<radec_.deldec << endl; 
  return;
}

/* --Methode-- */
int FPObjectListBuilder::ProcessFile(const char * filename, char filter, char camcol, int field, int run)
{
  map<int, int>::iterator it=filtindex_.find( FilterId(filter) );
  if (it==filtindex_.end())  {
    cout << " FPObjectListBuilder::ProcessFile()/ERROR filterId not found in filtindex_ map " << endl;
    return 0;
  }
  int jf = (*it).second;
  DataTable dts;
  string fitsflnm=filename;   fitsflnm += "[1]";
  fitsflnm += "[col id ; coord ; flux_psf ; refFlux ; objectId]";
  FitsInOutFile fis(fitsflnm, FitsInOutFile::Fits_RO);
  fis >> dts; 
  if (dts.NRows() < 1)  { 
  if (prtlev_>1) 
    cout << " FPObjectListBuilder::ProcessFile(" << run << "," << field << "," << camcol << "," << filter  << ") NRows=0" 
	 << " TotNbObjects=" << TotNbObjects() << endl;
    return 1;
  }
  DataTableRowPtr rsp=dts.EmptyRowPtr();
  double flx, refflx;
  double ra, dec;
  int_8 id, oid;
  for(size_t k=0; k<dts.NRows(); k++) {
    dts.GetCstRowPtr(k,rsp);
    id = rsp(0);   oid=rsp(4);
    ra=rsp(1)[0];   dec=rsp(1)[1];
    ra *= (180./M_PI);  dec *= (180./M_PI);  // convert radian to degree
    flx=rsp(2);  refflx=rsp(3);
    UpdateList(jf, id, oid, ra, dec, flx, refflx);
  }
  if (prtlev_>1) 
    cout << " FPObjectListBuilder::ProcessFile(" << run << "," << field << "," << camcol << "," << filter  << ") NRows= " 
	     << dts.NRows() << " TotNbObjects=" << TotNbObjects() << endl;
  return 0;
}

/* --Methode-- */
int FPObjectListBuilder::UpdateList(int jf, int_8 id, int_8 oid, double ra, double dec, double flx, double refflx)
{
  /*  DBG 
  if ((jf<0)||(jf>=nfilt_)) {
    cout << " !!!BUG!!! FPObjectListBuilder::UpdateList() out of range jf=" << jf << " nfilt=" << nfilt_ << endl;
    return -999;
  }
  */
  totnmes_++;
  if (!isfinite(flx)) { totbadmes_++; return -9; }
  // (kra,kdec) identifie la case en (ra,dec) dans lequel on se trouve
  int kdec=(dec-radec_.decmin)/radec_.deldec;
  int kra=(ra-radec_.ramin)/radec_.delra;
  // On verifie qu'on est dans la zone en alpha,delta selectionne 
  if ((kra<0)||(kdec<0)||(kra>=radec_.nbra)||(kdec>=radec_.nbdec))  return -1;
  totnmes_in_++;
  size_t rdidx = kdec*radec_.nbra+kra;
  FPMesList& fpml=vmeslist_[rdidx];
  FPMesList::iterator it = fpml.find(oid);
  if (it==fpml.end())  {  // adding a new source / object 
    FPMesures nsrc(id,oid,ra,dec,nfilt_);
    nsrc.AddMes(jf,flx);
    fpml[oid]=nsrc;
    //    fpml.insert( pair<int_8, FPMesures> (id,FPMesures(id,oid,ra,dec,nfilt_)) );
    totnobjs_++;
    return 1;
  }
  else {
    //  if (id!=(*it).second.id_) nerr_id_++;
    if ((fabs((*it).second.ra_-ra)>2./3600.)||(fabs((*it).second.dec_-dec)>2./3600.))  nerr_radec_++;
    (*it).second.AddMes(jf,flx);
    return 0;
  }
  return 0;
}

/* --Methode-- */
void FPObjectListBuilder::FinalizeScan()
{
  bool fg_path_only=false;
  if (outp_.length()>0) 
    if (outp_[outp_.length()-1]=='/')  fg_path_only=true;
  
  string flnm = outp_;
  if (fg_path_only)  flnm += "srclist.txt";
  cout << "ForcedPhotFileScan::FinalizeScan() : saving  object/source list to \n ->" << flnm << endl; 
  cout << " ... Total number of objects= " <<  TotNbObjects() << " NbMes=" << TotNbMes() << " NbMesIn=" << TotNbMesInRADecRegion() 
       << " BadMes=" << TotNbBadMes() << " NErrOId=" << nerr_oid_ << " NErrRADEC=" << nerr_radec_ << endl;
  ofstream of(flnm.c_str());
  of << "# id oid ra dec flxmean_1 flxsigma_1 nmes_1 flxmean_2 flxsigma_2 nmes_2 ... " << endl;
  size_t nsrc_wrt=0;
  // Boucle sur toutes les cellules en alpha, delta
  for(size_t i=0; i<vmeslist_.size(); i++) {
    FPMesList& fpml=vmeslist_[i];
    cout << " ra_dec_cell["<<i<<"] ra,dec="
	 << radec_.ramin+radec_.delra*((i%radec_.nbra)+0.5)<<","<<radec_.decmin+radec_.deldec*((i%radec_.nbdec)+0.5)
	 << " -> NbSrc=" << fpml.size() << endl; 
    // boucle sur les sources de chaque cellule
    for(FPMesList::iterator it=fpml.begin(); it!=fpml.end(); it++) {
      (*it).second.ComputeMean();
      if (((*it).second.vsflx_[0].sum_mean<minflxok_)||((*it).second.vsflx_[0].sum_mean>maxflxok_))  continue;
      (*it).second.Print(of);  // or (*it).second.Print(of," ; ");
      of << endl;   nsrc_wrt++;
    }
  }
  cout << "ForcedPhotFileScan::FinalizeScan() " << nsrc_wrt << " sources written to file out of total=" << TotNbObjects() << endl;
}


// }  Fin du namespace
