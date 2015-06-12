/*  -----------------------------------------
  class and functions to scan and analyze forced photometry output files - 
  R. Ansari - September 2013 
---------------------------------------------------- */
#include "analsstdc.h"

// some specific c includes (for stat function)
#include <unistd.h>
#include <sys/stat.h>

#include "strutilxx.h"   

using namespace std;
using namespace SOPHYA;
// namespace ANALSSTDC {
// -------------------------------------------------------
// ---- 
// -------------------------------------------------------


/* -- Fonction -- */
void Ana_LSST_DC_Usage()
{
  cout << " fpfscan (forced photometry output file scan/check) usage: \n"
       << " fpfscan ACT [-ra min,max,nbin] [-dec min,max,nbin] -flux min,max] \n"
       << "     [-outp path] [-basedir path] [-prtlev level] \n" 
       << "     -camcols c1,c2... -filters i,g...   \n"
       << "     -runlist r1,r2,r3...] [-fieldlist f1,f2,f3... \n"
       << "     [-runs first,last] [-fields first,last] \n"
       << "     [-runlistfile filename] [-fieldlistfile filename] \n"
       << " ACT = -h/help, LIST, SUMMFT , SRCLIST  ... \n"
       << "     list : list bad files and create [ run fiedlmin fieldmax ] list file \n"
       << "     summft : analyze each file and create a summary table (fits format) \n"
       << "     srclist : analyze each file and create a multi-filter source list \n"
       << " -ra  -dec : Right ascension, declination range definition and number of bins \n"
       << "     number of bins apply to sky maps,histos... \n"
       << " -flux : flux min-max cut \n" 
       << " -basedir : input fits files base directory path \n"
       << " -outp : output file name (without extension) or \n" 
       << "         output directory path (the last character should be /) \n"
       << " -prtlev : print level (0 1 2 ...) \n"
       << " -camcols -filters : camcol (1,2,3,4,5,6) and filter (u,g,r,i,z) list definition \n"
       << " -runlist -fieldlist : explicit list of runs / fields \n"
       << " -runs -fields : run and field list definition first<=run,field<=last \n"
       << " -runlistfile  : explicit list of [runs field_min field_max] from a file \n"
       << " -rfcclistfile  : explicit list of [runs field camcol] from a file \n"
       << "    file format: text file with one number per line with three numbers per line:\n"
       << "    [ run_number min_field_number max_field_number ]  OR  [ run field camcol ] \n"
       << "  NOTE: camcol, filter, runlist and fieldlist definition is MANDATORY \n" << endl;
  return;
}

int Ana_LSST_DC_read_listfile(const char* filename, vector<int>& list);

/* -- Fonction -- */
 int Ana_LSST_DC_DecodeArgs(int narg, char* arg[], string& action, string& basedir, 
			    vector<RunFieldMinMax>& runlist, vector<RunFieldCamCol>& runfcclist, 
			    vector<int>& fieldlist, vector<char>& camcols, vector<char>& filters,  
			    RADecLim& radec, string& outpath, int& prtlev, pair<double,double> & fluxcut)
{
  action="";
  runlist.clear(); 
  runfcclist.clear();  
  camcols.clear();
  filters.clear();
  fieldlist.clear();
  outpath=basedir="";
  radec.ramin=5.;  radec.ramax=55.;  radec.nbra=0;
  radec.decmin=-1.;  radec.decmax=1.;  radec.nbra=0;
  fluxcut.first=0.; fluxcut.second=-1.;

  prtlev=0;
  
  if ((narg>1)&&((strcmp(arg[1],"-h")==0)||(strcmp(arg[2],"help")==0)) )  { 
    Ana_LSST_DC_Usage();
    return 1;
  }
  if (narg<6) { 
    Ana_LSST_DC_Usage();
    return 2; 
  }
  action=arg[1];
  int ka=2;
  while (ka<(narg-1)) {
    //DBG cout << "*DBG*DecodeArgs() - ka="<<ka<<" -> "<<arg[ka]<<" "<<arg[ka+1]<<endl;
    if (strcmp(arg[ka],"-ra")==0) {
      sscanf(arg[ka+1],"%lg,%lg,%d",&radec.ramin,&radec.ramax,&radec.nbra);  
      ka+=2;
    }
    else if (strcmp(arg[ka],"-dec")==0) {
      sscanf(arg[ka+1],"%lg,%lg,%d",&radec.decmin,&radec.decmax,&radec.nbdec);  
      ka+=2;
    }
    else if (strcmp(arg[ka],"-flux")==0) {
      sscanf(arg[ka+1],"%lg,%lg",&(fluxcut.first),&(fluxcut.second));
      ka+=2;
    }
    else if (strcmp(arg[ka],"-basedir")==0) {
      basedir=arg[ka+1];
      ka+=2;
    }
    else if (strcmp(arg[ka],"-outp")==0) {
      outpath=arg[ka+1];
      ka+=2;
    }
    else if (strcmp(arg[ka],"-prtlev")==0) {
      prtlev=atoi(arg[ka+1]);
      ka+=2;
    }
    else if ((strcmp(arg[ka],"-camcols")==0)||(strcmp(arg[ka],"-filters")==0)) {
      string s = arg[ka+1];
      vector<string> vs;
      SplitStringToVString(s,vs,',');
      vector<char>* vc=&camcols;
      if (strcmp(arg[ka],"-filters")==0) vc=&filters;
      vc->clear();
      for(size_t i=0; i<vs.size(); i++) vc->push_back(vs[i][0]);
      ka+=2;
    }
    else if ((strcmp(arg[ka],"-runlist")==0)||(strcmp(arg[ka],"-fieldlist")==0)) {
      string s = arg[ka+1];
      vector<string> vs;
      SplitStringToVString(s,vs,',');
      vector<int> vi;
      for(size_t i=0; i<vs.size(); i++) vi.push_back(atoi(vs[i].c_str()));
      if (strcmp(arg[ka],"-fieldlist")==0) fieldlist=vi;
      else {
	runlist.clear();
	RunFieldMinMax runf; runf.run=runf.fieldmin=runf.fieldmax=-1;
	for(size_t i=0; i<vi.size(); i++)  {
	  runf.run=vi[i];  runlist.push_back(runf);
	}
      }
      ka+=2;
    }
    else if ((strcmp(arg[ka],"-runs")==0)||(strcmp(arg[ka],"-fields")==0)) {
      int first=0,last=0;
      sscanf(arg[ka+1],"%d,%d",&first,&last); 
      if (strcmp(arg[ka],"-fields")==0) {
	fieldlist.clear();
	for(int i=first; i<=last; i++)   fieldlist.push_back(i);
      }
      else {
	runlist.clear();
	RunFieldMinMax runf; runf.run=runf.fieldmin=runf.fieldmax=-1;
	for(int i=first; i<=last; i++)   {
	  runf.run=i;  runlist.push_back(runf);
	}
      }
      ka+=2;
    }
    else if (strcmp(arg[ka],"-runlistfile")==0) {
      Ana_LSST_DC_read_listfile(arg[ka+1], runlist);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-rfcclistfile")==0) {
      Ana_LSST_DC_read_listfile2(arg[ka+1], runfcclist);
      ka+=2;
    }
    else {
      cout << " fpfscan/DecodeArgs() Warning: bad argument: " << arg[ka] << endl;
      ka++;
    }
  }

  if (runfcclist.size()>0) {
    if (filters.size()<1)  { 
      cout << " Ana_LSST_DC_DecodeArgs/Bad args : filters.size()=0 with runfcclist " << endl;
      return 9;
    }
    else cout << " Ana_LSST_DC_DecodeArgs/Info: using [run field camcol] triplets " << endl;
  }
  else {
    if ((runlist.size()<1)||(fieldlist.size()<1)||(camcols.size()<1)||(filters.size()<1)) {
      cout << "  Ana_LSST_DC_DecodeArgs/Bad args - runlist.size()="<<runlist.size()<<" fieldlist.size()=" 
	   << fieldlist.size() << " \n camcols.size()="<<camcols.size()<<" filters.size()="<<filters.size()<<endl;
    return 9;
    }
    else cout << " Ana_LSST_DC_DecodeArgs/Info: using runlist fieldlist camcols filters " << endl;
  }
  return 0;
}

/* -- Fonction -- */
int Ana_LSST_DC_read_listfile(const char* filename, vector<RunFieldMinMax>& list)
// retourne le nombre de lignes lues
{
  list.clear();
  char line[128];  size_t kl=0;
  ifstream ifs(filename);  // std::ifstream::in);
  RunFieldMinMax  runf;
  while (ifs.is_open() && (!ifs.eof())) {
    ifs.getline(line,128);  line[127]='\0';
    if (ifs.eof())  break;
    runf.run=0; runf.fieldmin=runf.fieldmax=-1;
    sscanf(line,"%d %d %d",&runf.run, &runf.fieldmin, &runf.fieldmax);
    list.push_back(runf);
    // cout << " *DBG*Line= " << line << endl;
    kl++;
  }
  cout << "Ana_LSST_DC_read_listfile("<<filename<<") "<<kl<<" lines read from file"<<endl;
  return kl;
}

/* -- Fonction -- */
int Ana_LSST_DC_read_listfile2(const char* filename, vector<RunFieldCamCol>& list)
// retourne le nombre de lignes lues
{
  list.clear();
  char line[128];  size_t kl=0;
  ifstream ifs(filename);  // std::ifstream::in);
  RunFieldCamCol  runf;
  while (ifs.is_open() && (!ifs.eof())) {
    ifs.getline(line,128);  line[127]='\0';
    if (ifs.eof())  break;
    runf.run=0; runf.field=-1; runf.camcol=0;
    sscanf(line,"%d %d %d",&runf.run, &runf.field, &runf.camcol);
    list.push_back(runf);
    // cout << " *DBG*Line= " << line << endl;
    kl++;
  }
  cout << "Ana_LSST_DC_read_listfile2("<<filename<<") "<<kl<<" lines read from file"<<endl;
  return kl;
}

// }  Fin du namespace
