#  Example makefile to compile with SOPHYA and GSL 
include $(SOPHYABASE)/include/sophyamake.inc

#  Define our target list 
all : fpfscan mergeft adselectft

clean :
	rm -f Objs/*  

# list of include files 
MYINCLIST = analsstdc.h fpfscanb.h fpfscanmes.h
# list of obj files used by main program(s) 
MYOBJLIST = Objs/analsstdc.o Objs/fpfscanb.o Objs/fpfscanmes.o 

######
###  Forced photometry files scan program 
fpfscan : Objs/fpfscan.o $(MYOBJLIST)
	$(CXXLINK) -o Objs/fpfscan Objs/fpfscan.o $(MYOBJLIST) $(SOPHYAEXTSLBLIST) 

Objs/fpfscan.o : fpfscan.cc $(MYINCLIST)
	$(CXXCOMPILE) -o Objs/fpfscan.o  fpfscan.cc

Objs/analsstdc.o : analsstdc.cc  $(MYINCLIST)
	$(CXXCOMPILE) -o Objs/analsstdc.o  analsstdc.cc

Objs/fpfscanb.o : fpfscanb.cc  $(MYINCLIST)
	$(CXXCOMPILE) -o Objs/fpfscanb.o  fpfscanb.cc

Objs/fpfscanmes.o : fpfscanmes.cc  $(MYINCLIST)
	$(CXXCOMPILE) -o Objs/fpfscanmes.o  fpfscanmes.cc

######  stand alone programs (not depending on MYINCLIST  and MYOBJLIST (might use SOPHYA)
###  Merge FITS tables  
mergeft : Objs/mergeft.o 
	$(CXXLINK) -o Objs/mergeft Objs/mergeft.o  $(SOPHYAEXTSLBLIST) 

Objs/mergeft.o : mergeft.cc 
	$(CXXCOMPILE) -o Objs/mergeft.o  mergeft.cc

### alpha,delta select from FITS table
adselectft : Objs/adselectft.o 
	$(CXXLINK) -o Objs/adselectft Objs/adselectft.o  $(SOPHYAEXTSLBLIST) 

Objs/adselectft.o : adselectft.cc 
	$(CXXCOMPILE) -o Objs/adselectft.o  adselectft.cc

