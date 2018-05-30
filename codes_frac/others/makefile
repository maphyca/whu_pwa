# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ROOTSYS = ./
include $(ROOTSYS)/Makefile.arch

#------------------------------------------------------------------------------

ifeq ($(PLATFORM),win32)
LIBS        += "$(ROOTSYS)/lib/libMinuit2.lib"
else
LIBS        += -lMinuit2 
CXXFLAGS += -g      
endif

# for using with MPI
ifneq ($(USE_OPENMPI),)
CXX=mpic++
LD=mpic++
endif

PWASRC = test_minuit2.$(SrcSuf) DPFPWAPoint.$(SrcSuf) DPFAngular.$(SrcSuf) fit_parameter_interface.$(SrcSuf) whu_my_parameter.$(SrcSuf) data_obj.$(SrcSuf) whu_pwa_fcn.$(SrcSuf) whu_propogator.$(SrcSuf) calculate0p.$(SrcSuf)
PWAOBJ = test_minuit2.$(ObjSuf) DPFPWAPoint.$(ObjSuf) DPFAngular.$(ObjSuf) fit_parameter_interface.$(ObjSuf) whu_my_parameter.$(ObjSuf) data_obj.$(ObjSuf) whu_pwa_fcn.$(ObjSuf) whu_propogator.$(ObjSuf) calculate0p.$(ObjSuf)
PWA = test

OBJS          = $(PWAOBJ)
PROGRAMS      = $(PWA)
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) $(ExeSuf)


all:            $(PROGRAMS)

$(DEMOMINIMIZER): 	$(DEMOMINIMIZEROBJ) 
			$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
			@echo "$@ done"


$(PWAOBJ) : DPFPWAPoint.h DPFAngular.h fit_parameter_interface.h whu_constants_and_definitions.h whu_my_parameter.h whu_pwa_fcn.h DPFPropogator.h
$(PWA): $(PWAOBJ)
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
	@echo "$@ done"


clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS)  *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations

.SUFFIXES: .$(SrcSuf)


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
