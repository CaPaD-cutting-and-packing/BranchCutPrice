SYSTEM     = x86-64_linux#x86_sles10_4.1#x86_RHEL3.0_3.2
#SYSTEM     = x86-64_RHEL3.0_3.2
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Enterprise_Server1261/CPLEX_Studio/cplex/
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Enterprise_Server1261/CPLEX_Studio/concert/
#### CPLEX 12.6 is not as stable:
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio126/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio126/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++
CC  = gcc
JAVAC = javac

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -g -fomit-frame-pointer -pipe -Wall -fPIC -fno-strict-aliasing -Wno-sign-compare -Wno-write-strings -Wno-reorder -Wconversion -fexceptions -DDBG_ON -DIL_STD #-DILOSTRICTPOD. -Wconversion important!

CCOOPT = -O3 -fomit-frame-pointer -pipe -Wall -fPIC -fno-strict-aliasing -Wno-sign-compare -Wno-write-strings -Wno-reorder -Wconversion -fexceptions -DNDEBUG -DIL_STD # -march=core2 -mfpmath=sse #-DILOSTRICTPOD
COPT  = #-fPIC
JOPT  = -classpath ../../../lib/cplex.jar -O

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread #-static
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -lpthread 
JAVA      = java -d64 -Djava.library.path=$(CPLEXDIR)/bin/x86-64_RHEL3.0_3.2 -classpath $(CPLEXJARDIR):


all: bcp


CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = .
EXSRC         = $(EXDIR)/cpp
EXINC         = $(EXDIR)/inc
EXDATA        = $(EXDIR)/data

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR) -I$(EXINC)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(EXINC)
CCOFLAGS = $(CCOOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(EXINC)
JCFLAGS = $(JOPT)

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log

# -----

SOURCES = bb.cpp bb_mos.cpp cutgen.cpp main.cpp probl_csp1.cpp probl_cp22.cpp probl_pmp1.cpp timer.cpp bbcuts.cpp mydebug.cpp stdafx.cpp bcp.cpp bcp_branch.cpp bcp_cuts.cpp bcp_lp.cpp bcpstuff.cpp bcp2.cpp wrpcpx.cpp mytools.cpp subgr.cpp raster.cpp spprc_dp.cpp spprc_dp1.cpp bkp_dp.cpp

OBJECTS = $(SOURCES:.cpp=DBG.o)
OOBJECTS = $(SOURCES:.cpp=O3.o)

bcpDBG: $(OBJECTS)
	$(CCC) $(CCFLAGS) -o bcpDBG $(OBJECTS) $(CCLNFLAGS)

bcp: $(OOBJECTS)
	$(CCC) $(CCOFLAGS) -o bcp $(OOBJECTS) $(CCLNFLAGS)

%DBG.o: $(EXSRC)/%.cpp
	$(CCC) -c $(CCFLAGS)  $< -o $@

%O3.o: $(EXSRC)/%.cpp
	$(CCC) -c $(CCOFLAGS)  $< -o $@

gen: cpp/cspgen.cpp
	$(CCC) $(CCFLAGS) -o gen cpp/cspgen.cpp
