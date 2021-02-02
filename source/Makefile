#Makefile for parallel compiling  
#FC=h5pfc -r8 -O3 -fpp #-check bounds -traceback
#FC += -ip # -xAVX -axCORE-AVX2
#FC += -I$(HOME)/local/fftw-3.3.7/include

#-- Marenostrum
FC = h5pfc -r8 -O3 -fpp

# flow config
#FC+=-DMOVIE
FC+=-DPLANEMOVIE
FC+=-DSTRUCMOVIE

# Library
#LINKS=-lhdf5_fortran -lhdf5
#LINKS+=-L$(HOME)/local/fftw-3.3.7/lib -lfftw3 -lfftw3_omp
##LINKS+=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#LINKS+=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

#-- Marenostrum
LINKS=-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lsz -lz
LINKS+=-L${FFTW_DIR}/lib/ -lfftw3
LINKS+=-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

#============================================================================ 
# make PROGRAM   
#============================================================================
OBJECTS_RAW=inimov_hdf_ycut.o mkmov_hdf_ycut.o inimov_hdf_xcut.o mkmov_hdf_xcut.o inimov_hdf_zcut.o mkmov_hdf_zcut.o\
            initParams.o cfl.o cflr.o coetar.o cordin.o continua.o divg.o divgck.o \
            addbodyobj.o  wallBCs.o  \
            HdfAuxRoutines.o \
            divgloc.o gcurv.o hdnlq1.o hdnlq2.o hdnlq3.o hdnlsa.o hdnlte.o hdnlco2.o \
            indic.o inimov.o inirea.o initia.o inqpr.o invtrq1.o invtrq2.o \
            invtrq3.o invtrsa.o invtrte.o invtrco2.o matrix_transpose.o meshes.o \
            mpi_routines.o openfi.o param.o papero.o phcalc.o phini.o \
            prcalc.o solq1k.o solq2k.o solq3k.o solsak.o solxi_endwall.o solxj_endwall.o solxi_fs.o solxj_fs.o solxi_ns.o solxj_ns.o solxri.o solxrj.o solxi.o solxj.o soltek.o solco2k.o solxj_scafixval.o solxrj_scafixval.o  \
            solxj_ns_q1out.o solxj_ns_q3out.o \
            solxj_scafixval_co2inout.o solxj_scafixval_densinout.o solxrj_scafixval_dsalinout.o \
            tsch.o updvp.o vmaxv.o nusse_walls.o \
            coetarr.o indicr.o write_field.o  mkmov_dsal.o velbc.o  \
            mgrd_idc.o mgrd_mem.o mgrd_velitp.o mgrd_dsalc.o \
			readobj1.o readobj2.o obj_io.o geometry.o volavg.o\
			inoutflowBCs.o

PROGRAM=simexec
MODULES=obj/param.o
MODLOCS=mod/
OBJLOCS=obj/
OBJECTS=$(addprefix $(OBJLOCS),$(OBJECTS_RAW))
#============================================================================ 
# Linking    
#============================================================================
$(PROGRAM):objlocs modlocs $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(LINKS) -module $(MODLOCS) -o $@

#============================================================================
#  Dependencies
#============================================================================
obj/param.o: ./param.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/tsch.o: ./tsch.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/gcurv.o: ./gcurv.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/papero.o: ./papero.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/%.o: particle_ibm/%.f90 $(MODULES)
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/%.o: fluid_solver/%.f90 $(MODULES)
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/%.o: output/%.f90 $(MODULES)
	$(FC) -c -o $@ -module $(MODLOCS) $<
#============================================================================
#  Clean up
#============================================================================
clean :
	rm obj/*.o 
	rm mod/*.mod
	rm $(PROGRAM)

#============================================================================
#  Definitions
#============================================================================
.PHONY: objlocs
objlocs: ${OBJLOCS}
${OBJLOCS}:
	mkdir -p ${OBJLOCS}
.PHONY: modlocs
modlocs: ${MODLOCS}
${MODLOCS}:
	mkdir -p ${MODLOCS}
