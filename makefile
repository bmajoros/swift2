CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas

$(OBJ):
	mkdir $(OBJ)

all: \
	permissions \
	BOOM/libBOOM.a \
	swift2

BOOM/libBOOM.a:
	cd BOOM ; make all ; cd -

permissions:
	chmod a+x *.py

#---------------------------------------------------------
$(OBJ)/swift2.o:\
		swift2.C
	$(CC) $(CFLAGS) -o $(OBJ)/swift2.o -c \
		swift2.C
#---------------------------------------------------------
swift2: \
		$(OBJ) \
		$(OBJ)/InverseCDF.o \
		$(OBJ)/GridMap.o \
		$(OBJ)/PosteriorEstimator.o \
		$(OBJ)/DensityGrid.o \
		$(OBJ)/Trapezoids.o \
		$(OBJ)/PrefixSumArray.o \
		$(OBJ)/DensityWithPrior.o \
		$(OBJ)/DensityFunction.o \
		$(OBJ)/EmpiricalPvalues.o \
		$(OBJ)/Experiment.o \
		$(OBJ)/Swift.o \
		$(OBJ)/Simulator.o \
		$(OBJ)/Replicate.o \
		$(OBJ)/Replicates.o \
		$(OBJ)/SwiftSample.o \
		$(OBJ)/swift2.o
	$(CC) $(LDFLAGS) -o swift2 \
		$(OBJ)/InverseCDF.o \
		$(OBJ)/GridMap.o \
		$(OBJ)/PosteriorEstimator.o \
		$(OBJ)/DensityGrid.o \
		$(OBJ)/Trapezoids.o \
		$(OBJ)/PrefixSumArray.o \
		$(OBJ)/DensityWithPrior.o \
		$(OBJ)/DensityFunction.o \
		$(OBJ)/EmpiricalPvalues.o \
		$(OBJ)/Experiment.o \
		$(OBJ)/Swift.o \
		$(OBJ)/Simulator.o \
		$(OBJ)/Replicate.o \
		$(OBJ)/Replicates.o \
		$(OBJ)/SwiftSample.o \
		$(OBJ)/swift2.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/Replicate.o:\
		Replicate.C\
		Replicate.H
	$(CC) $(CFLAGS) -o $(OBJ)/Replicate.o -c \
		Replicate.C
#--------------------------------------------------------
$(OBJ)/Replicates.o:\
		Replicates.C\
		Replicates.H
	$(CC) $(CFLAGS) -o $(OBJ)/Replicates.o -c \
		Replicates.C
#--------------------------------------------------------
$(OBJ)/SwiftSample.o:\
		SwiftSample.C\
		SwiftSample.H
	$(CC) $(CFLAGS) -o $(OBJ)/SwiftSample.o -c \
		SwiftSample.C
#---------------------------------------------------------
$(OBJ)/Swift.o:\
		Swift.C\
		Swift.H
	$(CC) $(CFLAGS) -o $(OBJ)/Swift.o -c \
		Swift.C
#---------------------------------------------------------
$(OBJ)/Simulator.o:\
		Simulator.C\
		Simulator.H
	$(CC) $(CFLAGS) -o $(OBJ)/Simulator.o -c \
		Simulator.C
#--------------------------------------------------------
$(OBJ)/Experiment.o:\
		Experiment.C\
		Experiment.H
	$(CC) $(CFLAGS) -o $(OBJ)/Experiment.o -c \
		Experiment.C
#--------------------------------------------------------
$(OBJ)/EmpiricalPvalues.o:\
		EmpiricalPvalues.C\
		EmpiricalPvalues.H
	$(CC) $(CFLAGS) -o $(OBJ)/EmpiricalPvalues.o -c \
		EmpiricalPvalues.C
#---------------------------------------------------------
$(OBJ)/DensityFunction.o:\
		DensityFunction.C\
		DensityFunction.H
	$(CC) $(CFLAGS) -o $(OBJ)/DensityFunction.o -c \
		DensityFunction.C
#---------------------------------------------------------
$(OBJ)/DensityGrid.o:\
		DensityGrid.C\
		DensityGrid.H
	$(CC) $(CFLAGS) -o $(OBJ)/DensityGrid.o -c \
		DensityGrid.C
#--------------------------------------------------------
$(OBJ)/Trapezoids.o:\
		Trapezoids.C\
		Trapezoids.H
	$(CC) $(CFLAGS) -o $(OBJ)/Trapezoids.o -c \
		Trapezoids.C
#--------------------------------------------------------
$(OBJ)/PrefixSumArray.o:\
		PrefixSumArray.C\
		PrefixSumArray.H
	$(CC) $(CFLAGS) -o $(OBJ)/PrefixSumArray.o -c \
		PrefixSumArray.C
#--------------------------------------------------------
$(OBJ)/PosteriorEstimator.o:\
		PosteriorEstimator.C\
		PosteriorEstimator.H
	$(CC) $(CFLAGS) -o $(OBJ)/PosteriorEstimator.o -c \
		PosteriorEstimator.C
#---------------------------------------------------------
$(OBJ)/GridMap.o:\
		GridMap.C\
		GridMap.H
	$(CC) $(CFLAGS) -o $(OBJ)/GridMap.o -c \
		GridMap.C
#---------------------------------------------------------
$(OBJ)/InverseCDF.o:\
		InverseCDF.C\
		InverseCDF.H
	$(CC) $(CFLAGS) -o $(OBJ)/InverseCDF.o -c \
		InverseCDF.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/DensityWithPrior.o:\
		DensityWithPrior.C\
		DensityWithPrior.H
	$(CC) $(CFLAGS) -o $(OBJ)/DensityWithPrior.o -c \
		DensityWithPrior.C
#---------------------------------------------------------
