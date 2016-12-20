#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/cpp/BrownianMotion.o \
	${OBJECTDIR}/cpp/Init.o \
	${OBJECTDIR}/cpp/MonteCarlo.o \
	${OBJECTDIR}/cpp/Ode.o \
	${OBJECTDIR}/cpp/OdeProblems.o \
	${OBJECTDIR}/cpp/OdeRungeKuttaIntegrator.o \
	${OBJECTDIR}/cpp/OdeStabilizedIntegrators.o \
	${OBJECTDIR}/cpp/Sde.o \
	${OBJECTDIR}/cpp/SdeProblems.o \
	${OBJECTDIR}/cpp/SdeRungeKuttaIntegrator.o \
	${OBJECTDIR}/cpp/SdeStabilizedIntegrators.o \
	${OBJECTDIR}/cpp/TraditionalSdeRungeKuttaIntegrators.o \
	${OBJECTDIR}/cpp/main.o \
	${OBJECTDIR}/cpp/zufall.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-std=c++11 -fopenmp
CXXFLAGS=-std=c++11 -fopenmp

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/stochastic_adaptive_time_stepping

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/stochastic_adaptive_time_stepping: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/stochastic_adaptive_time_stepping ${OBJECTFILES} ${LDLIBSOPTIONS} -std=c++11 -fopenmp

${OBJECTDIR}/cpp/BrownianMotion.o: cpp/BrownianMotion.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/BrownianMotion.o cpp/BrownianMotion.cpp

${OBJECTDIR}/cpp/Init.o: cpp/Init.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/Init.o cpp/Init.cpp

${OBJECTDIR}/cpp/MonteCarlo.o: cpp/MonteCarlo.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/MonteCarlo.o cpp/MonteCarlo.cpp

${OBJECTDIR}/cpp/Ode.o: cpp/Ode.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/Ode.o cpp/Ode.cpp

${OBJECTDIR}/cpp/OdeProblems.o: cpp/OdeProblems.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/OdeProblems.o cpp/OdeProblems.cpp

${OBJECTDIR}/cpp/OdeRungeKuttaIntegrator.o: cpp/OdeRungeKuttaIntegrator.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/OdeRungeKuttaIntegrator.o cpp/OdeRungeKuttaIntegrator.cpp

${OBJECTDIR}/cpp/OdeStabilizedIntegrators.o: cpp/OdeStabilizedIntegrators.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/OdeStabilizedIntegrators.o cpp/OdeStabilizedIntegrators.cpp

${OBJECTDIR}/cpp/Sde.o: cpp/Sde.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/Sde.o cpp/Sde.cpp

${OBJECTDIR}/cpp/SdeProblems.o: cpp/SdeProblems.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/SdeProblems.o cpp/SdeProblems.cpp

${OBJECTDIR}/cpp/SdeRungeKuttaIntegrator.o: cpp/SdeRungeKuttaIntegrator.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/SdeRungeKuttaIntegrator.o cpp/SdeRungeKuttaIntegrator.cpp

${OBJECTDIR}/cpp/SdeStabilizedIntegrators.o: cpp/SdeStabilizedIntegrators.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/SdeStabilizedIntegrators.o cpp/SdeStabilizedIntegrators.cpp

${OBJECTDIR}/cpp/TraditionalSdeRungeKuttaIntegrators.o: cpp/TraditionalSdeRungeKuttaIntegrators.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/TraditionalSdeRungeKuttaIntegrators.o cpp/TraditionalSdeRungeKuttaIntegrators.cpp

${OBJECTDIR}/cpp/main.o: cpp/main.cpp 
	${MKDIR} -p ${OBJECTDIR}/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ihpp -I/u/anmc/rosilho/Programs/eigen-3.3.1 -I/u/anmc/rosilho/Programs/getpot-c++ -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cpp/main.o cpp/main.cpp

${OBJECTDIR}/cpp/zufall.o: cpp/zufall.f 
	${MKDIR} -p ${OBJECTDIR}/cpp
	$(COMPILE.f) -g -o ${OBJECTDIR}/cpp/zufall.o cpp/zufall.f

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/stochastic_adaptive_time_stepping
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
