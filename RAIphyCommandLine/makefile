#***************************************************************************/
#**    Author: Sam Way                                                     */
#**    Date: 02/22/2011                                                    */
#**    Description: Makefile for RAIphy.                                   */
#***************************************************************************/
CXX      = g++
CXXFLAGS = -O3 -funroll-loops -Wall
EXEC     = raiphy 
OBJS = \
    Main.o \
    Parms.o \
    RAIphyDatabase.o \
    RAIProfiler.o \
    FragmentClassifier.o
SOURCE = \
    Main.cpp \
    Parms.cpp \
    RAIphyDatabase.cpp \
    RAIProfiler.cpp \
    FragmentClassifier.cpp

${EXEC}: ${OBJS} ${SOURCE}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${EXEC} 

Main.o: Main.cpp Globals.h Parms.h RAIphyDatabase.h
	${CXX} -c ${CXXFLAGS} Main.cpp

Parms.o: Parms.cpp Globals.h Parms.h
	${CXX} -c ${CXXFLAGS} Parms.cpp

RAIphyDatabase.o: RAIphyDatabase.cpp Globals.h Parms.h RAIProfiler.h RAIphyDatabase.h
	${CXX} -c ${CXXFLAGS} RAIphyDatabase.cpp

RAIProfiler.o: RAIProfiler.cpp Globals.h Parms.h RAIProfiler.h 
	${CXX} -c ${CXXFLAGS} RAIProfiler.cpp

FragmentClassifier.o: FragmentClassifier.cpp Globals.h Parms.h RAIphyDatabase.h RAIProfiler.h FragmentClassifier.h
	${CXX} -c ${CXXFLAGS} FragmentClassifier.cpp

clean:
	rm -rf *.o 

