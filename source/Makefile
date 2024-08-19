
#  Detect the O/S and set a flag accordingly. This flag is used for
#  conditional compilation of some files
#  This is a simplification of stuff I found here:
#    http://stackoverflow.com/questions/714100/os-detecting-makefile
#
ifeq ($(OS),Windows_NT)
   OSFLAG = -D __windows__
	ThisOS := Windows
else
	UNAME_S := $(shell uname -s)
   ifeq ($(UNAME_S),linux)
      OSFLAG = -D __linux__
	endif
endif

#
#  Specify the compiler and any flags I want.
#
#  Notice that I am turning on all warnings (-Wall), but then turning OFF
#  a few categories. I hate seeing extra warnings that make it seem there might
#  be a problem when there is not.  The ones I turn off are:
#   1) unused-function
#      I have defined a few subroutines/functions that are currently unused.
#      That is ok.
#   2) target-lifetime
#      A few of my functions create pointer variables that are passed back to
#      the caller. Obviously, those variables in the caller could potentially
#      remain assigned even after the object gets destroyed. But I don't need
#      an ugly warning message about that.
#   3) surprising
#      This one is a bit more complicated.  When "surprising" is enabled,
#      a large number of warnings are generated that refer to the destructor
#      routines of the various data objects (TDlyData, TDlyDataOneLocation, etc).
#
#      The specific message is
#        Warning: Only array FINAL procedures declared for derived type 'tdlydata' defined at (1), suggest also scalar one [-Wsurprising]
#
#      I have researched this and tried several things to make the compiler "happy", 
#      but have had no success. Best I can tell, my code is valid and correct, though
#      I am certainly willing to accept additional correction if a smarter programmer
#      can provide it.  Given that my tests with memory leak detection software (valgrind)
#      have all come back clean, I am willing to proceed on the assumption that this issue
#      is either:
#        1) A lingering bug in the gfortran compiler (version 5.3.0)
#        2) Weak understanding of the issue on my part, that may "come back to bite me"
#           at some point in the future. But for now it seems to work.
#        3) Only an issue if we were to try using parallelization (this seems to be the
#           focus of the warning, from what I researched).
#      Whatever the case may be, I am moving on with the assumption that I am ok.  I am
#      turning off that class of warning messages.
#
#
CC = gfortran
ifeq ($(ThisOS), Windows)
   CFLAGS = -Wall -Wno-unused-function -Wno-target-lifetime -Wno-surprising -fcheck=all -g -O1 -static
   #CFLAGS = -Wall -Wmaybe-unitialized -g -O1 -static # jak used this for testing
#   CFLAGS = -Wall -Wno-unused-function -Wno-target-lifetime -Wno-surprising -static
else
   CFLAGS = -g -O1 -Wall
  # CFLAGS = -fcheck=all -g -O1
endif


#
#   The executable
#
lbrm:  lbrm.o glshfs_util.o cpp_util.o glshfs_global.o gl_constants.o glerldatatypesandunits.o 
	$(CC) $(CFLAGS) -o lbrm   lbrm.o glshfs_util.o cpp_util.o glshfs_global.o gl_constants.o      \
                      glerldatatypesandunits.o 

#
#   Specific blocks for object code
#
lbrm.o: lbrm.f90  glshfs_util.o cpp_util.o glshfs_global.o gl_constants.o glerldatatypesandunits.o
	$(CC) $(CFLAGS) -c lbrm.f90
   
#
#  The common/utility stuff (data structures, file read/write, etc.)
#  This stuff is mainly to fulfill various dependencies that LBRM has from GLSHFS.
#  Would be nice to reduce/eliminate LBRM dependencies in the future.
#
gl_constants.o: gl_constants.f90
	$(CC) $(CFLAGS) $(OSFLAG) -c gl_constants.f90

glerldatatypesandunits.o: glerldatatypesandunits.f90
	$(CC) $(CFLAGS) -c glerldatatypesandunits.f90						 
						 
#
#  Miscellaneous utility routines
#
#  Note that the source code for glshfs_util and glshfs_global include 
#  pre-processor directives that are conditional on the current O/S.  
#  The gfortran compiler will only respond to pre-processor directives if 
#  the both the file extension is F90 and the compile command specifies 
#  the file with extension of F90 (UPPERCASE F IS REQUIRED). 
#  I have seen discussions on the web that indicate this is a common
#  convention (Intel compiler, etc).
#
glshfs_util.o:  glshfs_util.F90 cpp_util.o
	$(CC) $(CFLAGS) $(OSFLAG) -c glshfs_util.F90

glshfs_global.o: glshfs_global.F90 glshfs_util.o gl_constants.o
	$(CC) $(CFLAGS) $(OSFLAG) -c glshfs_global.F90


#
#  Build special C++ routines that fill gaps in Fortran functionality.
#
cpp_util.o: cpp_util.cpp
	gcc $(OSFLAG) -c -g -O1 cpp_util.cpp

#
#  clean everything up
#         
clean:
	$(RM) *.mod
	$(RM) *.o
 
