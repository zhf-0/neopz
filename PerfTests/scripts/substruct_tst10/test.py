#! /usr/bin/env python2.7
#***************************************************************************
#*   Copyright (C) 2013 by Edson Borin                                     *
#*   edson@ic.unicamp.br                                                   *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 2 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#**************************************************************************/

# ---------------------------------------------
# Performance test module

# List of rdt files generated by the test
rdtfiles_l = [
        # short_name, option, filename, description
        ("ass", "-ass_rdt", "ass.rdt", "Assemble: dohrstruct->Assemble(...). Assemble element matrices and decompose matrices."),
        ("cre", "-cre_rdt", "cre.rdt", "Create: dohrstruct->Create()"),
        ("sol", "-sol_rdt", "sol.rdt", "Solver: cg.Solve(...)"),
        ("tot", "-tot_rdt", "tot.rdt", "Total: all the steps"),
        ("dohrass", "-tpz_dohr_ass", "tpzdohrass.rdt", "Assemble element matrices"),
        ("dohrdec", "-tpz_dohr_dec", "tpzdohrdec.rdt", "Decompose matrices")
]

def short_description() : return "substructure -- 8andares02.txt -- num_it 1000 - p2 - nsub 128 - serial"

def long_description():
	desc =  "Execute the substruct tool collecting statistics for the following steps:"
	for rdtarg in rdtfiles_l :
		desc = desc + '\n\t' + rdtarg[0] + ' (' + rdtarg[1] + ' ' + rdtarg [2] +  ') : ' + rdtarg[3]
	return desc
# ---------------------------------------------

import sys
import os.path
import shlex, subprocess
import resource

# Try to import rdt and stats modules, if available.
import sys

# Variables to be defined by cmake
builddir="@PERFTEST_BASE_DIR@"
datadir="@PERFTEST_SMALL_DATA_DIR@"

def error(message, status):
	sys.stderr.write('ERROR (test.py): '+message+'\n')
        sys.exit(status)

# Setup the command line
def setup_cmd():
	# Check build directory
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 1)
	# Check run directory
	rundir = os.path.join(builddir,'scripts','substruct_tst10')
	if not os.path.isdir(rundir) :
		error(rundir+' is an invalid run directory.', 1)
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 1)
	# Check executable
	executable=os.path.join(builddir,"progs","substruct", "substruct-perf")
	if not os.path.isfile(executable) :
		error(executable+' is an invalid executable file name.', 1)
	# Check input file
	inputfn = os.path.join(datadir,"substruct","inputs","8andares02.txt")
	if not os.path.isfile(inputfn) :
		error(inputfn+' is an invalid input file name.', 1)	
	# Put the arguments together
        arguments = ' -mp '+inputfn
	arguments = arguments + ' -num_it 1000'
	#NUMA aware Dohrman Assembly List thread work objects re-allocation.
	#arguments = arguments + ' -naDALora'
	#NUMA aware Dohrman Assembly List thread work objects re-allocation threshold.
	#arguments = arguments + ' -naDALorat 1835008' # 2/2MB(l2) + 6/8MB(l3)
	#NUMA aware (node round-robin) Dohrman Assembly List thread work scheduling.
	#arguments = arguments + ' -naDALtws' 
	arguments = arguments + ' -nsub 128'
	arguments = arguments + ' -nt_a 0' 
	arguments = arguments + ' -nt_d 0' 
	arguments = arguments + ' -nt_m 0' 
	arguments = arguments + ' -nt_sm 0' 
	arguments = arguments + ' -p 2' 
	for rdtarg in rdtfiles_l :
		arguments = arguments + ' ' + rdtarg[1] + ' ' + rdtarg[2]
	# TODO: Add arguments to enforce output checking!
	return rundir, executable+arguments

# Limits for this test
# 38400 = 64 (cores) * (60) * (10) = 10 minutes in 64 cores.
limits = { "cpu"   : (resource.RLIMIT_CPU,  38400, "Max CPU user time in seconds (not wall clock time)"), 
#	   "nofile": (resource.RLIMIT_NOFILE,   7, "The maximum number of open file descriptors for the current process."),
#	   "rss"   : (resource.RLIMIT_RSS,   1024, "The maximum resident set size that should be made available to the process"),
#	   "fsize" : (resource.RLIMIT_FSIZE,    1, "Max size of a file which the process may create"),
#	   "data"  : (resource.RLIMIT_DATA,  1024, "The maximum size (in bytes) of the process's heap"),
#	   "nproc" : (resource.RLIMIT_NPROC,    0, "The maximum number of processes the current process may create")
	 }

# Set the rlimits of the chidren process (see limits above)
# TODO: Improve the handling of sandboxing limits
def setlimits():
	print "Setting resource limit in child"
	for k, v in limits.iteritems() : 
		resource.setrlimit(v[0], (v[1],v[1])) 
		#print k, " : ", v[0], " => ", v[1]

# Sumarizes the RDT (Raw data table) files information
def sumarize_rdt_files(rundir) :
	results = {}
	for f in rdtfiles_l : 
		rdt_id  = f[0]   # Step name
		rdt_fn  = os.path.join(rundir,f[2]) # RDT file name
		rdt_dsc = f[3]   # Description
		results[rdt_id] = (rdt_fn, rdt_dsc)
	return results

# Execute the test.
def run_test(ntimes):
	rundir,cmd=setup_cmd()
	args = shlex.split(cmd)
	sout = None
	serr = None
	for i in range(ntimes) : 
		p = subprocess.Popen(args, preexec_fn=setlimits, stdout=sout, stderr=serr, cwd=rundir)
		p.wait()
		if (p.returncode != 0) : 
			return p.returncode, {}
	results = sumarize_rdt_files(rundir)
	return 0, results

