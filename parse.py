#!/usr/bin/env python

# Copyright (c) 2012 EPFL (Ecole Polytechnique federale de Lausanne)
# Laboratory for Biomolecular Modeling, School of Life Sciences
#
# POW is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# POW is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with POW ;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Author : Matteo Degiacomi, matteothomas.degiacomi@epfl.ch
# Web site : http://lbm.epfl.ch


import os, sys
import numpy as np

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


if rank == 0:
    print '\n>        Parallel Optimization Workbench (POW)         <'
    print   '> (c) 2012, Laboratory for Biomolecular Modeling, EPFL <'
    print   '>                                                      <'
    print   '>               -postprocess results -                 <\n'


#check input params consistency
if len(sys.argv)!=3 and len(sys.argv)!=4:
    if rank == 0:
        print "\nERROR: parameters are not what I expected!"
        print "USAGE: ./parse.py module input_file [logfile]\n"
    sys.exit(1)
 
#get input file
infile=str(sys.argv[2])
if os.path.isfile(infile)!=1 :
    if rank == 0:
        print "ERROR: input file not found!"
    sys.exit(1)

#get program installation location and declare an environment variable (needed by some modules)
h=os.path.dirname(sys.argv[0])
home_dir=os.path.abspath(str(h))
os.environ['POW']=home_dir
sys.path.append(home_dir)

#get local working directory and add to pythonpath
working_dir=os.getcwd()
sys.path.append(working_dir)


#add module path in path, in case it's not in POW home
p,mod=os.path.split(sys.argv[1])
if len(p)>0:
    sys.path.append(p)

#preprocessing performed only by master node
if rank == 0:

    #load appropriate module
    exec 'import %s as mode'%(mod.split('.')[0]) 
    
    print ">> parsing input file"
    #parse parameters
    params=mode.Parser()
    params.add_standard()
    params.set_default_values()
    params.parse(infile)
    params.check_standard_variables()
    params.check_variables()
    
    #checking logfile existence
    if len(sys.argv)==4:
        params.output_file=str(sys.argv[3])
    if os.path.isfile(params.output_file)!=1 :
        print "ERROR: input file not found!"
        sys.exit(1)
    #prepare datastructures
    print ">> loading data structures"
    data=mode.Data(params)

    #run postprocessing
    print ">> postprocessing logfile %s"%params.output_file
    post=mode.Postprocess(data,params)

else:
    params=None
    data=None
    post = None

#propagate parameters, data, space and fitness function to slaves
comm.Barrier()
params=comm.bcast(params,root=0)
data=comm.bcast(data,root=0)
post=comm.bcast( post,root=0)
comm.Barrier()

post.run()
