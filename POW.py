#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


from PSO import PSO

import wx
import os, sys
import time
import numpy as np
import ClusterAndDraw as CnD

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# import the post processor and tree builder


#check wether an appropriate number of parameters has been provided
if rank == 0:

    print '\n>        Parallel Optimization Workbench (POW)         <'
    print   '> (c) 2012, Laboratory for Biomolecular Modeling, EPFL <'

    # checking if input is correct and has 3 arguments
    if len(sys.argv)!=3:
        print "\nERROR: expected 2 parameters as input!"
        print "USAGE: mpiexec -n 2 POW.py module_name input_file\n"
        sys.exit(1)

#get input file
infile=str(sys.argv[2])

# ----- EXITING FROM ERROR
if rank == 0 and os.path.isfile(infile)!=1 :
    print "ERROR: input file not found!"
    sys.exit(1)
# -----

#get program installation location and declare an environment variable (needed by some modules)
h=os.path.dirname(sys.argv[0])
home_dir=os.path.abspath(str(h))
os.environ['POW']=home_dir
sys.path.append(home_dir)

#get local working directory and add to pythonpath
working_dir=os.getcwd()
sys.path.append(working_dir)

#preprocessing performed only by master node
if rank == 0:


    #init timer
    start=time.time()

    print '\n> PREPROCESSING %s\n'%sys.argv[1].split('.')[0]

    exec 'import %s as mode'%(sys.argv[1].split('.')[0]) # Xx mode is DockSymmCircle file! xX

    #parse input file
    print '>> parsing input file...'
    params=mode.Parser() #read user defined variables # Xx Params is class Parser of DockSymmCircle file xX
    params.add_standard() #add default variables _> in default as well <-> self.add_standard of Parser Class
    params.set_default_values() #set default values to all defined variables _> in default
    params.parse(infile) #parse input file _> in default
    params.check_standard_variables() #check consistency of standard variables _> in default
    params.check_variables() #check consistency of user defined variables _> in DockSymmCircle

    # ------------------------------------------------------------- TEST AREA -------------------------------------------------------------



    # -------------------------------------------------------------------------------------------------------------------------------------


    #load requested data structures
    print '>> importing data...'
    data=mode.Data(params) # Xx data is the class data of docksym with argument params from above,
    # -Xx so loaded into data all the parameters self.blabla of Parser class, in the class itself these are loaded in the
    # -Xx def __init__ (self, params)function
    # -Xx the Data class will then do appropriate modications, such checking whether there is conformational change etc

    #build search space
    print ">> generating search space..."
    space=mode.Space(params,data)
    print ">> "+str(len(space.low))+"D space generated (min_pos, max_pos):"
    for i in xrange(0,len(space.low),1):
        print "   %s, %s"%(space.low[i],space.high[i])

    #charge module fitness function or the one requestyed by user
    if params.fit=="NA":
        print ">> creating fitness function..."
        fitness=mode.Fitness(data,params)
    else:
        print ">>> importing user defined fitness function"
        try:
            exec 'import %s as user_fit'%(params.fit.split('.')[0])
        except ImportError, e:
            print "ERROR: load of user defined fitness function failed!"
            sys.exit(1)
        try:
            user_fit.Fitness
        except NameError:
            print 'ERROR: user defined fitness function file should contain a class Fitness with a function evaluate(id,values_array) '
        fitness=user_fit.Fitness(data,params)

else:
    params=None
    data=None
    space=None
    fitness=None
    post = None

#propagate parameters, data, space and fitness function to slaves
comm.Barrier()
params=comm.bcast(params,root=0)
space=comm.bcast(space,root=0)
fitness=comm.bcast(fitness,root=0)
data=comm.bcast(data,root=0)
comm.Barrier()

#prepare optimizer
search=PSO(params,space,fitness) # -> the init function of PSO

#init optimization timer
if rank==0:
    start2=time.time()

search.launch()

if rank==0:
    end2=time.time()
    one_step=end2-start2
    print(">> optimization time (sec): {t:<8.3f}".format(t=one_step))

    #postprocess
    print '\n> POSTPROCESSING...'
    post=mode.Postprocess(data,params)

comm.Barrier()
post=comm.bcast( post,root=0)
comm.Barrier()

post.create_distance_matrix()

    #post.run() # is now replaced by creating the distance matrix


if rank == 0:
    end=time.time()
    one_step=end-start
    print("\n>> total execution time (sec): {t:<8.3f}\n".format(t=one_step))

    # Code to create clusters, open interface and select PDBs
    app = wx.App(False)
    frame = CnD.MainFrame(None, "Clustering interface",post ,params, data)
    frame.RMSDPanel.computeMatrix()
    if params.cluster_threshold == "NA":
        frame.Show()
        app.MainLoop()
    else:
        frame.RMSDPanel.convertCoordsAndExportPDB(params.cluster_threshold)
