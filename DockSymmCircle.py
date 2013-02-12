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


from Default import Parser as R
from Default import Space as S
from Default import Postprocess as PP

import numpy as np
import os, sys
from scipy.cluster.vq import *
from copy import deepcopy

from Protein import Protein
import Multimer as M
import flexibility

class Parser(R):
	def __init__(self):

                #fitness function flags
                self.add('constraint','constraint','str',"NA")
		self.add('target','target','array float',np.array([]))
                self.add('detectClash','detect_clash','str',"on")
                self.add('mixingWeight','mix_weight','float', 0.2)

                #flags for search space and data structure
		self.add('degree','degree','int',-1)
		self.add('radius','radius','float',"NA")

                #postprocessing flags
		self.add('cluster_threshold','cluster_threshold','float',"2")
		self.add('output_folder','output_folder','str',"result")

                #define style (rigid, flexible)
		self.add('style','style','str',"rigid")

                #flexible assembly flags
		self.add('projection','proj_file','str',"NA")
		self.add('align','align','str',"no")
		self.add('ratio','ratio','float',0.9)
		self.add('trajectory','trajectory','str',"NA")
		self.add('trajSelection','trajselection','array str',np.array(["NA"]))
		self.add('topology','topology','str',"NA")

                #rigid assembly flags
		self.add('monomer','pdb_file_name','str',"NA")

                #receptor flags
                self.add('receptor','receptor','str',"NA")
		self.add('boundaryMaxReceptor','high_input_rec','array float',"NA")
		self.add('boundaryMinReceptor','low_input_rec','array float',"NA")
		self.add('z_padding','pad','int',10)


	def check_variables(self):

		#needed, set flag of clash detection
		if self.detect_clash == 'on' :
		    self.detect_clash=True
		elif self.detect_clash== 'off':
		    self.detect_clash=False
		else:
		    print 'ERROR: detectClash should either be "on" or "off"!'
		    sys.exit(1)

		#check if fitness function's energy vs geometry mixing weight is sane
		if self.mix_weight>1 or self.mix_weight<0:
		    print 'ERROR: mixingWeight should be a float within 0 and 1!'
		    sys.exit(1)

		#check if number of monomers is provided
		if self.degree==-1 :
		    print 'ERROR: number of monomer not provided!'
		    sys.exit(1)

		#check if number of target measures are provided
		if len(self.target) == 0 :
		    print 'ERROR: target measures not specified!'
		    sys.exit(1)
	
		#verify existence of constraint file (in search path)
		if self.constraint=="NA":
		    print "ERROR: constraint file not provided"
		    sys.exit(1)
 		if os.path.isfile(self.constraint)!=1:
		    print "ERROR: constraint file %s not found"%self.constraint
		    sys.exit(1)

		#add path of constraint file to search space
		a,b=os.path.split(self.constraint)
		self.constraint=b
		sys.path.append(a)

		#test if constraint file exists, and function constaint_check is available
		try:
			exec 'import %s as constraint'%(self.constraint.split('.')[0])
		except:
		    print "ERROR: load of user defined constraint function failed!"
		    print "ERROR: awww, there may be a syntax error in there..."
		    sys.exit(1)
		try:
		    constraint.constraint_check
		except AttributeError:
		    print 'ERROR: constraint_check function not found in file %s'%self.constraint
		    sys.exit(1)


                #check assembly style flag   
                if self.style!="rigid" and self.style!="flexible":
		    print 'ERROR: style should either be "rigid" or "flexible"!'
		    sys.exit(1)

                ###RIGID ASSEMBLY FLAGS###
                if self.pdb_file_name=="NA" and self.style=="rigid":
                        print "ERROR: in rigid style (default), a monomer should be provided!"
		        sys.exit(1)
                elif self.pdb_file_name!="NA" and self.style!="rigid":
                        print "ERROR: %s style requested"%self.style
                        print "ERROR: usage of monomer flag requires rigid style!"
		        sys.exit(1)
                elif os.path.isfile(self.pdb_file_name)!=1 and self.style=="rigid":
                        print "ERROR: pdb file %s not found"%self.pdb_file_name
		        sys.exit(1)


                ###FLEXIBLE ASSEMBLY FLAGS###
		#verify existence of topology
		if self.topology!="NA":
                    if self.style=="flexible":
		        tmp=os.path.abspath(self.topology)
		        self.topology=tmp
                    else:
                        print "ERROR: %s style requested"%self.style
                        print "ERROR: usage of topology flag requires flexible style!"
		        sys.exit(1)
		if os.path.isfile(self.topology)!=1 and self.style=="flexible":
		    print "ERROR: topology %s not found"%self.topology
		    sys.exit(1)

		#verify existence of traj
		if self.trajectory!="NA":
                    if self.style=="flexible":
                        tmp=os.path.abspath(self.trajectory)
                        self.trajectory=tmp
                    else:
                        print "ERROR: %s style requested"%self.style
                        print "ERROR: usage of trajectory flag requires flexible style!"
		        sys.exit(1)
		if os.path.isfile(self.trajectory)!=1 and self.style=="flexible":
		    print "ERROR: trajectory %s not found"%self.trajectory
		    sys.exit(1)

		#verify existence of projections file
		if self.proj_file!="NA":
                    if self.style=="flexible":
                        tmp=os.path.abspath(self.proj_file)
                        self.proj_file=tmp
                    else:
                        print "ERROR: %s style requested"%self.style
                        print "ERROR: usage of projection flag requires flexible style!"
		        sys.exit(1)
		if self.proj_file!="NA" and os.path.isfile(self.proj_file)!=1 and self.style=="flexible":
		    print "ERROR: projections file %s not found"%self.proj_file
		    sys.exit(1)

                #verify format trajectory selection (if needed)
                if self.trajselection[0]!="NA" and self.style!="flexible":
                        print "ERROR: %s style requested"%self.style
                        print "ERROR: usage of trajSelection flag requires flexible style!"
		        sys.exit(1)
                if self.style=="flexible" and self.trajselection[0]=="NA":
			tmp="protein"
                else:
                        tmp =' '.join(self.trajselection[0:len(self.trajselection)])                        
		self.trajselection=tmp


                #RECEPTOR FLAGS
		if self.receptor!="NA":

			#verify existence of receptor
                	tmp=os.path.abspath(self.receptor)
                	self.receptor=tmp
			if os.path.isfile(self.receptor)!=1:
	        	    print "ERROR: receptor file %s not found"%self.receptor
			    sys.exit(1)

		        #check receptor boundaries consistency
			if len(self.high_input_rec)!=len(self.low_input_rec):
			    print "ERROR: boundaryMinReceptor and boundaryMaxReceptor should have the same length!"
			    sys.exit(1)

			#if len(self.high_input_rec)!=len(self.boundary_type_rec) and self.boundary_type_rec!="NA":
			#    print "ERROR: boundaryTypeReceptor length inconsistent with low and high boundary length!"
			#    print 'ERROR: %s dimensions given, but %s needed!'%(len(self.boundary_type),len(self.low_input))
			#    sys.exit(1)

			if self.high_input_rec!="NA" and (self.low_input_rec>self.high_input_rec).any():
			    print 'ERROR: a receptor lower boundary condition is greated than a higher one!'
		            print self.low_input_rec
		            print self.high_input_rec
		            sys.exit(1)


class Data:

	index=[]
	index_receptor=[]

	def __init__(self,params):
		self.eigenspace_size=""
		
                if params.style=="flexible":
                        print ">> flexible docking requested, launching PCA..."
			#check align flag
			if params.align!="yes" and params.align!="no":
			    print "ERROR: align flag should be set either to \"yes\" or \"no\"!"
			    sys.exit(1)
			
			#load trajectory
			try:
			    self.traj=flexibility.loader(params.trajectory, params.topology, params.align, params.trajselection)
			except ImportError, e:
			    print "ERROR: Trajectory loading failed, aborting!"
			    sys.exit(1)
		    
			#if projection on eigenvectors has been already performed, use the user provided file
			if params.proj_file!="NA":
			    if os.path.isfile(params.proj_file)!=1 :
				print "ERROR: projections file not found!"
				sys.exit(1)
			    else:
				print ">> loading frames projections file..."
				p_file = open(params.proj_file, 'r+')
				p_list=[]
				line = p_file.readline()
				while line:
				    w = line.split()
				    p_list.append(w)
				    line = p_file.readline()
				self.proj=np.array(p_list).transpose().astype(float)
				p_file.close()
				if self.traj.shape[1]!=self.proj.shape[1]:
				    print "ERROR: trajectory has %s frames, projection has %s, aborting!"%(self.traj.shape[1],self.proj.shape[1])
				    sys.exit(1)       
			#else, peform PCA on trajectory
			else:
			    try:
				self.proj,e=flexibility.pca(params.trajectory,params.topology, params.ratio)
			    except ImportError, e:
				print "ERROR: PCA failed, aborting!"
				sys.exit(1)
		       
			#compute projections centroid
			try:
			    self.centroid=flexibility.clusterize(self.proj)

			except ImportError, e:
			    print "ERROR: Clustering failed, aborting!"
			    sys.exit(1)
		
			#suppose that data assigned to cluster can be fitted with a gaussian.
			#interval will be set at 95% of obtained distribution
			p=[]
			for i in xrange(0,len(self.proj),1):
			    p.append(np.std(self.proj[i])*2)
		
			self.eigenspace_size=np.array(p)
			print ">> Flexibility interval in 95% confidence..."
			print "   "+str(self.eigenspace_size)
			   
			#PCA automatically generates a reference pdb file called "protein.pdb"
		        #check its existence and load its name in local params object!
			params.pdb_file_name="protein.pdb"
			if os.path.isfile(params.pdb_file_name)!=1 :
			    print "ERROR: pdb file %s not found"%params.pdb_file_name
			    sys.exit(1)
                
                #load monomeric structure (both for rigid )
		self.structure = Protein()
                try:
		    self.structure.import_pdb(params.pdb_file_name)   
		except IOError:
		    sys.exit(1)
		
                #load receptor structure, if provided, and center it to the origin
                self.receptor_structure=[]
                if params.receptor!="NA":
                    print ">> loading receptor, and centering it to origin..."
		    self.receptor_structure = Protein()
                    try:
		        self.receptor_structure.import_pdb(params.receptor)   
		    except IOError,e:
                        print e
                        print "Man, a loading error!"
		        sys.exit(1)

                    coords=self.receptor_structure.get_xyz()
                    self.receptor_structure.set_xyz(coords-np.mean(coords,axis=0))


		#if clash avoidance is requested, compute indexes of relevant atom
		#this call currently just selects the CA atoms
		if params.detect_clash:
		    #self.index=self.get_index(["CA","CB"])
		    #self.index=self.get_index(["CA"])
		    self.get_index(["CA"])


	def get_index(self,atoms=["CA","CB"]):
        
        	#generate a dummy multimer
        	multimer = M.Multimer(self.structure)
        	multimer.create_multimer(3,10,np.array([0.0,0.0,0.0]))

                #extract the indexes where atoms of interest are located
        	for aname in atoms:
        		[m,p_index]=multimer.atomselect(1,"*","*",aname,True)
        		for i in p_index:
            			self.index.append(i)

                        if self.receptor_structure!=[]:
				[m,p_index]=self.receptor_structure.atomselect("*","*",aname,True)
				for i in p_index:
		    			self.index_receptor.append(i)
		#return self.index   


class Space(S):
	# n dim. space: 3 monomer rotations, 1 monomer radial translation, 1 multimer z translation, 3 receptor rotations, n monomer eigenvectors
	def __init__(self,params,data):

                #add receptor dimensions if needed
                rec_dim=0
                if params.receptor!="NA":
                    rec_dim=4

		self.low=np.zeros(4+rec_dim+len(data.eigenspace_size))
		self.high=np.zeros(4+rec_dim+len(data.eigenspace_size))
		self.cell_size=np.zeros(4+rec_dim+len(data.eigenspace_size))
		self.boundary_type=np.zeros(4+rec_dim+len(data.eigenspace_size))
		
		if len(params.high_input)!=len(params.low_input):
		    print "ERROR: seems that low and high boundaries have different length..."
		    sys.exit(1)
		
		#assign monomer low boundaries
		if params.low_input!="NA" :
		    if len(params.low_input)==3 or len(params.low_input)==4:
		        for i in xrange(0,len(params.low_input),1):
		            self.low[i]=params.low_input[i]
		    else:
		        print "ERROR: boundaryMin should contain within 3 or 4 values"
		        sys.exit(1)
		else:
		    print "WARNING: boundaryMin undefined, using default values for monomer rotation"
		
		#assign monomer high boundaries
		if params.high_input!="NA" :
		    if len(params.high_input)==3 or len(params.high_input)==4:
		        for i in xrange(0,len(params.high_input),1):
		            self.high[i]=params.high_input[i]
		    else:
		        print "ERROR: boundaryMax should contain 3 or 4 values"
		        sys.exit(1)
		else:
		    print "WARNING: boundaryMin undefined, using default values for monomer rotation"
		    self.high[0]=360.0
		    self.high[1]=180.0
		    self.high[2]=360.0
		
		#if needed assign radius boundaries
		if params.radius!="NA":
		    if params.high_input=="NA":
		        self.low[3]=params.radius
		        self.high[3]=params.radius
		    else:
		        print "ERROR: radius keyword used when also using boundary condition list!"
		        sys.exit(1)
		if params.radius=="NA" and params.high_input=="NA":
		    print "ERROR: translation undefined boundary or radius keywords needed!"
		    sys.exit(1)


                #if needed, add receptor dimensions
                if params.receptor!="NA":

			#multimer z-axis translation
			if params.high_input_rec=="NA":
		            d=data.structure.get_xyz()
		            self.low[4]=d[2].min()+params.pad
			    self.high[4]=d[2].max()+params.pad
		        else:
			    self.low[4]=params.low_input_rec[0]
			    self.high[4]=params.high_input_rec[0]

		            #receptor x,y,z axis rotations
			    if len(params.high_input_rec)!=4:
			        print "WARNING: Receptor rotations undefined in boundaryMinReceptor, using default values"
			       	self.low[5]=0.0
				self.low[6]=0.0
				self.low[7]=-360/(2*params.degree)
				self.high[5]=0
				self.high[6]=0
     				self.high[7]=360/(2*params.degree)
                            else:
			        self.low[5]=params.low_input_rec[1]
			        self.high[5]=params.high_input_rec[1]
			        self.low[6]=params.low_input_rec[2]
			        self.high[6]=params.high_input_rec[2]
			        self.low[7]=params.low_input_rec[3]
			        self.high[7]=params.high_input_rec[3]

   
		#add eigenvector fluctuations in search space
		for i in xrange(0,len(data.eigenspace_size),1):
		    self.low[4+rec_dim+i]=-data.eigenspace_size[i]
		    self.high[4+rec_dim+i]=data.eigenspace_size[i]
	    

		#final check for all boundary conditions consistency (cause we're paranoid)
		if (self.low>self.high).any():
		    print 'ERROR: a lower boundary condition is greated than a higher one'
		    sys.exit(1)

		#define cell size
		self.cell_size=self.high-self.low
		      
		#set boundary type
		if params.boundary_type=="NA":
		    for i in xrange(0,len(self.low),1):
		        self.boundary_type[i]=0
		elif params.boundary_type!="NA" and len(params.boundary_type)!=len(self.low):
		    print 'ERROR: boundaries type inconsistent with system dimensions'
		    print 'ERROR: %s dimensions given, but %s needed!'%(len(params.boundary_type),len(self.low))
		    sys.exit(1)
		else:
		    for i in xrange(0,len(self.low),1):
		        self.boundary_type[i]=params.boundary_type[i]


class Fitness:

    def __init__(self,data,params):

        self.style=params.style
        self.clash=params.detect_clash
        self.target=params.target
        
        self.constraint=params.constraint.split('.')[0]
        #test if constraint file exists, and function constaint_check is available
        try:
            exec 'import %s as constraint'%(self.constraint)
        except ImportError, e:
            print "ERROR: load of user defined constraint function failed!"
            sys.exit(1)
    
        try:
            constraint.constraint_check
        except AttributeError:
            print 'ERROR: constraint_check function not found'
            sys.exit(1)

        #data to manipulate
        self.data=data
        self.degree=params.degree

        #receptor space size
        self.rec_dim=0
        if params.receptor!="NA":
            self.rec_dim=4

        #energy vs geometry mixing weight
        self.c1=params.mix_weight

    
    def evaluate(self,num,pos):
        

        #import user provided constraint file and class needed to construct multimer
        exec 'import %s as constraint'%(self.constraint)
        import Multimer as M
        
        if self.style=="flexible":
		#pick monomeric structure from database
		deform_coeffs=pos[(4+self.rec_dim):len(pos)]
		pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
		code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
		target_frame=min_dist.argmin()
		coords=self.data.traj[:,target_frame]
		coords_reshaped=coords.reshape(len(coords)/3,3)
		self.data.structure.set_xyz(coords_reshaped)      
        
        #assemble multimer
        self.multimer = M.Multimer(self.data.structure)
        self.multimer.create_multimer(self.degree,pos[3],pos[0:3])

        measure=[]
        if self.data.receptor_structure!=[]:
                #shift multimer on z axis
		self.multimer.z_to_origin()
		self.multimer.z_shift(pos[4])

                #rotate receptor (backing up original position, and putting it back when measure is done)
                coords=self.data.receptor_structure.get_xyz()
		self.data.receptor_structure.rotation(pos[5],pos[6],pos[7])
                measure = constraint.constraint_check(self.multimer,self.data.receptor_structure)
                self.data.receptor_structure.set_xyz(coords)
        else:
                try:
        	    measure = constraint.constraint_check(self.multimer)
                except:
                    print "ERROR: syntax error in constraint file!"
                    raise

        if len(measure) != len(self.target) :
            print 'ERROR: constraint file produced %s measures, but %s target measures are provided!'%(len(measure),len(self.target))
            sys.exit(1)
        
        #measure distance within target values and obtained values
        diff=self.target-np.array(measure)
        distance=np.sqrt(np.dot(diff,diff))
  
        #detect clashes
        if self.clash:
            energy=self.interface_vdw()
            return self.c1*energy+(1-self.c1)*distance  
        else:
            return distance


    def interface_vdw(self):
    
        epsilon=1.0
        sigma=4.7
        cutoff=12.0
    
        energy=0

        #CLASH WITH NEIGHBORING MONOMER
        #extract coords of monomers 1 and 2 of multimeric structure according to desired atoms

        m1=self.multimer.get_multimer_uxyz()[0][self.data.index]
        m2=self.multimer.get_multimer_uxyz()[1][self.data.index]
        
        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
    
        dist=np.array(d)
    
        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues
        
        for i in xrange(0,len(couples[0]),1):    
            d=dist[couples[0,i],couples[1,i]]
            energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)
 
	#CLASH WITH RECEPTOR (if needed)
        if self.data.receptor_structure!=[]:
		m1=self.data.receptor_structure.get_xyz()[self.data.index_receptor] 
		
		#extract distances of every atom from all the others
		d=[]
		for i in xrange(0,len(m1),1):
		    d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
	    
		dist=np.array(d)
	    
		#detect interfacing atoms (atom couples at less than a certain cutoff distance
		couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues

		for i in xrange(0,len(couples[0]),1):    
		    d=dist[couples[0,i],couples[1,i]]
		    energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)
    
        return energy



class Postprocess(PP):

	def __init__(self,data,params):
		self.data=data
		self.params=params

		self.rec_dim=0
		if params.receptor!="NA":
		    self.rec_dim=4

		#load constraint file
		self.constraint=params.constraint.split('.')[0]
		try:
		    exec 'import %s as constraint'%(self.constraint)
		except ImportError, e:
		    print "ERROR: load of user defined constraint function failed!"
		    sys.exit(1)
	    
		try:
	            constraint.constraint_check
		except NameError:
		    print 'ERROR: constraint_check function not found'

	#clustering according to rmsd of solutions in search space
	#threshold2 = clustering threshold
	def run(self) :
	    
            exec 'import %s as constraint'%(self.constraint)

            #create output directory for generated PDB
            self.OUTPUT_DIRECTORY=self.params.output_folder
	    if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
		os.mkdir(self.OUTPUT_DIRECTORY)

	    clusters_file=open("%s/solutions.dat"%self.params.output_folder,"w")

	    #use superclass method to filter acceptable solutions
            self.log=self.select_solutions(self.params)
            print ">> %s solutions filtered"%len(self.log[:,1])
            if len(self.log[:,1])==0:
                return
 
	    #generate a dummy multimer and extract the indexes of C alpha
	    multimer = M.Multimer(self.data.structure)
	    multimer.create_multimer(2,10,np.array([0.0,0.0,0.0]))
	    [m,index]=multimer.atomselect(1,"*","*","CA",True)

            #if needed, extract receptor indexes too
            if self.data.receptor_structure!=[]:
	        [r1,index_receptor]=self.data.receptor_structure.atomselect("*","*","CA",True)

            #load the monomeric structure positions (needed for resetting atom position after displacement)
	    s = Protein()
	    s.import_pdb(self.params.pdb_file_name)   
            coords=s.get_xyz()

	    print ">> clustering best solutions..."
	    P=self.log[:,0:len(self.log[0,:])] #points list
	    V=self.log[:,-1] #values of best hits
	    C=[] #centroids array
	    P_C=np.zeros(len(P)) #points-to-cluster mapping
	    C_V=[] #centroids values
	    cnt=0 #centroids counter

	    #cluster accepted solutions
	    while(True) :
		#check if new clustering loop is needed
		k=np.nonzero(P_C==0)[0]
		if len(k)!=0 :
		    cnt=cnt+1
		    P_C[k[0]]=cnt
		    a=P[k[0]]
		    C.append(a)
		else :
		    break

                if self.params.style=="flexible":
			#build reference structure
			deform_coeffs=np.array(C[cnt-1])[(4+self.rec_dim):len(C[cnt-1])-1]

			pos=self.data.proj[:,self.data.centroid]+deform_coeffs
			code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
			target_frame1=min_dist.argmin()
			coords=self.data.traj[:,target_frame1]
			coords_reformat=coords.reshape(len(coords)/3,3)
			self.data.structure.set_xyz(coords_reformat)
                else:
			self.data.structure.set_xyz(coords)  

		#create multimer
		pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
		multimer1 = M.Multimer(self.data.structure)
       		multimer1.create_multimer(self.params.degree,pos[3],pos[0:3])
        	m1=multimer1.get_multimer_uxyz()[0][index]

		if self.data.receptor_structure!=[]:
		        #shift multimer on z axis
			multimer1.z_to_origin()
			multimer1.z_shift(pos[4])

		        #rotate receptor (backing up original position, and putting it back when measure is done),
                        #also write the rotated receptor structure
		        coords_tmp=self.data.receptor_structure.get_xyz()
			self.data.receptor_structure.rotation(pos[5],pos[6],pos[7])
                        r1=self.data.receptor_structure.get_xyz()[index_receptor]
                        self.data.receptor_structure.write_pdb("%s/receptor%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))
		        self.data.receptor_structure.set_xyz(coords_tmp)
                        m1=np.concatenate((m1,r1))

		#write multimer
		multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))

		#clustering loop
                cnt2=1
		for i in xrange(0,len(k),1) :

                    if self.params.style=="flexible":
			    deform_coeffs=np.array(P[k[i]])[(4+self.rec_dim):len(P[k[i]])-1]
			    pos=self.data.proj[:,self.data.centroid]+deform_coeffs
			    code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
			    target_frame=min_dist.argmin()
			    coords=self.data.traj[:,target_frame]
			    coords_reformat=coords.reshape(len(coords)/3,3)
			    self.data.structure.set_xyz(coords_reformat)
                    else:
			    self.data.structure.set_xyz(coords)  

		    multimer2 = M.Multimer(self.data.structure)
		    multimer2.create_multimer(self.params.degree,P[k[i]][3],np.array([P[k[i]][0],P[k[i]][1],P[k[i]][2]]))
		    m2=multimer2.get_multimer_uxyz()[0][index]

		    if self.data.receptor_structure!=[]:
		    	#shift multimer on z axis
		    	multimer2.z_to_origin()
			multimer2.z_shift(P[k[i]][4])

		        #rotate receptor (backing up original position, and putting it back when measure is done)
		        coords_tmp=self.data.receptor_structure.get_xyz()
			self.data.receptor_structure.rotation(P[k[i]][5],P[k[i]][6],P[k[i]][7])
                        r2=self.data.receptor_structure.get_xyz()[index_receptor]
		        self.data.receptor_structure.set_xyz(coords_tmp)
                        m2=np.concatenate((m2,r2))

                    #compute RMSD
		    rmsd=self.align(m1,m2)
                  		        
		    if rmsd<self.params.cluster_threshold :
                        cnt2+=1
		        P_C[k[i]]=cnt
		
                if self.params.style=="rigid":
		    print ">>> clustered %s solutions on multimer %s"%(cnt2-1,cnt)
                if self.params.style=="flexible":
		    print ">>> clustered %s solutions on multimer %s (from frame %s)"%(cnt2-1,cnt,target_frame1)
  
		#set centroid score with score of closes neighbor in set
		q=np.nonzero(P_C==cnt)[0]
		distance=10000
		targ=0
		for i in xrange(0,len(q),1) :
		    d=np.sqrt(np.dot(C[cnt-1]-P[q[i]],C[cnt-1]-P[q[i]]))         
		    if d<distance :
		        distance=d
		        targ=q[i]
		C_V.append(V[targ])
		
                #extract constraint values calculated for selected centroid
        	measure = constraint.constraint_check(multimer1)

		###generate output log (prepare data and formatting line, then dump in output file)###
                l=[]
                f=[]
		for item in C[cnt-1][0:len(C[cnt-1])-1]:
		    l.append(item)
                    f.append("%8.3f ")
                #write constraint values
                f.append("| ")
		for item in measure:
		    l.append(item)
                    f.append("%8.3f ")
                #write fitness
                f.append("| %8.3f\n")
                l.append(C_V[cnt-1])

                formatting=''.join(f)

		clusters_file.write(formatting%tuple(l))


	    clusters_file.close()

	    return









