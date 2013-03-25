# -*- coding: utf-8 -*-
import wx

#numpy
import numpy as np
from scipy.cluster.vq import *

import copy
import string

#import Multimer as M
#import AssemblyHeteroMultimer as AH # _CHANGED_
#from Protein import Protein
#import flexibility_new as F # _CHANGED_

import os, sys


class App(wx.App):

    def OnInit(self):
        self.frame = MainFrame(None, "RMSD Viewer")
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True

class MainFrame (wx.Frame):
    def __init__(self, parent, title, outputDir, params, data, post): # _CHANGED_
        wx.Frame.__init__(self, parent, title=title, size=(1000,600))

        self.post = post
        self.outputDir = outputDir # _CHANGED_
        self.params = params # _CHANGED_
        self.data = data # _CHANGED_


        # changing color of backgroundl
        self.SetBackgroundColour("black")

        # inserting the panels:
        self.treePanel = treeDisplay(self)
        self.RMSDPanel = MakeWork(self)

        # adding the panels to the main frame
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.RMSDPanel, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.treePanel, 0, wx.EXPAND | wx.ALL, 3 ) # the 3 is the padding compared to the frame!
        self.SetSizer(box)
        self.Centre()

# ------------------------------------------------------------- WORKING CLASS

class MakeWork (wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="Draw and Select Centroids")
        self.title.SetFont(font1)

        # ---------- INPUT FIELDS AND BUTTONS
        self.lblRMSD=wx.StaticText(self, label="RMSD threshold:")
        self.fieldRMSD=wx.TextCtrl(self, value="")
        self.buttonDrawTree=wx.Button(self, label="Draw Tree")
        self.buttonShowRMSD=wx.Button(self, label="Select RMSD")
        self.buttonExportFile=wx.Button(self, label="Export PDB")
        self.buttonSavePicture=wx.Button(self, label="Save as Png")
        self.Bind(wx.EVT_BUTTON, self.selectRMSD, self.buttonShowRMSD)
        self.Bind(wx.EVT_BUTTON, self.clickToConvertPDB, self.buttonExportFile)
        self.Bind(wx.EVT_BUTTON, self.DrawRMSDTree, self.buttonDrawTree)
        self.Bind(wx.EVT_BUTTON, self.saveSnapshot, self.buttonSavePicture)

        # ---------- MATRIX CONTAINER
        RMSDContainer = wx.FlexGridSizer(1, 6, 3, 3) # flexGridSzer (self, rows, cols, vgap, hgap) info at http://wxpython.org/docs/api/wx.FlexGridSizer-class.html
        RMSDContainer.SetFlexibleDirection(wx.HORIZONTAL) # specifies that rows are resized dynamically
        RMSDContainer.AddGrowableCol(4,1) # specifies that the column 1(starting from 0) can be regrown, ids specified below!
        RMSDContainer.AddMany([ # here insert them in order plaase
                    (self.buttonDrawTree, 0,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonSavePicture, 0,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonExportFile, 0,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.lblRMSD, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.fieldRMSD, 0,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonShowRMSD, 0,wx.ALIGN_LEFT | wx.EXPAND, 3),


                    ])

        # ----------- BOX CONTAINER
        BigRMSDContainer = wx.BoxSizer(wx.VERTICAL) # more info about box sizing at http://zetcode.com/wxpython/layout/
        BigRMSDContainer.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        BigRMSDContainer.Add(RMSDContainer, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(BigRMSDContainer)
        self.Centre()

    def DrawRMSDTree (self, event):
        self.Parent.treePanel.plotFigure(self.arrayReadyforDrawing)

    #def DrawFromPOW(self):
        #self.Parent.treePanel.plotFigure(self.arrayReadyforDrawing)
        #print "dendogram drawn"

    def selectRMSD (self, event):

        self.selectedRMSD = self.fieldRMSD.GetValue()
        self.Parent.treePanel.plotRMSDLine(self.selectedRMSD)

    def saveSnapshot(self, event):
        # based largely on code posted to wxpython-users by Andrea Gavana 2006-11-08
        dcSource = wx.PaintDC(self.Parent.treePanel)
        size = dcSource.Size

        # Create a Bitmap that will later on hold the screenshot image
        # Note that the Bitmap must have a size big enough to hold the screenshot
        # -1 means using the current default colour depth
        bmp = wx.EmptyBitmap(size.width, size.height)

        # Create a memory DC that will be used for actually taking the screenshot
        memDC = wx.MemoryDC()

        # Tell the memory DC to use our Bitmap
        # all drawing action on the memory DC will go to the Bitmap now
        memDC.SelectObject(bmp)

        # Blit (in this case copy) the actual screen on the memory DC
        # and thus the Bitmap
        memDC.Blit( 0, # Copy to this X coordinate
            0, # Copy to this Y coordinate
            size.width, # Copy this width
            size.height, # Copy this height
            dcSource, # From where do we copy?
            0, # What's the X offset in the original DC?
            0  # What's the Y offset in the original DC?
            )

        # Select the Bitmap out of the memory DC by selecting a new
        # uninitialized Bitmap
        memDC.SelectObject(wx.NullBitmap)

        img = bmp.ConvertToImage()
        img.SaveFile("%s/tree.png"%(self.Parent.outputDir), wx.BITMAP_TYPE_PNG)
        print ">> saved snapshot of dendogram tree"

    def clickToConvertPDB(self, event):
        self.convertCoordsAndExportPDB(self.selectedRMSD)

    def convertCoordsAndExportPDB (self, selectedRMSD):
        print ">> converting Coordinates of selected Centroids into PDB files"

        centroidArray = [] # this array will hold the indexes of the centroids
        average_RMSD_ARRAY = [] # this one will hold the average RMSD values of clusters in the same order as the centroidArray

        sortedRMSDArray =  copy.deepcopy(self.Parent.RMSDPanel.RMSDThresholdHash.keys())# this contains the sorted RMSDs of the RMSD threshold hash containing the RMDS and the corresponding centroid
        sortedRMSDArray.sort()

        

        # parse the sorted array and return an array of the centroids and within RMSD values with RMSDs below threshold
        for e in range(len(sortedRMSDArray)):
            if float(selectedRMSD) < float(sortedRMSDArray[e]):

                centroidArray = copy.deepcopy(self.Parent.RMSDPanel.RMSDThresholdHash[sortedRMSDArray[e-1]][0])
                average_RMSD_ARRAY = copy.deepcopy(self.Parent.RMSDPanel.RMSDThresholdHash[sortedRMSDArray[e-1]][1])
                break

        # in case the RMSD value is higher than the RMSDs in the sorted arrays of RMSDs:
        if len(centroidArray) == 0:
            centroidArray = copy.deepcopy(self.Parent.RMSDPanel.RMSDThresholdHash[sortedRMSDArray[-1]][0])
            average_RMSD_ARRAY = copy.deepcopy(self.Parent.RMSDPanel.RMSDThresholdHash[sortedRMSDArray[-1]][1])
            
        self.Parent.post.write_pdb( centroidArray,average_RMSD_ARRAY)


    def computeMatrix (self):

        #self.Parent.treePanel.Center()

        # import the files into arrays using numpy
        distanceMatrix = np.loadtxt('np_matrix.txt')
        coordinateArray = np.loadtxt('coordinateArray.txt')

#        distanceMatrix = self.Parent.post.distanceMatrix
#        coordinateArray = self.Parent.post.coordinateArray

        self.distanceMatrix = distanceMatrix#[0:50]
        self.coordinateArray = coordinateArray#[0:50,0:50]

        # create a dummy of the distanceMatrix
#        self.distanceMatrixDummy = self.distanceMatrix[:,0:len(self.distanceMatrix[0,:])]
        self.distanceMatrixDummy = copy.deepcopy(self.distanceMatrix)

        # create a hash to be sorted that will contain the RMSD values as well as the pairwise coodinates
        self.unSortedHash = {}
        self.clusterHash = {} # the first value of the self.clusterHash will be the centroid of the cluster and the second will be the distanceMatrix index of structures

        # keep track of all the index that have already been added in a cluster
        self.clusteredIndex = []

        # ------------------------------------------------------ CREATING CLUSTERS

#        self.sortedHashKeys = self.unSortedHash.keys()
#        self.sortedHashKeys.sort()
        foundStructureInCluster = False # a boolean that lets you know whether you have either of struture in one of the clusters
        NewAddedIndex = '' # the index of the structure found when iterating over all the clusters
        clusterCount = 1
        nonCentroids = [] # this array will be used to remove the non centroid elements in the sorted hash
        ArrayedRMSDIndexes = []
        valuesToBeRemoved = []
        self.centroidArray = [] # this will contain an array of indexes which are the centroids of all the clusters
        ClusterMergingSwitch = False
        centroidIndexForMerge = 0 # used to know which index is the centroid when one cluster merges with a smaller one
        self.maximalRMSD = 0.00 # was used to stop clustering at the certain point but is now used to compute the graduation bar in the dendogram

        self.iterant = 0 #### THIS IS THE ITERANT VARIABLE USED TO DRAW THE DENDOGRAM ####
        self.DrawingArray = [] # this array will host all the small arrays to be drawn
        arraylist1copy = []
        arraylist2copy = [] # this is used to transfer the data between the first and second merged array when drawing
        self.RMSDMutliplier = 0.00
        topRMSDBeforeRemoval = 0 # this is to store the RMSD value before being removed and using it to make functions
        self.withinBetween = True # this value will help you to know whether the within distance is always smaller than the between distance
        indexAway = 0

        # TAKING PDBs FROM RMSD THRESHOLD
        self.RMSDThresholdHash = {} # this hash will contain the centroids corresponding to each passed RMSD threshold

        # Numpy matrix shortcut variables:
        self.minimumIndex = 0
        self.firstMinimumIndex = 0
        self.secondMinimumIndex = 0


        # as long as the self.sortedHash.Keys is no empty continue assigning clusters
        
        while self.distanceMatrixDummy.min() != 100:
            
            self.updataMinimumValue()
#            self.distanceMatrixDummy[np.where(self.distanceMatrixDummy == self.distanceMatrixDummy.min())])
            #for lolo in range(30):

            # ======================================== WHEN CLUSTER IS EMPTY  ===================================================================

            # if there is nothing in the self.clusterHash, create your first cluster
            if len(self.clusterHash.keys()) == 0:
                self.clusterHash[clusterCount] = [None , [ self.firstMinimumIndex, self.secondMinimumIndex ]]
                # adding them to the clustered index hash
                self.clusteredIndex += [self.firstMinimumIndex , self.secondMinimumIndex]

                #print "removing first element "+ str(self.clusterHash[clusterCount][1]) +" from list "

#                self.DrawingArray += [ [clusterCount,[ self.makeSingleArray(self.distanceMatrixDummy.min())]] ]
                self.DrawingArray += [ [clusterCount, [["ini", self.makeSingleArray(self.distanceMatrixDummy.min())]] ] ]


                # remove value in the Matrix
                self.distanceMatrixDummy[self.minimumIndex] = 100

                #print "cluster so far -> "
                #print self.clusterHash[clusterCount]

                self.updataMinimumValue()
                # Xx DRAWING xX

                #print self.DrawingArray



            # ======================================= WHEN TWO CLUSTERS ARE CLOSE ==============================================================

            else:
                # iterate over all the cluster hash values and see whether the two next closest structures are already in one of the clusters
                for ki in self.clusterHash.keys():

                # ---------------------------------------- TWO BIG CLUSTERS MERGER --------------------------------------------------------------

                    # if both the index in the sorted hash belong to a cluster, merge the two
                    if ( ( (self.firstMinimumIndex in self.centroidArray) and (self.secondMinimumIndex in self.centroidArray) ) ):
#                        print "Two Big Clusters found at proximity! proceeding to merge"

                        topRMSDBeforeRemoval = self.distanceMatrixDummy.min()

                        # getting the cluster numbers
                        for cluster1No in self.clusterHash.keys():
                            if self.firstMinimumIndex == self.clusterHash[cluster1No][0]:
                                break
                        for cluster2No in self.clusterHash.keys():
                            if self.secondMinimumIndex == self.clusterHash[cluster2No][0]:
                                break

                        if True:

                            #print "Two clusters to be merged are "+str(cluster1No)+" with centroid "+str(self.clusterHash[cluster1No][0])+" and "+str(cluster2No)+" with centroid "+str(self.clusterHash[cluster2No][0])

                            # getting elements from first cluster and transfer it to next one
                            # THE FIRST CLUSTER IS THE ONE TO ADD THE INDEXES TO!!
                            for ClusIdx2 in self.clusterHash[cluster2No][1]:
                                self.clusterHash[cluster1No][1] += [ClusIdx2]


                            # ----------------------------- ASSIGN A CLUSTER

                            self.clusterHash[cluster1No][0] = self.getCentroid(self.clusterHash[cluster1No][1], self.clusterHash[cluster1No][0])

                            # ----------------------------- DELETING OLD CENTROID OF OTHER CLUSTER FROM LIST

                            #print ">> Deleting old cluster indexes from the Dummy matrix"

                            if self.clusterHash[cluster1No][0] != self.clusterHash[cluster2No][0]:

                                for ind in self.clusterHash[cluster2No][1]:

#                                    print "index "+str(ind)+" in dummy matrix deleted (replaced by 100s)"
                                    if ind != self.clusterHash[cluster1No][0]:
                                        self.distanceMatrixDummy[ind, :] = 100
                                        self.distanceMatrixDummy[:,ind] = 100


                                #print ">> Deleting centroid of old cluster: "+str(self.clusterHash[cluster2No][0])

                                # -------------------------- ADDING OLD CENTROID VALUE TO LIST OF INDEXES

                                self.clusteredIndex += [self.clusterHash[cluster2No][0]]

                                # --------------------------- REMOVING OLD CENTROID FROM CENTROID ARRAY

                                self.centroidArray.pop(self.centroidArray.index(self.clusterHash[cluster2No][0]))
                                #print "centroid array: "+str(self.centroidArray)

                                # --------------------------- ADDING NEW CENTROID TO CENTROID ARRAY

                                if self.clusterHash[cluster1No][0] not in self.centroidArray:
                                    self.centroidArray += [self.clusterHash[cluster1No][0]]

                            # --------------------------- DELETING SECOND CLUSTER

                            del self.clusterHash[cluster2No]
                            #print "deleted cluster number "+str(cluster2No)+" permanently!"


                            #self.CheckIfClusterFull(cluster1No, self.clusterHash[cluster1No][0], self.clusterHash[cluster1No][1])

                            ClusterMergingSwitch = True

                            # =================================== Xx DRAWING xX ===========================================

                            # those clusters are TWO BIG ONES SO THREAD CAREFULLY  #TODO: try see if height difference matters
                            # locate first cluster:
                            for arraylist1 in self.DrawingArray:
                                if arraylist1[0] == cluster1No:
                                    break
                            for arraylist2 in self.DrawingArray:
                                if arraylist2[0] == cluster2No:
                                    break
                            # add all the content of arraylist 2 into arraylist 1
                            arraylist1copy = copy.deepcopy(arraylist1[1])
                            arraylist2copy = copy.deepcopy(arraylist2[1])

                            # dealing with size now!

                            if len(arraylist1copy) > len(arraylist2copy):
                                indexAway = 0
                                for elem in arraylist2[1]:
                                    arraylist1[1] += [elem]
                                    indexAway += 1

                                # merge the two arrays
#                                arraylist1[1] += [self.joinTwoArrays(arraylist1copy[-1] ,arraylist2copy[-1],topRMSDBeforeRemoval )]
                                arraylist1[1].append( [ "self.joinTwoArrays",[ indexAway, topRMSDBeforeRemoval ] ])

                                # destroy the scond array permanently:
                                self.DrawingArray.pop(self.DrawingArray.index(arraylist2))

                            # in case the second one is bigger than the first one just reverse the addition of elements
                            else:
                                indexAway = 0
                                for elem in arraylist1[1]:
                                    arraylist2[1] += [elem]
                                    indexAway += 1
                                # merge the two arrays
#                                arraylist2[1] += [self.joinTwoArrays(arraylist2copy[-1] ,arraylist1copy[-1],topRMSDBeforeRemoval )]
                                arraylist2[1].append( [ "self.joinTwoArrays",[ indexAway, topRMSDBeforeRemoval ] ])

                                # destroy the first array permanently:
                                self.DrawingArray.pop(self.DrawingArray.index(arraylist1))
                                # rename the number label of the second array so that the others can still merge with it
                                arraylist2[0] = cluster1No



                            #print "MERGED 2 BIG CLUSTERS AND DELETED ELEMENT "+str(cluster2No)
#                            print self.DrawingArray

                            # save the centroids of that RMSD
                            self.SaveCentroidsAtRMSD(topRMSDBeforeRemoval)

                            self.updataMinimumValue()

                            break

                # ---------------------------------------- ONE SMALL AND ONE BIG CLUSTER MERGER -------------------------------------------------

                    elif ( ((self.firstMinimumIndex in self.centroidArray) and (self.secondMinimumIndex in self.clusterHash[ki][1]) ) \
                           or ( (self.secondMinimumIndex in self.centroidArray) and (self.firstMinimumIndex in self.clusterHash[ki][1]) ) ):
#                        print "Bigger cluster will absord smaller one"

                        # getting which index is centroid
                        if self.firstMinimumIndex in self.centroidArray:
                            centroidIndexForMerge = self.firstMinimumIndex # -> 21
                        else:
                            centroidIndexForMerge = self.secondMinimumIndex

                        topRMSDBeforeRemoval = self.distanceMatrixDummy.min()

                        # getting all the indexes from one cluster to another
                        # locate the cluster having for centroid  centroidIndexForMerge
                        for clusterKey in self.clusterHash.keys():
                            if centroidIndexForMerge == self.clusterHash[clusterKey][0]:
                                #self.checkWithinBetween(self.clusterHash[clusterKey], self.clusterHash[ki], topRMSDBeforeRemoval )
                                for ind in self.clusterHash[ki][1]:
                                    self.clusterHash[clusterKey][1] += [ind]

                                break



                        # ---------------- CREATING CENTROID

                        # create a new centroid index if there is more than 3 indexes in the cluster array
                        self.clusterHash[clusterKey][0] = self.getCentroid(self.clusterHash[clusterKey][1], self.clusterHash[clusterKey][0])
                        #print "Centroid assigned for newly merged cluster "+str(clusterKey)+" -> "+str(self.clusterHash[clusterKey][0])

                        # add the new centroid value to the centroid array:
                        if (self.clusterHash[clusterKey][0] not in self.centroidArray):
                            self.centroidArray += [self.clusterHash[clusterKey][0]]

                        # delete the centroid index from the self.clusteredIndex
                        if self.clusterHash[clusterKey][0] in self.clusteredIndex:
                            self.clusteredIndex.pop(self.clusteredIndex.index(self.clusterHash[clusterKey][0]))

                        # ---------------- REMOVING THE NON CENTROID ELEMENTS OF NEW CLUSTER FROM LIST

                        #print ">> deleting elements "+str(nonCentroids)+" from the RMSD lists "

                        for indexValues in self.clusterHash[ki][1]:
#                            print "index "+str(indexValues)+" in dummy matrix deleted (replaced by 100s)"
                            self.distanceMatrixDummy[indexValues, :] = 100
                            self.distanceMatrixDummy[:,indexValues] = 100


                        # ---------------- DELETING THE MERGED CLUSTER

                        #print "deleting cluster: "+str(ki)
                        del self.clusterHash[ki]

                        #self.CheckIfClusterFull(clusterKey, self.clusterHash[clusterKey][0], self.clusterHash[clusterKey][1])

#                        self.PrintSortedListState()


                        ClusterMergingSwitch = True

                        # Xx DRAWING xX
                        # those clusters are small and of the same height normally #TODO: try see if height difference matters
                        # locate first cluster:
                        for arraylist1 in self.DrawingArray:
                            if arraylist1[0] == clusterKey:
                                break
                        for arraylist2 in self.DrawingArray:
                            if arraylist2[0] == ki:
                                break


                        indexAway = 0
                        for elem in arraylist2[1]:
                            arraylist1[1] += [elem]
                            indexAway += 1

                        # merge the two arrays
#                        arraylist1[1] += [self.joinTwoArrays(arraylist1copy[-1] ,arraylist2copy[-1],topRMSDBeforeRemoval )]
                        arraylist1[1].append( ["self.joinTwoArrays", [indexAway, topRMSDBeforeRemoval ]] )

                        # destroy the scond array permanently:
                        self.DrawingArray.pop(self.DrawingArray.index(arraylist2))
                        #print "MERGED 1 SMALL AND 1 BIG CLUSTERS AND DELETED ELEMENT "
#                        print self.DrawingArray

                        # save the centroids of that RMSD
                        self.SaveCentroidsAtRMSD(topRMSDBeforeRemoval)
                        self.updataMinimumValue()

                        break

                # --------------------------------------------- MERGING 2 SMALL CLUSTERS ---------------------------------------

                    elif ((self.firstMinimumIndex in self.clusteredIndex) and (self.secondMinimumIndex in self.clusteredIndex)):

#                        print "merging two small clusters !!"

                        topRMSDBeforeRemoval = self.distanceMatrixDummy.min()

                        #getting the cluster numbers
                        for cluster1No in self.clusterHash.keys():
                            if self.firstMinimumIndex in self.clusterHash[cluster1No][1]:
                                break
                        for cluster2No in self.clusterHash.keys():
                            if self.secondMinimumIndex in self.clusterHash[cluster2No][1]:
                                break

#                        print "Two clusters to be merged are "+str(cluster1No)+" with centroid "+str(self.clusterHash[cluster1No][0])+" and "+str(cluster2No)+" with centroid "+str(self.clusterHash[cluster2No][0])

                        # getting elements from first cluster and transfer it to next one
                        # THE FIRST CLUSTER IS THE ONE TO ADD THE INDEXES TO!!
                        for ClusIdx2 in self.clusterHash[cluster2No][1]:
                            self.clusterHash[cluster1No][1] += [ClusIdx2]

                        # ----------------------------- ASSIGN A CLUSTER

                        self.clusterHash[cluster1No][0] = self.getCentroid(self.clusterHash[cluster1No][1], self.clusterHash[cluster1No][0])

                        # ----------------------------- DELETING FROM SORTED LIST ALL THE INDEXES 1ST CLUSTER ARRAY

                        # get the non centroid elements:
                        nonCentroids = copy.deepcopy(self.clusterHash[cluster1No][1])
                        nonCentroids.pop(nonCentroids.index(self.clusterHash[cluster1No][0]))
                        #print ">> deleting elements "+str(nonCentroids)+" from the RMSD lists, not get centroid function "

                        for indexValues in nonCentroids:
#                            print "index "+str(indexValues)+" in dummy matrix deleted (replaced by 100s)"
                            self.distanceMatrixDummy[indexValues, :] = 100
                            #print self.distanceMatrixDummy[indexValues]
                            self.distanceMatrixDummy[:,indexValues] = 100
                            #print self.distanceMatrixDummy[:,indexValues]


                        # ----------------------------- DELETE NEW CENTROID FROM INDEX ARRAY

                        if self.clusterHash[cluster1No][0] in self.clusteredIndex:
                            self.clusteredIndex.pop(self.clusteredIndex.index(self.clusterHash[cluster1No][0]))

                            #print "deleted "+str(self.clusterHash[cluster1No][0])+" from centroid array"

                        # ----------------------------- ADDING NEW CENTROID TO CENTROID ARRAY

                        if self.clusterHash[cluster1No][0] not in self.centroidArray:
                            self.centroidArray += [self.clusterHash[cluster1No][0]]

                        # --------------------------- DELETING SECOND CLUSTER

                        del self.clusterHash[cluster2No]
                        #print "deleted cluster number "+str(cluster2No)+" permanently!"

                        ClusterMergingSwitch = True

                        # Xx DRAWING xX
                        # those clusters are small and of the same height normally
                        # locate first cluster:
                        for arraylist1 in self.DrawingArray:
                            if arraylist1[0] == cluster1No:
                                break
                        for arraylist2 in self.DrawingArray:
                            if arraylist2[0] == cluster2No:
                                break

                        indexAway = 0 # counts the number of indexes away from arraylist1 arraylist2 is
                        for elem in arraylist2[1]:
                            arraylist1[1].append(elem)
                            indexAway += 1
                        # merge the two arrays
#                        arraylist1[1] += [self.joinTwoArrays(arraylist1copy[-1] ,arraylist2copy[-1], topRMSDBeforeRemoval )]

                        arraylist1[1] += [ ["self.joinTwoArrays", [ indexAway , topRMSDBeforeRemoval] ] ]

                        # destroy the scond array permanently:
                        self.DrawingArray.pop(self.DrawingArray.index(arraylist2))
#                        print "MERGED 2 SMALL CLUSTERS " #+str(cluster2No)
                        #print self.DrawingArray

                        # save the centroids of that RMSD
                        self.SaveCentroidsAtRMSD(topRMSDBeforeRemoval)
                        self.updataMinimumValue()

                        break
            # ======================================= WHEN TWO CLUSTERS ARE FAR AWAY BUT SINGLE INDEXES ARE CLOSE =====================================

            if ClusterMergingSwitch == False:

                for i in self.clusterHash.keys():

                    # if one of the structure index is already present in a cluster, add the other one
                    if ((self.firstMinimumIndex in self.clusterHash[i][1]) or (self.secondMinimumIndex in self.clusterHash[i][1])):
                        # adding either of the numbers in the array of the found cluster
                        if (self.firstMinimumIndex in self.clusterHash[i][1]):
                            self.clusterHash[i][1] += [self.secondMinimumIndex] # adding the index to appropriate cluster
                            self.clusteredIndex += [self.secondMinimumIndex] # adding the index to the list of added indexes
                            NewAddedIndex = self.secondMinimumIndex

                        elif (self.secondMinimumIndex in self.clusterHash[i][1]):
                            self.clusterHash[i][1] += [self.firstMinimumIndex] # adding the index to appropriate cluster
                            self.clusteredIndex +=  [self.firstMinimumIndex]  # adding the index to the list of added indexes
                            NewAddedIndex = self.firstMinimumIndex

#                        print "new index to be added to cluster "+str(i)+" : "+str(NewAddedIndex)

                        # Xx DRAWING xX
                        for arraylist in self.DrawingArray:
                            if arraylist[0] == i:
                                self.iterant += 1
#                                arraylist[1] += [self.addSIngleIndex (arraylist[1][-1], self.distanceMatrixDummy.min()) ]
                                arraylist[1] += [ ["addSingleIndex", "x"+str(self.iterant) ,self.distanceMatrixDummy.min()] ]
                                break

                        # save the centroids of that RMSD
                        self.SaveCentroidsAtRMSD(self.distanceMatrixDummy.min())

                        # remove value in the Matrix
                        self.distanceMatrixDummy[self.minimumIndex] = 100

                        # switch back the found in structure boolean
                        foundStructureInCluster = True

                        # -------------------------------- ASSIGNING NEW CENTROID

                        # create a new centroid index if there is more than 3 indexes in the cluster array
                        self.clusterHash[i][0] = self.getCentroid(self.clusterHash[i][1], self.clusterHash[i][0])
                        #print "Centroid assigned for cluster "+str(i)+" -> "+str(self.clusterHash[i][0])

                        # add the new centroid value to the centroid array:
                        if (self.clusterHash[i][0] not in self.centroidArray):
                            self.centroidArray += [self.clusterHash[i][0]]

                        # delete the centroid index from the self.clusteredIndex
                        if self.clusterHash[i][0] in self.clusteredIndex:
                            self.clusteredIndex.pop(self.clusteredIndex.index(self.clusterHash[i][0]))


                        # ---------------------------- REMOVING NON CENTROID PARTS
                        if len(self.clusterHash[i][1]) == 3:
                            # get the non centroid elements:
                            nonCentroids = copy.deepcopy(self.clusterHash[i][1])
                            nonCentroids.pop(nonCentroids.index(self.clusterHash[i][0]))
                            #print ">> deleting elements "+str(nonCentroids)+" from the RMSD lists "

                            for indexValues in nonCentroids:
#                                print "index "+str(indexValues)+" in dummy matrix deleted (replaced by 100s)"
                                self.distanceMatrixDummy[indexValues, :] = 100
                                self.distanceMatrixDummy[:,indexValues] = 100

                        elif len(self.clusterHash[i][1]) > 3:
                            # get the non centroid elements:
                            nonCentroids = copy.deepcopy(self.clusterHash[i][1])
                            nonCentroids.pop(nonCentroids.index(self.clusterHash[i][0]))

                            #print ">> deleting element "+str(NewAddedIndex)+" from the RMSD lists "

                            for indexValues in nonCentroids:
#                                print "index "+str(indexValues)+" in dummy matrix deleted (replaced by 100s)"
                                self.distanceMatrixDummy[indexValues, :] = 100
                                self.distanceMatrixDummy[:,indexValues] = 100


                        self.updataMinimumValue()

                        break

                # ---------------------------------- CREATING NEW CLUSTERS ---------------------------------

                if foundStructureInCluster == False:
                    clusterCount += 1
                    # create a new cluster with the following data
                    self.clusterHash[clusterCount] = [None , [ self.firstMinimumIndex , self.secondMinimumIndex ]]

                    # add the new elements to the list of already added elements
                    for newElem in self.clusterHash[clusterCount][1]:
                        self.clusteredIndex += [newElem]

                    self.DrawingArray.append( [clusterCount, [ [ "ini" , self.makeSingleArray(self.distanceMatrixDummy.min()) ] ]] )

                    # remove value in the Matrix
                    self.distanceMatrixDummy[self.minimumIndex] = 100

                    self.updataMinimumValue()

                    # Xx DRAWING xX

                    #print self.DrawingArray

            # switch back the foundStructureInCluster boolean:
            foundStructureInCluster = False
            self.withinBetween = True
            self.clusteredIndex.sort()
            #print "indexes assiged to clusters so far:"+ str(self.clusteredIndex)
            print " > RMSDs remaining: "+str(len(self.distanceMatrixDummy[self.distanceMatrixDummy < 100]))+" | number of clusters: "+str(len(self.clusterHash.keys()))


            # restart the clusterSwitch
            ClusterMergingSwitch = False

            self.PrintClusterState()
#            self.PrintSortedListState()
            #print self.clusteredIndex

# ==================================================================================================================================================== #
# =========================================================== END OF WHILE LOOP ====================================================================== #
# ==================================================================================================================================================== #

        print ">>> PREPARING TO DRAW NOW"
#        print self.DrawingArray

        #print ">>> RMSD mutliplier: "+str(self.RMSDMutliplier)

        # first need to extract the x data into a dictio to make the arrays real
        print ">>> initialising Drawing Array"
        self.unSeparatedArray = []
#        counter = 0.00
        for clusters in self.DrawingArray:
            for clusterElement in clusters[1]:
                self.unSeparatedArray += [clusterElement]
#            counter += 1.00
#            print str(float(counter/ len(self.DrawingArray))*100.00)+" %"

        xdictio = {}
        xiterator = 900.0/float(len(self.distanceMatrix)) # len matrix because there as many x as indexes

#        secondCounter = 0.00

        print ">>> Creating coordinate array for drawing"
        for coordArray in self.unSeparatedArray:


            if coordArray[0] == "ini":
                xdictio[coordArray[1][0][0]] = xiterator
                xiterator += 900.0/float(len(self.distanceMatrix))
                xdictio[coordArray[1][1][0]] = xiterator
                xiterator += 900.0/float(len(self.distanceMatrix))
            if coordArray[0] == "addSingleIndex":
                xdictio[coordArray[1]] = xiterator
                xiterator += 900.0/float(len(self.distanceMatrix))

#            secondCounter += 1.00
#            print str(float(secondCounter/ len(unSeparatedArray))*100.00)+" %"

        self.arrayReadyforDrawing = []
        thirdCounter = 0.00

        print ">>> Converting virtual coordinates into actual ones "
        # MAKE THE UNSEPARATED DRAWING ARRAY REAL!!!!!


        index = 0
        rmsd = 0
        indexOfFirst = 0
        indexOfSecond= 0
        for bluePrints in self.unSeparatedArray:
#            self.arrayReadyforDrawing += [self.makeArrayReal( bluePrints , xdictio)]

            if bluePrints[0] == "ini":
                self.arrayReadyforDrawing += [self.makeArrayReal( bluePrints[1] , xdictio)]
#                print self.arrayReadyforDrawing

            if bluePrints[0] == "addSingleIndex" :
                # getting the rmsd
                if self.unSeparatedArray[index - 1][0] == "ini" or self.unSeparatedArray[index - 1][0] == "self.joinTwoArrays":
                    rmsd = self.unSeparatedArray[index - 1][-1][-1]
                else:
                    rmsd = self.unSeparatedArray[index - 1][-1]

                self.arrayReadyforDrawing += [self.addSingleIndex(bluePrints[1], bluePrints[2], xdictio, rmsd)]

            if bluePrints[0] == "self.joinTwoArrays":
                # get the unseparated array index of the 2nd array to join:
                indexOfFirst = self.unSeparatedArray.index(bluePrints)
                indexOfSecond = indexOfFirst - bluePrints[1][0] - 1


                # getting the rmsd
                if self.unSeparatedArray[index - 1][0] == "ini" or self.unSeparatedArray[index - 1][0] == "self.joinTwoArrays":
                    rmsd = self.unSeparatedArray[index - 1][-1][-1]
                else:
                    rmsd = self.unSeparatedArray[index - 1][-1]


                self.arrayReadyforDrawing += [self.joinTwoArrays(indexOfSecond ,bluePrints[1][1],rmsd)]

            index += 1
            thirdCounter += 1.00
            percentage = float(thirdCounter/len(self.unSeparatedArray))*100.00
            if percentage == 20 or percentage == 40 or percentage == 60 or percentage == 80 or (float(thirdCounter/len(self.unSeparatedArray))*100.00) == 100:
                print ">>> Drawing %.2f"%(percentage)+ " %"
        #print arrayReadyforDrawing
        print ">>> Done!"


    # -------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                    FUNCTIONS FOR CLUSTER CREATION
    # -------------------------------------------------------------------------------------------------------------------------------------------------

    def updataMinimumValue(self):
        if self.distanceMatrixDummy.min() != 100:
            self.minimumIndex = np.where(self.distanceMatrixDummy == self.distanceMatrixDummy.min())
#            print self.distanceMatrixDummy.min()
#            print "error ?: "+str(self.minimumIndex[0])
            self.firstMinimumIndex = int(self.minimumIndex[0][0])
            self.secondMinimumIndex = int(self.minimumIndex[1][0])

            # saving the last RMSD value so as to compute the dendogram tree according to that height
            self.maximalRMSD = self.distanceMatrixDummy.min()
            self.RMSDMutliplier = (self.Parent.treePanel.GetSize()[1] - 50)/(self.distanceMatrixDummy.min())


    def getCentroid (self, array, oldCentroid):
        RMSDHash = {} # contains the average RMSD as key and index as value
        RMSDArray = [] # to hold the different RMSD to calculate mean value
        newCentroid = 0

        # first you need to get all rmsds relative of one to another
        for i in array: # -> the array contains the indexes of the proteins: 0,1,3, ...
            for j in array:
                if i == j:
                    pass # as it does not make sense to look for the RMSD value for self.distanceMatrix[0][0]

                elif self.distanceMatrix[i][j] == 100:
                    #print "RMSDArray of distance: "+str(i)+" "+str(j)+" -> "+str(self.distanceMatrix[j][i])
                    RMSDArray += [self.distanceMatrix[j][i]]
                    #print RMSDArray

                else:
                    #print "RMSDArray of distance: "+str(i)+" "+str(j)+" -> "+str(self.distanceMatrix[i][j])
                    RMSDArray += [self.distanceMatrix[i][j]]
                    #print RMSDArray

            # get the average value of the RMSDARRAY
            RMSDHash[sum(RMSDArray)/len(RMSDArray)] = i
#            print "RMSDHash: "
#            print RMSDHash
            # empty array
            RMSDArray = []

        newCentroid = RMSDHash[np.array(RMSDHash.keys())[np.array(RMSDHash.keys())<101].min()]


        # ------------------------ DELETING OLD CENTROID --------------------------------------

        # check whether the old centroid is the same as the new, in case it is not, delete the data of the old and put in the one of the new
        if ((newCentroid != oldCentroid) and (oldCentroid != None)):


            # -------------------------- DELETING OLD CENTROID FROM CENTROID ARRAY

            if oldCentroid in self.centroidArray:
                self.centroidArray.pop(self.centroidArray.index(oldCentroid))

            # -------------------------- ADDING OLD CENTROID VALUE TO LIST OF INDEXES ----------------------------

            self.clusteredIndex += [oldCentroid]

            # -------------------------- MAKING NEW CENTROID VALUE FREE TO INTERACT ----------------------------------------------

            self.distanceMatrixDummy[newCentroid] = copy.deepcopy(self.distanceMatrix[newCentroid])
            self.distanceMatrixDummy[:,newCentroid] = copy.deepcopy(self.distanceMatrix[:,newCentroid])

            # -------------------------- REPLACING INDEX VALUES SO THAT CAN INTERACT WITH NEW CENTROID -------------------------------

            for index in range(len(self.distanceMatrix)):
                if ((index != newCentroid) and (index not in self.clusteredIndex)):
#                    print "recovering values of indexes "+str(index)+" in Dummy Matrix "
                    self.distanceMatrixDummy[index] = copy.deepcopy(self.distanceMatrix[index])
                    self.distanceMatrixDummy[:,index] = copy.deepcopy(self.distanceMatrix[:,index])

            # -------------------------- IF PRESENT, DELETE THE NEW CENTROID FROM INDEX LIST ---------------------

            self.distanceMatrixDummy[oldCentroid, :] = 100
            self.distanceMatrixDummy[:,oldCentroid] = 100

#             !!! ---------------------- RE DELETE THE VERTICAL RMSDS OF INDEXES ONLY OF CLUSTERS CONTAINING CENTROIDS
            for array in self.clusterHash.values():
                if array[0] != None:
                    for indes in array[1]:
                        if (indes != array[0]) and (indes != newCentroid):
#                            print "re deleting vertical values of indexes: "+str(indes)
                            self.distanceMatrixDummy[indes, :] = 100
                            self.distanceMatrixDummy[:, indes] = 100

            # delete the centroid values from the clustered index
            if newCentroid in self.clusteredIndex:
                self.clusteredIndex.pop(self.clusteredIndex.index(newCentroid))


        self.updataMinimumValue()

        # return the Centroid
        return newCentroid # 21


    def PrintClusterState(self):
        #print " ---------------------------- CLUSTERS SO FAR ---------------------------------------"
        #for key in self.clusterHash.keys():
            #print "cluster "+str(key)+" "+str(self.clusterHash[key])
        #print " ------------------------------------------------------------------------------------"
        pass

    def PrintSortedListState(self):
        print "----------------------------- STATE OF SORTED ARRAY"
        if len(self.sortedHashKeys) > 10:
            for v in range(10):
                print str(self.unSortedHash[self.sortedHashKeys[v]])+" "+str(self.sortedHashKeys[v])
        else:
            for v in range(len(self.sortedHashKeys)):
                print str(self.unSortedHash[self.sortedHashKeys[v]])+" "+str(self.sortedHashKeys[v])
        print "---------------------------------------------------"

    def CheckIfClusterFull(self, clusterNumber, centroidIndex, clusterArray):
        '''This is used to check whether a cluster has exceeded its within RMSD value, in this case it cannot grow any more'''

        print "checking whether cluster has exceeded maximal RMSD value of cluster "+str(clusterNumber)+" centroid "+str(centroidIndex)+" "+str(clusterArray)
        distanceArray = [] # this will hold the distances between the centroids and the indexes
        ArrayedRMSDIndexes = []
        valuesToBeRemoved = []
        indexArray = []
        indexArray = copy.deepcopy(clusterArray)
        # removing the centroid from the indexArray
        indexArray.pop(indexArray.index(centroidIndex))

        # calculate the distance between the centroids and the distances:
        for i in range(len(indexArray)):
            if self.distanceMatrix[centroidIndex][indexArray[i]] == 100:
                distanceArray += [self.distanceMatrix[indexArray[i]][centroidIndex]]
            else:
                distanceArray += [self.distanceMatrix[centroidIndex][indexArray[i]]]

        # get the maximal element of the distance array and compare it with RMSD treshold
        if max(distanceArray) >= self.maximalRMSD:
            print "maximal RMSD in cluster "+str(clusterNumber)+" exceeded ("+ str(max(distanceArray)) +") proceeding to isolate cluster"
            print "removing centroid "+ str(centroidIndex) +" from RMSD sorted list"

            for RMSDval in self.sortedHashKeys:

                #print self.unSortedHash[RMSDval]
                ArrayedRMSDIndexes = self.unSortedHash[RMSDval]

                if (centroidIndex in ArrayedRMSDIndexes):
                    #print "found"+str(nonCentroids[0])+" "+str(nonCentroids[1])+" in "+str(ArrayedRMSDIndexes)
                    valuesToBeRemoved += [RMSDval]


            for valu in valuesToBeRemoved:
                del self.unSortedHash[valu]
                self.sortedHashKeys.pop(self.sortedHashKeys.index(valu))

            self.PrintSortedListState()

    def checkWithinBetween(self, array1, array2, betweenValue):
        '''This function is used to check that the within value of two merged clusters is always smaller than the between value of 2 clusters'''
        # get the centroids of both arrays:

        # ---------------------------------------- CASE OF 2 BIG CLUSTERS --------------------------
        if (len(array1[1]) > 2) and (len(array2[1]) > 2):

            centroid1 = array1[0]
            centroid2 = array2[0]
            indexArray1 = []
            indexArray2 = []

            # get all the rmsd values between the centroid and the indexes in the clusters:
            for value1 in array1[1]:
                if (self.distanceMatrix[centroid1][value1] == 100) and (value1 != centroid1):
                    indexArray1 += [self.distanceMatrix[value1][centroid1]]
                elif (self.distanceMatrix[centroid1][value1] != 100) and (value1 != centroid1):
                    indexArray1 += [self.distanceMatrix[centroid1][value1]]

            for value2 in array2[1]:
                if (self.distanceMatrix[centroid2][value2] == 100) and (value2 != centroid2):
                    indexArray2 += [self.distanceMatrix[value2][centroid2]]
                elif (self.distanceMatrix[centroid2][value2] != 100) and (value2 != centroid2):
                    indexArray2 += [self.distanceMatrix[centroid2][value2]]

            # get the maximal in between value of each cluster:
            max1 = max(indexArray1)
            max2 = max(indexArray2)

            # check those values against the in-between value:
            if max1 > betweenValue :
                self.withinBetween = False
                print "YOUR WITHIN VALUES BIGGER WHEN MERGING TWO BIG ARRAYS..."
                print "\n"
                print "max1 "+str(max1)+"  "+str(indexArray1)
                print "between value "+ str(betweenValue)
                #self.createNewCentroidWithLowerWithin(array1[1], centroid1, betweenValue)

            if max2 > betweenValue:
                self.withinBetween = False
                print "YOUR WITHIN VALUES BIGGER WHEN MERGING TWO BIG ARRAYS..."
                print "\n"
                print "max2 "+str(max2)+"  "+str(indexArray2)
                print "between value "+ str(betweenValue)
                #self.createNewCentroidWithLowerWithin(array2[1], centroid2, betweenValue)


        # ------------------------------------- CASE OF 1 BIG AND 1 SMALL ------------------------------------------------------------------

        elif (len(array1[1]) > 2) and (len(array2[1]) == 2):
            centroid1 = array1[0]
            indexArray1 = []
            indexArray2 = []
            indx1 = array2[1][0]
            indx2 = array2[1][1]

            # get all the rmsd values between the centroid and the indexes in the clusters:
            for value1 in array1[1]:
                if (self.distanceMatrix[centroid1][value1] == 100) and (value1 != centroid1):
                    indexArray1 += [self.distanceMatrix[value1][centroid1]]
                elif (self.distanceMatrix[centroid1][value1] != 100) and (value1 != centroid1):
                    indexArray1 += [self.distanceMatrix[centroid1][value1]]


            if self.distanceMatrix[indx1][indx2]  == 100:
                indexArray2 += [self.distanceMatrix[indx2][indx1] ]
            else:
                indexArray2 += [self.distanceMatrix[indx1][indx2] ]

            # get the maximal in between value of each cluster:
            max1 = max(indexArray1)
            max2 = indexArray2[0]

            # check those values against the in-between value:
            if max1 > betweenValue or max2 > betweenValue:
                self.withinBetween = False
                print "YOUR WITHIN VALUES ARE HIGHER THAN YOUR BETWEEN WHEN MERGING ONE BIG AND ONE SMALL..."
                print "\n"
                print "max1 "+str(max1)+"  "+str(indexArray1)
                print "max2 "+str(max2)+"  "+str(indexArray2)
                print "between value "+ str(betweenValue)

                # empty the self.sortedHashKeys and finish everythin
                #self.sortedHashKeys = []



    def getMaxRMSD(self):
        return self.RMSDMutliplier




    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                           FUNCTIONS FOR DRAWING DENDOGRAM
    # ----------------------------------------------------------------------------------------------------------------------------------------------------


    def makeSingleArray(self, RMSD):
        self.iterant += 1
        x1 = "x"+str(self.iterant)
        first = [ x1, "0.0"]

        self.iterant += 1
        x2 = "x"+str(self.iterant)
        second = [x2, "0.0"]

        height = str(RMSD)+"* self.RMSDMutliplier"

        third = ["((( "+str(x2)+" - "+str(x1)+" )/2.0)+"+str(x1)+")", str(first[1])+"+"+height]

        self.SaveCentroidsAtRMSD(RMSD)

        return [first, second, height, third, RMSD]

    def addSingleIndex (self, x, newRMSD, dictio, oldRMSD):
#        for xs in dictio.keys():
#            exec str(xs)+"="+str(dictio[x])

        x1 = self.arrayReadyforDrawing[-1][3][0]
        y1 = self.arrayReadyforDrawing[-1][3][1]
        first = [x1, y1]

        x2 = dictio[x]
        second = [x2, 0.0]

        height = (newRMSD - oldRMSD) * self.RMSDMutliplier

        third = [((x2-x1)/2.0)+x1, y1+height ]


        return [first, second, height, third]

    def joinTwoArrays (self, index2 ,newRMSD, oldRMSD):
        x1 = self.arrayReadyforDrawing[-1][3][0]
        y1 = self.arrayReadyforDrawing[-1][3][1]
        first = [x1, y1]

        x2 = self.arrayReadyforDrawing[index2][3][0]
        y2 = self.arrayReadyforDrawing[index2][3][1]
        second = [x2, y2]

        height = (newRMSD - oldRMSD) * self.RMSDMutliplier

        third = [((x2-x1)/2.0)+x1, y1+height ]

        # save the centroids of that RMSD
#        self.SaveCentroidsAtRMSD(RMSD)

        return [first, second, height, third]


    def makeArrayReal (self, array, dictio):
        for x in dictio.keys():
            exec str(x)+"="+str(dictio[x])

        exec "RMSDMutliplier = self.RMSDMutliplier"

        x1 = eval (array[0][0])
        y1 = eval(array[0][1])

        x2 = eval (array[1][0])
        y2 = eval(array[1][1])

        h = eval(array[2])

        x3 = eval (array[3][0])
        y3 = eval(array[3][1])

        return [[x1,y1],[x2,y2],h,[x3,y3]]



    def SaveCentroidsAtRMSD (self, rmsd):
        centroidArr = copy.deepcopy(self.centroidArray)

        # get the within cluster RSMD value, which is average value of centroid to all
        RMSD_from_centroid = []
        RMSDs_to_be_averaged = []
        array_of_indexes = []

        for centroid in centroidArr:
            RMSDs_to_be_averaged = []

            # first get the array of indexes for each centroid:
            for value in self.clusterHash.values():
                if centroid == value[0]:
                    array_of_indexes = copy.deepcopy(value[1])
                    break

            for index in array_of_indexes:
                if index == centroid:
                    pass
                elif self.distanceMatrix[centroid][index] == 100:
                    RMSDs_to_be_averaged.append(self.distanceMatrix[index][centroid])
                else:
                    RMSDs_to_be_averaged.append(self.distanceMatrix[centroid][index])


            RMSD_from_centroid.append( sum(RMSDs_to_be_averaged) / len(RMSDs_to_be_averaged) )


        self.RMSDThresholdHash[rmsd] = [centroidArr, RMSD_from_centroid]



# --------------------------------------------------------------- TREE PANEL

class treeDisplay (wx.Panel):

    def __init__(self, parent):
        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER, size=(1000,650))
        self.SetBackgroundColour("white")

    def plotFigure (self, array):
        dc = wx.PaintDC(self)
        dc.Clear()
        dc.SetPen(wx.Pen(wx.BLACK, 1))
        dc.SetDeviceOrigin(90, 620)
        dc.SetAxisOrientation(True, True)
#
#        dc.DrawRotatedText("mine", 40, 40, -90)
#        dc.DrawLine(25, 50, 250, 500) # (x1, y1, x2, y2)

        # creating the axis
        dc.DrawLine(-10, 0, -10, self.GetSize()[1]-50)
        # divide the RMSD axis in lines of 20

        height = self.GetSize()[1] - 50

        steps = 20
        graduationArray = []
        for grad in range(steps):
            graduationArray += [ ((height)/steps)*(grad) ] # + 1 because grad starts from 0

        # draw the graduations:
        for i in range(len(graduationArray)):
            dc.DrawLine(-10, int(graduationArray[i]), -15, int(graduationArray[i]))
            dc.DrawText("%.2f"%((self.Parent.RMSDPanel.maximalRMSD/steps)*i), -65, int(graduationArray[i])+10)

        # draw the last maximal gradutation
        dc.DrawLine(-10, self.GetSize()[1]-50, -15, self.GetSize()[1]-50)
        dc.DrawText("%.2f"%(self.Parent.RMSDPanel.maximalRMSD), -65, self.GetSize()[1]-40)

        #drawing in steps for each cluster creation/merging
        for element in array:
            # first vertical
            dc.DrawLine( element[0][0], element[0][1], element[0][0], (element[0][1] + element[2]))
            # then right
            dc.DrawLine (element[0][0], (element[0][1] + element[2]), element[1][0], (element[0][1] + element[2]) )
            # then down
            dc.DrawLine (element[1][0], (element[0][1] + element[2]), element[1][0], element[1][1] )

    def plotRMSDLine (self, RMSD):
        '''Plots the red line and all the centroid of clusters located directly beneath the line should be export as pdb bfiles'''

        # clear the screen and redraw everything:
        self.plotFigure(self.Parent.RMSDPanel.arrayReadyforDrawing)

        dc = wx.PaintDC(self)
        dc.SetPen(wx.Pen(wx.RED, 1))
        dc.SetDeviceOrigin(90, 620)
        dc.SetAxisOrientation(True, True)


        # transform the RMSD value into coordinates:
        widthScreen = self.GetSize()[0]

        multiplier = self.Parent.RMSDPanel.getMaxRMSD()
        #print str(multiplier*RMSD)
        #(self.distanceMatrix[self.distanceMatrix<100].max()

        dc.DrawLine(-9, int(multiplier*float(RMSD)), widthScreen -100,int(multiplier*float(RMSD)) )
