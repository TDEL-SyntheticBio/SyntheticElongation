from __future__ import division
from cc3d.core.PySteppables import *       
import math                      
import numpy                   
import sys
import random
RNG=random.SystemRandom()  #draw true random sequence, overkill but why not?
#Motility Variables
CtoM=52  #cell to adhesion value, self-consistent across simulations as 26 but x2 due to double interface formation
BASALY=28 #baseline motility for L929, due to their extremely motile behaviour, A, B
BASALB=40
BASALGR=31 #, A', B'' T0, change GR to 100 - see what happens 
SCF=0.5 #self-attenuator weighing basal motility vs loss of motility due to adhesion
#Self-Cutoff
ENDMCS=168010 #call runtime here directly
#Mitosis Variables
RADAVG=3 #average radius of the gaussian distribution to choose random radius
RADDEV=.425 #standard deviation of target radius, too low and division couples, too high and you'll lose cells at the start
MTFORCEMIN1=-3*10**(-4.362) #negative mitosis driving force fluctuation, usually only need to change the exponential part
MTFORCEMAX1=4*10**(-4.362)  #positive mitosis driving force fluctuation, usually change only the exponential part
MTFORCEMIN3=-3*10**(-4.541) #negative mitosis driving force fluctuation, usually only need to change the exponential part
MTFORCEMAX3=4*10**(-4.541)  #positive mitosis driving force fluctuation, usually change only the exponential part
MTFORCEMIN4=-3*10**(-4.673) #negative mitosis driving force fluctuation, usually only need to change the exponential part
MTFORCEMAX4=4*10**(-4.673)  #positive mitosis driving force fluctuation, usually change only the exponential part
#Signaling Variables
CONEXPSCF=10000 #Steady state expression of ligand expressed on a sender cell. This ligand is unaffected by signaling.
THETA=0 #time lag for expression of your constitutive, non-signaling affected ligand, start at 0 for simplicity, but can be adjusted depending on experiment results if known for generalizability
XI=1000 #controls how fast the sender cells reaches steady state for your constitutive, non-signaling affected ligand
FASTAPPROX=5000 #force approx for function of above variables at the time step, saves calling the mcs and doing the caluclation, purely computational speed effeciency

ALPHAYG=1 #controls how much your reporter synthesis magnitude due to signal S; can be set to 1 if decay is set properly
BETAYG=921.181 #threshold of signal required to generate a response in your cell due to signaling
EPSILONYG=526.389 #modulates how sharp the response is due to signaling, can turn synthesis to linear or heavi-side theta like if desired
KAPPAYG=25000 #general decay constant
THRESHOLDUPYG=17500 #activation threshold to change state, NOTE: line 38 was another set of parameters  
THRESHOLDDOYG=0 #deactivation threshold to revert state, as no clear deactivation/unsorting occured in reference experiments

ALPHABR=1 #controls how much your reporter synthesis magnitude due to signal S; can be set to 1 if decay is set properly
BETABR=921.181 #threshold of signal required to generate a response in your cell due to signaling
EPSILONBR=526.389 #modulates how sharp the response is due to signaling, can turn synthesis to linear or heavi-side theta like if desired
KAPPABR=25000 #general decay constant
THRESHOLDUPBR=17500 #activation threshold to change state, the choice of paramters in this paragraph render it similar to that of the previous paragraph
THRESHOLDDOBR=0 #deactivation threshold to revert state, as no clear deactivation/unsorting occured in reference experiments

#Single Cell Trace Variables
MARKEDCELLS=[11,112,247,277] # ID of cells to track if you desire single cell points tracked, change to fit setup

#Sampling and Comp Speed
RESOL=100 #Data sampling rate, choose to satisfy nyquist theorem if necessary
USEDNODES=8 #Choose a power of 2, otherwise the grids overlap and your simulation will eventually randomly crash, follow the recommendations given in the manual by developers

class ELUGMSteppable(SteppableBasePy):

    def __init__(self,_frequency=1):
        SteppableBasePy.__init__(self,_frequency)
    def start(self):

        self.pW1 = self.add_new_plot_window(
            title='Calibrate',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Psi',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True)                
        self.pW1.add_plot('IBR', style='Dots', color='green', size=3) 
        self.pW1.add_plot('IYG', style='Dots', color='red', size=3)
        self.pW1.add_plot('SBR', style='Dots', color='yellow', size=5)        
        self.pW1.add_plot('SYG', style='Dots', color='white', size=5)   #plots for morphospace/hetergenity measure
        
        self.pW2 = self.add_new_plot_window(
            title='Psi',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Psi',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True)                
        self.pW2.add_plot('Hamiltonian', style='Dots', color='white', size=3) #measure the system total energy
 
        self.pW3 = self.add_new_plot_window(
            title='Types',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Count',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True)                
        self.pW3.add_plot('Y', style='Dots', color='gray', size=3)
        self.pW3.add_plot('G', style='Dots', color='green', size=3)
        self.pW3.add_plot('B', style='Dots', color='blue', size=3)   
        self.pW3.add_plot('R', style='Dots', color='red', size=3)
        self.pW3.add_plot('AG', style='Dots', color='green', size=5)
        self.pW3.add_plot('AR', style='Dots', color='red', size=5) #count the types of cells in the simulation, along with how mny activate due to signaling
        
        self.pW4 = self.add_new_plot_window(
            title='Point System',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Count',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True)                
        self.pW4.add_plot('YG', style='Dots', color='green', size=3)
        self.pW4.add_plot('BR', style='Dots', color='blue', size=3) #measure average points per cell type, can be tied to average flurosence intensity
 
        self.pW5 = self.add_new_plot_window(
            title='Single Point System',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Count',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True)                
        self.pW5.add_plot('1', style='Dots', color='blue', size=3)
        self.pW5.add_plot('2', style='Dots', color='gray', size=3)
        self.pW5.add_plot('3', style='Dots', color='green', size=3)   
        self.pW5.add_plot('4', style='Dots', color='red', size=3)
        self.pW5.add_plot('1T', style='Dots', color='blue', size=3)
        self.pW5.add_plot('2T', style='Dots', color='gray', size=3)
        self.pW5.add_plot('3T', style='Dots', color='green', size=3)   
        self.pW5.add_plot('4T', style='Dots', color='red', size=3)  #single cell point traces 

        self.pW6 = self.add_new_plot_window(
            title='Sphericity',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Count',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True)                
        self.pW6.add_plot('BF', style='Dots', color='blue', size=3)
        self.pW6.add_plot('MR', style='Dots', color='green', size=3)  #bright field and merged field sphericity      
                
        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR #the adhesion matrix, call these values and store them for motility code
        YtoY=float(self.get_xml_element('YtoY').cdata)
        YtoG=GtoY=float(self.get_xml_element('YtoG').cdata)
        YtoB=BtoY=float(self.get_xml_element('YtoB').cdata)
        YtoR=RtoY=float(self.get_xml_element('YtoR').cdata)
        GtoG=float(self.get_xml_element('GtoG').cdata)
        GtoB=BtoG=float(self.get_xml_element('GtoB').cdata)
        GtoR=RtoG=float(self.get_xml_element('GtoR').cdata)
        BtoB=float(self.get_xml_element('BtoB').cdata)
        BtoR=RtoB=float(self.get_xml_element('BtoR').cdata)
        RtoR=float(self.get_xml_element('RtoR').cdata)              

        for cell in self.cellList:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #assign the cells a random target radius
            cell.lambdaSurface=2.5                    #temporary value, will be changed in later code
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2  #spherical surface area
            cell.lambdaVolume=2.5                     #temporary value, will be changed in later code
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            cell.dict["PTS"]=[0]                      #initial points for cell, feel free to randomize if desired
            cell.dict["P"]=[0,0]                      #activation counter, counts how many cells are active due to signaling at any given time

    def step(self,mcs):                 
        NUMTY=0 #number of type Y
        NUMTG=0 #number of type G
        NUMTB=0 #number of type B
        NUMTR=0 #number of type R
        
        YGPTS=0 #number of Y+G points
        BRPTS=0 #number of B+R points
        
        SYPSI=0 #system hamiltonian over all interaction over configuration

        CSAYGBR=0 #common surface area of YG to BR
        CSABRYG=0 #common surface area of BR to YG
        CSAYGYG=0 #common surface area of YG to YG
        CSABRBR=0 #common surface area of BR to BR
        YGCBR=0 #how many Y or G cells are in contact with B or R cell?
        BRCYG=0 #how many B or R cells are in contact with Y or G cells?
        YGCYG=0 #how many Y or G cells are in contact with Y or G cells?
        BRCBR=0 #how many B or R cells are in contact with B or R cells?
        
        SUMBFSF=0 #total bright field surface area
        SUMBFVL=0 #total bright field volume
        SUMMRSF=0 #total color field surface area
        SUMMRVL=0 #total color field volume
        
        NAR=0 #number of activated red cells due to signaling
        NAG=0 #number of activated green cells due to signaling
        
        if mcs==1:
            self.change_number_of_work_nodes(USEDNODES) #set to necessary computational nodes                  

        for cell in self.cellList: #iterate over cell list
            CSAY=0 #each cell detect how much sirface area it shares with Y cells
            CSAG=0 #each cell detect how much sirface area it shares with G cells
            CSAB=0 #each cell detect how much sirface area it shares with B cells
            CSAR=0 #each cell detect how much sirface area it shares with R cells
            CSAM=0 #each cell detect how much sirface area it shares with medium
            
            PTSY=0 #each cell gains points from neighbor type Y
            PTSG=0 #each cell gains points from neighbor type G
            PTSB=0 #each cell gains points from neighbor type B
            PTSR=0 #each cell gains points from neighbor type R
            DTRES=0 #change in reporter due to signal S
            SECLPTSR=0 #second signaling loop green points
            SECLPTSG=0 #second signaling loop red points

            for neighbor, commonSurfaceArea in self.get_cell_neighbor_data_list(cell):
                if neighbor is None:
                    continue
                if neighbor.type==1:
                    CSAY+=commonSurfaceArea
                    if mcs<FASTAPPROX:
                        PTSY+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else:
                        PTSY+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                if neighbor.type==2:
                    CSAG+=commonSurfaceArea
                    if mcs<FASTAPPROX:
                        PTSG+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else:
                        PTSG+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                    SECLPTSG+=commonSurfaceArea*neighbor.dict["PTS"][0]/(neighbor.surface)
                if neighbor.type==3:
                    CSAB+=commonSurfaceArea                    
                    if mcs<FASTAPPROX:
                        PTSB+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else:
                        PTSB+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                if neighbor.type==4:
                    CSAR+=commonSurfaceArea
                    if mcs<FASTAPPROX:
                        PTSR+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else:
                        PTSR+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                    SECLPTSR+=commonSurfaceArea*neighbor.dict["PTS"][0]/(neighbor.surface)
            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR) #alternative method to calculate common surface area with medium                 

# VETTING CODE                                        
#             if cell.id==3:
#                 for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell): #iterate for each cell its neighbors
#                     if neighbor:
#                         print "NID", neighbor.id, "TY", neighbor.type, "NCSA", commonSurfaceArea, "DICT", neighbor.dict["PTS"][0], "NS", neighbor.surface
#                 print "CID", cell.id, "CT",cell.type, "CS", cell.surface
#                 print CSAY, CSAG, PTSG, CSAB, PTSB, CSAR, PTSR, CSAM
           
            #if (cell.type==1 or cell.type==2):
                #DTRES=(1/(ALPHAYG+math.exp(-((PTSB+PTSR+SECLPTSG)-BETAYG)/EPSILONYG)))-(1/KAPPAYG)*cell.dict["PTS"][0]
                #cell.dict["PTS"][0]+=DTRES        
                #COMMENTED OUT TO ALLOW FOR MONODIRECTIONAL ELONGATION        
            if (cell.type==3 or cell.type==4):
                DTRES=(1/(ALPHABR+math.exp(-((PTSY+PTSG+SECLPTSR)-BETABR)/EPSILONBR)))-(1/KAPPABR)*cell.dict["PTS"][0]
                cell.dict["PTS"][0]+=DTRES  
                #(PTSB + PTSR + SECLPTSG) input for cell type 1 and 2
                # although negative sign, it is positive signaling 
            if cell.type==1: #change Y cell state
                if cell.dict["PTS"][0]>=THRESHOLDUPYG:
                    cell.type=2
                    cell.dict["P"][0]=1 #iterate one for activating
            if cell.type==2: #change G cell state
                if cell.dict["PTS"][0]<THRESHOLDDOYG:
                    cell.type=1
                    cell.dict["P"][0]=0 #set to 0 for deactivating
            if cell.type==3: #change B cell state
                if cell.dict["PTS"][0]>=THRESHOLDUPBR:
                    cell.type=4
                    cell.dict["P"][1]=1 #iterate one for activating
            if cell.type==4: #change R cell state
                if cell.dict["PTS"][0]<THRESHOLDDOBR:
                    cell.type=3
                    cell.dict["P"][1]=0 #set to 0 for deactivating

            if cell.type==1: #gray cells
                SUMBFSF+=CSAM #grays cell surface area count under bright field
                SUMBFVL+=cell.volume #gray cell volume count under bright field              
                CSAYGBR+=(CSAB+CSAR)/cell.surface #total common surface area of YG cells to BR cells normailzed to YG cell surface
                if (CSAB+CSAR)>0:                 # count the number of YG cells in contact with BR cells
                    YGCBR+=1
                CSAYGYG+=(CSAY+CSAG)/cell.surface # count the number of YG cells in contact with YG cells normailzed to YG cell surface Think of as proportions
                if (CSAY+CSAG)>0:                 # count the number of Yg cells in contact with YG cells
                    YGCYG+=1                    
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity  
                NUMTY+=1                          #count the number if gray cells
                cell.fluctAmpl=BASALY#+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted
                YGPTS+=cell.dict["PTS"][0] #count number of points YG cells have

            if cell.type==2: #green cells
                NAG+=cell.dict["P"][0] #number of activated green cells due to activation at any given time step
                SUMBFSF+=CSAM          #green cells surface area under bright field
                SUMBFVL+=cell.volume   #green cell volume under bright field
                SUMMRSF+=(CSAM+CSAY+CSAB) #green cells count under color field, medium blue and gray are invisible
                SUMMRVL+=cell.volume      #green cells contribute to color field volume         
                CSAYGBR+=(CSAB+CSAR)/cell.surface #total common surface area of YG cells to BR cells normailzed to YG cell surface
                if (CSAB+CSAR)>0:                 # count the number of YG cells in contact with BR cells
                    YGCBR+=1
                CSAYGYG+=(CSAY+CSAG)/cell.surface # count the number of YG cells in contact with YG cells normailzed to YG cell surface Think of as proportions
                if (CSAY+CSAG)>0:                 # count the number of Yg cells in contact with YG cells
                    YGCYG+=1
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity    
                NUMTG+=1                          #count the number of green cells
                cell.fluctAmpl=BASALGR#+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted
                YGPTS+=cell.dict["PTS"][0] #count number of points YG cells have
             
            if cell.type==3: #blue cells            
                SUMBFSF+=CSAM # blue surface area contributes to bright field surface area
                SUMBFVL+=cell.volume #blue volume contributes to bright field volume                
                CSABRYG+=(CSAY+CSAG)/cell.surface #total common surface area of BR cells to YG cells normalized
                if (CSAY+CSAG)>0:                 # count number of BR cells in contact with YG cels
                    BRCYG+=1
                CSABRBR+=(CSAB+CSAR)/cell.surface #total common surface area of BR cells to BR cells normalized
                if (CSAB+CSAR)>0:                 # count number of BR cells in contact with BR cels
                    BRCBR+=1                    
                cell.lambdaSurface=1.0           #change depending on cell adhesitivity
                cell.lambdaVolume=1.0            #change depending on cell adhesitivity      
                NUMTB+=1                         #count number of blue cells
                cell.fluctAmpl=BASALB#+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR)/cell.surface # corrected cell motility, vetted
                BRPTS+=cell.dict["PTS"][0] #count number of points BR cells have
                
            if cell.type==4: #red cells
                NAR+=cell.dict["P"][1] #number of activated red cells due to activation at any given time step
                SUMBFSF+=CSAM          #red cells are visible under bright field and thus conribute surface area
                SUMBFVL+=cell.volume   #red cell volume contributes to bright field                
                SUMMRSF+=(CSAM+CSAY+CSAB) #red cells contribute to color field, medium blue and gray are invisible
                SUMMRVL+=cell.volume      #red cells contribute to color field volume             
                CSABRYG+=(CSAY+CSAG)/cell.surface #common surface area of blue cells to gray cells
                if (CSAY+CSAG)>0:                 #count number of blue red cells that see y g cells
                    BRCYG+=1
                CSABRBR+=(CSAB+CSAR)/cell.surface #common surface area of BR cells that see BR cells
                if (CSAB+CSAR)>0:                 #count the number of BR cells that are in contact with BRBR cells
                    BRCBR+=1                      
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity      
                NUMTR+=1                          #count number of red cells
                cell.fluctAmpl=BASALGR#+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR)/cell.surface #corrected cell motility, vetted
                BRPTS+=cell.dict["PTS"][0] #count number of points of BR cells
                    
            SYPSI+=cell.fluctAmpl #sensitive measure to sorting events
            
#vetting code
            #print cell.id, cell.type, cell.fluctAmpl

            if mcs%RESOL==0: #record points for signle cells traces        
                if cell.id==MARKEDCELLS[0]:
                    self.pW5.add_data_point("1", mcs, cell.dict["PTS"][0])        
                    self.pW5.add_data_point("1T", mcs, cell.type) 
                if cell.id==MARKEDCELLS[1]:
                    self.pW5.add_data_point("2", mcs, cell.dict["PTS"][0])        
                    self.pW5.add_data_point("2T", mcs, cell.type)
                if cell.id==MARKEDCELLS[2]:
                    self.pW5.add_data_point("3", mcs, cell.dict["PTS"][0])        
                    self.pW5.add_data_point("3T", mcs, cell.type)                           
                if cell.id==MARKEDCELLS[3]:
                    self.pW5.add_data_point("4", mcs, cell.dict["PTS"][0])        
                    self.pW5.add_data_point("4T", mcs, cell.type)  
        
        if mcs%RESOL==0:
            
            self.pW2.add_data_point("Hamiltonian", mcs, SYPSI)
            
            self.pW3.add_data_point("Y", mcs, NUMTY)
            self.pW3.add_data_point("G", mcs, NUMTG) 
            self.pW3.add_data_point("B", mcs, NUMTB)   
            self.pW3.add_data_point("R", mcs, NUMTR)
            self.pW3.add_data_point("AR", mcs, NAR)
            self.pW3.add_data_point("AG", mcs, NAG)

            self.pW4.add_data_point("YG", mcs, YGPTS/(NUMTY+NUMTG))
            self.pW4.add_data_point("BR", mcs, BRPTS/(NUMTB+NUMTR))
            
            if YGCYG==0:
                self.pW1.add_data_point("IYG", mcs, 0)
            if YGCYG>0:
                self.pW1.add_data_point("IYG", mcs, CSAYGYG/(YGCYG))
            if BRCBR==0:
                self.pW1.add_data_point("IBR", mcs, 0)
            if BRCBR>0:
                self.pW1.add_data_point("IBR", mcs, CSABRBR/(BRCBR))
                
            if YGCBR==0:
                self.pW1.add_data_point("SYG", mcs, 0)
            if YGCBR>0:
                self.pW1.add_data_point("SYG", mcs, CSAYGBR/(YGCBR))
            if BRCYG==0:
                self.pW1.add_data_point("SBR", mcs, 0)
            if BRCYG>0:
                self.pW1.add_data_point("SBR", mcs, CSABRYG/(BRCYG))
            
            if SUMBFVL==0:
                self.pW6.add_data_point("BF", mcs, 0)                
            if SUMBFVL>0: 
                self.pW6.add_data_point("BF", mcs, (SUMBFSF/SUMBFVL))
            if SUMMRVL==0:
                self.pW6.add_data_point("MR", mcs, 0)
            if SUMMRVL>0:
                self.pW6.add_data_point("MR", mcs, (SUMMRSF/SUMMRVL))
             
            if mcs==ENDMCS:
                fileName = "Calibrate" + str(mcs) + ".txt"
                self.pW1.savePlotAsData(fileName)                
                fileName = "PSI" + str(mcs) + ".txt"
                self.pW2.savePlotAsData(fileName)
                fileName = "FOU" + str(mcs) + ".txt"
                self.pW3.savePlotAsData(fileName)
                fileName = "SIG" + str(mcs) + ".txt"
                self.pW4.savePlotAsData(fileName)
                fileName = "SCSIG" + str(mcs) + ".txt"
                self.pW5.savePlotAsData(fileName)                 
                fileName = "Sphericity" + str(mcs) + ".txt"
                self.pW6.savePlotAsData(fileName)                 
                self.stopSimulation()                 
    def finish(self):
        pass

#from PySteppablesExamples import MitosisSteppableBase

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_frequency=1):
        MitosisSteppableBase.__init__(self, _frequency)
        self.set_parent_child_position_flag(0) #randomize child cell position, see developer manual
    def step(self,mcs):        
        cells_to_divide=[]          #gen cells to divide list
        for cell in self.cellList:
            if cell.type==1: #this line decides which cells are proliferating; see powerpoint for classification           
                cell.dict["RDM"]+=RNG.uniform(MTFORCEMIN1,MTFORCEMAX1) #make cells grow in target radius by this much
                cell.targetSurface=4*math.pi*cell.dict["RDM"]**2 #spherical surface area
                cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
                if cell.volume>2*(4/3)*math.pi*RADAVG**3: #divide at two times the mean radius initialized with               
                    cells_to_divide.append(cell)           #add these cells to divide list
            if cell.type==3: #this line decides which cells are proliferating; see powerpoint for classification           
                cell.dict["RDM"]+=RNG.uniform(MTFORCEMIN3,MTFORCEMAX3) #make cells grow in target radius by this much
                cell.targetSurface=4*math.pi*cell.dict["RDM"]**2 #spherical surface area
                cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
                if cell.volume>2*(4/3)*math.pi*RADAVG**3: #divide at two times the mean radius initialized with               
                    cells_to_divide.append(cell)           #add these cells to divide list
            if cell.type==4: #this line decides which cells are proliferating; see powerpoint for classification           
                cell.dict["RDM"]+=RNG.uniform(MTFORCEMIN4,MTFORCEMAX4) #make cells grow in target radius by this much
                cell.targetSurface=4*math.pi*cell.dict["RDM"]**2 #spherical surface area
                cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
                if cell.volume>2*(4/3)*math.pi*RADAVG**3: #divide at two times the mean radius initialized with               
                    cells_to_divide.append(cell)           #add these cells to divide list
                
        for cell in cells_to_divide:
            self.divideCellRandomOrientation(cell)  #divide the cells

    def updateAttributes(self):
        self.parentCell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #reassign new target radius
        self.parentCell.targetVolume=(4/3)*math.pi*self.parentCell.dict["RDM"]**3 #new target volume
        self.parentCell.targetSurface=4*math.pi*self.parentCell.dict["RDM"]**2 #new target surface area
        self.cloneParent2Child()  #copy characterstics to child cell, indlucig signaling
        self.childCell.dict["P"][0]=0 #reset the activation counter, we dont care about cells from activated parent
        self.childCell.dict["P"][1]=0 #reset the activation counter, we dont care about cells from activated parent
        