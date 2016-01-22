import random
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool


class PSII(object):
    """
    Representation of a LHCII particle
    """   

    def __init__(self,state = "ground", Intensity = 75,timestep=2.5E-7):
        """

        Initialize a LHCII instance, saves all parameters as attributes of the instance.
        """

        self.timestep=timestep #in seconds
        self.Intensity = Intensity #in W/cm^2
        self.state = state
        self.absCrossection =7E-15 #cm^2
        self.absorptionrate = self.Intensity * self.absCrossection/(3.14*10**-19) #per second
        self.absorptionProbability=self.absorptionrate*self.timestep  #absorption probability per timestep
        self.probabilityDecay = 1-np.exp(-timestep/1.5E-9) #Probability for excited state decay during timestep without a car triplet
        self.probabilityDecayTriplet = 1-np.exp(-timestep/150E-12) #Probability for excited state decay during timestep in the presence of a car triplet
        self.ChlTriplet=0 # 0 equals no triplet; 1 equals a triplet
        self.CarTriplet=0 # 0 equals no triplet; 1 equals a triplet
        self.CarTripletDecay=1-np.exp(-timestep/9E-6) #Probability for car triplet state decay during timestep
        self.ChlTripletDecay=1-np.exp(-timestep/2E-3) #Probability for Chl triplet state decay during timestep
        self.ChlTripletYield=0.1
        self.CarTripletYield=0.1
        self.FlYield=0.180
        self.FlYieldTriplet=0.01
        



    def updatePhotonFlux(self, Intensity):
        """
        Updates the excitation intensity and dependent variables

        input: int
        """    
        self.absorptionrate = self.Intensity * self.absCrossection/(3.14*10**-19) #per second
        self.absorptionProbability=self.absorptionrate*self.timestep  

    def doesFluoresce(self):
        """
        Checks for fluorescence from the excited LHCII complex, the probablities are different for the two cases:
            - decay in the presence of a car triplet
            - decay in the absence of a car triplets

        returns boolean: True if photon is fluoresced False otherwise
        """          
        if self.state == "excited" and self.CarTriplet==0 and self.ChlTriplet==0:
            if random.random() <= self.probabilityDecay:
                self.state = "ground"
                if random.random() <= self.FlYield:	
                    return True
                if random.random()<= self.ChlTripletYield:
                    self.ChlTriplet=1
                    return False
                elif random.random()<= self.CarTripletYield:
                    self.CarTriplet+=1
                    return False
                else:
                    return False                                    
            else:
                return False
        if self.state == "excited" and (self.CarTriplet>=1 or self.ChlTriplet>=1):
            if random.random() <= self.probabilityDecayTriplet:
                self.state = "ground"
                if random.random()<= self.CarTripletYield/10.0:
                    self.CarTriplet+=1
                if self.CarTriplet+self.ChlTriplet>=2:
                    if random.random() <= self.FlYieldTriplet:
                        return True
                elif self.CarTriplet+self.ChlTriplet==1:
                    if random.random() <= self.FlYieldTriplet:
                        return True
                else:
                    return False
            else:
                return False
        
        else:
            return False

    def update(self, light):
        """
        Function describing the space of possible transition without or under illumination.
    

        Input:
            light: str "on" or "off" representing if the photon flux will be hitting the RCs during a timestep

        returns a pair of booleans: True/False if a photon is absorbed and True/False if a photon is fluoresced ??
        """            
        if light == "off":
            if self.ChlTriplet>=1:
                if random.random() <= self.ChlTripletDecay:
                    self.ChlTriplet-=1 
            if self.CarTriplet>=1:
                if random.random() <= self.CarTripletDecay:
                    self.CarTriplet-=1 
            if self.state == "excited":
                return False, self.doesFluoresce() 
            else:
                return False,False

        if light == "on":  
            if self.ChlTriplet>=1:
                if random.random() <= self.ChlTripletDecay:
                    self.ChlTriplet-=1
            if self.CarTriplet>=1:
                if random.random() <= self.CarTripletDecay:
                    self.CarTriplet-=1       
            Absorbed = False
            if random.random() <= self.absorptionProbability:
                Absorbed = True
                if self.state == "ground":
                    self.state = "excited"
                    return Absorbed, self.doesFluoresce()
                if self.state == "excited":
                    #print 'Singlet=Singlet annihilation!'
                    return Absorbed, self.doesFluoresce()
                    
            else:
                if self.state == "excited":
                    return False, self.doesFluoresce() 
                else:   
                    return False, False  
                    
def simulation(repetitions=1000000,Intensity=75,light='on'):
    complex=PSII(Intensity=Intensity)
    fluorescence=0
    SumChlTriplets=0
    SumCarTriplets=0
    for num in range(repetitions):
        SumChlTriplets+=complex.ChlTriplet
        SumCarTriplets+=complex.CarTriplet
        Abs,Fl= complex.update(light)
        if Fl==True:
            fluorescence+=1
    DetectionEfficiency=1
    fluorescence=fluorescence/float(repetitions*13.14E-9)*DetectionEfficiency #converted to counts per second and adjusted for the detection efficiency of our setup
    ChlTripletPro=SumChlTriplets/float(repetitions)
    CarTripletPro=SumCarTriplets/float(repetitions)
    return fluorescence,ChlTripletPro
            
def saturation(intensities): #Models the saturation curve for different excitation power
    Fl=[]
    Tr=[]
    for e in intensities:
        print e
        Fluo,Trip=simulation(Intensity=e)
        Fl.append(Fluo)
        Tr.append(Trip)
    plt.plot(intensities,Fl)
    plt.xlabel('Excitation intensity [W/cm^2]')
    plt.ylabel('Fluorescence Intensity [cps]')
    plt.title('Fluorescence saturation curve due to triplets')
    plt.figure(2)
    plt.plot(intensities,Tr)
    plt.xlabel('Excitation intensity [W/cm^2]')
    plt.ylabel('Average population of Car triplet states')
    plt.title('Amplitude ratio of S--T annihilation')
    plt.show()
    return Fl
    
#saturation([10,50,150,300,1000,2000])

def simulationAOM(numtrials=1,AOMtimes=[2.5E-3,10E-3],Intensity=75,ChlTripletYield=0.1,CarTripletYield=0.001,binning=2E-5):
    complex2=PSII(Intensity=Intensity)
    complex2.ChlTripletYield=ChlTripletYield
    complex2.CarTripletYield=CarTripletYield
    complex2.FlYield=0.15
    complex2.FlYieldTriplet=0.015
    timestep=float(complex2.timestep)
    SumChlTriplets=[]
    SumCarTriplets=[]
    fluorescence=[]
    Absorbed=[]
    Annihilation=[]
    binning=binning
    num_bins=int(AOMtimes[0]/binning)
    for i in range(num_bins):
        fluorescence.append(0)
        SumChlTriplets.append(0)
        SumCarTriplets.append(0)
        Absorbed.append(0)
        Annihilation.append(0)
    for e in range(numtrials):
        for num in range(int(AOMtimes[0]/timestep)):            
            Abs,Fl= complex2.update('on')
            SumChlTriplets[int(num*timestep/AOMtimes[0]*num_bins)]+=complex2.ChlTriplet
            SumCarTriplets[int(num*timestep/AOMtimes[0]*num_bins)]+=complex2.CarTriplet
            if Fl==True:
                fluorescence[int(num*timestep/AOMtimes[0]*num_bins)]+=1
                if (complex2.ChlTriplet+complex2.CarTriplet)>=1:
                    Annihilation[int(num*timestep/AOMtimes[0]*num_bins)]+=1
            if Abs==True:
                Absorbed[int(num*timestep/AOMtimes[0]*num_bins)]+=1
                
        complex2.CarTripletDecay=1-np.exp(-(float(AOMtimes[1])/3.0)/9.0E-6)
        complex2.ChlTripletDecay=1-np.exp(-(float(AOMtimes[1])/3.0)/2.0E-3)
        for num in range(3):
            Abs,Fl= complex2.update('off')
            #if Fl==True:
            #    fluorescence[int(num/binning)]+=1
        complex2.CarTripletDecay=1-np.exp(-timestep/9.0E-6)
        complex2.ChlTripletDecay=1-np.exp(-timestep/2.0E-3)
    return fluorescence, SumChlTriplets,SumCarTriplets,Absorbed,Annihilation

def AOM(numtrials=500,AOMtimes=[0.8E-3,0.1E-3],Intensities=[500]):
    #random.seed(1)
    ChlTripletYield=0.02
    CarTripletYield=0.15
    plt.figure(2)
    plt.figure(1)
    plt.figure(3)
    Offtimes=[10E-3,1.5E-3,0.5E-3,0.1E-3]
    plt.clf()
    colors=['k','g','darkkhaki','r']
    binning=2E-5
    for j in range(4):
        print j
        fluorescence,SumChlTriplets,SumCarTriplets,Absorbed,Annihilation=simulationAOM(numtrials,AOMtimes=[AOMtimes[0],Offtimes[j]],Intensity=Intensities[0],ChlTripletYield=ChlTripletYield,CarTripletYield=CarTripletYield,binning=binning)
        fluorescence=np.asarray(fluorescence)
        SumChlTriplets=np.asarray(SumChlTriplets)
        SumCarTriplets=np.asarray(SumCarTriplets)
        Absorbed=np.asarray(Absorbed)
        Annihilation=np.asarray(Annihilation)
        AvgTripletPop=Annihilation/0.15/(Annihilation/0.15+(fluorescence-Annihilation)/1.5)
        xaxis=[]
        for i in range(len(fluorescence)):
            xaxis.append(AOMtimes[0]/float(len(fluorescence))*i*1E3)  
        plt.figure(1)      
        plt.bar(xaxis,fluorescence,width=binning*1E3,color=colors[j])
        plt.xlabel('AOM time [ms]',size=15)
        plt.figure(2)
        plt.plot(xaxis,SumChlTriplets,color=colors[j])
        plt.xlabel('AOM time [ms]',size=15)
        plt.figure(3)
        #plt.plot(xaxis,SumCarTriplets,color=colors[j])
        #plt.xlabel('AOM time [ms]',size=15)
        plt.plot(xaxis,SumCarTriplets,color=colors[j])
        plt.xlabel('AOM time [ms]',size=15)
    plt.figure(1)    
    plt.ylabel('Fluorescence Intensity [a.u.]',size=15)
    #plt.ylim([0,700])
    plt.xlim([-0.1,0.9])
    plt.legend(('Offtime [ms]: ' + str(Offtimes[0]*1E3),'Offtime [ms]: ' + str(Offtimes[1]*1E3),'Offtime [ms]: ' + str(Offtimes[2]*1E3),'Offtime [ms]: ' + str(Offtimes[3]*1E3)))
    plt.title('Simulation of AOM kinetics on a C2S2 supercomplex',size=15)
    #plt.title(('ChlTripletYield: ' + str(ChlTripletYield) + ' \nCarTripletYield: ' + str(CarTripletYield)),size=15)
    #text1=plt.figtext(0.5, 0.6,r'$\tau_{Car}=9\mu s$' +'\n' + r'$\tau_{Chl}=2 ms$',size=15)

    plt.show()
    
AOM(Intensities=[75])



            