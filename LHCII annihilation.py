import random
import matplotlib.pyplot as plt
import numpy as np


class LHCII(object):
    """
    Representation of a LHCII particle
    """   

    def __init__(self,state = "ground", Intensity = 75,timestep=13.14E-9):
        """

        Initialize a LHCII instance, saves all parameters as attributes of the instance.
        """

        self.timestep=timestep #in seconds
        self.Intensity = Intensity #in W/cm^2
        self.state = state
        self.absCrossection = 1.4E-15 #cm^2
        self.absorptionrate = self.Intensity * self.absCrossection/(3.14*10**-19) #per second
        self.absorptionProbability=self.absorptionrate*self.timestep  #absorption probability per timestep
        self.probabilityDecay = 1-np.exp(-timestep/3.5E-9) #Probability for excited state decay during timestep without a car triplet
        self.probabilityDecayTriplet = 1-np.exp(-timestep/35E-12) #Probability for excited state decay during timestep in the presence of a car triplet
        self.triplet=0 # 0 equals no triplet; 1 equals a triplet
        self.TripletDecay=1-np.exp(-timestep/9E-6) #Probability for car triplet state decay during timestep



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
        if self.state == "excited" and self.triplet==0:
            if random.random() <= self.probabilityDecay:
                self.state = "ground"
                if random.random() <= 0.33:	
                    return True
                if random.random()<= 0.33:
                    self.triplet=1
                    return False
                else:
                    return False                                    
            else:
                return False
        if self.state == "excited" and self.triplet==1:
            if random.random() <= self.probabilityDecayTriplet:
                self.state = "ground"
                if random.random() <= 0.0033:
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
            if self.triplet==1:
                if random.random() <= self.TripletDecay:
                    self.triplet=0 
            if self.state == "excited":
                return False, self.doesFluoresce() 
            else:
                return False,False

        if light == "on":  
            if self.triplet==1:
                if random.random() <= self.TripletDecay:
                    self.triplet=0      
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
    complex=LHCII(Intensity=Intensity)
    fluorescence=0
    SumTriplets=0
    for num in range(repetitions):
        SumTriplets+=complex.triplet
        Abs,Fl= complex.update(light)
        if Fl==True:
            fluorescence+=1
    DetectionEfficiency=0.035
    fluorescence=fluorescence/float(repetitions*13.14E-9)*DetectionEfficiency #converted to counts per second and adjusted for the detection efficiency of our setup
    TripletPro=SumTriplets/float(repetitions)
    return fluorescence,TripletPro
            
def saturation(intensities):
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
    
saturation([10,50,150,300,1000,2000])

def simulationAOM(numtrials=100,AOMtimes=[3*760,3*760],Intensity=75):
    complex2=LHCII(Intensity=Intensity)
    fluorescence=[]
    binning=float(38)
    for i in range(int((AOMtimes[0]+AOMtimes[1])/binning)):
        fluorescence.append(0)
    for e in range(numtrials):
        for num in range(AOMtimes[0]):            
            Abs,Fl= complex2.update('on')
            if Fl==True:
                fluorescence[int(num/binning)]+=1
        for num in range(AOMtimes[0],AOMtimes[1]):
            Abs,Fl= complex2.update('off')
            if Fl==True:
                fluorescence[int(num/binning)]+=1

    return fluorescence

def AOM(numtrials=2000,AOMtimes=[2*760,3*760],Intensities=[500]):
    plt.figure(3)
    colors=['b','r','g']
    for j in range(len(Intensities)):
        print j
        fluorescence=simulationAOM(numtrials,AOMtimes,Intensity=Intensities[j])
        xaxis=[]
        for i in range(len(fluorescence)):
            xaxis.append(i*0.5)        
        plt.bar(xaxis,fluorescence,color=colors[j])
    plt.xlabel('AOM time [us]')
    plt.ylabel('Fluorescence Intensity [a.u.]')
    plt.ylim([0,700])
    plt.legend(('700 W/cm^2', '300 W/cm^2','100 W/cm^2'))
    plt.title('Pulse wave excitation: Triplet accumulation')
    plt.show()
    
AOM(Intensities=[700,300,100])



            