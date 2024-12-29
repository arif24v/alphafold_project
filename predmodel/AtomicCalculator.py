import numpy as np

class AtomicCalculator:
    #* Sigma (in Angstrom), Epsilon (in Kcal/mol)
    #* https://github.com/MolecularFlipbook/FlipbookApp/blob/master/mfb/MGLToolsPckgs/AutoDockTools/AD4.1_bound.dat
    LJ_PARAMETERS = {
    'C':(4.00,0.150),
    'N':(3.50, 0.160),
    'H':(2.00, 0.020),
    'CL':(4.09, 0.276),
    'OXT':(3.20,0.2),
    'F':(3.09,0.08),
    'S':(4.0,.2),
    'O':(3.20,0.200),
    'CA':(1.98,0.55),
    'BR':(4.33,0.389),
    'O1-':(3.20,0.2),
    'N1+':(3.5,0.160)
    }
    
    def __init__(self):
        pass  

    @staticmethod
    def lennard_jones_potential(self, particle1, particle2):
        sigma_1, epsilon_1 = self.LJ_PARAMETERS[particle1[0]]
        sigma_2, epsilon_2 = self.LJ_PARAMETERS[particle2[0]]
        r=np.sqrt((particle1[1]-particle2[1])**2+(particle1[2]-particle2[2])**2+(particle1[3]-particle2[3])**2)

        sigma = (sigma_1 + sigma_2) / 2
        epsilon = np.sqrt(epsilon_1 * epsilon_2)

        term1 = (sigma / r) ** 12
        term2 = (sigma / r) ** 6

        return 4 * epsilon * (term1 - term2) * 1.39 * 10**-20
    
    def force(self,particle1,particle2):
        sigma_1, epsilon_1 = self.LJ_PARAMETERS[particle1[0]]
        sigma_2, epsilon_2 = self.LJ_PARAMETERS[particle2[0]]
        r = np.array([particle1[1]-particle2[1],particle1[2]-particle2[2],particle1[3]-particle2[3]])
        distance = np.linalg.norm(r)
        sigma = (sigma_1 + sigma_2) / 2
        epsilon = np.sqrt(epsilon_1 * epsilon_2)
        
        #* kg/(m*mol)*6.9477*10^-21=N
        x = (sigma/distance)**6
        return r,distance,24*epsilon*(2 * (x**2)/distance-x/distance) * 6/9477*10**-21