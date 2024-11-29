import numpy as np
import os

def read_coordxvg(coordid, ANALYSIS_DIR):
    """Reads data from the XVG relative to the coordid given in input. It returns a numpy array N,2"""
    coordxvgname = f"cphmd-coord-{coordid}.xvg"
    coordxvg = os.path.join(ANALYSIS_DIR, coordxvgname)
    return np.loadtxt(coordxvg, comments=["@","#"])

def get_frac(array_xvg):
    """Count protonation and deprotonation states. It returns the total protonation or deprotonation fraction"""
    prot = array_xvg[np.where(array_xvg[:,1] < 0.2) ] #When lambda lower than 0.2
    deprot = array_xvg[np.where(array_xvg[:,1] > 0.8) ] #when lambda higher than 0.8

    prot_frac = len(prot) / (len(prot) + len(deprot))
    deprot_frac = len(deprot) / (len(prot) + len(deprot))

    return prot_frac, deprot_frac

def get_stats(coordid, dirs):
    """Read XVG and return statistics on total protonation fraction"""
    xvg_dict = {}
    prot_dict = {}
    deprot_dict = {}
    coordxvgname = f"cphmd-coord-{coordid}.xvg"

    for folder in dirs:
        xvg_dict[folder] = []

        coordxvg = os.path.join(folder, coordxvgname)
        
        xvg_dict[folder] = np.loadtxt(coordxvg, comments=["@","#"])

    for i, run in enumerate(xvg_dict.keys()):

        prot_dict[i], deprot_dict[i] = get_frac(xvg_dict[run])

    prot_avg = (sum(prot_dict.values())/len(prot_dict.values()))
    deprot_avg = (sum(deprot_dict.values())/len(deprot_dict.values()))

    prot_se = np.std(list(prot_dict.values()))/np.sqrt(len(prot_dict.values()))
    deprot_se = np.std(list(deprot_dict.values()))/np.sqrt(len(deprot_dict.values()))


    return prot_avg, deprot_avg, prot_se, deprot_se

def get_protfrac_ts(coordid, ANALYSIS_DIR):
    """Reads data from the XVG relative to the coordid given in input. It returns the time-series of the protonation fraction."""
    array_xvg = read_coordxvg(coordid, ANALYSIS_DIR)

    lista_prot = [] #Keep track of protonation
    lista_deprot = [] #Keep track of deprotonation

    for lambda_value in array_xvg[:,1]:
        #If lambda < 0.2 count as protonated
        if lambda_value < 0.2:
            lista_prot.append(1)
            lista_deprot.append(0)
        #If lambda > 0.8 count as deprotonated
        elif lambda_value > 0.8:
            lista_prot.append(0)
            lista_deprot.append(1)
        #If lambda > 0.8 count as deprotonated
        else:
            lista_prot.append(0)
            lista_deprot.append(0)

    #Cumlative sum along time
    protcumsum_array = np.cumsum(lista_prot)
    deprotcumsum_array = np.cumsum(lista_deprot)

    return protcumsum_array/(protcumsum_array + deprotcumsum_array) #Protonation fraction

def get_HISstats(coordids, dirs):
    """Read Histidines XVG and return statistics on total protonation fraction"""
    xvg_dict = {}
    prot_dict = {}


    for folder in dirs:
        xvg_dict[folder] = [] #Dict with a list for each directory

        for coordid in coordids:

            coordxvg = os.path.join(folder, f"cphmd-coord-{coordid}.xvg")
        
            xvg_dict[folder].append(np.loadtxt(coordxvg, comments=["@","#"])) #Each directory has a list of the three lambda states

    for i, run in enumerate(xvg_dict.keys()):

        lambda_list = xvg_dict[run]
        prot_dict[i] = get_HISfrac(lambda_list)

    prot_avg = (sum(prot_dict.values())/len(prot_dict.values()))

    prot_se = np.std(list(prot_dict.values()))/np.sqrt(len(prot_dict.values()))

    return prot_avg, prot_se

def get_HISfrac(array_xvgs):
    """Count protonation and deprotonation states of HISTIDINES. array_xvgs is the list containing the three HIS lambda states arrays. It returns the total protonation fraction"""
    nde = np.sum(array_xvgs[1][:,1]>0.8) + np.sum(array_xvgs[2][:,1]>0.8) #Sum neutral states HSD/HSE
    npr = np.sum(array_xvgs[0][:,1]>0.8) #Histidine protonated form HSP

    prot_frac = (0.+npr)/(nde+npr)

    return prot_frac

#Corrd id to lambda id
def coord2lambdaid(lambdareference):
    """Read lambda reference file and return a dictionary with coordid as keys and lambdaid as values"""
    
    coord2lambda = {}

    with open(lambdareference, "r") as fi:
        
        lambdaid = 0
        coordid = 0
        
        for line in fi.readlines():
            
            if line.startswith("resname"):
                continue
            
            elif line[5] != "1":
                lambdaid = lambdaid
                coordid += 1
            
            else:
                lambdaid +=1
                coordid +=1

        coord2lambda[coordid] = lambdaid
    return coord2lambda