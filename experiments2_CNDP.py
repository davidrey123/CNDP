#---modules
from src import Network
from src import OA_CNDP
from src import OA_CNDP_CG
#from src import HY_CNDP
from src import DuGP_CNDP
from src import CNDP_MILP
import math
#import polytope as pc

#import numpy as np
from src import OA_CNDP_CS
from src import OA_CNDP_CG_SB


CG = True

nets = ['SiouxFalls', 'EasternMassachusetts', 'Anaheim', 'BerlinMitteCenter']
inss = ['SF_CNDP_','EM_CNDP_','A_CNDP_','BMC_CNDP_']

'''
scale_dems = [
[1, 0.5],
[1, 2],
[1, 2],
[1, 2]
]
'''
scale_dems = [
[1],
[1],
[1],
[2]
]


inflate_costs = [
[1,100],
[1,100],
[1,100],
[0.1,1],
]


yvarss = [
[10,20],
[10,30],
[20,30, 60],
[20,30, 60]
]

#net = 'EasternMassachusetts'
#ins = 'EM_CNDP_'


b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':1}

#network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
#test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useLinkVF=True)


inflate_cost = 1
scale_dem = 1




for n in range (1, 4):
    net = nets[n]
    f = open("experiments_"+net+".txt", "w")
    
    for j in range (0, len(yvarss[n])):
        yvars = yvarss[n][j]


        for k in range (1, 3):
            if k == 1:
                inflate_cost = inflate_costs[n][0]
                scale_dem = scale_dems[n][0]
            elif k == 2:
                inflate_cost = inflate_costs[n][1]
                scale_dem = scale_dems[n][0]
            #elif k == 3:
            #    inflate_cost = inflate_costs[n][1]
            #    scale_dem = scale_dems[n][0]

            for i in range (1, 4):
                
                
                actual_ins = inss[n]+str(yvars)+"_"+str(i)
                print("\n\n\n", actual_ins, yvars, scale_dem, inflate_cost, CG, "\n\n\n")

                network = Network.Network(net,actual_ins,b_prop,1e-0,scal_flow[net],scale_dem * inflate_trips[net])
                test = OA_CNDP_CG_SB.OA_CNDP_CG_SB(network, inflate_cost)
                
                if net == "SiouxFalls":
                    test.params.min_gap = 1E-3
                elif net == "EasternMassachusetts":
                    test.params.min_gap = 4E-3
                    
                obj, avgy, tot_time, tap_time, iter, = test.solve()

                scientific_format = "{:.2e}".format(test.getAvgLinkCost())
                
                varlinks_str = str(len(test.varlinks))
                TD_str = str(round(network.TD,1))
                
                if not(k == 1 and i == 1):
                    varlinks_str = ""
                    
                if not(i == 1 and k==1):
                    TD_str = ""
                
                f.write("& "+varlinks_str +  " & " +str(i) + " & " + str(scientific_format) +  " & " + str("{:.1f}".format(obj))+ " & " + str("{:.3f}".format(avgy))+ " & "+ str("{:.2f}".format(100*test.gap))+ "\% & " + str("{:.1f}".format(test.tstt))  + " & " + str("{:.2f}".format(tot_time)) + "s & "+ str("{:.2f}".format(tap_time))+ "s & " + str(iter)+" \\\\")
                f.write("\n")
                
                if i == 3:
                    if k==2:
                        f.write("\cline{2-10}\n")
                    else:
                        f.write("\cline{4-10}\n")

                f.flush()
            
    f.close()