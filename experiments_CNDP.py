#---modules
from src import Network
from src import OA_CNDP
from src import OA_CNDP_CG
from src import HY_CNDP
from src import DuGP_CNDP
from src import CNDP_MILP
import math
#import polytope as pc

#import numpy as np
from src import OA_CNDP_CS


runMILP = True
runOA = True



net = 'HarkerFriesz'
ins = 'HF_CNDP_1'

b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':1}
print(net,ins)

network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])

inflate_cost = 1

#test = HY_CNDP.HY_CNDP(network)
#test.solve()


if runMILP:
    filename_milp = 'experiments_milp_'+net+'.txt'
    f_milp = open(filename_milp, "w")
    filename_milp_latex = 'experiments_milp_latex_'+net+'.txt'
    f_milp_latex = open(filename_milp_latex, "w")
    
if runOA:
    filename_oa = 'experiments_oa_'+net+'.txt'
    f_oa = open(filename_oa, "w")

    filename_oa_latex = 'experiments_oa_latex_'+net+'.txt'
    f_oa_latex = open(filename_oa_latex, "w")


header1_milp = "dem_scale\tcost_scale\t"
header1_milp_latex = "dem_scale & cost_scale & "
header2_milp = "\t\t"
header2_milp_latex = '&'


header0_milp_latex = "\begin{tabular}{cc"
header0_oa_latex = "\begin{tabular}{cc"
header1_oa = "dem_scale\tcost_scale\t"
header1_oa_latex = "dem_scale & cost_scale & "
header2_oa = "\t\t"
header2_oa_latex = '&'



for j in range(1, 4):
    pieces = j*5
    header1_milp += "\t"+str(pieces)+"\t"
    header1_milp_latex += "& \multicolumn{3}{|c}{"+str(pieces)+"}"
    header2_milp += "\tobj\trt (s)\tgap" # for 3 variations on num pieces for MILP
    header2_milp_latex += "obj & rt (s) & gap"
    header0_milp_latex += "|ccc"
    
if runOA:
    header1_oa += "\tB cut\t\t\t\tB cut, CG\t\t\t\tlink VF cut\t\t\t\tlink VF cut, CG"
    header1_oa_latex += "& \multicolumn{4}{|c|}{B cut} & \multicolumn{4}{|c|}{B cut, CG} & \multicolumn{4}{|c|}{link VF cut) & \multicolumn{4}{|c}{link VF cut, CG}"

    for j in range(0, 4):
        header2_oa += "\t obj \t rt (s) \t TAP (s) \t iter." # for 4 variations: with and without CG; with and without link-based VF cuts 
        header2_oa_latex += "& obj & rt (s) & TAP (s) & iter." 
        header0_oa_latex += "|cccc"

    header0_oa_latex += "}"
    header1_oa_latex += "\\\\"
    header2_oa_latex += "\\\\"

    

if runMILP:
    header1_milp_latex += "\\\\\\toprule"
    header2_milp_latex += "\\\\\\toprule"
    header0_milp_latex += "}"

if runMILP:
    f_milp_latex.write(header0_milp_latex+"\n"+header1_milp_latex+"\n"+header2_milp_latex+"\n")

    f_milp.write(header1_milp+"\n"+header2_milp+"\n")
    
    f_milp.flush()
    f_milp_latex.flush()

if runOA:
    f_oa_latex.write(header0_oa_latex+"\n"+header1_oa_latex+"\n"+header2_oa_latex+"\n")
    f_oa.write(header1_oa+"\n"+header2_oa+"\n")

    f_oa.flush()
    f_oa_latex.flush()
    

for i in range (1, 5):
    scale = i/4
    
    print("\n\n\ndemand scale", scale, "\n\n\n")
    
    if runMILP:
        f_milp_latex.write(str(scale * inflate_trips[net])+" & "+str(inflate_cost)+" ")
        f_milp.write(str(scale * inflate_trips[net])+" \t"+str(inflate_cost))
        
    if runOA:
        f_oa_latex.write(str(scale * inflate_trips[net])+" & "+str(inflate_cost)+" ")


        f_oa.write(str(scale * inflate_trips[net])+" \t"+str(inflate_cost))
    
    if runMILP:
        for j in range(1, 4):
            pieces = 3+j*2

            print("solving MILP with", pieces, "pieces")
            print("demand scale is ", scale * inflate_trips[net], "cost scale is ", inflate_cost, "scale flow is ", scal_flow[net])

            network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],scale * inflate_trips[net])
            test = CNDP_MILP.CNDP_MILP(network, pieces, pieces, 20, inflate_cost)
            obj, tot_time, gap = test.solve()

            time_id = ""

            if tot_time >= 4*3600:
                time_id = "\\tl"

            if obj == 1e100:
                f_milp.write("\tinfeas\t"+str(tot_time)+"\t")
                f_milp_latex.write(" & inf & "+str(round(tot_time, 2))+time_id+" & ")
            else:
                f_milp.write("\t"+str(obj)+"\t"+str(tot_time)+"\t"+str(gap))
                f_milp_latex.write(" & "+str(round(obj, 1))+" & "+str(round(tot_time, 2))+time_id+" & "+str(round(gap, 3)))
                
            if runMILP:
                f_milp.flush()
                f_milp_latex.flush()

            #print(obj, tot_time)
        
    if runOA:
        for j in range (2, 4):
            network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],scale * inflate_trips[net])
            test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useCG= (j%2 == 1), useLinkVF=(j >= 2))
            obj, tot_time, tap_time, iterations = test.solve()

            #print(obj, tot_time, tap_time, iterations)
            f_oa.write("\t"+str(obj)+"\t"+str(tot_time)+"\t"+str(tap_time)+"\t"+str(iterations))
            f_oa_latex.write(" & "+str(round(obj, 1))+" & "+str(round(tot_time, 2))+" & "+str(round(tap_time, 2))+" & "+str(iterations))

        f_oa.write("\n")
        f_oa_latex.write("\\\\ \n")
        
        if runOA:
            f_oa.flush()
            f_oa_latex.flush()

    
    if runMILP:
        f_milp.write("\n")
        f_milp_latex.write("\\\\ \n")
    
    if runOA:
        f_oa.flush()
        f_oa_latex.flush()
        
    if runMILP:
        f_milp.flush()
        f_milp_latex.flush()

if runOA:   
    f_oa_latex.write("\\hline") 
       
if runMILP:
    f_milp_latex.write("\\hline")

   
    
for i in range (0, 5):
    scale = math.pow(2, i-1)
    
    print("\n\n\nsolving MILP with", pieces, "pieces")
    print("cost scale is ", scale, "\n\n\n")
    
    if runMILP:
        f_milp_latex.write(str(inflate_trips[net])+" & "+str(scale*inflate_cost))
        f_milp.write(str(inflate_trips[net])+" \t"+str(scale*inflate_cost))
        
    if runOA:
        f_oa_latex.write(str(inflate_trips[net])+" & "+str(scale*inflate_cost))
        f_oa.write(str(inflate_trips[net])+" \t"+str(scale*inflate_cost))
    
    if runOA:
        for j in range (2, 4):
            network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
            test = OA_CNDP_CG.OA_CNDP_CG(network, scale*inflate_cost, useCG= (j%2 == 1), useLinkVF=(j >= 2))
            obj, tot_time, tap_time, iterations = test.solve()

            #print(obj, tot_time, tap_time, iterations)
            f_oa.write("\t"+str(obj)+"\t"+str(tot_time)+"\t"+str(tap_time)+"\t"+str(iterations))
            f_oa_latex.write(" & "+str(round(obj, 1))+" & "+str(round(tot_time, 2))+" & "+str(round(tap_time, 2))+" & "+str(iterations))

        f_oa.write("\n")
        f_oa_latex.write("\\\\ \n")
        
        if runOA:
            f_oa.flush()
            f_oa_latex.flush()

    
    
    if runMILP:
        for j in range(1, 4):
            pieces = 3+j*2
            network = Network.Network(net,ins,b_prop,1e-0,scale * scal_flow[net],inflate_trips[net])
            test = CNDP_MILP.CNDP_MILP(network, pieces, pieces, 20, inflate_cost)
            obj, tot_time, gap = test.solve()

            time_id = ""

            if tot_time >= 4*3600:
                time_id = "\\tl"

            if obj == 1e100:
                f_milp.write("\tinfeas\t"+str(tot_time)+"\t")
                f_milp_latex.write(" & inf & "+str(round(tot_time, 2))+time_id+" & ")
            else:
                f_milp.write("\t"+str(obj)+"\t"+str(tot_time)+"\t"+str(gap))
                f_milp_latex.write(" & "+str(round(obj, 1))+" & "+str(round(tot_time, 2))+time_id+" & "+str(round(gap, 3)))
            #print(obj, tot_time)
            
            if runMILP:
                f_milp.flush()
                f_milp_latex.flush()

        f_milp.write("\n")
        f_milp_latex.write("\\\\ \n")

    if runOA:
        f_oa.flush()
        f_oa_latex.flush()
        
    if runMILP:
        f_milp.flush()
        f_milp_latex.flush()


if runOA:
    f_oa_latex.write("\\bottomrule")  
    f_oa.close()
    f_oa_latex.close()  
    
if runMILP:
    f_milp_latex.write("\\bottomrule")
    f_milp.close()
    f_milp_latex.close()



