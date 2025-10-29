#---modules
from src import Network
from src import OA_dual_CG
from decimal import Decimal
#import polytope as pc
import sys

def main():

	args = sys.argv[1:]
	
	print(args)
	
	net = args[0]
	scenario = args[1]
	
	#net = 'Braess'
	if net == "Braess":
		ins = 'Braess_CNDP_1'
	#ins = 'Braess_CNDP_1'

	#net = 'SiouxFalls'
	if net == "SiouxFalls":
		ins = 'SF_CNDP_10_1'
	#ins = 'SF_CNDP_10_1'

	#net = 'EasternMassachusetts'
	if net == "EasternMassachusetts":
		ins = 'EM_CNDP_30_1'
	#ins = 'EM_CNDP_30_1'

	#net = 'HarkerFriesz'
	if net == "HarkerFriesz":
		ins = 'HF_CNDP_1'
	#ins = 'HF_CNDP_1'

	#net = 'NguyenDupuis'
	if net == "NguyenDupuis":
		ins = 'ND_CNDP_1'
	#ins = 'ND_CNDP_1'

	#net = 'Anaheim'
	if net == "Anaheim":
		ins = 'A_CNDP_30_1'
	
	if net == "BerlinMitteCenter":
		ins = 'BMC_CNDP_30_1'
		

	b_prop = 0.5
	scal_flow = {'SiouxFalls':1e-2,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1, 'NguyenDupuis':1}
	inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':0.25, 'NguyenDupuis':1}
	print(net,ins)

	inflate_cost = 5
	scale_dem = 1
	print(scale_dem * inflate_trips[net], inflate_cost, scal_flow[net])

	network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
	y = network.initCalcY()
	
	
	network.tapas("UE", y)
	#print("beckmann", network.getBeckmannOFV())
	#print("dual", network.getDualBeckmannOFV())
	#network.generateScenario("40_10_10_2", 0.4, 0.1, 0.1)
	network.generateScenario("20_20_20_2", 0.2, 0.2, 0.2)
	#network.generateScenario("20_10_10_2", 0.2, 0.1, 0.1)
	#network.generateScenario("20_10_10_1", 0.2, 0.1, 0.1)
	#network.generateScenario("20_10_10_2", 0.2, 0.1, 0.1)
	#network.generateScenario("20_10_10_3", 0.2, 0.1, 0.1)
	
	
	
	#print("starting tapas")

	#network.tapas('UE', y)

	#network.generateTarget(y)


	#for a in network.links:
	#    print(a, a.x, a.getTravelTime(a.x, "UE"))

	#print("TD", network.TD)

	#test = OA_CNDP_elastic.OA_CNDP_elastic(network)
	test = OA_dual_CG.OA_dual_CG(network, False, 0.5, scenario)

	#obj, tot_time, tap_time, iter, = test.solve()
	#test.solve()
	
if __name__=="__main__":
    main()