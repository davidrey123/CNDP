#---modules
from src import Network
from src import Benders
from decimal import Decimal
#import polytope as pc
import sys

def main():



	
	net = 'Anaheim'
	
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
		

	num_candidate = 40
	B = num_candidate/2
	max_cost = 10

	network = Network.Network(net, num_candidate, B)
	
	test = Benders.Benders(network, max_cost)


	#test.milp()
	test.compare()
	
if __name__=="__main__":
    main()