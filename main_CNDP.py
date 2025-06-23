#---modules
from src import Network
from src import OA_elastic_CG
from decimal import Decimal
#import polytope as pc

net = 'Braess'
ins = 'Braess_CNDP_1'

net = 'SiouxFalls'
ins = 'SF_CNDP_10_1'

#net = 'EasternMassachusetts'
#ins = 'EM_CNDP_30_1'

#net = 'HarkerFriesz'
#ins = 'HF_CNDP_1'

#net = 'NguyenDupuis'
#ins = 'ND_CNDP_1'

#net = 'Anaheim'
#ins = 'A_CNDP_30_1'

b_prop = 0.5
scal_flow = {'SiouxFalls':1e-2,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1, 'NguyenDupuis':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':0.25, 'NguyenDupuis':1}
print(net,ins)

inflate_cost = 5
scale_dem = 1
print(scale_dem * inflate_trips[net], inflate_cost, scal_flow[net])

network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])

#network.tapas("UE", None)
#network.checkDualBeckmann()
#network.generateScenarios(1, 5, 1.0/3, 0.1)

test = OA_elastic_CG.OA_elastic_CG(network, False, 0.5)

test.solve()
