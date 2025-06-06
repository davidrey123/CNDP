#---modules
from src import Network
from src import OA_CNDP
from src import OA_CNDP_CG
from src import OA_CNDP_CG_SB
#from src import HY_CNDP
#from src import DuGP_CNDP
#from src import CNDP_MILP
from decimal import Decimal
#import polytope as pc

#net = 'Braess'
#ins = 'Braess_CNDP_1'

net = 'SiouxFalls'
ins = 'SF_CNDP_20_2'

#net = 'EasternMassachusetts'
#ins = 'EM_CNDP_10_1'

#net = 'HarkerFriesz'
#ins = 'HF_CNDP_1'

#net = 'BerlinMitteCenter'
#ins = 'BMC_CNDP_60_3'

#net = 'Anaheim'
#ins = 'A_CNDP_30_1'

scale_dem = 1
b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':0.25}
print(net,ins)

inflate_cost = 1

print(inflate_trips[net]*scale_dem, inflate_cost, scal_flow[net])

network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net]*scale_dem)
print("TD", network.TD)

#test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost)
test = OA_CNDP_CG_SB.OA_CNDP_CG_SB(network, inflate_cost)
#test = HY_CNDP.HY_CNDP(network)
#test = CNDP_MILP.CNDP_MILP(network, 5, 5, 20, inflate_cost)

#obj, tot_time, tap_time, iter, = test.solve()
test.solve()