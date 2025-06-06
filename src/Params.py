INFTY = 1.0e9

class Params:   
    
    def resetPAS(self):
        # if these are too large, then we avoid flow shifting on PAS
        # if these are too small, then we use PAS that are ineffective
        self.pas_cost_mu = 0.01
        self.pas_flow_mu = 0.01
        
        # if this is too large, then we avoid flow shifting in PAS
        # if this is too small, we waste time shifting in PAS that is useless
        self.pas_cost_epsilon = 0.01
        self.line_search_gap = 1E-8
        
        
    def __init__(self):
        
        #---Numerical rounding precision (in decimals)
        self.rd = 6
        
        #---SP params
        self.SP_tol = 1e-6
        
        #---TAP params
        self.msa_max_iter = 500        
        self.tapas_max_iter = 100 
        self.min_gap = 1E-3
        self.warmstart = True        
        
        #---TAPAS params
        self.bush_gap = 1e-2
        self.pas_cost_mu = 1e-2     
        self.pas_cost_epsilon = 1e-2
        self.pas_flow_mu = 1e-2       
        self.line_search_gap = 1E-8
        
        self.resetPAS()
        
        self.flow_epsilon = 0.000001
        self.min_line_search_gap = 1E-8
        self.tapas_equilibrate_iter = 3
        
        #---used within TAPAS don't change        
        self.good_pas_cost_mu = 0
        self.good_pas_flow_mu = 0
        self.good_bush_gap = 0
        self.good_pas_cost_epsilon = 0
        
        #---BB params
        self.CPLEX_threads = 1
        self.BB_timelimit = 3600
        self.BB_tol = 1E-2
        
        #---BPC / BC params
        self.OAcut_tol = 0.05
        self.solveSO = False
        self.min_gap_SO_OA_cuts = 1E-1 
        self.runUEifCGIntegral = True
        self.useInterdictionCuts = True
        self.useValueFunctionCuts1 = False
        self.useValueFunctionCuts2 = False
        #self.initOAheuristic = 'kBestKNP'
        self.initOAheuristic = 'LocalSearchKNP'
        #self.initOAheuristic = 'LocalSearchY1'
        
        #---printing    
        self.DEBUG_CHECKS = False

        self.PRINT_PAS_INFO = False
        self.PRINT_BRANCH_INFO = False
        self.PRINT_TAPAS_INFO = False
        
        self.PRINT_TAP_ITER = False

        self.printBushEquilibrate = False
        self.printReducedCosts = False        
        self.PRINT_PARAM_ADJ = False
        self.PRINT_PAS_DEBUG = False
        
        self.PRINT_BB_INFO = False #---prints detailed BB info
        self.PRINT_BB_BASIC = False #---prints only basic BB info        
        self.PRINT_BB_NODE_DETAIL = False
        self.PRINT_BB_NODE = False

        self.PRINT_CG_INFO = False
        self.PRINT_LOG = False #---outputs instance log file - used in BPC only for now
        
        
        self.min_CNDP_gap = 1e-2
        
        self.OA_TOL_X = 2e-2
        self.OA_TOL_Y = 2e-2
        self.VF_TOL_X = 1e-2
        
        
        self.BACKTRACK = False
        self.BACKTRACK_INTERVAL = 10
        self.BACKTRACK_LEVEL_INTERVAL = 100