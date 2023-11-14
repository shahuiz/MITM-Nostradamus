deployment = {
    "BiD": "ON",
    "GnD": "OFF",
    "E_Time": 1e6,
    'mtype': 'MC', 
    'LB': 200,
    'Quantum': True
    }

config = {
     "tot": 6
    ,"enc": 2
    ,"mat": 0
    ,"ksc": 0
    }

folder = './Whirlpool/runs/'

import os
import time
from Whirlpool import solver_Whirlpool

if not os.path.exists(folder):
    os.mkdir(folder)


start_time = time.time()

ModelName, result = solver_Whirlpool(
    TOTAL=config["tot"], 
    ENCST=config["enc"], 
    MATCH=config["mat"], 
    KEYST=config["ksc"], 
    control_panel=deployment, 
    dir=folder
    )

time_cost = time.time() - start_time