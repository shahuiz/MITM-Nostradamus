deployment = {
    "BiD": "ON", 
    "GnD": "OFF", 
    "E_Time": 1e6, 
    'mtype': 'MC',
    'LB': 50, 
    'UB': 90, 
    'Quantum': True
    }

config = {
   "keylen": 128
    , "tot": 7
    , "enc": 2
    , "mat": 5
    , "ksc": 0
    }

if deployment["Quantum"]:
    folder = './AES/qruns/'
else: 
    folder = './AES/runs/'

import os
import time
from AES import solver_AES

if not os.path.exists(folder):
    os.mkdir(folder)


start_time = time.time()

ModelName, result = solver_AES(
    KeyLen=config["keylen"], 
    TOTAL=config["tot"], 
    ENCST=config["enc"], 
    MATCH=config["mat"], 
    KEYST=config["ksc"], 
    control_panel=deployment, 
    dir=folder
    )

time_cost = time.time() - start_time

        