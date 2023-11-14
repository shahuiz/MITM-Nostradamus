import numpy as np
import math
import gurobipy as gp
from gurobipy import GRB

# Whirlpool parameters
NROW = 8
Nb = 8
Nk = 8
LeftShift=[0,1,2,3,4,5,6,7]

NBRANCH = NROW + 1
NBYTE = 8
NWORD = 32
EPS = 1e-6
SCALE = 2

def XOR_rule(m: gp.Model, in1_b, in1_r, in1_w, in2_b, in2_r, in2_w, out_b, out_r, out_w, cost_fwd, cost_bwd):
    IEQ_A = np.asarray(
        [[0, 1, 0, 1, 0, -1, 0], [0, 0, 0, 0, -1, -1, -1], [0, -1, 0, 0, 0, 1, 0], [0, 0, 0, -1, 0, 1, 0], [-1, -1, 0, 0, 1, 1, 1], [0, 0, -1, -1, 1, 1, 1], [1, 0, 1, 0, -1, 0, -1]]
    )

    IEQ_B = np.asarray(
        [0, 1, 0, 0, 0, 0, 0]
    )
    
    enum_fwd = [in1_b, in1_w, in2_b, in2_w, out_b, out_w, cost_fwd]
    m.addMConstr(IEQ_A, list(enum_fwd), '>=', -IEQ_B)

    enum_bwd = [in1_r, in1_w, in2_r, in2_w, out_r, out_w, cost_bwd]
    m.addMConstr(IEQ_A, list(enum_bwd), '>=', -IEQ_B)
    pass

def SBox_rule(m: gp.Model, in_b, in_r, in_w, out_b, out_r, out_w):
    IEQ_A = np.asarray(
        [[0, 0, -1, 0, 0, 1], [0, 0, 0, -1, -1, -1]]
    )

    IEQ_B = np.asarray(
        [0, 1]
    )

    EQ_A = np.asarray(
        [[0, 1, 1, 0, -1, -1], [1, 0, 1, -1, 0, -1]]
    )
    EQ_B = np.asarray([0, 0])

    enum = [in_b, in_r, in_w, out_b, out_r, out_w]
    m.addMConstr(IEQ_A, list(enum), '>=', -IEQ_B)
    m.addMConstr(EQ_A, list(enum), '==', -EQ_B)

def SubBytes(m: gp.Model, Nb, xb, xr, xw, yb, yr, yw):
    for j in range(Nb):
        for i in range(NROW):
            SBox_rule(m, xb[j,i],xr[j,i],xw[j,i], yb[j,i],yr[j,i],yw[j,i])

def AddRoundKey(m: gp.Model, Nb, xb, xr, xw, kb, kr, kw, yb, yr, yw, cost_fwd, cost_bwd):
    for j in range(Nb):
        for i in range(NROW):
            XOR_rule(m, xb[j,i],xr[j,i],xw[j,i], kb[j,i],kr[j,i],kw[j,i], yb[j,i],yr[j,i],yw[j,i], cost_fwd[j,i],cost_bwd[j,i])

def MixColumns(m: gp.Model, Nb, xb, xr, xw, Eb, Er, Ew, yb, yr, yw, cost_fwd, cost_bwd):
    for j in range(Nb):
        # column-wise encoders: Eb, Er, Ew
        m.addConstr(Ew[j] == gp.max_(xw[j].tolist()))
        m.addConstr(Eb[j] == gp.max_(xb[j].tolist()))
        m.addConstr(Er[j] == gp.max_(xr[j].tolist()))
        # on yw
        m.addConstr(gp.quicksum(yw[j]) == NROW*Ew[j])
        # on yb
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(xb[j]) >= (NROW+1)*(Eb[j]-Ew[j]))
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(cost_fwd[j]) <= NROW*(1-Ew[j]))
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(cost_fwd[j]) <= NROW*Eb[j])
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(cost_fwd[j]) >= NROW*(Eb[j]-Ew[j]))
        # on yr
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(xr[j]) >= (NROW+1)*(Er[j]-Ew[j]))
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(cost_bwd[j]) <= NROW*(1-Ew[j]))
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(cost_bwd[j]) <= NROW*Er[j])
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(cost_bwd[j]) >= NROW*(Er[j]-Ew[j]))
    return

def GnDMixColumns(m: gp.Model, Nb, xb, xr, xw, Gb, Gr, Gw, Gbr, Eb, Er, Ew, yb, yr, yw, cost_fwd, cost_bwd, GnD_flag=False):
    if GnD_flag == False:
        MixColumns(m, Nb, xb, xr, xw, Eb, Er, Ew, yb, yr, yw, cost_fwd, cost_bwd)
        return
    
    for j in range(Nb):
        for i in range(NROW):
            m.addConstr(xw[j,i] == Gw[j,i]+Gb[j,i]+Gr[j,i]+Gbr[j,i])

        # column-wise encoders: Eb, Er, Ew
        m.addConstr(Ew[j] == gp.max_(Gw[j].tolist()))
        m.addConstr(Eb[j] == gp.max_(xb[j].tolist() + Gb[j].tolist() + Gbr[j].tolist()))
        m.addConstr(Er[j] == gp.max_(xr[j].tolist() + Gr[j].tolist() + Gbr[j].tolist()))

        # on yw
        m.addConstr(gp.quicksum(yw[j]) == NROW*Ew[j])
        # on yb
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(xb[j]) + gp.quicksum(Gb[j]) + gp.quicksum(Gbr[j]) >= (NROW+1)*(Eb[j]-Ew[j]))
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(cost_fwd[j]) <= NROW*(1-Ew[j]))
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(cost_fwd[j]) <= NROW*Eb[j])
        m.addConstr(gp.quicksum(yb[j]) + gp.quicksum(cost_fwd[j]) >= NROW*(Eb[j]-Ew[j]))
        # on yr
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(xr[j]) + gp.quicksum(Gr[j]) + gp.quicksum(Gbr[j]) >= (NROW+1)*(Er[j]-Ew[j]))
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(cost_bwd[j]) <= NROW*(1-Ew[j]))
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(cost_bwd[j]) <= NROW*Er[j])
        m.addConstr(gp.quicksum(yr[j]) + gp.quicksum(cost_bwd[j]) >= NROW*(Er[j]-Ew[j]))
    return

def MatchThruSB(m: gp.Model, Nb, xw, yw, meet, dom):
    for j in range(Nb):
        for i in range(NROW):
            m.addConstr(meet[j,i] + xw[j,i] <= 1)
            m.addConstr(meet[j,i] + yw[j,i] <= 1)
    m.addConstr(dom == gp.quicksum(meet.flatten()))
    return

def MatchThruARK(m: gp.Model, Nb, xw, yw, kw, meet, dom):
    for j in range(Nb):
        for i in range(NROW):
            m.addConstr(meet[j,i] + xw[j,i] <= 1)
            m.addConstr(meet[j,i] + yw[j,i] <= 1)
            m.addConstr(meet[j,i] + kw[j,i] <= 1)
    m.addConstr(dom == gp.quicksum(meet.flatten()))
    return

def MatchThruMC(m: gp.Model, Nb, MC_b, MC_r, MC_w, AK_b, AK_r, AK_w, meet, dom):
    meet_n = np.asarray(m.addVars(Nb, lb=-NROW, ub=NROW, vtype = GRB.INTEGER, name='match_case_1_signed').values())
    
    for j in range(Nb):
        m.addConstr(meet_n[j] == NROW - gp.quicksum(MC_w[j].flatten()) - gp.quicksum(AK_w[j].flatten()))
        m.addConstr(meet[j] == gp.max_(meet_n[j], 0))
    m.addConstr(dom == SCALE*gp.quicksum(meet.flatten()))
    return

def solver_Whirlpool(TOTAL:int, ENCST:int, MATCH:int, KEYST:int, control_panel, dir):
    #BlkLen = 128

    # define optimization model
    m = gp.Model('Whirlpool-%d-%d-%d-%d' % (TOTAL, ENCST, MATCH, KEYST))
    
    # Calculate Nk for key schedule
    KSTOTAL = math.ceil((TOTAL + 1) * Nb / Nk)
    
    # assign forward and backward rounds, excluding match round and last round
    if ENCST < MATCH:
        fwd = list(range(ENCST, MATCH))
        bwd = list(range(MATCH+1, TOTAL+1)) + list(range(0, ENCST))
    else:
        bwd = list(range(MATCH+1, ENCST))
        fwd = list(range(ENCST, TOTAL+1)) + list(range(0, MATCH))

    
    ### State Vars
    # define SR, MC, AK, SB states
    SR_b = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='SR_b').values()).reshape((TOTAL, Nb, NROW))
    SR_r = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='SR_r').values()).reshape((TOTAL, Nb, NROW))
    SR_w = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='SR_w').values()).reshape((TOTAL, Nb, NROW))

    MC_b = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='MC_b').values()).reshape((TOTAL, Nb, NROW))
    MC_r = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='MC_r').values()).reshape((TOTAL, Nb, NROW))
    MC_w = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='MC_w').values()).reshape((TOTAL, Nb, NROW))

    SB_b = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='SB_b').values()).reshape((TOTAL, Nb, NROW))
    SB_r = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='SB_r').values()).reshape((TOTAL, Nb, NROW))
    SB_w = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='SB_w').values()).reshape((TOTAL, Nb, NROW))

    # define Guess-and_Determine switches
    Gb = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='Gb').values()).reshape((TOTAL, Nb, NROW))
    Gr = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='Gr').values()).reshape((TOTAL, Nb, NROW))
    Gw = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='Gw').values()).reshape((TOTAL, Nb, NROW))   
    Gbr = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype= GRB.BINARY, name='Gbr').values()).reshape((TOTAL, Nb, NROW))

    # define MC column-wise indicators
    Eb = np.asarray(m.addVars(TOTAL, Nb, vtype= GRB.BINARY, name='Eb').values()).reshape((TOTAL, Nb))
    Er = np.asarray(m.addVars(TOTAL, Nb, vtype= GRB.BINARY, name='Er').values()).reshape((TOTAL, Nb))
    Ew = np.asarray(m.addVars(TOTAL, Nb, vtype= GRB.BINARY, name='Ew').values()).reshape((TOTAL, Nb))

    # define tag
    Tag_b = np.asarray(m.addVars(Nb, NROW, vtype= GRB.BINARY, name='Tag_b').values()).reshape((Nb, NROW))
    Tag_r = np.asarray(m.addVars(Nb, NROW, vtype= GRB.BINARY, name='Tag_r').values()).reshape((Nb, NROW))
    Tag_w = np.asarray(m.addVars(Nb, NROW, vtype= GRB.BINARY, name='Tag_w').values()).reshape((Nb, NROW))

    AT_b = np.asarray(m.addVars(Nb, NROW, vtype= GRB.BINARY, name='AT_b').values()).reshape((Nb, NROW))
    AT_r = np.asarray(m.addVars(Nb, NROW, vtype= GRB.BINARY, name='AT_r').values()).reshape((Nb, NROW))
    AT_w = np.asarray(m.addVars(Nb, NROW, vtype= GRB.BINARY, name='AT_w').values()).reshape((Nb, NROW))


    ### Cost Vars
    # define auxiliary vars tracking cost of df at MC operations
    mc_cost_fwd = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype=GRB.BINARY, name='mc_cost_fwd').values()).reshape((TOTAL, Nb, NROW))
    mc_cost_bwd = np.asarray(m.addVars(TOTAL, Nb, NROW, vtype=GRB.BINARY, name='mc_cost_bwd').values()).reshape((TOTAL, Nb, NROW))
    
    # define auxiliary vars tracking cost of df at Add Key operations in foward direction, at the LHS and the RHS of MC gate
    xor_cost_fwd = np.asarray(m.addVars(TOTAL+1, Nb, NROW, vtype= GRB.BINARY, name='xor_cost_fwd').values()).reshape((TOTAL+1, Nb, NROW))
    xor_cost_bwd = np.asarray(m.addVars(TOTAL+1, Nb, NROW, vtype= GRB.BINARY, name='xor_cost_bwd').values()).reshape((TOTAL+1, Nb, NROW))

    # define auxiliary vars trackin cost of df in the key expansion process, unpossible combinations are set to zeros
    key_cost_fwd = np.asarray(m.addVars(KSTOTAL, Nk, NROW, vtype= GRB.BINARY, name='key_cost_fwd').values()).reshape((KSTOTAL, Nk, NROW))
    key_cost_bwd = np.asarray(m.addVars(KSTOTAL, Nk, NROW, vtype= GRB.BINARY, name='key_cost_bwd').values()).reshape((KSTOTAL, Nk, NROW))
    
    # define match var
    MC_match = np.asarray(m.addVars(Nb, lb=0, ub=NROW, vtype=GRB.INTEGER, name='MC_match').values()).reshape(Nb)
    AK_match = np.asarray(m.addVars(Nb, NROW, vtype=GRB.BINARY, name='AK_match').values()).reshape((Nb, NROW))

    ### Obj vars
    DOF_b = m.addVar(lb=1, vtype=GRB.INTEGER, name="DF_b")
    DOF_r = m.addVar(lb=1, vtype=GRB.INTEGER, name="DF_r")
    DOM = m.addVar(lb=1, vtype=GRB.INTEGER, name="Match")
    
    GnD_b = m.addVar(lb=0, vtype=GRB.INTEGER, name="GND_b")
    GnD_r = m.addVar(lb=0, vtype=GRB.INTEGER, name="GND_r")
    GnD_br = m.addVar(lb=0, vtype=GRB.INTEGER, name="GND_br")

    d1 = m.addVar(lb=1, vtype=GRB.INTEGER, name="d1")
    d2 = m.addVar(lb=1, vtype=GRB.INTEGER, name="d2")
    dm = m.addVar(lb=1, vtype=GRB.INTEGER, name="dm") 
    dq = m.addVar(lb=1, vtype=GRB.INTEGER, name="dq")
    
    OBJ = m.addVar(lb=1, vtype=GRB.INTEGER, name="Obj")
    freeTag = m.addVar(lb=0, vtype=GRB.INTEGER, name="freeTag")
    blrdTag = m.addVar(lb=0, vtype=GRB.INTEGER, name="blrdTag")
    
    kw = m.addVar(lb=0, vtype=GRB.INTEGER, name="kw")
    kbr = m.addVar(lb=0, vtype=GRB.INTEGER, name="kbr")


    MITM_cplx = m.addVar(lb=4, vtype=GRB.CONTINUOUS, name="MITM_cplx")
    Diamond_cplx = m.addVar(lb=4, vtype=GRB.CONTINUOUS, name="Diamond_cplx")
    Cplx = m.addVar(lb=4, vtype=GRB.CONTINUOUS, name="Cplx")

    #### Main Procedure ####
    print('\n\n\n' + m.modelName)
    print('Parameters: Nb=%d, Nk=%d, keyexp_r=%d\n' %(Nb, Nk, KSTOTAL))
    print('Left Shift Cosntant: ', LeftShift)
    print('Round structure:')


    for r in range(TOTAL):
        for j in range(Nb):
            for i in range(NROW):
                m.addConstr(xor_cost_fwd[r,j,i] + xor_cost_bwd[r,j,i] == 0)

    # initialize start states and tag, allow no white cells
    for j in range(Nb):
        for i in range(NROW):
            m.addConstr(SR_w[ENCST,j,i] == 0)
            m.addConstr(SR_b[ENCST,j,i] + SR_r[ENCST,j,i] <= 1)

            # only allow grey or white
            m.addConstr(Tag_b[j,i] + Tag_r[j,i] == 0)

    # bulk perform SR
    for r in range(TOTAL):
        for j in range(Nb):
            for i in range(NROW):
                sj = (j + LeftShift[i]) % Nb
                m.addConstr(MC_b[r,j,i] == SR_b[r,sj,i])
                m.addConstr(MC_r[r,j,i] == SR_r[r,sj,i])
                m.addConstr(MC_w[r,j,i] == SR_w[r,sj,i])

    # fix no MixColumn cost in the last round
    for j in range(Nb):
        for i in range(NROW):
            m.addConstr(mc_cost_fwd[TOTAL-1,j,i] == 0)
            m.addConstr(mc_cost_bwd[TOTAL-1,j,i] == 0)

    # add constriants according to the encryption algorithm
    for r in range(TOTAL):
        nr, lr = r+1, r-1   # alias for next round
        
        # special matching: identity meet at last round, only allows same pure colored cells to match 
        if r == MATCH and MATCH == TOTAL - 1:
            print('mat lastr', r)
            SubBytes(m, Nb, SR_b[0], SR_r[0], SR_w[0], SB_b[r], SB_r[r], SB_w[r])
            GnDMixColumns(m, Nb, MC_b[r], MC_r[r], MC_w[r], Gb[r], Gr[r], Gw[r], Gbr[r], Eb[r], Er[r], Ew[r], AT_b, AT_r, AT_w, mc_cost_fwd[r], mc_cost_bwd[r], GnD_flag=True)    
            MatchThruARK(m, Nb, SB_w[r], AT_w[r], Tag_w, AK_match, DOM)

        elif r == MATCH:
            print('mat', r)
            SubBytes(m, Nb, SR_b[nr], SR_r[nr], SR_w[nr], SB_b[r], SB_r[r], SB_w[r])
            MatchThruMC(m, Nb, MC_b[r], MC_r[r], MC_w[r], SB_b[r], SB_r[r], SB_w[r], MC_match, DOM)

        # last round
        elif r == TOTAL - 1 and r in fwd:
            print('fwd lastr', r)
            GnDMixColumns(m, Nb, MC_b[r], MC_r[r], MC_w[r], Gb[r], Gr[r], Gw[r], Gbr[r], Eb[r], Er[r], Ew[r], AT_b, AT_r, AT_w, mc_cost_fwd[r], mc_cost_bwd[r], GnD_flag=True)    
            AddRoundKey(m, Nb, AT_b, AT_r, AT_w, Tag_b, Tag_r, Tag_w, SB_b[r], SB_r[r], SB_w[r], xor_cost_fwd[r], xor_cost_bwd[r])
            SubBytes(m, Nb, SB_b[r], SB_r[r], SB_w[r], SR_b[0], SR_r[0], SR_w[0])
        
        elif r == TOTAL - 1 and r in bwd:
            print('bwd lastr', r)
            SubBytes(m, Nb, SR_b[0], SR_r[0], SR_w[0], SB_b[r], SB_r[r], SB_w[r])
            AddRoundKey(m, Nb, SB_b[r], SB_r[r], SB_w[r], Tag_b, Tag_r, Tag_w, AT_b, AT_r, AT_w, xor_cost_fwd[r], xor_cost_bwd[r])
            GnDMixColumns(m, Nb, AT_b, AT_r, AT_w, Gb[r], Gr[r], Gw[r], Gbr[r], Eb[r], Er[r], Ew[r], MC_b[r], MC_r[r], MC_w[r], mc_cost_fwd[r], mc_cost_bwd[r], GnD_flag=True)    

        # ENC Propagation: forward direction
        elif r in fwd:
            print('fwd', r)
            GnDMixColumns(m, Nb, MC_b[r], MC_r[r], MC_w[r], Gb[r], Gr[r], Gw[r], Gbr[r], Eb[r], Er[r], Ew[r], SB_b[r], SB_r[r], SB_w[r], mc_cost_fwd[r], mc_cost_bwd[r], GnD_flag=True)    
            SubBytes(m, Nb, SB_b[r], SB_r[r], SB_w[r],  SR_b[nr], SR_r[nr], SR_w[nr])
            
        # ENC Propagation: backward direction
        elif r in bwd:
            print('bwd', r)
            SubBytes(m, Nb, SR_b[nr], SR_r[nr], SR_w[nr], SB_b[r], SB_r[r], SB_w[r])
            GnDMixColumns(m, Nb, SB_b[r], SB_r[r], SB_w[r], Gb[r], Gr[r], Gw[r], Gbr[r], Eb[r], Er[r], Ew[r], MC_b[r], MC_r[r], MC_w[r], mc_cost_fwd[r], mc_cost_bwd[r], GnD_flag=True)
    
    # setting objective
    m.addConstr(DOF_b == gp.quicksum(SR_b[ENCST].flatten()) + gp.quicksum(Tag_b.flatten()) - gp.quicksum(mc_cost_fwd.flatten()) - gp.quicksum(xor_cost_fwd.flatten()) - gp.quicksum(key_cost_fwd.flatten()))
    m.addConstr(DOF_r == gp.quicksum(SR_r[ENCST].flatten()) + gp.quicksum(Tag_r.flatten()) - gp.quicksum(mc_cost_bwd.flatten()) - gp.quicksum(xor_cost_bwd.flatten()) - gp.quicksum(key_cost_bwd.flatten()))
    
    m.addConstr(GnD_b == gp.quicksum(Gb.flatten()))
    m.addConstr(GnD_r == gp.quicksum(Gr.flatten()))
    m.addConstr(GnD_br == gp.quicksum(Gbr.flatten()))
    
    if control_panel["GnD"] == 'OFF':
        m.addConstr(GnD_b == 0)
        m.addConstr(GnD_r == 0)
        m.addConstr(GnD_br == 0)

    if control_panel["BiD"] == 'OFF':
        for r in range(TOTAL):
            if r in fwd:
                m.addConstr(gp.quicksum(mc_cost_fwd[r].flatten()) == 0)
            elif r in bwd:
                m.addConstr(gp.quicksum(mc_cost_bwd[r].flatten()) == 0)
        
        for r in range(TOTAL + 1):
            if r in fwd:
                m.addConstr(gp.quicksum(xor_cost_fwd[r].flatten()) == 0)
            elif r in bwd:
                m.addConstr(gp.quicksum(xor_cost_bwd[r].flatten()) == 0)

        for r in range(KSTOTAL):
            if r > KEYST:
                m.addConstr(gp.quicksum(key_cost_fwd[r].flatten()) == 0)
            elif r < KEYST:
                m.addConstr(gp.quicksum(key_cost_bwd[r].flatten()) == 0)

    m.addConstr(blrdTag == gp.quicksum(Tag_b.flatten()) + gp.quicksum(Tag_r.flatten()))
    m.addConstr(freeTag == gp.quicksum(Tag_w.flatten()))
    m.addConstr(freeTag >= 1)

    m.addConstr(kw <= freeTag*NBYTE)
    m.addConstr(kbr == kw + blrdTag*NBYTE)

    m.addConstr(d1 == SCALE*(DOF_b - GnD_r))
    m.addConstr(d2 == SCALE*(DOF_r - GnD_b))
    m.addConstr(dm == DOM - SCALE*(GnD_b + GnD_r + GnD_br))
    
    if control_panel["Quantum"]:
        # Cond for Quantum MITM
        m.addConstr(DOF_b >= DOF_r)
        m.addConstr(dq == DOF_b - DOF_r)

        m.addConstr(OBJ == gp.min_(dq, DOF_b, DOF_r, dm))
        #m.addConstr(OBJ == 6)
        #m.addConstr(kw == 6)
        m.addConstr(OBJ >= SCALE*blrdTag + 1)
        # MITM complexity
        m.addConstr(MITM_cplx == (512 - kw - NBYTE/SCALE*OBJ)/2)
        # Diamond complexity
        m.addConstr(Diamond_cplx == (kbr*2 + 512)/3 + 1)
        # Overall complexity
        m.addConstr(Cplx == gp.max_(MITM_cplx, Diamond_cplx))
    else: 
        m.addConstr(OBJ == gp.min_(d1, d2, dm))
        m.addConstr(OBJ >= SCALE*blrdTag + 1)
        # MITM complexity
        m.addConstr(MITM_cplx == 512 - kw - NBYTE/SCALE*OBJ)
        # Diamond complexity
        m.addConstr(Diamond_cplx == kbr/2 + 512/2)
        # Overall complexity
        m.addConstr(Cplx == gp.max_(MITM_cplx, Diamond_cplx))

    m.addConstr(Cplx >= control_panel["LB"])
    m.setObjective(Cplx, GRB.MINIMIZE)

    # set parameters
    m.setParam(GRB.Param.PoolSearchMode, 2)
    m.setParam(GRB.Param.PoolSolutions,  1)
    m.setParam(GRB.Param.Threads, 4)
    m.setParam(GRB.Param.TimeLimit, control_panel["E_Time"])

    # optimization
    m.optimize()
    
    m.write(dir + m.modelName + '.lp')

    if m.SolCount > 0:
        m.write(dir + m.modelName + '.sol')
        return m.ModelName, m.ObjVal
    else:
        return "", 0
