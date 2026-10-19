import math, random

CL_L, CL_R, CL_T, CL_B = 300.0, 740.0, 280.0, 440.0
DOCK = (520.0, CL_T + 10)
SERT_POS = [(380.0, CL_T + 5), (660.0, CL_T + 5)]
RECS = [dict(x=420.0, y=CL_B-5, aff=3.0), dict(x=540.0, y=CL_B-5, aff=0.8), dict(x=620.0, y=CL_B-5, aff=0.25)]
AUTO = (710.0, CL_T + 15)
DT = 1/60.0

def run(D, SERT_R, REC_R, k_off, EDGE_REFLECT, block, stress, desens, adapt,
        seed=1, t_end=400.0, burn=150.0):
    rng = random.Random(seed)
    caInt = (2.2 - 1.2*stress) * (1 + 0.8*block*(1-desens))
    sens = [(1.0-0.6*stress)*(1+0.5*adapt), (1.0+0.8*stress)*(1-0.45*adapt), (1-0.3*adapt)]
    autoSens = (1.0 - 0.8*stress) * (1 - desens)
    # per-SERT persistent blocked state, reroll every 10s (as in file)
    blocked = [False, False]; sertTimer = 0.0
    mol = []  # [x,y,boundT,recIdx] recIdx -1 free, -2 auto, 0..2 receptor
    occ = [None,None,None]; autoOccSlot = None; activeUntil = 0.0
    caT = 0.0; caActive=False; caLeft=0.0; fuse=None
    ves=[2,2,2]; rot=0; busy=-1
    samples=[]; relHz=[]; upHz=[]
    win_rel=0; win_up=0; winT=0.0
    sd = D*math.sqrt(DT)
    for i in range(int(t_end/DT)):
        t=i*DT
        sertTimer+=DT
        if sertTimer>10: 
            sertTimer=0
            blocked=[rng.random()<block for _ in range(2)]
        # release machinery (mirrors file): ca events every caInt, fuse when a full vesicle is active
        caT+=DT
        gate = caInt*(0.8+rng.random()*0.4)
        if not caActive and caT>gate:
            caT=0; caActive=True; caLeft=0.6
        if caActive:
            caLeft-=DT
            if caLeft<=0: caActive=False
        busyAvail = any(v>=4 for v in ves)
        if caActive and busy==-1 and busyAvail:
            busy=rng.randrange(3); fuse=0.0
        if busy!=-1:
            fuse+=DT
            if fuse>0.9:  # docking+fusing
                k=ves[busy]; ves[busy]=0; busy=-1
                for _ in range(k):
                    mol.append([DOCK[0]+(rng.random()-0.5)*40, DOCK[1]+rng.random()*12, None, -1])
                win_rel+=k
        # refill packaging ~ synthesis limited; approximate rate: 1 mol / (1/synthRate) with synthRate≈1.1/s calm
        # use stochastic packing: each ves gains 1 mol with prob p per sec (shared pool from SERT recycling implicit)
        if rng.random() < 1.1*DT:
            j=min(range(3), key=lambda x: ves[x])
            if ves[j]<4: ves[j]+=1
        newmol=[]
        for m in mol:
            if m[3]==-2:  # auto held
                m[2]-=DT
                if m[2]<=0:
                    autoOccSlot=None
                    activeUntil=t+3.0*autoSens
                    continue  # dies
                newmol.append(m); continue
            if m[3]>=0:   # bound to post-rec
                m[2]-=DT
                if m[2]<=0:
                    occ[m[3]]=None; m[3]=-1; m[2]=None
                newmol.append(m); continue
            m[0]+=rng.gauss(0,1)*sd; m[1]+=rng.gauss(0,1)*sd
            if m[1]<CL_T+4: m[1]=CL_T+4
            if m[1]>CL_B-4: m[1]=CL_B-4
            if m[0]<CL_L or m[0]>CL_R:
                if rng.random()<EDGE_REFLECT:
                    m[0]=CL_L+4 if m[0]<CL_L else CL_R-4
                else: continue
            cap=False
            for si,(sx,sy) in enumerate(SERT_POS):
                if not blocked[si] and math.hypot(m[0]-sx,m[1]-sy)<SERT_R:
                    win_up+=1; cap=True; break
            if cap: continue
            for ri,r in enumerate(RECS):
                if occ[ri] is not None: continue
                if math.hypot(m[0]-r['x'],m[1]-r['y'])<REC_R:
                    if rng.random() < 1-math.exp(-k_off*sens[ri]*DT):
                        occ[ri]=m; m[2]=1.0/k_off; m[3]=ri; break
            if m[3]==-1 and autoOccSlot is None:
                if math.hypot(m[0]-AUTO[0],m[1]-AUTO[1])<25:
                    if rng.random()<autoSens*0.5:
                        autoOccSlot=m; m[2]=2.0; m[3]=-2
        mol=newmol
        if t>burn:
            samples.append(sum(1 for m in mol if m[3]==-1))
        winT+=DT
        if winT>=1:
            winT-=1
            if t>burn: relHz.append(win_rel); upHz.append(win_up)
            win_rel=0; win_up=0
    ss=sorted(samples)
    f=lambda l: sum(l)/len(l)
    return dict(mean=f(samples), p10=ss[len(ss)//10], med=ss[len(ss)//2], p90=ss[9*len(ss)//10],
                rel=f(relHz), up=f(upHz))

if __name__=='__main__':
    import sys
    D=float(sys.argv[1]); SR=float(sys.argv[2]); RR=float(sys.argv[3]); KO=float(sys.argv[4]); ER=float(sys.argv[5])
    PH={'CTRL':(0,0,0),'F1':(0.92,0,0),'F2':(0.92,1,0),'F3':(0.92,1,1)}
    for st_name,st in (('calm',0.0),('strs',1.0)):
        out=[]
        for name,(b,d,a) in PH.items():
            rs=[run(D,SR,RR,KO,ER,b,st,d,a,seed=s) for s in (1,2,3)]
            m=sum(r['mean'] for r in rs)/3; p10=sum(r['p10'] for r in rs)/3; p90=sum(r['p90'] for r in rs)/3
            rel=sum(r['rel'] for r in rs)/3; up=sum(r['up'] for r in rs)/3
            out.append(f"{name}:{m:.0f}[{p10:.0f}-{p90:.0f}] rel={rel:.2f} up={up:.2f}")
        print(st_name, ' '.join(out))
