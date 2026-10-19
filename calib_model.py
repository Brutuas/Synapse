import math, random

# --- geometry & constants copied from synapse_v8.3.html ---
CL_L, CL_R, CL_T, CL_B = 300.0, 740.0, 280.0, 440.0
DOCK = (520.0, CL_T + 10)
SERT_POS = [(380.0, CL_T + 5), (660.0, CL_T + 5)]
SERT_R = 35.0
RECS = [dict(x=420.0, y=CL_B-5, aff=3.0), dict(x=540.0, y=CL_B-5, aff=0.8), dict(x=620.0, y=CL_B-5, aff=0.25)]
AUTO = (710.0, CL_T + 15, 25.0)   # x,y,radius ; capture prob prm.autoSens*0.5 per step
D = 50.0          # DIFFUSION_COEFF px/sqrt(s)
BIND_RATE = 0.5
EDGE_REFLECT = 0.99
DT = 1/60.0
VES_FULL = 4
TAU_BOUND_BASE = 1.0  # bound lifetime = 1.0/aff seconds (as in file: m.t > 1.0/target.aff -> but t resets; ok)

def params(block, stress, desens, adapt):
    p = {}
    p['caInterval'] = (2.2 - 1.2*stress) * (1 + 0.8*block*(1-desens))
    p['ht1a'] = (1.0 - 0.6*stress) * (1 + 0.5*adapt)
    p['ht2a'] = (1.0 + 0.8*stress) * (1 - 0.45*adapt)
    p['ht3']  = (1 - 0.3*adapt)
    p['sens'] = [p['ht1a'], p['ht2a'], p['ht3']]
    p['autoSens'] = (1.0 - 0.8*stress) * (1 - desens)
    return p

def run(block, stress, desens, adapt, seed=1, t_end=300.0, burn=100.0, ves_size=3.0):
    rng = random.Random(seed)
    p = params(block, stress, desens, adapt)
    # SERT blocked state, rerolled every 10 s per file
    blocked = [False, False]
    sertTimer = 0.0
    cleft = []  # [x, y, boundT(None), recIdx]
    occ = [None, None, None]
    autoOcc = None; autoT = 0.0
    activeUntil = 0.0
    caT = 0.0; caActive = False; caLeft = 0.0
    docked = None; fuseT = 0.0
    n_rel = 0; n_up = 0; n_diffloss = 0
    samples = []; rel_samples=[]; up_samples=[]
    steps = int(t_end/DT)
    win_rel = 0; win_up = 0; winT = 0.0
    rel_hist=[]; up_hist=[]
    for i in range(steps):
        t = i*DT
        sertTimer += DT
        if sertTimer > 10.0:
            sertTimer = 0
            blocked = [rng.random() < block for _ in range(2)]
        # release machinery: autoreceptor inhibits (active window doubles interval, as file delay 2.0 vs 0.3)
        caT += DT
        gate = p['caInterval'] * (0.8 + rng.random()*0.4)
        if t < activeUntil: gate *= 2.0
        if not caActive and caT > gate and docked is None:
            caT = 0; caActive = True; caLeft = 0.6
        if caActive:
            caLeft -= DT
            if caLeft <= 0:
                caActive = False
                docked = 'fusing'; fuseT = 0
        if docked == 'fusing':
            fuseT += DT
            if fuseT > 0.45:
                k = int(ves_size)
                for _ in range(k):
                    cleft.append([DOCK[0] + (rng.random()-0.5)*40, DOCK[1] + rng.random()*12, None, -1])
                n_rel += k; win_rel += k
                docked = None
        sd = D*math.sqrt(DT)
        keep = []
        for m in cleft:
            if m[2] is not None:  # bound to receptor
                m[2] -= DT
                if m[2] <= 0:
                    occ[m[3]] = None; m[2] = None; m[3] = -1
                keep.append(m); continue
            m[0] += rng.gauss(0,1)*sd
            m[1] += rng.gauss(0,1)*sd
            if m[1] < CL_T+4: m[1] = CL_T+4
            if m[1] > CL_B-4: m[1] = CL_B-4
            if m[0] < CL_L or m[0] > CL_R:
                if rng.random() < EDGE_REFLECT:
                    m[0] = CL_L+4 if m[0] < CL_L else CL_R-4
                else:
                    n_diffloss += 1; continue
            cap = False
            for si,(sx,sy) in enumerate(SERT_POS):
                if not blocked[si] and math.hypot(m[0]-sx, m[1]-sy) < SERT_R:
                    n_up += 1; win_up += 1; cap = True; break
            if cap: continue
            for ri,r in enumerate(RECS):
                if occ[ri] is not None: continue
                if math.hypot(m[0]-r['x'], m[1]-r['y']) < 20:
                    sens = p['sens'][ri]
                    prob = 1 - math.exp(-r['aff']*sens*BIND_RATE*DT)
                    if rng.random() < prob:
                        occ[ri] = id(m)
                        m[2] = 1.0/r['aff']  # bound lifetime per file
                        m[3] = ri
                        break
            if m[2] is None and autoOcc is None:
                if math.hypot(m[0]-AUTO[0], m[1]-AUTO[1]) < AUTO[2]:
                    if rng.random() < p['autoSens']*0.5:
                        autoOcc = m; m[2] = 2.0; m[3] = -2  # held 2s then triggers inhibition
        # autoreceptor handling
        if autoOcc is not None:
            autoOcc[2] -= DT
            if autoOcc[2] <= 0:
                autoOcc = None
                activeUntil = t + 3.0*p['autoSens']
        # sample counts
        if t > burn:
            samples.append(sum(1 for m in cleft if m[3]==-1 or m[2] is None) + len(cleft)*0 )
            samples[-1] = len(cleft)
        winT += DT
        if winT >= 1.0:
            winT -= 1.0
            if t > burn: rel_hist.append(win_rel); up_hist.append(win_up)
            win_rel = 0; win_up = 0
    ss = sorted(samples)
    mean = sum(samples)/len(samples)
    return dict(mean=mean, med=ss[len(ss)//2], p10=ss[len(ss)//10], p90=ss[9*len(ss)//10],
                relHz=sum(rel_hist)/len(rel_hist), upHz=sum(up_hist)/len(up_hist))

PHASES = {
 'CONTROL': (0.0, 0.0, 0.0),
 'F1':      (0.92, 0.0, 0.0),
 'F2':      (0.92, 1.0, 0.0),
 'F3':      (0.92, 1.0, 1.0),
}
if __name__ == '__main__':
    for name,(b,d,a) in PHASES.items():
        r = run(b, 0.0, d, a, seed=7)
        print(f"{name:8s} calm : mean={r['mean']:.1f} med={r['med']} p10-p90=[{r['p10']},{r['p90']}] rel={r['relHz']:.2f}/s up={r['upHz']:.2f}/s")
    for name,(b,d,a) in PHASES.items():
        r = run(b, 1.0, d, a, seed=7)
        print(f"{name:8s} stress: mean={r['mean']:.1f} med={r['med']} p10-p90=[{r['p10']},{r['p90']}] rel={r['relHz']:.2f}/s up={r['upHz']:.2f}/s")
