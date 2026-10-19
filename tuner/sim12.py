# -*- coding: utf-8 -*-
# Броуновская реплика synapse_v8.3.html (правка 12, новая логика: контроль = хроническая тревога)
import random, math, sys

CL,CR,CT,CB=300,740,280,440
D=75.0; BIND=0.49; EDGE=0.96; MAXPOP=28
PRE_CX,PRE_CY,PRE_RX,PRE_RY=440,150,420,130
VES_FULL=4
REV=dict(x=520,y=115); REV_R=70; SLOTS=[90,210,330]
DOCK=dict(x=520,y=CT-12)
SYN=dict(trp=(220,60),branch=(250,100),htp=(360,160),ht=(400,200))
KYN=dict(kyn=(140,180),kyna=(80,190),quin=(170,210))
MAO=dict(x=770,y=150)
SERT=[dict(x=380,y=CT+5,radius=19,blocked=False),dict(x=660,y=CT+5,radius=19,blocked=False)]
AUTO=dict(x=710,y=CT)
POST=[dict(x=420,y=CB-5,type='1A',aff=3.0,radius=21,occ=None),
      dict(x=540,y=CB-5,type='2A',aff=0.8,radius=21,occ=None),
      dict(x=620,y=CB-5,type='3', aff=0.25,radius=21,occ=None)]
FIXED_TARGETS=[[6,17],[17,28],[20,28],[20,28]]
PHASE_BIAS=[0,+3,+1,+2]
RELPROB_BASE=[0.75,0.85,0.9,0.9]
TPH=[dict(block=0,desens=0,adapt=0),dict(block=.92,desens=0,adapt=0),
     dict(block=.92,desens=1,adapt=0),dict(block=.92,desens=1,adapt=1)]
STRESS_BASE=0.85

class Mol:
    __slots__=('x','y','vx','vy','state','t','life','degradation','target','dead')
    def __init__(s,x,y,state):
        s.x,s.y,s.vx,s.vy=x,y,0,0; s.state=state; s.t=0; s.life=0
        s.degradation=False; s.target=None; s.dead=False

def barrelPos(i,rot):
    a=(SLOTS[i]+rot*120)*math.pi/180
    return dict(x=REV['x']+REV_R*math.cos(a), y=REV['y']+REV_R*math.sin(a))
def activeIdx(rotation): r=int(round(rotation))%3; return (3-r)%3
def insidePre(x,y,m_):
    dx=(x-PRE_CX)/(PRE_RX-m_); dy=(y-PRE_CY)/(PRE_RY-m_)
    return dx*dx+dy*dy<1
def moveToward(m,tx,ty,sp,dt):
    d=math.hypot(tx-m.x,ty-m.y)
    if d>0.1: m.x+=(tx-m.x)/d*sp*dt; m.y+=(ty-m.y)/d*sp*dt
def wander(m,dt,sp,rng):
    m.vx+=(rng.random()-0.5)*sp*6*dt; m.vy+=(rng.random()-0.5)*sp*6*dt
    v=math.hypot(m.vx,m.vy)
    if v>sp: m.vx=m.vx/v*sp; m.vy=m.vy/v*sp
    m.x+=m.vx*dt; m.y+=m.vy*dt

def run(phase_idx, seed, T=60, warm=60, dt=1/30):
    rng=random.Random(seed)
    tgt=TPH[phase_idx]
    th=dict(block=0,desens=0,adapt=0)
    vs=[dict(count=2,speed=1.0,t=0,tt=0,state='barrel'),
        dict(count=2,speed=0.8,t=0,tt=0,state='barrel'),
        dict(count=2,speed=0.5,t=0,tt=0,state='barrel')]
    mol=[Mol(rng.uniform(290,590),rng.uniform(120,200),'cytosol') for _ in range(8)]
    synTimer=caTimer=0.0; rot=rotdisp=0.0
    occInteg=0.0; occEma=None; relPending=0.0; relProb=RELPROB_BASE[phase_idx]
    autoOcc=None; autoActive=False; autoActiveT=0.0
    calcium=[]; ringRel=upEv=diffEv=maoEv=0.0
    win_n=0; win_sum=0.0; windows=[]; simTime=0.0
    bind={'1A':0,'2A':0,'3':0}
    def P():
        s=STRESS_BASE*(1-0.75*th['adapt'])
        return dict(synth=(0.85+0.15*s)*(1-0.3*th['adapt']),
            caInterval=(0.9-0.35*s)*(1+0.8*th['block']*(1-th['desens'])),
            mao=1.0+0.5*s, autoSens=(0.45+0.55*s)*(1-0.9*th['desens']),
            ht1a=(1.0-0.6*s)*(1+0.5*th['adapt']), ht2a=(1.0+0.8*s)*(1-0.45*th['adapt']),
            ht3=(1-0.3*th['adapt']), maoFrac=0.45+0.3*s, stress=s)
    def cleftCount(): return sum(1 for m in mol if m.state=='cleft' and not m.dead)
    for step in range(int((T+warm)/dt)):
        simTime+=dt
        for k in ('block','desens','adapt'):
            d=tgt[k]-th[k]; stp=dt/5
            th[k]+= d if abs(d)<stp else (stp if d>0 else -stp)
        prm=P()
        n=cleftCount()
        if occEma is None: occEma=n
        tauAU=10.0
        occEma+=(n-occEma)*min(1,dt/tauAU)
        relPending*=max(0.0,1-dt/tauAU)
        r=FIXED_TARGETS[phase_idx]
        center=max(r[0],min(r[1],0.5*(r[0]+r[1])+PHASE_BIAS[phase_idx]))
        err=(occEma+relPending)-center
        occInteg+=err*dt
        raw=1+0.35*err+0.25*occInteg
        if (raw>3.0 and occInteg>0) or (raw<0.25 and occInteg<0): occInteg-=err*dt
        occInteg=max(-3,min(3,occInteg))
        gain=max(0.25,min(3.0,raw))
        relProb=max(0.75,min(1.0,RELPROB_BASE[phase_idx]*(0.9+0.1*gain)))
        for s_ in SERT: s_['blocked']=rng.random()<th['block']
        synTimer+=dt
        if synTimer>0.5/prm['synth']:
            synTimer=0; mol.append(Mol(SYN['trp'][0],SYN['trp'][1],'trp'))
        caActive=any(e['active'] and e['t']<e['duration'] for e in calcium)
        caTimer+=dt
        if not caActive and caTimer>prm['caInterval']*(0.8+rng.random()*0.4)/gain:
            caTimer=0; calcium.append(dict(t=0,duration=0.6,active=True))
        rotating=abs(rotdisp-rot)>0.01
        if rotating:
            full=sum(1 for v in vs if v['count']>=VES_FULL)
            diff=rot-rotdisp
            rotdisp+= (1 if diff>0 else -1)*min(abs(diff),(2.5+full*1.5)*dt)
        busy=next((i for i,v in enumerate(vs) if v['state']!='barrel'),-1)
        for v in vs:
            v['t']+=dt
            if v['state']=='docking':
                v['tt']+=dt
                if v['tt']>0.25: v['state']='fusing'; v['t']=0; v['tt']=0
            elif v['state']=='fusing':
                if v['t']>0.45:
                    if rng.random()<relProb:
                        cyto=[mm for mm in mol if mm.state=='cytosol' and not mm.degradation and not mm.dead]
                        while v['count']<VES_FULL and cyto:
                            cyto.pop().dead=True; v['count']+=1
                        for k in range(v['count']):
                            mol.append(Mol(DOCK['x']+(rng.random()-0.5)*40, CT+10+rng.random()*12,'cleft'))
                        ringRel+=v['count']; relPending+=v['count']; v['count']=0
                    v['state']='returning'; v['t']=0; v['tt']=0
            elif v['state']=='returning':
                v['tt']+=dt
                if v['tt']>0.25: v['state']='barrel'; v['t']=0; v['tt']=0
        if busy==-1 and not rotating:
            ai=activeIdx(rot); av=vs[ai]
            delay=(2.0 if autoActive else 0.3)/av['speed']
            if av['count']>=VES_FULL and av['t']>delay and caActive:
                av['state']='docking'; av['t']=0; av['tt']=0
            elif av['count']<VES_FULL:
                if any(ix!=ai and vv['count']>=VES_FULL for ix,vv in enumerate(vs)):
                    rot=(rot+1)%3
        cytoN=sum(1 for m in mol if m.state=='cytosol' and not m.dead)
        for m in mol:
            if m.dead: continue
            m.t+=dt; m.life+=dt
            if m.state=='trp':
                moveToward(m,SYN['branch'][0],SYN['branch'][1],80,dt)
                if math.hypot(m.x-SYN['branch'][0],m.y-SYN['branch'][1])<8:
                    ramp=min(1,simTime/12)
                    m.state='kyn' if rng.random()<(0.03+0.32*prm['stress'])*ramp else 'htp'; m.t=0
            elif m.state=='htp':
                moveToward(m,SYN['ht'][0],SYN['ht'][1],80,dt)
                if math.hypot(m.x-SYN['ht'][0],m.y-SYN['ht'][1])<8: m.state='cytosol'; m.t=0; m.life=0
            elif m.state=='kyn':
                moveToward(m,KYN['kyn'][0],KYN['kyn'][1],70,dt)
                if math.hypot(m.x-KYN['kyn'][0],m.y-KYN['kyn'][1])<8:
                    m.state='quin' if rng.random()<0.6+0.25*prm['stress'] else 'kyna'; m.t=0
            elif m.state in ('kyna','quin'):
                p=KYN[m.state]; moveToward(m,p[0],p[1],60,dt)
                if math.hypot(m.x-p[0],m.y-p[1])<8: m.dead=True
            elif m.state=='cytosol':
                wander(m,dt,130,rng)
                ex=(m.x-PRE_CX)/PRE_RX; ey=(m.y-PRE_CY)/PRE_RY
                eN=math.hypot(ex,ey); wall=8.0; shell=wall/min(PRE_RX,PRE_RY)
                if eN>1-shell:
                    if eN>1: m.x=PRE_CX+ex/eN*PRE_RX; m.y=PRE_CY+ey/eN*PRE_RY
                    nx,ny=ex/PRE_RX,ey/PRE_RY; nl=math.hypot(nx,ny) or 1; nx/=nl; ny/=nl
                    vn=m.vx*nx+m.vy*ny
                    if vn>0: m.vx-=2*vn*nx; m.vy-=2*vn*ny; m.vx*=0.85; m.vy*=0.85
                    push=40*max(0.0,eN-(1-shell))/shell
                    m.vx-=nx*push*dt*10; m.vy-=ny*push*dt*10
                if not m.degradation and m.life>2.0:
                    rate=(0.15+0.05*max(0,cytoN-15))*prm['mao']
                    if rng.random()<rate*dt: m.degradation=True
                if m.life>8.0: m.dead=True; continue
                if m.degradation:
                    dx=MAO['x']-m.x; dy=MAO['y']-m.y; d=math.hypot(dx,dy) or 1
                    m.vx+=dx/d*500*dt; m.vy+=dy/d*500*dt
                    if d<45: m.dead=True; maoEv+=1
                    continue
                best=None; bc=VES_FULL+1; bi=-1
                for i,v in enumerate(vs):
                    if v['count']<VES_FULL and v['state']=='barrel' and v['count']<bc:
                        bc=v['count']; best=v; bi=i
                if best:
                    pos=barrelPos(bi,rotdisp)
                    dx=pos['x']-m.x; dy=pos['y']-m.y; d=math.hypot(dx,dy) or 1
                    if insidePre(pos['x'],pos['y'],30):
                        m.vx+=dx/d*180*dt; m.vy+=dy/d*180*dt
                        if d<26: best['count']+=1; m.dead=True; continue
            elif m.state=='cleft':
                st=D*math.sqrt(dt)
                m.x+=rng.gauss(0,st); m.y+=rng.gauss(0,st)
                m.y=max(CT+4,min(CB-4,m.y))
                if m.x<CL or m.x>CR:
                    if rng.random()<EDGE: m.x=CL+4 if m.x<CL else CR-4
                    else: m.dead=True; diffEv+=1; continue
                cap=False
                for s_ in SERT:
                    if not s_['blocked'] and math.hypot(m.x-s_['x'],m.y-s_['y'])<s_['radius']:
                        upEv+=1; m.degradation=rng.random()<prm['maoFrac']
                        m.state='cytosol'; m.t=0; m.life=0
                        m.x=s_['x']; m.y=s_['y']-40; m.vx=(rng.random()-0.5)*40; m.vy=-60
                        cap=True; break
                if cap: continue
                for r_ in POST:
                    sens={'1A':prm['ht1a'],'2A':prm['ht2a'],'3':prm['ht3']}[r_['type']]
                    if r_['occ'] is None and math.hypot(m.x-r_['x'],m.y-r_['y'])<r_['radius']:
                        if rng.random()<1-math.exp(-r_['aff']*sens*BIND*dt):
                            r_['occ']=m; bind[r_['type']]+=1
                            m.state='bound'; m.target=r_; m.t=0
                        break
                if m.state=='cleft' and autoOcc is None and math.hypot(m.x-AUTO['x'],m.y-(AUTO['y']+15))<25:
                    if rng.random()<min(1.0,prm['autoSens']*1.2):
                        autoOcc=m; m.state='auto'; m.t=0
            elif m.state=='bound':
                if m.t>1.0/m.target['aff']: m.target['occ']=None; m.state='cleft'; m.t=0
            elif m.state=='auto':
                if m.t>2.0: autoOcc=None; autoActive=True; autoActiveT=3.0*prm['autoSens']; m.dead=True
        frees=[m for m in mol if m.state=='cleft' and not m.dead]
        if len(frees)>MAXPOP:
            frees.sort(key=lambda a:a.life)
            for x_ in frees[:len(frees)-MAXPOP]: x_.dead=True; diffEv+=1
        mol=[m for m in mol if not m.dead]
        if autoActive:
            autoActiveT-=dt
            if autoActiveT<=0: autoActive=False
        for e in calcium: e['t']+=dt
        calcium=[e for e in calcium if e['t']<0.8]
        if simTime>=warm:
            win_n+=1; win_sum+=cleftCount()
            if win_n>=int(10/dt):
                windows.append(dict(t=round(simTime),mean=win_sum/win_n,
                                    relHz=ringRel/(win_n*dt),upHz=upEv/(win_n*dt)))
                win_n=0; win_sum=0.0; ringRel=upEv=0.0
    return windows, bind, dict(diff=diffEv,mao=maoEv)

if __name__=='__main__':
    phases=[int(x) for x in (sys.argv[1].split(',') if len(sys.argv)>1 else ['0','1','2','3'])]
    seeds=[int(x) for x in (sys.argv[2].split(',') if len(sys.argv)>2 else ['1','2','3'])]
    T=float(sys.argv[3]) if len(sys.argv)>3 else 60
    WARM=float(sys.argv[4]) if len(sys.argv)>4 else 60
    for ph in phases:
        allmeans=[]; rels=[]; ups=[]; btot={'1A':0,'2A':0,'3':0}; tb=0
        okcnt=wcnt=0; per_seed=[]
        for sd in seeds:
            wins,bind,loss=run(ph,sd,T=T,warm=WARM)
            ms=[w['mean'] for w in wins]
            per_seed.append((sd,[round(x,1) for x in ms]))
            allmeans+=ms; rels+=[w['relHz'] for w in wins]; ups+=[w['upHz'] for w in wins]
            for k in btot: btot[k]+=bind[k]
            tb+=sum(bind.values())
            r=FIXED_TARGETS[ph]
            for w in wins:
                okcnt+= 1 if r[0]<=w['mean']<=r[1] else 0; wcnt+=1
        mu=sum(allmeans)/len(allmeans)
        sd_=math.sqrt(sum((x-mu)**2 for x in allmeans)/len(allmeans))
        share={k:(btot[k]/tb*100 if tb else 0) for k in btot}
        print(f"F{ph}: mean={mu:.1f}+-{sd_:.1f} [{min(allmeans):.1f}..{max(allmeans):.1f}] цель={FIXED_TARGETS[ph]} "
              f"OK {okcnt}/{wcnt} | выброс={sum(rels)/len(rels):.2f}мол/с SERT={sum(ups)/len(ups):.2f}/с "
              f"| связ: 2A={share['2A']:.0f}% 1A={share['1A']:.0f}% 3={share['3']:.0f}%")
        for sd,ms in per_seed: print(f"   seed{sd}: {ms}")
