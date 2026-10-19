// Headless replica of synapse_v8.3.html update() — physics & controller verbatim.
'use strict';

// ---- tunable params (defaults = current HTML values) ----
const P_ = {
  D: 75, BIND: 0.49, EDGE: 0.96, MAXPOP: 28,
  SERT_R: 19, REC_R: 21, AUTO_R: 25, AUTO_K: 1.2,
  KP: 0.35, KI: 0.25, TAU: 10, SYNTH_MULT: 1.0, CA_MULT: 1.0, GAIN_MAX: 3.0, PH_SYNTH: [1,1,1,1], PH_CA: [1,1,1,1], MAO_BASE: 0.15, CYTO_LIFE: 8.0, TGT_MARGIN: [2,2,2,2], PAR_DOCK: false,
  REL_BASE: [0.75, 0.85, 0.9, 0.9],
  TARGETS_CALM: [[6,15],[19,28],[19,28],[19,28]],
  TARGETS_STR: [[6,17],[17,28],[20,28],[20,28]],
  BIAS: [0, 3, 1, 2],
};
if (process.env.TUNEP) Object.assign(P_, JSON.parse(process.env.TUNEP));

const PRE_CX = 440, PRE_CY = 150, PRE_RX = 420, PRE_RY = 130;
function yMem(x){ const t=(x-PRE_CX)/PRE_RX; return PRE_CY+PRE_RY*Math.sqrt(Math.max(0,1-t*t)); }
const CLEFT = { x: 300, y: 280, w: 440, h: 160 };
const CLEFT_LEFT = CLEFT.x, CLEFT_RIGHT = CLEFT.x + CLEFT.w, CLEFT_TOP = CLEFT.y, CLEFT_BOT = CLEFT.y + CLEFT.h;
const VES_FULL = 4;
const REV_CENTER = {x: 520, y: 115}, REV_RADIUS = 70;
const SLOT_ANGLES = [90, 210, 330];
const DOCK = {x: 520, y: CLEFT_TOP - 12};
const MAO_POS = {x: 770, y: 150, pulse: 0};
const CA_SPAWN = {x: 260, y: yMem(260) - 15};
const SYN = {trp: [220, 60], branch: [250, 100], htp: [360, 160], ht: [400, 200]};
const KYN_PATH = {tdo: [140, 130], kyn: [140, 180], kyna: [80, 190], quin: [170, 210]};

const THERAPY_PHASES = [
  {block:0, desens:0, adapt:0},
  {block:0.92, desens:0, adapt:0},
  {block:0.92, desens:1, adapt:0},
  {block:0.92, desens:1, adapt:1},
];

let molecules = [], calciumEvents = [], cortisolMols = [];
let time = 0, simTime = 0;
let stressLevel = 0, targetStress = 0, tdoRamp = 0;
let targetPhase = 0;
let therapy = {block:0, desens:0, adapt:0};
let revolverRotation = 0, revolverDisplay = 0;
let vesicles = [
  {count: 2, pool: 'RRP', speed: 1.0, t: 0, state: 'barrel', x: 0, y: 0},
  {count: 2, pool: 'Rec', speed: 0.8, t: 0, state: 'barrel', x: 0, y: 0},
  {count: 2, pool: 'Res', speed: 0.5, t: 0, state: 'barrel', x: 0, y: 0},
];
const serts = [
  {x: 380, y: CLEFT_TOP + 5, radius: 0, flash: 0, blocked: false},
  {x: 660, y: CLEFT_TOP + 5, radius: 0, flash: 0, blocked: false},
];
const autoRec = {x: 710, y: CLEFT_TOP, occupied: null, t: 0, active: false, activeT: 0, pulse: 0};
const postRecs = [
  {x: 420, y: CLEFT_BOT - 5, type: '5-HT1A', aff: 3.0, radius: 0, bindings: 0, flash: 0, occ: null, t: 0},
  {x: 540, y: CLEFT_BOT - 5, type: '5-HT2A', aff: 0.8, radius: 0, bindings: 0, flash: 0, occ: null, t: 0},
  {x: 620, y: CLEFT_BOT - 5, type: '5-HT3',  aff: 0.25, radius: 0, bindings: 0, flash: 0, occ: null, t: 0},
];

const RING_N = 10;
let ringIdx = 0, secTimer = 0;
let ring = { rel:Array(RING_N).fill(0), up:Array(RING_N).fill(0), diff:Array(RING_N).fill(0), mao:Array(RING_N).fill(0), lost:Array(RING_N).fill(0) };
function ringTick(dt){
  secTimer += dt;
  if (secTimer >= 1.0){ secTimer -= 1.0; ringIdx = (ringIdx+1)%RING_N; for (let k in ring) ring[k][ringIdx]=0; }
}
function ringAdd(key){ ring[key][ringIdx]++; }

function dist(x1,y1,x2,y2){ return Math.hypot(x2-x1, y2-y1); }
function insidePre(x, y, margin) {
  const dx = (x - PRE_CX) / (PRE_RX - margin), dy = (y - PRE_CY) / (PRE_RY - margin);
  return dx*dx + dy*dy < 1;
}
function barrelPos(i, rot){
  const a = (SLOT_ANGLES[i] + rot) * Math.PI/180;
  return {x: REV_CENTER.x + REV_RADIUS*Math.cos(a), y: REV_CENTER.y + REV_RADIUS*Math.sin(a)};
}
function activeVesicleIdx(rot){
  let best=-1, bd=1e9;
  for (let i=0;i<3;i++){ const p=barrelPos(i,rot); const d=Math.abs(p.y-REV_CENTER.y)+Math.abs(p.x-REV_CENTER.x)*0.001; if(d<bd){bd=d;best=i;} }
  return best;
}
function gaussRandom(){ let u=0,v=0; while(u===0)u=Math.random(); while(v===0)v=Math.random(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); }
class Mol { constructor(x,y,state){ this.x=x;this.y=y;this.state=state;this.vx=0;this.vy=0;this.t=0;this.target=null;this.life=0;this.recycled=false;this.degradation=false;this.trail=[];this.everBound=false;this.dead=false; } }
function moveToward(m,tx,ty,sp,dt){ const dx=tx-m.x,dy=ty-m.y,d=Math.hypot(dx,dy)||1; m.x+=dx/d*sp*dt; m.y+=dy/d*sp*dt; }
function wander(m,dt,speed){ m.vx+=(Math.random()-0.5)*speed*8*dt; m.vy+=(Math.random()-0.5)*speed*8*dt; const v=Math.hypot(m.vx,m.vy); if(v>speed){m.vx*=speed/v;m.vy*=speed/v;} m.x+=m.vx*dt; m.y+=m.vy*dt; }
let synTimer = 0, caTimer = 0, cortTimer = 0;
function spawnTrp(){ molecules.push(new Mol(SYN.trp[0], SYN.trp[1], 'trp')); }

function P(){
  const s = stressLevel;
  return {
    cortisol: 0.1+0.9*s,
    tdo: 0.03+0.32*s,
    synth: (0.85+0.15*s) * (1 - 0.3*therapy.adapt) * P_.SYNTH_MULT * (P_.PH_SYNTH[targetPhase]||1),
    caInterval: (0.9-0.35*s) * (1 + 0.8 * therapy.block * (1 - therapy.desens)) / P_.CA_MULT * (P_.PH_CA[targetPhase]||1),
    sert: 1.0+0.5*s,
    mao: 1.0+0.5*s,
    autoSens: (1.0-0.8*s) * (1 - therapy.desens),
    ht1a: (1.0-0.6*s) * (1 + 0.5*therapy.adapt),
    ht2a: (1.0+0.8*s) * (1 - 0.45*therapy.adapt),
    ht3:  (1 - 0.3*therapy.adapt),
    kynShift: 0.6+0.25*s,
    maoFraction: 0.45+0.3*s
  };
}

let occInteg = 0, occEma = null, relPending = 0, relProb = 0.75;
function cleftCount(){ let n=0; for (let m of molecules) if (m.state==='cleft') n++; return n; }
function normRange(){ return (stressLevel < 0.5 ? P_.TARGETS_CALM : P_.TARGETS_STR)[targetPhase]; }
function targetCenter(){ const r=normRange(); return r[0] + (P_.TGT_MARGIN[targetPhase]!=null?P_.TGT_MARGIN[targetPhase]:2); }
function occController(dt){
  const n = cleftCount();
  if (occEma === null) occEma = n;
  const tauAU = P_.TAU;
  occEma += (n - occEma) * Math.min(1, dt / tauAU);
  relPending *= Math.max(0, 1 - dt / tauAU);
  const tgt = targetCenter();
  const err = (occEma + relPending) - tgt;
  occInteg += err * dt;
  const raw = 1 + P_.KP * err + P_.KI * occInteg;
  if (raw > 3.0 && occInteg > 0 || raw < 0.25 && occInteg < 0) occInteg -= err * dt;
  occInteg = Math.max(-6, Math.min(6, occInteg));
  return Math.max(0.25, Math.min(P_.GAIN_MAX, raw));
}

function resetRadii(){
  serts[0].radius = serts[1].radius = P_.SERT_R;
  postRecs.forEach(r => r.radius = P_.REC_R);
}

function update(dt) {
  time += dt; simTime += dt;
  const tgt = THERAPY_PHASES[targetPhase];
  for (let k of ['block','desens','adapt']) {
    const d = tgt[k] - therapy[k];
    const step = dt / 5;
    therapy[k] += Math.abs(d) < step ? d : Math.sign(d) * step;
  }
  const prm = P();
  if (targetStress > stressLevel) stressLevel = Math.min(targetStress, stressLevel + 0.2*dt);
  else stressLevel = Math.max(targetStress, stressLevel - 0.1*dt);
  if (stressLevel > 0.5) tdoRamp = Math.min(1, tdoRamp + dt/12); else tdoRamp = Math.max(0, tdoRamp - dt/5);
  const tdoEff = prm.tdo * tdoRamp;
  ringTick(dt);
  for (let s of serts) s.blocked = Math.random() < therapy.block;
  let cytosolCount = 0;
  for (let m of molecules) if (m.state === 'cytosol') cytosolCount++;
  synTimer += dt;
  if (synTimer > 0.5 / prm.synth) { synTimer = 0; spawnTrp(); }
  const caActive = calciumEvents.some(e => e.active && e.t < e.duration);
  const occGain = occController(dt);
  relProb = Math.max(0.75, Math.min(1.0, P_.REL_BASE[targetPhase] * (0.9 + 0.1 * occGain)));
  caTimer += dt;
  if (!caActive && caTimer > prm.caInterval * (0.8 + Math.random()*0.4) / occGain) {
    caTimer = 0;
    calciumEvents.push({ x: CA_SPAWN.x+(Math.random()-0.5)*40, y: CA_SPAWN.y+(Math.random()-0.5)*8, t: 0, duration: 0.6, active: true });
  }
  cortTimer += dt;
  if (cortTimer > (3.0 - 2.2*stressLevel)) { cortTimer = 0;
    if (stressLevel > 0.15) cortisolMols.push({ x: CLEFT_LEFT+20+Math.random()*(CLEFT.w-40), y: CLEFT_TOP+20+Math.random()*(CLEFT.h-40), vx:0, vy:0, life:0, phase:'outside' });
  }
  for (let c of cortisolMols) {
    c.life += dt;
    if (c.phase === 'outside') { wander(c, dt, 20); c.y -= 70*dt;
      if (c.y > CLEFT_BOT-6) c.y = CLEFT_BOT-6;
      if (c.x < CLEFT_LEFT+6) c.x = CLEFT_LEFT+6;
      if (c.x > CLEFT_RIGHT-6) c.x = CLEFT_RIGHT-6;
      if (c.y < CLEFT_TOP+6) { c.phase='crossing'; c.t=0; }
    } else if (c.phase === 'crossing') { c.t += dt; c.y = CLEFT_TOP+5-c.t*30; if (c.t>0.5){c.phase='inside';c.t=0;} }
    else if (c.phase === 'inside') { wander(c, dt, 25); if (!insidePre(c.x,c.y,30)){c.x+=(PRE_CX-c.x)*0.05;c.y+=(PRE_CY-c.y)*0.05;} }
  }
  cortisolMols = cortisolMols.filter(c => c.life < 8);

  const rotating = Math.abs(revolverDisplay - revolverRotation) > 0.01;
  if (rotating) {
    const fullCount = vesicles.filter(v => v.count >= VES_FULL).length;
    const rotSpeed = 2.5 + fullCount*1.5;
    const diff = revolverRotation - revolverDisplay;
    revolverDisplay += Math.sign(diff)*Math.min(Math.abs(diff), rotSpeed*dt);
  }
  const nBusy = vesicles.filter(v => v.state !== 'barrel').length;
  const maxDock = (P_.PAR_DOCK && therapy.block > 0.5) ? 2 : 1;
  for (let i = 0; i < 3; i++) {
    const v = vesicles[i];
    v.t += dt;
    if (v.state === 'barrel') { const p = barrelPos(i, revolverDisplay); v.x = p.x; v.y = p.y; }
    else if (v.state === 'docking') { v.x += (DOCK.x-v.x)*0.12; v.y += (DOCK.y-v.y)*0.12; if (dist(v.x,v.y,DOCK.x,DOCK.y)<4){v.state='fusing';v.t=0;} }
    else if (v.state === 'fusing') {
      v.x = DOCK.x; v.y = DOCK.y;
      if (v.t > 0.45) {
        if (Math.random() < relProb) {
          let cyto = molecules.filter(mm => mm.state === 'cytosol' && !mm.degradation);
          while (v.count < VES_FULL && cyto.length > 0) { const pick = cyto.pop(); pick.dead = true; v.count++; }
          const spawnY = CLEFT_TOP + 10;
          for (let k = 0; k < v.count; k++) molecules.push(new Mol(DOCK.x+(Math.random()-0.5)*40, spawnY+Math.random()*12, 'cleft'));
          ringAdd('rel'); ring.rel[ringIdx] += v.count - 1;
          relPending += v.count;
          v.count = 0;
        }
        v.state = 'returning'; v.t = 0;
      }
    }
    else if (v.state === 'returning') { const home = barrelPos(i, revolverDisplay); v.x += (home.x-v.x)*0.12; v.y += (home.y-v.y)*0.12; if (dist(v.x,v.y,home.x,home.y)<4){v.state='barrel';v.t=0;} }
  }
  if (nBusy < maxDock && !rotating) {
    const activeIdx = activeVesicleIdx(revolverRotation);
    const activeVes = vesicles[activeIdx];
    const delay = (autoRec.active ? 2.0 : 0.3) / activeVes.speed;
    if (activeVes.count >= VES_FULL && activeVes.t > delay && caActive) { activeVes.state='docking'; activeVes.t=0; }
    else if (activeVes.count < VES_FULL) {
      const otherFull = vesicles.some((vv,idx)=>idx!==activeIdx && vv.count>=VES_FULL);
      if (otherFull) revolverRotation = (revolverRotation+1)%3;
    }
  }

  for (let m of molecules) {
    m.t += dt; m.life += dt;
    if (m.state === 'trp') { moveToward(m, SYN.branch[0], SYN.branch[1], 80, dt);
      if (dist(m.x,m.y,SYN.branch[0],SYN.branch[1])<8){ m.state = Math.random()<tdoEff?'kyn':'htp'; m.t=0; } }
    else if (m.state === 'htp') { moveToward(m, SYN.ht[0], SYN.ht[1], 80, dt);
      if (dist(m.x,m.y,SYN.ht[0],SYN.ht[1])<8){ m.state='cytosol'; m.t=0; m.life=0; } }
    else if (m.state === 'kyn') { moveToward(m, KYN_PATH.kyn[0], KYN_PATH.kyn[1], 70, dt);
      if (dist(m.x,m.y,KYN_PATH.kyn[0],KYN_PATH.kyn[1])<8){ m.state = Math.random()<prm.kynShift?'quin':'kyna'; m.t=0; } }
    else if (m.state === 'kyna') { moveToward(m, KYN_PATH.kyna[0], KYN_PATH.kyna[1], 60, dt);
      if (dist(m.x,m.y,KYN_PATH.kyna[0],KYN_PATH.kyna[1])<8){ m.dead=true; } }
    else if (m.state === 'quin') { moveToward(m, KYN_PATH.quin[0], KYN_PATH.quin[1], 60, dt);
      if (dist(m.x,m.y,KYN_PATH.quin[0],KYN_PATH.quin[1])<8){ m.dead=true; } }
    else if (m.state === 'cytosol') {
      wander(m, dt, 130);
      const wall = 8;
      const ex = (m.x-PRE_CX)/PRE_RX, ey = (m.y-PRE_CY)/PRE_RY;
      const eNorm = Math.sqrt(ex*ex+ey*ey);
      if (eNorm > 1 - wall/Math.min(PRE_RX,PRE_RY)) {
        let nx = ex/PRE_RX, ny = ey/PRE_RY; const nl = Math.hypot(nx,ny)||1; nx/=nl; ny/=nl;
        if (eNorm > 1) { m.x = PRE_CX + ex/eNorm*PRE_RX; m.y = PRE_CY + ey/eNorm*PRE_RY; }
        const vn = m.vx*nx + m.vy*ny;
        if (vn > 0) { m.vx -= 2*vn*nx; m.vy -= 2*vn*ny; m.vx*=0.85; m.vy*=0.85; }
        const push = 40*Math.max(0, eNorm-(1-wall/Math.min(PRE_RX,PRE_RY)))/(wall/Math.min(PRE_RX,PRE_RY));
        m.vx -= nx*push*dt*10; m.vy -= ny*push*dt*10;
      }
      if (m.degradation && !insidePre(m.x,m.y,10)) { m.x = PRE_CX+ex/Math.max(eNorm,1e-6)*PRE_RX*0.9; m.y = PRE_CY+ey/Math.max(eNorm,1e-6)*PRE_RY*0.9; }
      if (!m.degradation && m.life > 2.0 && cytosolCount > P_.CYTO_SET) {
        const overload = cytosolCount - P_.CYTO_SET;
        const rate = (P_.MAO_BASE+0.05*overload)*prm.mao;
        if (Math.random() < rate*dt) m.degradation = true;
      }
      if (m.life > P_.CYTO_LIFE) { m.dead = true; ringAdd('lost'); continue; }
      if (m.degradation) {
        const dx = MAO_POS.x-m.x, dy = MAO_POS.y-m.y, d = dist(m.x,m.y,MAO_POS.x,MAO_POS.y);
        m.vx += dx/d*500*dt; m.vy += dy/d*500*dt;
        if (d < 45) { m.dead = true; ringAdd('mao'); }
        continue;
      }
      let best=null, bc=VES_FULL+1, bi=-1;
      for (let i=0;i<3;i++){ if (vesicles[i].count<VES_FULL && vesicles[i].state==='barrel'){ if (vesicles[i].count<bc){bc=vesicles[i].count;best=vesicles[i];bi=i;} } }
      if (best) {
        const pos = barrelPos(bi, revolverDisplay);
        const dx=pos.x-m.x, dy=pos.y-m.y, d=dist(m.x,m.y,pos.x,pos.y)||1;
        if (insidePre(pos.x,pos.y,30)) {
          m.vx += dx/d*180*dt; m.vy += dy/d*180*dt;
          if (d < 26) { best.count++; m.dead = true; continue; }
        }
      }
    }
    else if (m.state === 'cleft') {
      const step = P_.D * Math.sqrt(dt);
      m.x += gaussRandom()*step; m.y += gaussRandom()*step;
      if (m.y < CLEFT_TOP+4) m.y = CLEFT_TOP+4;
      if (m.y > CLEFT_BOT-4) m.y = CLEFT_BOT-4;
      if (m.x < CLEFT_LEFT || m.x > CLEFT_RIGHT) {
        if (Math.random() < P_.EDGE) { m.x = m.x<CLEFT_LEFT ? CLEFT_LEFT+4 : CLEFT_RIGHT-4; }
        else { m.dead = true; ringAdd('diff'); continue; }
      }
      let captured = false;
      for (let s of serts) {
        if (!s.blocked && dist(m.x,m.y,s.x,s.y) < s.radius) {
          ringAdd('up'); m.recycled = true;
          m.degradation = Math.random() < prm.maoFraction;
          m.state='cytosol'; m.t=0; m.life=0; m.x=s.x; m.y=s.y-40; m.vx=(Math.random()-0.5)*40; m.vy=-60;
          captured = true; break;
        }
      }
      if (captured) continue;
      for (let r of postRecs) {
        const sens = r.type==='5-HT1A'?prm.ht1a:(r.type==='5-HT2A'?prm.ht2a:prm.ht3);
        if (dist(m.x,m.y,r.x,r.y) < r.radius) {
          const prob = 1 - Math.exp(-r.aff*sens*P_.BIND*dt);
          if (Math.random() < prob) { r.occ=m; r.bindings++; m.state='bound'; m.target=r; m.t=0; break; }
        }
      }
      if (m.state === 'cleft' && !autoRec.occupied && dist(m.x,m.y,autoRec.x,autoRec.y+15) < P_.AUTO_R) {
        if (Math.random() < Math.min(1.0, prm.autoSens * P_.AUTO_K)) { autoRec.occupied=m; m.state='auto'; m.t=0; }
      }
    }
    else if (m.state === 'bound') {
      if (m.t > 1.0/m.target.aff) { m.target.occ=null; m.state='cleft'; m.t=0; }
    }
    else if (m.state === 'auto') {
      if (m.t > 2.0) { autoRec.occupied=null; autoRec.active=true; autoRec.activeT=3.0*prm.autoSens; m.dead=true; }
    }
  }
  {
    const frees = molecules.filter(m => m.state === 'cleft');
    if (frees.length > P_.MAXPOP) {
      frees.sort((a,b)=>a.life-b.life);
      for (let i=0;i<frees.length-P_.MAXPOP;i++){ frees[i].dead=true; ringAdd('diff'); }
    }
  }
  molecules = molecules.filter(m => !m.dead);
  if (autoRec.active) { autoRec.activeT -= dt; if (autoRec.activeT<=0) autoRec.active=false; }
  for (let e of calciumEvents) e.t += dt;
  calciumEvents = calciumEvents.filter(e => e.t < 0.8);
}

// ---- harness ----
const SECONDS_PER_PHASE = parseFloat(process.env.SPP || '120');
const SEEDS = parseInt(process.env.SEEDS || '3');
const DT = 1/30;

function mulberry32(a){ return function(){ a|=0; a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; } }

function runSeed(seed, stress){
  Math.random = mulberry32(seed);
  molecules=[]; calciumEvents=[]; cortisolMols=[]; time=0; simTime=0;
  stressLevel=0; targetStress=stress; tdoRamp=0; targetPhase=0;
  therapy={block:0,desens:0,adapt:0};
  revolverRotation=0; revolverDisplay=0;
  vesicles=[{count:2,pool:'RRP',speed:1.0,t:0,state:'barrel',x:0,y:0},{count:2,pool:'Rec',speed:0.8,t:0,state:'barrel',x:0,y:0},{count:2,pool:'Res',speed:0.5,t:0,state:'barrel',x:0,y:0}];
  for (let s of serts){ s.blocked=false; }
  autoRec.occupied=null; autoRec.active=false; autoRec.activeT=0;
  for (let r of postRecs){ r.occ=null; r.bindings=0; }
  occInteg=0; occEma=null; relPending=0; relProb=P_.REL_BASE[0];
  synTimer=0; caTimer=0; cortTimer=0;
  ring={ rel:Array(RING_N).fill(0), up:Array(RING_N).fill(0), diff:Array(RING_N).fill(0), mao:Array(RING_N).fill(0), lost:Array(RING_N).fill(0) };
  ringIdx=0; secTimer=0;
  resetRadii();

  const results = {};
  let acc={n:0,sum:0}, timer=0;
  for (let ph=0; ph<4; ph++){
    targetPhase = ph;
    acc={n:0,sum:0}; timer=0;
    results[ph]=[];
    const steps = Math.round(SECONDS_PER_PHASE/DT);
    for (let i=0;i<steps;i++){
      update(DT);
      acc.n++; acc.sum+=cleftCount(); timer+=DT;
      if (timer>=10){ results[ph].push(acc.sum/acc.n); acc={n:0,sum:0}; timer=0; }
    }
  }
  return results;
}

const STRESS = parseFloat(process.env.STRESS || '0');
const all = [];
for (let s=0;s<SEEDS;s++) all.push(runSeed(1000+s*77, STRESS));
console.error(`Params(stress=${STRESS}): ${JSON.stringify(P_)}`);
let totalPass = 0, totalNeed = 0;
for (let ph=0; ph<4; ph++){
  const range = (STRESS<0.5?P_.TARGETS_CALM:P_.TARGETS_STR)[ph];
  const mid = ph>=1 ? 3 : 2;
  for (let s=0;s<SEEDS;s++){
    const rs = all[s][ph];
    console.log(`F${ph} seed${s}: [` + rs.map(v=>v.toFixed(1)).join(', ') + `] цель ${range[0]}-${range[1]}`);
  }
  let passSeeds=0;
  for (let s=0;s<SEEDS;s++){
    const tail = all[s][ph].slice(-mid);
    const okCount = tail.filter(v=>v>=range[0]&&v<=range[1]).length;
    if (okCount>=Math.ceil(mid*0.7)) passSeeds++;
  }
  totalPass += passSeeds; totalNeed += SEEDS;
  console.log(`== F${ph}: сидов в норме ${passSeeds}/${SEEDS}`);
}
console.log(`TOTAL: ${totalPass}/${totalNeed}`);
