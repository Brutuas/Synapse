#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Регрессионное тестирование synapse_v9.5.html в реальном браузере (Playwright).

Методика:
 - ускоряем часы страницы (performance.now), 1 фаза = 60 с симуляции;
 - после смены фазы 12 с переходного процесса не учитываем;
 - каждые 4 с читаем ДАШБОРД через DOM-мутацию (MutationObserver) — без canvas OCR;
 - метрики: щель(мол/нМ), выброс Гц, захват SERT Гц, MAO окно, распр 1A/2A/3,
   события 10с, соматика, Tph2, барабан, ПИ(EMA/gain/I).
Использование: python3 tuner/regtest.py [file.html] [seconds_per_phase] [warmup] [sample_every]
"""
import sys, re, json, statistics
from playwright.sync_api import sync_playwright

FILE = sys.argv[1] if len(sys.argv) > 1 else "/workspace/synapse_v9.5.html"
PHASE_SEC = float(sys.argv[2]) if len(sys.argv) > 2 else 60
WARMUP = float(sys.argv[3]) if len(sys.argv) > 3 else 12
SAMPLE_EVERY = float(sys.argv[4]) if len(sys.argv) > 4 else 4
SPEED = 30   # множитель виртуальных часов (dt капается в 0.05 — физика не меняется)

OBSERVER_JS = """
window.__snap = null;
const ta = document.getElementById('testReport');
ta.style.display = 'block';   // renderTestReport пишет только в видимый textarea
new MutationObserver(() => {
    const m = ta.value.match(/ДАШБОРД[\\s\\S]*?\\n([\\s\\S]*?)\\nСТАТУС/);
    window.__snap = m ? m[1] : null;
}).observe(ta, {childList:true, attributes:true});
"""

# ускорение стенда: патченные часы + отключённая отрисовка (draw=noop).
# update() идёт своими реальными шагами dt<=0.05 — динамика идентична браузеру.
# Тестер не парсит canvas: window.__dash() возвращает те же строки дашборда,
# что рисует drawDashboard (собираются через подменённый fillText).
INIT_JS = """
window.__t0 = Date.now();
const S = %d;
performance.now = () => (Date.now() - window.__t0) * S + 1000;
window.requestAnimationFrame = cb => setTimeout(() => cb(performance.now()), 2);
window.__dashLines = [];
window.addEventListener('load', () => {
    if (typeof draw === 'function') draw = function(){};
    const realDash = drawDashboard;
    window.__dashOn = false;
    drawDashboard = function(){
        if (!window.__dashOn) return;
        const proto = Object.getPrototypeOf(ctx);
        const realFill = proto.fillText;
        window.__dashLines = [];
        proto.fillText = function(txt){ if (String(txt).indexOf('ДАШБОРД')>=0 || String(txt).indexOf('СТАТУС')>=0) return; window.__dashLines.push(String(txt)); };
        try { realDash.call(this); } finally { proto.fillText = realFill; }
    };
    window.__dash = function(){ window.__dashOn = true; try { drawDashboard(); } finally { window.__dashOn = false; } return window.__dashLines.join('\\n'); };
});
""" % SPEED

def parse_snap(txt):
    d = {}
    if not txt: return d
    for line in txt.split("\n"):
        line = line.strip()
        m = re.match(r"5-HT в щели: (\d+) \(~(\d+) нМ\)", line)
        if m: d["cleft"] = int(m.group(1)); d["nm"] = int(m.group(2))
        # аддитивная телеметрия (v9.5-T): состояния, барабан, ПИ-ручка, задержка докинга
        m = re.match(r"T \[(.*)\]", line)
        if m:
            d["tele"] = {}
            for kv in m.group(1).split("|"):
                if "=" in kv:
                    k, v = kv.split("=", 1); d["tele"][k] = v
        m = re.match(r"норма фазы: (\d+)-(\d+) нМ", line)
        if m: d["lo_nm"], d["hi_nm"] = int(m.group(1)), int(m.group(2))
        m = re.match(r"Выброс: \d+ \(([\d.]+) Гц\)", line)
        if m: d["relHz"] = float(m.group(1))
        m = re.match(r"Захват SERT: \d+ \(([\d.]+) Гц\)", line)
        if m: d["sertHz"] = float(m.group(1))
        m = re.match(r"MAO окно: (\d+)", line)
        if m: d["mao10"] = int(m.group(1))
        m = re.match(r"Распр: (\d+)/(\d+)/(\d+)", line)
        if m: d["dist"] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        m = re.match(r"События 10с 1A/2A/3: (\d+)/(\d+)/(\d+)", line)
        if m: d["ev10"] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        m = re.match(r"Тревога тела \(соматика\): ([\d.]+)/10", line)
        if m: d["soma"] = float(m.group(1))
        m = re.match(r"Истощение Tph2: (\d+)%", line)
        if m: d["exh"] = int(m.group(1))
        m = re.match(r"Барабан: (\d+),(\d+),(\d+)", line)
        if m: d["bar"] = (m.group(1), m.group(2), m.group(3))
        m = re.match(r"Уровень тревоги: (\d+)%", line)
        if m: d["anx"] = int(m.group(1))
        m = re.match(r"ПИ-диаг: EMA=([\d.]+|—) gain=([\d.]+) I=([-\d.]+)", line)
        if m: d["pi"] = (m.group(1), m.group(2), m.group(3))
        m = re.match(r"Время: ([\d.]+) c", line)
        if m: d["t"] = float(m.group(1))
    return d

def main():
    results = {}
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1400, "height": 1100})
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.add_init_script(INIT_JS)   # часы ускорены ДО старта loop()
        page.goto("file://" + FILE)
        page.evaluate("() => { document.getElementById('testReport').style.display='block'; }")
        for phase in range(4):
            page.click(f'#therapyBar button[data-ph="{phase}"]')
            samples = []
            t_end = WARMUP + PHASE_SEC
            last_t = -1
            stall = 0
            while True:
                page.wait_for_timeout(int(SAMPLE_EVERY / SPEED * 1000 * 6))
                # прямой pull-замер дашборда из canvas-текста (надёжнее DOM-наблюдателя)
                snap = page.evaluate("() => window.__dash ? window.__dash() : null")
                d = parse_snap(snap)
                t = d.get("t", 0)
                if t <= last_t:
                    stall += 1
                    if stall > 15:
                        print(f"[F{phase}] ВРЕМЯ ЗАСТОРИЛО НА t={last_t} — стоп")
                        break
                else:
                    stall = 0; last_t = t
                if t >= WARMUP and "cleft" in d:
                    samples.append(d)
                if t >= t_end:
                    break
            # журнал testLog тоже полезнеем
            log_txt = page.evaluate("() => document.getElementById('testReport').value")
            ok_mimo = re.findall(r"\[(OK|МИМО)\]", log_txt or "")
            verdict = re.search(r"ВЕРДИКТ: (\d)/5", log_txt or "")
            results[f"F{phase}"] = {
                "n_samples": len(samples),
                "cleft_mean": round(statistics.mean([s["cleft"] for s in samples]), 1) if samples else None,
                "cleft_min": min([s["cleft"] for s in samples], default=None),
                "cleft_max": max([s["cleft"] for s in samples], default=None),
                "relHz_mean": round(statistics.mean([s.get("relHz", 0) for s in samples]), 2) if samples else None,
                "sertHz_mean": round(statistics.mean([s.get("sertHz", 0) for s in samples]), 2) if samples else None,
                "mao10_mean": round(statistics.mean([s.get("mao10", 0) for s in samples]), 1) if samples else None,
                "dist_last": samples[-1].get("dist") if samples else None,
                "soma_last": samples[-1].get("soma") if samples else None,
                "exh_last": samples[-1].get("exh") if samples else None,
                "pi_last": samples[-1].get("pi") if samples else None,
                "ok_windows": f"{sum(1 for x in ok_mimo if x=='OK')}/{len(ok_mimo)}",
                "verdict_tail": verdict.group(1) if verdict else None,
            }
            print(f"F{phase}: {json.dumps(results[f'F{phase}'], ensure_ascii=False)}")
        if errors:
            print("PAGE ERRORS:", errors[:5])
        else:
            print("JS errors: нет")
        browser.close()
    out = "/workspace/tuner/regtest_result.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print("saved ->", out)

if __name__ == "__main__":
    main()
