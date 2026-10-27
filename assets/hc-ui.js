/* hc-ui.js: canvas rendering and controls for the live horizontal-convection demo */
(function () {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const cv = $('sim-canvas');
  if (!cv || !window.HCSim) return;

  const sim = HCSim.create();
  const { NX, NY, W, L } = sim;
  const DT = 0.1;
  const ctx = cv.getContext('2d');
  const off = document.createElement('canvas');
  off.width = NX; off.height = NY;
  const octx = off.getContext('2d');
  const img = octx.createImageData(NX, NY);
  const reduced = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

  // colour map: dense (navy) -> neutral -> light (dark red)
  const STOPS = [
    [-1, [20, 28, 75]], [-0.6, [28, 80, 160]], [-0.25, [110, 170, 215]], [0, [238, 232, 228]],
    [0.25, [225, 150, 120]], [0.6, [175, 50, 40]], [1, [60, 10, 30]],
  ];
  const LUT = new Uint8ClampedArray(256 * 3);
  for (let n = 0; n < 256; n++) {
    const v = (n / 255) * 2 - 1;
    let col = STOPS[STOPS.length - 1][1];
    for (let i = 0; i < STOPS.length - 1; i++) {
      const [a, ca] = STOPS[i], [b, cb] = STOPS[i + 1];
      if (v <= b) { const f = (v - a) / (b - a); col = ca.map((c, j) => c + (cb[j] - c) * f); break; }
    }
    LUT.set(col, n * 3);
  }
  const rgb = (v) => { const n = Math.round((Math.max(-1, Math.min(1, v)) + 1) * 127.5) * 3; return `rgb(${LUT[n]},${LUT[n + 1]},${LUT[n + 2]})`; };

  // sizing: the real domain is 8:1, so stretch vertically to keep it readable
  let stretch = 2, cssW = 0, cssH = 0;
  function resize() {
    const r = cv.parentElement.getBoundingClientRect();
    cssW = r.width;
    stretch = cssW >= 900 ? 2 : cssW >= 560 ? 3 : 4;
    cssH = (cssW / 8) * stretch;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    cv.style.height = cssH + 'px';
    cv.width = Math.round(cssW * dpr);
    cv.height = Math.round(cssH * dpr);
    const s = $('stretch'); if (s) s.textContent = stretch;
    draw();
  }

  const bfill = new Float64Array(NX * NY);
  function draw() {
    const d = img.data, B = sim.params.B, b = sim.b, mask = sim.mask;
    for (let k = NY - 1; k >= 0; k--) {
      for (let i = 0; i < NX; i++) {
        const c = (k + 1) * W + (i + 1), q = (NY - 1 - k) * NX + i;
        const v = mask[c] ? b[c] : (k < NY - 1 ? bfill[q - NX] : 0);
        bfill[q] = v;
        const n = Math.round((Math.max(-1, Math.min(1, v / B)) + 1) * 127.5) * 3;
        const o = q * 4;
        d[o] = LUT[n]; d[o + 1] = LUT[n + 1]; d[o + 2] = LUT[n + 2]; d[o + 3] = 255;
      }
    }
    octx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(off, 0, 0, cv.width, cv.height);
    if (sim.params.ridges > 0) {
      ctx.fillStyle = '#3d362f';
      ctx.beginPath();
      ctx.moveTo(0, cv.height);
      for (let px = 0; px <= cv.width; px += 2) {
        const x = -L + (px / cv.width) * 2 * L;
        ctx.lineTo(px, (1 - sim.topo(x)) * cv.height);
      }
      ctx.lineTo(cv.width, cv.height);
      ctx.closePath();
      ctx.fill();
    }
    const t = $('time'); if (t) t.textContent = sim.time.toFixed(0);
  }

  // surface-forcing strip and legend
  function paintStrip() {
    const s = $('forcing-strip'); if (!s) return;
    const stops = [];
    for (let i = 0; i <= 24; i++) {
      const x = -L + (i / 24) * 2 * L;
      stops.push(`${rgb(sim.bs(x) / sim.params.B)} ${(i / 24) * 100}%`);
    }
    s.style.background = `linear-gradient(90deg, ${stops.join(',')})`;
    const lg = $('legend-bar');
    if (lg) {
      const ls = [];
      for (let i = 0; i <= 12; i++) ls.push(`${rgb(i / 6 - 1)} ${(i / 12) * 100}%`);
      lg.style.background = `linear-gradient(90deg, ${ls.join(',')})`;
    }
    const lo = $('legend-lo'), hi = $('legend-hi');
    if (lo) lo.textContent = (-sim.params.B).toFixed(2);
    if (hi) hi.textContent = '+' + sim.params.B.toFixed(2);
  }

  // loop
  let running = !reduced, visible = true, raf = 0, steps = 1;
  function frame() {
    raf = 0;
    if (!running || !visible || document.hidden) return;
    for (let s = 0; s < steps; s++) sim.step(DT);
    draw();
    raf = requestAnimationFrame(frame);
  }
  function kick() { if (!raf && running && visible) raf = requestAnimationFrame(frame); }

  const playBtn = $('play');
  function setRunning(v) {
    running = v;
    playBtn.textContent = v ? 'Pause' : 'Play';
    playBtn.setAttribute('aria-pressed', String(!v));
    kick();
  }
  playBtn.addEventListener('click', () => setRunning(!running));
  $('restart').addEventListener('click', () => { sim.reset(); draw(); kick(); });

  document.querySelectorAll('input[name="ridges"]').forEach((el) =>
    el.addEventListener('change', () => { sim.params.ridges = +el.value; sim.reset(); draw(); kick(); }));

  const forcing = $('forcing'), forcingOut = $('forcing-out');
  forcing.addEventListener('input', () => {
    sim.params.B = +forcing.value;
    forcingOut.textContent = (+forcing.value).toFixed(2);
    sim.reset(); paintStrip(); draw(); kick();
  });
  const speed = $('speed');
  speed.addEventListener('input', () => { steps = +speed.value; });
  $('noise').addEventListener('change', (e) => { sim.params.noise = e.target.checked ? 0.6 : 0; });

  // hover readout
  const ro = $('readout');
  cv.addEventListener('pointermove', (e) => {
    const r = cv.getBoundingClientRect();
    const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
    const i = Math.min(NX - 1, Math.max(0, Math.floor(fx * NX)));
    const k = Math.min(NY - 1, Math.max(0, NY - 1 - Math.floor(fy * NY)));
    const c = (k + 1) * W + (i + 1);
    const x = -L + fx * 2 * L, z = -fy;
    ro.textContent = sim.mask[c]
      ? `x ${x.toFixed(2)}   z ${z.toFixed(2)}   buoyancy ${sim.b[c] >= 0 ? '+' : '−'}${Math.abs(sim.b[c]).toFixed(3)}`
      : 'sea floor';
  });
  cv.addEventListener('pointerleave', () => { ro.textContent = 'Hover over the water to read the buoyancy.'; });

  new IntersectionObserver((es) => { visible = es[0].isIntersecting; kick(); }, { threshold: 0.05 }).observe(cv);
  document.addEventListener('visibilitychange', kick);
  window.addEventListener('resize', resize);

  paintStrip();
  resize();
  setRunning(running);
  if (reduced) $('motion-note').hidden = false;
})();

/* page-level: depth gauge, looping movie */
(function () {
  const gauge = document.getElementById('depth-mark');
  if (gauge) {
    const upd = () => {
      const max = document.documentElement.scrollHeight - innerHeight;
      gauge.style.top = (max > 0 ? (scrollY / max) * 100 : 0) + '%';
    };
    addEventListener('scroll', upd, { passive: true });
    addEventListener('resize', upd);
    upd();
  }
  const v = document.getElementById('movie');
  if (v && !(window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches)) {
    new IntersectionObserver((es) => {
      if (es[0].isIntersecting) v.play().catch(() => {}); else v.pause();
    }, { threshold: 0.4 }).observe(v);
  }
})();
