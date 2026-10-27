/*
 * hc-core.js: a small 2D horizontal-convection solver for the browser.
 *
 * Vorticity / streamfunction form of the Boussinesq equations in an x-z slice:
 *   D(omega)/Dt = db/dx + nu * lap(omega),   lap(psi) = omega,   u = -psi_z, w = psi_x
 *   D(b)/Dt     = kappa * lap(b),            b relaxed toward b_s(x) in the top cells
 * Free-slip walls and a rigid lid. Gaussian ridges are masked out as solid cells.
 * Semi-Lagrangian advection, conjugate-gradient Poisson solve (warm-started).
 *
 * This is a teaching toy: laminar, 2D, low Rayleigh number. It is NOT the
 * production Oceananigans.jl model.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.HCSim = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const L = 4;            // half-width of the domain, x in [-L, L]
  const NX = 384;
  const NY = 48;          // depth H = 1, so cells are square: h = 2L/NX = 1/NY
  const h = (2 * L) / NX;
  const W = NX + 2;       // one solid pad cell on every side
  const N = W * (NY + 2);

  const RIDGES = [
    { x: 0, a: 0.26 },
    { x: 2, a: 0.4 },
  ];
  const SIGMA = 0.24;

  function create() {
    const f64 = () => new Float64Array(N);
    const b = f64(), b2 = f64(), om = f64(), om2 = f64(), psi = f64();
    const u = f64(), w = f64();
    const rr = f64(), pp = f64(), Ap = f64();
    const mask = new Uint8Array(N);
    const clear = new Uint8Array(N);
    let fl = new Int32Array(0);
    let t = 0;

    const params = {
      B: 0.06,       // surface buoyancy amplitude
      ridges: 2,     // 0, 1 or 2
      nu: 1e-4,
      kappa: 5e-5,
      tau0: 0.1,     // surface relaxation time, top cell
      tau1: 0.4,     // second cell
      noise: 0,      // amplitude of random surface-forcing wobble (fraction of B)
    };
    const api = { params, NX, NY, W, N, L, h, b, mask, RIDGES, SIGMA, lastIters: 0 };

    function topo(x) {
      let m = 0;
      for (let r = 0; r < params.ridges; r++) {
        const R = RIDGES[r];
        const v = R.a * Math.exp(-((x - R.x) * (x - R.x)) / (2 * SIGMA * SIGMA));
        if (v > m) m = v;
      }
      return m;
    }
    api.topo = topo;

    const bsAt = (x) => -params.B * Math.cos((Math.PI * (x + L)) / (2 * L));
    api.bs = bsAt;

    function reset() {
      for (const a of [b, b2, om, om2, psi, u, w, rr, pp, Ap]) a.fill(0);
      mask.fill(0);
      const list = [];
      for (let k = 0; k < NY; k++) {
        const z = -1 + (k + 0.5) * h;
        for (let i = 0; i < NX; i++) {
          const x = -L + (i + 0.5) * h;
          if (z > -1 + topo(x)) {
            const c = (k + 1) * W + (i + 1);
            mask[c] = 1;
            list.push(c);
          }
        }
      }
      fl = Int32Array.from(list);
      // cells whose whole 4x4 neighbourhood (offsets -1..+2) is fluid can use cubic interpolation
      clear.fill(0);
      for (let k = 1; k < NY - 2; k++) {
        for (let i = 1; i < NX - 2; i++) {
          let ok = 1;
          for (let dk = -1; dk <= 2 && ok; dk++)
            for (let di = -1; di <= 2; di++)
              if (!mask[(k + dk + 1) * W + (i + di + 1)]) { ok = 0; break; }
          clear[(k + 1) * W + (i + 1)] = ok;
        }
      }
      t = 0;
    }

    // semi-Lagrangian advection of src into dst along (u, w)
    function advect(src, dst, dt, weighted) {
      const s = dt / h;
      for (let n = 0; n < fl.length; n++) {
        const c = fl[n];
        const i = (c % W) - 1;
        const k = ((c / W) | 0) - 1;
        let x = i - u[c] * s;
        let z = k - w[c] * s;
        if (x < 0) x = 0; else if (x > NX - 1) x = NX - 1;
        if (z < 0) z = 0; else if (z > NY - 1) z = NY - 1;
        const i0 = x | 0, k0 = z | 0;
        const fx = x - i0, fz = z - k0;
        const c00 = (k0 + 1) * W + i0 + 1, c10 = c00 + 1, c01 = c00 + W, c11 = c01 + 1;
        if (clear[c00]) {
          // Catmull-Rom cubic: much less numerical diffusion than bilinear
          const x2 = fx * fx, x3 = x2 * fx, z2 = fz * fz, z3 = z2 * fz;
          const ax0 = (-x3 + 2 * x2 - fx) * 0.5, ax1 = (3 * x3 - 5 * x2 + 2) * 0.5;
          const ax2 = (-3 * x3 + 4 * x2 + fx) * 0.5, ax3 = (x3 - x2) * 0.5;
          const az0 = (-z3 + 2 * z2 - fz) * 0.5, az1 = (3 * z3 - 5 * z2 + 2) * 0.5;
          const az2 = (-3 * z3 + 4 * z2 + fz) * 0.5, az3 = (z3 - z2) * 0.5;
          let v = 0;
          let r = c00 - W - 1;
          v += az0 * (ax0 * src[r] + ax1 * src[r + 1] + ax2 * src[r + 2] + ax3 * src[r + 3]);
          r += W;
          v += az1 * (ax0 * src[r] + ax1 * src[r + 1] + ax2 * src[r + 2] + ax3 * src[r + 3]);
          r += W;
          v += az2 * (ax0 * src[r] + ax1 * src[r + 1] + ax2 * src[r + 2] + ax3 * src[r + 3]);
          r += W;
          v += az3 * (ax0 * src[r] + ax1 * src[r + 1] + ax2 * src[r + 2] + ax3 * src[r + 3]);
          if (weighted) {
            // keep buoyancy monotone: no new extrema from the cubic
            const s00 = src[c00], s10 = src[c10], s01 = src[c01], s11 = src[c11];
            const lo = Math.min(s00, s10, s01, s11), hi = Math.max(s00, s10, s01, s11);
            v = v < lo ? lo : v > hi ? hi : v;
          }
          dst[c] = v;
          continue;
        }
        let w00 = (1 - fx) * (1 - fz), w10 = fx * (1 - fz), w01 = (1 - fx) * fz, w11 = fx * fz;
        if (weighted) {
          w00 *= mask[c00]; w10 *= mask[c10]; w01 *= mask[c01]; w11 *= mask[c11];
          const sm = w00 + w10 + w01 + w11;
          if (sm < 1e-9) { dst[c] = src[c]; continue; }
          dst[c] = (w00 * src[c00] + w10 * src[c10] + w01 * src[c01] + w11 * src[c11]) / sm;
        } else {
          dst[c] = w00 * src[c00] + w10 * src[c10] + w01 * src[c01] + w11 * src[c11];
        }
      }
    }

    function solvePsi(tol, maxIt) {
      const h2 = h * h;
      let rs = 0, bn = 0;
      for (let n = 0; n < fl.length; n++) {
        const c = fl[n];
        const ap = 4 * psi[c] - psi[c - 1] - psi[c + 1] - psi[c - W] - psi[c + W];
        const r = -h2 * om[c] - ap;
        rr[c] = r; pp[c] = r; rs += r * r; bn += om[c] * om[c];
      }
      bn = Math.sqrt(bn) * h2;
      if (bn < 1e-14) { api.lastIters = 0; return; }
      const stop = (tol * bn) * (tol * bn);
      let it = 0;
      for (; it < maxIt && rs > stop; it++) {
        let pAp = 0;
        for (let n = 0; n < fl.length; n++) {
          const c = fl[n];
          const a = 4 * pp[c] - pp[c - 1] - pp[c + 1] - pp[c - W] - pp[c + W];
          Ap[c] = a; pAp += pp[c] * a;
        }
        const alpha = rs / pAp;
        let rsn = 0;
        for (let n = 0; n < fl.length; n++) {
          const c = fl[n];
          psi[c] += alpha * pp[c];
          const r = rr[c] - alpha * Ap[c];
          rr[c] = r; rsn += r * r;
        }
        const beta = rsn / rs; rs = rsn;
        for (let n = 0; n < fl.length; n++) {
          const c = fl[n];
          pp[c] = rr[c] + beta * pp[c];
        }
      }
      api.lastIters = it;
    }

    function step(dt) {
      const { nu, kappa } = params;
      const inv2h = 1 / (2 * h), invh2 = 1 / (h * h);

      // 1. vorticity: advect, baroclinic torque, viscosity
      advect(om, om2, dt, false);
      for (let n = 0; n < fl.length; n++) {
        const c = fl[n];
        const bE = mask[c + 1] ? b[c + 1] : b[c];
        const bW = mask[c - 1] ? b[c - 1] : b[c];
        om2[c] += dt * (bE - bW) * inv2h;
      }
      for (let n = 0; n < fl.length; n++) {
        const c = fl[n];
        const lap = (om2[c - 1] + om2[c + 1] + om2[c - W] + om2[c + W] - 4 * om2[c]) * invh2;
        om[c] = om2[c] + dt * nu * lap;
      }

      // 2. streamfunction and velocities
      solvePsi(4e-3, 40);
      for (let n = 0; n < fl.length; n++) {
        const c = fl[n];
        u[c] = -(psi[c + W] - psi[c - W]) * inv2h;
        w[c] = (psi[c + 1] - psi[c - 1]) * inv2h;
      }

      // 3. buoyancy: advect, diffuse (no flux through solids), relax at the surface
      advect(b, b2, dt, true);
      for (let n = 0; n < fl.length; n++) {
        const c = fl[n];
        const v = b2[c];
        const E = mask[c + 1] ? b2[c + 1] : v;
        const Wd = mask[c - 1] ? b2[c - 1] : v;
        const Nn = mask[c + W] ? b2[c + W] : v;
        const S = mask[c - W] ? b2[c - W] : v;
        b[c] = v + dt * kappa * (E + Wd + Nn + S - 4 * v) * invh2;
      }
      const wob = params.noise ? 1 + params.noise * (Math.random() - 0.5) * 2 : 1;
      const a0 = 1 - Math.exp(-dt / params.tau0);
      const a1 = 1 - Math.exp(-dt / params.tau1);
      for (let i = 0; i < NX; i++) {
        const target = bsAt(-L + (i + 0.5) * h) * wob;
        const c0 = NY * W + i + 1;       // top row (k = NY-1)
        const c1 = c0 - W;
        if (mask[c0]) b[c0] += (target - b[c0]) * a0;
        if (mask[c1]) b[c1] += (target - b[c1]) * a1;
      }
      t += dt;
    }

    api.step = step;
    api.reset = reset;
    Object.defineProperty(api, 'time', { get: () => t });
    Object.defineProperty(api, 'fluid', { get: () => fl });
    reset();
    return api;
  }

  return { create, L, NX, NY, h };
});
