import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { runCoupled } from '../src/viscous/context.js';
import { naca4 } from '../src/geometry/airfoil.js';
import { createHash } from 'node:crypto';

const fixtures = JSON.parse(await readFile(new URL('./fixtures/fortran/kernels.json', import.meta.url)));
const state = ([s, ampl, ctau, theta, deltaStar, wakeGap, ue], turbulent) => ({ s, theta, deltaStar, wakeGap, ue, aux: turbulent ? ctau : ampl });
const close = (a, b, tolerance = 2e-11) => assert.ok(Math.abs(a - b) < tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

test('all integral blocks and compressible chain derivatives agree with original Fortran', async () => {
  const fixture = JSON.parse(await readFile(new URL('./fixtures/fortran/integral.json', import.meta.url)));
  for (const [path, hash] of Object.entries(fixture.provenance.sha256)) {
    assert.equal(createHash('sha256').update(await readFile(new URL(`../${path}`, import.meta.url))).digest('hex'), hash, path);
  }
  for (const { input: c, expected } of fixture.cases) {
    const kernel = createIntegralKernel(c.parameters);
    const half = c.matched && { ...c.states[0], ...c.matched, theta: c.matched.theta / 2, deltaStar: c.matched.deltaStar / 2 };
    const [upper,lower]=c.surfaceStates??[half,half];
    const actual = c.regime === 'te' ? kernel.trailingEdge(upper, lower, c.states[1],c.gap??0)
      : kernel.interval({ upstream: c.states[0], downstream: c.states[1], regime: c.regime, tripS: c.tripS });
    if(c.regime==='te')for(const key of ['aux','theta','deltaStar'])close(actual.matched[key],c.matched[key]);
    for (const key of ['residual','upstream','downstream','reynoldsDerivative','machSquaredDerivative','tripDerivative']) {
      actual[key].flat().forEach((v,i) => close(v, expected[key].flat()[i], 3e-10));
    }
    if (c.regime === 'transition') {
      close(actual.transition.s, expected.transition.s);
      assert.equal(actual.transition.forced, expected.transition.forced);
    }
  }
});

test('topology-independent BL blocks preserve native laminar/turbulent/wake residuals and derivatives', () => {
  for (const { input, expected } of fixtures.intervals) {
    const kernel = createIntegralKernel({ reynolds: input.reynolds, velocityConvention: 'xfoil' });
    const [upstream, downstream] = input.stations.map(s => state(s, input.type >= 2));
    const result = kernel.interval({ upstream, downstream, regime: ['', 'laminar', 'turbulent', 'wake'][input.type] });
    result.residual.forEach((v, i) => close(v, -expected.residual[i]));
    for (const key of ['upstream', 'downstream']) result[key].forEach((row, i) => row.forEach((v, j) => close(v, expected[key][i][j])));
    // The legacy fixture differentiates fixed REYBL and physical U inside
    // BLDIF. This adapter differentiates freestream Re and the selected edge
    // velocity convention; the new native integral fixture checks that chain.
  }
});

test('optional complete BL Jacobian preserves native residuals and matches independent state differences',async()=>{
  const fixture=JSON.parse(await readFile(new URL('./fixtures/fortran/integral.json',import.meta.url)));
  for(const {input:c} of fixture.cases){
    if(['te','transition'].includes(c.regime))continue;
    const kernel=createIntegralKernel({...c.parameters,exactJacobian:true}),native=createIntegralKernel(c.parameters);
    const input={upstream:c.states[0],downstream:c.states[1],regime:c.regime,tripS:c.tripS},block=kernel.interval(input);
    assert.deepEqual(block.residual,native.interval(input).residual);
    for(const side of ['upstream','downstream'])for(const [k,key] of ['aux','theta','deltaStar','ue','s'].entries()){
      const value=input[side][key],h=Math.cbrt(Number.EPSILON)*Math.max(Math.abs(value),k===0?.01:k<3?1e-7:1e-6);
      const p=kernel.interval({...input,[side]:{...input[side],[key]:value+h}}).residual;
      const m=kernel.interval({...input,[side]:{...input[side],[key]:value-h}}).residual;
      for(let row=0;row<3;row++)close(block[side][row][k],(p[row]-m[row])/(2*h),2e-6);
    }
  }
});

test('physical edge velocities obey independent perfect-gas and viscosity relations', () => {
  const mach = .4; const reynolds = 2e6; const gamma = 1.4;
  const kernel = createIntegralKernel({ mach, reynolds });
  const input = { s: .2, ue: 1.3, theta: .0006, deltaStar: .0016, aux: 0 };
  const p = kernel.station(input);
  const temperature = 1 + .5 * (gamma - 1) * mach ** 2 * (1 - input.ue ** 2);
  const rho = temperature ** (1 / (gamma - 1));
  close(p.ue, input.ue); close(p.rho, rho);
  close(p.machSquared, mach ** 2 * input.ue ** 2 / temperature);
  // XFOIL uses a stagnation-temperature-based Sutherland constant HVRAT=.35.
  const stagnationTemperature = 1 + .5 * (gamma - 1) * mach ** 2;
  const sutherland = .35 * stagnationTemperature;
  const viscosity = temperature ** 1.5 * (1 + sutherland) / (temperature + sutherland) / reynolds;
  close(p.viscosity, viscosity); close(p.reTheta, rho * input.ue * input.theta / viscosity);
});

test('full native single-element state satisfies independently assembled BL, transition and TE/wake blocks', () => {
  const { bl } = runCoupled(naca4('0012',160), { alpha: 4, reynolds: 1e6 });
  const kernel = createIntegralKernel({ reynolds: 1e6 });
  const at = (i, side) => ({ s: bl.XSSI[i][side], theta: bl.THET[i][side], deltaStar: bl.DSTR[i][side],
    ue: bl.UEDG[i][side], aux: bl.CTAU[i][side], wakeGap: i > bl.IBLTE[side] ? bl.WGAP[i - bl.IBLTE[side]] : 0 });
  const errors = { similarity: 0, laminar: 0, transition: 0, turbulent: 0, wake: 0, te: 0 };
  for (const side of [1,2]) for (let i = 2; i <= bl.NBL[side]; i++) {
    const downstream = at(i,side);
    if (i === bl.IBLTE[side] + 1) {
      const result = kernel.trailingEdge(at(bl.IBLTE[1],1), at(bl.IBLTE[2],2), downstream);
      errors.te = Math.max(errors.te, ...result.residual.map(Math.abs)); continue;
    }
    const regime = i === 2 ? 'similarity' : i > bl.IBLTE[side] ? 'wake' : i < bl.ITRAN[side] ? 'laminar' : i === bl.ITRAN[side] ? 'transition' : 'turbulent';
    const result = kernel.interval({ upstream: i === 2 ? downstream : at(i-1,side), downstream, regime });
    errors[regime] = Math.max(errors[regime], ...result.residual.map(Math.abs));
  }
  for (const [regime, error] of Object.entries(errors)) assert.ok(error < 3e-4, `${regime}: residual ${error}`);
});
