import { test, expect } from '@playwright/test';
import { vortexChannel, vortexErrors, directChannelConservation } from '../oracles/streamtube.js';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directBodyConservation } from '../oracles/streamtube-body.js';
import { stagnationPressureError } from '../oracles/streamtube-entropy.js';

test('two-body trust-region relaxation transfers a recovered lifting grid to conservative Euler in a browser WASM worker', async ({ page }) => {
  await page.goto('/');
  const input = { ...intrinsicBodyFixture({ elements: 2, tubes: 5, alpha: 2 }), flowModel: 'incompressible', mach: 0 };
  const stages = await page.evaluate(input => new Promise((resolve, reject) => {
    const moduleUrl = new URL('/src/euler/streamtube-body.js', location.href).href;
    const restartUrl = new URL('/src/euler/tests/streamtube-body-restart.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createStreamtubeBodySystem,solveStreamtubeBody}=await import(${JSON.stringify(moduleUrl)});
      const {initializeStreamtubeBodyFromGrid}=await import(${JSON.stringify(restartUrl)});
      const system=createStreamtubeBodySystem(data);
      const result=solveStreamtubeBody(system,{stepMethod:'dogleg',maxIterations:25,tolerance:1e-11});
      if(!result.converged)throw new Error(result.reason);
      const target=initializeStreamtubeBodyFromGrid({...data,flowModel:'compressible',mach:.2,streamwiseMode:'momentum'},result);
      const euler=solveStreamtubeBody(target.system,{initial:target.initial,maxIterations:20,tolerance:1e-11});
      self.postMessage([{result,conditions:system.conditions},{result:euler,conditions:target.system.conditions}]);
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = error => { close(); reject(new Error(error.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(stages[0].result.stepMethod).toBe('dogleg');
  expect(stages[0].result.history.slice(1).some(h => h.stepKind !== 'newton')).toBe(true);
  expect(stages.map(s => s.conditions.mach)).toEqual([0, .2]);
  for (const { result, conditions } of stages) {
    expect(result.converged).toBe(true); expect(result.diagnostics.residual).toBeLessThan(1e-11);
    expect(result.linearDiagnostics.maxRelativeResidual).toBeLessThan(1e-10);
    const c = directBodyConservation(result, input.bodies, conditions);
    for (const k of result.streamwiseMode === 'momentum' ? [0, 1, 2, 3] : [0, 3]) expect(Math.abs(c.balance[k])).toBeLessThan(2e-9);
    expect(Math.max(...c.cutTraction.map(Math.abs))).toBeLessThan(2e-9);
    expect(result.surfaces).toHaveLength(4); expect(result.forceStatus).toContain('Unvalidated');
  }
  expect(stages[1].result.streamwiseMode).toBe('momentum');
  expect(stages[1].result.diagnostics.maxMomentumResidual).toBeLessThan(1e-11);
});

test('lifting two-body incompressible relaxation transfers to Mach 0.2 Euler in a browser WASM worker', async ({ page }) => {
  await page.goto('/');
  const input = { ...intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 16, surfaceSpacing: 'cosine', tubes: 5, tubeGrowth: 3 }), flowModel: 'incompressible', mach: 0 };
  const { stages, bodies } = await page.evaluate(input => new Promise((resolve, reject) => {
    const solverUrl = new URL('/src/euler/streamtube-body.js', location.href).href;
    const restartUrl = new URL('/src/euler/tests/streamtube-body-restart.js', location.href).href;
    const initializerUrl = new URL('/src/euler/streamtube-body-initializer.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createStreamtubeBodySystem,solveStreamtubeBody}=await import(${JSON.stringify(solverUrl)});
      const {initializeStreamtubeBodyFromGrid}=await import(${JSON.stringify(restartUrl)});
      const {initializePanelStreamtubeBody}=await import(${JSON.stringify(initializerUrl)});
      const prepared=initializePanelStreamtubeBody(data,{crosslinePlacement:'potential',outerCrosslineSpread:.15,normalSpacing:'stagnation'});
      const {system}=prepared,gridResult=solveStreamtubeBody(system,{initial:prepared.initial,tolerance:1e-11});
      if(!gridResult.converged)throw new Error(gridResult.reason);
      const target={...prepared.input,flowModel:'compressible',mach:.2,streamwiseMode:'isentropic'};
      const restarted=initializeStreamtubeBodyFromGrid(target,gridResult);
      const result=solveStreamtubeBody(restarted.system,{initial:restarted.initial,tolerance:1e-11});
      self.postMessage({bodies:prepared.input.bodies,stages:[{result:gridResult,conditions:system.conditions,n:system.layout.n},
        {result,conditions:restarted.system.conditions,n:restarted.system.layout.n}]});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = error => { close(); reject(new Error(error.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(stages.map(s => s.result.flowModel)).toEqual(['incompressible', 'compressible']);
  expect(stages.map(s => s.conditions.mach)).toEqual([0, .2]);
  expect(stages.map(s => s.n)).toEqual([1524, 2979]);
  for (const { result, conditions } of stages) {
    expect(result.converged).toBe(true);
    expect(result.linearBackend).toBe('klu'); expect(result.jacobianBackend).toBe('analytic');
    expect(result.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
    expect(result.diagnostics.residual).toBeLessThanOrEqual(1e-11);
    const c = directBodyConservation(result, bodies, conditions);
    for (const k of [0, 3]) expect(Math.abs(c.balance[k])).toBeLessThan(2e-9);
    expect(Math.max(...c.cutTraction.map(Math.abs))).toBeLessThan(2e-9);
    expect(result.surfaces).toHaveLength(4); expect(result.forceStatus).toContain('Unvalidated');
  }
  const { result, conditions } = stages[1];
  const p = stagnationPressureError(result.sections.map(row => row.flat()), { gamma: conditions.gamma,
    referencePressure: conditions.pInf, freestreamMach: .2 });
  expect(p.maxRelativeError).toBeLessThan(2e-10);
});

test('intrinsic moving-quadrilateral Euler channel solves in a browser module worker using analytic derivatives and WASM', async ({ page }) => {
  await page.goto('/');
  const oracle = vortexChannel(32, 8);
  const { exact, ordinate, ...input } = oracle;
  const result = await page.evaluate(input => new Promise((resolve, reject) => {
    const moduleUrl = new URL('/src/euler/tests/streamtube-channel.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createStreamtubeChannel,solveStreamtubeChannel}=await import(${JSON.stringify(moduleUrl)});
      const system=createStreamtubeChannel(data),result=solveStreamtubeChannel(system,{tolerance:1e-11});
      self.postMessage({result:{...result,unknowns:system.n}});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = error => { close(); reject(new Error(error.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data.result); };
    worker.postMessage(input);
  }), input);
  expect(result.converged).toBe(true);
  expect(result.unknowns).toBe(487);
  expect(result.jacobianBackend).toBe('analytic'); expect(result.linearBackend).toBe('klu');
  expect(result.linearDiagnostics.solves).toBeGreaterThan(0);
  expect(result.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(Math.max(...Array.from(result.residual, Math.abs))).toBeLessThan(1e-11);
  const errors = vortexErrors(result, oracle);
  expect(errors.pressureMax).toBeLessThan(.00015); expect(errors.positionMax).toBeLessThan(8e-6);
  expect(Math.max(...directChannelConservation(result).external.map(Math.abs))).toBeLessThan(2e-9);
});

for (const streamwiseMode of ['momentum', 'isentropic']) test(`two-body moving-streamtube ${streamwiseMode} equations solve in a browser worker with the complete analytic KLU backend`, async ({ page }) => {
  test.setTimeout(90000);
  await page.goto('/');
  const input = { ...intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 16,
    surfaceSpacing: 'cosine', tubes: 5, tubeGrowth: 3 }), streamwiseMode };
  const { result, conditions, unknowns, bodies, normalAllocation } = await page.evaluate(input => new Promise((resolve, reject) => {
    const solverUrl = new URL('/src/euler/streamtube-body.js', location.href).href;
    const initializerUrl = new URL('/src/euler/streamtube-body-initializer.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {solveStreamtubeBody}=await import(${JSON.stringify(solverUrl)});
      const {initializePanelStreamtubeBody}=await import(${JSON.stringify(initializerUrl)});
      const controls=data.streamwiseMode==='isentropic'
        ? {crosslinePlacement:'potential',outerCrosslineSpread:.15,normalSpacing:'stagnation'} : {};
      const {system,initial,input,diagnostics}=initializePanelStreamtubeBody(data,controls);
      const result=solveStreamtubeBody(system,{initial,tolerance:1e-10});
      self.postMessage({result,conditions:system.conditions,unknowns:system.layout.n,bodies:input.bodies,normalAllocation:diagnostics.normalAllocation});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = error => { close(); reject(new Error(error.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(result.converged).toBe(true);
  expect(unknowns).toBe(streamwiseMode === 'momentum' ? 1260 : 2979);
  expect(result.streamwiseMode).toBe(streamwiseMode);
  expect(result.jacobianBackend).toBe('analytic'); expect(result.linearBackend).toBe('klu');
  expect(result.linearDiagnostics.solves).toBeGreaterThan(0);
  expect(result.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(Math.max(...Array.from(result.residual, Math.abs))).toBeLessThan(1e-10);
  expect(result.surfaces).toHaveLength(4); expect(result.diagnosticForces).toHaveLength(2);
  expect(result.forceStatus).toContain('Unvalidated');
  const conservation = directBodyConservation(result, bodies, conditions);
  if (streamwiseMode === 'momentum') expect(Math.max(...conservation.balance.map(Math.abs))).toBeLessThan(2e-9);
  else {
    const targets = normalAllocation.groups.flatMap(g => Object.values(g.targets));
    expect(targets).toHaveLength(4);
    targets.forEach(t => expect(Math.abs(t.actualAspect / 2.5 - 1)).toBeLessThan(1.1e-4));
    for (const k of [0, 3]) expect(Math.abs(conservation.balance[k])).toBeLessThan(2e-9);
    expect(Math.max(...conservation.cutTraction.map(Math.abs))).toBeLessThan(2e-9);
    const pressure = stagnationPressureError(result.sections.map(row => row.flat()), {
      gamma: conditions.gamma, referencePressure: conditions.pInf, freestreamMach: input.mach });
    expect(pressure.maxRelativeError).toBeLessThan(2e-10);
  }
});
