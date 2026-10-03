import { test, expect } from '@playwright/test';
import { conformalPolynomialGridFixture } from '../fixtures/polynomial-grid.js';

test('certified curved harmonic grid and independent WASM reference execute in a browser module worker', async ({ page }) => {
  // An empty same-origin harness isolates module compatibility from the
  // application startup solve; this test does not claim GUI integration.
  await page.route('**/curved-kernel-check', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Curved kernel check</title>' }));
  await page.goto('/curved-kernel-check');
  const data = conformalPolynomialGridFixture({ nx: 4, nt: 3 });
  const input = { nodes: data.nodes, controlPoints: data.controlPoints, directions: data.directions, massFlows: data.massFlows };
  const result = await page.evaluate(input => new Promise((resolve, reject) => {
    const geometryURL = new URL('/src/geometry/polynomial-grid-geometry.js', location.href).href;
    const solverURL = new URL('/src/geometry/transverse-harmonic-grid.js', location.href).href;
    const referenceURL = new URL('/src/geometry/tests/curved-harmonic-reference.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createPolynomialGridGeometry}=await import(${JSON.stringify(geometryURL)});
      const {createTransverseHarmonicGrid,smoothTransverseHarmonicGrid}=await import(${JSON.stringify(solverURL)});
      const {solveCurvedHarmonicReference}=await import(${JSON.stringify(referenceURL)});
      const geometry=createPolynomialGridGeometry(data),nx=data.nodes.length-1,nt=data.massFlows.length;
      const nodes=data.nodes.map((row,i)=>row.map((p,j)=>({...p,y:p.y+(!i||!j||i===nx||j===nt?0:.025*Math.sin(Math.PI*i/nx)*Math.sin(Math.PI*j/nt))})));
      const result=smoothTransverseHarmonicGrid(createTransverseHarmonicGrid({...data,nodes,curvedGeometry:geometry}));
      if(!result.converged)throw new Error(result.reason);
      const reference=solveCurvedHarmonicReference({...data,nodes:result.nodes,geometry:geometry.onGrid(result.nodes)},{refinement:2});
      self.postMessage({nodes:result.nodes,converged:result.converged,valid:result.quality.valid,physicsValidated:result.physicsValidated,
        minimumJacobian:result.quality.minimumJacobian,minimumTransversality:result.quality.minimumTransversality,
        referenceError:reference.maximumTubeIntervals,linear:reference.linear});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(result.converged).toBe(true); expect(result.valid).toBe(true); expect(result.physicsValidated).toBe(false);
  expect(result.minimumJacobian).toBeGreaterThan(0); expect(result.minimumTransversality).toBeGreaterThan(.8);
  expect(result.referenceError).toBeLessThan(1e-9); expect(result.linear.relativeResidual).toBeLessThan(1e-10);
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    expect(Math.abs(data.psiAt(p) - j / data.nt)).toBeLessThan(1e-10);
    expect(p.x).toBe(data.nodes[i][j].x);
    if (!i || !j || i === data.nx || j === data.nt) expect(p).toEqual(data.nodes[i][j]);
  }));
});
