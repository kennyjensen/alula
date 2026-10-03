// SPDX-License-Identifier: GPL-2.0-or-later
self.onmessage = async ({ data: { id, controls } }) => {
  try {
    const [{buildEulerLabCase},{solveEuler},{solveStreamlineChannel}]=await Promise.all([
      import('../euler/cases.js'),import('../euler/solve.js'),import('../euler/streamline.js')]);
    const start = performance.now(); const problem = buildEulerLabCase(controls);
    const options = { ...problem.options,
      onIteration: iteration => self.postMessage({ id, type: 'iteration', iteration }),
      onMesh: state => self.postMessage({ id, type: 'mesh', ...state }) };
    const result = (problem.moving ? solveStreamlineChannel : solveEuler)(problem.mesh, options);
    self.postMessage({ id, type: 'result', result, elapsed: performance.now() - start });
  } catch (error) { self.postMessage({ id, type: 'error', message: error.message }); }
};
