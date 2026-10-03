// SPDX-License-Identifier: GPL-2.0-or-later
import { Worker } from 'node:worker_threads';
export function runBoundedWorker(workerUrl, data, { timeoutMs = 180000, onProgress } = {}) {
  if (!(Number.isFinite(timeoutMs) && timeoutMs > 0)) throw new Error('Supply a positive worker time limit.');
  return new Promise(resolve => {
    const started = performance.now(), worker = new Worker(workerUrl, { workerData: data });
    const history = [], retained = {};
    let settled = false;
    const finish = async terminal => {
      if (settled) return;
      settled = true; clearTimeout(timer);
      await worker.terminate();
      resolve({ ...terminal, seconds: (performance.now() - started) / 1000, history, retained });
    };
    const timer = setTimeout(() => void finish({ type: 'timeout', message: `Exceeded ${timeoutMs / 1000} seconds.`,
      stage: retained.stage ?? data.stage }), timeoutMs);
    worker.on('message', message => {
      if (settled) return;
      if (['result', 'error', 'geometry-ready', 'mesh-ready'].includes(message.type)) { void finish(message); return; }
      if (message.type === 'mesh') retained.mesh = message.mesh;
      if (message.type === 'accepted-flow') retained.acceptedFlow = message.parentResult;
      if (message.type === 'coefficients') retained.coefficients = message.coefficients;
      if (message.type === 'flow-stage') retained.stage = message.stage;
      if (message.type === 'iteration') history.push({ ...message.iteration, seconds: (performance.now() - started) / 1000 });
      if (['flow-stage', 'iteration', 'mesh-unavailable', 'coefficients-unavailable'].includes(message.type)) {
        try { onProgress?.(message); }
        catch (error) { void finish({ type: 'error', code: 'PROGRESS_CAPTURE_FAILED', message: error.message, stage: retained.stage ?? data.stage }); }
      }
    });
    worker.on('error', error => void finish({ type: 'error', message: error.message, stack: error.stack, stage: retained.stage ?? data.stage }));
    worker.on('exit', code => { if (!settled) void finish({ type: 'error', message: `Worker exited without a terminal result (${code}).` }); });
  });
}
