import { parentPort } from 'node:worker_threads';
parentPort.postMessage({ type: 'iteration', iteration: { iteration: 1, residual: 3 } });
for (;;) {} // Intentional: an iteration callback cannot interrupt this worker.
