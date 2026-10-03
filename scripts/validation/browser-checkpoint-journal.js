// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import path from 'node:path';

const serial = value => JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v) + '\n';

// Atomic replacement leaves either the previous complete JSON or the next
// complete JSON after interruption. A flushed checkpoint precedes its index.
export function writeAtomicJson(directory, name, value) {
  const target = path.join(directory, name + '.json'), temporary = target + '.pending';
  const descriptor = fs.openSync(temporary, 'w');
  try { fs.writeFileSync(descriptor, serial(value)); fs.fsyncSync(descriptor); }
  finally { fs.closeSync(descriptor); }
  fs.renameSync(temporary, target);
}

export function createBrowserCheckpointJournal(directory, { archiveEvery = 10 } = {}) {
  if (!Number.isInteger(archiveEvery) || archiveEvery < 1) throw new Error('Invalid archive observation interval.');
  fs.mkdirSync(path.join(directory, 'checkpoints'), { recursive: true });
  let lastSequence = 0, checkpointCount = 0, terminal, latest;
  const record = event => {
    if (!Number.isInteger(event.sequence) || event.sequence !== lastSequence + 1)
      throw new Error(`Missing or reordered diagnostic event: ${event.sequence} after ${lastSequence}.`);
    const sequence = String(event.sequence).padStart(6, '0');
    if (event.type === 'coupled-checkpoint') {
      if (!event.checkpoint?.restart || !event.checkpoint?.continuation)
        throw new Error('Diagnostic coupled checkpoint is incomplete.');
      const final = event.details?.history?.at(-1), iteration = final?.iteration;
      const packet = { checkpoint: event.checkpoint, details: event.details,
        observation: { sequence: event.sequence, workerSeconds: event.workerSeconds } };
      writeAtomicJson(directory, 'latest-checkpoint', packet);
      if (!checkpointCount++) writeAtomicJson(directory, 'initial-checkpoint', packet);
      if (iteration === 0 || Number.isInteger(iteration) && iteration % archiveEvery === 0 || event.details?.retained)
        writeAtomicJson(path.join(directory, 'checkpoints'), sequence, packet);
      latest = { sequence: event.sequence, stage: event.details?.stage, startupAttempt: event.details?.startupAttempt,
        iteration, iterationOffset: event.details?.iterationOffset, actualNcrit: event.checkpoint.restart.options?.ncrit,
        families: event.checkpoint.families, workerSeconds: event.workerSeconds };
    } else if (event.type === 'assembly-checkpoint') {
      writeAtomicJson(directory, `assembly-${sequence}`, event);
      writeAtomicJson(directory, 'latest-assembly-checkpoint', event);
    } else if (event.type === 'prepared-euler') {
      writeAtomicJson(directory, `prepared-euler-${sequence}`, event);
    } else if (event.type === 'mesh') {
      writeAtomicJson(directory, 'latest-mesh', event);
    } else if (event.type === 'result' || event.type === 'error') {
      writeAtomicJson(directory, 'terminal', event); terminal = event;
      if (event.type === 'result') writeAtomicJson(directory, 'result', event.result);
    }
    // Large state arrays are separate, complete JSON files. This compact log
    // records every event and every accepted iteration without duplicating them.
    const compact = event.type === 'iteration' ? event : { sequence: event.sequence, type: event.type,
      workerSeconds: event.workerSeconds, stage: event.stage, startupAttempt: event.startupAttempt,
      ...(event.type === 'coupled-checkpoint' ? { checkpoint: latest } : {}),
      ...(event.type === 'mesh' ? { cells: event.mesh.cells.length, quality: event.mesh.quality } : {}),
      ...(event.message ? { message: event.message } : {}) };
    fs.appendFileSync(path.join(directory, 'progress.jsonl'), serial(compact));
    lastSequence = event.sequence;
    writeAtomicJson(directory, 'journal', { lastSequence, checkpointCount, latest,
      terminal: terminal?.type ?? null, scope: 'Diagnostic observation, not GUI or physical acceptance.' });
    return { lastSequence, checkpointCount, latest, terminal };
  };
  return { record, status: () => ({ lastSequence, checkpointCount, latest, terminal }) };
}
