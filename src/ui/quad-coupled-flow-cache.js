// SPDX-License-Identifier: GPL-2.0-or-later
// One accepted coupled checkpoint survives Mach edits. This module imports
// no PDE, meshing or WASM code and never retains a displayed flow solution.

const positive = value => Number.isFinite(value) && value > 0;
const object = value => value !== null && typeof value === 'object' && !Array.isArray(value)
  && (Object.getPrototypeOf(value) === Object.prototype || Object.getPrototypeOf(value) === null);
const vector = value => Array.isArray(value) || ArrayBuffer.isView(value) && !(value instanceof DataView);
const finiteVector = value => vector(value) && value.length > 0 && Array.from(value).every(Number.isFinite);
const families = ['euler', 'boundaryLayer', 'edgeMatching'];

function equal(a, b) {
  if (Object.is(a, b)) return true;
  if (typeof a !== typeof b || a === null || b === null || typeof a !== 'object') return false;
  if (Array.isArray(a) || Array.isArray(b)) {
    if (!Array.isArray(a) || !Array.isArray(b) || a.length !== b.length) return false;
  } else if (!object(a) || !object(b)) return false;
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && equal(a[k], b[k]));
}

function matchingCase(a, b) {
  if (!object(a) || !object(b)) return false;
  const ka = Object.keys(a).filter(k => k !== 'mach'), kb = Object.keys(b).filter(k => k !== 'mach');
  return ka.length === kb.length && ka.every(k => Object.hasOwn(b, k) && equal(a[k], b[k]));
}

function nodesMatch(a, b) {
  return Array.isArray(a) && Array.isArray(b) && a.length > 0 && a.length === b.length
    && a.every((group, g) => Array.isArray(group) && Array.isArray(b[g]) && group.length > 0 && group.length === b[g].length
      && group.every((row, i) => Array.isArray(row) && Array.isArray(b[g][i]) && row.length > 0 && row.length === b[g][i].length
        && row.every((point, j) => [point?.x, point?.y, b[g][i][j]?.x, b[g][i][j]?.y].every(Number.isFinite))));
}

// Initialization/refinement records can themselves contain full histories
// and geometry. Retain descriptive scalar records and short scalar lists,
// not those secondary copies of solved data. The checkpoint stays complete.
const bulky = new Set(['flow', 'x', 'residual', 'profiles', 'stations', 'history', 'nodes', 'undisplacedNodes',
  'vertices', 'checkpoint', 'restart', 'initialEuler', 'initialBL']);
function metadata(value, key = '') {
  if (value === null || ['string', 'boolean', 'undefined'].includes(typeof value)) return value;
  if (typeof value === 'number') return Number.isFinite(value) ? value : undefined;
  if (bulky.has(key)) return undefined;
  if (Array.isArray(value)) return value.length <= 32 && value.every(v => v === null || ['string', 'boolean'].includes(typeof v)
    || typeof v === 'number' && Number.isFinite(v)) ? value.slice() : undefined;
  if (!object(value)) return undefined;
  return Object.fromEntries(Object.entries(value).flatMap(([k, v]) => {
    const copy = metadata(v, k); return copy === undefined ? [] : [[k, copy]];
  }));
}

function eligible(result) {
  const sourceCase = result?.sourceCase, checkpoint = result?.checkpoint, f = checkpoint?.restart;
  const sequence = result?.gridSequence
    ?? (result?.automaticRefinement?.kind === 'coarse-to-fine' ? result.automaticRefinement : undefined);
  const tolerance = result?.solverSettings?.tolerance;
  return result?.model === 'research-streamtube-euler-bl' && result.converged === true
    && result.mesh?.quality?.valid === true && object(sourceCase) && sourceCase.flowModel === 'streamtube-grid'
    && sourceCase.quadBoundaryLayers === true && positive(sourceCase.mach) && sourceCase.mach < 1
    && checkpoint?.version === 1 && object(f?.input) && object(f?.options) && object(checkpoint.continuation)
    && sourceCase.mach === f.input.mach && ['isentropic', 'hybrid'].includes(f.input.streamwiseMode)
    && (sourceCase.ncrit ?? 9) === (f.options.ncrit ?? 9)
    && result.ncritContinuation?.reachedTarget !== false
    && result.gridSequence?.reachedTarget !== false
    && !(result.automaticRefinement?.kind === 'coarse-to-fine' && result.automaticRefinement.reachedTarget === false)
    && (!Number.isInteger(sequence?.actualGridIntervals)
      || sequence.actualGridIntervals === sourceCase.gridIntervals
        && (sequence.requestedGridIntervals === undefined
          || sequence.actualGridIntervals === sequence.requestedGridIntervals))
    && (f.input.streamwiseMode !== 'hybrid' || f.options.blThermodynamics === 'historical-common-isentrope')
    && f.options.edgeMatching === 'section-velocity'
    && finiteVector(f.initialEuler?.x) && finiteVector(f.initialBL) && f.initialBL.length % 4 === 0
    && nodesMatch(f.initialEuler.nodes, f.initialEuler.undisplacedNodes)
    && object(checkpoint.families) && object(result.families)
    && Object.keys(checkpoint.families).length === 3 && Object.keys(result.families).length === 3
    && positive(tolerance) && families.every(k => Number.isFinite(checkpoint.families[k])
      && checkpoint.families[k] >= 0 && checkpoint.families[k] <= tolerance
      && checkpoint.families[k] === result.families[k])
    && ['referenceChord', 'referenceReynolds', 'solverLength', 'kernelReynolds'].every(k => positive(result[k]));
}

export function createQuadCoupledFlowCache() {
  let parent = null;
  return {
    remember(result) {
      try {
        if (!eligible(result)) return false;
        const next = { model: result.model, converged: true, mesh: { quality: structuredClone(result.mesh.quality) },
          sourceCase: structuredClone(result.sourceCase), checkpoint: structuredClone(result.checkpoint),
          families: structuredClone(result.families), referenceChord: result.referenceChord,
          referenceReynolds: result.referenceReynolds, solverLength: result.solverLength, kernelReynolds: result.kernelReynolds,
          solverSettings: metadata(result.solverSettings), initialization: metadata(result.initialization) ?? {},
          ...(result.refinement ? { refinement: metadata(result.refinement) ?? {} } : {}) };
        parent = next; return true;
      } catch {
        // A malformed/uncloneable new result must not displace a good seed.
        return false;
      }
    },
    forCase(caseData) {
      if (!parent || !positive(caseData?.mach) || caseData.mach >= 1 || !matchingCase(parent.sourceCase, caseData)) return null;
      return structuredClone(parent);
    },
    clear() { parent = null; },
  };
}
