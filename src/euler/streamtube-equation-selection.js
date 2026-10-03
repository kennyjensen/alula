// SPDX-License-Identifier: GPL-2.0-or-later
// MSES 3.05 manual §§1.2.4–1.2.5: selectable S-momentum/entropy rows.
// The existing `hybrid` envelope keeps the shared upwind gas, historical BL
// adapter and complete restart format. Omitted ismom retains the old blend.
// ISMOM3 uses an explicitly documented logical interpretation of the manual's
// approximate region: from the inlet through leadingIndex+10 inclusive, and
// four tube slots on either side of each body dividing streamline. Overlapping
// body regions form a union. The unpublished SETUP DATA bounds are not claimed
// identical; this region is fixed within a grid and has no state derivative.
// ISMOM4 retains our existing documented blend, not an exact native switch.

export const ISMOM3_LEADING_REGION = Object.freeze({ downstreamCells: 10, transverseCells: 4 });

// Explicit inherited regions are a refinement research control. Their logical
// bounds and source checkpoint travel with the equations; the default 10/4
// interpretation above remains unchanged. Bounds do not follow the moving
// Cartesian grid and therefore have no state derivative.
const exactKeys = (object, keys) => object && typeof object === 'object' && !Array.isArray(object)
  && Object.keys(object).length === keys.length && keys.every(key => Object.hasOwn(object, key));
const positiveInteger = value => Number.isSafeInteger(value) && value > 0;

export function normalizeStreamtubeEquationRegions(value) {
  if (!exactKeys(value, ['version', 'parentCheckpointSha256', 'topology', 'regions']) || value.version !== 1
    || typeof value.parentCheckpointSha256 !== 'string' || !/^[0-9a-f]{64}$/.test(value.parentCheckpointSha256))
    throw new Error('ISMOM3 inherited entropy regions require version 1 and a parent checkpoint SHA-256.');
  const t = value.topology;
  if (!exactKeys(t, ['nx', 'tubes', 'bodies']) || !Number.isSafeInteger(t.nx) || t.nx < 4
    || !Array.isArray(t.bodies) || !t.bodies.length || !Array.isArray(t.tubes)
    || t.tubes.length !== t.bodies.length + 1 || !Array.from(t.tubes).every(positiveInteger)
    || Array.from(t.bodies).some(b => !exactKeys(b, ['leadingIndex', 'trailingIndex'])
      || !positiveInteger(b.leadingIndex) || !positiveInteger(b.trailingIndex)
      || b.leadingIndex >= b.trailingIndex || b.trailingIndex >= t.nx))
    throw new Error('Invalid ISMOM3 inherited entropy-region topology.');
  const total = t.tubes.reduce((sum, n) => sum + n, 0);
  let cut = 0;
  if (!Number.isSafeInteger(total) || !Array.isArray(value.regions) || value.regions.length !== t.bodies.length
    || Array.from(value.regions).some((r, b) => {
      cut += t.tubes[b];
      return !exactKeys(r, ['body', 'throughRow', 'lowerTube', 'upperTube']) || r.body !== b
        || !positiveInteger(r.throughRow) || r.throughRow < t.bodies[b].leadingIndex || r.throughRow >= t.nx
        || !Number.isSafeInteger(r.lowerTube) || r.lowerTube < 0 || r.lowerTube >= cut
        || !Number.isSafeInteger(r.upperTube) || r.upperTube <= cut || r.upperTube > total;
    })) throw new Error('Invalid ISMOM3 inherited entropy-region bounds.');
  return Object.freeze({ version: 1, parentCheckpointSha256: value.parentCheckpointSha256,
    topology: Object.freeze({ nx: t.nx, tubes: Object.freeze(t.tubes.slice()),
      bodies: Object.freeze(t.bodies.map(b => Object.freeze({ ...b }))) }),
    regions: Object.freeze(value.regions.map(r => Object.freeze({ ...r }))) });
}

export function validateStreamtubeEquationRegionTopology(value, { nx, tubes, bodies }) {
  const t = value.topology;
  if (nx !== t.nx || !Array.isArray(tubes) || tubes.length !== t.tubes.length
    || Array.from(tubes).some((n, g) => n !== t.tubes[g]) || !Array.isArray(bodies) || bodies.length !== t.bodies.length
    || Array.from(bodies).some((b, k) => !b || b.leadingIndex !== t.bodies[k].leadingIndex || b.trailingIndex !== t.bodies[k].trailingIndex))
    throw new Error('ISMOM3 inherited entropy-region topology does not match this grid; explicitly remap it before refinement.');
}

export function normalizeStreamtubeEquationSelection(hybrid) {
  if (!hybrid || typeof hybrid !== 'object' || !Number.isFinite(hybrid.epsilonP) || hybrid.epsilonP <= 0)
    throw new Error('Hybrid body rows require explicit upwinding and positive epsilonP.');
  if (hybrid.ismom !== undefined && (!Number.isInteger(hybrid.ismom) || hybrid.ismom < 1 || hybrid.ismom > 4))
    throw new Error('ISMOM must be an integer from 1 to 4.');
  if (hybrid.entropyRegions !== undefined && hybrid.ismom !== 3)
    throw new Error('Inherited entropy regions require explicit ISMOM3.');
  return Object.freeze({ epsilonP: hybrid.epsilonP, ...(hybrid.ismom === undefined ? {} : { ismom: hybrid.ismom }),
    ...(hybrid.entropyRegions === undefined ? {} : { entropyRegions: normalizeStreamtubeEquationRegions(hybrid.entropyRegions) }) });
}

export function streamtubeEquationAt({ hybrid, bodies, tubes, nx, i, group, tube }) {
  const ismom = hybrid?.ismom ?? 4;
  if (ismom === 1) return 'momentum';
  if (ismom === 2) return 'isentropic';
  if (ismom === 4) return 'hybrid';
  if (ismom !== 3) throw new Error('ISMOM must be an integer from 1 to 4.');
  if (!Array.isArray(bodies) || bodies.length === 0 || !Array.isArray(tubes) || tubes.length !== bodies.length + 1
    || !tubes.every(n => Number.isInteger(n) && n > 0) || !Number.isInteger(i) || i < 1
    || !Number.isInteger(group) || group < 0 || group >= tubes.length
    || !Number.isInteger(tube) || tube < 0 || tube >= tubes[group]
    || bodies.some(b => !Number.isInteger(b.leadingIndex) || b.leadingIndex < 1))
    throw new Error('ISMOM3 requires valid logical body/cell indices.');
  let globalTube = tube;
  for (let g = 0; g < group; g++) globalTube += tubes[g];
  if (hybrid.entropyRegions !== undefined) {
    validateStreamtubeEquationRegionTopology(hybrid.entropyRegions, { nx, tubes, bodies });
    if (i >= nx) throw new Error('Invalid ISMOM3 inherited entropy-region row.');
    return hybrid.entropyRegions.regions.some(r => i <= r.throughRow && globalTube >= r.lowerTube && globalTube < r.upperTube)
      ? 'isentropic' : 'momentum';
  }
  let cut = 0;
  for (let b = 0; b < bodies.length; b++) {
    cut += tubes[b];
    if (i <= bodies[b].leadingIndex + ISMOM3_LEADING_REGION.downstreamCells
      && globalTube >= cut - ISMOM3_LEADING_REGION.transverseCells
      && globalTube < cut + ISMOM3_LEADING_REGION.transverseCells) return 'isentropic';
  }
  return 'momentum';
}

// Public explicit ISMOM controls. Omission preserves each adapter's defaults.
export function streamtubeEquationControls(ismom) {
  if (ismom === undefined) return {};
  if (!Number.isInteger(ismom) || ismom < 1 || ismom > 4)
    throw new Error('Euler ISMOM must be an integer from 1 to 4, or omitted for the default equations.');
  return { streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5, ismom },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
}
