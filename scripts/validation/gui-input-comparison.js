// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';

// Node and Chrome can differ in the last bits of Math.cos/sin when generating
// NACA coordinates. Only generated coordinates receive a roundoff allowance;
// physical controls, metadata and supplied benchmark coordinates remain exact.
export function compareGuiInput(actual, expected) {
  const supplied = Boolean(expected.geometrySource);
  let maxCoordinateDifference = 0, comparedCoordinates = 0;
  const withoutCoordinates = input => ({ ...input, elements: input.elements.map(({ points, ...element }) => element) });
  assert.deepEqual(withoutCoordinates(actual), withoutCoordinates(expected));
  for (let e = 0; e < expected.elements.length; e++) {
    const a = actual.elements[e].points, b = expected.elements[e].points;
    assert.equal(a.length, b.length);
    for (let i = 0; i < b.length; i++) {
      assert.deepEqual(Object.keys(a[i]).sort(), Object.keys(b[i]).sort());
      for (const key of ['x', 'y']) {
        const difference = Math.abs(a[i][key] - b[i][key]);
        const tolerance = supplied ? 0 : 8 * Number.EPSILON * Math.max(1, Math.abs(b[i][key]));
        assert.ok(Number.isFinite(difference) && difference <= tolerance,
          `Element ${e}, point ${i}, ${key}: coordinate difference ${difference} exceeds roundoff allowance ${tolerance}.`);
        maxCoordinateDifference = Math.max(maxCoordinateDifference, difference);
        comparedCoordinates++;
      }
    }
  }
  return { suppliedCoordinatesExact: supplied, comparedCoordinates, maxCoordinateDifference,
    generatedCoordinateRelativeAllowance: supplied ? 0 : 8 * Number.EPSILON, controlsExact: true };
}
