// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { compareGuiInput } from '../scripts/validation/gui-input-comparison.js';
import { buildGuiDefaultCoupledCase } from '../scripts/validation/gui-default-coupled-cases.js';

test('generated browser coordinates admit only platform roundoff; controls remain exact', () => {
  const expected = buildGuiDefaultCoupledCase('single'), actual = structuredClone(expected);
  actual.elements[0].points[17].x += Number.EPSILON;
  const before = structuredClone(actual);
  assert.equal(compareGuiInput(actual, expected).maxCoordinateDifference, Number.EPSILON);
  assert.deepEqual(actual, before);
  actual.mach += Number.EPSILON;
  assert.throws(() => compareGuiInput(actual, expected));
  actual.mach = expected.mach;
  actual.elements[0].points[17].x += 1e-12;
  assert.throws(() => compareGuiInput(actual, expected), /roundoff allowance/);
});

test('benchmark coordinates and point counts remain exact across browser and Node', () => {
  const expected = buildGuiDefaultCoupledCase('nlr7301'), actual = structuredClone(expected);
  assert.equal(compareGuiInput(actual, expected).suppliedCoordinatesExact, true);
  actual.elements[0].points[17].x += Number.EPSILON;
  assert.throws(() => compareGuiInput(actual, expected), /roundoff allowance/);
  actual.elements[0].points = expected.elements[0].points.slice(1);
  assert.throws(() => compareGuiInput(actual, expected));
});
