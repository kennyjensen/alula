import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4, prepareContour, signedArea, transform, validateAssembly } from '../src/geometry/airfoil.js';
import { parseCoordinates } from '../src/geometry/parse.js';

test('NACA geometry is closed, CCW, symmetric and uses a sharp trailing edge', () => {
  const p = naca4('0012', 80);
  assert.equal(p.length, 81); assert.deepEqual(p[0], p.at(-1));
  assert.ok(signedArea(p) > 0);
  for (let i = 0; i <= 80; i++) { assert.equal(p[i].x, p[80 - i].x); assert.ok(Math.abs(p[i].y + p[80 - i].y) < 1e-14); }
  assert.deepEqual(prepareContour(p.toReversed()), prepareContour(p));
});
test('invalid geometry is rejected, including collisions and containment', () => {
  assert.throws(() => naca4('2012'), /Invalid/);
  const p = naca4();
  const open = p.map(x => ({ ...x })); open.at(-1).y = 0.001;
  assert.throws(() => prepareContour(open), /sharp trailing/);
  const duplicate = [...p]; duplicate[5] = duplicate[4];
  assert.throws(() => prepareContour(duplicate), /duplicate|intersects/);
  assert.throws(() => validateAssembly([p, transform(p, { x: 0.5 })]), /intersect/);
  assert.throws(() => validateAssembly([p, transform(p, { chord: 0.1, x: 0.4 })]), /inside/);
});
test('Selig and MSES parsing preserves shared coordinates, domain and element order', () => {
  const points = naca4('0012', 20);
  const rows = points.map(p => `${p.x} ${p.y}`).join('\n');
  const parsed = parseCoordinates(`Two elements\n-2 4 -3 3\n${rows}\n999. 999.\n${rows}\n`);
  assert.equal(parsed.elements.length, 2); assert.equal(parsed.domain.xMin, -2);
  assert.deepEqual(parsed.elements[0].points, points);
  assert.equal(parseCoordinates(rows).elements.length, 1);
  assert.throws(() => parseCoordinates('Name\n1 2 3'), /Invalid coordinate/);
});
