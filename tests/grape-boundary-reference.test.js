import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { reconstructOrthogonalBoundary } from '../src/geometry/orthogonal-boundary-control.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/grape-boundary.json', import.meta.url)));
const close = (a, b, label, tolerance = 2e-10) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(b)), `${label}: ${a} != ${b}`);

test('orthogonal Q=0 reconstruction agrees with original native GRAPE inner and outer boundary equations', t => {
  const source = fs.readFileSync(new URL('../third_party/grape/original.src', import.meta.url), 'utf8');
  assert.equal(createHash('sha256').update(source).digest('hex'), fixture.sourceSha256);
  for (const block of Object.values(fixture.extractedBlocks)) assert.ok(source.includes(block));
  assert.equal(fixture.cases.length, 40);
  let maximumRelativePoissonError = 0, maximumNativeQ = 0;
  for (const c of fixture.cases) {
    const result = reconstructOrthogonalBoundary(c.arguments), expected = c.expected;
    close(result.normalSpeed, c.normalSpeed, `${c.name}: inward speed`);
    close(result.stretch / result.gamma, expected.poisson, `${c.name}: P`);
    close(expected.transversePoisson, 0, `${c.name}: native Q at selected speed`, 2e-11);
    for (const key of ['x', 'y']) close(result.normalSecond[key], expected.normalSecond[key], `${c.name}: normal second ${key}`);
    maximumNativeQ = Math.max(maximumNativeQ, Math.abs(expected.transversePoisson));
    maximumRelativePoissonError = Math.max(maximumRelativePoissonError,
      Math.abs(result.stretch / result.gamma - expected.poisson) / Math.max(1, Math.abs(expected.poisson)));
  }
  t.diagnostic(JSON.stringify({ nativeCases: fixture.cases.length, maximumRelativePoissonError, maximumNativeQ }));
});

test('native boundary fixture covers unequal interior distances and both normal orientations in the failing passage', () => {
  const manufactured = fixture.cases.filter(c => c.name.startsWith('cubic'));
  const passage = fixture.cases.filter(c => c.name.startsWith('lower passage'));
  assert.equal(manufactured.length, 18); assert.equal(passage.length, 22);
  for (const cases of [manufactured, passage]) for (const sign of [-1, 1]) {
    assert.ok(cases.some(c => c.arguments.normalSign === sign));
    assert.ok(cases.some(c => c.arguments.normalSign === sign
      && Math.abs(c.arguments.secondDistance / c.arguments.firstDistance - 2) > .1));
  }
});
