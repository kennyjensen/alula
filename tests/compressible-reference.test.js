import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';

test('finite-Mach whole-flow reference retains original Fortran provenance and converged native outputs',()=>{
  const fixture=JSON.parse(readFileSync(new URL('./fixtures/fortran/compressible-flow.json',import.meta.url)));
  for(const [path,hash] of Object.entries(fixture.provenance.sha256))assert.equal(createHash('sha256').update(readFileSync(path)).digest('hex'),hash,path);
  assert.deepEqual(fixture.cases.map(c=>c.options.mach),[0,.1,.2]);
  for(const c of fixture.cases){
    assert.equal(c.expected.converged,true);assert.equal(c.panels,320);
    for(const key of ['cl','cd','cm'])assert.ok(Number.isFinite(c.expected[key]),key);
    assert.ok(c.expected.cp.length>=c.panels);
  }
});
