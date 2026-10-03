// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { parseCoordinates } from '../src/geometry/parse.js';
import { serializeCoordinates } from '../src/geometry/serialize.js';
import { prepareAirfoilElement } from '../src/geometry/airfoil-element.js';
import { createContourTopology } from '../src/geometry/contour-topology.js';
import { naca4, prepareContour, transform } from '../src/geometry/airfoil.js';

const rowText = points => points.map(p => `${p.x} ${p.y}`).join('\n');
const open = () => naca4('0012', 20).map((p, i) => ({ x: p.x, y: p.y + (i < 10 ? 1 : -1) * .002 * p.x }));
const domain = { xMin: -2, xMax: 3, yMin: -3, yMax: 3.5 };

test('plain/labeled XFOIL parsing preserves orientation, open TE endpoints and coordinate frame', () => {
  const points = transform(open(), { chord: 2.3, angle: 17, x: 3.4, y: -1.7 }).reverse();
  for (const prefix of ['', '0012 NACA example\n']) {
    const input = parseCoordinates(prefix + rowText(points));
    assert.deepEqual(input.elements[0].points, points);
    assert.equal(input.format, prefix ? 'xfoil-labeled' : 'xfoil-plain');
    assert.equal(input.domain, null);
    assert.equal(Object.hasOwn(input.elements[0], 'trailingEdge'), false);
    assert.deepEqual(parseCoordinates(serializeCoordinates(input)).elements[0].points, points);
  }
});

test('BOM, CRLF, Fortran D exponents, commas and comments retain numerical coordinates', () => {
  const points = open(), rows = points.map(p => `${p.x.toExponential(17).replace('e', 'D')}, ${p.y.toExponential(17).replace('e', 'd')} ! ordinate`);
  const input = parseCoordinates('\uFEFF# leading comment\r\nFoil # label text\r\n' + rows.join('\r\n') + '\r\n');
  assert.deepEqual(input.elements[0].points, points);
  assert.equal(input.name, 'Foil # label text');
  assert.deepEqual(parseCoordinates(serializeCoordinates(input)).elements[0].points, points);
});

test('MSES header, separators, element order and legacy header values round-trip', () => {
  const a = open(), b = transform(open(), { chord: .3, x: 1.1, y: -.08 });
  const body = `${rowText(a)}\n999.0 999.0\n${rowText(b)}\n999 999\n`;
  const current = parseCoordinates('Assembly\n-2 3 -3 3.5\n' + body);
  assert.deepEqual(current.domain, domain); assert.equal(current.format, 'mses');
  assert.deepEqual(current.elements.map(e => e.points), [a, b]);
  const replay = parseCoordinates(serializeCoordinates(current));
  assert.deepEqual(replay.domain, current.domain); assert.deepEqual(replay.elements, current.elements);
  const legacy = parseCoordinates('Old ISES\n-2 3 -2.5 3 0.2 7\n' + rowText(a));
  assert.equal(legacy.domain, null); assert.deepEqual(legacy.header, { kind: 'legacy-ises-mses', values: [-2, 3, -2.5, 3, .2, 7] });
  assert.match(legacy.warnings.join(' '), /not converted/);
  assert.deepEqual(parseCoordinates(serializeCoordinates(legacy)).header, legacy.header);
  const headerless = parseCoordinates(body);
  assert.equal(headerless.labeled, false);
  assert.deepEqual(parseCoordinates(serializeCoordinates(headerless)).elements, headerless.elements);
});

test('standard open TE chains become one explicit straight base without moving a source point', () => {
  for (const points of [open(), open().reverse()]) {
    const element = { name: 'Finite section', points, trips: [.2, .4] }, before = structuredClone(element);
    const prepared = prepareAirfoilElement(element), topology = createContourTopology(prepared.points, prepared);
    assert.deepEqual(prepared.trailingEdge, { kind: 'finite-base', upperIndex: 0, lowerIndex: points.length - 1 });
    assert.equal(topology.base.panels.length, 1);
    assert.deepEqual(prepared.sourcePoints, points);
    assert.deepEqual(topology.surface.points, open());
    assert.deepEqual(prepared.points.at(-1), prepared.points[0]);
    assert.deepEqual(prepareAirfoilElement(prepared), prepared);
    assert.deepEqual(element, before); assert.notEqual(prepared.sourcePoints, points);
    const text = serializeCoordinates({ name: element.name, elements: [prepared] }, { format: 'xfoil-labeled' });
    assert.deepEqual(parseCoordinates(text).elements[0].points, points);
  }
});

test('closed sharp preparation keeps the existing schema, coordinates and canonical orientation', () => {
  for (const points of [naca4('0012', 20), naca4('0012', 20).reverse()]) {
    const element = { name: 'Sharp', points, trips: [.1, .1] }, expected = { ...element, points: prepareContour(points) };
    assert.deepEqual(prepareAirfoilElement(element), expected);
    assert.equal(Object.hasOwn(prepareAirfoilElement(element), 'sourcePoints'), false);
    assert.equal(Object.hasOwn(prepareAirfoilElement(element), 'trailingEdge'), false);
  }
});

test('explicit measured NLR bases are retained; standard wetted export is explicit and JSON remains lossless', () => {
  const data = JSON.parse(fs.readFileSync('third_party/airfoils/nlr7301/geometry.json', 'utf8'));
  const elements = data.elements.map((e, i) => ({ name: e.name, points: e.points.map(([x, y]) => ({ x, y })),
    trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: i === 0 ? 424 : 216 } }));
  const prepared = elements.map(prepareAirfoilElement);
  assert.deepEqual(prepared, elements);
  prepared.forEach(e => assert.equal(createContourTopology(e.points, e).base.panels.length, 32));
  assert.deepEqual(JSON.parse(JSON.stringify({ elements: prepared })).elements, elements);
  const source = parseCoordinates(serializeCoordinates({ name: 'NLR source coordinates', elements: prepared, domain }, { format: 'mses' }));
  assert.deepEqual(source.elements.map(e => e.points), elements.map(e => e.points));
  const wetted = parseCoordinates(serializeCoordinates({ name: 'NLR wetted surfaces', elements: prepared, domain }, { format: 'mses', coordinates: 'wetted' }));
  assert.deepEqual(wetted.elements.map(e => e.points.length), [425, 217]);
  wetted.elements.forEach((e, i) => assert.deepEqual(e.points, elements[i].points.slice(0, elements[i].trailingEdge.lowerIndex + 1)));
  // A standard open file can imply only a straight base; it cannot carry the
  // measured 32-segment base and explicit corner identity. No hidden guess.
  wetted.elements.map(prepareAirfoilElement).forEach(e => assert.equal(createContourTopology(e.points, e).base.panels.length, 1));
});

test('included XFOIL and MSES coordinate files round-trip without numeric or frame changes', () => {
  const files = fs.readdirSync('third_party/Xfoil/runs').filter(f => f.endsWith('.dat')).map(f => `third_party/Xfoil/runs/${f}`);
  files.push('third_party/airfoils/rae2822-mses/blade.rae', 'public/examples/blade.nlr7301-wind');
  for (const file of files) {
    const parsed = parseCoordinates(fs.readFileSync(file, 'utf8')), replay = parseCoordinates(serializeCoordinates(parsed));
    assert.deepEqual(replay.elements, parsed.elements, file); assert.deepEqual(replay.domain, parsed.domain, file);
    assert.deepEqual(replay.header, parsed.header, file); assert.equal(replay.name, parsed.name, file);
  }
});

test('corner markers remain in parsed data but unsupported solver preparation rejects them explicitly', () => {
  const points = open(); points.splice(5, 0, { ...points[5] });
  const parsed = parseCoordinates('Cove markers\n' + rowText(points));
  assert.deepEqual(parsed.elements[0].points, points); assert.match(parsed.warnings.join(' '), /corner markers/);
  assert.deepEqual(parseCoordinates(serializeCoordinates(parsed)).elements[0].points, points);
  assert.throws(() => prepareAirfoilElement(parsed.elements[0]), /corner\/split-spline/);
});

test('export selection, signed zero and malformed inputs have explicit behavior', () => {
  const a = { name: 'A', points: open() }, b = { name: 'B', points: transform(open(), { x: 2 }) };
  const data = { name: 'Pair', elements: [a, b], domain };
  assert.throws(() => serializeCoordinates(data, { format: 'xfoil-plain' }), /one element/);
  assert.deepEqual(parseCoordinates(serializeCoordinates(data, { format: 'xfoil-plain', elementIndex: 1 })).elements[0].points, b.points);
  assert.throws(() => serializeCoordinates({ name: 'One', elements: [a] }, { format: 'mses' }), /domain bounds/);
  const signed = open(); signed[10].x = -0;
  const parsed = parseCoordinates(serializeCoordinates({ elements: [{ points: signed }] }, { format: 'xfoil-plain' }));
  assert.equal(Object.is(parsed.elements[0].points[10].x, -0), true);
  for (const text of ['', '# comment only', '999 999\n' + rowText(open()), 'Bad\n-2 3 5 4\n' + rowText(open()),
    'Name\n1 2 3', '00 12 NACA\n' + rowText(open()), rowText(open()) + '\n999 999\n999 999', rowText(open()) + '\nNaN 1'])
    assert.throws(() => parseCoordinates(text));
  assert.throws(() => prepareAirfoilElement({ ...a, trailingEdge: { kind: 'sharp' } }), /not snapped/);
});
