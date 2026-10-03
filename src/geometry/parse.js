// SPDX-License-Identifier: GPL-2.0-or-later
// Plain/labeled XFOIL and ISES/MSES blade files. Parsing retains the supplied
// coordinate order, open endpoints and common frame. Solver preparation is
// deliberately separate: see prepareAirfoilElement.
export function parseCoordinates(text) {
  if (typeof text !== 'string') throw new Error('Coordinate file contents must be text.');
  const lines = text.replace(/^\uFEFF/, '').split(/\r\n|\n|\r/).map(line => line.trim()).filter(line => line && !/^[#!]/.test(line));
  if (!lines.length) throw new Error('The coordinate file is empty.');
  const token = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:[dDeE][+-]?\d+)?$/;
  const parts = line => line.split(/[#!]/)[0].trim().split(/[\s,]+/);
  const numeric = value => token.test(value) ? Number(value.replace(/[dD]/, 'e')) : NaN;
  const numbers = line => parts(line).map(numeric);
  let name = 'Imported airfoil';
  const first = numbers(lines[0]), labeled = !(first.length >= 2 && first.slice(0, 2).every(Number.isFinite));
  if (labeled) name = lines.shift();
  let domain = null, header = null;
  const warnings = [];
  if (lines.length && numbers(lines[0]).length >= 4) {
    const values = numbers(lines.shift());
    if (!values.every(Number.isFinite)) throw new Error('Invalid ISES/MSES grid header.');
    if (values.length === 4) {
      const [xMin, xMax, yMin, yMax] = values;
      if (xMin >= xMax || yMin >= yMax) throw new Error('Invalid MSES domain bounds.');
      domain = { xMin, xMax, yMin, yMax };
      header = { kind: 'mses-domain', values };
    } else {
      // The old MSES format uses different domain semantics. Preserve the
      // complete header rather than treating its first four values as bounds.
      header = { kind: 'legacy-ises-mses', values };
      warnings.push('Legacy ISES/MSES header values are preserved; their domain semantics are not converted to current MSES bounds.');
    }
  }
  const elements = []; let points = [], separated = false;
  const flush = () => {
    if (!points.length) throw new Error('An MSES separator must follow a nonempty element.');
    if (points.length < 9) throw new Error('Each imported element needs at least 9 coordinates. Lednicer count headers are not supported.');
    elements.push({ name: `${name}${elements.length ? ` ${elements.length + 1}` : ''}`, points });
    points = [];
  };
  for (const line of lines) {
    const pair = numbers(line);
    if (pair.length !== 2 || !pair.every(Number.isFinite)) throw new Error(`Invalid coordinate row: ${line.slice(0, 80)}`);
    if (pair[0] === 999 && pair[1] === 999) { flush(); separated = true; continue; }
    points.push({ x: pair[0], y: pair[1] });
  }
  if (points.length) flush();
  if (!elements.length) throw new Error('No airfoil coordinates found.');
  if (elements.some(e => e.points.some((p, i) => i > 0 && p.x === e.points[i - 1].x && p.y === e.points[i - 1].y)))
    warnings.push('Successive duplicate coordinates are retained as file corner markers; current solver preparation rejects unsupported corner/split-spline conditions.');
  return { name, elements, domain, format: header || separated ? 'mses' : labeled ? 'xfoil-labeled' : 'xfoil-plain',
    header, labeled, warnings };
}
