// SPDX-License-Identifier: GPL-2.0-or-later
// Standard coordinate text stores geometry, not solver/checkpoint metadata.
// "source" retains input order/points; "wetted" explicitly exports only the
// TE-to-TE surface of a finite-base solid. Use JSON to retain a measured base
// polyline together with its corner indices and all other metadata.
import { prepareAirfoilElement } from './airfoil-element.js';
import { createContourTopology } from './contour-topology.js';

export function serializeCoordinates(data, { format = 'source', coordinates = 'source', elementIndex } = {}) {
  if (!data || !Array.isArray(data.elements) || !data.elements.length)
    throw new Error('Coordinate export requires at least one element.');
  if (!['source', 'xfoil-plain', 'xfoil-labeled', 'mses'].includes(format)
    || !['source', 'wetted'].includes(coordinates)) throw new Error('Unknown coordinate export format or point selection.');
  if (elementIndex !== undefined && (!Number.isInteger(elementIndex) || elementIndex < 0 || elementIndex >= data.elements.length))
    throw new Error('Invalid coordinate export element index.');
  const elements = elementIndex === undefined ? data.elements : [data.elements[elementIndex]];
  const sourceFormat = format === 'source';
  if (sourceFormat) format = data.format ?? (elements.length > 1 || data.domain || data.header ? 'mses' : 'xfoil-labeled');
  if (!['xfoil-plain', 'xfoil-labeled', 'mses'].includes(format)) throw new Error('Unknown original coordinate format.');
  if (format !== 'mses' && elements.length !== 1) throw new Error('XFOIL coordinate files contain one element; select an element or export MSES.');
  if (format === 'mses' && !sourceFormat && !data.domain && !data.header)
    throw new Error('MSES export requires explicit domain bounds or a preserved ISES/MSES header; source export can retain a headerless input.');
  const number = value => {
    if (!Number.isFinite(value)) throw new Error('Coordinate export requires finite numbers.');
    return Object.is(value, -0) ? '-0' : String(value);
  };
  const name = data.name ?? (elements.length === 1 ? elements[0].name : 'Airfoil assembly') ?? 'Imported airfoil';
  if (typeof name !== 'string' || !name.trim() || /[\r\n]/.test(name)) throw new Error('Coordinate name must occupy one nonempty line.');
  const rows = [];
  if (format === 'xfoil-labeled' || format === 'mses' && (!sourceFormat || data.labeled !== false)) {
    const first = name.trim().split(/[\s,]+/).slice(0, 2);
    if (/^[#!]/.test(name.trim()) || first.length === 2 && first.every(v => /^[+-]?(?:\d+\.?\d*|\.\d+)(?:[dDeE][+-]?\d+)?$/.test(v)))
      throw new Error('An airfoil label must not look like a coordinate pair or comment.');
    rows.push(name);
  }
  if (format === 'mses') {
    if (data.header?.kind === 'legacy-ises-mses') {
      if (!Array.isArray(data.header.values) || data.header.values.length < 5) throw new Error('Invalid legacy ISES/MSES header.');
      rows.push(data.header.values.map(number).join(' '));
    } else if (data.domain) {
      const { xMin, xMax, yMin, yMax } = data.domain;
      if (!(xMin < xMax && yMin < yMax)) throw new Error('Invalid MSES domain bounds.');
      rows.push([xMin, xMax, yMin, yMax].map(number).join(' '));
    } else if (data.header) throw new Error('MSES current-domain export requires its complete domain bounds.');
  }
  elements.forEach((element, index) => {
    let points = element.sourcePoints ?? element.points;
    if (coordinates === 'wetted') {
      const prepared = prepareAirfoilElement(element);
      points = prepared.trailingEdge?.kind === 'finite-base'
        ? createContourTopology(prepared.points, prepared).surface.points : prepared.points;
    }
    if (!Array.isArray(points) || points.length < 9) throw new Error('Each exported element needs at least nine coordinates.');
    if (index) rows.push('999.0 999.0');
    for (const p of points) {
      if (!p || p.x === 999 && p.y === 999) throw new Error('The 999/999 pair is reserved for MSES element separation.');
      rows.push(`${number(p.x)} ${number(p.y)}`);
    }
  });
  return rows.join('\n') + '\n';
}
