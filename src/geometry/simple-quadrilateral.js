// SPDX-License-Identifier: GPL-2.0-or-later
// Domain check for an intermediate finite-volume polygon. Concavity is
// permitted; crossed/touching opposite edges, reversals and zero area are not.
// This does not certify a bilinear map or a physically accepted final mesh.
import { requirePositiveSimplePolygon } from './simple-polygon.js';
export function requirePositiveSimpleQuad(vertices) {
  if (!Array.isArray(vertices) || vertices.length !== 4 || vertices.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))
    throw new Error('A simple quadrilateral requires four finite vertices.');
  return requirePositiveSimplePolygon(vertices, 'quadrilateral');
}
