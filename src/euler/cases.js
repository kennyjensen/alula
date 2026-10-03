// SPDX-License-Identifier: GPL-2.0-or-later
import { naca4, transform } from '../geometry/airfoil.js';
import { channelMesh } from './mesh.js';
import { multielementMesh } from './multielement-mesh.js';

export function buildEulerLabCase(input = {}) {
  const controls = { configuration: 'two', panels: 20, rows: 2, padding: 1, mach: .3, alpha: 0, maxIterations: 30, ...input };
  const { configuration, panels, rows, padding, mach, alpha, maxIterations } = controls;
  if (!['single', 'two', 'three', 'channel'].includes(configuration) || ![20, 40].includes(panels) || ![2, 3, 4].includes(rows)
    || ![1, 2, 4].includes(padding) || !Number.isFinite(mach) || mach < .15 || mach > .5 || !Number.isFinite(alpha) || Math.abs(alpha) > 6
    || !Number.isInteger(maxIterations) || maxIterations < 1 || maxIterations > 40) throw new Error('Invalid Euler laboratory controls.');
  const moving = configuration === 'channel';
  if (moving && alpha !== 0) throw new Error('The channel verification case requires zero inlet angle.');
  const main = { name: 'Main · NACA 0012', points: naca4('0012', panels) };
  const flap = { name: 'Flap · c = 0.4', points: transform(naca4('0012', panels), { chord: .4, x: 1.05, y: -.2, angle: -10 }) };
  const slat = { name: 'Slat · c = 0.25', points: transform(naca4('0012', panels), { chord: .25, x: -.4, y: .22 }) };
  const elements = moving ? [] : configuration === 'single' ? [main] : configuration === 'two' ? [main, flap] : [slat, main, flap];
  const mesh = moving ? channelMesh({ nx: panels / 2, ny: rows, upper: x => 1 - .08 * Math.sin(Math.PI * x / 2) ** 2 })
    : multielementMesh(elements.map(e => e.points), { rows, padding });
  if (mesh.cells.length > 600) throw new Error(`This grid has ${mesh.cells.length} cells; the dense reference limit is 600. Reduce panels or rows.`);
  return { controls, moving, elements, mesh, options: { mach, alpha, maxIterations },
    title: moving ? 'Moving streamlines' : configuration === 'single' ? 'Single element' : configuration === 'two' ? 'Main element + flap' : 'Slat + main + flap' };
}
