// SPDX-License-Identifier: GPL-2.0-or-later
// Shared geometry definitions. Clone before modifying a case.
export const builtinAirfoils = {
  single: [{ name: 'Main element', code: '0012', chord: 1, x: 0, y: 0, deflection: 0 }],
  flap: [{ name: 'Main element', code: '2412', chord: 1, x: 0, y: 0, deflection: 0 },
    { name: 'Flap', code: '0012', chord: 0.3, x: 0.94, y: -0.08, deflection: 15 }],
  three: [{ name: 'Slat', code: '0012', chord: 0.2, x: -0.22, y: 0.05, deflection: -20 },
    { name: 'Main element', code: '2412', chord: 1, x: 0, y: 0, deflection: 0 },
    { name: 'Flap', code: '0012', chord: 0.3, x: 0.94, y: -0.08, deflection: 15 }],
};
