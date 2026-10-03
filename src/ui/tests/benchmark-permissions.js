// SPDX-License-Identifier: GPL-2.0-or-later
// Compatibility API for external inventories. Preset metadata never locks a
// solver: all airfoils may attempt all modes. Numerical/input validation and
// convergence checks remain in the solvers and their result handling.
export function benchmarkRestriction() { return ''; }
