// SPDX-License-Identifier: GPL-2.0-or-later
// Numeric-only display rows. Never label an active temporary residual as final.
export function quadResidualProgressRows(iteration) {
  const rows = [], names = { euler: 'Euler equation error', boundaryLayer: 'BL equation error', edgeMatching: 'Edge matching error' };
  for (const [key, label] of Object.entries(names))
    if (Number.isFinite(iteration[key])) rows.push({ label, value: iteration[key].toExponential(2) });
  if (Number.isFinite(iteration.residualContext?.mcrit)) rows.push({ label: 'Active MCRIT', value: iteration.residualContext.mcrit.toFixed(4) });
  if (iteration.prescribedResidual) {
    for (const [key, label] of Object.entries(names))
      if (Number.isFinite(iteration.prescribedResidual[key])) rows.push({ label: `Prescribed ${label.toLowerCase()}`, value: iteration.prescribedResidual[key].toExponential(2) });
    if (iteration.prescribedResidual.unavailable) rows.push({ label: 'Prescribed-equation check', value: 'State inadmissible under prescribed equations' });
  } else if (iteration.residualContext) rows.push({ label: 'Prescribed-equation check', value: 'Not sampled this iteration' });
  if (Number.isFinite(iteration.progress?.stateChange)) rows.push({ label: 'Relative state change', value: iteration.progress.stateChange.toExponential(2) });
  return rows;
}
