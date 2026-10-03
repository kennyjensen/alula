// Independent perfect-gas stagnation-pressure check from static p, rho, q.
// No reduced streamtube equations, enthalpy values or entropy-row imports.
export function stagnationPressureError(sections, { gamma, referencePressure, freestreamMach }) {
  const p0 = referencePressure * (1 + .5 * (gamma - 1) * freestreamMach ** 2) ** (gamma / (gamma - 1));
  let maxRelativeError = 0, maxStreamwiseChange = 0;
  const values = sections.map(row => row.map(s => {
    const machSquared = s.rho * s.q ** 2 / (gamma * s.p);
    const pressure = s.p * (1 + .5 * (gamma - 1) * machSquared) ** (gamma / (gamma - 1));
    maxRelativeError = Math.max(maxRelativeError, Math.abs(pressure / p0 - 1)); return pressure;
  }));
  for (let i = 1; i < values.length; i++) for (let j = 0; j < values[i].length; j++)
    maxStreamwiseChange = Math.max(maxStreamwiseChange, Math.abs((values[i][j] - values[i - 1][j]) / p0));
  return { maxRelativeError, maxStreamwiseChange };
}
