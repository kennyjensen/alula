// SPDX-License-Identifier: GPL-2.0-or-later
// Independent ODE oracle. f''' + f f''/2 = 0, f(0)=f'(0)=0, f'(infinity)=1.
// Fourth-order Runge–Kutta plus scalar shooting; not a BL closure correlation.
export function blasius(step = 0.01, etaMax = 12) {
  const rhs = ([f, u, shear]) => [u, shear, -f * shear / 2, u * (1 - u)];
  function integrate(shear) {
    let state = [0, 0, shear, 0];
    for (let i = 0; i < Math.round(etaMax / step); i++) {
      const k1 = rhs(state);
      const k2 = rhs(state.map((v, j) => v + step * k1[j] / 2));
      const k3 = rhs(state.map((v, j) => v + step * k2[j] / 2));
      const k4 = rhs(state.map((v, j) => v + step * k3[j]));
      state = state.map((v, j) => v + step * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6);
    }
    return state;
  }
  let low = 0.3; let high = 0.4;
  for (let i = 0; i < 45; i++) { const mid = (low + high) / 2; if (integrate(mid)[1] > 1) high = mid; else low = mid; }
  const shear = (low + high) / 2; const solution = integrate(shear);
  return { theta: solution[3], deltaStar: etaMax - solution[0], h: (etaMax - solution[0]) / solution[3], cf: 2 * shear };
}
