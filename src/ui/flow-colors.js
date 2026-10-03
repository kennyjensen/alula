// SPDX-License-Identifier: GPL-2.0-or-later
// A fixed scale keeps an unchanged color tied to the same speed throughout
// the solve. Values above the legend's upper end saturate, not rescale.
export const flowSpeedMaximum = 3;
const stops = [[94, 147, 250], [74, 219, 199], [244, 209, 111], [245, 128, 93]];
export function flowSpeedColor(speed) {
  const value = Math.max(0, Math.min(flowSpeedMaximum, speed));
  const i = Math.min(stops.length - 2, Math.floor(value)), t = value - i;
  return `rgb(${stops[i].map((a, k) => Math.round(a + t * (stops[i + 1][k] - a))).join(',')})`;
}
export const flowSpeedGradient = `linear-gradient(90deg, ${stops.map((_, i) => flowSpeedColor(i)).join(', ')})`;
export const flowChangeMaximum = .1;
export function flowChangeColor(delta) {
  const center = [140, 170, 181], end = delta < 0 ? [112, 162, 255] : [255, 175, 96];
  const t = Math.min(1, Math.abs(delta) / flowChangeMaximum);
  return `rgb(${center.map((a, k) => Math.round(a + t * (end[k] - a))).join(',')})`;
}
export const flowChangeGradient = `linear-gradient(90deg, ${[-flowChangeMaximum, 0, flowChangeMaximum].map(flowChangeColor).join(', ')})`;
