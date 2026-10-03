// SPDX-License-Identifier: GPL-2.0-or-later
// Spacing references for a shared finite-base wake use the center between
// its two banks. The actual solid corner samples and branch potentials stay
// in their separate guide arrays; this is only the cut's geometric metric.
export function streamtubeCutEdge(profile, body, end) {
  const incoming = end === 'upstream';
  const upper = profile.curve.branch('upper', incoming ? 0 : 1, profile.stag).point;
  if (incoming || body.trailingEdge?.kind !== 'finite-base') return upper;
  const lower = profile.curve.branch('lower', 1, profile.stag).point;
  return { x: .5 * (upper.x + lower.x), y: .5 * (upper.y + lower.y) };
}

export function streamtubeCutReferencePoints(body, guides, end, first, last) {
  const lower = guides.lower.slice(first, last + 1);
  if (end !== 'wake' || body.trailingEdge?.kind !== 'finite-base') return lower;
  return lower.map((p, k) => ({ ...p, x: .5 * (p.x + guides.upper[first + k].x),
    y: .5 * (p.y + guides.upper[first + k].y) }));
}
