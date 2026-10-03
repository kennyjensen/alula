// SPDX-License-Identifier: GPL-2.0-or-later
import { auditStreamtubeShocks } from '../streamtube-shock-audit.js';
// A mesh-continuation trigger, not a governing-equation switch or shock certificate.
export function shockRefinementRequest(result) {
  if (!result?.converged || !result.mesh?.quality?.valid || !(result.flow?.diagnostics?.maxMach > 1)) return null;
  const audit = auditStreamtubeShocks(result.flow);
  const candidate = audit.candidates.find(c => c.sonicCrossing && c.pressureRise >= .01);
  return candidate ? { reason: 'Resolved supersonic-to-subsonic compression; refine before further continuation',
    trigger: 'sonic-compression', candidate, maximumSkewDegrees: audit.maximumSkewDegrees } : null;
}
