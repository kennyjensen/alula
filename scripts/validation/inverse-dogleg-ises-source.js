// SPDX-License-Identifier: GPL-2.0-or-later
// Build a temporary research variant of the unchanged ISES update driver.
// Every geometry/transition operation stays in its original order.
export function inverseDoglegIsesSource(original) {
  const modifications = [
    ['onMesh, onCheckpoint, resume, dissipationEnhancement, iterationRecovery, maxProgressExtraIterations = 0 } = {})',
      'onMesh, onCheckpoint, resume, dissipationEnhancement, iterationRecovery, maxProgressExtraIterations = 0, researchModel, inversePrepare, trustRadius = resume?.inverseDogleg?.radius } = {})'],
    ['state.set(candidate.euler.adoptGeometry(state.subarray(0, candidate.ne), nodes));',
      `state.set(candidate.euler.adoptGeometry(state.subarray(0, candidate.ne), nodes));
    const inverse = inversePrepare(candidate, state);
    if (JSON.stringify(inverse.phase) !== JSON.stringify(candidate.bl.snapshotActive())) throw new Error('Grid maintenance changed the prepared inverse transition interval.');
    state.set(inverse.x);`],
    ['const checkpoint = () => ({ version: 1, families:',
      'const checkpoint = () => ({ version: 1, inverseDogleg: { radius: trustRadius }, families:'],
    ['let ordering, linear, stationPolicy;', 'let ordering, linear, stationPolicy, trustModel;'],
    ['if (!ordering) preferredOrdering = linear.ordering ?? preferredOrdering;',
      `trustModel = researchModel(system, state, value, { preferredOrdering, pivotTolerance });
      linear = trustModel.linear;
      if (trustRadius === undefined) trustRadius = trustModel.initialRadius;
      if (!ordering) preferredOrdering = linear.ordering ?? preferredOrdering;`],
    ['let maximumStep = recovery ? .5 : 1, redistributionBacktracks = 0, geometryRedistribution = false;',
      'let maximumStep = recovery ? .5 : 1, dogleg, redistributionBacktracks = 0, geometryRedistribution = false;'],
    ["stage = 'coupled density proposal'; proposal = propose(system, state, linear.x, { stagnationLimiter, maximumStep });",
      `stage = 'reduced dogleg proposal'; dogleg = trustModel.propose(trustRadius);
          stage = 'coupled density proposal'; proposal = propose(system, state, dogleg.direction, { stagnationLimiter, maximumStep });
          proposal.stepKind = 'coupled-inverse-' + dogleg.kind;`],
    ["} else {\n            stage = 'material-trip transfer'; event = events.prepare(proposal.x, state);",
      `} else {
            stage = 'inverse BL preparation'; const inverse = inversePrepare(system, proposal.x); proposal.x.set(inverse.x);
            stage = 'material-trip transfer'; event = events.prepare(proposal.x, state);
            if (JSON.stringify(inverse.phase) !== JSON.stringify(system.bl.snapshotActive())) throw new Error('Guarded inverse transition interval mismatch.');
            proposal.x.set(inverse.x);`],
    ['accepted = { candidate, event, inlet, dekink, triggeredBodies, r, redistributionAttempted,',
      `stage = 'complete trust-region reduction';
          proposal.trustRegion = trustModel.assess(dogleg.reducedDirection.map(v => proposal.step * v), candidate.value.residual, trustRadius);
          if (!proposal.trustRegion.accepted) throw new Error('Complete inverse-BL residual missed the trust-region acceptance gate.');
          trustRadius = proposal.trustRegion.nextRadius;
          accepted = { candidate, event, inlet, dekink, triggeredBodies, r, redistributionAttempted,`],
    ['rejections.push({ step, stage, message: error.message,',
      'rejections.push({ step, radius: trustRadius, kind: dogleg?.kind, trustRegion: proposal?.trustRegion, stage, message: error.message,'],
    ["maximumStep = (error.code === 'streamtube-grid-step' ? error.stepFraction : .5) * step;",
      'trustRadius *= .25;'],
    ['activeChange: true, changes: event.changes, meritComparable: false',
      'activeChange: true, changes: event.changes, meritComparable: true'],
    ["stepMethod: 'ises-density-newton'", "stepMethod: 'ises-inverse-dogleg'"],
  ];
  let code = original;
  for (const [a, b] of modifications) {
    if (code.split(a).length !== 2) throw new Error(`Ambiguous ISES research replacement: ${a}`);
    code = code.replace(a, b);
  }
  return { code, modifications };
}
