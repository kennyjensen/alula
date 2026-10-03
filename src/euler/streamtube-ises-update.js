// SPDX-License-Identifier: GPL-2.0-or-later
// Research integration of additive-density Newton with ISES grid maintenance.
// Optional common-step backtracking checks the grid/gas domain and, with
// Armijo selection, residual decrease after the complete grid maintenance. Fresh charts validate complete updates before commitment.
import { tangentialNewton } from './streamtube-tangential-newton.js';
import { msesTemporaryMcrit } from './streamtube-shock-audit.js';
import { requireCoupledResidualDecrease as requireResidualDecrease, repeatedPassageStepLimit, repeatedGridStepLimit } from './streamtube-iteration-progress.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from './streamtube-body.js';
import { proposeDensityNewton } from './streamtube-density-newton.js';
import { solveSparseDirectAligned } from '../numerics/klu.js';
import { streamtubeEquationOrder } from './streamtube-body-layout.js';
import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';
import { captureStreamtubeInletFractions, adjustStreamtubeInlets, reparameterizeStreamtubeInlets, dekinkStreamtubeInteriors } from '../geometry/streamtube-grid-maintenance.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { limitStreamtubeGridStep, assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';
import { requireConvexGridUpdate, requireConvexPublishedGrid } from './streamtube-grid-update.js';

export function solveStreamtubeIses(input, { initialEuler, maxIterations = 12, tolerance = 1e-10, stagnationLimiter = 'listing', iterationGeometry = 'convex', stepAcceptance = 'listing', maxBacktracks = 12, onIteration, onMesh, onCheckpoint, retainCheckpoint = false, retainBestCheckpoint = false, traceUpdates = false, tangentialGridRecovery = false, initialRedistributionMode = 'smove', adaptiveMcrit = false, targetMcrit, firstOrderStartup = false, targetMucon, broadShockStartup = false, broadShockMucon = 2, stopOnGridStagnation = false, gridCorrectionBacktracking = 'after-search', resume } = {}) {
  if (resume !== undefined) {
    if (input !== undefined || initialEuler !== undefined || resume?.version !== 1
      || !resume.input || !resume.initialEuler || !resume.continuation)
      throw new Error('Resume requires one complete Euler ISES checkpoint and no separate initial state.');
    const continuation = resume.continuation;
    if (tangentialGridRecovery !== (continuation.tangentialGridRecovery ?? false)
      || initialRedistributionMode !== (continuation.initialRedistributionMode ?? 'smove')
      || iterationGeometry !== continuation.iterationGeometry || stepAcceptance !== continuation.stepAcceptance
      || stagnationLimiter !== continuation.stagnationLimiter
      || adaptiveMcrit !== (continuation.adaptiveMcrit ?? false)
      || targetMcrit !== continuation.targetMcrit
      || firstOrderStartup !== (continuation.firstOrderStartup ?? false)
      || targetMucon !== continuation.targetMucon
      || broadShockStartup !== (continuation.broadShockStartup ?? false)
      || broadShockMucon !== (continuation.broadShockMucon ?? 2)
      || stopOnGridStagnation !== (continuation.stopOnGridStagnation ?? false)
      || gridCorrectionBacktracking !== (continuation.gridCorrectionBacktracking ?? 'after-search')) throw new Error('Euler ISES checkpoint update controls do not match.');
    resume = structuredClone(resume);
    ({ input, initialEuler } = resume);
  }
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0 || !['listing', 'prose'].includes(stagnationLimiter)
    || !['convex', 'ises-sampled'].includes(iterationGeometry) || !['listing', 'admissible', 'armijo'].includes(stepAcceptance)
    || !Number.isInteger(maxBacktracks) || maxBacktracks < 0 || maxBacktracks > 20
    || tangentialGridRecovery && stepAcceptance !== 'armijo'
    || typeof tangentialGridRecovery !== 'boolean' || !['smove', 'supplied'].includes(initialRedistributionMode)
    || !resume && initialRedistributionMode === 'supplied' && !initialEuler
    || typeof traceUpdates !== 'boolean' || typeof adaptiveMcrit !== 'boolean' || typeof firstOrderStartup !== 'boolean'
    || !adaptiveMcrit && targetMcrit !== undefined
    || !firstOrderStartup && targetMucon !== undefined
    || typeof stopOnGridStagnation !== 'boolean' || !Number.isFinite(broadShockMucon) || broadShockMucon < 2 || broadShockMucon > 8
    || typeof broadShockStartup !== 'boolean' || broadShockStartup && (!firstOrderStartup || !adaptiveMcrit)
    || !['after-search', 'before-damping'].includes(gridCorrectionBacktracking)
    || typeof retainCheckpoint !== 'boolean' || typeof retainBestCheckpoint !== 'boolean'
    || retainBestCheckpoint && stepAcceptance === 'listing'
    || (onCheckpoint !== undefined && typeof onCheckpoint !== 'function'))
    throw new Error('Invalid ISES research iteration controls.');
  if (adaptiveMcrit) {
    if (!input.upwind) throw new Error('Adaptive MCRIT requires speed upwinding.');
    targetMcrit ??= input.upwind.mcrit;
    msesTemporaryMcrit(targetMcrit, 0);
  }
  if (firstOrderStartup) {
    targetMucon ??= input.upwind?.mucon;
    if (!Number.isFinite(targetMucon) || targetMucon <= 0
      || !input.upwind || (resume && broadShockStartup
        ? ![targetMucon, -Math.max(broadShockMucon, targetMucon)].includes(input.upwind.mucon)
        : (resume ? Math.abs(input.upwind.mucon) : input.upwind.mucon) !== targetMucon))
      throw new Error('First-order Euler startup requires a positive requested MUCON and a matching checkpoint law.');
    // Start on the same grid and gas state. The sign disables only the
    // second-order correction; it does not change ISMOM or the operating point.
    if (!resume) input = { ...input, upwind: { ...input.upwind,
      mucon: broadShockStartup ? -Math.max(broadShockMucon, targetMucon) : -targetMucon,
      ...(broadShockStartup ? { mcrit: .75 } : {}) } };
  }
  if (input.stagnationMotion !== undefined && !['walls-only', 'interpolated'].includes(input.stagnationMotion)) throw new Error('Unknown ISES stagnation motion.');
  const geometryDomain = iterationGeometry === 'convex' ? 'convex' : 'positive-simple';
  // Sampled equation geometry is useful for derivatives, but does not grant
  // permission to accept a concave mesh. All public solves use this policy.
  const preserveConvexity = stepAcceptance !== 'listing';
  if (input.geometryDomain !== undefined && input.geometryDomain !== geometryDomain) throw new Error('Conflicting ISES iteration geometry domains.');
  input = { ...input, stagnationMotion: input.stagnationMotion ?? 'walls-only', normalStencil: input.normalStencil ?? 'body-stations', geometryDomain };
  let system = createStreamtubeBodySystem(input);
  // Governing pressure rows use pInf to condition the Newton matrix. At low
  // Mach that scaling hides pressure imbalance from the line-search merit
  // (pInf is O(1/M²)), while coordinate rows remain O(1). Compare pressure
  // in the same dynamic units as the velocities. Keep the governing rows,
  // Jacobian and original convergence tolerance unchanged.
  const pressureRows = new Set(['streamwise', 'internalPressure', 'farfieldPressure',
    'cutPressure', 'trailingKutta', 'leadingKutta']);
  const meritWeights = Float64Array.from(system.layout.rows,
    row => pressureRows.has(row.kind) ? system.conditions.pressureScale : 1);
  const requireCoupledResidualDecrease = (current, candidate, controls) =>
    requireResidualDecrease(current, candidate, { ...controls, weights: meritWeights });
  if (!system.layout.densityCount) throw new Error('ISES density update requires compressible equations.');
  let state = initialEuler ? system.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes) : system.initial.slice();
  let value = system.evaluate(state), reason = 'iteration limit', failed = false, lastRejectedStep = null;
  if (preserveConvexity) {
    assertConvexStreamtubeGrid(value.nodes);
    requireConvexPublishedGrid(system, value.nodes);
  }
  if (resume) {
    if ((!Array.isArray(resume.residual) && !ArrayBuffer.isView(resume.residual))
      || resume.residual.length !== value.residual.length
      || value.residual.some((v, k) => v !== resume.residual[k]))
      throw new Error('Euler ISES checkpoint residual does not replay exactly.');
    const expected = initialEuler.nodes;
    if (!Array.isArray(expected) || expected.length !== value.nodes.length
      || value.nodes.some((group, g) => expected[g]?.length !== group.length
        || group.some((line, i) => expected[g][i]?.length !== line.length
          || line.some((p, j) => p.x !== expected[g][i][j]?.x || p.y !== expected[g][i][j]?.y))))
      throw new Error('Euler ISES checkpoint physical nodes do not replay exactly.');
  }
  let fractions = resume ? resume.continuation.fractions : captureStreamtubeInletFractions(value.nodes, system.layout.bodies);
  if (resume && (!Array.isArray(fractions) || fractions.length !== input.bodies.length || fractions.some((f, body) =>
    !Array.isArray(f) || f.length !== input.bodies[body].leadingIndex + 1 || f[0] !== 0 || f.at(-1) !== 1
    || !f.every(Number.isFinite) || f.some((v, i) => i && v <= f[i - 1]))))
    throw new Error('Invalid Euler ISES checkpoint inlet fractions.');
  if (resume && (!Array.isArray(resume.continuation.lastRedistributedStagnation)
    || resume.continuation.lastRedistributedStagnation.length !== input.bodies.length
    || !resume.continuation.lastRedistributedStagnation.every((s, b) => Number.isFinite(s) && s > 0 && s < system.curves[b].length)))
    throw new Error('Invalid Euler ISES checkpoint redistribution history.');
  const adopt = (x, nodes) => {
    if (preserveConvexity) assertConvexStreamtubeGrid(nodes);
    const candidate = createStreamtubeBodySystem(input), state = candidate.adoptGeometry(x, nodes);
    const value = candidate.evaluate(state);
    if (preserveConvexity) {
      assertConvexStreamtubeGrid(value.nodes);
      requireConvexPublishedGrid(candidate, value.nodes);
    }
    return { system: candidate, state, value };
  };
  const changeUpwind = changes => {
    if (Object.entries(changes).every(([key, value]) => value === input.upwind[key])) return;
    const previousInput = input;
    input = { ...input, upwind: { ...input.upwind, ...changes } };
    try { ({ system, state, value } = adopt(state, value.nodes)); }
    catch (error) { input = previousInput; throw error; }
  };
  const requestedDissipation = () => (!adaptiveMcrit || input.upwind.mcrit === targetMcrit)
    && (!firstOrderStartup || input.upwind.mucon === targetMucon);
  const restoreRootDissipation = () => {
    // Converge a broad first-order shock at fixed MCRIT, then switch both
    // parameters together. Only the restored equations can certify a result.
    if (broadShockStartup && input.upwind.mucon < 0) {
      if (value.diagnostics.residual <= tolerance) changeUpwind({ mucon: targetMucon, mcrit: targetMcrit });
      return;
    }
    if (adaptiveMcrit && value.diagnostics.residual <= tolerance) changeUpwind({ mcrit: targetMcrit });
    // First converge the broadened shock, then solve the requested law.
    // Both phases share this Newton loop, its budget and checkpoint history.
    if (firstOrderStartup && value.diagnostics.residual <= tolerance) changeUpwind({ mucon: targetMucon });
  };
  // Optional read-only diagnosis. Every stage is evaluated on its own chart;
  // the actual proposal, residual checks and accepted state are untouched.
  const updateTrace = [];
  let activeTrace;
  const traceStage = (name, x, nodes) => {
    if (!activeTrace) return;
    let maximumNodeMotion = 0;
    nodes.forEach((group,g)=>group.forEach((row,i)=>row.forEach((p,j)=>{
      const q=value.nodes[g][i][j];
      maximumNodeMotion=Math.max(maximumNodeMotion,Math.hypot(p.x-q.x,p.y-q.y));
    })));
    try {
      const candidate=adopt(x,nodes), r=candidate.value.residual;
      let squaredNorm=0,maximum=0,worstRow=-1;
      for(let k=0;k<r.length;k++){squaredNorm+=r[k]*r[k];if(Math.abs(r[k])>maximum){maximum=Math.abs(r[k]);worstRow=k;}}
      activeTrace.stages.push({name,admissible:true,squaredNorm,maximumResidual:maximum,
        worstRow,row:structuredClone(candidate.system.layout.rows[worstRow]),
        residualByFamily:{...candidate.value.diagnostics.residualByFamily},maximumNodeMotion});
    } catch(error) {
      activeTrace.stages.push({name,admissible:false,maximumNodeMotion,message:error.message,
        code:error.code,diagnostics:structuredClone(error.diagnostics)});
    }
  };
  const redistribute = (nodes, correctionScale = 1) => {
    const passages = [], moved = nodes.map((group, g) => {
      const r = redistributeStreamtubeTangentially(group, { referenceBank: g === 0 ? group[0].length - 1 : 0,
        fixedBanks: [g !== 0, g !== nodes.length - 1], correctionScale, quadratureDomain: iterationGeometry === 'convex' ? 'convex' : 'sampled-positive' });
      passages.push({ group: g, referenceBank: r.referenceBank, fixedBanks: r.fixedBanks, pairs: r.solution.pairs,
        coordinateRelativeResidual: r.solution.relativeResidual, maxDisplacement: r.maxDisplacement,
        ...(correctionScale === 1 ? {} : { correctionScale }) });
      return r.nodes;
    });
    return { nodes: moved, passages };
  };
  const initialRedistribution = { beforeResidual: value.diagnostics.residual, accepted: false };
  try {
    if (resume) Object.assign(initialRedistribution, { accepted: true, resumed: true, afterResidual: value.diagnostics.residual, passages: [] });
    else if (initialRedistributionMode === 'supplied') {
      Object.assign(initialRedistribution, { accepted: true, supplied: true, afterResidual: value.diagnostics.residual, passages: [] });
    } else {
      // The source SMOVE map is tried unchanged first. Under the existing
      // admissible policy, globalize only its coordinate correction, keeping
      // the same original state and one common scale for every passage.
      // This is a domain safeguard, not a residual-descent condition.
      let correctionScale = 1;
      const rejections = [];
      for (let trial = 0; ; trial++) {
        let stage = 'SMOVE';
        try {
          const r = redistribute(value.nodes, correctionScale);
          stage = 'admissibility';
          const candidate = adopt(state, r.nodes);
          ({ system, state, value } = candidate);
          Object.assign(initialRedistribution, { accepted: true, afterResidual: value.diagnostics.residual, passages: r.passages,
            ...(rejections.length ? { correctionScale, backtracks: rejections.length, rejections } : {}) });
          break;
        } catch (error) {
          if (stepAcceptance === 'listing') throw error;
          rejections.push({ correctionScale, stage, message: error.message });
          Object.assign(initialRedistribution, { correctionScale, backtracks: rejections.length - 1, rejections });
          if (trial >= maxBacktracks) throw error;
          correctionScale *= .5;
        }
      }
    }
  } catch (error) {
    reason = `ISES initial redistribution rejected: ${error.message}`; failed = true; initialRedistribution.rejection = error.message;
    // This is a rejected initialization movement, before any Newton update.
    // Use the existing failure envelope so the precursor/Worker can retain
    // the original local pressure or geometry diagnosis without extra flow data.
    lastRejectedStep = { stage: 'initial redistribution', message: error.message,
      ...(error.code === undefined ? {} : { code: error.code }),
      ...(error.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }) };
  }
  let lastRedistributedStagnation = resume ? resume.continuation.lastRedistributedStagnation.slice() : value.stagnation.slice();
  const history = [{ iteration: 0, residual: value.diagnostics.residual, step: 0,
    ...(firstOrderStartup ? { dissipation: { mucon: input.upwind.mucon, targetMucon,
      mcrit: input.upwind.mcrit, targetMcrit: targetMcrit ?? input.upwind.mcrit } } : {}) }];
  const equationOrder = streamtubeEquationOrder(system.layout);
  const linearDiagnostics = { solves: 0, maxRelativeResidual: 0, refinements: 0, maxFactorNonzeros: 0, orderingFallbacks: 0 };
  const checkpoint = () => structuredClone({ version: 1, input,
    initialEuler: { x: Array.from(state), nodes: value.nodes }, residual: Array.from(value.residual),
    continuation: { ...(tangentialGridRecovery ? { tangentialGridRecovery } : {}),
      ...(initialRedistributionMode === 'supplied' ? { initialRedistributionMode } : {}),
      fractions, lastRedistributedStagnation, stagnationLimiter, iterationGeometry, stepAcceptance,
      ...(adaptiveMcrit ? { adaptiveMcrit, targetMcrit } : {}),
      ...(broadShockStartup ? { broadShockStartup, broadShockMucon } : {}),
      ...(stopOnGridStagnation ? { stopOnGridStagnation } : {}),
      ...(gridCorrectionBacktracking !== 'after-search' ? { gridCorrectionBacktracking } : {}),
      ...(firstOrderStartup ? { firstOrderStartup, targetMucon } : {}) } });
  let bestCheckpoint;
  const emitCheckpoint = () => {
    // Optional startup evidence only. Preserve the terminal Euler result and
    // its honest convergence status. Copy a complete admitted state only on
    // improvement; a residual vector alone cannot restore its geometry chart.
    if (retainBestCheckpoint && initialRedistribution.accepted
      && Number.isFinite(value.diagnostics.residual)
      && (!bestCheckpoint || value.diagnostics.residual < bestCheckpoint.residual))
      bestCheckpoint = { iteration: history.length - 1, residual: value.diagnostics.residual, checkpoint: checkpoint() };
    // User cancellation must propagate, not become a rejected Newton step.
    // Both state and details are detached from the active solver.
    if (initialRedistribution.accepted && onCheckpoint)
      onCheckpoint(checkpoint(), structuredClone({ history, linearDiagnostics, initialRedistribution }));
  };
  onIteration?.(history[0]); onMesh?.({ system, nodes: value.nodes, flow: value, iteration: { ...history[0] } });
  emitCheckpoint();
  const rejectUpdate = (error, stage, proposal) => {
    const { x, ...details } = proposal ?? {};
    lastRejectedStep = { ...details, stage, message: error.message,
      ...(error.code === undefined ? {} : { code: error.code }),
      ...(error.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }),
      ...(error.relativeResidual === undefined ? {} : { linearRelativeResidual: error.relativeResidual }),
      ...(error.attempts ? { linearAttempts: structuredClone(error.attempts) } : {}) };
    reason = `ISES research update rejected during ${stage}: ${error.message}`; failed = true;
  };
  for (let iteration = 1; !failed && iteration <= maxIterations
    && (value.diagnostics.residual > tolerance || !requestedDissipation()); iteration++) {
    let proposal, stage = 'linear solve';
    // A checkpoint may itself be a temporary-law root. Publish the restored
    // operator even when no additional Newton step is needed. Cancellation
    // from the checkpoint observer must propagate outside rejection handling.
    if (!requestedDissipation() && value.diagnostics.residual <= tolerance) {
      try { restoreRootDissipation(); }
      catch (error) { rejectUpdate(error, 'dissipation update'); break; }
      const entry = history.at(-1);
      entry.residual = value.diagnostics.residual;
      entry.dissipation = { ...entry.dissipation, mcrit: input.upwind.mcrit,
        targetMcrit: targetMcrit ?? input.upwind.mcrit,
        ...(firstOrderStartup ? { mucon: input.upwind.mucon, targetMucon } : {}) };
      onIteration?.(entry); onMesh?.({ system, nodes: value.nodes, flow: value, iteration: { ...entry } });
      emitCheckpoint();
      if (value.diagnostics.residual <= tolerance) break;
    }
    try {
      const jacobian = system.jacobian(state, { sparse: true });
      const linear = solveSparseDirectAligned(jacobian, value.residual.map(v => -v), equationOrder);
      const recordLinear = result => {
        linearDiagnostics.maxFactorNonzeros = Math.max(linearDiagnostics.maxFactorNonzeros, result.factorNonzeros);
        linearDiagnostics.orderingFallbacks += result.equationOrdering === 'original' ? 1 : 0;
        linearDiagnostics.solves++;
        linearDiagnostics.maxRelativeResidual = Math.max(linearDiagnostics.maxRelativeResidual, result.relativeResidual);
        linearDiagnostics.refinements += result.refinements;
      };
      recordLinear(linear);
      const normalDirection = linear.x;
      let coordinateDirection = null, coordinateRecoveryFailure;
      if (tangentialGridRecovery) {
        try {
          const raw = proposeDensityNewton(system, state, linear.x, { stagnationLimiter });
          const limit = limitStreamtubeGridStep(value.nodes, system.decode(raw.x).nodes);
          if (limit.limited && limit.step * raw.step < 1e-3) {
            const candidate = tangentialNewton(system, state, value, jacobian, linear.x, equationOrder);
            candidate?.linear?.forEach(recordLinear);
            if (candidate?.direction) { coordinateDirection = candidate; linear.x = candidate.direction; }
            else coordinateRecoveryFailure = candidate?.diagnostics;
          }
        } catch (error) { coordinateRecoveryFailure = { message: error.message, code: error.code }; }
      }
      const rejections = []; let maximumStep = 1, accepted, redistributionBacktracks = 0;
      let geometryRedistribution = false, redistributionRecovery = false, gridRepairRecovery = false, searchStart = 0;
      let gridCorrectionBacktracks = 0;
      for (let trial = 0; ; trial++) {
        let redistributionAttempted = false, dekinkAttempted = false, rawGridRepair = false, inletGridRepair;
        if(traceUpdates){activeTrace={iteration,trial,mcrit:input.upwind?.mcrit,mucon:input.upwind?.mucon,
          beforeSquaredNorm:value.residual.reduce((a,v)=>a+v*v,0),stages:[]};updateTrace.push(activeTrace);}

        try {
          stage = 'density proposal'; proposal = proposeDensityNewton(system, state, linear.x, { stagnationLimiter, maximumStep });
          const decoded = system.decode(proposal.x);
          if (coordinateDirection) decoded.nodes = coordinateDirection.mix(decoded.nodes, proposal.step);
          if(activeTrace){activeTrace.step=proposal.step;activeTrace.limiter=structuredClone(proposal.limiter);
            activeTrace.undampedUpdate=structuredClone(proposal.undampedUpdate);}
          traceStage('raw-newton',proposal.x,decoded.nodes);
          stage = 'Newton grid step';
          if (preserveConvexity) {
            try { requireConvexGridUpdate(value.nodes, decoded.nodes); }
            catch (error) {
              const cell = error.diagnostics?.cell;
              const inletLimited = !coordinateDirection && error.code === 'streamtube-grid-step' && system.layout.bodies.some((b, body) =>
                cell.i < b.leadingIndex && (cell.group === body && cell.tube === decoded.nodes[body][0].length - 2
                  || cell.group === body + 1 && cell.tube === 0));
              if (inletLimited) {
                // Normal motion and inlet re-spacing form one complete
                // trial. The raw corner can fold while their combined path
                // stays convex. Do not shrink to that intermediate boundary
                // before testing the actual proposed inlet coordinates.
                stage = 'inlet grid repair';
                const repaired = reparameterizeStreamtubeInlets(decoded.nodes, system.layout.bodies, fractions);
                try { requireConvexGridUpdate(value.nodes, repaired.nodes); }
                catch (repairError) {
                  // Re-spacing is nonlinear in the Newton step: use ordinary
                  // halving, not an affine root of the unmaintained movement.
                  if (repairError.code === 'streamtube-grid-step') repairError.code = 'streamtube-inlet-grid-step';
                  throw repairError;
                }
                repaired.reparameterization.beforeGridGate = true;
                inletGridRepair = repaired; decoded.nodes = repaired.nodes;
              } else {
                if (!gridRepairRecovery || error.code !== 'streamtube-grid-step') throw error;
                // A trial is not an accepted state. After the ordinary search
                // fails on DEKINK, let that same repair run before deciding
                // whether the complete update has a valid physical grid.
                rawGridRepair = true;
              }
            }
          }
          if (coordinateDirection && !rawGridRepair) geometryRedistribution = false;
          stage = 'inlet adjustment';
          // A compensated tangent proposal already includes the flow response
          // to its new cut spacing. Reapplying the old fractions would undo
          // that coordinate freedom. Capture new fractions only on acceptance.
          const inlet = coordinateDirection ? { nodes: decoded.nodes, maxDisplacement: 0 }
            : inletGridRepair ?? adjustStreamtubeInlets(decoded.nodes, system.layout.bodies, fractions,
              { preserveConvexity: preserveConvexity && !rawGridRepair });
          traceStage('inlet-adjustment',proposal.x,inlet.nodes);
          stage = 'DEKINK'; const dekink = dekinkStreamtubeInteriors(inlet.nodes, { preserveConvexity: preserveConvexity && !rawGridRepair });
          dekinkAttempted = dekink.repairs.length > 0;
          if (rawGridRepair) {
            geometryRedistribution = false;
            try { assertConvexStreamtubeGrid(dekink.nodes); }
            catch { geometryRedistribution = true; }
          }
          traceStage('dekink',proposal.x,dekink.nodes);
          const triggeredBodies = coordinateDirection ? [] : system.layout.bodies.flatMap((b, body) => {
            const lower = dekink.nodes[body][b.leadingIndex + 1].at(-1), upper = dekink.nodes[body + 1][b.leadingIndex + 1][0];
            const spacing = .5 * Math.hypot(upper.x - lower.x, upper.y - lower.y);
            return Math.abs(decoded.stagnation[body] - lastRedistributedStagnation[body]) > .5 * spacing ? [body] : [];
          });
          // Crossing the cumulative stagnation threshold activates a finite
          // coordinate correction. A convexity-rejected trial can also need
          // SMOVE before it can be evaluated.  That is a trial-local repair:
          // it must not be applied unchanged to every shorter Armijo trial.
          // Otherwise the Newton part tends to zero while the finite SMOVE
          // motion remains, making the accepted-state residual discontinuous
          // at a zero line-search step.
          stage = 'SMOVE'; redistributionAttempted = triggeredBodies.length > 0 || geometryRedistribution;
          let r = redistributionAttempted ? redistribute(dekink.nodes, 2 ** -redistributionBacktracks) : { nodes: dekink.nodes, passages: [] };
          if ((redistributionRecovery || gridCorrectionBacktracking === 'before-damping' && stepAcceptance === 'armijo') && redistributionAttempted) {
            let preRedistributionDecreases = false;
            try {
              const repaired = adopt(proposal.x, dekink.nodes);
              requireCoupledResidualDecrease(value.residual, repaired.value.residual, { step: proposal.step, tolerance });
              preRedistributionDecreases = true;
            } catch { /* An invalid or nondecreasing repair still needs the ordinary common-step search. */ }
            // If DEKINK already gave a valid decreasing step, backtrack only
            // the optional tangential correction before shrinking Newton.
            // Otherwise a finite SMOVE jump can hide every usable repair.
            while (preRedistributionDecreases) {
              try {
                const moved = adopt(proposal.x, r.nodes);
                requireCoupledResidualDecrease(value.residual, moved.value.residual, { step: proposal.step, tolerance });
                break;
              } catch (error) {
                if (gridCorrectionBacktracks >= maxBacktracks) break;
                rejections.push({ step: proposal.step, stage: 'grid correction',
                  correctionScale: 2 ** -redistributionBacktracks, message: error.message, code: error.code });
                gridCorrectionBacktracks++; redistributionBacktracks++;
                r = redistribute(dekink.nodes, 2 ** -redistributionBacktracks);
              }
            }
          }
          if (activeTrace && gridRepairRecovery) activeTrace.gridRepairRecovery = true;
          if(activeTrace){activeTrace.redistributionAttempted=redistributionAttempted;activeTrace.gridCorrectionScale=redistributionAttempted?2 ** -redistributionBacktracks:0;}
          traceStage('tangential-redistribution',proposal.x,r.nodes);
          let discardedGridRepair, candidate;
          stage = 'admissibility';
          try { candidate = adopt(proposal.x, r.nodes); }
          catch (error) {
            // Preserve useful SMOVE repairs. If a request inherited from a
            // larger rejected trial instead spoils this valid shorter grid,
            // test the current trial without that obsolete repair. A real
            // stagnation trigger still requires its coordinate correction.
            if (stepAcceptance !== 'armijo' || !geometryRedistribution
              || !redistributionAttempted || triggeredBodies.length) throw error;
            try {
              const unchangedGrid = adopt(proposal.x, dekink.nodes);
              requireCoupledResidualDecrease(value.residual, unchangedGrid.value.residual,
                { step: proposal.step, tolerance });
              candidate = unchangedGrid;
            } catch { throw error; }
            discardedGridRepair = { code: error.code, message: error.message };
            r = { nodes: dekink.nodes, passages: [] };
            geometryRedistribution = false; redistributionAttempted = false;
            if(activeTrace){activeTrace.redistributionAttempted=false;activeTrace.gridCorrectionScale=0;
              activeTrace.discardedGridRepair=discardedGridRepair;}
            traceStage('retained-pre-redistribution',proposal.x,r.nodes);
          }
          let residualDecrease;
          if (stepAcceptance === 'armijo') {
            stage = 'residual decrease';
            try { residualDecrease = requireCoupledResidualDecrease(value.residual, candidate.value.residual,
              { step: proposal.step, tolerance }); }
            catch (error) {
              if (error.code === 'COUPLED_RESIDUAL_DECREASE') {
                error.code = 'EULER_RESIDUAL_DECREASE';
                error.message = error.message.replace('coupled', 'Euler');
              }
              throw error;
            }
          }
          if(activeTrace)activeTrace.accepted=true;
          accepted = { candidate, inlet, dekink, triggeredBodies, r, redistributionAttempted, residualDecrease, discardedGridRepair };
          break;
        } catch (error) {
          if(activeTrace){activeTrace.accepted=false;activeTrace.rejection={stage,message:error.message,code:error.code};}
          rejections.push({ step: proposal?.step, stage, message: error.message,
            ...(error.code === undefined ? {} : { code: error.code }),
            ...(error.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }) });
          if (coordinateDirection && (trial - searchStart >= maxBacktracks || stage === 'density proposal')) {
            coordinateRecoveryFailure = { ...coordinateDirection.diagnostics, accepted: false,
              stage, message: error.message, code: error.code };
            coordinateDirection = null; linear.x = normalDirection;
            maximumStep = 1; searchStart = trial + 1; redistributionBacktracks = 0;
            geometryRedistribution = false; redistributionRecovery = false; gridRepairRecovery = false;
            continue;
          }
          // First preserve the ordinary ISES search. A different tangential
          // grid can change the subsequent viscous startup even when both
          // Euler grids converge. Only after that search exhausts itself on
          // a maintained merit rejection, retry this same Newton direction
          // with a backtracked coordinate correction as well. This second
          // search is bounded by the same budget and retains every final
          // physical, grid and residual acceptance gate.
          if (stepAcceptance === 'armijo' && maxBacktracks > 0 && trial >= maxBacktracks
            && !redistributionRecovery && (redistributionAttempted || dekinkAttempted) && stage === 'residual decrease') {
            redistributionRecovery = true; searchStart = trial + 1;
            gridRepairRecovery = dekinkAttempted;
            maximumStep = 1; redistributionBacktracks = 1; geometryRedistribution = false;
            continue;
          }
          if (stepAcceptance === 'listing' || trial - searchStart >= maxBacktracks || stage === 'density proposal') throw error;
          if (stage === 'Newton grid step' && error.code === 'streamtube-grid-step') geometryRedistribution = true;
          // A triggered redistribution is part of the complete trial too.
          // Leaving it at full strength after a merit rejection creates a
          // finite grid/residual jump even as the Newton step tends to zero.
          if (redistributionAttempted && (stage === 'SMOVE' || stage === 'admissibility'
            || redistributionRecovery && stage === 'residual decrease')) redistributionBacktracks++;
          // A grid-step failure asks the next trial to repair its geometry.
          // Once that repaired trial is admissible but fails Armijo, retry the
          // smaller Newton step without carrying its finite SMOVE correction.
          // A genuinely triggered stagnation redistribution is recomputed
          // from the new trial and remains active.
          if (stage === 'residual decrease') geometryRedistribution = false;
          // Keep the Newton direction and every equation unchanged. Retry
          // from the last accepted state with a smaller common scalar.
          maximumStep = (error.code === 'streamtube-grid-step' ? error.stepFraction : .5) * proposal.step;
        }
      }
      const { candidate, inlet, dekink, triggeredBodies, r, redistributionAttempted } = accepted;
      let densityChange = 0;
      if (adaptiveMcrit) for (let k = 0; k < system.layout.densityCount; k++) densityChange = Math.max(densityChange,
        Math.abs(Math.expm1(candidate.state[k] - state[k])));
      ({ system, state, value } = candidate);
      if (coordinateDirection) fractions = captureStreamtubeInletFractions(value.nodes, system.layout.bodies);
      if (redistributionAttempted || coordinateDirection) lastRedistributedStagnation = value.stagnation.slice();
      let dissipation;
      const linearizedMucon = input.upwind?.mucon;
      if (adaptiveMcrit && !(broadShockStartup && input.upwind.mucon < 0)) {
        // Freeze the law for the entire Jacobian/line search. Only the accepted
        // physical density update schedules the NEXT iteration's broadening.
        // Re-evaluate on the identical state/grid; do not compare Armijo norms
        // belonging to different dissipation laws or accept a temporary root.
        stage = 'dissipation update';
        const linearizedMcrit = input.upwind.mcrit, acceptedResidual = value.diagnostics.residual;
        changeUpwind({ mcrit: acceptedResidual <= tolerance ? targetMcrit : msesTemporaryMcrit(targetMcrit, densityChange) });
        if (value.diagnostics.residual <= tolerance) changeUpwind({ mcrit: targetMcrit });
        dissipation = { mcrit: input.upwind.mcrit, targetMcrit, linearizedMcrit, densityChange, acceptedResidual };
      }
      if (firstOrderStartup) {
        stage = 'dissipation update';
        const firstOrderResidual = value.diagnostics.residual;
        restoreRootDissipation();
        dissipation = { ...dissipation, mucon: input.upwind.mucon, targetMucon, linearizedMucon,
          mcrit: input.upwind.mcrit, targetMcrit: targetMcrit ?? input.upwind.mcrit,
          ...(linearizedMucon !== input.upwind.mucon ? { restoredSecondOrder: true, firstOrderResidual } : {}) };
      }
      const { x, ...step } = proposal;
      const entry = { ...step, ...(coordinateDirection ? { coordinateRecovery: coordinateDirection.diagnostics } : {}),
        ...(coordinateRecoveryFailure ? { coordinateRecoveryFailure } : {}), iteration, residual: value.diagnostics.residual,
        ...(dissipation ? { dissipation } : {}),
        ...(preserveConvexity ? { backtracks: rejections.length, rejections } : {}),
        ...(accepted.residualDecrease ? { residualDecrease: accepted.residualDecrease } : {}),
        maintenance: { inlet: { maxDisplacement: inlet.maxDisplacement, maxArcErrorBefore: inlet.maxArcErrorBefore, maxArcErrorAfter: inlet.maxArcErrorAfter,
            ...(inlet.reparameterization ? { reparameterization: inlet.reparameterization } : {}),
            ...(inlet.convexity ? { convexity: inlet.convexity } : {}) },
          dekinkRepairs: dekink.repairs, ...(dekink.convexity ? { dekinkConvexity: dekink.convexity } : {}),
          ...(redistributionRecovery ? { residualBacktracking: true } : {}),
          ...(gridRepairRecovery ? { gridRepairRecovery: true } : {}),
          ...(gridCorrectionBacktracks ? { gridCorrectionBacktracks } : {}),
          ...(accepted.discardedGridRepair ? { discardedGridRepair: accepted.discardedGridRepair } : {}),
          triggeredBodies, ...(geometryRedistribution ? { geometryRedistribution: true } : {}), passages: r.passages } };
      history.push(entry); onIteration?.(entry); onMesh?.({ system, nodes: value.nodes, flow: value, iteration: { ...entry } });
    } catch (error) {
      rejectUpdate(error, stage, proposal);
    }
    if (!failed) {
      emitCheckpoint();
      const stalled = stopOnGridStagnation ? repeatedGridStepLimit(history, tolerance) : null;
      const passageStall = stopOnGridStagnation && system.layout.bodies.length > 1
        && system.conditions.mach <= .3 && value.diagnostics.maxMach < 1
        ? repeatedPassageStepLimit(history, tolerance) : null;
      if (!stalled && passageStall) rejectUpdate(Object.assign(new Error('Repeated small passage updates made no useful residual progress.'), {
        code: 'EULER_PASSAGE_STAGNATION', diagnostics: passageStall,
      }), 'passage progress', proposal);
      if (stalled) rejectUpdate(Object.assign(new Error('Repeated convexity-limited steps made no useful residual progress.'), {
        code: 'EULER_GRID_STAGNATION', diagnostics: stalled,
      }), 'grid progress', proposal);
    }
  }
  const result = solveStreamtubeBody(system, { initial: state, maxIterations: 0, tolerance });
  const finalQuality = streamtubeMeshSnapshot({ system, nodes: result.nodes }).quality;
  const requestedLaw = requestedDissipation();
  const converged = result.converged && requestedLaw && finalQuality.valid;
  return { ...result, converged, residualConverged: result.converged && requestedLaw, finalQuality, history,
    reason: converged ? 'residual' : !requestedLaw && result.converged ? 'Temporary dissipation residual met; requested MCRIT or second-order dissipation has not been restored.'
      : result.converged ? 'Euler residual met; final grid fails the convexity/quality gate.' : reason, linearDiagnostics, initialRedistribution,
    lastRejectedStep, solverInput: input, stagnationLimiter, iterationGeometry, stepAcceptance,
    ...(adaptiveMcrit ? { adaptiveMcrit, targetMcrit } : {}),
    ...(firstOrderStartup ? { firstOrderStartup, targetMucon } : {}),
    ...(broadShockStartup ? { broadShockStartup, broadShockMucon } : {}),
      ...(stopOnGridStagnation ? { stopOnGridStagnation } : {}),
    ...(gridCorrectionBacktracking !== 'after-search' ? { gridCorrectionBacktracking } : {}),
    ...(preserveConvexity ? { gridAcceptance: 'convex' } : {}),
    ...((retainCheckpoint || onCheckpoint || resume) && initialRedistribution.accepted ? { checkpoint: checkpoint() } : {}),
    ...(bestCheckpoint ? { bestCheckpoint } : {}),
    ...(traceUpdates ? { updateTrace } : {}),
    formulation: `Research subcritical simultaneous Euler with additive-density Newton, common scalar damping, inlet arc adjustment, DEKINK and initial/triggered five-pair tangential redistribution; iteration geometry ${iterationGeometry}, step acceptance ${stepAcceptance}, final convexity/quality required; explicit multielement passage extension; no evolving BL or shock dissipation.` };
}

// Solve at unchanged physical conditions and resolution: first-order startup,
// then the prescribed dissipation. The caller chooses when to use this sequence.
// MCRIT is frozen within each Newton solve and updated only from its accepted
// physical density increment. maxIterations is a per-phase budget; the total
// is at most twice that value for a requested second-order solution. This
// helper does not select an app default or certify a temporary-law root.
export function continueEulerDissipationFromCheckpoint(source, {
  maxIterations = 40, tolerance = 1e-10, onIteration, onCheckpoint,
} = {}) {
  const target = structuredClone(source?.input?.upwind);
  if (!target || !Number.isFinite(target.mucon) || target.mucon === 0
    || !Number.isInteger(maxIterations) || maxIterations < 1 || !Number.isFinite(tolerance) || tolerance <= 0
    || onIteration !== undefined && typeof onIteration !== 'function'
    || onCheckpoint !== undefined && typeof onCheckpoint !== 'function')
    throw new Error('Invalid Euler dissipation startup controls.');
  msesTemporaryMcrit(target.mcrit, 0);
  // Validate the saved residual, physical nodes and iteration controls before
  // changing the equation parameters. Rebuilding a corrupt residual first
  // would silently bypass the ordinary checkpoint replay guard.
  const initial = solveStreamtubeIses(undefined, { resume: source, ...source.continuation,
    maxIterations: 0, tolerance, retainCheckpoint: true });
  let checkpoint = initial.checkpoint, reason = 'iteration limit', lastResult, converged = false;
  const history = [], started = performance.now();
  function rebuild() {
    const system = createStreamtubeBodySystem(checkpoint.input);
    const x = system.adoptGeometry(Float64Array.from(checkpoint.initialEuler.x), checkpoint.initialEuler.nodes);
    const value = system.evaluate(x);
    // Equation changes must leave the physical state untouched.
    if (x.some((v, i) => v !== checkpoint.initialEuler.x[i])
      || JSON.stringify(value.nodes) !== JSON.stringify(checkpoint.initialEuler.nodes))
      throw new Error('Dissipation change moved the source physical state.');
    checkpoint.residual = Array.from(value.residual);
    return { system, value };
  }
  const phases = target.mucon > 0 ? [-target.mucon, target.mucon] : [target.mucon];
  for (const mucon of phases) {
    checkpoint.input.upwind = { ...target, mucon };
    let phaseConverged = rebuild().value.diagnostics.residual <= tolerance;
    for (let iteration = 1; !phaseConverged && iteration <= maxIterations; iteration++) {
      const before = checkpoint, mcrit = checkpoint.input.upwind.mcrit;
      lastResult = solveStreamtubeIses(undefined, { resume: checkpoint, ...checkpoint.continuation,
        maxIterations: 1, tolerance, retainCheckpoint: true });
      if (lastResult.history.length < 2) {
        reason = lastResult.reason;
        history.push({ mucon, iteration, mcrit, residual: lastResult.diagnostics.residual, step: 0, reason });
        break;
      }
      checkpoint = lastResult.checkpoint;
      const system = createStreamtubeBodySystem(checkpoint.input);
      let densityChange = 0;
      for (let k = 0; k < system.layout.densityCount; k++) densityChange = Math.max(densityChange,
        Math.abs(Math.expm1(checkpoint.initialEuler.x[k] - before.initialEuler.x[k])));
      const nextMcrit = msesTemporaryMcrit(target.mcrit, densityChange);
      checkpoint.input.upwind.mcrit = nextMcrit;
      const { value } = rebuild(), h = lastResult.history.at(-1);
      const entry = { mucon, iteration, mcrit, nextMcrit, densityChange, step: h.step, backtracks: h.backtracks,
        residual: lastResult.diagnostics.residual, nextResidual: value.diagnostics.residual, quality: lastResult.finalQuality };
      history.push(entry);
      onIteration?.(structuredClone(entry)); onCheckpoint?.(structuredClone(checkpoint));
      if (value.diagnostics.residual <= tolerance && nextMcrit === target.mcrit) phaseConverged = true;
    }
    if (!phaseConverged) break;
    if (mucon === target.mucon) { converged = true; reason = 'prescribed residual'; }
  }
  return { converged, reason, history, checkpoint, lastRejectedStep: lastResult?.lastRejectedStep,
    seconds: (performance.now() - started) / 1000 };
}
