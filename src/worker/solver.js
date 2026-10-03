// SPDX-License-Identifier: GPL-2.0-or-later
import { solveInviscid } from '../inviscid/linear-vortex.js';
import { traceStreamlines } from '../inviscid/streamlines.js';
import { solveViscousAssembly } from '../viscous/result.js';
import { quadCoupledIterationPressure } from '../ui/quad-coupled-iteration-pressure.js';
import { quadCoupledNcrit, quadCoupledNcritResult } from '../ui/quad-coupled-ncrit.js';
import { quadCoupledIterationCoefficients as iterationLoads } from '../ui/quad-coupled-iteration-coefficients.js';
import { quadCoupledTransonicCoefficients as transonicLoads } from '../ui/quad-coupled-transonic-coefficients.js';

self.onmessage = async ({ data }) => {
  const { id, caseData, bounds, meshOnly = false } = data;
  try {
    const start = performance.now();
    const pressureFields = frame => ({ iteration: frame.iteration?.iteration ?? frame.iteration,
      mach: frame.checkpoint?.restart.input.mach ?? frame.mach, targetMach: frame.targetMach ?? caseData?.mach,
      actualAlpha: frame.checkpoint?.restart.input.alpha ?? frame.actualAlpha, targetAlpha: frame.targetAlpha ?? caseData?.alpha,
      stage: frame.stage, coarseStage: frame.coarseStage, startupAttempt: frame.startupAttempt,
      ...(frame.checkpoint?.continuation ? { shearCoordinate: frame.checkpoint.continuation.shearCoordinate ?? 'linear' } : {}),
      ...(frame.checkpoint?.restart?.options ? { hkFloorLinearization: frame.checkpoint.restart.options.hkFloorLinearization ?? 'exact' } : {}),
      ...(frame.gridLevel === undefined ? {} : { gridLevel: frame.gridLevel }),
      ...(frame.requestedGridIntervals === undefined ? {} : { requestedGridIntervals: frame.requestedGridIntervals }),
      ...(frame.retained === undefined ? {} : { retained: frame.retained }),
      ...quadCoupledNcrit(frame, caseData?.ncrit ?? 9) });
    const pressureInput = frame => ({ ...frame, elements: caseData.elements,
      referenceChord: frame.referenceChord ?? frame.normalization?.referenceChord ?? caseData.referenceChord });
    // Loads and Cp are independent observables. A failed wake-drag estimate
    // must not hide usable surface pressure or abort the numerical solve.
    const publishCoefficients = (frame, calculate) => {
      let pressure, pressureError;
      try { pressure = quadCoupledIterationPressure(pressureInput(frame)); } catch (error) { pressureError = error; }
      try { self.postMessage({ id, type: 'coefficients', coefficients: calculate(), ...(pressure ? { pressure } : {}), ...pressureFields(frame) }); }
      catch (error) {
        self.postMessage({ id, type: 'coefficients-unavailable', message: error.message, ...pressureFields(frame) });
        if (pressure) self.postMessage({ id, type: 'pressure', pressure, ...pressureFields(frame) });
      }
      if (pressureError) self.postMessage({ id, type: 'pressure-unavailable', message: pressureError.message, ...pressureFields(frame) });
    };
    const publishObservedFlow = observed => {
      const frame = { ...observed,
        solverLength: observed.solverLength ?? observed.normalization?.solverLength,
        referenceChord: observed.referenceChord ?? observed.normalization?.referenceChord ?? caseData.referenceChord,
        momentReference: caseData.momentReference ?? { x: (caseData.referenceChord ?? 1) / 4, y: 0 },
        alpha: observed.checkpoint.restart.input.alpha ?? 0,
        gamma: observed.checkpoint.restart.input.gamma ?? 1.4 };
      publishCoefficients(frame, () => frame.checkpoint.restart.options.blThermodynamics === 'historical-common-isentrope'
        ? transonicLoads(frame) : iterationLoads(frame));
    };
    const publishEulerCoefficients = mesh => {
      const frame = mesh.coefficientProgress;
      if (!frame) return; // Geometry-only initialization has no flow loads.
      self.postMessage({ id, type: frame.coefficients ? 'coefficients' : 'coefficients-unavailable',
        ...frame, ...pressureFields({ ...mesh.iteration, ...frame,
          stage: mesh.iteration?.stage ?? 'quad-euler' }) });
    };
    if (data.task === 'audit-grid') {
      const { auditStreamtubeMesh } = await import('../geometry/streamtube-grid-audit.js');
      const audit = auditStreamtubeMesh(data.mesh, { onProgress: progress => self.postMessage({ id, type: 'grid-audit-progress', progress }) });
      self.postMessage({ id, type: 'grid-audit', audit, elapsed: performance.now() - start });
      return;
    }
    if (data.task === 'refine-coupled') {
      const { solveCoupledStreamtubeRefinement } = await import('../euler/streamtube-coupled-refinement-assembly.js');
      const { quadCoupledResultForDisplay } = await import('../ui/quad-coupled-result.js');
      const { createQuadMeshProgress } = await import('../ui/quad-mesh-progress.js');
      const progress = createQuadMeshProgress(caseData.referenceChord); let latestMesh, prepared, checkpoint;
      const raw = solveCoupledStreamtubeRefinement(caseData, data.parentResult, {
        onPrepared: value => { prepared = value; checkpoint = value.checkpoint; },
        onCheckpoint: value => { checkpoint = value; },
        onStage: stage => { progress.stage(stage); self.postMessage({ id, type: 'flow-stage', ...stage }); },
        onMesh: (mesh, stage, state) => {
          latestMesh = progress.mesh(mesh); self.postMessage({ id, type: 'mesh', mesh: latestMesh, stage });
          if (state && checkpoint) {
            const { system, settings } = prepared;
            const frame = { checkpoint, flow: state.flow, bl: system.bl,
              bodies: system.euler.layout.bodies, solverLength: settings.normalization.solverLength,
              referenceChord: settings.normalization.referenceChord, momentReference: caseData.momentReference,
              alpha: caseData.alpha ?? 0, mach: checkpoint.restart.input.mach, gamma: system.euler.conditions.gamma,
              iteration: state.iteration, stage: 'coupled-refinement' };
            publishObservedFlow(frame);
          }
        },
        onIteration: iteration => self.postMessage({ id, type: 'iteration', iteration }),
      });
      const result = quadCoupledResultForDisplay(raw, caseData, latestMesh);
      self.postMessage({ id, type: 'result', result, elapsed: performance.now() - start }); return;
    }
    if(caseData.flowModel==='streamtube-grid'){
      if(caseData.quadBoundaryLayers&&!meshOnly){
        const {solveCoupledStreamtubeAssembly}=await import('../euler/streamtube-coupled-assembly.js');
        const {quadCoupledResultForDisplay}=await import('../ui/quad-coupled-result.js');
        const {createQuadMeshProgress}=await import('../ui/quad-mesh-progress.js');
        const progress=createQuadMeshProgress(caseData.referenceChord);let latestMesh;
        const raw=solveCoupledStreamtubeAssembly(caseData,{
          direct:true,
          convergence:'mses',
          // Natural-transition interval changes can use most of the startup
          // budget before the final full Newton steps. Explicit caps still win.
          maxIterations:caseData.maxIterations??(caseData.transitionMode==='automatic'?40:20),eulerMaxIterations:caseData.eulerMaxIterations??40,
          onStage:stage=>{progress.stage(stage);self.postMessage({id,type:'flow-stage',...stage});},
          onFlow:publishObservedFlow,
          onMesh:(mesh,stage)=>{latestMesh=progress.mesh(mesh);self.postMessage({id,type:'mesh',mesh:latestMesh,stage});publishEulerCoefficients(latestMesh);},
          onIteration:iteration=>self.postMessage({id,type:'iteration',iteration})});
        // The assembly retains a structured startup failure so callers can
        // inspect its last admissible mesh. The Worker error channel is the
        // public terminal result for that state; preserve its typed cause.
        if(raw.failure){
          const failure=raw.failure, envelope=raw.initialization!==undefined||failure.failure?raw:failure, detail=envelope.failure??envelope;
          throw Object.assign(new Error(envelope.reason??raw.reason??detail.reason??'Coupled solve failed.'),{
            code:detail.code,stage:detail.stage,diagnostics:{...envelope.diagnostics,...detail.diagnostics,
              ...(envelope.initialization?.coldMachStartup?{coldMachStartup:envelope.initialization.coldMachStartup}:{})},
            cause:detail.cause,initialization:envelope.initialization??raw.initialization,
            actualMach:envelope.actualMach??raw.actualMach,targetMach:envelope.targetMach??raw.targetMach});
        }
        const displayed=quadCoupledResultForDisplay(raw,caseData,latestMesh);
        const result=quadCoupledNcritResult({...displayed,
          sourceCase:structuredClone(displayed.sourceCase??raw.sourceCase??caseData)},caseData);
        self.postMessage({id,type:'result',result,elapsed:performance.now()-start});return;
      }
      const {solveStreamtubeAssembly}=await import('../euler/streamtube-result.js');
      const result=solveStreamtubeAssembly(caseData,{meshOnly,maxIterations:caseData.maxIterations,
        onMesh:(mesh,stage)=>{self.postMessage({id,type:'mesh',mesh,stage});publishEulerCoefficients(mesh);},
        onIteration:iteration=>self.postMessage({id,type:'iteration',iteration})});
      if(meshOnly)self.postMessage({id,type:'mesh-ready',elapsed:performance.now()-start});
      else self.postMessage({id,type:'result',result,elapsed:performance.now()-start});
      return;
    }
    const subcritical=caseData.flowModel==='subcritical',coupled=subcritical||caseData.flowModel==='coupled';
    // Install the message handler before loading the WASM-backed module.
    // A top-level await in the static import graph can otherwise discard the
    // first posted message while the module worker is still initializing.
    const controls={onIteration:iteration=>self.postMessage({id,type:'iteration',iteration}),
      onMesh:mesh=>self.postMessage({id,type:'mesh',mesh,stage:'solving'})};
    if(meshOnly&&!subcritical)throw new Error('The selected panel flow model has no volume mesh.');
    if(subcritical){
      const {buildSubcriticalMeshPreview}=await import('../potential/mesh-preview.js');
      try{
        self.postMessage({id,type:'mesh',mesh:buildSubcriticalMeshPreview(caseData,{mesh:{surfaceScale:caseData.outerMeshScale??.3}}),stage:'initial'});
      }catch(error){
        // A cold wake guide can fail even if the assembly initialization
        // later gives a valid mesh. A display preview must not veto that solve.
        if(meshOnly)throw error;
        self.postMessage({id,type:'mesh-unavailable',message:error.message});
      }
      if(meshOnly){self.postMessage({id,type:'mesh-ready',elapsed:performance.now()-start});return;}
    }
    const solveSubcriticalAssembly=subcritical?(await import('../potential/result.js')).solveSubcriticalAssembly:null;
    const result=subcritical?solveSubcriticalAssembly(caseData,{...controls,mesh:{surfaceScale:caseData.outerMeshScale??.3}})
      :coupled?solveViscousAssembly(caseData,controls):solveInviscid(caseData);
    const elapsed = performance.now() - start;
    self.postMessage({ id, type: 'result', result, elapsed });
    if (coupled) return;
    try {
      const lines = traceStreamlines(result, bounds);
      self.postMessage({ id, type: 'streamlines', lines });
    } catch (error) { self.postMessage({ id, type: 'streamline-error', message: error.message }); }
  } catch (error) {
    const precursor = error?.initialization?.result;
    self.postMessage({ id, type: 'error', message: error?.message ?? String(error), code: error?.code,
      diagnostics: error?.diagnostics, stage: error?.stage, stack: error?.stack,
      actualMach: error?.actualMach ?? precursor?.conditions?.mach ?? precursor?.mach ?? error?.diagnostics?.mach,
      targetMach: error?.targetMach ?? caseData?.mach,
      ...quadCoupledNcrit({ ...error, checkpoint: precursor?.checkpoint,
        conditions: precursor?.conditions, actualNcrit: error?.actualNcrit ?? error?.diagnostics?.actualNcrit }, caseData?.ncrit ?? 9),
      ...(precursor ? { precursor: { reason: precursor.diagnostics?.reason,
        solverStopReason: precursor.diagnostics?.solverStopReason ?? precursor.flow?.reason,
        lastRejectedStep: precursor.flow?.lastRejectedStep, quality: precursor.mesh?.quality } } : {}),
      ...(error?.cause ? { cause: { message: error.cause.message, code: error.cause.code,
        stage: error.cause.stage, diagnostics: error.cause.diagnostics } } : {}) });
  }
};
