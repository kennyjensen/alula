// SPDX-License-Identifier: GPL-2.0-or-later
// Diagnostic entry point only. This function is serialized into a Chrome
// module Worker by capture-cold-browser-checkpoints.js; production is unedited.
export function coldBrowserCheckpointObservers(emit) {
  return {
    onEulerPrepared: packet => {
      // This handoff is process-local and includes a live system with methods.
      // capturePreparedStreamtubeAssembly has already structured-cloned all
      // remaining data; preserve those values without its uncloneable identity.
      const { system, ...prepared } = packet.prepared;
      emit({ type: 'prepared-euler', packet: { version: packet.version, caseData: packet.caseData, prepared },
        processLocalSystemOmitted: true });
    },
    onCheckpoint: observation => emit({ type: 'assembly-checkpoint', observation }),
    onIterationCheckpoint: (checkpoint, details) => emit({ type: 'coupled-checkpoint', checkpoint, details }),
  };
}

export function coldBrowserCheckpointWorker(checkpointObservers) {
  self.onmessage = async ({ data: { caseData, origin } }) => {
    let sequence = 0;
    const started = performance.now(), before = JSON.stringify(caseData);
    const emit = data => {
      const fields = { sequence: ++sequence, workerSeconds: (performance.now() - started) / 1000 };
      try { self.postMessage({ ...fields, ...data }); }
      catch (error) {
        // An observer serialization failure must not reject a numerical step.
        // Report the missing capture explicitly; its receipt cannot qualify.
        self.postMessage({ ...fields, type: 'observation-error', failedType: data.type, message: error.message });
      }
    };
    const errorFields = error => ({ message: error?.message ?? String(error), stack: error?.stack,
      code: error?.code, stage: error?.stage, diagnostics: error?.diagnostics,
      initialization: error?.initialization, cause: error?.cause && { message: error.cause.message,
        code: error.cause.code, stage: error.cause.stage, diagnostics: error.cause.diagnostics } });
    try {
      if (!caseData?.quadBoundaryLayers || caseData.flowModel !== 'streamtube-grid' || !(caseData.mach <= .3))
        throw new Error('This diagnostic captures the ordinary low-Mach cold coupled assembly route only.');
      const { solveCoupledStreamtubeAssembly } = await import(`${origin}/src/euler/streamtube-coupled-assembly.js`);
      const { quadCoupledResultForDisplay } = await import(`${origin}/src/ui/quad-coupled-result.js`);
      const { quadCoupledIterationPressure } = await import(`${origin}/src/ui/quad-coupled-iteration-pressure.js`);
      const { quadCoupledNcrit, quadCoupledNcritResult } = await import(`${origin}/src/ui/quad-coupled-ncrit.js`);
      const { createQuadMeshProgress } = await import(`${origin}/src/ui/quad-mesh-progress.js`);
      const progress = createQuadMeshProgress(caseData.referenceChord);
      let latestMesh;
      // These budgets and observable calculations match the ordinary low-Mach
      // branch in src/worker/solver.js. Extra callbacks only serialize states
      // already produced by the public adapter. Never evaluate a checkpoint.
      const raw = solveCoupledStreamtubeAssembly(caseData, {
        maxIterations: caseData.maxIterations ?? (caseData.transitionMode === 'automatic' ? 40 : 20),
        eulerMaxIterations: caseData.eulerMaxIterations ?? 20,
        ...checkpointObservers(emit),
        onStage: stage => { progress.stage(stage); emit({ type: 'flow-stage', ...stage }); },
        onIteration: iteration => emit({ type: 'iteration', iteration }),
        onMesh: (mesh, stage) => { latestMesh = progress.mesh(mesh); emit({ type: 'mesh', mesh: latestMesh, stage }); },
        onFlow: frame => {
          const fields = { iteration: frame.iteration?.iteration ?? frame.iteration,
            mach: frame.checkpoint?.restart.input.mach ?? frame.mach, targetMach: frame.targetMach ?? caseData.mach,
            stage: frame.stage, coarseStage: frame.coarseStage, startupAttempt: frame.startupAttempt,
            ...(frame.gridLevel === undefined ? {} : { gridLevel: frame.gridLevel }),
            ...(frame.requestedGridIntervals === undefined ? {} : { requestedGridIntervals: frame.requestedGridIntervals }),
            ...(frame.retained === undefined ? {} : { retained: frame.retained }),
            ...quadCoupledNcrit(frame, caseData.ncrit ?? 9) };
          try {
            const pressure = quadCoupledIterationPressure({ ...frame, elements: caseData.elements,
              referenceChord: frame.referenceChord ?? frame.normalization?.referenceChord ?? caseData.referenceChord });
            emit({ type: 'pressure', pressure, ...fields });
          } catch (error) { emit({ type: 'pressure-unavailable', message: error.message, ...fields }); }
        },
      });
      const displayed = quadCoupledResultForDisplay(raw, caseData, latestMesh);
      const result = quadCoupledNcritResult({ ...displayed,
        sourceCase: structuredClone(displayed.sourceCase ?? raw.sourceCase ?? caseData) }, caseData);
      emit({ type: 'result', result, inputUnchanged: before === JSON.stringify(caseData) });
    } catch (error) { emit({ type: 'error', ...errorFields(error), inputUnchanged: before === JSON.stringify(caseData) }); }
  };
}
