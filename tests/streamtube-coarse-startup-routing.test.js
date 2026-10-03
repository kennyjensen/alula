import { coupledResolutionRecoveryCase } from '../src/euler/streamtube-coupled-resolution-recovery.js';
import { coupledMachRecoveryCase } from '../src/euler/streamtube-coupled-mach-recovery.js';
import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';
import { planCoupledLogarithmicShearRecovery } from '../src/euler/streamtube-coupled-log-shear-recovery.js';
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {streamtubeEquationControls} from '../src/euler/streamtube-equation-selection.js';
import {retainBestCoupledCheckpoint} from '../src/euler/streamtube-transition-recovery.js';
import {coupledStartupExtensionPlan} from '../src/euler/streamtube-coupled-startup.js';
import {coupledNcritStartupRecoveryPlan} from '../src/euler/streamtube-coupled-ncrit-recovery.js';
const draftPath=new URL('../src/euler/streamtube-coupled-assembly.js',import.meta.url);
const draft=fs.readFileSync(draftPath,'utf8'), baseline=fs.readFileSync(new URL('../docs/rae2822/coarse-startup-integration/assembly-before.js.txt',import.meta.url),'utf8');
console.log(JSON.stringify({draftSha256:createHash('sha256').update(draft).digest('hex'),operations:'All numerical dependencies replaced by controlled stubs; zero flow/geometry/BL evaluations.'}));
let sequence=0;
const clone=structuredClone;
// The new cold assembly explicitly selects the tested TE-following chart.
// Public metadata now spells out the legacy linear/exact policies.
// Archived routing controls otherwise remain exact; no checkpoint replay is
// exercised by these cold-only stubs.
const withFreshWakeChart=value=>{
  const next=clone(value);
  const visit=v=>{
    if(!v||typeof v!=='object'||ArrayBuffer.isView(v))return;
    if(v.solverSettings){v.solverSettings.convergence ??= 'residual';v.solverSettings.shearCoordinate ??= 'linear';v.solverSettings.hkFloorLinearization ??= 'exact';v.solverSettings.automaticResolutionRecovery ??= true;v.solverSettings.automaticMachRecovery ??= true;}
    if(v.wakeGeometry==='independent-banks'&&v.wakeOutlet==='banks'&&!Object.hasOwn(v,'stepMethod'))
      v.wakeDisplacementMotion='te-center';
    for(const child of Object.values(v))visit(child);
  };visit(next);return next;
};
const baseCase=()=>({elements:[{name:'test',points:[{x:1,y:0},{x:0,y:0}]}],gridIntervals:32,gridTubes:11,
  mach:.2,alpha:2.68,reynolds:2.7e6,ncrit:4,referenceChord:1,transitionMode:'automatic',materialTrips:[[1,1]],gridEllipticSmoothing:true});
const finiteCase=()=>{const input=baseCase();input.elements[0].trailingEdge={kind:'finite-base',upperIndex:0,lowerIndex:1};return input;};
const meshFor=(grid,kind)=>({tag:`${grid}:${kind}`,quality:{valid:true},cells:Array(grid).fill([0,1,2,3]),initialization:{gridSmoothing:{converged:true}}});
async function harness(config={},source=draft){
  const calls={euler:[],eulerBindings:[],initializers:[],solves:[],solveConditions:[],maps:[],sequences:[],sequenceConditions:[],continuations:[],evaluations:[],checkpoints:[],geometry:[],captures:[],restores:[],recoveries:[],levels:[],levelSolves:[],profilePolicies:[],continuationGridFrames:[],ncritRecoveries:[]};
  const nodes=grid=>[[[{x:grid,y:0},{x:grid,y:1}],[{x:grid+1,y:0},{x:grid+1,y:1}]]];
  const flow=grid=>({nodes:nodes(grid),undisplacedNodes:nodes(grid),diagnostics:{maxMach:.2},cells:[[[{interfacePressure:{lower:grid,upper:grid+1}}]]]});
  const makePrepared=(input,options,mapped=false)=>{
    const grid=input._grid,factor=mapped?1:(config.factor??(grid>16?.25:1)),m=meshFor(grid,mapped?'mapped':'prepared');
    const initial=Float64Array.from([grid,0,1,2,3]), ne=1;
    const bl={transitionMode:options.transitionMode??'fixed-trip',hasFiniteBase:!!config.finite,
      trips:clone(options.tripFractions),snapshotActive:()=>[2,3],
      stations:[{id:grid,kind:'surface',body:0,side:'upper',i:grid}],surfaces:[{body:0,side:'upper',ids:[grid]}],wakes:[]};
    return{initialization:{thicknessFactor:factor,...(mapped?{transfer:{mapped:true}}:{})},mesh:m,
      system:{ne,n:config.pairedFlow?grid+5:5,initial,bl,conditions:{ncrit:options.ncrit,hkFloorLinearization:options.hkFloorLinearization},
      euler:{layout:{nx:grid,tubes:[3,4],bodies:input.bodies},decode:()=>({nodes:nodes(grid)})},evaluate:()=>{
        calls.evaluations.push({grid,mapped});
        if(grid>16&&config.failFineEvaluation)throw new Error('fine seed rejected');
        return{families:{euler:0,boundaryLayer:.2,edgeMatching:0},outer:flow(grid)};
      }}};
  };
  const eulerInput=input=>({_grid:input.gridIntervals,bodies:input.elements.map((e,element)=>({element,...(e.trailingEdge?{trailingEdge:clone(e.trailingEdge)}:{})})),
    mach:input.mach,streamwiseMode:'isentropic',...streamtubeEquationControls(input.eulerIsmom)});
  const geometry=input=>({input:eulerInput(input),mesh:meshFor(input.gridIntervals,'geometry'),
    system:{baseGeometry:config.finite?[{}]:[],conditions:{lengthScale:1}}});
  const stubs={ coupledResolutionRecoveryCase, coupledMachRecoveryCase,
    retryCoupledEulerSeed: ({ onIteration }) => {
      if (!config.failedSeedRecovery) return null;
      onIteration?.({ iteration: 1, residual: 99, euler: 99, boundaryLayer: 98, edgeMatching: 97 });
      return { diagnostics: { accepted: false, reason: 'controlled seed failure', iterations: 1 } };
    }, stalledBoundaryLayerResolutionPlan: () => null,
    ...shearPolicy, planCoupledLogarithmicShearRecovery, observableFlow, observableBL,
    continueCoupledWakeCoordinates: result => result,
    streamtubeEquationControls,
    retainBestCoupledCheckpoint,
    coupledNcritStartupRecoveryPlan,
    coupledStartupExtensionPlan,
    finishCoupledNcritStartup(source,options){
      calls.ncritRecoveries.push({source:clone(source),plan:clone(options.plan),maxIterations:options.maxIterations});
      const result=clone(source),ncrit=config.recoveryActualNcrit??options.plan.targetNcrit;
      assert.equal(source.converged,true);assert.equal(source.checkpoint.restart.options.ncrit,4);
      options.onStage?.({stage:'coupled',actualNcrit:ncrit,targetNcrit:options.plan.targetNcrit});
      options.onIteration?.({iteration:1,actualNcrit:ncrit,targetNcrit:options.plan.targetNcrit});
      result.conditions.ncrit=ncrit;result.checkpoint.restart.options.ncrit=ncrit;result.coupledOptions.ncrit=ncrit;
      result.automaticRefinement={kind:'transition-local',refinedNx:23};
      result.mesh.initialization.gridRefinement=result.automaticRefinement;
      options.onIterationCheckpoint?.(result.checkpoint,{stage:'coupled',actualNcrit:ncrit,targetNcrit:options.plan.targetNcrit});
      options.onFlow?.({checkpoint:result.checkpoint,flow:result.flow,actualNcrit:ncrit,targetNcrit:options.plan.targetNcrit});
      return{result,diagnostics:{attempted:true,reachedTarget:ncrit===options.plan.targetNcrit,actualNcrit:ncrit,
        targetNcrit:options.plan.targetNcrit,stateConverged:true}};
    },
    continueCoupledNcrit(source,options){
      calls.continuations.push({phase:options.phase,sourceNcrit:source.checkpoint.restart.options.ncrit,
        targetNcrit:options.targetNcrit,sourceConverged:source.converged,sourceGrid:source.solverInput._grid,
        sourceHkPolicy:source.checkpoint.restart.options.hkFloorLinearization});
      assert.equal(source.converged,true,'Continuation must start from a complete converged source.');
      const result=clone(source);
      const actualNcrit=(options.phase==='coarse'?config.coarseRetainedNcrit:config.fineRetainedNcrit)??options.targetNcrit;
      result.checkpoint.restart.options.ncrit=actualNcrit;
      result.coupledOptions.ncrit=actualNcrit;
      result.conditions={...result.conditions,ncrit:actualNcrit};
      if(options.phase==='fine'&&config.distinctContinuationRoot){
        // Keep the actual observer-provenance object: final reporting must
        // replace it rather than retroactively relabel the earlier root.
        calls.continuationGridFrames.push(source.gridSequence);
        result.families={euler:8e-12,boundaryLayer:7e-12,edgeMatching:6e-12};
        result.checkpoint.families=clone(result.families);
        result.history=Array.from({length:5},(_,iteration)=>({iteration,residual:iteration===4?8e-12:1}));
        result.reason='continued root residual';
      }
      options.onStage?.({stage:options.phase==='coarse'?'coupled-coarse-initialization':'coupled',
        actualNcrit,targetNcrit:options.targetNcrit,ncritStartup:true});
      return {result,diagnostics:{attempted:true,phase:options.phase,actualNcrit,targetNcrit:options.targetNcrit,
        reachedTarget:actualNcrit===options.targetNcrit,stateConverged:true,attempts:[]}};
    },
    prepareStreamtubeAssembly(input,options){calls.geometry.push(clone(input));const prepared=geometry(input);
      options.onMesh?.(prepared.mesh,'initial');return prepared;},
    capturePreparedStreamtubeAssembly(prepared,input){calls.captures.push(input.gridIntervals);return{prepared,caseData:clone(input)};},
    restorePreparedStreamtubeAssembly(packet,input){calls.restores.push(input.gridIntervals);return packet.prepared;},
    solveStreamtubeAssembly(input,options){
      calls.eulerBindings.push({grid:input.gridIntervals,preparedEuler:options.preparedEuler,onEulerPrepared:options.onEulerPrepared});
      calls.euler.push(clone(input));const grid=input.gridIntervals,m=meshFor(grid,'euler');
      options.onMesh?.(m,'initial');
      if(grid===16&&config.translateSmoothingObserver){
        try{options.onMesh?.(m,'smoothing');}catch(error){throw new Error(`Elliptic grid initialization failed: ${String(error?.message)}`);}
      }
      if(grid>16&&config.fineGasError)throw config.fineGasError;
      if(grid===16&&config.coarseGasError)throw config.coarseGasError;
      options.onIteration?.({iteration:0,residual:1});
      const solverInput=eulerInput(input);
      return{status:'research-converged',mach:input.mach,solverInput,flow:{x:Float64Array.of(grid),nodes:nodes(grid),...(config.subsonicEuler?{diagnostics:{maxMach:.6}}:{})},mesh:m,
        diagnostics:{iterations:1,equationResidual:1e-12,cells:grid}};
    },
    createStreamtubeBodySystem(input){return{layout:{bodies:input.bodies},baseGeometry:config.finite?[{}]:[],conditions:{lengthScale:1},
      adoptGeometry:x=>x.slice(),decode:()=>({nodes:nodes(input._grid),undisplacedNodes:nodes(input._grid)})};},
    transferStreamtubeGeometry(_source,state){return state.slice();},initialStreamtubeDisplacement(){return{};},
    initializeCoupledStreamtubeBody(input,options,controls){calls.initializers.push({grid:input._grid,ncrit:options.ncrit,controls:clone(controls)});
      if(config.initializationError)throw config.initializationError;
      if(config.secondInitializationError&&options.ncrit===4)throw config.secondInitializationError;
      if(config.failFineInitialization&&input._grid>16)throw new Error('cold fine BL failed');
      const prepared=makePrepared(input,options);
      if(config.ncritRecovery)prepared.initialization.thicknessFactor=controls.initialThicknessFactor;
      if(config.initializationReceipt)Object.assign(prepared.initialization,config.initializationReceipt);
      return prepared;},
    createCoupledStreamtubeBody(input,options){return makePrepared(input,options).system;},
    prepareCoupledCoarseProfile(args,options){calls.maps.push({sourceGrid:args.sourceResult.solverInput._grid,targetGrid:args.input._grid,tolerance:options.tolerance});
      if(config.sourceGeometryReplay) assert.equal(args.options.geometryReplay, config.sourceGeometryReplay);
      if(config.mappingError)throw config.mappingError;return makePrepared(args.input,args.options,true);},
    prepareCoupledGridSequence(args,options){calls.sequences.push({sourceGrid:args.sourceResult.solverInput._grid,
      targetGrid:args.input._grid,requestedGridIntervals:args.requestedGridIntervals,sourceGridIntervals:args.sourceGridIntervals,tolerance:options.tolerance});
      calls.profilePolicies.push(options.blPredictor);
      calls.sequenceConditions.push({sourceNcrit:args.sourceResult.checkpoint.restart.options.ncrit,
        sourceSystemNcrit:args.sourceSystem.conditions.ncrit,targetNcrit:args.options.ncrit,
        sourceHkPolicy:args.sourceResult.checkpoint.restart.options.hkFloorLinearization,
        sourceSystemHkPolicy:args.sourceSystem.conditions.hkFloorLinearization,targetHkPolicy:args.options.hkFloorLinearization});
      if(config.coarseRetainedNcrit!==undefined){
        assert.equal(args.sourceResult.checkpoint.restart.options.ncrit,args.options.ncrit,'Refinement must retain the active coarse criterion.');
        assert.equal(args.sourceSystem.conditions.ncrit,args.options.ncrit);
        assert.equal(args.sourceResult.checkpoint.restart.options.hkFloorLinearization,args.options.hkFloorLinearization);
        assert.equal(args.sourceSystem.conditions.hkFloorLinearization,args.options.hkFloorLinearization);
      }
      if(config.sourceGeometryReplay) assert.equal(args.options.geometryReplay, config.sourceGeometryReplay);
      if(config.mappingError)throw config.mappingError;
      const grid=config.multilevel?args.requestedGridIntervals+4:36;
      const input={...args.input,_grid:grid},prepared=makePrepared(input,args.options,true);
      return {...prepared,input,options:clone(args.options),initialEuler:{x:[grid],...flow(grid)},initialBL:[0,1,2,3],
        transfer:{method:'complete-state-grid-sequence',requestedGridRetained:false,
          ...(config.multilevel?{nominalGridIntervals:args.requestedGridIntervals}: {})}};},
    prepareCoupledGridLevels(args,options){
      calls.levels.push({requested:args.requestedGridIntervals,source:args.sourceGridIntervals});
      let retained=args.sourceResult,sourceSystem=args.sourceSystem,sourceGridIntervals=args.sourceGridIntervals;
      const levels=[];
      const finish=(prepared,failure)=>({prepared,sourceResult:retained,sourceSystem,sourceGridIntervals,
        diagnostics:{actualGridIntervals:sourceGridIntervals,requestedGridIntervals:args.requestedGridIntervals,
          finalSeedPrepared:prepared!==null,reachedTarget:false,stateConverged:true,levels,...(failure?{failure}: {})}});
      for(const [index,nominal] of [32,64,128].entries()){
        if(config.rejectLevel===nominal)return finish(null,{reason:`controlled ${nominal} transfer rejected`,stage:'preparation'});
        const prepared=stubs.prepareCoupledGridSequence({...args,sourceResult:retained,sourceSystem,
          sourceGridIntervals,requestedGridIntervals:nominal},options);
        const level={index,final:nominal===128,sourceGridIntervals,requestedGridIntervals:nominal,nominalGridIntervals:nominal};
        levels.push(level);if(level.final)return finish(prepared);
        const value=options.solveIntermediateLevel(prepared,level),candidate=value.result??value;
        if(value.failure||!candidate.converged)return finish(null,{reason:value.failure?.reason??candidate.reason,stage:'intermediate-solve'});
        retained=candidate;sourceGridIntervals=nominal;const saved=retained.checkpoint.restart;
        sourceSystem=stubs.createCoupledStreamtubeBody(saved.input,{...saved.options,initialEuler:saved.initialEuler,initialBL:saved.initialBL});
      }
      throw new Error('controlled multilevel fallthrough');
    },
    solveCoupledGridLevel(prepared,level,options){
      calls.levelSolves.push({level:clone(level),options:clone(prepared.options)});
      const labels={stage:'coupled-grid-refinement',startupAttempt:options.startupAttempt,gridLevel:level.nominalGridIntervals,
        requestedGridIntervals:level.targetGridIntervals};
      options.onStage?.(labels);
      let checkpoint;
      return stubs.solveCoupledStreamtubeIses(prepared.input,{...prepared.options,maxIterations:options.maxIterations,
        onIteration:h=>options.onIteration?.({...h,...labels}),
        onCheckpoint:(cp,details)=>{checkpoint=cp;options.onIterationCheckpoint?.(cp,{...details,...labels});},
        onMesh:state=>{
          options.onMesh?.(state.mesh,'solving',labels.stage);
          if(options.onFlow&&checkpoint)options.onFlow({checkpoint,flow:state.flow,
            bl:{stations:clone(prepared.system.bl.stations),surfaces:clone(prepared.system.bl.surfaces)},
            bodies:prepared.input.bodies,normalization:options.normalization,iteration:state.iteration,mach:prepared.input.mach,...labels});
        }});
    },
    solveCoupledStreamtubeIses(input,options){
      input??=options.resume.restart.input;const grid=input._grid;calls.solves.push({grid,maxIterations:options.maxIterations});
      calls.solveConditions.push({grid,ncrit:options.ncrit,hkFloorLinearization:options.hkFloorLinearization});
      if(config.earnedExtension){
        const offset=options.resume?.testCompletedIterations??0, controls=options.resume?.restart.options??options;
        const count=offset>=40&&!config.extensionFailed?Math.min(5,options.maxIterations):options.maxIterations;
        const baseOptions={transitionMode:'automatic',tripFractions:[[1,1]],ncrit:controls.ncrit,
          reynolds:controls.reynolds,transitionState:[2,3]};
        const continuation=clone(options.resume?.continuation??{fractions:[[0,.5,1]],lastRedistributedStagnation:[.2],
          iterationGeometry:'ises-sampled',stepAcceptance:'admissible',stagnationLimiter:'listing',
          linearOrdering:'station-auto',stationFallback:true,preferredOrdering:'amd',pivotTolerance:.001,
          projectionGeometry:'boundary-increment',blUpdate:'xfoil'});
        if(offset>=40)for(const name of ['iterationGeometry','stepAcceptance','stagnationLimiter'])
          assert.equal(options[name],continuation[name],'Extension preserves the exact resume controls.');
        const history=[],linearDiagnostics={solves:count,refinements:0,pivotRecoveries:0,maxRelativeResidual:1e-15,
          iterations:Array.from({length:count},(_,i)=>({iteration:i+1,pivotTolerance:continuation.pivotTolerance}))};
        let cp,x,families,currentFlow,currentMesh,converged=false;
        for(let k=0;k<=count;k++){
          const total=offset+k;
          converged=offset>=40&&!config.extensionFailed&&k===count;
          families=converged?{euler:1e-12,boundaryLayer:2e-12,edgeMatching:3e-12}
            :{euler:10/(total+1),boundaryLayer:30/(total+1),edgeMatching:20/(total+1)};
          x=Float64Array.of(grid,total,1,2,3);currentFlow=flow(grid);
          currentFlow.cells[0][0][0].interfacePressure.lower=total;
          cp={version:1,testCompletedIterations:total,families:clone(families),restart:{input,options:baseOptions,
            initialEuler:{x:x.slice(0,1),nodes:clone(currentFlow.nodes),undisplacedNodes:clone(currentFlow.undisplacedNodes)},initialBL:x.slice(1)},
            continuation:clone(continuation)};
          const h={iteration:k,step:k?1:0,residual:Math.max(...Object.values(families)),...families,
            backtracks:0,rejections:[],maintenance:{triggeredBodies:[],passages:[],dekinkRepairs:[]}};
          history.push(h);options.onIteration?.(h);
          options.onCheckpoint?.(cp,{history,linearDiagnostics,initialRedistribution:{accepted:true,resumed:offset>0}});
          currentMesh={...meshFor(grid,'extension'),nodes:clone(currentFlow.nodes),iteration:h};currentMesh.quality.invalidCells=[];
          options.onMesh?.({system:{layout:{bodies:input.bodies}},mesh:currentMesh,nodes:currentFlow.nodes,flow:currentFlow,iteration:h});
        }
        return{converged,reason:converged?'residual':'iteration limit',history,solverInput:input,coupledOptions:baseOptions,
          conditions:{transitionMode:'automatic',ncrit:baseOptions.ncrit,reynolds:baseOptions.reynolds,mach:input.mach},
          initialRedistribution:{accepted:true,resumed:offset>0},lastRejectedStep:null,flow:currentFlow,mesh:currentMesh,x,
          residual:Float64Array.of(families.euler,-families.boundaryLayer,families.edgeMatching,0,0),
          checkpoint:cp,families,boundaryLayer:{stations:[{id:0}]},linearDiagnostics};
      }
      if(config.stallRecovery){
        const offset=options.resume&&options.maxIterations?20:0;
        const controls=options.resume?.restart.options??options;
        const history=[],baseOptions={transitionMode:'automatic',tripFractions:[[1,1]],ncrit:controls.ncrit??4,
          reynolds:controls.reynolds??2.7e6,transitionState:[2,3]};
        const makeCp=iteration=>({version:1,
          families:options.maxIterations===0&&options.resume?clone(options.resume.families):{euler:.001,edgeMatching:.001,
            boundaryLayer:iteration===26?.01:iteration>26?.2:.3},
          restart:{input,options:baseOptions,initialEuler:{x:[grid],nodes:nodes(grid),undisplacedNodes:nodes(grid)},initialBL:[1,2,3,4]},
          continuation:{iterationGeometry:'ises-sampled',stepAcceptance:'admissible',stagnationLimiter:'listing',
            linearOrdering:'station-auto',projectionGeometry:'boundary-increment',blUpdate:'xfoil'}});
        let cp;
        for(let k=0;k<=options.maxIterations;k++){
          cp=makeCp(k+offset);const h={iteration:k,step:k?1:0,...cp.families,residual:cp.families.boundaryLayer};history.push(h);
          options.onIteration?.(h);options.onCheckpoint?.(cp,{history,initialRedistribution:{accepted:true}});
        }
        const value={converged:false,reason:'iteration limit',history,solverInput:input,coupledOptions:baseOptions,
          conditions:{transitionMode:'automatic',ncrit:baseOptions.ncrit},flow:flow(grid),mesh:meshFor(grid,'solved'),
          checkpoint:cp,families:cp.families,boundaryLayer:{stations:[{id:0}]},
          linearDiagnostics:{solves:options.maxIterations,refinements:0,pivotRecoveries:0,maxRelativeResidual:0,iterations:[]}};
        return value;
      }
      if(config.logRecovery){
        // Complete synthetic checkpoints exercise the real recovery planner and
        // assembly observers; no Euler, BL, mesh or Jacobian calculation runs.
        const controls=options.resume?.restart.options??options;
        const coordinate=options.resume?.continuation.shearCoordinate??options.shearCoordinate??'linear';
        const converged=coordinate==='logarithmic',count=converged?2:options.maxIterations;
        const cpInput={...input,bodies:input.bodies.map(body=>({...body,leadingIndex:1}))};
        const cpOptions={transitionMode:'automatic',ncrit:controls.ncrit,reynolds:controls.reynolds,
          tripFractions:[[1,1]],transitionState:[2,3],
          ...(controls.hkFloorLinearization?{hkFloorLinearization:controls.hkFloorLinearization}:{})};
        const continuation={fractions:[[0,1]],lastRedistributedStagnation:[.2],
          iterationGeometry:'ises-sampled',stepAcceptance:'admissible',stagnationLimiter:'listing',
          blUpdate:'xfoil',projectionGeometry:'boundary-increment',
          ...(coordinate==='logarithmic'?{shearCoordinate:coordinate}:{})};
        const x=Float64Array.of(grid,1,2,3,4),history=[];
        const families=converged?{euler:1e-12,boundaryLayer:2e-12,edgeMatching:3e-12}
          :{euler:.4,boundaryLayer:.2,edgeMatching:.1};
        const checkpoint={version:1,families:clone(families),restart:{input:cpInput,options:cpOptions,
          initialEuler:{x:x.slice(0,1),nodes:nodes(grid),undisplacedNodes:nodes(grid)},initialBL:x.slice(1)},continuation};
        const m=meshFor(grid,'log-recovery');m.quality.invalidCells=[];
        const linearDiagnostics={solves:count,refinements:0,pivotRecoveries:0,maxRelativeResidual:0,
          iterations:Array.from({length:count},(_,i)=>({iteration:i+1}))};
        for(let iteration=0;iteration<=count;iteration++){
          const h={iteration,step:iteration?1:0,residual:Math.max(...Object.values(families)),...families};
          history.push(h);options.onIteration?.(h);
          options.onCheckpoint?.(checkpoint,{history,linearDiagnostics,initialRedistribution:{accepted:true}});
        }
        return{converged,reason:converged?'residual':controls.ncrit===4?'iteration limit':'line search',
          history,solverInput:cpInput,coupledOptions:cpOptions,
          conditions:{transitionMode:'automatic',ncrit:cpOptions.ncrit,reynolds:cpOptions.reynolds,mach:cpInput.mach},
          flow:flow(grid),mesh:m,x,checkpoint,families,initialRedistribution:{accepted:true},
          residual:Float64Array.of(families.euler,families.boundaryLayer,families.edgeMatching,0,0),
          boundaryLayer:{stations:[{id:0}]},linearDiagnostics};
      }
      const m=meshFor(grid,'solved'),x=Float64Array.from([grid,0,1,2,3]),history=[{iteration:0,residual:1},{iteration:1,residual:1e-12}];
      options.onIteration?.(history[0]);if(!config.pairedFlow)options.onMesh?.({mesh:m,iteration:history[0]});
      const converged=config.ncritRecovery ? options.ncrit===4&&!config.sourceNcritFailed
        : grid===16?!config.coarseFailed:!config.fineFailed&&grid!==config.failLevelActualGrid;
      const f={input,options:{transitionMode:options.transitionMode,tripFractions:options.tripFractions,
        edgeMatching:options.edgeMatching,blThermodynamics:options.blThermodynamics,
        ...(config.sourceGeometryReplay?{geometryReplay:config.sourceGeometryReplay}:{}),
        ...(options.ncrit===undefined?{}:{ncrit:options.ncrit}),
        ...(options.hkFloorLinearization===undefined?{}:{hkFloorLinearization:options.hkFloorLinearization})},
        initialEuler:{x:x.slice(0,1),nodes:nodes(grid),undisplacedNodes:nodes(grid)},initialBL:x.slice(1)};
      const checkpoint={version:1,restart:f,continuation:{},families:{euler:1e-12,boundaryLayer:1e-12,edgeMatching:1e-12}};
      options.onCheckpoint?.(checkpoint,{kind:'accepted'});
      if(config.pairedFlow)options.onMesh?.({mesh:m,system:{layout:{bodies:input.bodies}},nodes:nodes(grid),flow:flow(grid),iteration:history[1]});
      return{converged,reason:converged?'residual':'line search',history,solverInput:input,coupledOptions:f.options,
        conditions:{ncrit:options.ncrit,...(options.hkFloorLinearization===undefined?{}:{hkFloorLinearization:options.hkFloorLinearization})},
        flow:flow(grid),mesh:m,x,checkpoint,families:checkpoint.families,
        boundaryLayer:{stations:[{id:0}]},linearDiagnostics:{solves:1,refinements:0,pivotRecoveries:0,maxRelativeResidual:0,iterations:[]}};
    },
    streamtubeMeshSnapshot(state){return config.earnedExtension
      ? {...(state.mesh??meshFor(32,'snapshot')),nodes:state.nodes,iteration:clone(state.iteration)}
      : state.mesh??meshFor(32,'snapshot');},
    transitionRecoveryPlan(result,options){
      if(config.stallRecovery&&options?.bestState&&result.history.at(-1).iteration===40){
        assert.equal(options.bestState.iteration,26);assert.equal(options.bestState.checkpoint.families.boundaryLayer,.01);
        return {reason:'boundary-layer-transition-stall',bestIteration:26};
      }
      return config.transitionRecovery&&result.solverInput._grid>16?{kind:'cycle'}:null;
    },
    recoverCoupledTransition(result,options){
      calls.recoveries.push({grid:result.solverInput._grid,gridSmoothing:clone(result.mesh.initialization.gridSmoothing)});
      if(config.stallRecovery){
        assert.equal(result.families.boundaryLayer,.01,'Only the accepted best complete state may seed recovery.');
        assert.equal(options.maxIterations,20);assert.equal(options.tolerance,1e-10);
        assert.equal(result.checkpoint.continuation.linearOrdering,'station-auto');
        if(config.recoveryError)throw new Error('controlled refinement failure');
        options.onIterationCheckpoint?.(result.checkpoint,{kind:'refined'});
        options.onFlow?.({checkpoint:result.checkpoint,bl:{stations:[{id:44}],surfaces:[]},stage:'coupled',refined:true,iteration:{iteration:0}});
        if(config.recoveryUnconverged)return {...result,history:[{iteration:0},{iteration:1}],reason:'refined iteration limit'};
        return {...result,converged:true,reason:'residual',automaticRefinement:{kind:'transition-local'},
          history:Array.from({length:8},(_,iteration)=>({iteration}))};
      }
      return {...result,automaticRefinement:{kind:'transition-recovery'}};
    },
  };
  const key=`__coarseReview${++sequence}`;globalThis[key]=stubs;
  const rewritten=source.replace(/import \{([^}]+)\} from '[^']+';/g,(_,names)=>`const {${names.replace(/\s+as\s+/g, ": ")}} = globalThis[${JSON.stringify(key)}];`);
  const module=await import('data:text/javascript;base64,'+Buffer.from(rewritten+'\n//# sourceURL=coarse-startup-routing-'+key+'.js').toString('base64'));
  return{...module,calls,packetFor:input=>({prepared:geometry(input),caseData:clone(input)}),release:()=>delete globalThis[key]};
}

test('planner preserves all physical controls, detaches input, and excludes legacy paths',async()=>{
  const h=await harness(),input=baseCase(),before=clone(input);
  const prepared={initialization:{thicknessFactor:.25},system:{bl:{transitionMode:'automatic',hasFiniteBase:false,trips:[[1,1]]}}};
  const target=h.coupledCoarseStartupCase(input,prepared,{maxIterations:20});assert.deepEqual(target,{...before,gridIntervals:16,gridTubes:7});
  target.elements[0].points[0].x=-1;target.materialTrips[0][0]=.5;assert.deepEqual(input,before);
  assert.deepEqual(h.coupledCoarseStartupCase(input,{...prepared,initialization:{thicknessFactor:1}},{maxIterations:20}),
    {...before,gridIntervals:16,gridTubes:7});
  for(const [c,p,o] of [[{...input,gridIntervals:16},prepared,{}],
    [input,{...prepared,system:{bl:{...prepared.system.bl,transitionMode:'fixed-trip'}}},{}],
    [input,{...prepared,system:{bl:{...prepared.system.bl,trips:[[.05,1]]}}},{}],
    [input,prepared,{enabled:false}],[input,prepared,{maxIterations:0}]])
    assert.equal(h.coupledCoarseStartupCase(c,p,{maxIterations:20,...o}),null);
  h.release();
});
test('finite-base sequencing bypasses a failing unused cold fine BL guess and publishes its actual topology',async()=>{
  const h=await harness({finite:true,failFineInitialization:true}),input=finiteCase(),before=clone(input),checkpoints=[];
  const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,onCheckpoint:c=>checkpoints.push(c)});
  assert.deepEqual(h.calls.solves.map(c=>c.grid),[16,36]);
  assert.equal(h.calls.maps.length,0);
  assert.equal(h.calls.levels.length,0);assert.equal(h.calls.levelSolves.length,0);
  assert.deepEqual(h.calls.profilePolicies,['interpolate']);
  assert.deepEqual(h.calls.sequences,[{sourceGrid:16,targetGrid:32,requestedGridIntervals:32,sourceGridIntervals:16,tolerance:1e-10}]);
  assert.equal(result.converged,true);assert.equal(result.restart.input._grid,36);
  assert.equal(result.mesh.cells.length,36);
  assert.equal(result.initialization.coarseInitialization.requestedGridRetained,false);
  assert.equal(result.initialization.coarseInitialization.directFineBLInitialization,'deferred');
  assert.deepEqual(h.calls.initializers.map(c=>c.grid),[16]);
  assert.equal(result.automaticRefinement.kind,'coarse-to-fine');
  assert.equal(result.mesh.initialization.gridRefinement.kind,'coarse-to-fine');
  assert.ok(checkpoints.filter(c=>c.restart).every(c=>c.restart.input._grid===36));
  assert.deepEqual(input,before);h.release();
});
test('sharp profiles also bypass unused fine gas and BL guesses',async()=>{
  const h=await harness({fineGasError:new Error('unused fine gas'),failFineInitialization:true});
  const input=baseCase(),before=clone(input),meshes=[],iterations=[];
  const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,
    onMesh:(mesh,phase,stage)=>meshes.push({grid:mesh.cells.length,stage}),onIteration:h=>iterations.push(h)});
  assert.deepEqual(h.calls.euler.map(i=>[i.gridIntervals,i.gridTubes]),[[16,7]]);
  assert.deepEqual(h.calls.initializers.map(c=>c.grid),[16]);
  assert.deepEqual(h.calls.solves.map(c=>c.grid),[16,36]);
  assert.equal(h.calls.maps.length,0);assert.equal(h.calls.sequences.length,1);
  assert.deepEqual(h.calls.profilePolicies,['xfoil-mrchdu']);
  assert.equal(result.converged,true);assert.equal(result.mesh.cells.length,36);
  assert.equal(result.initialization.euler.requestedFineEulerPerformed,false);
  assert.ok(meshes.filter(m=>m.grid===16).every(m=>m.stage==='coupled-coarse-initialization'));
  assert.ok(iterations.some(i=>i.coarseStage==='euler'));
  assert.deepEqual(input,before);h.release();
});
test('rejected refinement retains the converged parent without a cold fine retry',async()=>{
  for(const finite of [false,true]){
    const h=await harness({finite,mappingError:new Error('refined grid domain'),fineGasError:new Error('must not run')});
    const result=h.solveCoupledStreamtubeAssembly(finite?finiteCase():baseCase(),{maxIterations:3});
    assert.deepEqual(h.calls.euler.map(c=>c.gridIntervals),[16]);
    assert.deepEqual(h.calls.solves.map(c=>c.grid),[16]);
    assert.equal(result.converged,false);assert.equal(result.stateConverged,true);
    assert.equal(result.restart.input._grid,16);assert.equal(result.mesh.cells.length,16);
    assert.equal(result.gridSequence.reachedTarget,false);
    assert.equal(result.gridSequence.requestedGridIntervals,32);
    assert.match(result.reason,/retained converged 16-interval level/);
    assert.match(result.reason,/refined grid domain/);h.release();
  }
});
test('prepared geometry belongs to requested preview; coarse precursor receives no unrelated packet',async()=>{
  const h=await harness(),input=baseCase(),packet=h.packetFor(input),handoffs=[];
  try {
    h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,preparedEuler:packet,onEulerPrepared:p=>handoffs.push(p)});
    assert.equal(h.calls.eulerBindings.length,1);
    assert.deepEqual(h.calls.eulerBindings[0],{grid:16,preparedEuler:undefined,onEulerPrepared:undefined});
    assert.deepEqual(h.calls.restores,[32]);assert.equal(handoffs.length,1);
    assert.equal(handoffs[0].caseData.gridIntervals,32);
  } finally {h.release();}
});
test('a failed coarse solve can still try the original fine startup once',async()=>{
  const h=await harness({coarseFailed:true}),r=h.solveCoupledStreamtubeAssembly(baseCase(),{maxIterations:3,maxStartupAttempts:2});
  assert.deepEqual(h.calls.euler.map(i=>i.gridIntervals),[16,32]);assert.deepEqual(h.calls.solves.map(c=>c.grid),[16,16,32]);
  assert.equal(h.calls.sequences.length,0);assert.equal(r.initialization.coarseInitialization.accepted,false);assert.equal(r.mesh.cells.length,32);h.release();
});
test('coarse stage/iteration/mesh Error cancellation aborts before mapping or fine solve',async()=>{
  for(const hook of ['onStage','onIteration','onMesh']){
    const h=await harness(),cancel=new Error(`cancel ${hook}`);
    const callbacks={onStage:e=>{if(e.stage==='coupled-coarse-initialization')throw cancel;},
      onIteration:e=>{if(e.stage==='coupled-coarse-initialization')throw cancel;},onMesh:(m,p,s)=>{if(s==='coupled-coarse-initialization')throw cancel;}};
    assert.throws(()=>h.solveCoupledStreamtubeAssembly(baseCase(),{maxIterations:3,[hook]:callbacks[hook]}),e=>e===cancel);
    assert.equal(h.calls.maps.length,0);assert.equal(h.calls.solves.some(c=>c.grid===32),false);h.release();
  }
});
test('non-Error coarse cancellation preserves the thrown value',async()=>{
  for(const cancel of [undefined,null,Object.freeze(new Error('frozen cancellation'))]){
    const h=await harness();let thrown=false,value='sentinel';
    try{h.solveCoupledStreamtubeAssembly(baseCase(),{maxIterations:3,onStage:e=>{if(e.stage==='coupled-coarse-initialization')throw cancel;}});}
    catch(error){thrown=true;value=error;}assert.equal(thrown,true);assert.equal(value,cancel);h.release();
  }
});
test('callback cancellation survives lower-layer smoothing error translation',async()=>{
  const h=await harness({translateSmoothingObserver:true}),cancel=new Error('stop during coarse smoothing');
  assert.throws(()=>h.solveCoupledStreamtubeAssembly(baseCase(),{maxIterations:3,onMesh:(mesh,phase,stage)=>{
    if(stage==='coupled-coarse-initialization'&&phase==='smoothing')throw cancel;
  }}),error=>error===cancel);assert.equal(h.calls.maps.length,0);h.release();
});
test('legacy16/fixed/disabled paths preserve baseline stub outputs exactly',async()=>{
  for(const [casePatch,config,options] of [[{gridIntervals:16},{},{}],
    [{transitionMode:'fixed-trip',materialTrips:[[.05,.05]]},{},{}],[{}, {},{coarseStartup:false}]]){
    const a=await harness(config),b=await harness(config,baseline),input={...baseCase(),...casePatch};
    const ra=a.solveCoupledStreamtubeAssembly(input,{maxIterations:3,...options});
    const rb=b.solveCoupledStreamtubeAssembly(input,{maxIterations:3,...options});
    assert.equal(ra.solverInput.wakeDisplacementMotion,'te-center');
    assert.deepEqual(ra,withFreshWakeChart(rb));assert.deepEqual(a.calls,withFreshWakeChart(b.calls));a.release();b.release();
  }
});

test('finite-base fixed trips, disabled sequencing and zero iterations retain their original fine precursor route',async()=>{
  for(const [patch,options] of [
    [{transitionMode:'fixed-trip',materialTrips:[[.05,.05]]},{}],
    [{},{coarseStartup:false}],
    [{},{maxIterations:0}],
  ]){
    const current=await harness({finite:true}),old=await harness({finite:true},baseline);
    const input={...finiteCase(),...patch},controls={maxIterations:3,...options};
    try{
      const actual=current.solveCoupledStreamtubeAssembly(input,controls);
      assert.equal(actual.solverInput.wakeDisplacementMotion,'te-center');
      assert.deepEqual(actual,withFreshWakeChart(old.solveCoupledStreamtubeAssembly(input,controls)));
      assert.deepEqual(current.calls,withFreshWakeChart(old.calls));
      assert.equal(current.calls.geometry.length,0);
      assert.deepEqual(current.calls.euler.map(i=>i.gridIntervals),[32]);
      assert.equal(current.calls.sequences.length,0);
    }finally{current.release();old.release();}
  }
});

test('public exhausted BL-stall recovery selects the best accepted checkpoint across resumed chunks',async()=>{
  const h=await harness({stallRecovery:true}),frames=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,eulerIsmom:4},
      {maxIterations:40,maxStartupAttempts:1,onFlow:frame=>frames.push(frame)});
    assert.equal(result.converged,true);assert.equal(result.automaticRefinement.kind,'transition-local');
    assert.deepEqual(h.calls.solves.map(c=>c.maxIterations),[20,20,0]);
    assert.equal(h.calls.recoveries.length,1);assert.equal(result.initialization.attempts[0].iterations,47);
    assert.equal(frames.at(-1).refined,true);assert.deepEqual(frames.at(-1).bl.stations,[{id:44}]);
    assert.equal(result.solverInput.hybrid.ismom,4);
  }finally{h.release();}
});

test('public BL-stall recovery failure retains the best complete state and publishes matching coarse flow',async()=>{
  for(const config of [{recoveryError:true},{recoveryUnconverged:true}]){
    const h=await harness({stallRecovery:true,...config}),frames=[],meshes=[];
    try{
      const result=h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,eulerIsmom:4},
        {maxIterations:40,maxStartupAttempts:1,onFlow:frame=>frames.push(frame),onMesh:mesh=>meshes.push(mesh)});
      assert.equal(result.converged,false);assert.equal(result.families.boundaryLayer,.01);
      assert.equal(result.recoveryAttempt.retained,'best accepted unrefined state');
      assert.equal(result.automaticRefinement,undefined);assert.match(result.reason,/Retained the best accepted state/);
      assert.equal(frames.at(-1).checkpoint.families.boundaryLayer,.01);
      assert.equal(frames.at(-1).bl.stations[0].id,16,'Abandoned refined topology must not describe the retained flow.');
      assert.equal(h.calls.recoveries.length,1);
    }finally{h.release();}
  }
});

test('public BL-stall recovery propagates cancellation even when the observer throws undefined',async()=>{
  const h=await harness({stallRecovery:true});let caught=false;
  try{
    try{h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16},
      {maxIterations:40,onFlow:frame=>{if(frame.refined)throw undefined;}});}
    catch(error){caught=true;assert.equal(error,undefined);}
    assert.equal(caught,true);assert.equal(h.calls.initializers.length,1);
  }finally{h.release();}
});

test('zero or short public update budgets do not trigger an exhausted-history recovery',async()=>{
  for(const maxIterations of [0,3,20]){
    const h=await harness({stallRecovery:true});
    try{
      h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16},{maxIterations,maxStartupAttempts:1});
      assert.equal(h.calls.recoveries.length,0);
      assert.deepEqual(h.calls.solves.map(c=>c.maxIterations),[maxIterations]);
    }finally{h.release();}
  }
});


test('complete finite-base sequencing does not require a valid fine gas guess or fine Euler solve',async()=>{
  const gas=Object.assign(new Error('fine interface pressure is negative'),{code:'streamtube-interface-pressure',stage:'gas-initialization'});
  const h=await harness({finite:true,fineGasError:gas,failFineInitialization:true});
  const input={...finiteCase(),mach:.185,alpha:6,reynolds:2.7e6,ncrit:9,eulerIsmom:4,gridEllipticSmoothing:false};
  const before=clone(input),handoffs=[],meshes=[],stages=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,onEulerPrepared:p=>handoffs.push(p),
      onMesh:(mesh,phase,stage)=>meshes.push({grid:mesh.cells.length,stage}),onStage:s=>stages.push(s)});
    assert.deepEqual(h.calls.euler.map(i=>i.gridIntervals),[16]);
    assert.deepEqual(h.calls.geometry.map(i=>i.gridIntervals),[32]);
    assert.deepEqual(h.calls.initializers.map(i=>i.grid),[16]);
    assert.deepEqual(h.calls.solves.map(i=>i.grid),[16,36]);
    assert.deepEqual({...h.calls.euler[0],gridIntervals:32,gridTubes:11,ncrit:9,coupledNativeHk:undefined},
      {...before,coupledNativeHk:undefined});
    assert.equal(h.calls.euler[0].ncrit,4);
    assert.equal(h.calls.euler[0].coupledNativeHk,true);
    assert.equal(handoffs.length,1);assert.equal(handoffs[0].caseData.gridIntervals,32);
    assert.equal(h.calls.eulerBindings[0].preparedEuler,undefined);assert.equal(h.calls.eulerBindings[0].onEulerPrepared,undefined);
    assert.equal(result.mach,.185);assert.equal(result.solverSettings.eulerIsmom,4);
    assert.equal(result.solverSettings.hybrid.ismom,4);assert.equal(result.solverSettings.blThermodynamics,'historical-common-isentrope');
    assert.equal(result.initialization.euler.source,'coarse-grid');assert.equal(result.initialization.euler.requestedFineEulerPerformed,false);
    assert.equal(result.initialization.coarseInitialization.directFineEulerInitialization,'skipped');
    assert.equal(result.initialization.coarseInitialization.requestedGridRetained,false);
    assert.ok(meshes.filter(m=>m.grid===16).every(m=>m.stage==='coupled-coarse-initialization'));
    assert.equal(meshes.at(-1).grid,36);assert.deepEqual(input,before);
  }finally{h.release();}
});

test('retained coarse Ncrit7.5 maps unchanged into the fine root before continuation reaches the requested Ncrit9',async()=>{
  const h=await harness({finite:true,coarseRetainedNcrit:7.5,
    fineGasError:new Error('unused requested-grid Euler must not run'),failFineInitialization:true});
  const input={...finiteCase(),mach:.185,alpha:6,reynolds:2.7e6,ncrit:9,eulerIsmom:4,gridTubes:9},before=clone(input);
  const stages=[],iterations=[],checkpoints=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,maxStartupAttempts:2,
      onStage:e=>stages.push(e),onIteration:e=>iterations.push(e),onCheckpoint:e=>checkpoints.push(e)});
    assert.deepEqual(input,before);
    assert.deepEqual(h.calls.euler.map(c=>[c.gridIntervals,c.ncrit]),[[16,4]]);
    assert.deepEqual(h.calls.initializers.map(c=>c.grid),[16]);
    assert.deepEqual(h.calls.solveConditions,[{grid:16,ncrit:4,hkFloorLinearization:'native'},
      {grid:36,ncrit:7.5,hkFloorLinearization:'native'}]);
    assert.deepEqual(h.calls.sequenceConditions,[{sourceNcrit:7.5,sourceSystemNcrit:7.5,targetNcrit:7.5,
      sourceHkPolicy:'native',sourceSystemHkPolicy:'native',targetHkPolicy:'native'}]);
    assert.deepEqual(h.calls.continuations,[
      {phase:'coarse',sourceNcrit:4,targetNcrit:9,sourceConverged:true,sourceGrid:16,sourceHkPolicy:'native'},
      {phase:'fine',sourceNcrit:7.5,targetNcrit:9,sourceConverged:true,sourceGrid:36,sourceHkPolicy:'native'},
    ]);
    assert.equal(result.converged,true);assert.equal(result.conditions.ncrit,9);
    assert.equal(result.actualNcrit,9);assert.equal(result.targetNcrit,9);
    assert.equal(result.ncritContinuation.reachedTarget,true);assert.equal(result.ncritContinuation.stateConverged,true);
    assert.equal(result.ncritContinuation.coarse.actualNcrit,7.5);assert.equal(result.ncritContinuation.coarse.reachedTarget,false);
    assert.equal(result.ncritContinuation.fine.actualNcrit,9);assert.equal(result.ncritContinuation.fine.reachedTarget,true);
    assert.equal(result.restart.options.ncrit,9);assert.equal(result.restart.options.hkFloorLinearization,'native');
    assert.equal(result.restart.input._grid,36);assert.equal(result.mesh.cells.length,36);
    assert.equal(result.initialization.euler.requestedFineEulerPerformed,false);
    assert.equal(result.initialization.coarseInitialization.directFineEulerInitialization,'skipped');
    assert.equal(result.initialization.attempts.length,1);
    const fineInitial=checkpoints.find(e=>e.restart);
    assert.equal(fineInitial.restart.options.ncrit,7.5);assert.equal(checkpoints.at(-1).restart.options.ncrit,9);
    assert.ok(iterations.some(e=>e.stage==='coupled'&&e.actualNcrit===7.5&&e.targetNcrit===9));
    assert.ok(stages.some(e=>e.stage==='coupled'&&e.actualNcrit===9&&e.targetNcrit===9));
  }finally{h.release();}
});

test('an unreached Ncrit9 retains the converged Ncrit8 fine state but marks the overall requested solve unconverged',async()=>{
  const h=await harness({finite:true,coarseRetainedNcrit:7.5,fineRetainedNcrit:8,
    fineGasError:new Error('unused requested-grid Euler must not run'),failFineInitialization:true});
  const input={...finiteCase(),mach:.185,alpha:6,reynolds:2.7e6,ncrit:9,eulerIsmom:4,gridTubes:9},before=clone(input);
  const checkpoints=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,maxStartupAttempts:2,onCheckpoint:e=>checkpoints.push(e)});
    assert.deepEqual(input,before);assert.deepEqual(h.calls.euler.map(c=>c.gridIntervals),[16]);
    assert.deepEqual(h.calls.solves.map(c=>c.grid),[16,36]);assert.deepEqual(h.calls.initializers.map(c=>c.grid),[16]);
    assert.deepEqual(h.calls.continuations.map(c=>[c.phase,c.sourceNcrit,c.targetNcrit]),[['coarse',4,9],['fine',7.5,9]]);
    assert.equal(result.converged,false);assert.equal(result.stateConverged,true);assert.equal(result.status,'unconverged');
    assert.equal(result.conditions.ncrit,8);assert.equal(result.actualNcrit,8);assert.equal(result.targetNcrit,9);
    assert.equal(result.ncritContinuation.reachedTarget,false);assert.equal(result.ncritContinuation.stateConverged,true);
    assert.equal(result.ncritContinuation.fine.actualNcrit,8);assert.equal(result.ncritContinuation.fine.reachedTarget,false);
    assert.match(result.reason,/Requested Ncrit 9 was not reached; retained Ncrit 8/);
    assert.equal(result.restart.options.ncrit,8);assert.equal(result.checkpoint.restart.options.ncrit,8);
    assert.equal(result.coupledOptions.ncrit,8);assert.equal(result.restart.options.hkFloorLinearization,'native');
    assert.equal(result.restart.input._grid,36);assert.equal(result.mesh.cells.length,36);
    assert.deepEqual(result.restart.initialEuler.nodes,result.flow.nodes);
    assert.equal(result.initialization.attempts.length,1);assert.equal(result.initialization.attempts[0].converged,false);
    assert.equal(result.initialization.euler.requestedFineEulerPerformed,false);
    assert.equal(result.initialization.coarseInitialization.directFineEulerInitialization,'skipped');
    assert.equal(checkpoints.at(-1).restart.options.ncrit,8);
  }finally{h.release();}
});

test('final grid-solve metadata follows the retained Ncrit root while earlier continuation provenance remains unchanged',async()=>{
  for(const actualNcrit of [9,8]){
    const h=await harness({finite:true,coarseRetainedNcrit:7.5,fineRetainedNcrit:actualNcrit,distinctContinuationRoot:true});
    try{
      const result=h.solveCoupledStreamtubeAssembly({...finiteCase(),ncrit:9,eulerIsmom:4},{maxIterations:3});
      const sequence=result.gridSequence,final=sequence.finalSolve;
      assert.equal(result.converged,actualNcrit===9);
      assert.equal(sequence.reachedTarget,true,'The retained root reached the grid even if requested Ncrit is unfinished.');
      assert.equal(final.converged,true);assert.equal(final.actualNcrit,actualNcrit);
      assert.equal(final.reason,'continued root residual');assert.equal(final.iterations,4);
      assert.equal(final.cells,result.mesh.cells.length);assert.deepEqual(final.families,result.checkpoint.families);
      assert.deepEqual(final.families,result.families);
      assert.equal(result.automaticRefinement,sequence);assert.equal(result.mesh.initialization.gridRefinement,sequence);
      const earlier=h.calls.continuationGridFrames[0];
      assert.notEqual(earlier,sequence);assert.equal(earlier.reachedTarget,true);
      assert.equal(earlier.finalSolve.actualNcrit,7.5);assert.equal(earlier.finalSolve.iterations,1);
      assert.equal(earlier.finalSolve.reason,'residual');
      assert.deepEqual(earlier.finalSolve.families,{euler:1e-12,boundaryLayer:1e-12,edgeMatching:1e-12});
    }finally{h.release();}
  }
});

test('sequenced flow can enter transition recovery without referring to an unperformed fine Euler solve',async()=>{
  const h=await harness({finite:true,fineGasError:new Error('unused fine gas'),transitionRecovery:true});
  try{
    const result=h.solveCoupledStreamtubeAssembly(finiteCase(),{maxIterations:3});
    assert.deepEqual(h.calls.euler.map(i=>i.gridIntervals),[16]);
    assert.deepEqual(h.calls.recoveries,[{grid:36,gridSmoothing:{converged:true}}]);
    assert.equal(result.initialization.euler.source,'coarse-grid');
    assert.equal(result.initialization.euler.requestedFineEulerPerformed,false);
    assert.equal(result.automaticRefinement.kind,'transition-recovery');
    assert.equal(result.initialization.attempts.length,1);
  }finally{h.release();}
});

test('failed coarse then failed fine gas restores requested preview and preserves typed error plus coarse reason',async()=>{
  const coarse=Object.assign(new Error('coarse gas rejected'),{stage:'gas-initialization'});
  const fine=Object.assign(new Error('fine gas rejected'),{stage:'gas-initialization',code:'pressure',diagnostics:{group:0,i:167,tube:2}});
  const h=await harness({finite:true,coarseGasError:coarse,fineGasError:fine}),meshes=[];
  try{
    assert.throws(()=>h.solveCoupledStreamtubeAssembly(finiteCase(),{maxIterations:3,onMesh:m=>meshes.push(m)}),error=>error===fine);
    assert.deepEqual(h.calls.euler.map(i=>i.gridIntervals),[16,32]);assert.equal(h.calls.geometry.length,1);
    assert.equal(meshes.at(-1).cells.length,32);assert.equal(fine.stage,'gas-initialization');
    assert.deepEqual(fine.diagnostics,{group:0,i:167,tube:2});
    assert.equal(fine.initialization.coarseInitialization.reason,'coarse gas rejected');
    assert.equal(h.calls.solves.length,0);
  }finally{h.release();}
});

test('early geometry handoff and coarse observers preserve cancellation without a fine gas fallback',async()=>{
  for(const hook of ['onEulerPrepared','onStage','onIteration','onMesh'])for(const cancel of [new Error('cancel'),undefined,Object.freeze({cancelled:true})]){
    const h=await harness({finite:true}),callbacks={onEulerPrepared:()=>{throw cancel;},
      onStage:e=>{if(e.stage==='coupled-coarse-initialization')throw cancel;},
      onIteration:e=>{if(e.stage==='coupled-coarse-initialization')throw cancel;},
      onMesh:(_m,_p,stage)=>{if(stage==='coupled-coarse-initialization')throw cancel;}};
    try{
      let caught=false;try{h.solveCoupledStreamtubeAssembly(finiteCase(),{maxIterations:3,[hook]:callbacks[hook]});}
      catch(error){caught=true;assert.equal(error,cancel);}
      assert.equal(caught,true);assert.equal(h.calls.euler.some(i=>i.gridIntervals===32),false);
      assert.equal(h.calls.sequences.length,0);
    }finally{h.release();}
  }
});

test('128-interval sharp startup solves accepted intermediate levels and preserves every requested physical control',async()=>{
  const h=await harness({multilevel:true,pairedFlow:true,failFineInitialization:true,fineGasError:new Error('unused cold fine gas')});
  const input={...baseCase(),gridIntervals:128,gridTubes:7,reynolds:2.51e6,eulerIsmom:4},before=clone(input),flows=[],meshes=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,
      onMesh:(mesh,phase,stage)=>meshes.push({mesh:clone(mesh),phase,stage}),onFlow:flow=>flows.push(clone(flow))});
    assert.deepEqual(h.calls.geometry.map(c=>c.gridIntervals),[128]);
    assert.deepEqual(h.calls.euler.map(c=>c.gridIntervals),[16]);
    assert.deepEqual(h.calls.initializers.map(c=>c.grid),[16]);
    assert.deepEqual(h.calls.solves.map(c=>c.grid),[16,36,68,132]);
    assert.deepEqual(h.calls.levels,[{requested:128,source:16}]);
    assert.deepEqual(h.calls.levelSolves.map(c=>c.level.nominalGridIntervals),[32,64]);
    assert.deepEqual(h.calls.sequences.map(c=>[c.sourceGridIntervals,c.requestedGridIntervals]),[[16,32],[32,64],[64,128]]);
    assert.deepEqual(h.calls.profilePolicies,['xfoil-mrchdu','xfoil-mrchdu','xfoil-mrchdu']);
    for(const call of h.calls.euler)assert.deepEqual({...call,gridIntervals:128,gridTubes:7},before);
    for(const call of h.calls.levelSolves){
      assert.equal(call.options.reynolds,2.51e6);assert.equal(call.options.ncrit,4);
      assert.equal(call.level.targetGridIntervals,128);assert.deepEqual(call.options.tripFractions,[[1,1]]);
    }
    assert.equal(result.converged,true);assert.equal(result.restart.input._grid,132);assert.equal(result.mesh.cells.length,132);
    assert.equal(result.gridSequence.parentUnknowns,73,'The last64-level parent owns the final128 transfer.');
    assert.equal(result.gridSequence.initialization.nominalGridIntervals,128);
    assert.equal(result.initialization.coarseInitialization.directFineEulerInitialization,'skipped');
    assert.equal(result.solverSettings.eulerIsmom,4);assert.equal(result.restart.input.hybrid.ismom,4);
    const intermediate=flows.filter(f=>f.stage==='coupled-grid-refinement');
    assert.deepEqual(intermediate.map(f=>[f.gridLevel,f.checkpoint.restart.input._grid]),[[32,36],[64,68]]);
    for(const f of flows){
      assert.deepEqual(f.flow.nodes,f.checkpoint.restart.initialEuler.nodes);
      assert.equal(f.bl.stations[0].i,f.checkpoint.restart.input._grid);
    }
    assert.equal(meshes.at(-1).mesh.cells.length,132);assert.equal(flows.at(-1).checkpoint.restart.input._grid,132);
    assert.deepEqual(input,before);
  }finally{h.release();}
});

test('failed128 hierarchy retains its last converged resolution with matching Cp frame and never retries cold fine gas',async()=>{
  for(const config of [{rejectLevel:64,retained:32,actual:36},{rejectLevel:128,retained:64,actual:68},
    {failLevelActualGrid:68,retained:32,actual:36}]){
    const h=await harness({multilevel:true,pairedFlow:true,fineGasError:new Error('unused cold fine gas'),...config});
    const input={...baseCase(),gridIntervals:128,gridTubes:7,eulerIsmom:4},flows=[],meshes=[],checkpoints=[];
    try{
      const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,
        onMesh:mesh=>meshes.push(clone(mesh)),onFlow:flow=>flows.push(clone(flow)),onCheckpoint:cp=>checkpoints.push(clone(cp))});
      assert.equal(result.converged,false);assert.equal(result.stateConverged,true);
      assert.equal(result.gridSequence.reachedTarget,false);assert.equal(result.gridSequence.actualGridIntervals,config.retained);
      assert.equal(result.gridSequence.requestedGridIntervals,128);
      assert.ok(result.initialization.euler, 'Retained intermediate levels preserve Euler startup provenance.');
      assert.equal(result.restart.input._grid,config.actual);assert.equal(result.mesh.cells.length,config.actual);
      assert.match(result.reason,new RegExp(`retained converged ${config.retained}-interval level`));
      assert.deepEqual(h.calls.euler.map(c=>c.gridIntervals),[16]);assert.deepEqual(h.calls.initializers.map(c=>c.grid),[16]);
      assert.equal(h.calls.solves.some(c=>c.grid===132||c.grid===128),false);
      assert.equal(meshes.at(-1).cells.length,config.actual);
      const cp=flows.at(-1).checkpoint;
      assert.equal(cp.restart.input._grid,config.actual);
      assert.deepEqual(flows.at(-1).flow.nodes,cp.restart.initialEuler.nodes);
      assert.equal(flows.at(-1).bl.stations[0].i,config.actual);
      assert.equal(checkpoints.at(-1).restart.input._grid,config.actual);
    }finally{h.release();}
  }
});

test('128 hierarchy observer cancellation reaches the caller unchanged before a cold fallback or later level',async()=>{
  for(const hook of ['onStage','onIteration','onMesh','onFlow','onIterationCheckpoint'])
    for(const stop of [undefined,'stop',Object.freeze({cancelled:true})]){
      const h=await harness({multilevel:true,pairedFlow:true});
      const callbacks={
        onStage:e=>{if(e.stage==='coupled-grid-refinement')throw stop;},
        onIteration:e=>{if(e.stage==='coupled-grid-refinement')throw stop;},
        onMesh:(_m,_p,stage)=>{if(stage==='coupled-grid-refinement')throw stop;},
        onFlow:e=>{if(e.stage==='coupled-grid-refinement')throw stop;},
        onIterationCheckpoint:(_cp,e)=>{if(e.stage==='coupled-grid-refinement')throw stop;},
      };
      try{
        assert.throws(()=>h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:128},{maxIterations:3,[hook]:callbacks[hook]}),value=>Object.is(value,stop));
        assert.deepEqual(h.calls.euler.map(c=>c.gridIntervals),[16]);
        assert.ok(h.calls.solves.every(c=>c.grid===16||c.grid===36));
        assert.ok(h.calls.sequences.length<=1);
      }finally{h.release();}
    }
});

test('failed default automatic startup reuses its Euler state for one thin Ncrit source and returns the requested refined root',async()=>{
  const h=await harness({ncritRecovery:true,pairedFlow:true}),input={...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4};
  const before=clone(input),frames=[],iterations=[],checkpoints=[],stages=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:3,onFlow:f=>frames.push(f),
      onIteration:h=>iterations.push(h),onCheckpoint:c=>checkpoints.push(c),onStage:s=>stages.push(s)});
    assert.deepEqual(input,before);assert.equal(h.calls.euler.length,1);
    assert.deepEqual(h.calls.initializers.map(c=>[c.ncrit,c.controls.initialThicknessFactor]),[[9,1],[4,.25]]);
    assert.deepEqual(h.calls.solveConditions.map(c=>[c.ncrit,c.hkFloorLinearization]),[[9,undefined],[4,'native']]);
    assert.equal(h.calls.ncritRecoveries.length,1);assert.equal(h.calls.ncritRecoveries[0].plan.coarseAdvanceNcrit,5.5);
    assert.equal(h.calls.ncritRecoveries[0].plan.shearCoordinate,'linear');
    assert.equal(h.calls.ncritRecoveries[0].plan.hkFloorLinearization,'native');
    for(const [attempt,hk] of [[1,'exact'],[2,'native']]){
      const startup=stages.find(s=>s.stage==='coupled'&&s.startupAttempt===attempt);
      assert.equal(startup.shearCoordinate,'linear');assert.equal(startup.hkFloorLinearization,hk);
      const update=iterations.find(h=>h.stage==='coupled'&&h.startupAttempt===attempt);
      assert.equal(update.shearCoordinate,'linear');assert.equal(update.hkFloorLinearization,hk);
    }
    assert.equal(result.converged,true);assert.equal(result.actualNcrit,9);assert.equal(result.ncritContinuation.reachedTarget,true);
    assert.equal(result.checkpoint.restart.options.ncrit,9);assert.equal(result.initialization.attempts.length,2);
    assert.ok(iterations.some(h=>h.stage==='coupled'&&h.actualNcrit===4&&h.targetNcrit===9));
    assert.ok(frames.some(f=>f.actualNcrit===4&&f.targetNcrit===9));
    assert.ok(checkpoints.some(c=>c.stage==='boundary-layer-initialization'&&c.actualNcrit===4&&c.targetNcrit===9));
    assert.equal(result.initialization.ncritStartupRecovery.attempted,true);
  }finally{h.release();}
});

test('public method observations report the selected seed, linear startup and native logarithmic update 41',async()=>{
  const receipt={method:'MRCHUE surfaces and ISET-style wake guess',wakeInitialization:'iset-linear-shape'};
  const h=await harness({ncritRecovery:true,logRecovery:true,initializationReceipt:receipt});
  const stages=[],iterations=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4},
      {maxIterations:40,onStage:s=>stages.push(s),onIteration:h=>iterations.push(h)});
    const selected=stages.find(s=>s.startupAttempt===1&&s.boundaryLayerInitialization);
    assert.deepEqual(selected.boundaryLayerInitialization,{...receipt,thicknessFactor:1,accepted:true});
    assert.equal(selected.shearCoordinate,'linear');assert.equal(selected.hkFloorLinearization,'exact');
    const ordinary=iterations.filter(h=>h.stage==='coupled'&&h.startupAttempt===1);
    assert.ok(ordinary.length>0);assert.ok(ordinary.every(h=>h.shearCoordinate==='linear'&&h.hkFloorLinearization==='exact'));
    const lastLinear=iterations.find(h=>h.stage==='coupled'&&h.startupAttempt===2&&h.iteration===40);
    assert.equal(lastLinear.shearCoordinate,'linear');assert.equal(lastLinear.hkFloorLinearization,'native');
    const recovery=stages.find(s=>s.shearRecovery);
    assert.equal(recovery.shearCoordinate,'logarithmic');assert.equal(recovery.hkFloorLinearization,'native');
    assert.equal(recovery.shearRecovery.originalIterations,40);assert.equal(recovery.shearRecovery.additionalIterations,40);
    const firstLog=iterations.find(h=>h.stage==='coupled'&&h.startupAttempt===2&&h.iteration===41);
    assert.equal(firstLog.shearCoordinate,'logarithmic');assert.equal(firstLog.hkFloorLinearization,'native');
    assert.equal(firstLog.shearRecovery.shearCoordinate,'logarithmic');
    assert.equal(result.converged,true);assert.equal(result.checkpoint.continuation.shearCoordinate,'logarithmic');
    assert.equal(h.calls.ncritRecoveries[0].plan.shearCoordinate,'logarithmic');
  }finally{h.release();}
});

test('incomplete startup continuation retains its actual intermediate Ncrit and marks the requested case unconverged',async()=>{
  const h=await harness({ncritRecovery:true,recoveryActualNcrit:7.5,pairedFlow:true});
  try{
    const result=h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4},{maxIterations:3});
    assert.equal(result.converged,false);assert.equal(result.stateConverged,true);assert.equal(result.actualNcrit,7.5);
    assert.equal(result.conditions.ncrit,7.5);assert.equal(result.restart.options.ncrit,7.5);
    assert.equal(result.ncritContinuation.reachedTarget,false);assert.match(result.reason,/Ncrit 9.*Ncrit 7.5/);
    assert.equal(h.calls.initializers.length,2);assert.equal(h.calls.euler.length,1);
  }finally{h.release();}
});

test('successful first startups, disabled retry budgets and initialization exceptions never launch the Ncrit recovery',async()=>{
  for(const config of [{},{ncritRecovery:true}])for(const options of [{maxStartupAttempts:1},{maxIterations:0}]){
    const h=await harness(config);
    try{h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9},{maxIterations:3,...options});
      assert.equal(h.calls.ncritRecoveries.length,0);assert.equal(h.calls.initializers.length,1);
    }finally{h.release();}
  }
  const success=await harness();
  try{success.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9},{maxIterations:3});
    assert.equal(success.calls.ncritRecoveries.length,0);assert.equal(success.calls.initializers.length,1);
  }finally{success.release();}
  const error=Object.assign(new Error('Invalid displacement initialization'),{code:'coupled-initialization'});
  const failed=await harness({initializationError:error});
  try{assert.throws(()=>failed.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9},{maxIterations:3}),e=>e===error);
    assert.equal(failed.calls.ncritRecoveries.length,0);assert.equal(failed.calls.initializers.length,1);
  }finally{failed.release();}
});

test('Ncrit startup completion propagates cancellation before another initialization or false target report',async()=>{
  const h=await harness({ncritRecovery:true});let caught=false;
  try{
    try{h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9},{maxIterations:3,
      onStage:event=>{if(event.actualNcrit===9&&h.calls.ncritRecoveries.length)throw undefined;}});}
    catch(error){caught=true;assert.equal(error,undefined);}
    assert.equal(caught,true);assert.equal(h.calls.initializers.length,2);assert.equal(h.calls.ncritRecoveries.length,1);
  }finally{h.release();}
});

test('a rejected second initializer restores the previous complete checkpoint, pressure and geometry',async()=>{
  const error=Object.assign(new Error('thin lower-Ncrit grid rejected'),{code:'coupled-initialization'});
  const h=await harness({ncritRecovery:true,pairedFlow:true,secondInitializationError:error}),frames=[],checkpoints=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly({...baseCase(),gridIntervals:16,ncrit:9},{maxIterations:3,
      onFlow:f=>frames.push(f),onIterationCheckpoint:cp=>checkpoints.push(cp)});
    assert.equal(result.converged,false);assert.equal(result.conditions.ncrit,9);assert.equal(result.actualNcrit,9);
    assert.equal(result.ncritContinuation.reachedTarget,true);assert.equal(result.ncritContinuation.stateConverged,false);
    assert.equal(result.ncritContinuation.startup.reachedTarget,false);
    assert.equal(result.ncritContinuation.startup.failure.code,'coupled-initialization');
    assert.equal(result.initialization.attempts.length,2);assert.equal(result.initialization.attempts[1].actualNcrit,4);
    assert.equal(h.calls.ncritRecoveries.length,0);assert.equal(h.calls.solves.length,1);
    assert.deepEqual(frames.at(-1).checkpoint,checkpoints.at(-1));
    assert.deepEqual(frames.at(-1).flow.nodes,result.flow.nodes);
    assert.equal(frames.at(-1).actualNcrit,9);
  }finally{h.release();}
});


test('earned extension resumes the complete first requested state once and offsets every history/pressure/checkpoint consistently',async()=>{
  const h=await harness({earnedExtension:true,pairedFlow:true}),input={...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4};
  const iterations=[],frames=[],checkpoints=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:40,onIteration:x=>iterations.push(x),
      onIterationCheckpoint:(cp,d)=>checkpoints.push({cp:clone(cp),d:clone(d)}),onFlow:x=>frames.push(x)});
    assert.deepEqual(h.calls.solves.map(x=>x.maxIterations),[20,20,20]);
    assert.equal(h.calls.initializers.length,1);assert.equal(h.calls.ncritRecoveries.length,0);
    assert.equal(result.converged,true);assert.equal(result.checkpoint.restart.options.ncrit,9);
    assert.deepEqual(result.history.map(h=>h.iteration),Array.from({length:46},(_,i)=>i));
    assert.deepEqual(iterations.filter(x=>x.stage==='coupled').map(h=>h.iteration),Array.from({length:46},(_,i)=>i));
    assert.equal(result.linearDiagnostics.solves,45);
    assert.deepEqual(result.linearDiagnostics.iterations.map(h=>h.iteration),Array.from({length:45},(_,i)=>i+1));
    assert.equal(result.initialization.attempts[0].iterations,45);
    assert.equal(result.startupExtension.original.history.length,41);assert.equal(result.startupExtension.iterations,5);
    assert.equal(result.startupExtension.original.checkpoint.testCompletedIterations,40);
    assert.equal(result.startupExtension.original.linearDiagnostics.solves,40);
    const last=checkpoints.at(-1);assert.equal(last.d.iterationOffset,40);assert.equal(last.d.history.at(-1).iteration,5);
    assert.equal(last.cp.testCompletedIterations,45);
    for(const frame of frames.filter(x=>x.iteration.iteration>40)){
      assert.equal(frame.iteration.iteration,frame.checkpoint.testCompletedIterations);
      assert.equal(frame.checkpoint.restart.options.ncrit,9);
      assert.equal(frame.flow.cells[0][0][0].interfacePressure.lower,frame.checkpoint.restart.initialBL[0]);
      assert.deepEqual(frame.flow.nodes,frame.checkpoint.restart.initialEuler.nodes);
    }
    assert.equal(frames.at(-1).iteration.iteration,45);
  }finally{h.release();}
});

test('an exhausted extension reserve retains its last accepted state and complete original evidence',async()=>{
  const h=await harness({earnedExtension:true,extensionFailed:true,pairedFlow:true}),input={...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4};
  const iterations=[],frames=[],checkpoints=[];
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:40,maxStartupAttempts:1,
      onIteration:x=>iterations.push(x),onFlow:x=>frames.push(x),
      onIterationCheckpoint:(cp,d)=>checkpoints.push({cp:clone(cp),d:clone(d)})});
    assert.equal(result.converged,false);assert.equal(result.checkpoint.testCompletedIterations,120);
    assert.equal(result.history.length,121);assert.equal(result.linearDiagnostics.solves,120);
    assert.equal(result.initialization.attempts[0].iterations,120);
    assert.equal(result.startupExtension.original.checkpoint.testCompletedIterations,40);
    assert.equal(result.startupExtension.original.history.length,41);
    assert.equal(result.startupExtension.iterations,80);assert.equal(result.startupExtension.converged,false);
    assert.equal(result.initialization.attempts[0].startupExtension.original.checkpoint.testCompletedIterations,40);
    assert.deepEqual(h.calls.solves.map(x=>x.maxIterations),[20,20,20,20,20,20]);assert.equal(h.calls.initializers.length,1);
    assert.deepEqual(iterations.filter(x=>x.stage==='coupled').map(x=>x.iteration),Array.from({length:121},(_,i)=>i));
    assert.deepEqual(result.linearDiagnostics.iterations.map(x=>x.iteration),Array.from({length:120},(_,i)=>i+1));
    assert.deepEqual(result.startupExtension.chunks.map(x=>x.offset),[40,60,80,100]);
    assert.equal(checkpoints.at(-1).d.iterationOffset,100);
    assert.equal(checkpoints.at(-1).cp.testCompletedIterations,120);
    for(const frame of frames.filter(x=>x.iteration.iteration>40))
      assert.equal(frame.iteration.iteration,frame.checkpoint.testCompletedIterations);
  }finally{h.release();}
});

test('extension observer cancellation propagates every thrown value without a lower-Ncrit retry',async()=>{
  for(const hook of ['onStage','onIteration','onIterationCheckpoint','onMesh','onFlow'])for(const error of [undefined,new Error('stop')]){
    const h=await harness({earnedExtension:true,pairedFlow:true}),input={...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4};
    let caught=false;
    const callback=(value,details)=>{
      const selected=hook==='onStage'?value.startupExtension:hook==='onIterationCheckpoint'?details.iterationOffset===40
        :hook==='onIteration'?value.iteration>40:hook==='onMesh'?value.iteration?.iteration>40:value.iteration.iteration>40;
      if(selected)throw error;
    };
    try{h.solveCoupledStreamtubeAssembly(input,{maxIterations:40,[hook]:callback});}
    catch(value){caught=true;assert.equal(value,error);}
    finally{h.release();}
    assert.equal(caught,true);assert.equal(h.calls.initializers.length,1);assert.equal(h.calls.ncritRecoveries.length,0);
    assert.ok(h.calls.solves.length<=3);
  }
});

test('discarding a failed seed restores the retained residual progress with its mesh and checkpoint', async () => {
  const h = await harness({ earnedExtension: true, extensionFailed: true, failedSeedRecovery: true, pairedFlow: true });
  const iterations = [], frames = [];
  try {
    const r = h.solveCoupledStreamtubeAssembly({ ...baseCase(), gridIntervals: 16, eulerIsmom: 4 }, {
      direct: true, maxIterations: 40, maxStartupAttempts: 1, machRecovery: false, resolutionRecovery: false,
      onIteration: row => iterations.push(row), onFlow: frame => frames.push(frame),
    });
    assert.equal(r.converged, false);
    assert.equal(r.seedRecovery.accepted, false);
    assert.equal(iterations.at(-2).residual, 99);
    const last = iterations.at(-1), frame = frames.at(-1);
    assert.equal(last.retained, true);
    assert.equal(last.iteration, r.history.at(-1).iteration);
    assert.equal(last.residual, Math.max(...Object.values(r.families)));
    assert.equal(frame.iteration.iteration, last.iteration);
    assert.deepEqual(frame.checkpoint.families, r.families);
    assert.equal(r.history.length, 121, 'Restoration is a display event, not another numerical update.');
  } finally { h.release(); }
});


test('a failed earned extension can still select lower-Ncrit recovery and keeps its original and extended histories',async()=>{
  const h=await harness({earnedExtension:true,extensionFailed:true,pairedFlow:true}),input={...baseCase(),gridIntervals:16,ncrit:9,eulerIsmom:4};
  try{
    const result=h.solveCoupledStreamtubeAssembly(input,{maxIterations:40});
    assert.deepEqual(h.calls.initializers.map(c=>[c.ncrit,c.controls.initialThicknessFactor]),[[9,1],[4,.25]]);
    assert.equal(result.initialization.attempts.length,2);
    const first=result.initialization.attempts[0];assert.equal(first.iterations,120);assert.equal(first.history.length,121);
    assert.equal(first.startupExtension.original.history.length,41);
    assert.equal(first.startupExtension.original.checkpoint.testCompletedIterations,40);
    assert.equal(result.checkpoint.restart.options.ncrit,4);assert.equal(result.actualNcrit,4);assert.equal(result.targetNcrit,9);
    assert.equal(result.converged,false);assert.equal(result.ncritContinuation.reachedTarget,false);
    assert.deepEqual(h.calls.solves.map(x=>x.maxIterations),[20,20,20,20,20,20,20,20]);
  }finally{h.release();}
});

test('direct app policy reuses its final-grid precursor for one earned coupled extension', async () => {
  const h = await harness({ earnedExtension: true }), input = { ...baseCase(), mach: .74 };
  const stages = [];
  try {
    const r = h.solveCoupledStreamtubeAssembly(input, {
      direct: true, coarseStartup: false, maxIterations: 40, onStage: event => stages.push(event),
    });
    assert.deepEqual(h.calls.solves, [{ grid: 32, maxIterations: 40 }, { grid: 32, maxIterations: 20 }]);
    assert.equal(h.calls.euler.length, 1);
    assert.equal(h.calls.sequences.length, 0);
    assert.equal(h.calls.recoveries.length, 0);
    assert.equal(h.calls.ncritRecoveries.length, 0);
    assert.equal(r.solverSettings.maxStartupAttempts, 1);
    assert.equal(r.solverSettings.automaticTransitionRefinement, false);
    assert.equal(r.converged, true, 'An earned extension must still solve the requested equations.');
    assert.ok(stages.every(s => s.mach === .74 && s.actualAlpha === 2.68 && s.gridLevel === 32));
    assert.ok(stages.every(s => ['euler', 'boundary-layer-initialization', 'coupled'].includes(s.stage)));
  } finally { h.release(); }
});

test('direct viscous stall does not trigger transition refinement or a second startup', async () => {
  const h = await harness({ stallRecovery: true }), stages=[];
  try {
    const r=h.solveCoupledStreamtubeAssembly(baseCase(), {
      direct:true,coarseStartup:false,maxIterations:40,onStage:e=>stages.push(e.stage),
    });
    assert.equal(h.calls.solves.length,1);
    assert.equal(h.calls.recoveries.length,0);
    assert.equal(h.calls.sequences.length,0);
    assert.equal(h.calls.ncritRecoveries.length,0);
    assert.equal(r.initialization.attempts.length,1);
    assert.ok(stages.every(s=>['euler','boundary-layer-initialization','coupled'].includes(s)));
  } finally {h.release();}
});

 test('direct automatic sequencing preserves Ncrit9 and all physical conditions through the fine solve', async () => {
  const h = await harness({finite:true,subsonicEuler:true,sourceGeometryReplay:'native'}), input = {...finiteCase(), ncrit:9, gridTubes:24};
  const before = clone(input), stages=[];
  try {
    const r=h.solveCoupledStreamtubeAssembly(input,{direct:true,maxIterations:40,onStage:s=>stages.push(s)});
    assert.equal(r.converged,true);
    assert.deepEqual(input,before);
    assert.deepEqual(h.calls.euler.map(c=>[c.gridIntervals,c.gridTubes]),[[16,7]]);
    assert.deepEqual(h.calls.solves.map(c=>c.grid),[16,36]);
    assert.equal(h.calls.continuations.length,0);
    assert.ok(h.calls.initializers.every(c=>c.ncrit===9));
    assert.equal(r.actualNcrit??r.conditions.ncrit,9);
    assert.equal(r.restart.options.edgeMatching,'pressure');
    assert.equal(r.restart.options.blThermodynamics,undefined);
    assert.equal(r.restart.options.geometryReplay,'native');
    assert.ok(stages.filter(s=>s.coarseStage).every(s=>s.gridLevel===16&&s.requestedGridIntervals===32));
  } finally {h.release();}
 });
