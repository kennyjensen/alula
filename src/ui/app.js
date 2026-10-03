// SPDX-License-Identifier: GPL-2.0-or-later
import { naca4Standard as naca4, transform } from '../geometry/airfoil.js';
import { parseCoordinates } from '../geometry/parse.js';
import { prepareAirfoilElement } from '../geometry/airfoil-element.js';
import { serializeCoordinates } from '../geometry/serialize.js';
import { flowModelRoute } from './flow-model.js';
import { getBenchmarkAirfoil } from '../geometry/benchmark-airfoils.js';
import { builtinAirfoils as presets } from '../geometry/builtin-airfoils.js';
import { colors, assemblyBounds, drawGeometry, drawPressure, drawBoundaryLayer } from './plots.js';
import { bindPlotNavigation } from './plot-navigation.js';
import { bindWorkspaceTabs } from './workspace-tabs.js';
import { flowSpeedGradient } from './flow-colors.js';
import { createQuadCoupledFlowCache } from './quad-coupled-flow-cache.js';
import { quadCoupledDisplayedStateLabel } from './quad-coupled-failure.js';
import { quadCoupledNcrit, quadCoupledNcritLabel, quadCoupledNcritResult } from './quad-coupled-ncrit.js';
import { quadCoupledGrid, quadCoupledGridLabel, quadCoupledGridResult } from './quad-coupled-grid.js';
import { buildSolverErrorContext } from './solver-error-context.js';
import { quadSmoothingLabel } from './quad-smoothing-status.js';
import { createQuadSolveProgress, quadStartingGridLabel } from './quad-solve-progress.js';

const $ = id => document.getElementById(id);
const escapeHTML = value => String(value).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);

let definitions = structuredClone(presets.flap);
let caseName = 'Main element + flap';
let selectedBenchmark = null;
let coordinateDocument = null;
let elements = []; let result = null; let lines = []; let bounds; let worker; let jobId = 0; let stale = true;
let zoom = 1; let pan = { x: 0, y: 0 }; let plotPoints = []; let hover = null;
let livePressure = null;
let overlayCoefficients = null;
let lastCase = null; let solverBusy = false;
let activeRequest = null, lastErrorReport = null;
let pendingMesh = null; let meshPhase = null;
let gridAudit = null, auditingGrid = false;
const quadFlowCache = createQuadCoupledFlowCache();
const workspaceTabs = bindWorkspaceTabs({ shell: document.querySelector('.app-shell'),
  tabs: $('workspace-tabs'), settings: $('settings-pane'), visualization: $('visualization-pane'),
  actions: $('analysis-actions'), actionSlot: $('analysis-actions-slot') });
const geometryCanvas = $('geometry-canvas'); const pressureCanvas = $('pressure-canvas');
const blCanvas = $('bl-canvas');
$('flow-speed-bar').style.background = flowSpeedGradient;
const isStreamtubeGrid=()=>['streamtube-grid','streamtube-bl'].includes($('flow-model').value);
const isQuadCoupled=()=>$('flow-model').value==='streamtube-bl';
const quadAutomatic=()=>$('quad-transition').value==='automatic';
const quadTripsByMode={'fixed-trip':['0.05','0.05'],automatic:['1','1']};
let previousQuadTransition=$('quad-transition').value;
const quadSolveLabel=()=>isQuadCoupled()?'Run quad Euler/BL':'Run quad Euler';
const quadIterationLabel=mesh=>mesh?.iteration?.stage==='coupled-coarse-initialization'?'Coarse startup':mesh?.iteration?.stage?.startsWith('coupled')?'Euler/BL':'Euler';
const plottedResult=()=>!stale||isQuadCoupled()&&['failed','rejected'].includes(meshPhase)?result:null;
const plottedPressure=()=>isStreamtubeGrid()&&!isQuadCoupled()?null:livePressure??plottedResult();
const isExperimentalMesh=mesh=>mesh?.initialization?.gridSmoothing?.experimental===true;
const experimentalMeshTitle=mesh=>mesh.initialization.gridSmoothing.retainedOriginal?'Experimental smoothing · original mesh retained':mesh.initialization.gridSmoothing.converged?'Experimental smoothed mesh':'Experimental grid · smoothing incomplete';
const experimentalAuditExplanation='Check grid is unavailable for this preview: its polygonal boundaries and uniform cross-line target differ from the curved harmonic smoother. The preview accuracy report below does not establish physical acceptance.';
const hasVolumeMesh=isStreamtubeGrid;
const isCoupled=()=>$('flow-model').value==='coupled';
// The two quad modes share controls on screen but retain independent Mach,
// equation and mesh choices. BL-specific controls already have their own fields.
const quadSettingIds = ['quad-mach','euler-ismom','euler-startup','grid-intervals','grid-tubes','grid-inlet','grid-outlet',
  'grid-chord-exponent','grid-upper-streamlines','grid-lower-streamlines','grid-gap-streamlines',
  'grid-surface-spacing','grid-le-ratio','grid-te-ratio','grid-curvature-exponent',
  'grid-elliptic','grid-match-aspect','grid-aspect-ratio'];
const readQuadSettings = () => Object.fromEntries(quadSettingIds.map(id => [id,
  $(id).type === 'checkbox' ? $(id).checked : $(id).value]));
const quadSettings = { 'streamtube-grid': readQuadSettings(), 'streamtube-bl': readQuadSettings() };
const eulerFormulations = { 1: 'Momentum', 2: 'Entropy', 3: 'Leading-edge hybrid', 4: 'Automatic hybrid' };
let previousFlowMode = $('flow-model').value;
const panelLevels=[40,60,80,120,160,240,320,400,480,600];
const panelCount=(d,defaultCount=Number($('resolution').value))=>d.points?d.points.length-1:d.panels??defaultCount;
const analysisPanelCount=(d,defaultCount)=>{
  const n=panelCount(d,defaultCount), te=d.trailingEdge;
  if(!d.points||te?.kind!=='finite-base')return n;
  if(![te.upperIndex,te.lowerIndex].every(i=>Number.isInteger(i)&&i>=0&&i<n)||te.upperIndex===te.lowerIndex)
    throw new Error('Invalid finite trailing-edge corner indices.');
  return (te.lowerIndex-te.upperIndex+n)%n;
};
const canRefineCoupled = () => isQuadCoupled() && !solverBusy && !stale && result?.converged === true
  && result.model === 'research-streamtube-euler-bl' && result.checkpoint?.version === 1
  && result.checkpoint.restart.options.blThermodynamics !== 'historical-common-isentrope'
  && result.sourceCase && JSON.stringify(result.sourceCase) === JSON.stringify(lastCase);
function syncCoupledRefinement() {
  $('refine-coupled-button').hidden = $('refine-coupled-help').hidden = !isQuadCoupled();
  $('refine-coupled-button').disabled = !canRefineCoupled();
  $('refine-coupled-button').title = result?.checkpoint?.restart?.options?.blThermodynamics === 'historical-common-isentrope'
    ? 'Refinement of the transonic model is not validated yet.' : '';
}

function syncResolution(){
  for(const option of $('resolution').options)option.disabled=definitions.reduce((n,d)=>n+panelCount(d,Number(option.value)),0)>700;
  if($('resolution').selectedOptions[0].disabled){
    const allowed=[...$('resolution').options].filter(o=>!o.disabled);
    if(allowed.length)$('resolution').value=allowed.at(-1).value;
  }
  $('resolution').disabled=definitions.every(d=>d.points||d.panels!==undefined);
  document.querySelectorAll('[data-key="panels"]').forEach(select=>{
    const index=Number(select.dataset.element);
    for(const option of select.options){
      const count=option.value===''?Number($('resolution').value):Number(option.value);
      option.disabled=count+definitions.reduce((n,d,i)=>n+(i===index?0:panelCount(d)),0)>700;
    }
  });
}

function buildElements({ includeTrips = true } = {}) {
  const panels = Number($('resolution').value);
  return definitions.map(d => prepareAirfoilElement({ name: d.name, ...(d.sourcePoints?{sourcePoints:d.sourcePoints}:{}), points: d.points ?? transform(naca4(d.code, d.panels??panels), {
    chord: Number(d.chord), x: Number(d.x), y: Number(d.y), angle: -Number(d.deflection),
  }),...(d.trailingEdge?{trailingEdge:structuredClone(d.trailingEdge)}:{}),...(includeTrips&&isCoupled()&&(d.tripUpper!==undefined||d.tripLower!==undefined)?{trips:[d.tripUpper??Number($('trip-upper').value),d.tripLower??Number($('trip-lower').value)]}:{}) }));
}
function refreshCoordinateExport() {
  const select = $('coordinate-element'), previous = select.value;
  select.innerHTML = '<option value="all">All elements</option>' + definitions.map((d, i) =>
    `<option value="${i}">${escapeHTML(d.name)}</option>`).join('');
  if ([...select.options].some(o => o.value === previous)) select.value = previous;
  const format = $('coordinate-format').value, dat = format === 'xfoil-labeled';
  select.options[0].disabled = dat;
  if (dat && select.value === 'all') select.value = '0';
  $('coordinate-points').disabled = format === 'json';
  $('coordinate-domain').hidden = format !== 'mses' || Boolean(coordinateDocument?.header || coordinateDocument?.domain);
  $('coordinate-export-note').textContent = format === 'json'
    ? 'JSON retains original coordinates, measured base corners and all geometry metadata.'
    : 'Standard text files cannot retain measured base corner labels. Original points retain the source polyline; wetted surface explicitly omits the base. Use JSON for a lossless geometry record.';
}
function exportCoordinates() {
  try {
    const format = $('coordinate-format').value, selected = $('coordinate-element').value;
    const elementIndex = selected === 'all' ? undefined : Number(selected);
    const data = { ...coordinateDocument, name: coordinateDocument?.name ?? caseName,
      elements: buildElements({ includeTrips: false }), referenceChord: Number($('reference').value),
      ...(selectedBenchmark ? { geometrySource: structuredClone(selectedBenchmark.provenance) } : {}) };
    if (format === 'mses' && !data.header && !data.domain) {
      const value = id => {
        if (!$(id).value || !$(id).validity.valid) throw new Error('Enter all four MSES domain bounds before exporting a new blade file.');
        return Number($(id).value);
      };
      data.domain = { xMin: value('coordinate-x-min'), xMax: value('coordinate-x-max'),
        yMin: value('coordinate-y-min'), yMax: value('coordinate-y-max') };
    }
    const content = format === 'json' ? JSON.stringify({ schemaVersion: 1, application: 'alula geometry', ...data,
      elements: elementIndex === undefined ? data.elements : [data.elements[elementIndex]] }, null, 2) + '\n'
      : serializeCoordinates(data, { format, coordinates: $('coordinate-points').value, elementIndex });
    const url = URL.createObjectURL(new Blob([content], { type: format === 'json' ? 'application/json' : 'text/plain' }));
    const link = document.createElement('a'); link.href = url;
    link.download = format === 'json' ? 'alula-geometry.json' : format === 'mses' || format === 'source' && data.format === 'mses' ? 'blade.export' : 'airfoil.dat';
    link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
    $('error-message').hidden = true;
  } catch (error) { showError(error.message); }
}
function renderElements() {
  $('element-list').innerHTML = definitions.map((d, i) => `<details class="element-card" ${i === definitions.length - 1 ? 'open' : ''}>
    <summary><span class="element-swatch" style="background:${colors[i]}"></span>${escapeHTML(d.name)}<small>${d.points ? 'DAT' : `NACA ${escapeHTML(d.code)}`}</small></summary>
    <div class="element-fields">${d.points ? `<div class="imported-description">${d.points.length - 1} panels · original coordinates.<br>Imported contours retain their original spacing.</div>` : [
      ['code', 'NACA', 'text', ''], ['chord', 'Chord', 'number', 'min="0.01" step="0.01"'],
      ['deflection', 'Deflection ↓ °', 'number', 'step="1"'], ['x', 'LE x', 'number', 'step="0.01"'], ['y', 'LE y', 'number', 'step="0.01"'],
    ].map(([key, label, type, attributes]) => `<label>${label}<input data-element="${i}" data-key="${key}" aria-label="${escapeHTML(d.name)} ${label}" type="${type}" value="${escapeHTML(d[key])}" ${attributes} ${key === 'code' ? 'maxlength="4" pattern="[0-9]{4}"' : ''}></label>`).join('')+`<label>Panels<select data-element="${i}" data-key="panels" aria-label="${escapeHTML(d.name)} panels"><option value="" ${d.panels===undefined?'selected':''}>Use default</option>${panelLevels.map(n=>`<option value="${n}" ${d.panels===n?'selected':''}>${n}</option>`).join('')}</select></label>`}</div>
    <div class="element-trips" ${isCoupled()?'':'hidden'}><small>Transition trips · blank uses global</small><div class="two-fields">${[['tripUpper','Upper'],['tripLower','Lower']].map(([key,label])=>`<label>${label} x/c<input data-element="${i}" data-key="${key}" aria-label="${escapeHTML(d.name)} ${label} trip" type="number" min="0.000001" max="1" step="any" placeholder="Global" value="${d[key]??''}"></label>`).join('')}</div></div>
    <div class="quad-element-trips" ${isQuadCoupled()?'':'hidden'}><small>Trip limits · surface fractions · blank uses global</small><div class="two-fields">${[['quadTripUpper','Upper'],['quadTripLower','Lower']].map(([key,label])=>`<label>${label}<input data-element="${i}" data-key="${key}" aria-label="${escapeHTML(d.name)} ${label} quad trip" type="number" min="0.000001" max="${quadAutomatic()?'1':'0.999999'}" step="any" placeholder="Global" value="${d[key]??''}"></label>`).join('')}</div></div>
    ${definitions.length > 1 ? `<div class="element-actions"><button class="remove-element" data-remove="${i}">Remove element</button></div>` : ''}</details>`).join('');
  $('element-count').textContent = `${definitions.length} ELEMENT${definitions.length === 1 ? '' : 'S'}`;
  $('add-element').disabled = definitions.length >= 6;
  syncResolution();
  refreshCoordinateExport();
  $('geometry-legend').innerHTML = definitions.map((d, i) => `<span class="legend-item" style="--line-color:${colors[i]}">${escapeHTML(d.name)}</span>`).join('');
}
function setStatus(text, state = '') { $('status').textContent = text; $('status').className = `status-pill ${state}`; }
function experimentalAccuracyText(mesh) {
  const smoothing=mesh.initialization.gridSmoothing, accuracy=smoothing.accuracy;
  const measured=Number.isFinite(accuracy?.maximumLocalIntervals);
  const mismatch=measured?`${accuracy.maximumLocalIntervals.toFixed(3)} local intervals`:'unresolved';
  const target=Number.isFinite(accuracy?.limit)?accuracy.limit:0.05;
  const invalidCells=smoothing.regions?.reduce((n,r)=>n+(r.quality?.invalidCells?.length??0),0)??0;
  return `Streamline mismatch: ${mismatch}; target ≤ ${target} ${accuracy?.passed===true?'met':'not met or unresolved'}. Numerical convergence: ${smoothing.converged?'yes':'no'}. Curved cell geometry: ${smoothing.geometryCertified?'locally certified':invalidCells?`${invalidCells} rejected cells`:'unresolved'}. Physical acceptance pending.${smoothing.reason?` ${smoothing.reason}`:''}`;
}
function setBusy(busy, meshOnly = false) {
  solverBusy = busy;
  $('solve-button').disabled = busy; $('stop-button').hidden = !busy;
  $('solve-button').title = '';
  $('solve-button').firstElementChild.textContent = busy && !meshOnly ? 'Solving…' : 'Run';
  syncCoupledRefinement();
}
function syncModelUI() {
  syncCoupledRefinement();
  const coupled = isCoupled();
  for (const id of ['grid-le-ratio', 'grid-te-ratio', 'grid-curvature-exponent']) $(id).disabled = $('grid-surface-spacing').value !== 'curvature';
  $('grid-aspect-ratio').disabled = !$('grid-match-aspect').checked;
  for (const id of ['grid-upper-streamlines', 'grid-lower-streamlines', 'grid-gap-streamlines']) $(id).placeholder = `Auto · ${Number($('grid-tubes').value) + 1}`;
  for (const id of ['grid-inlet', 'grid-outlet'])
    $(id).options[0].textContent = `Auto · ${Math.ceil(Number($('grid-intervals').value) / 2)}`;
  document.body.classList.toggle('grid-preview', isStreamtubeGrid());
  document.querySelector('.pressure-panel').hidden = isStreamtubeGrid()&&!isQuadCoupled();
  $('show-streamlines').closest('label').hidden = false;
  $('show-streamlines').title = isStreamtubeGrid() ? 'Current Euler streamtube centerlines, colored by speed relative to freestream' : '';
  $('flow-mesh-note').textContent = isStreamtubeGrid()
    ? isQuadCoupled() ? `Euler flow, quad nodes, all surface boundary layers and wakes solve together. ${quadAutomatic()?'Automatic-transition':'Fixed-trip'} research model; physical accuracy remains unvalidated.` : 'Euler flow and quad node positions solve together. Select Quad Euler + boundary layers for viscous flow.'
    : coupled ? 'Surface panels and boundary layers solve together. Incompressible flow; no volume mesh.' : 'Incompressible surface-panel flow; no volume mesh.';
  $('viscous-conditions').hidden = !coupled;
  document.querySelectorAll('.element-trips').forEach(node=>{node.hidden=!coupled;});
  document.querySelectorAll('.quad-element-trips').forEach(node=>{node.hidden=!isQuadCoupled();});
  $('quad-viscous-conditions').hidden=!isQuadCoupled();
  for(const id of ['quad-trip-upper','quad-trip-lower'])$(id).max=quadAutomatic()?'1':'0.999999';
  document.querySelectorAll('.quad-element-trips input').forEach(node=>{node.max=quadAutomatic()?'1':'0.999999';});
  $('quad-trip-note').textContent='Trips are fractions of each initial surface: 0 at stagnation, 1 at the trailing edge, not chord x/c. Element overrides take precedence. '
    +(quadAutomatic()?'Use 1 to allow natural transition or a laminar surface to the TE. A smaller value forces transition if reached first.':'Fixed-trip mode prescribes transition at each trip. Select Automatic (eᴺ) to solve for natural transition.');
  $('fixed-mach').hidden=hasVolumeMesh();
  $('quad-mach-help').hidden=!isQuadCoupled();
  $('euler-startup-field').hidden=isQuadCoupled();
  $('streamtube-grid-conditions').hidden=!isStreamtubeGrid();
  $('solve-button').hidden=false;
  $('grid-audit').hidden = true;
  $('grid-build-help').textContent = `${quadSolveLabel()} builds the initial grid and updates the grid and speed-colored streamlines as it solves; Stop keeps the latest view.`;
  $('flow-tag').textContent = isStreamtubeGrid()?(isQuadCoupled()?'EULER + BL':'EULER · NO BL'):coupled ? 'VISCOUS' : 'INVISCID';
  $('flow-toggle-label').textContent = isQuadCoupled() ? 'Flow & BL' : coupled ? 'BL & wake' : 'Streamlines';
  const finiteBasePanel = $('flow-model').value === 'inviscid' && definitions.some(d => d.trailingEdge?.kind === 'finite-base');
  $('residual-label').textContent = isQuadCoupled() ? 'Euler/BL equation residual' : isStreamtubeGrid() ? 'Euler equation residual' : coupled ? 'Coupled equation residual' : finiteBasePanel ? 'Panel equation residual' : 'Boundary residual ‖Vn‖∞';
  $('residual-description').textContent = isStreamtubeGrid() ? 'Maximum equation error · target 10⁻¹⁰' : coupled ? 'Maximum scaled equation error · target 10⁻⁸' : finiteBasePanel ? 'Nodal streamfunction and Kutta equation error' : 'No-penetration equation error';
  $('cp-description').textContent = isQuadCoupled() ? 'BL-edge pressure from solved edge speed · research, unvalidated' : coupled ? 'Solid: viscous · dashed: inviscid' : 'Negative Cp upwards';
}
function updateAirfoilOverlay() {
  const state = overlayCoefficients;
  const conditions = state?.conditions;
  const input = activeRequest?.input ?? (!stale ? lastCase : null);
  const viscous = isCoupled() || isQuadCoupled();
  const numeric = (...values) => values.find(Number.isFinite);
  const mach = numeric(conditions?.coefficientConditions?.mach, conditions?.conditions?.mach,
    conditions?.actualMach, conditions?.mach, input?.mach, isStreamtubeGrid() ? Number($('quad-mach').value) : 0);
  const alpha = numeric(conditions?.actualAlpha, conditions?.alpha, conditions?.sourceCase?.alpha,
    input?.alpha, Number($('alpha-number').value));
  const re = viscous ? numeric(conditions?.conditions?.reynolds, conditions?.sourceCase?.reynolds,
    input?.reynolds, Number($(isQuadCoupled() ? 'quad-reynolds' : 'reynolds').value)) : undefined;
  const ncr = viscous ? numeric(quadCoupledNcrit(conditions).actualNcrit, conditions?.conditions?.ncrit,
    conditions?.sourceCase?.ncrit, input?.ncrit, Number($(isQuadCoupled() ? 'quad-ncrit' : 'ncrit').value)) : undefined;
  const { cl, cm, cd } = state?.values ?? {};
  const values = { mach, alpha, re, ncr, cl, cm, cd,
    ld: Number.isFinite(cl) && Number.isFinite(cd) && cd > 0 ? cl / cd : undefined };
  $('overlay-name').textContent = caseName;
  $('airfoil-overlay').classList.toggle('provisional', !!state?.provisional);
  for (const [key, value] of Object.entries(values)) {
    $('overlay-' + key).textContent = !Number.isFinite(value) ? '—'
      : key === 're' ? value.toExponential(2)
        : value.toFixed(({ mach: 3, alpha: 2, ncr: 1, cl: 4, cm: 4, cd: 6, ld: 1 })[key]);
  }
}
function draw() {
  updateAirfoilOverlay();
  const mesh = stale ? pendingMesh : result?.mesh, showMesh = $('show-mesh').checked;
  const liveFlow = isStreamtubeGrid() ? mesh?.flow : null;
  $('flow-legend').hidden = !liveFlow || !$('show-streamlines').checked;
  $('streamline-status').hidden = !!liveFlow;
  $('flow-iteration').hidden = !liveFlow;
  $('flow-change-readout').hidden = !liveFlow;
  if (liveFlow) {
    $('flow-speed-bar').style.background = flowSpeedGradient;
    $('flow-speed-ticks').textContent = '';
    for (const tick of ['0', '1', '2', '≥3']) {
      const label = document.createElement('span'); label.textContent = tick; $('flow-speed-ticks').append(label);
    }
    const phase = meshPhase === 'stopped' ? 'stopped · provisional'
      : meshPhase === 'failed' ? 'unconverged · provisional'
        : meshPhase === 'rejected' ? 'equations converged · grid rejected · provisional'
        : !stale ? 'equations converged · research, unvalidated'
          : 'solving · provisional';
    const progress = liveFlow.iteration === 0 ? 'Initial flow estimate' : `${quadIterationLabel(mesh)} iteration ${liveFlow.iteration}`;
    $('flow-iteration').textContent = `${progress} · ${phase} · residual ${(mesh.iteration?.residual??liveFlow.residual).toExponential(2)}`;
    $('flow-legend').dataset.iteration = String(liveFlow.iteration);
    $('flow-legend').classList.toggle('provisional', stale);
    $('flow-change-readout').textContent = `Largest speed change this step: ${(100 * liveFlow.maximumSpeedChangeFromPrevious).toFixed(3)}% of U∞`;
    $('streamline-status').textContent = `${progress} · ${phase}`;
  }
  const experimental = isExperimentalMesh(mesh);
  const smoothing = mesh?.initialization?.gridSmoothing;
  const slorLabel = mesh?.topology === 'intrinsic-quadrilateral-streamtubes' && !experimental
    ? quadSmoothingLabel(smoothing) : '';
  const meshLabel = { initial: 'initial · unsolved', smoothing: 'SLOR smoothing · unsolved', experimental: 'unsolved · physical acceptance pending', solving: 'solving', failed: 'unconverged', rejected: 'equations converged · grid rejected', stopped: 'stopped · unsolved' }[meshPhase];
  const movement = mesh?.iteration?.maximumNodeMovement;
  const meshIteration = mesh?.iteration ? ` · ${quadIterationLabel(mesh)} iteration ${mesh.iteration.iteration}`
    + (Number.isFinite(movement) ? ` · max move ${(movement / mesh.iteration.movementReferenceChord).toExponential(1)} c_ref` : '') : '';
  const meshName = experimental ? experimentalMeshTitle(mesh) : mesh?.topology === 'intrinsic-quadrilateral-streamtubes' ? 'MSES-style quadrilateral mesh'
    : mesh?.cells.every(cell => cell.length === 3) ? 'Triangular mesh' : 'Flow mesh';
  const incompatibleAudit=experimental;
  $('grid-audit-overlay').disabled = incompatibleAudit || !gridAudit;
  $('grid-audit-export').disabled = !gridAudit && !experimental;
  if(incompatibleAudit) $('grid-audit-summary').textContent = `${experimentalAuditExplanation}${experimental?` ${experimentalAccuracyText(mesh)}`:''}`;
  $('mesh-status').hidden = !hasVolumeMesh();
  $('mesh-update').hidden = !isStreamtubeGrid() || !mesh;
  $('mesh-update').dataset.iteration = String(mesh?.iteration?.iteration ?? 0);
  $('mesh-update').textContent = mesh?.iteration?.iteration > 0
    ? `Grid update ${mesh.iteration.iteration}`
      + (Number.isFinite(movement) ? ` · largest move ${(movement / mesh.iteration.movementReferenceChord).toExponential(2)} c_ref this step` : '')
      + (Number.isFinite(mesh.iteration.maximumNodeMovementFromInitial) ? ` · ${(mesh.iteration.maximumNodeMovementFromInitial / mesh.iteration.movementReferenceChord).toExponential(2)} c_ref from start` : '')
    : quadStartingGridLabel(mesh, meshPhase);
  $('mesh-status').textContent = mesh
    ? `${meshName} · ${mesh.cells.length.toLocaleString()} cells${stale && meshLabel ? ` · ${meshLabel}` : !stale && isStreamtubeGrid() ? ' · equations converged · research' : ''}${meshIteration}${slorLabel ? ` · ${slorLabel}` : ''}${meshPhase === 'initial' && mesh.initialization?.gridRepair?.converged ? ' · repaired' : ''}${!experimental && meshPhase === 'initial' && !smoothing?.converged && !mesh.initialization?.gridRepair?.converged ? ' · unsmoothed' : ''}${mesh.initialization?.surfaceDistributions?.some(d => d.artificialLeadingCurvature) ? ' · local LE refinement' : ''}${mesh.initialization?.normalAllocation?.groups.some(g => g.addedTubes > 0) ? ' · extra streamlines for aspect' : ''}`
    : hasVolumeMesh() ? meshPhase === 'building' ? 'Building initial mesh…' : meshPhase === 'waiting' ? 'Waiting for a mesh from flow initialization…' : isStreamtubeGrid()?`${quadSolveLabel()} to view the grid.`:'Run analysis to view the grid.' : 'This panel flow model has no volume mesh.';
  if (mesh?.initialization?.tubes) $('mesh-status').textContent +=
    ` · Actual streamlines, bottom → top passages [${mesh.initialization.tubes.map(n => n + 1).join(', ')}]`
    + (mesh.initialization.surfaceIntervals ? ` · Surface points/side: ${mesh.initialization.surfaceIntervals
      .slice().sort((a, b) => a.element - b.element)
      .map(b => `${elements[b.element]?.name ?? `Element ${b.element + 1}`} ${b.intervals + 1}`).join('; ')}` : '');
  $('mach-contour-status').hidden = !$('show-mach-contours').checked;
  $('mach-contour-status').textContent = mesh?.flow?.lines?.some(line => line.machNumbers?.some(Number.isFinite))
    ? `Mach contours · ΔM = 0.05 · gold: M = 1 · ${stale ? 'current iterate · unconverged' : 'retained solution'} · interpolated between tube centers`
    : 'Mach contours available after an Euler flow update.';
  const displayed = plottedResult();
  drawGeometry(geometryCanvas, { elements, lines, bounds, zoom, pan,
    showMachContours: $('show-mach-contours').checked,
    showStreamlines: $('show-streamlines').checked, showPanels: $('show-panels').checked,
    showMesh, mesh, gridAudit, gridAuditOverlay: $('grid-audit-overlay').value,
    flowColorMode: 'speed',
    boundaryLayer: displayed?.boundaryLayer,
    referenceChord: Number($('reference').value) > 0 ? Number($('reference').value) : 1 });
  drawPressurePanel();
  if (!$('bl-panel').hidden) drawBoundaryLayer(blCanvas, displayed, $('bl-quantity').value);
}
function drawPressurePanel() {
  const pressure = plottedPressure();
  const provisional = Boolean(livePressure || pressure?.coefficientStatus === 'unconverged');
  const mach = livePressure?.mach ?? pressure?.mach;
  const targetMach = livePressure?.targetMach ?? pressure?.machContinuation?.targetMach ?? pressure?.continuation?.targetMach;
  const machLabel = Number.isFinite(mach) ? ` · Mach ${mach.toFixed(3)}`
    + (Number.isFinite(targetMach) && targetMach !== mach ? ` → target ${targetMach.toFixed(3)}` : '') : '';
  const ncrit = quadCoupledNcrit(pressure), ncritLabel = quadCoupledNcritLabel(pressure);
  const grid = quadCoupledGrid(pressure), gridLabel = quadCoupledGridLabel(pressure);
  $('cp-status').hidden = !livePressure && !provisional;
  $('cp-status').textContent = livePressure ? `Live · iteration ${livePressure.iteration}${machLabel}${ncritLabel}${gridLabel} · provisional`
    : provisional ? `Unconverged pressure${machLabel}${ncritLabel}${gridLabel} · provisional` : '';
  document.querySelector('.pressure-panel').classList.toggle('provisional', provisional);
  for (const [key, value] of Object.entries({ iteration: livePressure?.iteration, mach, ncrit: ncrit.actualNcrit,
    gridIntervals: grid.actualGridIntervals })) {
    if (Number.isFinite(value)) pressureCanvas.dataset[key] = String(value);
    else delete pressureCanvas.dataset[key];
  }
  plotPoints = drawPressure(pressureCanvas, pressure, hover, solverBusy ? 'Waiting for the next pressure update…' : undefined);
}
function clearResults() {
  result = null;
  overlayCoefficients = null;
  clearLivePressure();
  $('residual').textContent = '—';
  $('coefficient-warning').hidden = true;
  $('bl-panel').hidden = true;
  $('bl-panel').classList.remove('provisional');
  $('export-button').disabled = true; $('refine-coupled-button').disabled = true;
  $('cp-tooltip').hidden = true; $('cp-hover').textContent = 'Hover to inspect'; hover = null;
}
function clearLivePressure() {
  livePressure = null; hover = null; plotPoints = [];
  $('cp-tooltip').hidden = true; $('cp-hover').textContent = 'Hover to inspect';
}
function showLivePressure(data) {
  clearLivePressure();
  if (data.pressure) {
    livePressure = { ...data.pressure, iteration: data.iteration, mach: data.mach, targetMach: data.targetMach,
      ...quadCoupledNcrit(data), stage: data.stage, gridLevel: data.gridLevel,
      requestedGridIntervals: data.requestedGridIntervals, retained: data.retained };
    $('cp-description').textContent = `${data.pressure.pressureKind} · unvalidated`;
  }
  drawPressurePanel();
}
function markStale() {
  worker?.terminate(); jobId++; stale = true; lines = [];
  activeRequest = null; lastErrorReport = null;
  $('solve-methods').hidden = true;
  pendingMesh = null; meshPhase = null; setBusy(false);
  resetGridAudit();
  $('error-message').hidden = true;
  $('streamline-status').textContent = 'Geometry preview';
  clearResults(); syncModelUI(); setStatus('Inputs changed', 'stale');
  $('run-meta').textContent = isQuadCoupled()?`Simultaneous Euler/grid/BL/wakes · ${quadAutomatic()?'automatic transition':'fixed trips'} · research`:isStreamtubeGrid()?'Moving quadrilateral Euler · research · no boundary layers':isCoupled() ? 'Coupled boundary layers · native Fortran reference' : 'Linear vortex panel baseline';
  try { elements = buildElements(); bounds = assemblyBounds(elements); }
  catch (error) { showError(error.message); elements = []; }
  $('flow-angle').textContent = `${Number($('alpha-number').value).toFixed(2)}°`;
  document.querySelector('.flow-arrow').style.transform = `rotate(${-Number($('alpha-number').value)}deg)`;
  draw();
}
function resetGridAudit() {
  gridAudit = null; auditingGrid = false; $('grid-audit-overlay').value = 'none';
  $('grid-audit-summary').textContent = 'Build mesh, then Check grid. The harmonic comparison tests the SLOR target, including on an unsmoothed mesh.';
  $('grid-audit-details').hidden = true; $('grid-audit-legend').hidden = true;
}
function showGridAudit(audit, elapsed) {
  gridAudit = audit; auditingGrid = false; setBusy(false);
  const tube = Math.max(...audit.harmonic.map(r => r.levels.at(-1).maximum.tubeIntervals));
  const resolved = audit.harmonic.every(r => r.referenceResolved);
  $('grid-audit-summary').textContent = `${audit.status}. Largest spacing jump: ${audit.spacing.maximumAdjacentRatio.toFixed(2)}×; streamline mismatch: ${tube.toFixed(3)} local tube spacings. Reference refinement ${resolved ? 'resolved to the stated target' : 'unresolved'}. Checked in ${(elapsed / 1000).toFixed(2)} s; no flow solve.`;
  $('grid-audit-table').innerHTML = '<table><thead><tr><th>Region</th><th>Streamline mismatch</th><th>Cross-line mismatch</th><th>Last reference change (streamline / cross-line)</th></tr></thead><tbody>' + audit.harmonic.map((r, g) => {
    const m = r.levels.at(-1).maximum;
    return `<tr><td>${g + 1}</td><td>${m.tubeIntervals.toFixed(3)}</td><td>${m.crosslineIntervals.toFixed(3)}</td><td>${r.referenceChange.tubeIntervals.toFixed(3)} / ${r.referenceChange.crosslineIntervals.toFixed(3)}</td></tr>`;
  }).join('') + '</tbody></table><p>Errors are measured in local grid intervals; they are not percentages of aerodynamic error.</p>';
  $('grid-audit-details').hidden = false; $('grid-audit-overlay').value = 'spacing'; $('show-mesh').checked = true;
  setStatus('Grid checked · see report', audit.status.includes('attention') ? 'failed' : 'stale');
  $('run-meta').textContent = 'Initial grid checked independently · physical validation pending'; updateGridAuditOverlay();
}
function updateGridAuditOverlay() {
  const mode = $('grid-audit-overlay').value;
  $('grid-audit-legend').hidden = !gridAudit || mode === 'none';
  $('grid-audit-legend').textContent = mode === 'spacing'
    ? 'Orange: adjacent streamwise spacing ratio > 1.5× or corner sine < 0.1. These are spacing/shear screening targets.'
    : 'Red: independent harmonic-coordinate mismatch > 0.05 local intervals. Darker cells have larger discrepancies; reference refinement may still be unresolved.';
  draw();
}
function errorFormValues() {
  return Object.fromEntries([...document.querySelectorAll('aside input, aside select')].filter(el => el.id)
    .map(el => [el.id, el.type === 'checkbox' ? el.checked : el.value]));
}
function showError(message, event = {}) {
  const details = buildSolverErrorContext({ message, input: activeRequest?.input,
    event, progress: activeRequest?.progress, mesh: event.mesh ?? pendingMesh,
    label: activeRequest?.label ?? caseName });
  lastErrorReport = { ...details.report, reportedAt: new Date().toISOString(),
    request: activeRequest ? { id: activeRequest.id, task: activeRequest.task, submittedAt: activeRequest.submittedAt } : null,
    formValues: structuredClone(activeRequest?.formValues ?? errorFormValues()),
    environment: { url: location.origin + location.pathname, userAgent: navigator.userAgent } };
  $('error-description').textContent = message;
  $('error-context').textContent = details.summary;
  $('error-context').hidden = !details.summary;
  $('error-copy-status').textContent = '';
  $('error-debug-text').hidden = true; $('error-debug-text').value = '';
  $('error-message').hidden = false; setStatus('Check inputs', 'failed');
}
function showCoefficients(next, unconverged = false, conditions = next) {
  overlayCoefficients = { values: next, conditions, provisional: unconverged && [next.cl, next.cm, next.cd].some(Number.isFinite) };
  updateAirfoilOverlay();
  const available = ['cl', 'cm', 'cd'].some(id => Number.isFinite(next[id]));
  $('coefficient-warning').hidden = !unconverged || !available;
  $('coefficient-warning').textContent = 'Unconverged coefficients · amber values are provisional and may change substantially.'
    + (next.coefficientConditions && next.coefficientConditions.mach !== next.mach
      ? ` Last iterate: Mach ${next.coefficientConditions.mach}; requested: ${next.mach}.` : '');
}
function showUnconvergedResult(next, elapsed) {
  worker.terminate(); clearResults(); setBusy(false);
  result = next; stale = true; pendingMesh = next.mesh ?? pendingMesh; meshPhase = 'failed';
  showCoefficients(next, true);
  const d = next.diagnostics, residual = d.equationResidual ?? d.linearResidual;
  $('residual').textContent = Number.isFinite(residual) ? residual.toExponential(1) : '—';
  showError(next.status === 'unconverged'
    ? `Coupled flow did not converge after ${d.iterations} iterations (equation error ${Number.isFinite(residual) ? residual.toExponential(2) : 'unresolved'}): ${d.reason}. Any displayed coefficients are provisional.`
    : next.status === 'outside-model' ? next.warnings.join(' ') : 'The linear system did not meet its residual tolerance. Any displayed coefficients are provisional.', next);
  setStatus(next.status === 'outside-model' ? 'Outside model' : 'Flow unconverged', 'failed');
  $('run-meta').textContent = `${(elapsed / 1000).toFixed(2)} s · ${next.status}`;
  $('streamline-status').textContent = 'Flow unconverged';
  $('export-button').disabled = false;
  draw();
}
function showResult(next, elapsed) {
  result = next; stale = false; pendingMesh = null; meshPhase = null;
  const coupled = Boolean(result.boundaryLayer); const d = result.diagnostics;
  showCoefficients(result);
  $('residual').textContent = (coupled ? d.equationResidual : d.normalVelocityResidual ?? d.nodalStreamfunctionResidual).toExponential(1);
  setStatus(result.warnings.length ? 'Solved · review' : 'Solved');
  $('run-meta').textContent = `${result.panelCount} panels · ${(elapsed/1000).toFixed(2)} s solve · ${coupled ? 'Simultaneous multielement boundary layers' : 'Linear vortex panels'}`;
  $('export-button').disabled = false;
  setBusy(false);
  if(result.model==='multielement-coupled-subcritical'){
    $('run-meta').textContent+=` · Mach ${result.mach.toFixed(2)}`;
  }
  $('bl-panel').hidden = !coupled;
  if (coupled) {
    $('streamline-status').textContent = `Computed wakes (${result.boundaryLayer.wakes.length}) & δ* · true scale`;
    $('bl-legend').innerHTML=result.elements.map((e,i)=>`<span class="legend-item" style="--line-color:${colors[i]}">${escapeHTML(e.name)}</span>`).join('');
  }
  syncResolution();
  draw();
}
function showQuadResult(next, elapsed) {
  worker.terminate(); setBusy(false);
  const d = next.diagnostics, converged = next.status === 'research-converged';
  const residualConverged = d.residualConverged ?? next.flow.residualConverged;
  result = next; stale = !converged; pendingMesh = converged ? null : next.mesh;
  meshPhase = converged ? null : residualConverged ? 'rejected' : 'failed';
  showCoefficients(next, !converged);
  $('residual').textContent = d.equationResidual.toExponential(1);
  $('run-meta').textContent = `Mach ${next.mach.toFixed(2)} · ${d.iterations} Euler iterations · ${(elapsed / 1000).toFixed(2)} s · research, no BL`;
  $('streamline-status').textContent = converged ? 'Euler equations converged · forces unvalidated' : 'Euler flow unconverged';
  $('export-button').disabled = false;
  if (converged) setStatus('Euler converged · research', 'stale');
  else {
    const outcome = residualConverged ? 'Euler equations converged, but the final grid was rejected' : 'Quad Euler did not converge';
    showError(`${outcome} after ${d.iterations} iterations (equation error ${d.equationResidual.toExponential(2)}): ${d.reason} Displayed coefficients are provisional; drag is pressure-only.`, next);
    setStatus(residualConverged ? 'Euler converged · grid rejected' : 'Euler unconverged', 'failed');
  }
  draw();
}
function showQuadCoupledResult(next, elapsed) {
  worker.terminate(); setBusy(false);
  next = quadCoupledGridResult(quadCoupledNcritResult(next, activeRequest?.input ?? lastCase), activeRequest?.input ?? lastCase);
  const grid = quadCoupledGrid(next, activeRequest?.input ?? lastCase);
  const d = next.diagnostics, converged = next.converged;
  const ncrit = quadCoupledNcrit(next), ncritNotReached = ncrit.actualNcrit !== ncrit.targetNcrit
    || next.ncritContinuation?.reachedTarget === false;
  result = next; stale = !converged; pendingMesh = converged ? null : next.mesh;
  quadFlowCache.remember(next);
  meshPhase = converged ? null : 'failed';
  showCoefficients(next, !converged);
  $('coefficient-warning').hidden = false;
  const transonic = next.checkpoint?.restart?.options?.blThermodynamics === 'historical-common-isentrope';
  $('cp-description').textContent = transonic ? 'Physical Euler pressure on the solid contour · research, unvalidated'
    : 'BL-edge pressure from solved edge speed · research, unvalidated';
  const continuation = next.continuation ?? (next.machContinuation ? {
    ...next.machContinuation, currentMach: next.machContinuation.actualMach,
    attempts: next.machContinuation.attempts ?? [],
  } : null);
  const targetNotReached = continuation && !continuation.reachedTarget;
  $('coefficient-warning').textContent = `${next.alphaContinuation?.reachedTarget === false ? `Requested alpha ${next.targetAlpha}° was not reached; these coefficients belong to alpha ${next.actualAlpha}°. ` : ''}${converged ? 'Unvalidated Euler/BL coefficient estimates' : 'Unconverged Euler/BL coefficients · amber values are provisional'}. CL and Cm use pressure on the solid contour; CD uses ${transonic ? 'viscous wake and Euler wave momentum loss' : 'the summed wake momentum loss'}. Physical accuracy remains unvalidated.${targetNotReached ? ` Requested Mach ${continuation.targetMach.toFixed(3)} was not reached; these values belong to Mach ${next.mach.toFixed(3)}.` : ''}${ncritNotReached ? ` Requested Ncrit ${ncrit.targetNcrit} was not reached; these values belong to Ncrit ${ncrit.actualNcrit ?? 'unknown'}.` : ''}${grid.differentGrid ? ` Requested grid (${grid.requestedGridIntervals} intervals/side) was not reached; these values belong to the retained ${grid.actualGridIntervals}-interval grid.` : ''}`;
  $('residual').textContent = d.equationResidual.toExponential(1);
  const attempts = continuation?.attempts ?? next.initialization.attempts ?? [];
  $('run-meta').textContent = `Mach ${next.mach.toFixed(2)} · Re ${next.referenceReynolds.toExponential(2)}${quadCoupledNcritLabel(next)}${quadCoupledGridLabel(next)} · ${next.boundaryLayer.surfaces.length} BLs / ${next.boundaryLayer.wakes.length} wakes · ${(elapsed/1000).toFixed(2)} s · research`;
  $('export-button').disabled = false; $('bl-panel').hidden = false;
  $('bl-panel').classList.toggle('provisional', !converged);
  $('bl-legend').innerHTML = next.elements.map((e,i)=>`<span class="legend-item" style="--line-color:${colors[i]}">${escapeHTML(e.name)}</span>`).join('');
  if(converged)setStatus(next.refinement ? 'Refined Euler/BL converged · research' : 'Euler/BL converged · research','stale');
  else { showError(next.alphaContinuation?.reachedTarget === false
    ? `Requested alpha ${next.targetAlpha}° was not reached: ${next.reason} Showing the retained solution at alpha ${next.actualAlpha}°; its coefficients do not represent the requested angle.`
    : grid.differentGrid
    ? `Requested grid (${grid.requestedGridIntervals} surface intervals per side) was not reached: ${next.reason}. Showing the ${quadCoupledDisplayedStateLabel(next)} on the ${grid.actualGridIntervals}-interval grid; its coefficients do not represent the requested grid.`
    : ncritNotReached
    ? `Requested Ncrit ${ncrit.targetNcrit} was not reached: ${next.reason}. Showing the ${quadCoupledDisplayedStateLabel(next)} at Ncrit ${ncrit.actualNcrit ?? 'unknown'}; its coefficients do not represent the requested condition.`
    : targetNotReached
    ? `Requested Mach ${continuation.targetMach.toFixed(3)} was not reached: ${next.reason}. Showing the ${quadCoupledDisplayedStateLabel(next)} at Mach ${next.mach.toFixed(3)}; its coefficients do not represent the requested condition.`
    : `Coupled Euler/BL did not converge ${next.refinement ? 'during refinement' : `after ${attempts.length} startup attempt(s)`}: ${next.reason}. Displayed BL and pressure profiles are provisional.`, next);setStatus(next.alphaContinuation?.reachedTarget === false?'Target alpha not reached':grid.differentGrid?'Target grid not reached':ncritNotReached?'Target Ncrit not reached':targetNotReached?'Target Mach not reached':'Euler/BL unconverged','failed'); }
  syncCoupledRefinement();
  draw();
}
function run({ meshOnly = false, refineCoupled = false } = {}) {
  if (refineCoupled && !canRefineCoupled()) return;
  const parentResult = refineCoupled ? result : null;
  worker?.terminate(); const id = ++jobId;
  activeRequest = { id, input: null, progress: {}, label: caseName,
    task: refineCoupled ? 'refine-coupled' : meshOnly ? 'mesh' : 'solve',
    submittedAt: new Date().toISOString(), formValues: errorFormValues() };
  $('solve-methods').hidden = true;
  $('solve-method-history').replaceChildren();
  $('solve-methods').querySelector('details').open = false;
  lastErrorReport = null;
  resetGridAudit();
  stale = true; lines = []; clearResults(); $('error-message').hidden = true;
  pendingMesh = parentResult?.mesh ?? null; meshPhase = hasVolumeMesh() ? 'building' : null;

  if (meshOnly || isStreamtubeGrid()) $('show-mesh').checked = true;
  try {
    if(definitions.reduce((n,d)=>n+analysisPanelCount(d,Number($('resolution').value)),0)>700)throw new Error('Use at most 700 surface panels across the assembly. Reduce an element’s panel count.');
    const referenceChord = Number($('reference').value);
    if (!refineCoupled) {
      elements = buildElements(); bounds = assemblyBounds(elements);
      const alpha = Number($('alpha-number').value);
      if (!$('alpha-number').validity.valid || !$('alpha-number').value || !$('reference').validity.valid || !$('reference').value) throw new Error('Enter a valid angle and positive reference chord.');
      lastCase = { elements, alpha, referenceChord, momentReference: { x: referenceChord / 4, y: 0 }, ...flowModelRoute($('flow-model').value, Number($('quad-mach').value)) };
      if (selectedBenchmark) lastCase.geometrySource = { id: selectedBenchmark.id, referenceChord: selectedBenchmark.referenceChord,
        ...structuredClone(selectedBenchmark.provenance) };
      if (isCoupled() && !meshOnly) {
        for (const id of ['reynolds', 'ncrit', 'trip-upper', 'trip-lower']) if (!$(id).value || !$(id).validity.valid) throw new Error('Enter valid Reynolds, Ncrit and transition trip values.');
        Object.assign(lastCase, { reynolds: Number($('reynolds').value), ncrit: Number($('ncrit').value), trips: [Number($('trip-upper').value), Number($('trip-lower').value)] });
      }
      if(isStreamtubeGrid()){
        if(!meshOnly&&(!$('quad-mach').value||!$('quad-mach').validity.valid))throw new Error('Enter a valid freestream Mach number for quad Euler.');
        Object.assign(lastCase,{mach:Number($('quad-mach').value),gridIntervals:Number($('grid-intervals').value),gridTubes:Number($('grid-tubes').value),gridCrosslinePlacement:'potential'});
        const eulerIsmom = Number($('euler-ismom').value);
        if (!Number.isInteger(eulerIsmom) || eulerIsmom < 1 || eulerIsmom > 4) throw new Error('Select a valid Euler formulation (ISMOM 1–4).');
        lastCase.eulerIsmom = eulerIsmom;
        if (!isQuadCoupled()) lastCase.eulerStartup = $('euler-startup').value;
        for (const [id, key] of [['grid-inlet', 'gridInletIntervals'], ['grid-outlet', 'gridOutletIntervals']])
          if ($(id).value !== 'auto') lastCase[key] = Number($(id).value);
        const number = id => {
          if (!$(id).value || !$(id).validity.valid) throw new Error(`Enter a valid ${document.querySelector(`label[for="${id}"]`).textContent}.`);
          return Number($(id).value);
        };
        lastCase.gridChordExponent = number('grid-chord-exponent');
        for (const [id, key] of [['grid-upper-streamlines', 'gridUpperTubes'], ['grid-lower-streamlines', 'gridLowerTubes'], ['grid-gap-streamlines', 'gridGapTubes']])
          if ($(id).value !== '') lastCase[key] = number(id) - 1;
        lastCase.gridSurfaceSpacing = $('grid-surface-spacing').value;
        if (lastCase.gridSurfaceSpacing === 'curvature') lastCase.gridCurvatureSpacing = {
          exponent: number('grid-curvature-exponent'), leadingSpacingRatio: number('grid-le-ratio'), trailingSpacingRatio: number('grid-te-ratio') };
        lastCase.gridEllipticSmoothing = $('grid-elliptic').checked;
        lastCase.gridSmoothingMethod = 'elliptic';
        if ($('grid-match-aspect').checked) lastCase.gridStagnationAspectRatio = number('grid-aspect-ratio');
        lastCase.quadBoundaryLayers = isQuadCoupled();
        if(isQuadCoupled()&&!meshOnly){
          Object.assign(lastCase,{reynolds:number('quad-reynolds'),ncrit:number('quad-ncrit')});
          if(quadAutomatic())lastCase.transitionMode='automatic';
          const globalTrips=[number('quad-trip-upper'),number('quad-trip-lower')];
          lastCase.materialTrips=definitions.map(d=>[d.quadTripUpper??globalTrips[0],d.quadTripLower??globalTrips[1]]);
          if(lastCase.materialTrips.some(pair=>pair.some(v=>!Number.isFinite(v)||v<=0||(quadAutomatic()?v>1:v>=1))))throw new Error(quadAutomatic()?'Every quad trip must be greater than 0 and at most 1.':'Every fixed quad trip must be a surface fraction strictly between 0 and 1.');
        }
      }
    }
    // Run always starts at the requested conditions. Cached roots are not
    // automatic continuation seeds for the direct inviscid → viscous route.
    activeRequest.input = structuredClone(lastCase);
    if (isQuadCoupled() && !meshOnly) activeRequest.methodProgress = createQuadSolveProgress(lastCase, parentResult);
    setBusy(true, meshOnly);
    $('streamline-status').textContent = refineCoupled ? 'Refining solved grid…' : hasVolumeMesh() ? 'Building requested mesh…' : 'Computing flow field…';
    $('run-meta').textContent = refineCoupled ? 'Transferring the coupled flow to a finer grid'
      : isQuadCoupled() && !meshOnly ? 'Step 1/2: inviscid at requested grid, alpha and Mach'
      : hasVolumeMesh() ? 'Preparing geometry and initial wake guides' : 'Initializing flow';
    setStatus(refineCoupled ? 'Refining coupled solution' : meshOnly ? 'Building mesh' : 'Solving');
    worker = new Worker(new URL('../worker/solver.js', import.meta.url), { type: 'module' });
    worker.onmessage = ({ data }) => {
      if (data.id !== jobId) return;
      if (activeRequest?.methodProgress) {
        const event = data.type === 'iteration' ? data.iteration : ['flow-stage', 'pressure', 'pressure-unavailable', 'coefficients', 'coefficients-unavailable'].includes(data.type) ? data : null;
        if (event) {
          const progress = activeRequest.methodProgress.update(event, { stageChange: data.type === 'flow-stage' });
          const labels = JSON.stringify(progress);
          if (progress.current && labels !== activeRequest.methodLabels) {
            activeRequest.methodLabels = labels;
            $('solve-methods').hidden = false;
            $('solve-method-current').textContent = `Method: ${progress.current}`;
            $('solve-method-summary').textContent = `Methods attempted (${progress.history.length})`;
            $('solve-method-history').replaceChildren(...progress.history.map(label => {
              const item = document.createElement('li'); item.textContent = label; return item;
            }));
          }
        }
        if (['result', 'error'].includes(data.type)) $('solve-method-current').textContent = $('solve-method-current').textContent.replace(/^Method:/, 'Last method:');
      }
      if (activeRequest?.id === data.id) {
        const progress = data.type === 'iteration' ? data.iteration : data.type === 'flow-stage' ? data : null;
        if (progress) {
          const sameStage = progress.stage === activeRequest.progress.stage
            && (progress.actualNcrit === undefined || progress.actualNcrit === activeRequest.progress.actualNcrit);
          activeRequest.progress = { ...(sameStage ? activeRequest.progress : {}), ...progress };
          if (Number.isFinite(progress.mach)) activeRequest.progress.actualMach = progress.actualMach ?? progress.mach;
        }
      }
      if(data.type==='flow-stage'){
        clearLivePressure(); drawPressurePanel();
        if (Object.keys(quadCoupledNcrit(data)).length || data.stage==='coupled-grid-refinement') showCoefficients({ cl: null, cd: null, cm: null }, true, data);
        const label=data.stage==='coupled-grid-refinement'?`${data.retained?'Retained':'Preparing'} Euler/BL grid ${data.gridLevel ?? 'refinement'} → ${data.requestedGridIntervals ?? lastCase.gridIntervals} intervals/side`:data.stage?.includes('ncrit')?'Preparing automatic transition startup':data.stage==='coupled-coarse-initialization'?'Preparing a coarse-grid BL starting solution':data.stage==='transition-refinement'?'Refining the boundary-layer grid':data.stage==='coupled-wake-grid'?'Rebuilding wake grid':data.stage==='coupled-refinement'?'Refining solved grid':data.stage==='hybrid-certification'?'Preparing saved flow':data.stage==='coupled-cold-recovery'?'Restarting from lower Mach':data.stage==='euler-cold-startup'?'Initializing existing mesh at lower Mach':data.stage==='coupled-mach'?'Advancing Mach':data.stage==='euler'?'Solving Euler precursor':data.stage==='boundary-layer-initialization'
          ?`Initializing ${elements.length*2} boundary layers and ${elements.length} wakes`:refineCoupled?'Solving refined Euler/BL':`Solving Euler/BL · attempt ${data.startupAttempt}`;
        const machLabel = Number.isFinite(data.mach) && Number.isFinite(data.targetMach)
          ? ` · Mach ${data.mach.toFixed(3)} → ${data.targetMach.toFixed(3)}` : '';
        $('run-meta').textContent=label+machLabel+quadCoupledNcritLabel(data)+(data.retryReason?' · recovering with a thinner initial guess':'');
        $('status').dataset.stage=data.stage;
      }
      if (data.type === 'mesh') {
        pendingMesh = data.mesh; meshPhase = data.stage;
        $('streamline-status').textContent = data.mesh.iteration ? `${quadIterationLabel(data.mesh)} iteration ${data.mesh.iteration.iteration} · moving grid · solving`
          : data.stage === 'smoothing' ? 'Relaxing the grid · no Euler iteration' : isExperimentalMesh(data.mesh) ? `${experimentalMeshTitle(data.mesh)} · no flow solution` : meshOnly ? 'Initial mesh · no flow solution' : 'Mesh available · computing flow';
        draw();
      }
      if (data.type === 'mesh-ready') {
        worker.terminate(); setBusy(false); const partialSlor = pendingMesh?.initialization?.gridSmoothing?.initialGuessAccepted === true && pendingMesh.initialization.gridSmoothing.converged === false; setStatus(isExperimentalMesh(pendingMesh) ? experimentalMeshTitle(pendingMesh) : partialSlor ? 'Mesh ready · SLOR incomplete' : 'Mesh ready', isExperimentalMesh(pendingMesh) || partialSlor ? 'stale' : '');
        $('run-meta').textContent = `${pendingMesh.cells.length.toLocaleString()} cells · ${(data.elapsed / 1000).toFixed(2)} s · ${isExperimentalMesh(pendingMesh) ? 'experimental mesh, physical acceptance pending' : 'initial mesh, flow not solved'}`;
        draw();
      }
      if (data.type === 'mesh-unavailable') { meshPhase = 'waiting'; draw(); }
      if (data.type === 'accepted-flow') quadFlowCache.remember(data.parentResult);
      if (data.type === 'pressure' && isQuadCoupled()) showLivePressure(data);
      if (data.type === 'pressure-unavailable') { clearLivePressure(); drawPressurePanel(); }
      if (data.type === 'coefficients-unavailable') {
        clearLivePressure(); drawPressurePanel();
        showCoefficients({ cl: null, cd: null, cm: null }, true, data);
        $('coefficient-warning').hidden = false;
        $('coefficient-warning').textContent = `Current coefficients unavailable at Mach ${data.mach.toFixed(3)}${quadCoupledNcritLabel(data)}: ${data.message}. The flow solve continues.`;
      }
      if (data.type === 'coefficients' && isStreamtubeGrid()) {
        const eulerPressure = data.kind === 'euler-pressure';
        showCoefficients(data.coefficients, true, data);
        showLivePressure(data);
        $('coefficient-warning').hidden = false;
        $('coefficient-warning').textContent = eulerPressure
          ? 'Current Euler iterate · provisional pressure loads. CD excludes viscous drag.'
          : 'Current Euler/BL iterate · unconverged estimates. Amber values may change substantially.';
        if (Number.isFinite(data.mach) && Number.isFinite(data.targetMach) && data.mach !== data.targetMach) {
          $('coefficient-warning').hidden = false;
          $('coefficient-warning').textContent += ` Current Mach ${data.mach.toFixed(3)} → target ${data.targetMach.toFixed(3)}.`;
        }
        if (!Number.isFinite(data.coefficients.cd)) {
          $('coefficient-warning').textContent += ' CD is unavailable for the current wake state; CL and Cm remain provisional.';
        }
        if (Object.keys(quadCoupledNcrit(data)).length) {
          $('coefficient-warning').hidden = false;
          $('coefficient-warning').textContent += `${quadCoupledNcritLabel(data)}. These values belong to the current startup condition.`;
        }
      }
      if (data.type === 'iteration') {
        const h=data.iteration;
        if(isQuadCoupled()&&(['euler','coupled','coupled-mach','coupled-wake-grid','coupled-coarse-initialization','coupled-grid-refinement'].includes(h.stage)||h.stage?.includes('ncrit'))){
          const machLabel = (Number.isFinite(h.mach) && Number.isFinite(h.targetMach) ? ` · Mach ${h.mach.toFixed(3)} → ${h.targetMach.toFixed(3)}` : '')+quadCoupledNcritLabel(h);
          $('run-meta').textContent=h.stage==='euler'?`Euler precursor${machLabel} · iteration ${h.iteration} · error ${h.residual.toExponential(2)}`
            :h.stage==='coupled-grid-refinement'?`Euler/BL grid ${h.gridLevel ?? 'refinement'} → ${h.requestedGridIntervals ?? lastCase.gridIntervals} intervals/side${machLabel} · iteration ${h.iteration} · error ${h.residual.toExponential(2)}`
            :h.stage==='coupled-coarse-initialization'||h.stage?.includes('ncrit')?`${(h.gridLevel??h.ncritContinuation?.phase)==='fine'?'Fine-grid ':(h.gridLevel??h.ncritContinuation?.phase)==='coarse'||h.stage==='coupled-coarse-initialization'?'Coarse-grid ':''}${h.coarseStage==='euler'?'Euler':'Euler/BL'} startup${machLabel} · iteration ${h.iteration} · error ${h.residual.toExponential(2)}`
            :h.stage==='coupled-wake-grid'?`Rebuilding wake grid${machLabel} · iteration ${h.iteration} · error ${h.residual.toExponential(2)}`
            :`${h.ncritContinuation?.phase==='fine'?'Fine-grid Euler/BL startup':refineCoupled || h.transitionRecovery ? 'Refined Euler/BL' : `Euler/BL · attempt ${h.startupAttempt ?? h.attempt ?? 1}`}${machLabel} · iteration ${h.iteration} · error ${h.residual.toExponential(2)}`;
          $('residual').textContent=h.residual.toExponential(1);
          return;
        }
        if(h.stage==='quad-euler'){
          const retry = h.startupStrategy === 'broader-shock-retry' ? ' · restarting with broader shock smoothing (same grid)' : '';
          const order = h.dissipation?.targetMucon > 0
            ? h.dissipation.mucon < 0 ? ' · first-order shock startup' : ' · requested second-order dissipation' : '';
          const dissipation = h.dissipation ? `${order} · MCRIT ${h.dissipation.mcrit.toFixed(3)}`
            + (h.dissipation.mcrit !== h.dissipation.targetMcrit ? ` → ${h.dissipation.targetMcrit.toFixed(3)} (shock smoothing)` : '') : '';
          $('run-meta').textContent = `Quad Euler · Mach ${h.mach.toFixed(2)} · iteration ${h.iteration} · equation error ${h.residual.toExponential(2)}${retry}${dissipation}`;
          return;
        }
        const stage=h.stage==='assembly-initialization'?'Incompressible initialization'
          :h.stage==='outer-initialization'||h.stage==='incompressible-initialization'||h.edgeVelocityFraction<1?'Flow initialization':null;
        const progress=h.kind==='remesh'?'rebuilding mesh':h.kind==='wake'?`normal velocity ${h.wakeResidual.toExponential(2)}`
          :`Newton ${h.iteration} · equation error ${h.residual.toExponential(2)}`;
        $('run-meta').textContent = [stage,Number.isFinite(h.continuationMach)?`Mach ${h.continuationMach.toFixed(2)}`:null,
          `Wake ${(h.wakeIteration??0)+1}`,progress].filter(Boolean).join(' · ');
      }
      if (data.type === 'result') {
        clearLivePressure();
        if(data.result.model==='research-streamtube-euler-bl'){showQuadCoupledResult(data.result,data.elapsed);return;}
        if(data.result.model==='research-streamtube-euler'){showQuadResult(data.result,data.elapsed);return;}
        if (data.result.status !== 'solved') {
          showUnconvergedResult(data.result, data.elapsed); return;
        }
        showResult(data.result, data.elapsed);
      }
      if (data.type === 'streamlines') { lines = data.lines; $('streamline-status').textContent = 'Potential-flow streamlines'; draw(); }
      if (data.type === 'streamline-error') { $('streamline-status').textContent = `Flow traces unavailable: ${data.message}`; }
      if (data.type === 'error') {
        worker.terminate(); stale = true; meshPhase = 'failed'; clearResults(); showError(data.message, data); setBusy(false);
        if(isStreamtubeGrid())setStatus(data.code==='streamtube-sonic-capacity'?'Euler initialization blocked':isQuadCoupled()?'Euler/BL solve failed':'Quad Euler / grid failed','failed');
        $('streamline-status').textContent = 'No flow solution'; draw();
      }
    };
    worker.onerror = event => { if (id !== jobId) return; worker.terminate(); stale = true; meshPhase = 'failed'; clearResults(); showError(event.message || 'Solver worker failed.', { stage: 'worker', diagnostics: { filename: event.filename, lineno: event.lineno, colno: event.colno } }); setBusy(false); draw(); };
    worker.postMessage({ id, caseData: lastCase, meshOnly, ...(refineCoupled ? { task: 'refine-coupled', parentResult }
      : {}), bounds: {
      xMin: bounds.xMin - 0.5 * referenceChord, xMax: bounds.xMax + 0.5 * referenceChord,
      yMin: bounds.yMin - 0.2 * referenceChord, yMax: bounds.yMax + 0.2 * referenceChord,
    } });
    draw();
  } catch (error) { meshPhase = 'failed'; showError(error.message, error); setBusy(false); draw(); }
}

$('element-list').addEventListener('input', event => {
  const { element, key } = event.target.dataset;
  if (element === undefined) return;
  if(['tripUpper','tripLower','quadTripUpper','quadTripLower','panels'].includes(key)&&event.target.value==='')delete definitions[Number(element)][key];
  else definitions[Number(element)][key] = key === 'code' ? event.target.value : (event.target.value === '' ? NaN : Number(event.target.value));
  if (key === 'code') event.target.closest('details').querySelector('summary small').textContent = `NACA ${event.target.value}`;
  // Benchmark coordinates are read-only; their editable trips are flow inputs.
  if (!selectedBenchmark) { $('preset').value = 'custom'; caseName = 'Custom assembly'; }
  if(key==='panels')syncResolution();markStale();
});
$('element-list').addEventListener('click', event => {
  if (event.target.dataset.remove === undefined) return;
  selectedBenchmark = null;
  definitions.splice(Number(event.target.dataset.remove), 1); $('preset').value = 'custom'; caseName = 'Custom assembly'; renderElements(); markStale();
});
$('preset').addEventListener('change', () => {
  coordinateDocument = null; $('coordinate-import-note').hidden = true;
  selectedBenchmark = getBenchmarkAirfoil($('preset').value);
  definitions = selectedBenchmark ? selectedBenchmark.elements : structuredClone(presets[$('preset').value]);
  if (selectedBenchmark) $('reference').value = String(selectedBenchmark.referenceChord);
  if (['single', 'flap', 'three'].includes($('preset').value)) {
    for (const id of ['quad-mach', 'quad-reynolds', 'reynolds', 'quad-ncrit', 'ncrit',
      'alpha-number', 'alpha', 'trip-upper', 'trip-lower']) $(id).value = $(id).defaultValue;
    for (const settings of Object.values(quadSettings)) settings['quad-mach'] = $('quad-mach').defaultValue;
    quadTripsByMode.automatic = ['1', '1'];
    quadTripsByMode['fixed-trip'] = ['0.05', '0.05'];
    ['quad-trip-upper', 'quad-trip-lower'].forEach((id, k) => {
      $(id).value = quadTripsByMode[previousQuadTransition][k];
    });
  }
  if (selectedBenchmark?.id === 'rae2822-mses') {
    for (const [id, value] of Object.entries({
      'quad-mach': '0.74', 'quad-reynolds': '2700000', 'reynolds': '2700000',
      'quad-ncrit': '4', 'ncrit': '4',
      'alpha-number': '2.68', 'alpha': '2.68',
      'quad-trip-upper': '0.03', 'quad-trip-lower': '0.07',
      'trip-upper': '0.03', 'trip-lower': '0.07',
    })) $(id).value = value;
    for (const settings of Object.values(quadSettings)) settings['quad-mach'] = '0.74';
    quadTripsByMode[previousQuadTransition] = ['0.03', '0.07'];
  }
  if (selectedBenchmark?.id === 'nlr7301') {
    for (const [id, value] of Object.entries({
      'quad-mach': '0.185', 'quad-reynolds': '2510000', 'reynolds': '2510000',
      'alpha-number': '6', 'alpha': '6',
    })) $(id).value = value;
    for (const settings of Object.values(quadSettings)) settings['quad-mach'] = '0.185';
  }
  caseName = selectedBenchmark?.name ?? $('preset').selectedOptions[0].text;
  zoom = 1; pan = { x: 0, y: 0 }; renderElements(); markStale();
  if(!isStreamtubeGrid())run();
});
$('add-element').addEventListener('click', () => {
  if (definitions.length >= 6) return;
  selectedBenchmark = null;
  definitions.push({ name: `Element ${definitions.length + 1}`, code: '0012', chord: 0.2, x: 1.4 + (definitions.length - 2) * 0.3, y: -0.1, deflection: 0 });
  $('preset').value = 'custom'; caseName = 'Custom assembly'; renderElements(); markStale();
});
$('alpha').addEventListener('input', () => { $('alpha-number').value = $('alpha').value; markStale(); });
$('alpha-number').addEventListener('input', () => { $('alpha').value = $('alpha-number').value; markStale(); });
$('flow-model').addEventListener('change', () => {
  if (quadSettings[previousFlowMode]) quadSettings[previousFlowMode] = readQuadSettings();
  previousFlowMode = $('flow-model').value;
  if (quadSettings[previousFlowMode]) for (const [id, value] of Object.entries(quadSettings[previousFlowMode])) {
    if ($(id).type === 'checkbox') $(id).checked = value;
    else $(id).value = value;
  }
  markStale();
});
$('quad-transition').addEventListener('change',()=>{
  quadTripsByMode[previousQuadTransition]=['quad-trip-upper','quad-trip-lower'].map(id=>$(id).value);
  previousQuadTransition=$('quad-transition').value;
  ['quad-trip-upper','quad-trip-lower'].forEach((id,k)=>{$(id).value=quadTripsByMode[previousQuadTransition][k];});
  markStale();
});
for (const id of ['reynolds', 'ncrit', 'trip-upper', 'trip-lower','quad-mach','quad-reynolds','quad-ncrit','quad-trip-upper','quad-trip-lower']) $(id).addEventListener('input', markStale);
for(const id of ['euler-ismom','euler-startup','grid-intervals','grid-tubes','grid-inlet','grid-outlet'])$(id).addEventListener('change',markStale);
for (const id of ['grid-chord-exponent','grid-upper-streamlines','grid-lower-streamlines','grid-gap-streamlines','grid-le-ratio','grid-te-ratio','grid-curvature-exponent','grid-aspect-ratio']) $(id).addEventListener('input',markStale);
for (const id of ['grid-elliptic','grid-match-aspect']) $(id).addEventListener('change',markStale);
$('grid-surface-spacing').addEventListener('change', markStale);
$('bl-quantity').addEventListener('change', draw);
$('reference').addEventListener('input', markStale); $('resolution').addEventListener('change',()=>{syncResolution();markStale();});
$('solve-button').addEventListener('click', run);
$('refine-coupled-button').addEventListener('click', () => run({ refineCoupled: true }));
$('grid-audit-overlay').addEventListener('change', updateGridAuditOverlay);
$('grid-audit-export').addEventListener('click', () => {
  const mesh = stale ? pendingMesh : result?.mesh; if (!mesh || (!gridAudit && !isExperimentalMesh(mesh))) return;
  const payload = isExperimentalMesh(mesh)
    ? { schemaVersion: 1, application: 'alula experimental grid preview', exportedAt: new Date().toISOString(), input: lastCase, mesh, smoothing: mesh.initialization.gridSmoothing, physicalAcceptance: false }
    : { schemaVersion: 1, application: 'alula initial-grid audit', exportedAt: new Date().toISOString(), mesh, audit: gridAudit };
  const url = URL.createObjectURL(new Blob([JSON.stringify(payload)], { type: 'application/json' }));
  const link = document.createElement('a'); link.href = url; link.download = isExperimentalMesh(mesh) ? 'alula-experimental-grid.json' : 'alula-grid-audit.json'; link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
});
$('stop-button').addEventListener('click', () => {
  if (auditingGrid) {
    worker?.terminate(); jobId++; auditingGrid = false; setBusy(false);
    $('grid-audit-summary').textContent = 'Grid check stopped; no verification result.'; setStatus('Grid check stopped'); draw(); return;
  }
  worker?.terminate(); jobId++; stale = true; meshPhase = 'stopped'; lines = [];
  clearResults(); setBusy(false); setStatus('Stopped');
  $('solve-method-current').textContent = $('solve-method-current').textContent.replace(/^Method:/, 'Last method:');
  $('run-meta').textContent = 'Analysis stopped · coefficients withheld';
  $('streamline-status').textContent = pendingMesh ? 'Mesh only · no flow solution' : 'Geometry preview'; draw();
});
$('show-panels').addEventListener('change', draw); $('show-streamlines').addEventListener('change', draw);
$('show-mesh').addEventListener('change', draw);
$('show-mach-contours').addEventListener('change', draw);
$('reset-view').addEventListener('click', () => { zoom = 1; pan = { x: 0, y: 0 }; draw(); });
$('import-file').addEventListener('change', async event => {
  const file = event.target.files[0]; if (!file) return;
  try {
    if (file.size > 2e6) throw new Error('Coordinate files must be smaller than 2 MB.');
    const parsed = parseCoordinates(await file.text());
    if (parsed.elements.length > 6) throw new Error('The workbench supports up to six elements.');
    selectedBenchmark = null;
    definitions = parsed.elements.map(prepareAirfoilElement);
    coordinateDocument = { ...parsed, elements: undefined };
    $('coordinate-import-note').textContent = [...parsed.warnings, ...(parsed.header ? ['The imported domain/header is retained for coordinate export; it does not override the selected solver mesh settings.'] : [])].join(' ');
    $('coordinate-import-note').hidden = !$('coordinate-import-note').textContent;
    $('preset').value = 'custom'; caseName = parsed.name;
    zoom = 1; pan = { x: 0, y: 0 }; renderElements(); markStale(); if(!isStreamtubeGrid())run();
  } catch (error) { markStale(); showError(error.message); }
  event.target.value = '';
});
$('coordinate-format').addEventListener('change', refreshCoordinateExport);
$('coordinate-export-button').addEventListener('click', exportCoordinates);
$('copy-error-debug').addEventListener('click', async () => {
  if (!lastErrorReport) return;
  const report = lastErrorReport;
  const text = JSON.stringify(report, (_, value) => ArrayBuffer.isView(value) ? Array.from(value) : value, 2);
  try {
    if (!navigator.clipboard?.writeText) throw new Error('Clipboard unavailable.');
    await navigator.clipboard.writeText(text);
    if (lastErrorReport === report) $('error-copy-status').textContent = 'Diagnostics copied';
  } catch {
    if (lastErrorReport !== report) return;
    $('error-debug-text').value = text; $('error-debug-text').hidden = false;
    $('error-debug-text').focus(); $('error-debug-text').select();
    $('error-copy-status').textContent = 'Select and copy the diagnostic text below.';
  }
});
$('export-button').addEventListener('click', () => {
  if (!result || $('export-button').disabled) return;
  const { field, ...output } = result;
  const payload = { schemaVersion: 3, application: 'alula 0.5.0-preview.1', exportedAt: new Date().toISOString(),
    input: lastCase, result: output, ...(lastErrorReport ? { failureDiagnostics: lastErrorReport } : {}),
    limitations: result.limitations ?? 'Incompressible inviscid baseline. CD is unavailable. No boundary layers, transition, wakes, separation or shocks.' };
  const url = URL.createObjectURL(new Blob([JSON.stringify(payload,
    (_, value) => ArrayBuffer.isView(value) && typeof value.length === 'number' ? Array.from(value) : value, 2)], { type: 'application/json' }));
  const link = document.createElement('a'); link.href = url; link.download = 'alula-result.json'; link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
});
bindPlotNavigation(geometryCanvas, { getView: () => ({ zoom, pan }),
  onSwipe: direction => workspaceTabs.swipe(direction),
  swipeEnabled: () => workspaceTabs.mobile.matches && zoom === 1,
  setView: view => { zoom = view.zoom; pan = view.pan; draw(); } });
pressureCanvas.addEventListener('pointermove', event => {
  const pressure = plottedPressure();
  if (!pressure) return;
  const box = pressureCanvas.getBoundingClientRect(); const x = event.clientX - box.left; const y = event.clientY - box.top;
  hover = plotPoints.reduce((best, p) => Math.hypot(p.screenX - x, p.screenY - y) < Math.hypot(best.screenX - x, best.screenY - y) ? p : best, plotPoints[0]);
  if (!hover) return;
  const tooltip = $('cp-tooltip'); tooltip.hidden = false;
  tooltip.textContent = `${hover.name} · x/c ${(hover.x / pressure.referenceChord).toFixed(3)} · Cp ${hover.cp.toFixed(3)}`;
  tooltip.style.left = `${Math.max(4, Math.min(x + 12, box.width - tooltip.offsetWidth - 5))}px`;
  tooltip.style.top = `${Math.max(4, y - 35)}px`;
  $('cp-hover').textContent = `Cp ${hover.cp.toFixed(3)}`; drawPressure(pressureCanvas, pressure, hover);
});
pressureCanvas.addEventListener('pointerleave', () => { hover = null; $('cp-tooltip').hidden = true; $('cp-hover').textContent = 'Hover to inspect'; drawPressurePanel(); });
const observer = new ResizeObserver(draw); observer.observe(geometryCanvas); observer.observe(pressureCanvas); observer.observe(blCanvas);
renderElements(); markStale(); run();
