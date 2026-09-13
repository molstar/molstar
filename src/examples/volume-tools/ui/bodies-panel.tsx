/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import * as React from 'react';
import { useCallback, useEffect, useRef, useState } from 'react';
import * as ReactDOM from 'react-dom';
import { Volume } from '../../../mol-model/volume';
import { OpenFiles } from '../../../mol-plugin-state/actions/file';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../../mol-plugin/context';
import { Asset } from '../../../mol-util/assets';
import { Color } from '../../../mol-util/color';
import { UUID } from '../../../mol-util/uuid';
import { isBodyMaskCell, VolumeSegmentorManager, VolumeSegmentorState } from '../../../extensions/volume-tools/segmentor';
import type { BodyInfo, Point2D } from '../../../extensions/volume-tools/segmentor';
import { DrawingCanvas } from './drawing-canvas';
import { SectionDivider, btnStyle, hintStyle, labelStyle, numInputStyle, rowStyle, smallBtnStyle } from './panel-styles';

interface Props {
    plugin: PluginContext;
    manager: VolumeSegmentorManager;
}

const DEFERRED_THRESHOLD_PREVIEW_DIMENSION = 280;
const FALLBACK_COLOR = '#FF6B00';

export function BodiesPanel({ plugin, manager }: Props) {
    const [state, setState] = useState<VolumeSegmentorState>(manager.state);
    const [previewViewId, setPreviewViewId] = useState<string | null>(null);
    const [busy, setBusy] = useState<string | null>(null);
    const [dustMinVoxels, setDustMinVoxels] = useState(100);
    const [dustResult, setDustResult] = useState<number | null>(null);
    const [exportOptions, setExportOptions] = useState({ bundleZip: true, compress: false });
    const [thresholdDraft, setThresholdDraft] = useState<number | null>(null);
    const [, forceUpdate] = useState(0);
    const volumeInputRef = useRef<HTMLInputElement>(null);

    useEffect(() => {
        const sub = manager.behaviors.state.subscribe(setState);
        return () => sub.unsubscribe();
    }, [manager]);

    // Re-render when plugin state changes (e.g. visibility toggled via Molstar UI)
    useEffect(() => {
        const sub = plugin.state.data.events.changed.subscribe(() => forceUpdate(n => n + 1));
        return () => sub.unsubscribe();
    }, [plugin]);

    // Auto-select the first (source) volume on load
    useEffect(() => {
        const sub = plugin.managers.volume.hierarchy.behaviors.selection.subscribe(() => {
            const volumes = plugin.managers.volume.hierarchy.current.volumes.filter(v => !isBodyMaskCell(v.cell));
            const target = manager.state.targetVolumeRef;
            const targetExists = !!target && volumes.some(v => v.cell.transform.ref === target);
            if (volumes.length > 0 && !targetExists) void manager.setTargetVolume(volumes[0].cell.transform.ref);
            else if (volumes.length === 0 && target) void manager.setTargetVolume(undefined);
        });
        return () => sub.unsubscribe();
    }, [plugin, manager]);

    // Ctrl/Cmd+Z undoes the last body change unless a text field has focus
    useEffect(() => {
        const onKeyDown = (e: KeyboardEvent) => {
            if (!(e.ctrlKey || e.metaKey) || e.key.toLowerCase() !== 'z' || e.shiftKey) return;
            const target = e.target as HTMLElement | null;
            if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA')) return;
            if (manager.state.busy || manager.state.undoDepth === 0) return;
            e.preventDefault();
            void manager.undo();
        };
        window.addEventListener('keydown', onKeyDown);
        return () => window.removeEventListener('keydown', onKeyDown);
    }, [manager]);

    // Forget the previewed view when the active body changes
    useEffect(() => { setPreviewViewId(null); }, [state.activeBodyId]);

    const run = useCallback(async (key: string, action: () => Promise<unknown>) => {
        setBusy(key);
        try { await action(); } catch (e) { alert(String(e)); } finally { setBusy(null); }
    }, []);

    const handleOpenVolume = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        e.target.value = '';
        await plugin.clear();
        await manager.setTargetVolume(undefined);
        try {
            await plugin.runTask(plugin.state.data.applyAction(OpenFiles, {
                files: [Asset.File(file)],
                format: { name: 'auto', params: {} },
                visuals: true,
            }));
        } catch (err) {
            alert(String(err));
        }
    }, [plugin, manager]);

    const volumeObj = state.targetVolumeRef
        ? plugin.state.data.select(state.targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined
        : undefined;
    const volumeLabel = volumeObj?.label ?? 'Volume';
    const volumeHidden = state.targetVolumeRef ? !!plugin.state.data.cells.get(state.targetVolumeRef)?.state.isHidden : false;
    const volStats = volumeObj?.data.grid.stats;
    const dimensions = volumeObj?.data.grid.cells.space.dimensions as [number, number, number] | undefined;
    const deferThresholdSliderPreview = !!dimensions && Math.max(...dimensions) > DEFERRED_THRESHOLD_PREVIEW_DIMENSION;

    const thresholdAbs = volStats ? Volume.IsoValue.toAbsolute(state.threshold, volStats).absoluteValue : null;
    const thresholdShown = thresholdDraft ?? thresholdAbs;
    const thresholdSigma = (volStats && thresholdShown !== null && volStats.sigma > 0) ? (thresholdShown - volStats.mean) / volStats.sigma : null;

    const activeBody = state.bodies.find(b => b.id === state.activeBodyId);
    const activeColor = activeBody ? Color.toHexStyle(activeBody.color) : FALLBACK_COLOR;
    const previewView = previewViewId ? activeBody?.views.find(v => v.id === previewViewId) : undefined;
    const exportable = state.bodies.filter(b => b.voxelCount > 0).length;
    const hasRemainder = state.bodies.some(b => b.remainder);
    const anyBusy = !!busy || state.busy;

    const viewportEl = (
        (plugin.layout.root as HTMLElement | undefined)?.querySelector('.msp-viewport') ??
        document.querySelector('.msp-viewport')
    ) as HTMLElement | null;

    const handlePolygonComplete = useCallback((polygon: Point2D[], canvasW: number, canvasH: number) => {
        const cam = plugin.canvas3d?.camera;
        const body = manager.activeBody;
        if (!cam || !body) return;
        const vp = cam.viewport;
        manager.setDrawing(false);
        void manager.addView(body.id, {
            id: UUID.create22(),
            label: `View ${body.views.length + 1}`,
            polygon, canvasWidth: canvasW, canvasHeight: canvasH,
            viewportWidth: vp.width, viewportHeight: vp.height,
            cameraSnapshot: cam.getSnapshot(),
        });
    }, [plugin, manager]);

    const commitThreshold = useCallback((value: number) => {
        setThresholdDraft(null);
        manager.setThreshold(Volume.IsoValue.absolute(value));
    }, [manager]);

    const handleThresholdSliderChange = useCallback((value: number) => {
        if (deferThresholdSliderPreview) setThresholdDraft(value);
        else commitThreshold(value);
    }, [deferThresholdSliderPreview, commitThreshold]);

    const flushThresholdDraft = useCallback(() => {
        if (thresholdDraft !== null) commitThreshold(thresholdDraft);
    }, [thresholdDraft, commitThreshold]);

    return (
        <div style={{ padding: '8px', fontFamily: 'sans-serif', fontSize: '13px', color: '#ccc' }}>

            {/* ── Data ─────────────────────────────────────────── */}
            <input ref={volumeInputRef} type='file'
                accept='.mrc,.ccp4,.map,.dsn6,.brix,.dx,.dxbin,.cub,.cube,.gz'
                style={{ display: 'none' }} onChange={handleOpenVolume} />

            <div style={{ display: 'flex', gap: '4px', marginBottom: '2px' }}>
                <button style={{ ...btnStyle, flex: 1, background: '#1e5230' }} onClick={() => volumeInputRef.current?.click()}>
                    Open Volume
                </button>
                {state.targetVolumeRef && (
                    <button style={{ ...btnStyle, background: volumeHidden ? '#444' : '#225' }}
                        title={volumeHidden ? 'Show volume' : 'Hide volume'}
                        onClick={() => setSubtreeVisibility(plugin.state.data, state.targetVolumeRef!, !volumeHidden)}>
                        {volumeHidden ? '○' : '●'}
                    </button>
                )}
            </div>
            {state.targetVolumeRef && (
                <div style={{ fontSize: '11px', marginBottom: '8px', color: '#7c7', wordBreak: 'break-all' }}>
                    {volumeLabel}{dimensions && <span style={{ color: '#666' }}> · {dimensions.join('×')}</span>}
                </div>
            )}

            {state.targetVolumeRef && volStats && (<>
                {/* ── Source volume ─────────────────────────────── */}
                <SectionDivider label='Source volume' />

                <label style={labelStyle}>
                    Threshold:&nbsp;
                    <input type='number' step={volStats.sigma / 20 || 0.001} style={numInputStyle}
                        value={thresholdShown !== null ? parseFloat(thresholdShown.toPrecision(5)) : ''}
                        onChange={e => { const v = parseFloat(e.target.value); if (!isNaN(v)) commitThreshold(v); }} />
                    {thresholdSigma !== null && <span style={{ color: '#888' }}>&nbsp;({thresholdSigma.toFixed(1)}σ)</span>}
                </label>
                <input type='range' min={volStats.min} max={volStats.max} step={volStats.sigma / 20 || 0.001} style={{ width: '100%' }}
                    value={thresholdShown ?? volStats.mean}
                    onChange={e => handleThresholdSliderChange(Number(e.target.value))}
                    onPointerUp={flushThresholdDraft}
                    onKeyUp={flushThresholdDraft} />
                {deferThresholdSliderPreview && (
                    <p style={hintStyle}>Large volume: the slider updates the numeric threshold immediately and refreshes the surface on release.</p>
                )}

                <div style={{ display: 'flex', gap: '4px', margin: '8px 0 4px' }}>
                    <button style={{ ...btnStyle, flex: 1, background: '#1a1a2e' }}
                        disabled={anyBusy}
                        title='Mirror the source volume along X (cannot be undone)'
                        onClick={() => void run('flip', () => manager.flipHandedness())}>
                        {busy === 'flip' ? '⏳ Flipping…' : '↔ Flip Handedness'}
                    </button>
                    <button style={{ ...btnStyle, flex: 1, background: '#152a1a' }}
                        disabled={anyBusy}
                        title='Download the source volume as MRC, including dust removal and flips'
                        onClick={() => manager.saveVolume()}>
                        ⬇ Save Volume MRC
                    </button>
                </div>

                <label style={labelStyle}>Remove dust — min size: {dustMinVoxels} voxels (uses threshold above)</label>
                <div style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
                    <input type='range' min={1} max={10000} step={1} style={{ flex: 1 }}
                        value={dustMinVoxels}
                        onChange={e => setDustMinVoxels(Number(e.target.value))} />
                    <button style={{ ...btnStyle, background: '#2a1a0a', minWidth: '120px' }}
                        disabled={anyBusy}
                        title='Zero isolated blobs smaller than the minimum in the source volume (cannot be undone)'
                        onClick={() => void run('dust', async () => setDustResult(await manager.removeDust(dustMinVoxels)))}>
                        {busy === 'dust' ? '⏳ Removing…' : '✦ Remove Dust'}
                    </button>
                </div>
                {dustResult !== null && <p style={hintStyle}>{dustResult > 0 ? `Removed ${dustResult.toLocaleString()} dust voxels.` : 'No dust found at this threshold.'}</p>}

                {/* ── Bodies ────────────────────────────────────── */}
                <SectionDivider label={`Bodies (${state.bodies.length})`} />

                {state.bodies.length > 1 && (
                    <p style={{ ...hintStyle, margin: '0 0 4px' }}>Resolved top to bottom: where polygons overlap, the upper body wins.</p>
                )}
                {state.bodies.map((body, index) => (
                    <BodyRow key={body.id} body={body} index={index} state={state} manager={manager} busy={anyBusy} />
                ))}
                <div style={{ display: 'flex', gap: '4px', marginTop: '4px' }}>
                    <button style={{ ...btnStyle, flex: 1, background: '#1a3a6a' }} disabled={anyBusy}
                        onClick={() => manager.addBody()}>
                        + Add body
                    </button>
                    <button style={{ ...btnStyle, flex: 1, background: hasRemainder ? '#222' : '#3a2608' }} disabled={anyBusy || hasRemainder}
                        title={hasRemainder ? 'A remainder body already exists' : 'Add a body that takes every voxel above the threshold not claimed by other bodies'}
                        onClick={() => void run('remainder', () => manager.addRemainderBody())}>
                        ✦ Add remainder body
                    </button>
                    <button style={{ ...btnStyle, background: '#2a2a3a' }} disabled={anyBusy || state.undoDepth === 0}
                        title='Undo the last body change (Ctrl/Cmd+Z)'
                        onClick={() => void run('undo', () => manager.undo())}>
                        ↶
                    </button>
                </div>

                {activeBody && (
                    <BodySettings body={activeBody} state={state} manager={manager} busy={anyBusy}
                        previewViewId={previewViewId} setPreviewViewId={setPreviewViewId} />
                )}

                {/* ── Mask ──────────────────────────────────────── */}
                <SectionDivider label='Mask' />

                <div style={{ display: 'flex', gap: '16px' }}>
                    <div style={{ flex: 1 }}>
                        <label style={labelStyle}>Extend: {state.defaults.extend} vx</label>
                        <input type='range' min='0' max='30' step='1' style={{ width: '100%' }}
                            value={state.defaults.extend}
                            onChange={e => manager.setDefaults({ ...state.defaults, extend: Number(e.target.value) })} />
                    </div>
                    <div style={{ flex: 1 }}>
                        <label style={labelStyle}>Soft edge: {state.defaults.softEdge} vx</label>
                        <input type='range' min='0' max='30' step='1' style={{ width: '100%' }}
                            value={state.defaults.softEdge}
                            onChange={e => manager.setDefaults({ ...state.defaults, softEdge: Number(e.target.value) })} />
                    </div>
                </div>
                <label style={labelStyle}>Preview mask</label>
                <div style={{ display: 'flex', gap: '4px' }}>
                    {(['none', 'active', 'all'] as const).map(mode => (
                        <button key={mode}
                            style={{ ...btnStyle, flex: 1, fontSize: '11px',
                                background: state.preview === mode ? '#334' : '#282828',
                                outline: state.preview === mode ? '1px solid #667' : 'none' }}
                            onClick={() => manager.setPreview(mode)}>
                            {mode === 'none' ? 'None' : mode === 'active' ? 'Active body' : 'All bodies'}
                        </button>
                    ))}
                </div>
                {state.message && <p style={hintStyle}>{state.message}</p>}

                {/* ── Export ────────────────────────────────────── */}
                <SectionDivider label='Export' />

                <div style={{ display: 'flex', gap: '12px', fontSize: '11px', marginBottom: '6px' }}>
                    <label style={{ cursor: 'pointer' }}>
                        <input type='checkbox' checked={exportOptions.bundleZip}
                            onChange={e => setExportOptions({ ...exportOptions, bundleZip: e.target.checked })} /> Bundle as zip
                    </label>
                    <label style={{ cursor: 'pointer', color: exportOptions.bundleZip ? '#ccc' : '#666' }}>
                        <input type='checkbox' checked={exportOptions.compress} disabled={!exportOptions.bundleZip}
                            onChange={e => setExportOptions({ ...exportOptions, compress: e.target.checked })} /> Compress (slow)
                    </label>
                </div>
                <button style={{ ...btnStyle, width: '100%', padding: '8px', fontSize: '13px', background: exportable > 0 ? '#1a3a22' : '#222' }}
                    disabled={anyBusy || exportable === 0}
                    title='One MRC mask per body on the source grid, largest body first'
                    onClick={() => void run('export', () => manager.exportMasks(exportOptions))}>
                    {busy === 'export' ? '⏳ Exporting…' : `⬇ Export ${exportable} mask${exportable === 1 ? '' : 's'}`}
                </button>
                {exportable === 1 && <p style={hintStyle}>Multi-body refinement needs at least two bodies.</p>}

                <StatusText state={state} />
            </>)}

            {/* Drawing / preview canvas portal */}
            {(state.isDrawing || previewView) && viewportEl && ReactDOM.createPortal(
                <DrawingCanvas
                    color={activeColor}
                    previewMask={state.isDrawing ? undefined : previewView}
                    onPolygonComplete={state.isDrawing ? handlePolygonComplete : undefined}
                />,
                viewportEl
            )}
        </div>
    );
}

function BodyRow({ body, index, state, manager, busy }: { body: BodyInfo, index: number, state: VolumeSegmentorState, manager: VolumeSegmentorManager, busy: boolean }) {
    const isActive = body.id === state.activeBodyId;
    const previewed = state.preview === 'all' || (state.preview === 'active' && isActive);
    const total = Math.max(1, state.stats.candidates);
    const detail = body.remainder ? 'rest' : `${body.views.length} view${body.views.length === 1 ? '' : 's'}`;
    return (
        <div style={{ ...rowStyle, cursor: 'pointer', outline: isActive ? `1px solid ${Color.toHexStyle(body.color)}` : 'none' }}
            onClick={() => manager.setActiveBody(body.id)}
            title={isActive ? 'Active body' : 'Click to make this the active body'}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '6px', minWidth: 0 }}>
                <input type='color' value={Color.toHexStyle(body.color)} title='Body color'
                    style={{ width: '18px', height: '18px', padding: 0, border: 'none', background: 'none', cursor: 'pointer', flexShrink: 0 }}
                    onClick={e => e.stopPropagation()}
                    onChange={e => void manager.setBodyColor(body.id, Color.fromHexStyle(e.target.value))} />
                <div style={{ display: 'flex', flexDirection: 'column', minWidth: 0 }}>
                    <span style={{ fontSize: '12px', fontWeight: isActive ? 600 : 400, color: isActive ? '#fff' : '#ccc', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
                        {body.name}
                    </span>
                    <span style={{ fontSize: '10px', color: '#888', whiteSpace: 'nowrap' }}>
                        {detail} · {body.voxelCount.toLocaleString()} vx · {(100 * body.voxelCount / total).toFixed(1)}%
                    </span>
                </div>
            </div>
            <div style={{ display: 'flex', gap: '2px', flexShrink: 0 }}>
                <button style={{ ...smallBtnStyle, background: '#2a2a3a' }} disabled={busy || index === 0} title='Move up (higher priority)'
                    onClick={e => { e.stopPropagation(); void manager.moveBody(body.id, -1); }}>▲</button>
                <button style={{ ...smallBtnStyle, background: '#2a2a3a' }} disabled={busy || index === state.bodies.length - 1} title='Move down (lower priority)'
                    onClick={e => { e.stopPropagation(); void manager.moveBody(body.id, 1); }}>▼</button>
                <button style={{ ...smallBtnStyle, background: previewed ? '#1a4f7a' : '#2a2a3a' }}
                    title='Preview this body mask'
                    onClick={e => {
                        e.stopPropagation();
                        if (previewed && state.preview === 'active') void manager.setPreview('none');
                        else { manager.setActiveBody(body.id); void manager.setPreview('active'); }
                    }}>👁</button>
                <button style={{ ...smallBtnStyle, background: '#4a1a1a' }} disabled={busy}
                    title='Remove body'
                    onClick={e => { e.stopPropagation(); void manager.removeBody(body.id); }}>✕</button>
            </div>
        </div>
    );
}

function BodySettings({ body, state, manager, busy, previewViewId, setPreviewViewId }: {
    body: BodyInfo, state: VolumeSegmentorState, manager: VolumeSegmentorManager, busy: boolean,
    previewViewId: string | null, setPreviewViewId: (id: string | null) => void
}) {
    const custom = body.extend !== undefined || body.softEdge !== undefined;
    const extend = body.extend ?? state.defaults.extend;
    const softEdge = body.softEdge ?? state.defaults.softEdge;
    const color = Color.toHexStyle(body.color);
    return (
        <div style={{ border: `1px solid ${color}`, borderRadius: '4px', padding: '6px 8px', background: '#1a1a1a', marginTop: '6px' }}>
            <label style={labelStyle}>
                Name:&nbsp;
                <input type='text' style={{ ...numInputStyle, width: '140px' }} value={body.name}
                    onChange={e => manager.renameBody(body.id, e.target.value)} />
            </label>

            {body.remainder ? (
                <p style={hintStyle}>Takes every voxel above the threshold that no other body claims.</p>
            ) : (<>
                <label style={{ ...labelStyle, marginTop: '8px' }}>Views ({body.views.length}) — voxels inside all views belong to this body</label>
                {body.views.map((v, idx) => (
                    <div key={v.id}
                        style={{ ...rowStyle, cursor: 'pointer', outline: previewViewId === v.id ? `1px solid ${color}` : 'none' }}
                        onClick={() => { setPreviewViewId(previewViewId === v.id ? null : v.id); manager.flyTo(v); }}
                        title='Click to show this view'>
                        <span style={{ fontSize: '11px' }}>#{idx + 1} {v.label}{v.inverted ? ' ⊘' : ''}</span>
                        <div style={{ display: 'flex', gap: '2px' }}>
                            <button style={{ ...smallBtnStyle, background: v.inverted ? '#664400' : '#2a2a3a' }} disabled={busy}
                                title={v.inverted ? 'Inverted — click to select inside again' : 'Invert: select outside this polygon'}
                                onClick={e => { e.stopPropagation(); void manager.invertView(body.id, v.id); }}>⊘</button>
                            <button style={{ ...smallBtnStyle, background: '#4a1a1a' }} disabled={busy}
                                title='Remove view'
                                onClick={e => { e.stopPropagation(); void manager.removeView(body.id, v.id); if (previewViewId === v.id) setPreviewViewId(null); }}>✕</button>
                        </div>
                    </div>
                ))}
                <button
                    style={{ ...btnStyle, background: state.isDrawing ? '#b85000' : '#333', width: '100%', marginTop: '4px' }}
                    disabled={busy}
                    title='Draw a polygon over the viewport; rotate and add more views to intersect them'
                    onClick={() => { manager.setDrawing(!state.isDrawing); setPreviewViewId(null); }}>
                    {state.isDrawing ? '✏ Drawing active — click to stop' : `✏ Draw polygon for ${body.name}`}
                </button>
                {state.isDrawing && (
                    <p style={{ color: '#888', fontSize: '11px', margin: '4px 0 0' }}>
                        Click to add vertices · Double-click or near start to close · Esc cancels
                    </p>
                )}
            </>)}

            <label style={{ ...labelStyle, cursor: 'pointer', marginTop: '8px' }}>
                <input type='checkbox' checked={custom}
                    onChange={e => void manager.setBodyMaskParams(body.id, e.target.checked ? { extend, softEdge } : {})} />
                &nbsp;Custom extend / soft edge for this body
            </label>
            {custom && (
                <div style={{ display: 'flex', gap: '16px' }}>
                    <div style={{ flex: 1 }}>
                        <label style={labelStyle}>Extend: {extend} vx</label>
                        <input type='range' min='0' max='30' step='1' style={{ width: '100%' }} value={extend}
                            onChange={e => void manager.setBodyMaskParams(body.id, { extend: Number(e.target.value), softEdge })} />
                    </div>
                    <div style={{ flex: 1 }}>
                        <label style={labelStyle}>Soft edge: {softEdge} vx</label>
                        <input type='range' min='0' max='30' step='1' style={{ width: '100%' }} value={softEdge}
                            onChange={e => void manager.setBodyMaskParams(body.id, { extend, softEdge: Number(e.target.value) })} />
                    </div>
                </div>
            )}
        </div>
    );
}

function StatusText({ state }: { state: VolumeSegmentorState }) {
    const { candidates, unassigned } = state.stats;
    const pct = candidates > 0 ? (100 * unassigned / candidates).toFixed(1) : '0';
    const empty = state.bodies.filter(b => b.voxelCount === 0);
    return (
        <p style={{ ...hintStyle, marginTop: '10px' }}>
            {state.busy && <>⏳ Updating bodies… </>}
            {unassigned.toLocaleString()} of {candidates.toLocaleString()} voxels above threshold unassigned ({pct}%).
            {empty.length > 0 && <> Empty: {empty.map(b => b.name).join(', ')}.</>}
        </p>
    );
}
