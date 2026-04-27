/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import * as React from 'react';
import { useEffect, useRef, useState, useCallback } from 'react';
import * as ReactDOM from 'react-dom';
import { PluginContext } from '../../../mol-plugin/context';
import { Volume } from '../../../mol-model/volume';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { OpenFiles } from '../../../mol-plugin-state/actions/file';
import { Asset } from '../../../mol-util/assets';
import { VolumeMaskController } from '../behavior';
import { DrawingCanvas } from './drawing-canvas';
import type { MaskCreatorState, MaskSource, ViewMask, Point2D } from '../types';

interface Props {
    plugin: PluginContext;
    controller: VolumeMaskController;
}

const DEFERRED_THRESHOLD_PREVIEW_DIMENSION = 280;

export function MaskCreatorPanel({ plugin, controller }: Props) {
    const [state, setState] = useState<MaskCreatorState>(controller.current);
    const [previewMaskId, setPreviewMaskId] = useState<string | null>(null);
    const [busy, setBusy] = useState<string | null>(null);
    const [, forceUpdate] = useState(0);
    const [dustMinVoxels, setDustMinVoxels] = useState(100);
    const volumeInputRef = useRef<HTMLInputElement>(null);
    const structureInputRef = useRef<HTMLInputElement>(null);

    useEffect(() => {
        const sub = controller.state.subscribe(setState);
        return () => sub.unsubscribe();
    }, [controller]);

    // Re-render when plugin state changes (e.g. visibility toggled via Molstar UI)
    useEffect(() => {
        const sub = plugin.state.data.events.changed.subscribe(() => forceUpdate(n => n + 1));
        return () => sub.unsubscribe();
    }, [plugin]);

    // Auto-select first volume on load
    useEffect(() => {
        const sub = plugin.managers.volume.hierarchy.behaviors.selection.subscribe(() => {
            const volumes = plugin.managers.volume.hierarchy.current.volumes;
            if (volumes.length > 0 && !controller.current.targetVolumeRef) {
                controller.setTargetVolume(volumes[0].cell.transform.ref);
            }
        });
        return () => sub.unsubscribe();
    }, [plugin, controller]);

    // Auto-select first structure on load
    useEffect(() => {
        const sub = plugin.managers.structure.hierarchy.behaviors.selection.subscribe(() => {
            const structures = plugin.managers.structure.hierarchy.current.structures;
            if (structures.length > 0 && !controller.current.targetStructureRef) {
                controller.setTargetStructure(structures[0].cell.transform.ref);
            }
        });
        return () => sub.unsubscribe();
    }, [plugin, controller]);

    const openFile = useCallback(async (
        e: React.ChangeEvent<HTMLInputElement>,
        onLoaded?: () => void
    ) => {
        const file = e.target.files?.[0];
        if (!file) return;
        e.target.value = '';
        try {
            await plugin.runTask(
                plugin.state.data.applyAction(OpenFiles, {
                    files: [Asset.File(file)],
                    format: { name: 'auto', params: {} },
                    visuals: true,
                })
            );
            onLoaded?.();
        } catch (err) {
            alert(String(err));
        }
    }, [plugin]);

    const handleOpenVolume = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
        await plugin.clear();
        controller.setTargetVolume(undefined);
        controller.setTargetStructure(undefined);
        await openFile(e);
    }, [plugin, controller, openFile]);

    const handleOpenStructure = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
        return openFile(e);
    }, [openFile]);

    // Labels
    const volumes = plugin.managers.volume.hierarchy.current.volumes;
    const structures = plugin.managers.structure.hierarchy.current.structures;

    const currentVolumeLabel = state.targetVolumeRef
        ? (volumes.find(v => v.cell.transform.ref === state.targetVolumeRef)?.cell.obj?.label ?? 'Volume')
        : null;
    const currentStructureLabel = state.targetStructureRef
        ? (structures.find(s => s.cell.transform.ref === state.targetStructureRef)?.cell.obj?.label ?? 'Structure')
        : null;

    // Visibility state (read from state tree)
    const volumeHidden = state.targetVolumeRef
        ? !!plugin.state.data.cells.get(state.targetVolumeRef)?.state.isHidden
        : false;
    const structureHidden = state.targetStructureRef
        ? !!plugin.state.data.cells.get(state.targetStructureRef)?.state.isHidden
        : false;

    const handlePolygonComplete = useCallback((polygon: Point2D[], canvasW: number, canvasH: number) => {
        const cam = plugin.canvas3d?.camera;
        if (!cam) return;
        const vp = cam.viewport;
        const mask: ViewMask = {
            id: `mask-${Date.now()}`,
            label: `View ${state.viewMasks.length + 1}`,
            polygon, canvasWidth: canvasW, canvasHeight: canvasH,
            viewportWidth: vp.width, viewportHeight: vp.height,
            cameraSnapshot: cam.getSnapshot(),
        };
        controller.addViewMask(mask);
        controller.setDrawing(false);
    }, [plugin, controller, state.viewMasks.length]);

    const volumeObj = state.targetVolumeRef
        ? plugin.state.data.select(state.targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined
        : undefined;

    const volStats = volumeObj?.data.grid.stats;
    const volumeDimensions = volumeObj?.data.grid.cells.space.dimensions as [number, number, number] | undefined;
    const deferThresholdSliderPreview = !!volumeDimensions && Math.max(...volumeDimensions) > DEFERRED_THRESHOLD_PREVIEW_DIMENSION;

    const thresholdAbs = volStats
        ? Volume.IsoValue.toAbsolute(state.threshold, volStats).absoluteValue
        : null;

    const thresholdSigma = (volStats && thresholdAbs !== null && volStats.sigma > 0)
        ? (thresholdAbs - volStats.mean) / volStats.sigma
        : null;

    const viewportEl = (
        (plugin.layout.root as HTMLElement | undefined)?.querySelector('.msp-viewport') ??
        document.querySelector('.msp-viewport')
    ) as HTMLElement | null;

    const needsStructure = state.maskSource === 'structure';
    const canCreate = !!(state.targetVolumeRef && (!needsStructure || state.targetStructureRef));

    const availableChains = state.availableChainIds;

    const handleThresholdInputChange = useCallback((value: number) => {
        controller.previewThreshold(Volume.IsoValue.absolute(value));
    }, [controller]);

    const handleThresholdSliderChange = useCallback((value: number) => {
        const isoValue = Volume.IsoValue.absolute(value);
        if (deferThresholdSliderPreview) controller.setThreshold(isoValue);
        else controller.previewThreshold(isoValue);
    }, [controller, deferThresholdSliderPreview]);

    const commitThresholdSliderPreview = useCallback(() => {
        if (deferThresholdSliderPreview) controller.commitThresholdPreview();
    }, [controller, deferThresholdSliderPreview]);

    return (
        <div style={{ padding: '8px', fontFamily: 'sans-serif', fontSize: '13px', color: '#ccc' }}>
            <h3 style={{ margin: '0 0 8px', color: '#fff' }}>Volume Mask Creator</h3>

            {/* ── Data ─────────────────────────────────────────── */}
            <input ref={volumeInputRef} type='file'
                accept='.mrc,.ccp4,.map,.dsn6,.brix,.dx,.dxbin,.cub,.cube,.gz'
                style={{ display: 'none' }} onChange={handleOpenVolume} />
            <input ref={structureInputRef} type='file'
                accept='.pdb,.ent,.cif,.mmcif,.mcif,.bcif,.mol2,.sdf,.mol,.pdbqt,.gro,.xyz,.cub,.cube'
                style={{ display: 'none' }} onChange={handleOpenStructure} />

            <div style={{ display: 'flex', gap: '4px', marginBottom: '2px' }}>
                <button style={{ ...btnStyle, flex: 1, background: '#1e5230' }} onClick={() => volumeInputRef.current?.click()}>
                    Open Volume
                </button>
                {currentVolumeLabel && (<>
                    <button style={{ ...btnStyle, background: volumeHidden ? '#444' : '#225' }}
                        title={volumeHidden ? 'Show volume' : 'Hide volume'}
                        onClick={() => controller.toggleVolumeVisibility()}>
                        {volumeHidden ? '○' : '●'}
                    </button>
                    <button style={{ ...btnStyle, flex: 1, background: '#1e3a6a' }} onClick={() => structureInputRef.current?.click()}>
                        Open Structure
                    </button>
                    {currentStructureLabel && (
                        <button style={{ ...btnStyle, background: structureHidden ? '#444' : '#225' }}
                            title={structureHidden ? 'Show structure' : 'Hide structure'}
                            onClick={() => controller.toggleStructureVisibility()}>
                            {structureHidden ? '○' : '●'}
                        </button>
                    )}
                </>)}
            </div>
            {currentVolumeLabel && (
                <div style={{ display: 'flex', gap: '8px', fontSize: '11px', marginBottom: '8px' }}>
                    <span style={{ flex: 1, color: '#7c7', wordBreak: 'break-all' }}>
                        {currentVolumeLabel}
                    </span>
                    <span style={{ flex: 1, color: currentStructureLabel ? '#88f' : '#555', wordBreak: 'break-all' }}>
                        {currentStructureLabel ?? 'No structure'}
                    </span>
                </div>
            )}

            {/* ── Source volume ─────────────────────────────────── */}
            {state.targetVolumeRef && (<>
                <SectionDivider label='Source volume' />

                {volStats && (<>
                    <label style={labelStyle}>
                        Threshold:&nbsp;
                        <input type='number' step={volStats.sigma / 20 || 0.001} style={numInputStyle}
                            value={thresholdAbs !== null ? parseFloat(thresholdAbs.toPrecision(5)) : ''}
                            onChange={e => {
                                const v = parseFloat(e.target.value);
                                if (!isNaN(v)) handleThresholdInputChange(v);
                            }} />
                        {thresholdSigma !== null && <span style={{ color: '#888' }}>&nbsp;({thresholdSigma.toFixed(1)}σ)</span>}
                    </label>
                    <input type='range' min={volStats.min} max={volStats.max} step={volStats.sigma / 20 || 0.001} style={{ width: '100%' }}
                        value={thresholdAbs ?? volStats.mean}
                        onChange={e => handleThresholdSliderChange(Number(e.target.value))}
                        onPointerUp={commitThresholdSliderPreview}
                        onKeyUp={commitThresholdSliderPreview} />
                    {deferThresholdSliderPreview && (
                        <p style={{ color: '#666', fontSize: '11px', margin: '4px 0 0', fontStyle: 'italic' }}>
                            Large volume: the slider updates the numeric threshold immediately and refreshes the surface on release.
                        </p>
                    )}
                </>)}

                {currentStructureLabel && (<>
                    <label style={labelStyle}>Opacity: {(state.volumeOpacity * 100).toFixed(0)}%</label>
                    <input type='range' min={0} max={1} step={0.01} style={{ width: '100%' }}
                        value={state.volumeOpacity}
                        onChange={e => controller.setVolumeOpacity(Number(e.target.value))} />
                </>)}

                <div style={{ display: 'flex', gap: '4px', marginTop: '8px', marginBottom: '6px' }}>
                    <button style={{ ...btnStyle, flex: 1, background: '#2a1a40' }}
                        disabled={!!busy}
                        onClick={async () => {
                            setBusy('flip');
                            try { await controller.flipHandedness(); } catch (e) { alert(String(e)); } finally { setBusy(null); }
                        }}>
                        {busy === 'flip' ? '⏳ Flipping…' : '↔ Flip Handedness'}
                    </button>
                    <button style={{ ...btnStyle, flex: 1, background: '#152a1a' }}
                        disabled={!!busy}
                        onClick={() => { try { controller.saveVolume(); } catch (e) { alert(String(e)); } }}>
                        ⬇ Save Volume MRC
                    </button>
                </div>

                <label style={labelStyle}>Remove dust — min size: {dustMinVoxels} voxels (uses threshold above)</label>
                <div style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
                    <input type='range' min={1} max={10000} step={1} style={{ flex: 1 }}
                        value={dustMinVoxels}
                        onChange={e => setDustMinVoxels(Number(e.target.value))} />
                    <button style={{ ...btnStyle, background: '#2a1a0a', minWidth: '120px' }}
                        disabled={!!busy}
                        onClick={async () => {
                            setBusy('dust');
                            try { await controller.removeDust(dustMinVoxels); } catch (e) { alert(String(e)); } finally { setBusy(null); }
                        }}>
                        {busy === 'dust' ? '⏳ Removing…' : '✦ Remove Dust'}
                    </button>
                </div>
            </>)}

            {/* ── Mask creation ─────────────────────────────────── */}
            {state.targetVolumeRef && (<>
            <SectionDivider label='Mask creation' />

            {currentStructureLabel && (<>
                <label style={labelStyle}>Mask source</label>
                <div style={{ display: 'flex', gap: '4px', marginBottom: '6px' }}>
                    {(['volume', 'structure'] as MaskSource[]).map(src => (
                        <button key={src}
                            style={{ ...btnStyle, flex: 1, fontSize: '11px',
                                background: state.maskSource === src ? '#334' : '#282828',
                                outline: state.maskSource === src ? '1px solid #667' : 'none' }}
                            onClick={() => controller.setMaskSource(src)}>
                            {src === 'volume' ? 'Volume threshold' : 'Structure proximity'}
                        </button>
                    ))}
                </div>
            </>)}

            {state.maskSource === 'structure' && currentStructureLabel ? (<>
                {availableChains.length > 0 && (<>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline' }}>
                        <label style={labelStyle}>
                            Chains ({state.selectedChainIds.length}/{availableChains.length})
                        </label>
                        <div style={{ display: 'flex', gap: '6px' }}>
                            <button style={{ ...smallBtnStyle, background: '#2a2a2a', fontSize: '10px' }}
                                onClick={() => controller.setSelectedChainIds(availableChains)}>All</button>
                            <button style={{ ...smallBtnStyle, background: '#2a2a2a', fontSize: '10px' }}
                                onClick={() => controller.setSelectedChainIds([])}>None</button>
                        </div>
                    </div>
                    <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px', marginBottom: '6px' }}>
                        {availableChains.map(c => {
                            const selected = state.selectedChainIds.includes(c);
                            return (
                                <button key={c}
                                    style={{ ...smallBtnStyle, fontFamily: 'monospace', minWidth: '28px',
                                        background: selected ? '#1e3a6a' : '#282828',
                                        outline: selected ? '1px solid #4a7ac0' : 'none',
                                        color: selected ? '#fff' : '#888' }}
                                    onClick={() => controller.toggleChainId(c)}>
                                    {c}
                                </button>
                            );
                        })}
                    </div>
                    <p style={{ color: '#666', fontSize: '11px', margin: '0 0 6px', fontStyle: 'italic' }}>
                        Selected chains are highlighted in the structure view.
                    </p>
                </>)}
                <label style={labelStyle}>Protein radius: {state.proteinRadius.toFixed(1)} Å</label>
                <input type='range' min='0.5' max='20' step='0.5' style={{ width: '100%' }}
                    value={state.proteinRadius}
                    onChange={e => controller.setProteinRadius(Number(e.target.value))} />
            </>) : (
                <p style={{ color: '#666', fontSize: '11px', margin: '0 0 6px', fontStyle: 'italic' }}>
                    Uses the threshold from the source volume above.
                </p>
            )}

            <div style={{ display: 'flex', gap: '16px' }}>
                <div style={{ flex: 1 }}>
                    <label style={labelStyle}>Dilation: {state.dilation} vx</label>
                    <input type='range' min='0' max='20' step='1' style={{ width: '100%' }}
                        value={state.dilation}
                        onChange={e => controller.setDilation(Number(e.target.value))} />
                </div>
                <div style={{ flex: 1 }}>
                    <label style={labelStyle}>Soft edge: {state.softEdge} vx</label>
                    <input type='range' min='0' max='20' step='1' style={{ width: '100%' }}
                        value={state.softEdge}
                        onChange={e => controller.setSoftEdge(Number(e.target.value))} />
                </div>
            </div>

            <button
                style={{ ...btnStyle, background: state.isDrawing ? '#b85000' : '#333', width: '100%', marginTop: '8px' }}
                onClick={() => { controller.setDrawing(!state.isDrawing); setPreviewMaskId(null); }}
                disabled={!state.targetVolumeRef}>
                {state.isDrawing ? '✏ Drawing active — click to stop' : '✏ Draw polygon regions'}
            </button>
            {state.isDrawing && (
                <p style={{ color: '#888', fontSize: '11px', margin: '4px 0 0' }}>
                    Click to add vertices · Double-click or near start to close
                </p>
            )}

            {state.viewMasks.length > 0 && (
                <div style={{ marginTop: '6px' }}>
                    <label style={{ ...labelStyle, marginTop: '4px' }}>Polygons ({state.viewMasks.length})</label>
                    {state.viewMasks.map((m, idx) => (
                        <div key={m.id}
                            style={{ ...rowStyle, cursor: 'pointer', outline: previewMaskId === m.id ? '1px solid #FF6B00' : 'none' }}
                            onClick={() => { setPreviewMaskId(id => id === m.id ? null : m.id); plugin.canvas3d?.camera.setState(m.cameraSnapshot, 400); }}
                            title='Click to preview this view'>
                            <span style={{ fontSize: '11px' }}>#{idx + 1} {m.label}{m.inverted ? ' ⊘' : ''}</span>
                            <div style={{ display: 'flex', gap: '2px' }}>
                                <button style={{ ...smallBtnStyle, background: m.inverted ? '#664400' : '#2a2a3a' }}
                                    title={m.inverted ? 'Inverted — click to revert' : 'Invert selection'}
                                    onClick={e => { e.stopPropagation(); controller.invertViewMask(m.id); }}>⊘</button>
                                <button style={{ ...smallBtnStyle, background: '#4a1a1a' }}
                                    onClick={e => { e.stopPropagation(); controller.removeViewMask(m.id); if (previewMaskId === m.id) setPreviewMaskId(null); }}>✕</button>
                            </div>
                        </div>
                    ))}
                </div>
            )}

            <button style={{ ...btnStyle, width: '100%', marginTop: '10px', padding: '8px', background: canCreate ? '#1a4f7a' : '#222', fontSize: '13px' }}
                disabled={!canCreate || !!busy}
                onClick={async () => {
                    setBusy('create');
                    try { await controller.createMask(); } catch (e) { alert(String(e)); } finally { setBusy(null); }
                }}>
                {busy === 'create' ? '⏳ Creating mask…' : '⬡ Create Mask'}
            </button>
            </>)}

            {/* ── Mask result ───────────────────────────────────── */}
            {state.maskVolumeRef && (<>
                <SectionDivider label='Mask' />
                <div style={{ border: '1px solid #333', borderRadius: '4px', padding: '8px', background: '#1a1a1a' }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '8px' }}>
                        <span style={{ color: '#e08020', fontWeight: 600, fontSize: '12px' }}>
                            ● Mask{state.maskInverted ? ' (inverted)' : ''}
                        </span>
                        <button style={{ ...smallBtnStyle, background: '#3a1010' }}
                            title='Delete mask'
                            disabled={!!busy}
                            onClick={() => controller.clearMask()}>
                            ✕
                        </button>
                    </div>
                    <div style={{ display: 'flex', gap: '4px', marginBottom: '6px' }}>
                        <button style={{ ...btnStyle, flex: 1, background: state.maskInverted ? '#5a3300' : '#252525', outline: state.maskInverted ? '1px solid #e08020' : 'none' }}
                            disabled={!!busy}
                            onClick={async () => {
                                setBusy('invert');
                                try { await controller.invertMask(); } catch (e) { alert(String(e)); } finally { setBusy(null); }
                            }}>
                            {busy === 'invert' ? '⏳' : '⊘ Invert'}
                        </button>
                        <button style={{ ...btnStyle, flex: 1, background: '#1a3a22' }}
                            disabled={!state.maskData}
                            onClick={() => { try { controller.saveMask(); } catch (e) { alert(String(e)); } }}>
                            ⬇ Save Mask
                        </button>
                    </div>
                    <button style={{ ...btnStyle, width: '100%', background: '#3a2608', padding: '7px' }}
                        disabled={!!busy}
                        onClick={async () => {
                            setBusy('apply');
                            try { await controller.applyMaskToVolume(); } catch (e) { alert(String(e)); } finally { setBusy(null); }
                        }}>
                        {busy === 'apply' ? '⏳ Applying mask to volume…' : '✱ Apply Mask to Volume'}
                    </button>
                </div>
            </>)}

            {/* Drawing / preview canvas portal */}
            {(state.isDrawing || previewMaskId) && viewportEl && ReactDOM.createPortal(
                <DrawingCanvas
                    previewMask={previewMaskId ? state.viewMasks.find(m => m.id === previewMaskId) : undefined}
                    onPolygonComplete={state.isDrawing ? handlePolygonComplete : undefined}
                />,
                viewportEl
            )}
        </div>
    );
}

function SectionDivider({ label }: { label: string }) {
    return (
        <div style={{ display: 'flex', alignItems: 'center', gap: '8px', margin: '14px 0 6px' }}>
            <span style={{ color: '#999', fontSize: '12px', fontWeight: 600, textTransform: 'uppercase', letterSpacing: '0.08em', whiteSpace: 'nowrap' }}>{label}</span>
            <div style={{ flex: 1, height: '1px', background: '#333' }} />
        </div>
    );
}

const labelStyle: React.CSSProperties = { display: 'block', color: '#888', marginTop: '4px', marginBottom: '2px', fontSize: '11px' };
const numInputStyle: React.CSSProperties = { width: '60px', background: '#252525', border: '1px solid #444', color: '#ddd', borderRadius: '2px', padding: '1px 4px', fontSize: '12px' };
const btnStyle: React.CSSProperties = { border: 'none', color: '#ddd', padding: '6px 8px', cursor: 'pointer', borderRadius: '3px', fontSize: '12px' };
const smallBtnStyle: React.CSSProperties = { ...btnStyle, padding: '2px 7px', fontSize: '11px' };
const rowStyle: React.CSSProperties = { display: 'flex', justifyContent: 'space-between', alignItems: 'center', background: '#1e1e1e', padding: '3px 6px', marginBottom: '2px', borderRadius: '2px' };
