/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { BehaviorSubject } from 'rxjs';
import { OrderedSet } from '../../mol-data/int';
import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { Volume } from '../../mol-model/volume';
import { Structure, StructureElement, StructureProperties, Unit } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';
import { MarkerAction } from '../../mol-util/marker-action';
import { Vec3 } from '../../mol-math/linear-algebra';
import { StateTransform } from '../../mol-state';
import type { UnitIndex } from '../../mol-model/structure/structure/element/util';
import { MaskVolumeFromSource } from './transformers';
import { downloadMrc } from './internal/mrc-export';
import { flipVolumeX, removeDust } from './internal/volume-ops';
import type { MaskCreatorState, MaskSource, ViewMask } from './types';

export const MASK_OVERLAY_COLOR = Color(0xFF6B00);
const THRESHOLD_PREVIEW_DEBOUNCE_MS = 75;

const defaultState: MaskCreatorState = {
    isDrawing: false,
    targetVolumeRef: undefined,
    targetStructureRef: undefined,
    availableChainIds: [],
    viewMasks: [],
    threshold: Volume.IsoValue.relative(1),
    dilation: 4,
    softEdge: 2,
    maskSource: 'volume',
    proteinRadius: 1,
    maskVolumeRef: undefined,
    maskData: undefined,
    maskInverted: false,
    volumeOpacity: 1,
    selectedChainIds: [],
};

export class VolumeMaskController {
    readonly state = new BehaviorSubject<MaskCreatorState>({ ...defaultState });
    private highlightedChains?: StructureElement.Loci;
    private thresholdPreviewTimer: ReturnType<typeof setTimeout> | undefined;
    private thresholdPreviewInFlight = false;
    private thresholdPreviewSession = 0;
    private pendingThresholdPreview: { value: Volume.IsoValue, targetVolumeRef: StateTransform.Ref, session: number } | undefined;

    constructor(private plugin: PluginContext) {}

    get current() { return this.state.value; }

    private update(patch: Partial<MaskCreatorState>) {
        this.state.next({ ...this.current, ...patch });
    }

    setTargetVolume(ref: StateTransform.Ref | undefined) {
        this.clearPendingThresholdPreview();
        this.thresholdPreviewSession++;

        let threshold: Volume.IsoValue = Volume.IsoValue.relative(3);
        if (ref) {
            const stats = (this.plugin.state.data.select(ref)[0]?.obj as SO.Volume.Data | undefined)?.data.grid.stats;
            if (stats) threshold = Volume.IsoValue.absolute(stats.mean + 3 * stats.sigma);
        }
        this.update({ targetVolumeRef: ref, viewMasks: [], maskVolumeRef: undefined, maskData: undefined, threshold });
    }

    setTargetStructure(ref: StateTransform.Ref | undefined) {
        let availableChainIds: string[] = [];
        let selectedChainIds: string[] = [];
        if (ref) {
            const structObj = this.plugin.state.data.select(ref)[0]?.obj as SO.Molecule.Structure | undefined;
            if (structObj) {
                availableChainIds = getChainIds(structObj.data);
                selectedChainIds = availableChainIds;
            }
        }
        this.update({ targetStructureRef: ref, maskVolumeRef: undefined, maskData: undefined, availableChainIds, selectedChainIds });
        this.syncSelectedChainHighlight();
    }

    setSelectedChainIds(ids: string[]) {
        this.update({ selectedChainIds: ids });
        this.syncSelectedChainHighlight();
    }

    toggleChainId(id: string) {
        const ids = this.current.selectedChainIds;
        const next = ids.includes(id) ? ids.filter(c => c !== id) : [...ids, id];
        this.setSelectedChainIds(next);
    }

    setDrawing(active: boolean) { this.update({ isDrawing: active }); }

    addViewMask(mask: ViewMask) {
        this.update({ viewMasks: [...this.current.viewMasks, mask] });
    }

    removeViewMask(id: string) {
        this.update({ viewMasks: this.current.viewMasks.filter(m => m.id !== id) });
    }

    invertViewMask(id: string) {
        this.update({ viewMasks: this.current.viewMasks.map(m => m.id === id ? { ...m, inverted: !m.inverted } : m) });
    }

    setDilation(value: number) { this.update({ dilation: value }); }
    setSoftEdge(value: number) { this.update({ softEdge: value }); }
    setMaskSource(value: MaskSource) {
        this.update({ maskSource: value });
        this.syncSelectedChainHighlight();
    }
    setProteinRadius(value: number) { this.update({ proteinRadius: value }); }

    setThreshold(value: Volume.IsoValue) {
        this.update({ threshold: value });
    }

    /** Update threshold and schedule a coalesced live preview on the source isosurface. */
    previewThreshold(value: Volume.IsoValue) {
        this.setThreshold(value);

        const { targetVolumeRef } = this.current;
        if (!targetVolumeRef) return;

        this.pendingThresholdPreview = { value, targetVolumeRef, session: this.thresholdPreviewSession };
        this.clearThresholdPreviewTimer();
        this.thresholdPreviewTimer = setTimeout(() => {
            this.thresholdPreviewTimer = undefined;
            void this.flushThresholdPreview();
        }, THRESHOLD_PREVIEW_DEBOUNCE_MS);
    }

    /** Immediately preview the latest threshold value on the source isosurface. */
    commitThresholdPreview() {
        const { targetVolumeRef, threshold } = this.current;
        if (!targetVolumeRef) return;

        this.pendingThresholdPreview = { value: threshold, targetVolumeRef, session: this.thresholdPreviewSession };
        this.clearThresholdPreviewTimer();
        void this.flushThresholdPreview();
    }

    private clearThresholdPreviewTimer() {
        if (this.thresholdPreviewTimer !== undefined) {
            clearTimeout(this.thresholdPreviewTimer);
            this.thresholdPreviewTimer = undefined;
        }
    }

    private clearPendingThresholdPreview() {
        this.clearThresholdPreviewTimer();
        this.pendingThresholdPreview = undefined;
    }

    private async flushThresholdPreview() {
        if (this.thresholdPreviewInFlight) return;

        this.thresholdPreviewInFlight = true;
        try {
            while (this.pendingThresholdPreview) {
                const request = this.pendingThresholdPreview;
                this.pendingThresholdPreview = undefined;

                if (request.session !== this.thresholdPreviewSession || this.current.targetVolumeRef !== request.targetVolumeRef) continue;
                await this.applyThresholdPreview(request.targetVolumeRef, request.value);
            }
        } finally {
            this.thresholdPreviewInFlight = false;
            if (this.pendingThresholdPreview) void this.flushThresholdPreview();
        }
    }

    private async applyThresholdPreview(targetVolumeRef: StateTransform.Ref, value: Volume.IsoValue) {
        const volItem = this.plugin.managers.volume.hierarchy.current.volumes
            .find(v => v.cell.transform.ref === targetVolumeRef);
        if (!volItem?.cell.obj) return;

        const builder = this.plugin.build();
        let changed = false;
        for (const repr of volItem.representations) {
            const params = repr.cell.transform.params as any;
            if (params?.type?.name !== 'isosurface') continue;
            builder.to(repr.cell).update(
                StateTransforms.Representation.VolumeRepresentation3D,
                old => ({ ...old, type: { name: 'isosurface', params: { ...old.type.params, isoValue: value } } })
            );
            changed = true;
        }
        if (changed) await builder.commit({ canUndo: false });
    }

    async setVolumeOpacity(value: number) {
        this.update({ volumeOpacity: value });

        const { targetVolumeRef, maskVolumeRef } = this.current;
        if (!targetVolumeRef) return;

        const volItem = this.plugin.managers.volume.hierarchy.current.volumes
            .find(v => v.cell.transform.ref === targetVolumeRef);
        if (!volItem) return;

        const directReprs = volItem.representations.filter(r => r.cell.transform.parent !== maskVolumeRef);
        if (directReprs.length === 0) return;

        const builder = this.plugin.build();
        for (const repr of directReprs) {
            builder.to(repr.cell).update(
                StateTransforms.Representation.VolumeRepresentation3D,
                old => ({ ...old, type: { ...old.type, params: { ...old.type.params, alpha: value } } })
            );
        }
        await builder.commit({ canUndo: false });
    }

    toggleVolumeVisibility() {
        const { targetVolumeRef } = this.current;
        if (!targetVolumeRef) return;
        const cell = this.plugin.state.data.cells.get(targetVolumeRef);
        if (!cell) return;
        setSubtreeVisibility(this.plugin.state.data, targetVolumeRef, !cell.state.isHidden);
    }

    toggleStructureVisibility() {
        const { targetStructureRef } = this.current;
        if (!targetStructureRef) return;
        const cell = this.plugin.state.data.cells.get(targetStructureRef);
        if (!cell) return;
        setSubtreeVisibility(this.plugin.state.data, targetStructureRef, !cell.state.isHidden);
    }

    async createMask() {
        const { targetVolumeRef, targetStructureRef, viewMasks, threshold, dilation, softEdge, maskSource, proteinRadius } = this.current;
        if (!targetVolumeRef) throw new Error('No target volume loaded');

        if (maskSource === 'structure' && !targetStructureRef) {
            throw new Error('No structure loaded');
        }

        // Extract atom positions if needed
        let atomPositions: Float32Array | undefined;
        if (targetStructureRef && maskSource === 'structure') {
            const structObj = this.plugin.state.data.select(targetStructureRef)[0]?.obj as SO.Molecule.Structure | undefined;
            if (structObj) atomPositions = extractAtomPositions(structObj.data, this.current.selectedChainIds);
        }

        await this.clearMask();

        const tree = this.plugin.build();
        const maskNode = tree
            .to(targetVolumeRef)
            .apply(MaskVolumeFromSource, { viewMasks, threshold, dilation, softEdge, maskSource, atomPositions, proteinRadius, invertOutput: false });

        const maskSelector = await maskNode.commit();
        const maskRef = maskSelector.ref;
        this.update({ maskVolumeRef: maskRef });

        const maskVolume = maskSelector.data;
        const reprParams = createVolumeRepresentationParams(this.plugin, maskVolume, {
            type: 'isosurface',
            typeParams: { isoValue: Volume.IsoValue.absolute(0.5), alpha: 0.45 },
            color: 'uniform',
            colorParams: { value: MASK_OVERLAY_COLOR },
        });
        await this.plugin.build()
            .to(maskRef)
            .apply(StateTransforms.Representation.VolumeRepresentation3D, reprParams)
            .commit();

        const maskVolObj = this.plugin.state.data.select(maskRef)[0]?.obj as SO.Volume.Data | undefined;
        if (maskVolObj) {
            const rawData = maskVolObj.data.grid.cells.data as unknown as ArrayLike<number>;
            const maskData = new Float32Array(rawData.length);
            for (let i = 0; i < maskData.length; i++) maskData[i] = rawData[i];
            this.update({ maskData, maskInverted: false });
        }
    }

    async clearMask() {
        const { maskVolumeRef } = this.current;
        if (!maskVolumeRef) return;
        const state = this.plugin.state.data;
        if (state.select(maskVolumeRef).length > 0) {
            await state.build().delete(maskVolumeRef).commit();
        }
        this.update({ maskVolumeRef: undefined, maskData: undefined, maskInverted: false });
    }

    async invertMask() {
        const { maskVolumeRef, maskInverted } = this.current;
        if (!maskVolumeRef) return;

        this.update({ maskInverted: !maskInverted });

        await this.plugin.build()
            .to(maskVolumeRef)
            .update(MaskVolumeFromSource, old => ({ ...old, invertOutput: !old.invertOutput }))
            .commit({ canUndo: false });

        // Pull fresh maskData from the recomputed volume (for MRC export)
        const maskVolObj = this.plugin.state.data.select(maskVolumeRef)[0]?.obj as SO.Volume.Data | undefined;
        if (maskVolObj) {
            const rawData = maskVolObj.data.grid.cells.data as unknown as ArrayLike<number>;
            const maskData = new Float32Array(rawData.length);
            for (let i = 0; i < maskData.length; i++) maskData[i] = rawData[i];
            this.update({ maskData });
        }
    }

    async applyMaskToVolume() {
        const { targetVolumeRef, maskData } = this.current;
        if (!targetVolumeRef) throw new Error('No target volume loaded');
        if (!maskData) throw new Error('No mask generated yet');

        const volObj = this.plugin.state.data.select(targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined;
        if (!volObj) throw new Error('Target volume not found');

        const cells = volObj.data.grid.cells.data as unknown as { [i: number]: number, length: number };
        if (cells.length !== maskData.length) {
            throw new Error(`Mask/volume size mismatch (${maskData.length} vs ${cells.length})`);
        }

        for (let i = 0; i < cells.length; i++) cells[i] = cells[i] * maskData[i];

        await this.finalizeVolumeEdit(targetVolumeRef, volObj);
    }

    async flipHandedness() {
        const { targetVolumeRef } = this.current;
        if (!targetVolumeRef) throw new Error('No target volume loaded');

        const volObj = this.plugin.state.data.select(targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined;
        if (!volObj) throw new Error('Target volume not found');

        const grid = volObj.data.grid;
        const [nx, ny, nz] = grid.cells.space.dimensions as [number, number, number];
        const cells = grid.cells.data as unknown as { [i: number]: number, length: number };

        flipVolumeX(cells, nx, ny, nz, grid.cells.space);

        await this.finalizeVolumeEdit(targetVolumeRef, volObj);
    }

    async removeDust(minVoxels: number) {
        const { targetVolumeRef, threshold } = this.current;
        if (!targetVolumeRef) throw new Error('No target volume loaded');

        const volObj = this.plugin.state.data.select(targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined;
        if (!volObj) throw new Error('Target volume not found');

        const grid = volObj.data.grid;
        const [nx, ny, nz] = grid.cells.space.dimensions as [number, number, number];
        const cells = grid.cells.data as unknown as { [i: number]: number, length: number };

        const thresholdAbs = Volume.IsoValue.toAbsolute(threshold, grid.stats).absoluteValue;
        removeDust(cells, nx, ny, nz, grid.cells.space, minVoxels, thresholdAbs);

        await this.finalizeVolumeEdit(targetVolumeRef, volObj);
    }

    /** Recompute stats, invalidate GPU caches, and refresh representations after an in-place edit. */
    private async finalizeVolumeEdit(targetVolumeRef: string, volObj: SO.Volume.Data) {
        const cells = volObj.data.grid.cells.data as unknown as { [i: number]: number, length: number };

        let min = Infinity, max = -Infinity, sum = 0;
        for (let i = 0; i < cells.length; i++) {
            const v = cells[i];
            if (v < min) min = v;
            if (v > max) max = v;
            sum += v;
        }
        const mean = sum / cells.length;
        let sqSum = 0;
        for (let i = 0; i < cells.length; i++) { const d = cells[i] - mean; sqSum += d * d; }
        const sigma = Math.sqrt(sqSum / cells.length);
        const stats = volObj.data.grid.stats as { min: number, max: number, mean: number, sigma: number };
        stats.min = min; stats.max = max; stats.mean = mean; stats.sigma = sigma;

        // Invalidate per-volume GPU texture caches (keyed by Volume identity, not data contents).
        volObj.data.customProperties.dispose();
        for (const k of Object.keys(volObj.data._propertyData)) delete volObj.data._propertyData[k];
        for (const k of Object.keys(volObj.data._localPropertyData)) delete volObj.data._localPropertyData[k];

        // Force representations to re-apply by delete+recreate. Exclude the mask's
        // orange overlay repr (its parent is maskVolumeRef, not the volume root).
        const { maskVolumeRef } = this.current;
        const volItem = this.plugin.managers.volume.hierarchy.current.volumes
            .find(v => v.cell.transform.ref === targetVolumeRef);
        const directReprs = (volItem?.representations ?? []).filter(
            r => r.cell.transform.parent !== maskVolumeRef
        );
        if (directReprs.length > 0) {
            const saved = directReprs.map(r => ({
                transformer: r.cell.transform.transformer,
                params: r.cell.transform.params,
            }));
            const delBuilder = this.plugin.build();
            for (const r of directReprs) delBuilder.delete(r.cell.transform.ref);
            await delBuilder.commit({ canUndo: false });

            const addBuilder = this.plugin.build();
            for (const s of saved) addBuilder.to(targetVolumeRef).apply(s.transformer as any, s.params);
            await addBuilder.commit({ canUndo: false });
        }
    }

    saveMask() {
        const { targetVolumeRef, maskData } = this.current;
        if (!maskData) throw new Error('No mask generated yet');
        const volObj = targetVolumeRef
            ? (this.plugin.state.data.select(targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined)
            : undefined;
        if (!volObj) throw new Error('Source volume not found');
        downloadMrc(volObj.data.grid, maskData, deriveName(volObj.label, 'mask'));
    }

    saveVolume() {
        const { targetVolumeRef } = this.current;
        if (!targetVolumeRef) throw new Error('No target volume loaded');
        const volObj = this.plugin.state.data.select(targetVolumeRef)[0]?.obj as SO.Volume.Data | undefined;
        if (!volObj) throw new Error('Target volume not found');
        const raw = volObj.data.grid.cells.data as unknown as ArrayLike<number>;
        const data = raw instanceof Float32Array ? raw : Float32Array.from(raw as any);
        downloadMrc(volObj.data.grid, data, deriveName(volObj.label, 'crop'));
    }

    dispose() {
        this.clearPendingThresholdPreview();
        this.clearSelectedChainHighlight();
        this.clearMask().catch(() => {});
    }

    private clearSelectedChainHighlight() {
        if (this.highlightedChains && !StructureElement.Loci.isEmpty(this.highlightedChains)) {
            this.plugin.canvas3d?.mark({ loci: this.highlightedChains }, MarkerAction.RemoveHighlight);
        }
        this.highlightedChains = undefined;
    }

    private syncSelectedChainHighlight() {
        this.clearSelectedChainHighlight();

        const { maskSource, targetStructureRef, selectedChainIds } = this.current;
        if (maskSource !== 'structure' || !targetStructureRef || selectedChainIds.length === 0) return;

        const structObj = this.plugin.state.data.select(targetStructureRef)[0]?.obj as SO.Molecule.Structure | undefined;
        if (!structObj) return;

        const loci = getChainLoci(structObj.data, selectedChainIds);
        if (StructureElement.Loci.isEmpty(loci)) return;

        this.plugin.canvas3d?.mark({ loci }, MarkerAction.Highlight);
        this.highlightedChains = loci;
    }
}

function deriveName(label: string | undefined, suffix: string): string {
    const base = (label ?? '').replace(/\.[^.]+$/, '').trim() || 'volume';
    return `${base}_${suffix}.mrc`;
}

/** All unique auth chain IDs in the structure, sorted. */
export function getChainIds(structure: Structure): string[] {
    const set = new Set<string>();
    const loc = StructureElement.Location.create(structure);
    for (const unit of structure.units) {
        loc.unit = unit;
        for (let i = 0; i < unit.elements.length; i++) {
            loc.element = unit.elements[i];
            set.add(StructureProperties.chain.auth_asym_id(loc));
        }
    }
    return Array.from(set).sort();
}

function extractAtomPositions(structure: Structure, selectedChainIds: string[]): Float32Array {
    const elementsByUnit = getChainElements(structure, selectedChainIds);
    let count = 0;
    for (const { indices } of elementsByUnit) count += OrderedSet.size(indices);

    const positions = new Float32Array(count * 3);
    let idx = 0;
    const pos = Vec3();
    for (const { unit, indices } of elementsByUnit) {
        const { elements, conformation: c } = unit;
        for (let i = 0, il = OrderedSet.size(indices); i < il; i++) {
            const element = elements[OrderedSet.getAt(indices, i)];
            // Use position() (not invariantPosition()) to apply the unit's operator matrix,
            // giving correct world-space coords for assemblies and symmetry-expanded structures.
            c.position(element, pos);
            positions[idx++] = pos[0];
            positions[idx++] = pos[1];
            positions[idx++] = pos[2];
        }
    }
    return positions;
}

function getChainLoci(structure: Structure, selectedChainIds: string[]) {
    return StructureElement.Loci(structure, getChainElements(structure, selectedChainIds));
}

function getChainElements(structure: Structure, selectedChainIds: string[]) {
    const filter = new Set(selectedChainIds);
    const loc = StructureElement.Location.create(structure);
    const elements: Array<{ unit: Unit, indices: OrderedSet<UnitIndex> }> = [];

    for (const unit of structure.units) {
        if (unit.elements.length === 0) continue;

        const indices: UnitIndex[] = [];
        loc.unit = unit;
        for (let i = 0; i < unit.elements.length; i++) {
            loc.element = unit.elements[i];
            if (filter.has(StructureProperties.chain.auth_asym_id(loc))) {
                indices.push(i as UnitIndex);
            }
        }

        if (indices.length > 0) {
            elements.push({ unit, indices: OrderedSet.ofSortedArray(indices) });
        }
    }

    return elements;
}

/** PluginBehavior that marks the volume-mask extension as active in the plugin. */
export const VolumeMaskBehavior = PluginBehavior.create({
    name: 'volume-mask-behavior',
    category: 'misc',
    display: { name: 'Volume Mask Creator' },
    ctor: class extends PluginBehavior.Handler {
        register() { /* MaskVolumeFromSource is a BuiltIn transformer — registered at module load */ }
        unregister() {}
    },
    params: () => ({}),
});
