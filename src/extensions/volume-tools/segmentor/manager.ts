/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';
import { StatefulPluginComponent } from '../../../mol-plugin-state/component';
import { createVolumeRepresentationParams } from '../../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectCell, StateTransform } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { Color } from '../../../mol-util/color';
import { ColorLists } from '../../../mol-util/color/lists';
import { computeCandidates } from '../candidates';
import { downloadVolumeMrc, exportBodyMasks, maskBaseName, resolveBodyMaskParams } from './internal/export';
import { assignPolygons, assignRemainder, countUnassigned, countVoxels } from './internal/label-ops';
import { padGridBox, voxelBBox } from '../soft-mask';
import { flipHandednessInPlace, removeDustInPlace } from './internal/volume-edit';
import { BodyLabels } from './labels';
import { BodyLabelColorThemeProvider } from './theme';
import { BodyMaskFromLabels, BodyMaskFromLabelsTag } from './transformers';
import { BodyId, BodyInfo, BodyMaskParams, LabelStore, MaxBodyId, ViewMask } from './types';

export type PreviewMode = 'none' | 'active' | 'all';

export interface VolumeSegmentorStats {
    /** Voxels above the threshold. */
    candidates: number
    /** Voxels above the threshold not assigned to any body. */
    unassigned: number
    total: number
}

export interface VolumeSegmentorState {
    targetVolumeRef: StateTransform.Ref | undefined
    threshold: Volume.IsoValue
    isDrawing: boolean
    activeBodyId: BodyId | undefined
    /** Snapshot of the body definitions (new array on every change). */
    bodies: BodyInfo[]
    defaults: BodyMaskParams
    preview: PreviewMode
    undoDepth: number
    stats: VolumeSegmentorStats
    busy: boolean
    /** Mirrors `LabelStore.version`. */
    version: number
    /** Last non-fatal message for the UI (e.g. why a preview was skipped). */
    message: string | undefined
}

const ThresholdPreviewDebounceMs = 75;
const MaxUndoDepth = 50;
const MaxPreviewAllBytes = 512 * 1024 * 1024;
const UnassignedColor = Color(0x9a9a9a);

const BodyPalette: Color[] = ColorLists['many-distinct'].list.map(e => Array.isArray(e) ? e[0] : e);

const managers = new WeakMap<PluginContext, VolumeSegmentorManager>();

function defaultState(): VolumeSegmentorState {
    return {
        targetVolumeRef: undefined,
        threshold: Volume.IsoValue.relative(3),
        isDrawing: false,
        activeBodyId: undefined,
        bodies: [],
        defaults: { extend: 6, softEdge: 2, pruneBelowThreshold: true },
        preview: 'none',
        undoDepth: 0,
        stats: { candidates: 0, unassigned: 0, total: 0 },
        busy: false,
        version: 0,
        message: undefined,
    };
}

/** Sentinel for "no call in flight", so clearing the target is not mistaken for it. */
const NoPendingTarget = Symbol('no-pending-target');

/**
 * Drives interactive body segmentation of one volume. Bodies are defined by their view
 * polygons (plus an optional remainder body); voxel labels are recomputed from these
 * definitions after every change and kept in sync with the plugin state (colored source
 * surface, per-body mask previews).
 */
export class VolumeSegmentorManager extends StatefulPluginComponent<VolumeSegmentorState> {
    static get(plugin: PluginContext): VolumeSegmentorManager | undefined {
        return managers.get(plugin);
    }
    static register(plugin: PluginContext, manager: VolumeSegmentorManager) {
        managers.set(plugin, manager);
    }
    static unregister(plugin: PluginContext) {
        managers.delete(plugin);
    }

    readonly behaviors = {
        state: this.ev.behavior<VolumeSegmentorState>(this.state),
    };

    private candidates: Int32Array = new Int32Array(0);
    private undoStack: BodyInfo[][] = [];
    private previewRefs = new Map<BodyId, { mask: StateTransform.Ref, repr: StateTransform.Ref }>();
    private themeSyncInFlight = false;
    private pendingTargetRef: StateTransform.Ref | undefined | typeof NoPendingTarget = NoPendingTarget;
    private thresholdTimer: ReturnType<typeof setTimeout> | undefined;
    private thresholdInFlight = false;
    private pendingThreshold: Volume.IsoValue | undefined;
    private recomputing = false;
    private recomputeQueued = false;

    constructor(readonly plugin: PluginContext) {
        super(defaultState());
        // Representations can be created after the volume node; keep every isosurface of the
        // target colored by body label.
        this.subscribe(plugin.state.data.events.changed, () => {
            if (this.state.targetVolumeRef && !this.themeSyncInFlight) void this.applyThemeToVolume();
        });
    }

    private update(patch: Partial<VolumeSegmentorState>) {
        if (this.updateState(patch)) this.behaviors.state.next(this.state);
    }

    // --- data access ---

    get volume(): Volume | undefined {
        const ref = this.state.targetVolumeRef;
        if (!ref) return undefined;
        return (this.plugin.state.data.select(ref)[0]?.obj as SO.Volume.Data | undefined)?.data;
    }

    get store(): LabelStore | undefined {
        const volume = this.volume;
        return volume && BodyLabels.get(volume);
    }

    get thresholdAbs(): number {
        const volume = this.volume;
        if (!volume) return 0;
        return Volume.IsoValue.toAbsolute(this.state.threshold, volume.grid.stats).absoluteValue;
    }

    get activeBody(): BodyInfo | undefined {
        return this.state.bodies.find(b => b.id === this.state.activeBodyId);
    }

    /** Direct isosurface representations of the target volume (previews excluded). */
    private get volumeIsosurfaceReprs(): StateObjectCell[] {
        const ref = this.state.targetVolumeRef;
        return ref ? this.volumeIsosurfaceReprsOf(ref) : [];
    }

    private volumeIsosurfaceReprsOf(ref: StateTransform.Ref): StateObjectCell[] {
        return this.volumeReprsOf(ref).filter(c => (c.transform.params as any)?.type?.name === 'isosurface');
    }

    /** All representations of the volume at `ref` (previews excluded). */
    private volumeReprsOf(ref: StateTransform.Ref): StateObjectCell[] {
        const state = this.plugin.state.data;
        return state
            .selectQ(q => q.ofTransformer(StateTransforms.Representation.VolumeRepresentation3D, ref))
            .filter(c => {
                // The nearest (non-decorator) volume ancestor must be the target itself:
                // representations of the body mask previews sit under their own derived volumes.
                let parent = state.cells.get(c.transform.parent);
                while (parent && parent.transform.ref !== ref) {
                    if (SO.Volume.Data.is(parent.obj) && !parent.transform.transformer.definition.isDecorator) return false;
                    parent = state.cells.get(parent.transform.parent);
                }
                return !!parent;
            });
    }

    // --- target volume ---

    async setTargetVolume(ref: StateTransform.Ref | undefined) {
        if (this.pendingTargetRef === ref) return;
        this.pendingTargetRef = ref;
        try {
            await this.applyTargetVolume(ref);
        } finally {
            this.pendingTargetRef = NoPendingTarget;
        }
    }

    private async applyTargetVolume(ref: StateTransform.Ref | undefined) {
        this.clearThresholdTimer();
        await this.removePreviews();
        this.undoStack = [];

        if (!ref) {
            this.candidates = new Int32Array(0);
            this.update({ ...defaultState() });
            return;
        }

        const volume = (this.plugin.state.data.select(ref)[0]?.obj as SO.Volume.Data | undefined)?.data;
        if (!volume) return;

        const store = BodyLabels.ensure(volume);
        const { mean, sigma } = volume.grid.stats;
        // Adopt the iso value of an existing isosurface so the surface and the threshold agree.
        const existingIso = (this.volumeIsosurfaceReprsOf(ref)[0]?.transform.params as any)?.type?.params?.isoValue as Volume.IsoValue | undefined;
        const threshold = existingIso ?? Volume.IsoValue.absolute(mean + 3 * sigma);
        this.candidates = computeCandidates(volume, Volume.IsoValue.toAbsolute(threshold, volume.grid.stats).absoluteValue);

        this.update({
            ...defaultState(),
            targetVolumeRef: ref,
            threshold,
            bodies: store.bodies.slice(),
            activeBodyId: store.bodies[0]?.id,
            version: store.version,
            stats: this.computeStats(store),
        });
        await this.applyThemeToVolume();
        if (store.bodies.length > 0) await this.recompute();
    }

    private computeStats(store: LabelStore): VolumeSegmentorStats {
        return {
            candidates: this.candidates.length,
            unassigned: countUnassigned(store.labels, this.candidates),
            total: store.labels.length,
        };
    }

    /**
     * Colors the source isosurfaces by body label and disables wrapping across box faces.
     * Idempotent: only representations that are not set up yet are touched.
     */
    private async applyThemeToVolume() {
        const version = this.store?.version ?? 0;
        const reprs = this.volumeIsosurfaceReprs.filter(cell => {
            const params = cell.transform.params as any;
            return params?.colorTheme?.name !== BodyLabelColorThemeProvider.name || params?.type?.params?.wrap !== 'off';
        });
        if (reprs.length === 0) return;
        const isoValue = this.state.threshold;
        this.themeSyncInFlight = true;
        try {
            const builder = this.plugin.build();
            for (const cell of reprs) {
                builder.to(cell).update(StateTransforms.Representation.VolumeRepresentation3D, old => ({
                    ...old,
                    type: { name: 'isosurface', params: { ...old.type.params, wrap: 'off', isoValue } },
                    colorTheme: { name: BodyLabelColorThemeProvider.name, params: { unassignedColor: UnassignedColor, version } },
                }));
            }
            await builder.commit({ canUndo: false });
        } finally {
            this.themeSyncInFlight = false;
        }
    }

    private async refreshTheme() {
        const reprs = this.volumeIsosurfaceReprs;
        if (reprs.length === 0) return;
        const version = this.store?.version ?? 0;
        const builder = this.plugin.build();
        for (const cell of reprs) {
            builder.to(cell).update(StateTransforms.Representation.VolumeRepresentation3D, old => ({
                ...old,
                colorTheme: old.colorTheme.name === BodyLabelColorThemeProvider.name
                    ? { name: old.colorTheme.name, params: { ...old.colorTheme.params, version } }
                    : { name: BodyLabelColorThemeProvider.name, params: { unassignedColor: UnassignedColor, version } },
            }));
        }
        await builder.commit({ canUndo: false });
    }

    // --- threshold ---

    /** Sets the threshold; the isosurface follows with a short debounce and bodies are recomputed. */
    setThreshold(value: Volume.IsoValue) {
        if (!this.state.targetVolumeRef) return;
        this.update({ threshold: value });
        this.pendingThreshold = value;
        this.clearThresholdTimer();
        this.thresholdTimer = setTimeout(() => {
            this.thresholdTimer = undefined;
            void this.flushThreshold();
        }, ThresholdPreviewDebounceMs);
    }

    private clearThresholdTimer() {
        if (this.thresholdTimer !== undefined) {
            clearTimeout(this.thresholdTimer);
            this.thresholdTimer = undefined;
        }
    }

    private async flushThreshold() {
        if (this.thresholdInFlight) return;
        this.thresholdInFlight = true;
        try {
            while (this.pendingThreshold) {
                const value = this.pendingThreshold;
                this.pendingThreshold = undefined;
                await this.applyThreshold(value);
            }
        } finally {
            this.thresholdInFlight = false;
            if (this.pendingThreshold) void this.flushThreshold();
        }
    }

    private async applyThreshold(value: Volume.IsoValue) {
        const volume = this.volume;
        if (!volume) return;

        const reprs = this.volumeIsosurfaceReprs;
        if (reprs.length > 0) {
            const builder = this.plugin.build();
            for (const cell of reprs) {
                builder.to(cell).update(StateTransforms.Representation.VolumeRepresentation3D, old => ({
                    ...old,
                    type: { name: 'isosurface', params: { ...old.type.params, isoValue: value } },
                }));
            }
            await builder.commit({ canUndo: false });
        }

        this.candidates = computeCandidates(volume, Volume.IsoValue.toAbsolute(value, volume.grid.stats).absoluteValue);
        await this.recompute();
    }

    /**
     * Zeroes isolated blobs smaller than `minVoxels` (at the current threshold) in the source
     * volume itself, then rebuilds its representations and recomputes the bodies. Returns the
     * number of voxels zeroed. This edit cannot be undone.
     */
    async removeDust(minVoxels: number): Promise<number> {
        const ref = this.state.targetVolumeRef;
        const volume = this.volume;
        if (!ref || !volume) return 0;

        const thresholdAbs = this.thresholdAbs;
        const zeroed = await this.runBusy(Task.create('Remove dust', async () => removeDustInPlace(volume, minVoxels, thresholdAbs)));
        if (zeroed === 0) return 0;

        this.candidates = computeCandidates(volume, thresholdAbs);
        await this.rebuildVolumeRepresentations(ref);
        await this.recompute();
        return zeroed;
    }

    /**
     * Mirrors the source volume along X to fix a map stored with the opposite handedness, then
     * rebuilds its representations and recomputes the bodies. View polygons are kept, so bodies
     * are relabelled against the mirrored data. This edit cannot be undone.
     */
    async flipHandedness() {
        const ref = this.state.targetVolumeRef;
        const volume = this.volume;
        if (!ref || !volume) return;

        await this.runBusy(Task.create('Flip handedness', async () => flipHandednessInPlace(volume)));
        this.candidates = computeCandidates(volume, this.thresholdAbs);
        await this.rebuildVolumeRepresentations(ref);
        await this.recompute();
    }

    /** Downloads the source volume as MRC, including dust removal and handedness flips. */
    saveVolume() {
        const volume = this.volume;
        if (!volume) return;
        downloadVolumeMrc(volume, maskBaseName(volume.label));
    }

    /** Representations compare volumes by grid reference, so in-place edits need a rebuild. */
    private async rebuildVolumeRepresentations(ref: StateTransform.Ref) {
        const reprs = this.volumeReprsOf(ref);
        if (reprs.length === 0) return;
        const saved = reprs.map(c => ({ transformer: c.transform.transformer, params: c.transform.params }));

        const remove = this.plugin.build();
        for (const c of reprs) remove.delete(c.transform.ref);
        await remove.commit({ canUndo: false });

        const add = this.plugin.build();
        for (const s of saved) add.to(ref).apply(s.transformer as any, s.params);
        await add.commit({ canUndo: false });
    }

    // --- body definitions ---

    private snapshot(store: LabelStore) {
        this.undoStack.push(store.bodies.slice());
        if (this.undoStack.length > MaxUndoDepth) this.undoStack.shift();
    }

    private async setBodies(store: LabelStore, bodies: BodyInfo[], recompute = true) {
        store.bodies = bodies;
        this.update({ bodies: bodies.slice(), undoDepth: this.undoStack.length });
        if (recompute) await this.recompute();
    }

    addBody(name?: string, color?: Color): BodyInfo | undefined {
        const store = this.store;
        if (!store) return undefined;
        if (store.nextId > MaxBodyId) {
            this.update({ message: `At most ${MaxBodyId} bodies are supported.` });
            return undefined;
        }
        this.snapshot(store);
        const id = store.nextId++;
        const body: BodyInfo = {
            id,
            name: name ?? `Body ${id}`,
            color: color ?? BodyPalette[(id - 1) % BodyPalette.length],
            views: [],
            remainder: false,
            voxelCount: 0,
        };
        void this.setBodies(store, [...store.bodies, body], false);
        this.update({ activeBodyId: id, message: undefined });
        return body;
    }

    /** Adds (or activates) the body that collects every voxel no other body claims. */
    async addRemainderBody(name?: string) {
        const store = this.store;
        if (!store) return;
        const existing = store.bodies.find(b => b.remainder);
        if (existing) {
            this.update({ activeBodyId: existing.id });
            return;
        }
        const body = this.addBody(name ?? 'Remainder');
        if (!body) return;
        await this.setBodies(store, store.bodies.map(b => b.id === body.id ? { ...b, remainder: true } : b));
    }

    async removeBody(id: BodyId) {
        const store = this.store;
        if (!store) return;
        this.snapshot(store);
        await this.removePreview(id);
        const bodies = store.bodies.filter(b => b.id !== id);
        if (this.state.activeBodyId === id) this.update({ activeBodyId: bodies[0]?.id });
        await this.setBodies(store, bodies);
    }

    /** Moves a body up (`-1`) or down (`+1`) in the priority order. */
    async moveBody(id: BodyId, delta: -1 | 1) {
        const store = this.store;
        if (!store) return;
        const index = store.bodies.findIndex(b => b.id === id);
        const target = index + delta;
        if (index < 0 || target < 0 || target >= store.bodies.length) return;
        this.snapshot(store);
        const bodies = store.bodies.slice();
        [bodies[index], bodies[target]] = [bodies[target], bodies[index]];
        await this.setBodies(store, bodies);
    }

    setActiveBody(id: BodyId | undefined) {
        this.update({ activeBodyId: id });
        void this.syncPreviews();
    }

    renameBody(id: BodyId, name: string) {
        const store = this.store;
        if (!store) return;
        void this.setBodies(store, store.bodies.map(b => b.id === id ? { ...b, name } : b), false);
    }

    async setBodyColor(id: BodyId, color: Color) {
        const store = this.store;
        if (!store) return;
        store.bodies = store.bodies.map(b => b.id === id ? { ...b, color } : b);
        BodyLabels.bump(store);
        this.update({ bodies: store.bodies.slice(), version: store.version });
        await this.refreshTheme();
        await this.syncPreviews();
    }

    /** Per-body extend / soft edge; `undefined` values fall back to the defaults. */
    async setBodyMaskParams(id: BodyId, params: { extend?: number, softEdge?: number }) {
        const store = this.store;
        if (!store) return;
        store.bodies = store.bodies.map(b => b.id === id ? { ...b, extend: params.extend, softEdge: params.softEdge } : b);
        this.update({ bodies: store.bodies.slice() });
        await this.syncPreviews();
    }

    async setDefaults(defaults: BodyMaskParams) {
        this.update({ defaults });
        await this.syncPreviews();
    }

    // --- views ---

    setDrawing(isDrawing: boolean) {
        this.update({ isDrawing });
    }

    async addView(bodyId: BodyId, view: ViewMask) {
        const store = this.store;
        if (!store) return;
        this.snapshot(store);
        await this.setBodies(store, store.bodies.map(b => b.id === bodyId ? { ...b, views: [...b.views, view] } : b));
    }

    async removeView(bodyId: BodyId, viewId: string) {
        const store = this.store;
        if (!store) return;
        this.snapshot(store);
        await this.setBodies(store, store.bodies.map(b => b.id === bodyId ? { ...b, views: b.views.filter(v => v.id !== viewId) } : b));
    }

    /** Toggles a view between selecting inside and outside its polygon. */
    async invertView(bodyId: BodyId, viewId: string) {
        const store = this.store;
        if (!store) return;
        this.snapshot(store);
        await this.setBodies(store, store.bodies.map(b => b.id === bodyId
            ? { ...b, views: b.views.map(v => v.id === viewId ? { ...v, inverted: !v.inverted } : v) }
            : b));
    }

    flyTo(view: ViewMask, durationMs = 400) {
        this.plugin.canvas3d?.camera.setState(view.cameraSnapshot, durationMs);
    }

    async undo() {
        const store = this.store;
        const bodies = this.undoStack.pop();
        if (!store || !bodies) return;
        const activeBodyId = bodies.some(b => b.id === this.state.activeBodyId) ? this.state.activeBodyId : bodies[0]?.id;
        this.update({ activeBodyId });
        await this.setBodies(store, bodies);
    }

    /** Removes all bodies. */
    async reset() {
        const store = this.store;
        if (!store) return;
        this.snapshot(store);
        await this.removePreviews();
        this.update({ activeBodyId: undefined, isDrawing: false });
        await this.setBodies(store, []);
    }

    // --- labels ---

    /** Recomputes all voxel labels from the body definitions (list order = priority). */
    async recompute() {
        if (this.recomputing) {
            this.recomputeQueued = true;
            return;
        }
        this.recomputing = true;
        try {
            do {
                this.recomputeQueued = false;
                await this.recomputeOnce();
            } while (this.recomputeQueued);
        } finally {
            this.recomputing = false;
        }
    }

    private async recomputeOnce() {
        const volume = this.volume;
        const store = this.store;
        if (!volume || !store) return;
        const { labels } = store;
        const candidates = this.candidates;
        const bodies = store.bodies;

        await this.runBusy(Task.create('Update bodies', async ctx => {
            labels.fill(0);
            for (const body of bodies) {
                if (body.remainder || body.views.length === 0) continue;
                await assignPolygons(labels, candidates, volume, body.views, body.id, 'unassigned-only', ctx);
            }
            const remainder = bodies.find(b => b.remainder);
            if (remainder) assignRemainder(labels, candidates, remainder.id);
        }));

        const counts = countVoxels(labels);
        store.bodies = store.bodies.map(b => b.voxelCount === counts[b.id] ? b : { ...b, voxelCount: counts[b.id] });
        BodyLabels.bump(store);
        this.update({
            bodies: store.bodies.slice(),
            version: store.version,
            undoDepth: this.undoStack.length,
            stats: this.computeStats(store),
        });
        await this.refreshTheme();
        await this.syncPreviews();
    }

    private async runBusy<T>(task: Task<T>): Promise<T> {
        this.update({ busy: true });
        try {
            return await this.plugin.runTask(task);
        } finally {
            this.update({ busy: false });
        }
    }

    // --- mask previews ---

    async setPreview(preview: PreviewMode) {
        this.update({ preview, message: undefined });
        await this.syncPreviews();
    }

    private previewTargets(): BodyId[] {
        const { preview, bodies, activeBodyId } = this.state;
        if (preview === 'none') return [];
        if (preview === 'active') {
            const active = bodies.find(b => b.id === activeBodyId);
            return active && active.voxelCount > 0 ? [active.id] : [];
        }
        return bodies.filter(b => b.voxelCount > 0).map(b => b.id);
    }

    private estimatePreviewBytes(ids: BodyId[]): number {
        const volume = this.volume;
        const store = this.store;
        if (!volume || !store) return 0;
        const { space } = volume.grid.cells;
        let bytes = 0;
        for (const id of ids) {
            const body = store.bodies.find(b => b.id === id);
            const tight = voxelBBox(store.labels, space, id);
            if (!body || !tight) continue;
            const params = resolveBodyMaskParams(body, this.state.defaults);
            const box = padGridBox(tight, params.extend + params.softEdge + 1, space.dimensions);
            bytes += box.dims[0] * box.dims[1] * box.dims[2] * 4;
        }
        return bytes;
    }

    /** Creates, updates or removes per-body mask preview nodes to match the current state. */
    private async syncPreviews() {
        const ref = this.state.targetVolumeRef;
        const store = this.store;
        if (!ref || !store) return;

        let targets = this.previewTargets();
        if (this.state.preview === 'all' && this.estimatePreviewBytes(targets) > MaxPreviewAllBytes) {
            this.update({ message: 'Previewing all bodies would need too much memory; showing the active body instead.' });
            const active = this.activeBody;
            targets = active && active.voxelCount > 0 ? [active.id] : [];
        }

        for (const id of Array.from(this.previewRefs.keys())) {
            if (!targets.includes(id)) await this.removePreview(id);
        }

        for (const id of targets) {
            const body = store.bodies.find(b => b.id === id);
            if (!body) continue;
            const params = resolveBodyMaskParams(body, this.state.defaults);
            const transformParams = {
                bodyId: id,
                extend: params.extend,
                softEdge: params.softEdge,
                pruneBelowThreshold: params.pruneBelowThreshold,
                threshold: this.state.threshold,
                version: store.version,
            };

            const existing = this.previewRefs.get(id);
            if (existing && this.plugin.state.data.cells.has(existing.mask)) {
                const cell = this.plugin.state.data.cells.get(existing.mask)!;
                const old = cell.transform.params as typeof transformParams;
                const paramsChanged = old.version !== transformParams.version || old.extend !== transformParams.extend
                    || old.softEdge !== transformParams.softEdge || old.pruneBelowThreshold !== transformParams.pruneBelowThreshold
                    || !Volume.IsoValue.areSame(old.threshold, transformParams.threshold, this.volume!.grid.stats);
                const builder = this.plugin.build();
                if (paramsChanged) builder.to(existing.mask).update(BodyMaskFromLabels, () => transformParams);
                builder.to(existing.repr).update(StateTransforms.Representation.VolumeRepresentation3D, o => ({
                    ...o,
                    colorTheme: { name: 'uniform', params: { value: body.color } },
                }));
                await builder.commit({ canUndo: false });
                continue;
            }

            const mask = await this.plugin.build().to(ref)
                .apply(BodyMaskFromLabels, transformParams, { tags: BodyMaskFromLabelsTag })
                .commit({ canUndo: false });
            if (!mask.data) continue;

            const reprParams = createVolumeRepresentationParams(this.plugin, mask.data, {
                type: 'isosurface',
                typeParams: { isoValue: Volume.IsoValue.absolute(0.5), alpha: 0.45, wrap: 'off' },
                color: 'uniform',
                colorParams: { value: body.color },
            });
            const repr = await this.plugin.build().to(mask.ref)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, reprParams, { tags: BodyMaskFromLabelsTag })
                .commit({ canUndo: false });
            this.previewRefs.set(id, { mask: mask.ref, repr: repr.ref });
        }
    }

    private async removePreview(id: BodyId) {
        const refs = this.previewRefs.get(id);
        this.previewRefs.delete(id);
        if (!refs) return;
        if (this.plugin.state.data.cells.has(refs.mask)) {
            await this.plugin.build().delete(refs.mask).commit({ canUndo: false });
        }
    }

    private async removePreviews() {
        for (const id of Array.from(this.previewRefs.keys())) await this.removePreview(id);
        // Preview nodes may also exist from an earlier manager instance on the same volume.
        const ref = this.state.targetVolumeRef;
        if (!ref) return;
        const leftovers = this.plugin.state.data.selectQ(q => q.ofTransformer(BodyMaskFromLabels, ref));
        if (leftovers.length === 0) return;
        const builder = this.plugin.build();
        for (const cell of leftovers) builder.delete(cell.transform.ref);
        await builder.commit({ canUndo: false });
    }

    // --- export ---

    /** Downloads one MRC mask per body (largest first), optionally bundled as a zip. */
    async exportMasks(options: { bundleZip: boolean, compress: boolean }) {
        const volume = this.volume;
        const store = this.store;
        if (!volume || !store) return [];
        const baseName = maskBaseName(volume.label);
        return this.runBusy(Task.create('Export body masks', ctx =>
            exportBodyMasks(volume, store, this.state.defaults, this.thresholdAbs, { baseName, ...options }, ctx)
        ));
    }

    dispose() {
        this.clearThresholdTimer();
        this.pendingThreshold = undefined;
        void this.removePreviews();
        super.dispose();
    }
}

/** True for state cells created by the body mask preview, so volume pickers can skip them. */
export function isBodyMaskCell(cell: StateObjectCell) {
    return cell.transform.transformer === BodyMaskFromLabels;
}
