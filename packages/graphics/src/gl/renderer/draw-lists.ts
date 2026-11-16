/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The renderer's draw lists: a pass of a scene's table recorded in chunks (WebGPU render bundles) and replayed while
// what it draws is unchanged, recorded again where it changed (`GraphicsContext.drawLists`)

import type { GraphicsContext } from '../context.js';
import { DrawFact, DrawRow, type DrawTable } from '../draw-table.js';
import type { DrawList, GlobalBindings, GraphicsRenderVariant } from '../render-item.js';
import type { RenderState } from '../render-state.js';
import { DrawKey, type GraphicsRenderable } from '../renderable.js';
import type { Scene } from '../scene.js';
import { getObjectState, isPassDrawn, type PassDraws, type SolidInteriorMode } from './pass-draws.js';

/**
 * A part of a pass over a table's order: the objects of `[begin, end)` with their inline caps (`objects`), the caps
 * drawn after all objects of the range (`caps`); a chunk of a draw list's recordings
 */
export type PassPart = {
  readonly begin: number;
  readonly end: number;
  readonly objects: boolean;
  readonly caps: boolean;
};
export const WholePass: PassPart = { begin: 0, end: Infinity, objects: true, caps: true };

/**
 * What the draw lists read of the renderer: how it draws a part of a pass (`renderPassRows`) and culls, and the state
 * of the frame a recording depends on besides the table (getters: read when drawn)
 */
export interface DrawListHost {
  readonly ctx: GraphicsContext;
  readonly globals: GlobalBindings;
  /** the culling view's planes, as the renderer compares them (`isSameView`) */
  readonly viewPlanes: Float64Array;
  /** per cap of the last `getSolidInteriorCaps` the segments it cuts */
  readonly capSegments: readonly (readonly number[])[];
  readonly modelScale: number;
  readonly near: number;
  /** advances where GPU culling's records change */
  readonly cullRecords: number;
  /** objects culled on the GPU */
  readonly gpuCulledCount: number;
  /** an occlusion test set (Hi-Z, changing between frames) */
  readonly occlusionTest: boolean;
  /** drawing the marked variant (`renderMarked`) */
  readonly drawingMarked: boolean;
  renderPassRows: (pass: PassDraws, table: DrawTable, part?: PassPart) => void;
  getRowVisibility: (table: DrawTable) => ArrayLike<number>;
  setRowBoundingSphere: (table: DrawTable, slot: number) => void;
  cullObject: (r: GraphicsRenderable) => void;
  /** whether the sphere raster draws `r` instead of the pass (`sphereRaster` 'only') */
  isRasterDrawn: (r: GraphicsRenderable) => boolean;
  /** whether solid interior caps are drawn */
  isCapsDrawn: () => boolean;
  hasSolidInteriorCap: (facts: number) => boolean;
  getSolidInteriorCaps: (r: GraphicsRenderable, mode: SolidInteriorMode, out: number[]) => readonly number[];
  isSameView: (planes: Float64Array) => boolean;
  /** set while a pass draws from a draw list (its rows unculled while the view moves) */
  setDrawingFromList: (value: boolean) => void;
}

export interface RendererDrawLists {
  /** Draws `pass` of `group` from its draw list, recording what changed; false where it can't (no table or lists) */
  render: (pass: PassDraws, group: Scene.Group) => boolean;
  /** Drops all recordings */
  reset: () => void;
}

export function createRendererDrawLists(host: DrawListHost): RendererDrawLists {
  const {
    ctx,
    globals,
    viewPlanes,
    capSegments,
    renderPassRows,
    getRowVisibility,
    setRowBoundingSphere,
    cullObject,
    isCapsDrawn,
    hasSolidInteriorCap,
    getSolidInteriorCaps,
    isSameView,
  } = host;

  /**
   * What a recording drew, by the rows of its table (`walkPass`): per draw its slot, object, render state and the
   * versions of its program, values (`DrawTable.valuesVersions`) and draws (`Renderable.drawVersion`), and what it
   * binds (`digests`, for data changes); the caps as a flat list (per object drawing caps: slot, program and values
   * versions, count, caps); and what the walk read besides the rows, which unchanged with the table means the same
   * draws (`isSnapshotCurrent`). Written in place, so that comparing a frame allocates nothing
   */
  type Snapshot = {
    valid: boolean;
    tableVersion: number;
    valuesEpoch: number;
    readonly view: Float64Array;
    modelScale: number;
    near: number;
    caps: boolean;
    /** `DrawLists.resourcesVersion` */
    drawResources: number;
    cullRecords: number;
    count: number;
    slots: Uint32Array;
    programVersions: Uint32Array;
    valuesVersions: Uint32Array;
    drawVersions: Uint32Array;
    readonly objects: GraphicsRenderable[];
    readonly states: RenderState[];
    readonly digests: (string | undefined)[];
    readonly capKey: DrawKey;
  };
  const createSnapshot = (): Snapshot => ({
    valid: false,
    tableVersion: -1,
    valuesEpoch: -1,
    view: new Float64Array(28).fill(NaN),
    modelScale: NaN,
    near: NaN,
    caps: false,
    drawResources: -1,
    cullRecords: 0,
    count: 0,
    slots: new Uint32Array(0),
    programVersions: new Uint32Array(0),
    valuesVersions: new Uint32Array(0),
    drawVersions: new Uint32Array(0),
    objects: [],
    states: [],
    digests: [],
    capKey: DrawKey(),
  });
  const ensureSnapshotDraws = (snap: Snapshot, n: number) => {
    if (snap.slots.length > n) return;
    const grow = (a: Uint32Array) => {
      const b = new Uint32Array(Math.max(16, n * 2));
      b.set(a);
      return b;
    };
    snap.slots = grow(snap.slots);
    snap.programVersions = grow(snap.programVersions);
    snap.valuesVersions = grow(snap.valuesVersions);
    snap.drawVersions = grow(snap.drawVersions);
  };
  const setSnapshotInputs = (snap: Snapshot, table: DrawTable) => {
    snap.tableVersion = table.version;
    snap.valuesEpoch = table.valuesEpoch;
    snap.view.set(viewPlanes);
    snap.modelScale = host.modelScale;
    snap.near = host.near;
    snap.caps = isCapsDrawn();
    snap.drawResources = ctx.drawLists?.resourcesVersion ?? 0;
    snap.cullRecords = host.cullRecords;
  };
  /**
   * Whether nothing the walk read changed (the rows and order, values, programs, the view and model scale, the near
   * plane, whether caps are drawn): the snapshot holds without walking. Not with an occlusion test (Hi-Z changes
   * between frames)
   */
  const isSnapshotCurrent = (snap: Snapshot, table: DrawTable) => {
    if (!snap.valid || host.occlusionTest) return false;
    // volatile cells compared
    table.refresh();
    return snap.tableVersion === table.version && snap.valuesEpoch === table.valuesEpoch && isSameInputs(snap);
  };
  /** Whether what the walk read besides the table is unchanged (see `isSnapshotCurrent`) */
  const isSameInputs = (snap: Snapshot) =>
    isSameView(snap.view) &&
    snap.modelScale === host.modelScale &&
    snap.near === host.near &&
    snap.caps === isCapsDrawn() &&
    snap.drawResources === (ctx.drawLists?.resourcesVersion ?? 0) &&
    snap.cullRecords === host.cullRecords;

  const walkCapList: number[] = [];
  /**
   * The caps `r` (in `slot`) draws, at `c` of the snapshot's cap list: written, or compared (-1 where they differ);
   * nothing where it draws none (one starting or stopping to draw caps changes the list all the same). Assumes `r`
   * is visible, `boundingSphere` set; returns the position after its entry
   */
  const walkCaps = (
    r: GraphicsRenderable,
    slot: number,
    valuesVersion: number,
    mode: SolidInteriorMode,
    key: DrawKey,
    c: number,
    write: boolean,
  ) => {
    const caps = getSolidInteriorCaps(r, mode, walkCapList);
    if (caps.length === 0) return c;
    if (write) {
      DrawKey.add(key, slot);
      DrawKey.add(key, r.programVersion);
      DrawKey.add(key, valuesVersion);
      DrawKey.add(key, caps.length);
      for (let k = 0, kl = caps.length; k < kl; ++k) DrawKey.add(key, caps[k]);
      // the segments each cap draws (`capSegments`): a recording replays while the same are cut
      for (let k = 0, kl = caps.length; k < kl; ++k) {
        const segments = capSegments[k];
        DrawKey.add(key, segments.length);
        for (let i = 0, il = segments.length; i < il; ++i) DrawKey.add(key, segments[i]);
      }
      return key.length;
    }
    const v = key.values;
    if (c + 4 + caps.length > key.length) return -1;
    if (v[c] !== slot || v[c + 1] !== r.programVersion || v[c + 2] !== valuesVersion || v[c + 3] !== caps.length) {
      return -1;
    }
    for (let k = 0, kl = caps.length; k < kl; ++k) {
      if (v[c + 4 + k] !== caps[k]) return -1;
    }
    let o = c + 4 + caps.length;
    for (let k = 0, kl = caps.length; k < kl; ++k) {
      const segments = capSegments[k];
      if (o + 1 + segments.length > key.length || v[o] !== segments.length) return -1;
      for (let i = 0, il = segments.length; i < il; ++i) if (v[o + 1 + i] !== segments[i]) return -1;
      o += 1 + segments.length;
    }
    return o;
  };

  /**
   * The draws of `pass` as `renderPassRows` draws them (culled the same way) with the caps they draw (the near plane
   * and clip objects cutting objects, camera dependent), written into `snap` or compared with it: 'same'; 'data'
   * where only values of draws changed (their data may have, see `flushDataChanges`); 'different' at the first draw,
   * state, program, range or cap that differs
   */
  const walkPass = (
    pass: PassDraws,
    table: DrawTable,
    snap: Snapshot,
    write: boolean,
    part: PassPart = WholePass,
  ): 'same' | 'data' | 'different' => {
    const { variant, select, check, draws, caps } = pass;
    const withCaps = !!caps && isCapsDrawn();
    const visible = getRowVisibility(table);
    const { order, bits, objects, valuesVersions } = table;
    const begin = part.begin;
    const end = Math.min(part.end, table.count);
    if (write) snap.capKey.length = 0;
    let n = 0;
    let c = 0;
    let data = false;
    for (let k = part.objects ? begin : end; k < end; ++k) {
      const slot = order[k];
      if (!visible[slot]) continue;
      const o = slot * DrawRow.Stride;
      const facts = bits[o + DrawRow.Facts];
      const state = bits[o + DrawRow.State];
      if (!isPassDrawn(pass, state) || (select && !select(facts, state)) || !check(facts, state)) continue;
      const r = objects[slot]!;
      setRowBoundingSphere(table, slot);
      const objectDraws = draws(facts);
      // as `drawObject`: culled, the raster's objects left out (`sphereRaster` 'only')
      if (objectDraws.length > 0) cullObject(r);
      if (objectDraws.length > 0 && !host.isRasterDrawn(r)) {
        const { programVersion } = r;
        const drawVersion = host.drawingMarked && r.renderMarked ? r.markedDrawVersion! : r.drawVersion;
        const valuesVersion = valuesVersions[slot];
        for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
          const renderState = getObjectState(facts, variant, objectDraws[j].flag, objectDraws[j].passState);
          if (!renderState) continue;
          if (write) {
            ensureSnapshotDraws(snap, n);
            snap.slots[n] = slot;
            snap.objects[n] = r;
            snap.states[n] = renderState;
            snap.programVersions[n] = programVersion;
            snap.valuesVersions[n] = valuesVersion;
            snap.drawVersions[n] = drawVersion;
            snap.digests[n] = undefined;
          } else {
            if (
              n >= snap.count ||
              snap.slots[n] !== slot ||
              snap.objects[n] !== r ||
              snap.states[n] !== renderState ||
              snap.programVersions[n] !== programVersion ||
              snap.drawVersions[n] !== drawVersion
            ) {
              return 'different';
            }
            if (snap.valuesVersions[n] !== valuesVersion) data = true;
          }
          n += 1;
        }
      }
      if (withCaps && caps!.inline && hasSolidInteriorCap(facts) && caps!.check(facts, state)) {
        c = walkCaps(r, slot, valuesVersions[slot], caps!.mode, snap.capKey, c, write);
        if (c < 0) return 'different';
      }
    }
    if (part.caps && withCaps && !caps!.inline) {
      for (let k = begin; k < end; ++k) {
        const slot = order[k];
        if (!visible[slot]) continue;
        const o = slot * DrawRow.Stride;
        const facts = bits[o + DrawRow.Facts];
        const state = bits[o + DrawRow.State];
        if (!isPassDrawn(pass, state) || !hasSolidInteriorCap(facts) || !caps!.check(facts, state)) continue;
        const r = objects[slot]!;
        setRowBoundingSphere(table, slot);
        c = walkCaps(r, slot, valuesVersions[slot], caps!.mode, snap.capKey, c, write);
        if (c < 0) return 'different';
      }
    }
    if (write) {
      snap.count = n;
      // shorter only where draws were dropped: no references kept to their objects
      if (snap.objects.length > n) {
        snap.objects.length = n;
        snap.states.length = n;
        snap.digests.length = n;
      }
      return 'same';
    }
    if (n !== snap.count || c !== snap.capKey.length) return 'different';
    return data ? 'data' : 'same';
  };

  /** What the `i`th draw of `snap` binds, its data written as by drawing it (`DrawLists.capture`) */
  const captureDraw = (snap: Snapshot, i: number, variant: GraphicsRenderVariant) => {
    return ctx.drawLists!.capture(() => {
      ctx.setRenderState(snap.states[i]);
      const r = snap.objects[i];
      if (host.drawingMarked && r.renderMarked) r.renderMarked(variant, globals);
      else r.render(variant, globals);
    });
  };

  /**
   * The draws of `snap` whose values changed written as by drawing them, without encoding them: true where each binds
   * what it bound when recorded (only data changed, e.g. markers on hover), so the recording replays with it; false
   * where one binds otherwise (a texture or buffer replaced, a draw count changed): record again
   */
  const flushDataChanges = (snap: Snapshot, table: DrawTable, variant: GraphicsRenderVariant) => {
    if (!ctx.drawLists) return false;
    const { valuesVersions } = table;
    for (let i = 0, il = snap.count; i < il; ++i) {
      const version = valuesVersions[snap.slots[i]];
      if (snap.valuesVersions[i] === version) continue;
      const digest = snap.digests[i];
      if (digest === undefined || captureDraw(snap, i, variant) !== digest) return false;
      snap.valuesVersions[i] = version;
    }
    return true;
  };

  /**
   * The draw lists of a group's pass, one per time it is drawn between two submissions (stereo eyes, sub-renders, the
   * pick types, marking samples): the n-th draw of the pass after a submission binds the n-th global slot each time
   * (`DrawLists.submitCount`), whether a frame or a marking redraw without a new frame; per pass, so that
   * submissions drawing other passes (a pick, a marking redraw) don't shift the lists of a pass
   */
  type PassDrawLists = {
    submit: number;
    next: number;
    lists: PassDrawList[];
  };
  /**
   * A recording of a chunk (`part`) of a pass: its draw list and what it drew (`snapshot`, `valid` once recorded);
   * for a chunk of objects, the rows its walk read (`rows`: per position of its range the slot, its row and values
   * versions; `invalidations` of the table): unchanged, it holds while other objects change (`isChunkCurrent`)
   */
  type Chunk = {
    readonly list: DrawList;
    readonly snapshot: Snapshot;
    readonly part: PassPart;
    rows: Uint32Array;
    rowCount: number;
    invalidations: number;
    /** the visibility of its rows per position (`getRowVisibility`) when drawn */
    readonly visible: Uint8Array;
    /** whether a visible row's draws follow the view besides its visibility (`isViewDependent`) */
    viewDependent: boolean;
    /** the submission its inputs last changed in (`drawChunk`: drawn directly while changing) */
    changedAt: number;
    /** how it was drawn last (`drawChunk`) */
    outcome: ChunkOutcome;
  };
  type ChunkOutcome = 'unplanned' | 'replayed' | 'flushed' | 'empty' | 'recorded' | 'direct';
  /**
   * The recordings of a pass in chunks of `ChunkSize` objects of the table's order, then one of the caps drawn after
   * all objects: a change of what the pass draws (an object marked, hidden, entering the view) records the chunks it
   * touches again, the others replay (a hover sweeping across many objects, a zoomed view rotating)
   */
  type PassDrawList = {
    readonly chunks: Chunk[];
    caps: Chunk | undefined;
    /** what the walks read when every chunk was last drawn: unchanged, all replay without a check of their own */
    readonly inputs: Snapshot;
  };
  const ChunkSize = 64;
  /**
   * Submissions a chunk's inputs stay unchanged before it is recorded again after changing in consecutive ones: drawn
   * directly meanwhile, its data or view changing every frame (a hover, a recolor, a rotation) would else be walked,
   * flushed or recorded each frame, which cost more than drawing it (merged hover 19-21%, unmerged recolor 17-21%)
   */
  const StableSubmissions = 3;
  /**
   * Rows of a chunk drawn directly as cheaply as walked and flushed (a merged object is one row): more rows changing
   * every frame (a hover sweeping across separate objects, a rotation) are flushed and replayed, cheaper than drawn
   */
  const FewChunkRows = 4;
  /** The caps drawn after all objects, a chunk of their own */
  const CapsPart: PassPart = { begin: 0, end: Infinity, objects: false, caps: true };
  let groupDrawLists = new WeakMap<Scene.Group, Map<PassDraws, PassDrawLists>>();
  const resetDrawLists = () => {
    groupDrawLists = new WeakMap();
  };
  const getDrawList = (pass: PassDraws, group: Scene.Group) => {
    if (!ctx.drawLists) return;
    let passLists = groupDrawLists.get(group);
    if (!passLists) groupDrawLists.set(group, (passLists = new Map()));
    let lists = passLists.get(pass);
    if (!lists) passLists.set(pass, (lists = { submit: -1, next: 0, lists: [] }));
    const submit = ctx.drawLists?.submitCount ?? 0;
    if (lists.submit !== submit) {
      lists.submit = submit;
      lists.next = 0;
    }
    const i = lists.next++;
    return (lists.lists[i] ??= { chunks: [], caps: undefined, inputs: createSnapshot() });
  };
  const createChunk = (part: PassPart): Chunk => ({
    list: ctx.drawLists!.create(),
    snapshot: createSnapshot(),
    part,
    rows: new Uint32Array(part.objects ? 3 * ChunkSize : 0),
    rowCount: -1,
    invalidations: -1,
    visible: new Uint8Array(part.objects ? ChunkSize : 0),
    viewDependent: true,
    changedAt: -Infinity,
    outcome: 'empty',
  });
  /**
   * Whether the draws of a row (`facts`) follow the view besides whether it is visible: instances culled
   * (`cullObject`), levels of detail, solid interior caps (cut by the near plane)
   */
  /**
   * Whether an object's draws in a pass change with the view: its levels and culled instances; its solid interior's
   * caps only where the pass draws them with the object (`inlineCaps`), else the caps chunk compares all it read
   */
  const isViewDependent = (facts: number, inlineCaps: boolean) =>
    (facts & DrawFact.LodLevels) !== 0 ||
    (inlineCaps && (facts & DrawFact.SolidInterior) !== 0) ||
    (facts & (DrawFact.InstanceGrid | DrawFact.MultipleInstances)) ===
      (DrawFact.InstanceGrid | DrawFact.MultipleInstances);
  /**
   * Whether nothing the walk of `chunk` read changed: for a chunk of objects its rows (slots in order, row and values
   * versions) and the table's invalidations, so that a change of other objects leaves it current; the caps chunk
   * reads all objects (`isSnapshotCurrent`). Assumes the table refreshed
   */
  const isChunkCurrent = (chunk: Chunk, table: DrawTable) =>
    chunk.part.objects
      ? chunk.snapshot.valid && isChunkInputsCurrent(chunk, table)
      : isSnapshotCurrent(chunk.snapshot, table);
  /** `isChunkCurrent` of a chunk of objects but its recording: whether what it would draw is as when last drawn */
  const isChunkInputsCurrent = (chunk: Chunk, table: DrawTable) => {
    const { snapshot, part, rows } = chunk;
    if (chunk.rowCount < 0 || host.occlusionTest) return false;
    const end = Math.min(part.end, table.count);
    if (snapshot.tableVersion !== table.version || snapshot.valuesEpoch !== table.valuesEpoch) {
      if (chunk.invalidations !== table.invalidations || chunk.rowCount !== end - part.begin) return false;
      const { order, rowVersions, valuesVersions } = table;
      for (let k = part.begin, i = 0; k < end; ++k, i += 3) {
        const slot = order[k];
        if (rows[i] !== slot || rows[i + 1] !== rowVersions[slot] || rows[i + 2] !== valuesVersions[slot]) {
          return false;
        }
      }
    }
    if (isSameInputs(snapshot)) return true;
    // a recording binding a replaced resource is not current, whatever the view did
    if (snapshot.drawResources !== (ctx.drawLists?.resourcesVersion ?? 0)) return false;
    // the view changed (a camera rotating): the same draws where its rows are as visible as they were and none of them
    // draws by the view otherwise
    if (chunk.viewDependent || snapshot.modelScale !== host.modelScale || snapshot.caps !== isCapsDrawn()) return false;
    const visible = getRowVisibility(table);
    const { order } = table;
    for (let k = part.begin, i = 0; k < end; ++k, ++i) {
      if (visible[order[k]] !== chunk.visible[i]) return false;
    }
    return true;
  };
  /** `chunk` drawn as the table is now (`isChunkCurrent`) */
  const setChunkInputs = (chunk: Chunk, table: DrawTable, inlineCaps: boolean) => {
    const { snapshot, part, rows } = chunk;
    setSnapshotInputs(snapshot, table);
    if (!part.objects) return;
    const end = Math.min(part.end, table.count);
    const { order, rowVersions, valuesVersions } = table;
    for (let k = part.begin, i = 0; k < end; ++k, i += 3) {
      const slot = order[k];
      rows[i] = slot;
      rows[i + 1] = rowVersions[slot];
      rows[i + 2] = valuesVersions[slot];
    }
    chunk.rowCount = Math.max(0, end - part.begin);
    chunk.invalidations = table.invalidations;
    const visible = getRowVisibility(table);
    const { bits } = table;
    let viewDependent = false;
    for (let k = part.begin, i = 0; k < end; ++k, ++i) {
      const slot = order[k];
      chunk.visible[i] = visible[slot];
      if (visible[slot] && isViewDependent(bits[slot * DrawRow.Stride + DrawRow.Facts], inlineCaps)) {
        viewDependent = true;
      }
    }
    chunk.viewDependent = viewDependent;
  };
  const isEmptySnapshot = (snapshot: Snapshot) => snapshot.count === 0 && snapshot.capKey.length === 0;

  /**
   * Draws `part` of `pass` from `chunk`: 'unplanned' replayed without walking (nothing the walk read changed),
   * 'replayed' or 'flushed' (its data written first) after walking, 'empty' drawing nothing, 'recorded', 'direct'
   * drawn without its recording while its inputs change in consecutive submissions (`StableSubmissions`)
   */
  const drawChunk = (pass: PassDraws, table: DrawTable, chunk: Chunk): ChunkOutcome => {
    const outcome = drawChunkOnce(pass, table, chunk);
    chunk.outcome = outcome;
    return outcome;
  };
  /**
   * Whether `chunk`, changed again in the submission after its last change, is drawn directly: where its recording
   * was made again (draws changing every frame, a recolor), or it was drawn directly, or it has few rows
   */
  const isDrawnDirectly = (chunk: Chunk) =>
    chunk.outcome === 'recorded' || chunk.outcome === 'direct' || chunk.rowCount <= FewChunkRows;
  /**
   * Whether `part` holds an object culled on the GPU while it moves: the records its lists hold change from frame to
   * frame (a recording per frame would cost more than it saves), and multi-draw indirect, where the device has it, is
   * not in bundles; so drawn directly. At rest they hold those with instances (`recordsCompact`), recorded as other
   * draws
   */
  const isMultiDrawn = (table: DrawTable, part: PassPart) => {
    if (host.gpuCulledCount === 0) return false;
    const { order, objects } = table;
    for (let k = part.begin, kl = Math.min(part.end, table.count); k < kl; ++k) {
      if (objects[order[k]]?.recordsCompact === false) return true;
    }
    return false;
  };
  const drawChunkOnce = (pass: PassDraws, table: DrawTable, chunk: Chunk): ChunkOutcome => {
    const { snapshot, list, part } = chunk;
    const inlineCaps = !!pass.caps?.inline;
    if (part.objects && isMultiDrawn(table, part)) {
      renderPassRows(pass, table, part);
      if (snapshot.valid) list.reset();
      snapshot.valid = false;
      setChunkInputs(chunk, table, inlineCaps);
      return 'direct';
    }
    if (part.objects) {
      const submit = ctx.drawLists?.submitCount ?? 0;
      const current = isChunkInputsCurrent(chunk, table);
      if (current && snapshot.valid) {
        if (isEmptySnapshot(snapshot)) return 'empty';
        if (list.replay(globals)) return 'unplanned';
      } else if (
        current
          ? submit - chunk.changedAt < StableSubmissions && chunk.outcome === 'direct'
          : submit - chunk.changedAt <= 1 && isDrawnDirectly(chunk)
      ) {
        // changing every frame, or not stable for long enough since: drawn as without draw lists
        if (!current) chunk.changedAt = submit;
        renderPassRows(pass, table, part);
        // its recording no longer of the draws, released once
        if (snapshot.valid) list.reset();
        snapshot.valid = false;
        setChunkInputs(chunk, table, inlineCaps);
        return 'direct';
      } else if (!current) {
        chunk.changedAt = submit;
      }
    } else if (isChunkCurrent(chunk, table)) {
      if (isEmptySnapshot(snapshot)) return 'empty';
      if (list.replay(globals)) return 'unplanned';
    }
    // a recording binding a replaced resource (`DrawLists.resourcesVersion`) does not replay, whatever the
    // walk finds
    const bound = snapshot.drawResources === (ctx.drawLists?.resourcesVersion ?? 0);
    const match = snapshot.valid && bound ? walkPass(pass, table, snapshot, false, part) : 'different';
    if (match === 'same' || (match === 'data' && flushDataChanges(snapshot, table, pass.variant))) {
      if (isEmptySnapshot(snapshot)) {
        setChunkInputs(chunk, table, inlineCaps);
        return 'empty';
      }
      if (list.replay(globals)) {
        setChunkInputs(chunk, table, inlineCaps);
        return match === 'data' ? 'flushed' : 'replayed';
      }
    } else {
      // the draws changed: so are the recordings for other targets
      list.reset();
    }
    // as drawn directly, the caps (and their clears) included
    list.record(globals, () => renderPassRows(pass, table, part));
    // drawing may change values (e.g. what culling writes): the snapshot after it is the recorded one
    walkPass(pass, table, snapshot, true, part);
    if (ctx.drawLists) {
      for (let i = 0, il = snapshot.count; i < il; ++i) snapshot.digests[i] = captureDraw(snapshot, i, pass.variant);
    }
    setChunkInputs(chunk, table, inlineCaps);
    snapshot.valid = true;
    // a chunk drawing nothing keeps no recording
    if (isEmptySnapshot(snapshot)) list.reset();
    return 'recorded';
  };

  /**
   * The chunks of `drawList` replayed as recorded (those drawing nothing skipped), the caps chunk last (index
   * `chunks.length`): the index of the first that can't replay (no recording for this target, globals bound
   * otherwise), from which the chunks are drawn one by one; -1 where all replayed
   */
  const replayChunks = (drawList: PassDrawList) => {
    const { chunks, caps } = drawList;
    for (let c = 0, il = chunks.length; c < il; ++c) {
      const { snapshot, list } = chunks[c];
      // drawn directly last (`drawChunk`): no recording, its snapshot not of the draws
      if (!snapshot.valid) return c;
      if (!isEmptySnapshot(snapshot) && !list.replay(globals)) return c;
    }
    if (caps && !isEmptySnapshot(caps.snapshot) && !caps.list.replay(globals)) return chunks.length;
    return -1;
  };

  /**
   * Draws `pass` from a draw list (`drawLists`) chunk by chunk (`PassDrawList`): replays each chunk's recording while
   * it draws the recorded draws with unchanged values, caps and globals, records it otherwise (as `renderPass` draws
   * it); false without draw lists
   */
  const renderPassFromDrawList = (pass: PassDraws, group: Scene.Group) => {
    const { table } = group;
    if (!table) return false;
    host.setDrawingFromList(true);
    try {
      return drawPassFromList(pass, group, table);
    } finally {
      host.setDrawingFromList(false);
    }
  };
  const drawPassFromList = (pass: PassDraws, group: Scene.Group, table: DrawTable) => {
    const drawList = getDrawList(pass, group);
    if (!drawList) return false;
    const { chunks, inputs } = drawList;
    // volatile cells compared, once for all chunks
    table.refresh();
    // nothing any chunk's walk read changed: each replays its recording, from the first that can't one by one
    let first = 0;
    if (isSnapshotCurrent(inputs, table)) {
      first = replayChunks(drawList);
      if (first < 0) {
        ctx.stats.drawLists.replayed += 1;
        ctx.stats.drawLists.unplanned += 1;
        return true;
      }
    }
    inputs.valid = false;
    const chunkCount = Math.ceil(table.count / ChunkSize);
    while (chunks.length < chunkCount) {
      const begin = chunks.length * ChunkSize;
      chunks.push(createChunk({ begin, end: begin + ChunkSize, objects: true, caps: false }));
    }
    // objects removed: their chunks' recordings go
    if (chunks.length > chunkCount) chunks.length = chunkCount;
    let recorded = 0;
    let direct = 0;
    let unplanned = true;
    let flushed = false;
    const count = (result: ReturnType<typeof drawChunk>) => {
      if (result === 'recorded') recorded += 1;
      if (result === 'direct') direct += 1;
      if (result !== 'unplanned' && result !== 'empty') unplanned = false;
      if (result === 'flushed') flushed = true;
    };
    for (let c = first; c < chunkCount; ++c) count(drawChunk(pass, table, chunks[c]));
    if (pass.caps && !pass.caps.inline) {
      drawList.caps ??= createChunk(CapsPart);
      count(drawChunk(pass, table, drawList.caps));
    }
    setSnapshotInputs(inputs, table);
    inputs.valid = true;
    const { drawLists } = ctx.stats;
    drawLists.chunksDirect += direct;
    if (recorded > 0) {
      drawLists.recorded += 1;
      drawLists.chunksRecorded += recorded;
    } else {
      drawLists.replayed += 1;
      if (unplanned) drawLists.unplanned += 1;
      if (flushed) drawLists.flushed += 1;
    }
    return true;
  };

  return { render: renderPassFromDrawList, reset: resetDrawLists };
}
