/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { type GraphicsRenderable, RenderableState } from './renderable.js';

/**
 * What the renderer reads from an object's values, as bits of a row (`DrawRow.Facts`); `getDrawFacts` derives them,
 * from the cells in `DrawFactCells`: what its pass checks read (`getCheckFacts`) and what culling reads
 */
export enum DrawFact {
  /** `drawCount` > 0 */
  DrawCount = 1 << 0,
  /** `uAlpha` === 1 */
  AlphaOne = 1 << 1,
  /** `uAlpha` === 0 */
  AlphaZero = 1 << 2,
  /** `transparencyAverage` > 0 */
  TransparencyAny = 1 << 3,
  /** `transparencyAverage` === 1 */
  TransparencyFull = 1 << 4,
  /** `markerAverage` > 0 */
  MarkerAny = 1 << 5,
  /** `markerAverage` === 1 */
  MarkerFull = 1 << 6,
  /** `emissiveAverage` + `uEmissive` > 0 */
  Emissive = 1 << 7,
  /** `dGeometryType` 'directVolume' */
  DirectVolume = 1 << 8,
  /** `dGeometryType` 'text' or 'image' */
  TextOrImage = 1 << 9,
  /** `dPointStyle` 'fuzzy' */
  FuzzyPoints = 1 << 10,
  /** `dXrayShaded` 'on' or 'inverted' */
  Xray = 1 << 11,
  /** `dTransparentBackfaces` 'opaque' */
  BackfacesOpaque = 1 << 12,
  /** `dTransparentBackfaces` 'off' */
  BackfacesOff = 1 << 13,
  /** `uDoubleSided` */
  DoubleSided = 1 << 14,
  /** `dGeometryType` 'mesh' or 'textureMesh' with `dSolidInterior` */
  SolidInterior = 1 << 15,
  /** `instanceGrid` with cells (`cellSize` > 0) */
  InstanceGrid = 1 << 16,
  /** `uInstanceCount` > 1 */
  MultipleInstances = 1 << 17,
  /** `lodLevels` */
  LodLevels = 1 << 18,
  /** `uAlpha` < 1 */
  AlphaBelowOne = 1 << 19,
  /** `dFlipSided` */
  FlipSided = 1 << 20,
  /** `hasReflection` (reflected instances flip the winding) */
  HasReflection = 1 << 21,
  /** `uDoubleSided` given and false (without it: neither single nor double sided, nothing culled) */
  SingleSided = 1 << 22,
  /** wiggling (`uWiggleAmplitude` or `wiggleAverage`, with speed and frequency) or tumbling: drawn anew per frame */
  Animated = 1 << 23,
}

/** The cells `getDrawFacts` and the row's bounds read: a change of any refreshes the row */
export const DrawFactCells = [
  'drawCount',
  'uAlpha',
  'transparencyAverage',
  'markerAverage',
  'emissiveAverage',
  'uEmissive',
  'dGeometryType',
  'dPointStyle',
  'dXrayShaded',
  'dTransparentBackfaces',
  'uDoubleSided',
  'dFlipSided',
  'hasReflection',
  'dSolidInterior',
  'instanceGrid',
  'uInstanceCount',
  'lodLevels',
  'boundingSphere',
  'uLod',
  // the scene's averages (`DrawTable.contributions`)
  'wiggleAverage',
  'transparencyMin',
  'uWiggleAmplitude',
  'uWiggleSpeed',
  'uWiggleFrequency',
  'uTumbleAmplitude',
  'uTumbleSpeed',
  'uTumbleFrequency',
] as const;

type Values = GraphicsRenderable['values'];

function getValue(values: Values, name: (typeof DrawFactCells)[number]): unknown {
  return (values as unknown as { [k: string]: ValueCell<unknown> | undefined })[name]?.ref.value;
}

/** The `DrawFact` bits of `values` */
export function getDrawFacts(values: Values): number {
  let f = getCheckFacts(values) | getCullFacts(values);
  if (isAnimated(values)) f |= DrawFact.Animated;
  return f;
}

function isAnimated(values: Values) {
  const n = (name: string) => (getValue(values, name as (typeof DrawFactCells)[number]) as number | undefined) ?? 0;
  const wiggle =
    (n('uWiggleAmplitude') > 0 || n('wiggleAverage') > 0) && n('uWiggleSpeed') > 0 && n('uWiggleFrequency') > 0;
  const tumble = n('uTumbleAmplitude') > 0 && n('uTumbleSpeed') > 0 && n('uTumbleFrequency') > 0;
  return wiggle || tumble;
}

/**
 * What an object adds to the scene's averages (`Scene.markerAverage`, ...), per slot (`DrawTable.contributions`):
 * marker, emission and wiggle averages, opacity, the smallest transparency of its transparent parts (1 where none),
 * whether its values make parts of it opaque (1; its state may as well, `DrawState.Opaque`)
 */
export type DrawContributions = {
  readonly marker: Float32Array;
  readonly emissive: Float32Array;
  readonly wiggle: Float32Array;
  readonly opacity: Float32Array;
  readonly transparencyMin: Float32Array;
  readonly opaque: Uint8Array;
};

function writeContributions(c: DrawContributions, slot: number, values: Values, facts: number) {
  c.marker[slot] = values.markerAverage.ref.value;
  c.emissive[slot] = values.emissiveAverage.ref.value + values.uEmissive.ref.value;
  c.wiggle[slot] = (getValue(values, 'wiggleAverage') as number | undefined) ?? 0;
  const alpha = values.uAlpha.ref.value;
  // x-ray shaded, fuzzy points, text and images are half transparent
  const half = (facts & (DrawFact.Xray | DrawFact.FuzzyPoints | DrawFact.TextOrImage)) !== 0;
  const xray = (facts & DrawFact.Xray) !== 0 ? 0.5 : 1;
  const fuzzy = (facts & DrawFact.FuzzyPoints) !== 0 ? 0.5 : 1;
  const textOrImage = (facts & DrawFact.TextOrImage) !== 0 ? 0.5 : 1;
  c.opacity[slot] = (1 - values.transparencyAverage.ref.value) * alpha * xray * fuzzy * textOrImage;
  let min = 1;
  if (alpha < 1) min = Math.min(min, 1 - alpha);
  if (half) min = Math.min(min, 0.5);
  const transparencyMin = (getValue(values, 'transparencyMin') as number | undefined) ?? 0;
  if (transparencyMin > 0) min = Math.min(min, transparencyMin);
  c.transparencyMin[slot] = min;
  c.opaque[slot] =
    (alpha === 1 && values.transparencyAverage.ref.value !== 1) || (facts & DrawFact.BackfacesOpaque) !== 0 ? 1 : 0;
}

/** The `DrawFact` bits of the draw count and of what culling reads: the instance grid and count, lod levels */
export function getCullFacts(values: Values): number {
  let f = 0;
  if (values.drawCount.ref.value) f |= DrawFact.DrawCount;
  if (values.instanceGrid.ref.value.cellSize > 0) f |= DrawFact.InstanceGrid;
  if (values.uInstanceCount.ref.value > 1) f |= DrawFact.MultipleInstances;
  if (getValue(values, 'lodLevels')) f |= DrawFact.LodLevels;
  return f;
}

/**
 * The `DrawFact` bits the renderer's pass checks and render states read: of alpha, transparency, markers, emission, the
 * geometry type, its style and faces; also of objects standing in for a kind (only those values, see the renderer's
 * pipeline requests)
 */
export function getCheckFacts(values: Values): number {
  let f = 0;
  const alpha = values.uAlpha.ref.value;
  if (alpha === 1) f |= DrawFact.AlphaOne;
  if (alpha === 0) f |= DrawFact.AlphaZero;
  if (alpha < 1) f |= DrawFact.AlphaBelowOne;
  const transparency = values.transparencyAverage.ref.value;
  if (transparency > 0) f |= DrawFact.TransparencyAny;
  if (transparency === 1) f |= DrawFact.TransparencyFull;
  const marker = values.markerAverage.ref.value;
  if (marker > 0) f |= DrawFact.MarkerAny;
  if (marker === 1) f |= DrawFact.MarkerFull;
  if (values.emissiveAverage.ref.value + values.uEmissive.ref.value > 0) f |= DrawFact.Emissive;
  const geometryType = values.dGeometryType.ref.value;
  if (geometryType === 'directVolume') f |= DrawFact.DirectVolume;
  if (geometryType === 'text' || geometryType === 'image') f |= DrawFact.TextOrImage;
  if (getValue(values, 'dPointStyle') === 'fuzzy') f |= DrawFact.FuzzyPoints;
  const xray = getValue(values, 'dXrayShaded');
  if (xray === 'on' || xray === 'inverted') f |= DrawFact.Xray;
  const backfaces = getValue(values, 'dTransparentBackfaces');
  if (backfaces === 'opaque') f |= DrawFact.BackfacesOpaque;
  if (backfaces === 'off') f |= DrawFact.BackfacesOff;
  const doubleSided = getValue(values, 'uDoubleSided');
  if (doubleSided) f |= DrawFact.DoubleSided;
  else if (doubleSided !== undefined) f |= DrawFact.SingleSided;
  if (getValue(values, 'dFlipSided')) f |= DrawFact.FlipSided;
  if (getValue(values, 'hasReflection')) f |= DrawFact.HasReflection;
  if ((geometryType === 'mesh' || geometryType === 'textureMesh') && getValue(values, 'dSolidInterior')) {
    f |= DrawFact.SolidInterior;
  }
  return f;
}

/** An object's state as bits of a row (`DrawRow.State`) */
export enum DrawState {
  Disposed = 1 << 0,
  Visible = 1 << 1,
  Pickable = 1 << 2,
  ColorOnly = 1 << 3,
  Opaque = 1 << 4,
  WriteDepth = 1 << 5,
}

export function getDrawState(state: RenderableState): number {
  return (
    (state.disposed ? DrawState.Disposed : 0) |
    (state.visible ? DrawState.Visible : 0) |
    (state.pickable ? DrawState.Pickable : 0) |
    (state.colorOnly ? DrawState.ColorOnly : 0) |
    (state.opaque ? DrawState.Opaque : 0) |
    (state.writeDepth ? DrawState.WriteDepth : 0)
  );
}

/**
 * The layout of a row, in 32-bit words (as a WGSL struct with std430 layout): the unscaled bounding sphere, the lod
 * distances (`uLod`, min and max), the facts and the state; the first six read as floats, the last two as bits
 */
export enum DrawRow {
  CenterX = 0,
  CenterY = 1,
  CenterZ = 2,
  Radius = 3,
  LodMin = 4,
  LodMax = 5,
  Facts = 6,
  State = 7,
  Stride = 8,
}

/**
 * The objects of a scene group as rows of flat data: what decides whether and how they draw, without reading their
 * values. Each object holds a slot (stable while it is in the group, reused after its removal); `order` lists the slots
 * in the group's order. Rows follow the objects' values and state: their changes queue the row, `refresh` writes the
 * queued rows (and reorders), advancing the row's version and the table's. Any change of an object's values advances
 * its values version and the values epoch (what draws read, not only the facts: data changes such as markers)
 */
export interface DrawTable {
  /** the rows by slot, as floats (`DrawRow`); replaced when the table grows */
  readonly floats: Float32Array;
  /** the rows by slot, as bits: the same buffer as `floats` */
  readonly bits: Uint32Array;
  /** the object of each slot, undefined for free slots */
  readonly objects: readonly (GraphicsRenderable | undefined)[];
  /** advances per slot when its row is written */
  readonly rowVersions: Uint32Array;
  /** the slots in the group's order, the first `count` */
  readonly order: Uint32Array;
  readonly count: number;
  /** slots in use and free: the length of the rows */
  readonly capacity: number;
  /** advances when a row is written or the order changes */
  readonly version: number;
  /** advances when a row's visibility (`DrawState.Visible`, not disposed) changes: objects shown, hidden, added, removed */
  readonly visibleVersion: number;
  /** advances with `invalidate` only: with the rows' and values' versions, what a part of the table derived from */
  readonly invalidations: number;
  /** what each slot's object adds to the scene's averages, written with its row (by slot; replaced when growing) */
  readonly contributions: DrawContributions;
  /**
   * advances per slot when any value of its object changes (`ValueCell.observe`), at once; volatile cells (changing
   * without notifying) are compared on `refresh`
   */
  readonly valuesVersions: Uint32Array;
  /** advances when any value of an object changes (`valuesVersions`), or an object is added or removed */
  readonly valuesEpoch: number;
  /**
   * the slots whose values changed (or that were added) since the last call, each once: the objects to update
   * (`Renderable.update`) after changes of a few of many (a hover); some may be free by now
   */
  takeChanged(): number[];
  /**
   * the slot's object is updated outside `takeChanged` (`Scene.update` of given objects): not taken until its values
   * change again
   */
  dropChanged(slot: number): void;
  /** the slot of `r`, -1 if not in the table */
  getSlot(r: GraphicsRenderable): number;
  add(r: GraphicsRenderable): void;
  remove(r: GraphicsRenderable): void;
  /** the group's objects were reordered: `order` follows on `refresh` */
  invalidateOrder(): void;
  /**
   * something of every object the rows don't hold changed (their programs, `Renderable.setTransparency`): advances
   * `version`, so that what was derived from the table is derived again
   */
  invalidate(): void;
  /** writes the queued rows and the order where they changed, compares volatile cells */
  refresh(): void;
  dispose(): void;
}

const InitialCapacity = 64;

/** Queues its slot when an observed cell or state changes (`ValueCell.Observer`) */
class SlotObserver implements ValueCell.Observer {
  constructor(
    readonly slot: number,
    private readonly queue: (slot: number) => void,
  ) {}
  get changed() {
    return false;
  }
  set changed(v: boolean) {
    if (v) this.queue(this.slot);
  }
}

/** Advances its slot's values version when an observed cell changes (`ValueCell.Observer`) */
class ValuesObserver implements ValueCell.Observer {
  constructor(
    readonly slot: number,
    private readonly changedValues: (slot: number) => void,
  ) {}
  get changed() {
    return false;
  }
  set changed(v: boolean) {
    if (v) this.changedValues(this.slot);
  }
}

/** What a slot observes: its row's cells and state, all its values, and the volatile ones among them */
type SlotObservers = {
  readonly row: SlotObserver;
  readonly values: ValuesObserver;
  readonly cells: ValueCell<unknown>[];
  readonly volatile: { readonly cell: ValueCell<unknown>; version: number }[];
};

/** The value cells of `values` (the entries holding a `ref`) */
function getValueCells(values: GraphicsRenderable['values']) {
  const cells: ValueCell<unknown>[] = [];
  for (const cell of Object.values(values) as (ValueCell<unknown> | undefined)[]) {
    if (cell && typeof cell === 'object' && 'ref' in cell) cells.push(cell);
  }
  return cells;
}

/** A table of the objects of `renderables` (a scene group's, kept in sync by its owner through `add` and `remove`) */
export function createDrawTable(renderables: ReadonlyArray<GraphicsRenderable>): DrawTable {
  let capacity = 0;
  let buffer = new ArrayBuffer(0);
  let floats = new Float32Array(buffer);
  let bits = new Uint32Array(buffer);
  let rowVersions = new Uint32Array(0);
  let valuesVersions = new Uint32Array(0);
  let contributions: DrawContributions = {
    marker: new Float32Array(0),
    emissive: new Float32Array(0),
    wiggle: new Float32Array(0),
    opacity: new Float32Array(0),
    transparencyMin: new Float32Array(0),
    opaque: new Uint8Array(0),
  };
  let queued = new Uint8Array(0);
  let order = new Uint32Array(0);
  const objects: (GraphicsRenderable | undefined)[] = [];
  const observers: (SlotObservers | undefined)[] = [];
  /** slots with volatile cells, compared when refreshing */
  const withVolatile = new Set<number>();
  const slots = new Map<GraphicsRenderable, number>();
  const free: number[] = [];
  const queue: number[] = [];
  /** slots whose state is not observable (`RenderableState.isObservable`): compared when refreshing */
  const polled: number[] = [];
  let orderDirty = false;
  let version = 0;
  let visibleVersion = 0;
  let invalidations = 0;
  let valuesEpoch = 0;

  const enqueue = (slot: number) => {
    if (queued[slot]) return;
    queued[slot] = 1;
    queue.push(slot);
  };

  /**
   * the slots of `takeChanged`, flagged in `inChanged`: 1 changed, 2 in `changed` but dropped (`dropChanged`, left in
   * the list so that dropping is O(1); skipped when taken), 0 not in the list
   */
  let changed: number[] = [];
  let inChanged = new Uint8Array(0);
  const changedValues = (slot: number) => {
    valuesVersions[slot] += 1;
    valuesEpoch += 1;
    if (inChanged[slot] === 1) return;
    if (inChanged[slot] === 0) changed.push(slot);
    inChanged[slot] = 1;
  };

  const grow = (min: number) => {
    const next = Math.max(min, capacity * 2, InitialCapacity);
    const nextBuffer = new ArrayBuffer(next * DrawRow.Stride * 4);
    new Uint32Array(nextBuffer).set(bits);
    buffer = nextBuffer;
    floats = new Float32Array(buffer);
    bits = new Uint32Array(buffer);
    const nextVersions = new Uint32Array(next);
    nextVersions.set(rowVersions);
    rowVersions = nextVersions;
    const nextValuesVersions = new Uint32Array(next);
    nextValuesVersions.set(valuesVersions);
    valuesVersions = nextValuesVersions;
    const grown = <T extends Float32Array | Uint8Array>(a: T): T => {
      const b = new (a.constructor as new (n: number) => T)(next);
      b.set(a);
      return b;
    };
    contributions = {
      marker: grown(contributions.marker),
      emissive: grown(contributions.emissive),
      wiggle: grown(contributions.wiggle),
      opacity: grown(contributions.opacity),
      transparencyMin: grown(contributions.transparencyMin),
      opaque: grown(contributions.opaque),
    };
    const nextQueued = new Uint8Array(next);
    nextQueued.set(queued);
    queued = nextQueued;
    const nextInChanged = new Uint8Array(next);
    nextInChanged.set(inChanged);
    inChanged = nextInChanged;
    capacity = next;
  };

  const isVisibleState = (state: number) => (state & (DrawState.Visible | DrawState.Disposed)) === DrawState.Visible;
  const write = (slot: number) => {
    const r = objects[slot];
    const o = slot * DrawRow.Stride;
    const wasVisible = isVisibleState(bits[o + DrawRow.State]);
    if (r) {
      const { values } = r;
      const { center, radius } = values.boundingSphere.ref.value;
      floats[o + DrawRow.CenterX] = center[0];
      floats[o + DrawRow.CenterY] = center[1];
      floats[o + DrawRow.CenterZ] = center[2];
      floats[o + DrawRow.Radius] = radius;
      const lod = values.uLod.ref.value;
      floats[o + DrawRow.LodMin] = lod[0];
      floats[o + DrawRow.LodMax] = lod[1];
      const facts = getDrawFacts(values);
      bits[o + DrawRow.Facts] = facts;
      bits[o + DrawRow.State] = getDrawState(r.state);
      writeContributions(contributions, slot, values, facts);
    } else {
      for (let k = 0; k < DrawRow.Stride; ++k) bits[o + k] = 0;
      bits[o + DrawRow.State] = DrawState.Disposed;
    }
    if (isVisibleState(bits[o + DrawRow.State]) !== wasVisible) visibleVersion += 1;
    rowVersions[slot] += 1;
    version += 1;
  };

  const observe = (r: GraphicsRenderable, slot: number): SlotObservers => {
    const row = new SlotObserver(slot, enqueue);
    const valuesObserver = new ValuesObserver(slot, changedValues);
    const values = r.values as unknown as { [k: string]: ValueCell<unknown> | undefined };
    for (const name of DrawFactCells) {
      const cell = values[name];
      if (cell) ValueCell.observe(cell, row);
    }
    const cells = getValueCells(r.values);
    const volatile: SlotObservers['volatile'] = [];
    for (const cell of cells) {
      ValueCell.observe(cell, valuesObserver);
      if (cell.volatile) volatile.push({ cell, version: cell.ref.version });
    }
    if (volatile.length > 0) withVolatile.add(slot);
    if (RenderableState.isObservable(r.state)) RenderableState.observe(r.state, row);
    else polled.push(slot);
    return { row, values: valuesObserver, cells, volatile };
  };

  const unobserve = (r: GraphicsRenderable, slot: number, o: SlotObservers) => {
    const values = r.values as unknown as { [k: string]: ValueCell<unknown> | undefined };
    for (const name of DrawFactCells) {
      const cell = values[name];
      if (cell) ValueCell.unobserve(cell, o.row);
    }
    for (const cell of o.cells) ValueCell.unobserve(cell, o.values);
    withVolatile.delete(slot);
    RenderableState.unobserve(r.state, o.row);
    const p = polled.indexOf(slot);
    if (p >= 0) polled.splice(p, 1);
  };

  const refreshOrder = () => {
    if (order.length < renderables.length) order = new Uint32Array(Math.max(renderables.length, order.length * 2));
    for (let i = 0, il = renderables.length; i < il; ++i) {
      const slot = slots.get(renderables[i]);
      if (slot === undefined) throw new Error('object of the group not in its draw table');
      order[i] = slot;
    }
    orderDirty = false;
    version += 1;
  };

  return {
    get floats() {
      return floats;
    },
    get bits() {
      return bits;
    },
    objects,
    get rowVersions() {
      return rowVersions;
    },
    get order() {
      return order;
    },
    get count() {
      return renderables.length;
    },
    get capacity() {
      return capacity;
    },
    get version() {
      return version;
    },
    get visibleVersion() {
      return visibleVersion;
    },
    get invalidations() {
      return invalidations;
    },
    get valuesVersions() {
      return valuesVersions;
    },
    get contributions() {
      return contributions;
    },
    get valuesEpoch() {
      return valuesEpoch;
    },
    takeChanged() {
      const taken = changed;
      changed = [];
      let n = 0;
      for (let i = 0, il = taken.length; i < il; ++i) {
        const slot = taken[i];
        if (inChanged[slot] === 1) taken[n++] = slot;
        inChanged[slot] = 0;
      }
      taken.length = n;
      return taken;
    },
    dropChanged(slot) {
      if (inChanged[slot] === 1) inChanged[slot] = 2;
    },
    getSlot(r) {
      return slots.get(r) ?? -1;
    },
    add(r) {
      if (slots.has(r)) return;
      let slot = free.pop();
      if (slot === undefined) {
        slot = objects.length;
        if (slot >= capacity) grow(slot + 1);
        objects.push(undefined);
        observers.push(undefined);
      }
      objects[slot] = r;
      slots.set(r, slot);
      observers[slot] = observe(r, slot);
      changedValues(slot);
      enqueue(slot);
      orderDirty = true;
    },
    remove(r) {
      const slot = slots.get(r);
      if (slot === undefined) return;
      unobserve(r, slot, observers[slot]!);
      slots.delete(r);
      objects[slot] = undefined;
      observers[slot] = undefined;
      free.push(slot);
      changedValues(slot);
      // its row cleared on refresh: a free slot reads as disposed
      enqueue(slot);
      orderDirty = true;
    },
    invalidateOrder() {
      orderDirty = true;
    },
    invalidate() {
      version += 1;
      invalidations += 1;
    },
    refresh() {
      if (withVolatile.size > 0) {
        for (const slot of withVolatile) {
          for (const v of observers[slot]!.volatile) {
            if (v.cell.ref.version === v.version) continue;
            v.version = v.cell.ref.version;
            changedValues(slot);
          }
        }
      }
      for (let i = 0, il = polled.length; i < il; ++i) {
        const slot = polled[i];
        if (bits[slot * DrawRow.Stride + DrawRow.State] !== getDrawState(objects[slot]!.state)) enqueue(slot);
      }
      for (let i = 0, il = queue.length; i < il; ++i) {
        const slot = queue[i];
        queued[slot] = 0;
        write(slot);
      }
      queue.length = 0;
      if (orderDirty) refreshOrder();
    },
    dispose() {
      for (let slot = 0, il = objects.length; slot < il; ++slot) {
        const r = objects[slot];
        if (r) unobserve(r, slot, observers[slot]!);
      }
      objects.length = 0;
      observers.length = 0;
      slots.clear();
      free.length = 0;
      queue.length = 0;
      polled.length = 0;
      withVolatile.clear();
    },
  };
}
