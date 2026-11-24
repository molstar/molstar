/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval, OrderedSet } from '@molstar/core/data/int';
import {
  CopyLayerBlock,
  CopyLayerKinds,
  CopyLayerValue,
  writeCopyLayerImage,
} from '@molstar/graphics/gl/instance-expansion';
import type { TextureImage } from '@molstar/graphics/gl/renderable/util';
import { applyClippingGroups } from '@molstar/graphics/geo/geometry/clipping-data';
import { applyEmissiveValue } from '@molstar/graphics/geo/geometry/emissive-data';
import { applyOverpaintColor, clearOverpaint } from '@molstar/graphics/geo/geometry/overpaint-data';
import { applySubstanceMaterial, clearSubstance } from '@molstar/graphics/geo/geometry/substance-data';
import { applyTransparencyValue } from '@molstar/graphics/geo/geometry/transparency-data';
import { applyWiggleValue } from '@molstar/graphics/geo/geometry/wiggle-data';
import { type Loci as ModelLoci, isEveryLoci } from '@molstar/model/model/loci';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';
import { DefaultCopyBlockBudget } from './copy-markers.js';
import type { TemplateCopies } from './expanded.js';

/** The layers a copy has its own of, in the order of a copy's words (`CopyLayerKinds`) */
export const CopyLayerKeys = ['overpaint', 'transparency', 'clipping', 'emissive', 'substance', 'wiggle'] as const;
export type CopyLayerKey = (typeof CopyLayerKeys)[number];

/** a layer of any of them: its loci and what it writes */
export type AnyLayer = { readonly loci: ModelLoci } & Record<string, any>;
export type AnyLayers = { readonly layers: ReadonlyArray<AnyLayer> };

/** a kind's bytes an item, and a layer's write of items `start` to `end` (as `Visual.setOverpaint` and the others) */
const Kinds: Record<
  CopyLayerKey,
  { size: number; apply(array: Uint8Array, start: number, end: number, layer: AnyLayer): void }
> = {
  overpaint: {
    size: 4,
    apply: (a, s, e, l) => (l.clear ? clearOverpaint(a, s, e) : applyOverpaintColor(a, s, e, l.color)),
  },
  transparency: { size: 1, apply: (a, s, e, l) => applyTransparencyValue(a, s, e, l.value) },
  clipping: { size: 1, apply: (a, s, e, l) => applyClippingGroups(a, s, e, l.groups) },
  emissive: { size: 1, apply: (a, s, e, l) => applyEmissiveValue(a, s, e, l.value) },
  substance: {
    size: 4,
    apply: (a, s, e, l) => (l.clear ? clearSubstance(a, s, e) : applySubstanceMaterial(a, s, e, l.material)),
  },
  wiggle: { size: 1, apply: (a, s, e, l) => applyWiggleValue(a, s, e, l.value) },
};

/**
 * a copy's own layer: its value (the same for all its template particles), or a block by template particle; `coarse`
 * where a block did not fit the budget (`CopyLayerRoom`): its layers applied to its value whole
 */
export type CopyState = { value?: Uint8Array; block?: Uint8Array; coarse?: boolean };

/**
 * The blocks' bytes left of a budget (`getCopyBlockBudget`), shared by the kinds and templates of a layout's layers,
 * and the copies past it (taking a layer whole, `CopyState.coarse`). Why: a layer differing by copy makes a block each
 * (a template's particles times the kind's bytes: the isg cell's 1.8 MB of overpaint a copy), past the device's limit
 * the buffer not made
 */
export interface CopyLayerRoom {
  bytes: number;
  coarse: number;
}

/** `room` past its budget: warned of once (a layout's), whether now */
export function warnCopyLayerRoom(room: CopyLayerRoom, budget: number) {
  if (room.coarse === 0) return false;
  console.warn(
    `particle copies: their own layers past ${Math.round(budget / 2 ** 20)} MiB, ${room.coarse} copies layered in part layered whole`,
  );
  return true;
}

/** A kind's bytes an item */
export function getCopyLayerSize(key: CopyLayerKey) {
  return Kinds[key].size;
}

/** `layer`'s write of a kind's items `start` to `end` of `array` (as `Visual.setOverpaint` and the others write them) */
export function writeLayerItems(key: CopyLayerKey, array: Uint8Array, start: number, end: number, layer: AnyLayer) {
  Kinds[key].apply(array, start, end, layer);
}

/** `set`'s runs (`OrderedSet`) as half-open intervals */
function forEachRun(set: OrderedSet<number>, f: (start: number, end: number) => void) {
  if (OrderedSet.isInterval(set)) {
    if (Interval.size(set) > 0) f(Interval.start(set), Interval.end(set));
    return;
  }
  let start = -1;
  let end = -1;
  OrderedSet.forEach(set, (i) => {
    if (i === end) {
      end += 1;
      return;
    }
    if (start >= 0) f(start, end);
    start = i;
    end = i + 1;
  });
  if (start >= 0) f(start, end);
}

/**
 * The own layers of `template`'s copies `tops` (of `list`) of a kind by `layers` in order: a copy has one once a layer
 * reaches it by copy (whole: a value; in part: a block, from its value or from the template's layers so far); a layer
 * on the same particles of every copy the template's (its objects' own, `getTemplateWideLocals`), and those copies'
 * that have one. The copies without one their objects' (none here)
 */
export function getCopyLayerStates(
  list: ParticleList,
  template: ParticleList,
  tops: ArrayLike<number>,
  key: CopyLayerKey,
  layers: AnyLayers,
  room?: CopyLayerRoom,
  blockBytes?: number,
): Map<number, CopyState> {
  const { size, apply } = Kinds[key];
  const n = template.count;
  /** a block's bytes as the GPU holds it (the paged store's by its pages' span, its levels' too) */
  const cost = blockBytes ?? n * size;
  const all = Interval.ofBounds(0, n);
  const states = new Map<number, CopyState>();
  /** the template's layers so far: a block made later starts from them */
  const wide: [OrderedSet<number>, AnyLayer][] = [];
  const applyTo = (s: CopyState, locals: OrderedSet<number>, layer: AnyLayer) => {
    if (OrderedSet.size(locals) === n) {
      s.value = new Uint8Array(size);
      apply(s.value, 0, 1, layer);
      s.block = undefined;
      return;
    }
    if (s.coarse) {
      apply(s.value!, 0, 1, layer);
      return;
    }
    if (!s.block && room) {
      if (room.bytes < cost) {
        // no block within the budget: the copy's value, this layer and the ones after whole
        s.coarse = true;
        room.coarse += 1;
        if (!s.value) s.value = new Uint8Array(size);
        apply(s.value, 0, 1, layer);
        return;
      }
      room.bytes -= cost;
    }
    if (!s.block) {
      const block = new Uint8Array(n * size);
      if (s.value) for (let i = 0; i < n; ++i) block.set(s.value, i * size);
      else for (const [w, l] of wide) forEachRun(w, (a, b) => apply(block, a, b, l));
      s.block = block;
      s.value = undefined;
    }
    const block = s.block;
    forEachRun(locals, (a, b) => apply(block, a, b, layer));
  };
  const topIndex = new Map<number, number>();
  for (let c = 0; c < tops.length; ++c) topIndex.set(tops[c], c);
  for (const layer of layers.layers) {
    const loci = layer.loci;
    let common: OrderedSet<number> | undefined;
    const reached: [number, OrderedSet<number>][] = [];
    if (isEveryLoci(loci)) {
      common = all;
    } else if (Particle.isLoci(loci) && loci.particles === list) {
      // per copy what it holds: whole by its top particle, in part by the loci's partial
      const held = new Map<number, OrderedSet<number>>();
      OrderedSet.forEach(loci.indices, (top) => {
        const c = topIndex.get(top);
        if (c !== undefined) held.set(c, all);
      });
      const partial = loci.partial;
      if (partial) {
        for (let i = 0; i < partial.tops.length; ++i) {
          const c = topIndex.get(partial.tops[i]);
          if (c !== undefined && !held.has(c)) held.set(c, partial.locals[i]);
        }
      }
      if (held.size === 0) continue;
      if (held.size === tops.length) {
        for (const locals of held.values()) {
          common = common ? OrderedSet.intersect(common, locals) : locals;
          if (OrderedSet.isEmpty(common)) break;
        }
        if (common && OrderedSet.isEmpty(common)) common = undefined;
      }
      for (const [c, locals] of held) {
        const rest = common ? OrderedSet.subtract(locals, common) : locals;
        if (!OrderedSet.isEmpty(rest)) reached.push([c, rest]);
      }
    } else {
      continue;
    }
    if (common) {
      wide.push([common, layer]);
      for (const s of states.values()) applyTo(s, common, layer);
    }
    for (const [c, locals] of reached) {
      let s = states.get(c);
      if (!s) states.set(c, (s = {}));
      applyTo(s, locals, layer);
    }
  }
  return states;
}

/**
 * A list with templates' copies' own layers (its templates' copies one after another, as `TemplateCopyMarkers`):
 * `CopyLayerKinds` words a copy (none, `CopyLayerValue` or `CopyLayerBlock` and its first byte), then their bytes;
 * what its objects' vertex stage reads by an instance's copy and template particle (`tCopyLayers`)
 */
export interface TemplateCopyLayers {
  readonly image: TextureImage<Uint8Array>;
  /** whether a copy of `template` has a layer of its own (of `key`, or of any) */
  has(template: ParticleList, key?: CopyLayerKey): boolean;
  /** the copies of `templates` (in order), no layer */
  layout(list: ParticleList, templates: readonly (readonly [ParticleList, TemplateCopies])[]): void;
  /** the copies' layers by `layers` (each given its kind's again), the image written; whether a copy's changed */
  set(list: ParticleList, layers: Partial<Record<CopyLayerKey, AnyLayers>>): boolean;
}

/** `getBudget` the blocks' bytes at most (`getCopyBlockBudget`), read at each `set` (a representation's graphics set later) */
export function createTemplateCopyLayers(getBudget: () => number = () => DefaultCopyBlockBudget): TemplateCopyLayers {
  /** warned of copies layered whole since the layout */
  let warned = false;
  let templates: (readonly [ParticleList, TemplateCopies])[] = [];
  /** per kind per template its copies' states (by copy) */
  const states = new Map<CopyLayerKey, Map<ParticleList, Map<number, CopyState>>>();
  let image = writeCopyLayerImage(undefined, new Uint8Array(4));
  const write = () => {
    let copyCount = 0;
    for (const [, copies] of templates) copyCount += copies.tops.length;
    const wordBytes = copyCount * CopyLayerKinds * 4;
    let byteLength = wordBytes;
    for (const byTemplate of states.values()) {
      for (const byCopy of byTemplate.values()) {
        for (const s of byCopy.values()) byteLength += ((s.value ?? s.block)!.length + 3) & ~3;
      }
    }
    const buffer = new Uint8Array(Math.max(4, byteLength));
    const words = new Uint32Array(buffer.buffer, 0, copyCount * CopyLayerKinds);
    let at = wordBytes;
    let base = 0;
    for (const [template, copies] of templates) {
      CopyLayerKeys.forEach((key, k) => {
        const byCopy = states.get(key)?.get(template);
        if (!byCopy) return;
        for (const [c, s] of byCopy) {
          const bytes = (s.value ?? s.block)!;
          buffer.set(bytes, at);
          words[(base + c) * CopyLayerKinds + k] = ((s.value ? CopyLayerValue : CopyLayerBlock) | at) >>> 0;
          at += (bytes.length + 3) & ~3;
        }
      });
      base += copies.tops.length;
    }
    image = writeCopyLayerImage(image, buffer);
  };
  return {
    get image() {
      return image;
    },
    has(template, key) {
      if (key) return (states.get(key)?.get(template)?.size ?? 0) > 0;
      for (const k of CopyLayerKeys) if ((states.get(k)?.get(template)?.size ?? 0) > 0) return true;
      return false;
    },
    layout(_list, next) {
      templates = [...next];
      states.clear();
      warned = false;
      write();
    },
    set(list, layers) {
      let changed = false;
      const budget = getBudget();
      const room: CopyLayerRoom = { bytes: budget, coarse: 0 };
      // the kinds kept (not given) hold their blocks
      for (const [key, byTemplate] of states) {
        if (layers[key]) continue;
        for (const byCopy of byTemplate.values()) for (const s of byCopy.values()) room.bytes -= s.block?.length ?? 0;
      }
      for (const key of CopyLayerKeys) {
        const value = layers[key];
        if (!value) continue;
        const had = states.get(key);
        const byTemplate = new Map<ParticleList, Map<number, CopyState>>();
        for (const [template, copies] of templates) {
          const byCopy = getCopyLayerStates(list, template, copies.tops, key, value, room);
          if (byCopy.size > 0) byTemplate.set(template, byCopy);
        }
        if ((had?.size ?? 0) > 0 || byTemplate.size > 0) changed = true;
        states.set(key, byTemplate);
      }
      if (!warned) warned = warnCopyLayerRoom(room, budget);
      if (changed) write();
      return changed;
    },
  };
}
