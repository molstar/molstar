/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval, OrderedSet } from '@molstar/core/data/int';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util/value-cell';
import {
  createExpansionCopies,
  createInstanceExpansion,
  createInstanceExpansionValues,
  type ExpansionCopies,
  getExpansionBoundingSphere,
  type InstanceExpansionValues,
  writeCopyMarkerImage,
} from '@molstar/graphics/gl/instance-expansion';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { TextureImage } from '@molstar/graphics/gl/renderable/util';
import { calculateTransformDataBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import type { TransformData } from '@molstar/graphics/geo/geometry/transform-data';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { type Loci as ModelLoci, isEveryLoci } from '@molstar/model/model/loci';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { createCopyMarkers } from './copy-markers.js';
import { StorePageSize } from './page-store.js';
import {
  Particle,
  type ParticleList,
  ParticleRotations,
  ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';

// A list with templates drawn where the vertex stage composes the instances (`instanceExpansion`, WebGL2): a template's
// leaf target an object of its particles in the template's frame, drawn once per copy the CPU's cull sees (its poses
// `instance-expansion.ts`); its own leaves drawn as a flat list's. Its copies' particles never expanded on the CPU

/** Whether `graphics` draws a list with templates by its copies' poses (`instanceExpansion`), not flattened */
export function canExpandTemplates(graphics: GraphicsContext | undefined) {
  return !!graphics?.capabilities.instanceExpansion;
}

/** The templates of `list` (its id space's order) */
export function getTemplates(list: ParticleList): ParticleList[] {
  return [...ParticleNesting.getIdSpace(list).moves.keys()];
}

/** A template's copies: their top particles (ascending), their poses and coarse levels, their bounds */
export interface TemplateCopies {
  readonly tops: readonly number[];
  readonly copies: ExpansionCopies;
  /** the corners of the template's box placed by each (`Particle.getBoundary`'s of a list with templates) */
  readonly boundingSphere: Sphere3D;
}

const _p = Vec3();
const _c = Vec3();
const _q = Quat();
const _m = Mat4();

/**
 * The copies of `template` in `list`, placed by `transform` (the representation's): a copy's pose
 * `transform T(t) R(q) T(-c)` (`c` its target's center); its template's coarse levels (`getTemplateLevels`) and its
 * particles' spacing (the cube root of its sphere's volume a particle, as the page store's)
 */
export function getTemplateCopies(list: ParticleList, template: ParticleList, transform: Mat4): TemplateCopies {
  // its particles drawn by a visual a target (WebGL): read, not by range (`Particle.DeferredParticles`)
  Particle.assertRead(template, "drawing a template's copies without the GPU's cull");
  const tops = ParticleNesting.getCopiesOf(list, template);
  const transforms = new Float32Array(tops.length * 16);
  const identity = Mat4.isIdentity(transform);
  for (let i = 0; i < tops.length; ++i) {
    const top = tops[i];
    const target = ParticleNesting.getTemplate(list, ParticleValues.get(list.targets, top))!;
    const center = ParticleTarget.getCenter(target) ?? Vec3.set(_c, 0, 0, 0);
    GridPositions.get(list.coordinates, top, _p);
    ParticleRotations.get(list.rotations, top, _q);
    // q (x - c) + t
    Mat4.fromQuat(_m, _q);
    Vec3.transformQuat(_c, Vec3.negate(_c, center), _q);
    Mat4.setTranslation(_m, Vec3.add(_c, _c, _p));
    if (!identity) Mat4.mul(_m, transform, _m);
    Mat4.toArray(_m, transforms, i * 16);
  }
  const { sphere, box } = Particle.getBoundary(template);
  const r = sphere.radius;
  const spacing = Math.cbrt(((4 / 3) * Math.PI * r * r * r) / Math.max(1, template.count));
  const levels = ParticleNesting.getTemplateLevels(list, template).length;
  const copies = createExpansionCopies(transforms, tops.length, sphere, spacing, levels);
  return { tops, copies, boundingSphere: getExpansionBoundingSphere(copies, box.min, box.max, 0) };
}

/**
 * `renderObject` (a template's leaf target at its coarse level `level`, 0 the template whole; its instances in the
 * template's frame) drawn by `copies`' poses: the expansion's values set (in place where it has them), its bounds the
 * copies' padded by its instances' reach
 */
export function setTargetExpansion(
  renderObject: GraphicsRenderObject,
  copies: TemplateCopies,
  level: number,
  /** the copies' markers' image (`writeCopyMarkerImage`) and the template's first copy in them */
  markers: { image: TextureImage<Uint8Array>; base: number },
  /** an instance's template particle by its instance id */
  toLocal: (instance: number) => number,
  /** the copies' own layers' image (`TemplateCopyLayers`, by the markers' copies) and whether a copy has one */
  layers: { image: TextureImage<Uint8Array>; has: boolean },
) {
  const values = renderObject.values as unknown as TransformData &
    Partial<InstanceExpansionValues> & {
      aInstance: ValueCell<ArrayLike<number>>;
      instanceCount: ValueCell<number>;
      invariantBoundingSphere: ValueCell<Sphere3D>;
      boundingSphere: ValueCell<Sphere3D>;
    };
  // the instances' sphere in the template's frame (a copy's culled by it placed)
  const sphere = calculateTransformDataBoundingSphere(values.invariantBoundingSphere.ref.value, values as never);
  const expansion = createInstanceExpansion(copies.copies, sphere, level);
  const n = values.instanceCount.ref.value;
  // per instance as stored (the grid's order) its template particle
  const ids = values.aInstance.ref.value;
  const locals = new Uint32Array(n);
  for (let i = 0; i < n; ++i) locals[i] = toLocal(ids[i]);
  if (values.instanceExpansion) {
    ValueCell.update(values.instanceExpansion, expansion);
    ValueCell.updateIfChanged(values.uExpansionCount!, n);
    ValueCell.updateIfChanged(values.tCopyPoses!, copies.copies.image);
    ValueCell.updateIfChanged(values.tCopyMarkers!, markers.image);
    ValueCell.updateIfChanged(values.uCopyMarkerBase!, markers.base);
    ValueCell.update(values.expansionLocals!, locals);
    setExpansionLayers(renderObject, layers);
  } else {
    Object.assign(
      values,
      createInstanceExpansionValues(expansion, n, markers.image, markers.base, locals, layers.image, layers.has),
    );
  }
  const reach = values.invariantBoundingSphere.ref.value.radius;
  ValueCell.update(values.boundingSphere, Sphere3D.expand(Sphere3D(), copies.boundingSphere, reach));
}

/** `renderObject`'s copies' own layers: their image (another where it grew), read where a copy has one */
export function setExpansionLayers(
  renderObject: GraphicsRenderObject,
  layers: { image: TextureImage<Uint8Array>; has: boolean },
) {
  const values = renderObject.values as unknown as Partial<InstanceExpansionValues>;
  if (!values.tCopyLayers) return;
  ValueCell.updateIfChanged(values.tCopyLayers, layers.image);
  ValueCell.updateIfChanged(values.dCopyLayers!, layers.has);
}

/**
 * The particles of `template` that `loci` (of `list`) holds in every one of its copies `tops` (a copy whole all of
 * them, in part its locals): what its objects' own layers take (the rest a copy's, not its template's); none where a
 * copy holds none
 */
export function getTemplateWideLocals(
  list: ParticleList,
  template: ParticleList,
  tops: ArrayLike<number>,
  loci: ModelLoci,
): OrderedSet<number> | undefined {
  if (!Particle.isLoci(loci) || loci.particles !== list || tops.length === 0) return undefined;
  const partial = loci.partial;
  let byTop: Map<number, OrderedSet<number>> | undefined;
  if (partial) {
    byTop = new Map();
    for (let i = 0; i < partial.tops.length; ++i) byTop.set(partial.tops[i], partial.locals[i]);
  }
  let locals: OrderedSet<number> | undefined;
  for (let c = 0; c < tops.length; ++c) {
    if (OrderedSet.has(loci.indices, tops[c])) continue;
    const held = byTop?.get(tops[c]);
    if (!held) return undefined;
    locals = locals ? OrderedSet.intersect(locals, held) : held;
    if (OrderedSet.isEmpty(locals)) return undefined;
  }
  return locals ?? Interval.ofBounds(0, template.count);
}

/**
 * Template particles `locals` of `template` (of `list`) as its leaf list's (`getLeafList`), or of its coarse level
 * `level` (those in it)
 */
export function getLevelLoci(
  list: ParticleList,
  template: ParticleList,
  locals: OrderedSet<number>,
  level?: ParticleNesting.TemplateLevel,
): Particle.Loci | undefined {
  if (!level) return Particle.Loci(ParticleNesting.getLeafList(list, template), locals);
  // a coarse level's: its particles of those template particles (by its inverse, not a pass over the level a visual)
  if (OrderedSet.size(locals) === template.count) {
    return Particle.Loci(level.list, Interval.ofBounds(0, level.list.count));
  }
  const inverse = getLevelInverse(template, level);
  const sample: number[] = [];
  OrderedSet.forEach(locals, (t) => {
    const k = inverse[t];
    if (k >= 0) sample.push(k);
  });
  if (sample.length === 0) return undefined;
  sample.sort((a, b) => a - b);
  return Particle.Loci(level.list, OrderedSet.ofSortedArray(sample));
}

/**
 * A list with templates' copies' markers (`copy-markers.ts`, its templates' copies one after another, a template's
 * particles by their index in it): what its objects' vertex stage reads by an instance's copy (`tCopyMarkers`), so
 * that a copy is marked apart from the others; a template's particles in every copy (the same in each) its objects'
 * own markers instead (no block a copy)
 */
export interface TemplateCopyMarkers {
  /** the image the objects read (`writeCopyMarkerImage`), and a template's first copy in it */
  readonly image: TextureImage<Uint8Array>;
  base(template: ParticleList): number;
  /**
   * whether an object of `template` drawing target `targetId` (its leaf list's id) shows a copy's marker: a
   * copy of it marked whole (or where not known), or in part at a particle of that target
   */
  isMarked(template: ParticleList, targetId: number): boolean;
  /** the copies of `templates` (in order), every marker cleared */
  layout(list: ParticleList, templates: readonly (readonly [ParticleList, TemplateCopies])[]): void;
  /**
   * `action` on `loci`'s copies (whole or in part, by copy; every copy for an action that removes on every loci),
   * the image written where changed; the templates whose every copy it holds the same particles of returned, with
   * them, for their objects' own markers
   */
  mark(
    list: ParticleList,
    loci: ModelLoci,
    action: MarkerAction,
  ): {
    changed: boolean;
    templateWide: Map<ParticleList, OrderedSet<number>>;
  };
}

export function createTemplateCopyMarkers(getBudget?: () => number): TemplateCopyMarkers {
  const markers = createCopyMarkers(getBudget);
  let image = writeCopyMarkerImage(undefined, markers.words, markers.bytes, markers.byteLength);
  const bases = new Map<ParticleList, number>();
  /** per top particle its copy in `markers` (-1: a leaf), and its template */
  let copyOfTop = new Int32Array(0);
  const templateOfCopy: ParticleList[] = [];
  const copyCounts = new Map<ParticleList, number>();
  /** per template its leaf list's targets (`getLeafList`: the ids its objects draw) */
  const leafTargets = new Map<ParticleList, ParticleList['targets']>();
  /**
   * per template the particles marked in every copy by its objects' own markers (`templateWide`), by kind: what a
   * removal from some of its copies moves into every copy's markers first (the others keep it), what a removal from
   * every copy takes from the copies' markers too (a copy's own marks of them)
   */
  const wide = new Map<ParticleList, Partial<Record<MarkKind, OrderedSet<number>>>>();
  /** advanced as the markers change: what a template's summary is of */
  let version = 0;
  const write = () => {
    image = writeCopyMarkerImage(image, markers.words, markers.bytes, markers.byteLength);
    version += 1;
  };
  /** per template its marked copies' summary: one marked whole (or not known where), the targets of the others' */
  const summaries = new Map<ParticleList, { version: number; whole: boolean; targets: Set<number> }>();
  const summarize = (template: ParticleList) => {
    const base = bases.get(template) ?? 0;
    const count = copyCounts.get(template) ?? 0;
    const targets = new Set<number>();
    let whole = false;
    for (let c = base; c < base + count && !whole; ++c) {
      if (!markers.isMarked(c)) continue;
      const positions = markers.getMarkedPositions(c);
      if (!positions) whole = true;
      else for (const i of positions) targets.add(ParticleValues.get(leafTargets.get(template)!, i));
    }
    return { version, whole, targets };
  };
  return {
    get image() {
      return image;
    },
    base: (template) => bases.get(template) ?? 0,
    isMarked(template, targetId) {
      if (markers.markedCount === 0) return false;
      let e = summaries.get(template);
      if (!e || e.version !== version) {
        e = summarize(template);
        summaries.set(template, e);
      }
      return e.whole || e.targets.has(targetId);
    },
    layout(list, templates) {
      bases.clear();
      copyCounts.clear();
      leafTargets.clear();
      summaries.clear();
      wide.clear();
      templateOfCopy.length = 0;
      copyOfTop = new Int32Array(list.count).fill(-1);
      let base = 0;
      const layout = templates.map(([template, copies]) => {
        bases.set(template, base);
        copyCounts.set(template, copies.tops.length);
        leafTargets.set(template, ParticleNesting.getLeafList(list, template).targets);
        copies.tops.forEach((top, c) => {
          copyOfTop[top] = base + c;
          templateOfCopy[base + c] = template;
        });
        base += copies.tops.length;
        // its particles by their index: pages from 0
        const pages = Int32Array.from({ length: Math.max(1, Math.ceil(template.count / StorePageSize)) }, (_, k) => k);
        return { pages, copyCount: copies.tops.length };
      });
      markers.layout(layout);
      write();
    },
    mark(list, loci, action) {
      const templateWide = new Map<ParticleList, OrderedSet<number>>();
      let changed = false;
      const kinds = getMarkKinds(action);
      const removing = isRemovingAction(action);
      if (isEveryLoci(loci)) {
        // every copy's own markers cleared by an action that removes (its objects' own take every particle)
        if (removing) {
          changed = markers.markAll(action);
          for (const w of wide.values()) for (const k of kinds) delete w[k];
        }
      } else if (Particle.isLoci(loci) && loci.particles === list && copyOfTop.length === list.count) {
        const partial = loci.partial;
        if (partial) {
          // a template's particles in every copy (the same set each, `ofTemplateParticles`): its objects' own
          const byTemplate = new Map<ParticleList, { count: number; locals: OrderedSet<number> | undefined }>();
          for (let i = 0; i < partial.tops.length; ++i) {
            const c = copyOfTop[partial.tops[i]];
            if (c < 0) continue;
            const t = templateOfCopy[c];
            const e = byTemplate.get(t);
            const locals = partial.locals[i];
            if (!e) byTemplate.set(t, { count: 1, locals });
            else {
              e.count += 1;
              if (e.locals && !OrderedSet.areEqual(e.locals, locals)) e.locals = undefined;
            }
          }
          for (const [t, e] of byTemplate) if (e.locals && e.count === copyCounts.get(t)) templateWide.set(t, e.locals);
        }
        const forEachCopy = (t: ParticleList, f: (c: number) => void) => {
          const base = bases.get(t) ?? 0;
          for (let c = base, ce = base + (copyCounts.get(t) ?? 0); c < ce; ++c) f(c);
        };
        // kept by kind: added to, or taken from in every copy (and from every copy's own markers)
        for (const [t, locals] of templateWide) {
          let w = wide.get(t);
          if (!w) wide.set(t, (w = {}));
          for (const k of kinds) {
            const was = w[k];
            const next = removing
              ? was && OrderedSet.subtract(was, locals)
              : was
                ? OrderedSet.union(was, locals)
                : locals;
            if (next && OrderedSet.size(next) > 0) w[k] = next;
            else delete w[k];
          }
          if (removing) {
            forEachCopy(t, (c) => {
              if (markers.markIndices(c, locals, action)) changed = true;
            });
          }
        }
        // a removal from some copies of a template marked in every copy: that mark moved into every copy's markers
        // first (its objects' own let go of it: `templateWide`), the removal from those copies' after
        const pushed = new Set<ParticleList>();
        const pushDown = (c: number) => {
          const t = templateOfCopy[c];
          const w = wide.get(t);
          if (!w || pushed.has(t) || templateWide.has(t)) return;
          let all: OrderedSet<number> | undefined;
          for (const k of kinds) {
            const marked = w[k];
            if (!marked) continue;
            const add = k === 'select' ? MarkerAction.Select : MarkerAction.Highlight;
            forEachCopy(t, (q) => markers.markIndices(q, marked, add));
            all = all ? OrderedSet.union(all, marked) : marked;
            delete w[k];
          }
          if (!all) return;
          pushed.add(t);
          templateWide.set(t, all);
          changed = true;
        };
        if (removing) {
          OrderedSet.forEach(loci.indices, (top) => {
            if (copyOfTop[top] >= 0) pushDown(copyOfTop[top]);
          });
          if (partial) for (const top of partial.tops) if (copyOfTop[top] >= 0) pushDown(copyOfTop[top]);
        }
        OrderedSet.forEach(loci.indices, (top) => {
          const c = copyOfTop[top];
          if (c >= 0 && markers.markCopy(c, action)) changed = true;
        });
        if (partial) {
          for (let i = 0; i < partial.tops.length; ++i) {
            const c = copyOfTop[partial.tops[i]];
            if (c < 0) continue;
            const t = templateOfCopy[c];
            if (templateWide.has(t) && !pushed.has(t)) continue;
            if (markers.markIndices(c, partial.locals[i], action)) changed = true;
          }
        }
      }
      if (changed) write();
      return { changed, templateWide };
    },
  };
}

/** what a template-wide mark is kept by (`TemplateCopyMarkers`): selected or highlighted */
type MarkKind = 'select' | 'highlight';

/** the kinds `action` adds to or takes from (none: a toggle, which a template-wide mark is not kept through) */
function getMarkKinds(action: MarkerAction): readonly MarkKind[] {
  switch (action) {
    case MarkerAction.Select:
    case MarkerAction.Deselect:
      return ['select'];
    case MarkerAction.Highlight:
    case MarkerAction.RemoveHighlight:
      return ['highlight'];
    case MarkerAction.Clear:
      return ['select', 'highlight'];
    default:
      return [];
  }
}

/** whether `action` only takes markers away */
function isRemovingAction(action: MarkerAction) {
  return action === MarkerAction.RemoveHighlight || action === MarkerAction.Deselect || action === MarkerAction.Clear;
}

type MarkerValues = {
  uMarker: ValueCell<number>;
  markerAverage: ValueCell<number>;
  markerStatus: ValueCell<number>;
  tMarker: ValueCell<TextureImage<Uint8Array>>;
};

/**
 * An expanded object's marker state with its copies' (`TemplateCopyMarkers.isMarked`): its own (its template
 * particles', what `Visual.mark` keeps) in `own` by object id, its uniform the texture's and its average at least a
 * part's where a copy is marked (its marking passes draw it, the vertex stage or-ing the copies' in). `restore` before
 * the visual marks (its own state as it left it), `apply` after
 */
export function createExpandedMarkerState() {
  const own = new Map<number, number>();
  return {
    restore(renderObject: GraphicsRenderObject) {
      const values = renderObject.values as unknown as MarkerValues;
      const average = own.get(renderObject.id);
      if (average !== undefined) ValueCell.updateIfChanged(values.markerAverage, average);
    },
    apply(renderObject: GraphicsRenderObject, copiesMarked: boolean) {
      const values = renderObject.values as unknown as MarkerValues;
      const average = values.markerAverage.ref.value;
      own.set(renderObject.id, average);
      // its uniform its own status where no copy is marked
      const status = values.markerStatus.ref.value;
      if (!copiesMarked) {
        ValueCell.updateIfChanged(values.uMarker, status);
        return;
      }
      // its layer read where its uniform stood for it (an upload the visual deferred: set again from its status)
      if (values.uMarker.ref.value !== -1) {
        if (status !== -1) ValueCell.update(values.tMarker, values.tMarker.ref.value);
        ValueCell.update(values.uMarker, -1);
      }
      ValueCell.updateIfChanged(values.markerAverage, average === 0 ? 0.5 : average);
    },
    clear() {
      own.clear();
    },
  };
}

const levelInverses = new WeakMap<ParticleNesting.TemplateLevel, Int32Array>();

/** per template particle its index in `level` (-1: not in it) */
function getLevelInverse(template: ParticleList, level: ParticleNesting.TemplateLevel) {
  let inverse = levelInverses.get(level);
  if (!inverse) {
    inverse = new Int32Array(template.count).fill(-1);
    const { indices } = level;
    for (let k = 0; k < indices.length; ++k) inverse[indices[k]] = k;
    levelInverses.set(level, inverse);
  }
  return inverse;
}

/**
 * The instance of a pick of an expanded object's draw slot `slot` (`v n + k`: pose `v` of those its last cull saw,
 * its instance `k`), `indices` the template particles its instances are: the copy's top particle and the template
 * particle; undefined where the slot is past what the cull saw
 */
export function getExpandedInstance(
  renderObject: GraphicsRenderObject,
  tops: readonly number[],
  indices: OrderedSet<number>,
  slot: number,
  /** a coarse level's: its particles' in the template (`TemplateLevel.indices`) */
  locals?: Int32Array,
) {
  const values = renderObject.values as unknown as InstanceExpansionValues & { instanceCount: ValueCell<number> };
  const n = values.instanceCount.ref.value;
  if (!values.instanceExpansion || n === 0 || slot < 0) return undefined;
  const { seen } = values.instanceExpansion.ref.value;
  const v = Math.floor(slot / n);
  const k = slot - v * n;
  if (v >= seen.count || k >= OrderedSet.size(indices)) return undefined;
  const i = OrderedSet.getAt(indices, k);
  return { top: tops[seen.order[v]], local: locals ? locals[i] : i };
}
