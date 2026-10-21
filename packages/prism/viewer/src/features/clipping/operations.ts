/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { StateObjectCell } from '@molstar/core/state';
import { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import type { PluginContext } from '@molstar/plugin/context';
import { selectNodes } from '../../api.js';
import { watchedState } from '../../sources/watched.js';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult, type ViewerApi } from '../../api/area.js';
import { NodeParam, OptionalParam, TargetParam } from '../../api/params.js';
import { ClipGizmos } from '../../sources/clipping/gizmo/controller.js';
import {
  type Cut,
  type CutObject,
  type CutShape,
  CutShapeLabels,
  CutShapes,
  clippableCells,
  cutCandidates,
  type CutTarget,
  findCuts,
  nearlyEqual,
  MinPlacedRadius,
  placeCut,
  representationLabel,
  representationType,
  representationsSphere,
  reshape,
  rotationTo,
  shapeAxis,
  withObject,
} from '../../sources/clipping/cuts.js';

// The `clipping` area (clipping.md §2): cuts through the representations a target or a node selector reaches,
// written into those representations' own params (the state tree is the record), read back as cuts. Each run is one
// step of the history.

const CutParam = () => PD.Numeric(1, { min: 1, step: 1 }, { description: 'The cut by its number (clipping.list)' });
const Vec3Param = (description: string) => OptionalParam(PD.Value<readonly number[]>([0, 0, 0], { description }));

export class ClippingApi extends ApiArea {
  readonly id = 'clipping';

  /** The cuts as the views follow them: read again when the state changes, while something watches */
  readonly cuts: ReadonlyState<readonly Cut[]> = watchedState(
    () => findCuts(this.ctx.plugin),
    () => [
      this.ctx.plugin.state.data.events.changed,
      this.ctx.plugin.managers.structure.component.events.optionsUpdated,
    ],
    { compare: sameCuts },
  );

  /** The cuts' gizmos in the canvas (the section's Show handles): their drags commit through the person's API */
  readonly gizmos = new ClipGizmos(this.ctx.plugin, this.cuts);

  /** Shows the gizmos, a drag committing through `api` (the person's: one `clipping.update`, one history step) */
  showGizmos(api: ViewerApi): boolean {
    return this.gizmos.show(async (from, next) => {
      // by its object: another cut added or removed during the drag shifts the numbers
      const cut = findCuts(this.ctx.plugin).find((c) => RepresentationClip.objectsEqual(c.object, from));
      if (!cut) return false;
      const r = await api.run('clipping.update', {
        cut: cut.index,
        position: [...next.position],
        rotation: { axis: [...next.rotation.axis], angle: next.rotation.angle },
        scale: [...next.scale],
      });
      return r.ok;
    });
  }

  /**
   * The section's filter (which of the selection's representation types a new cut goes through): every type but
   * `excluded`, or, set by Clip this… on a representation, only `only` (of the types the selection has by then)
   */
  readonly filter = observableState<{ readonly excluded: readonly string[]; readonly only?: readonly string[] }>({
    excluded: [],
  });

  /**
   * Clip this… on a row. A structure's component or representation: its atoms selected (the person's step), the filter
   * narrowed to the representation's type (a component's: all its types); the section's Add makes the cut. Anything
   * else that can be cut (a volume's or a shape's representation, a volume, particles: no atoms to select): a plane
   * through it now, and that cut. The section opening is the caller's.
   */
  async clipThis(api: ViewerApi, ref: string): Promise<false | { readonly cut?: Cut }> {
    const state = this.ctx.plugin.state.data;
    const cell = state.cells.get(ref);
    if (!cell) return false;
    const representation = RepresentationClip.isClippable(cell);
    const component = representation ? state.cells.get(cell.transform.parent) : cell;
    const data = component?.obj?.data;
    if (data instanceof Structure) {
      const r = await api.run('select.set', { target: { loci: Structure.toStructureElementLoci(data) } });
      if (!r.ok) return false;
      this.filter.set(representation ? { excluded: [], only: [representationType(cell)] } : { excluded: [] });
      return {};
    }
    const r = await api.run('clipping.add', { shape: 'plane', node: `#${ref}` });
    const cut = (r as { value?: Cut }).value;
    return r.ok && cut ? { cut } : false;
  }

  /** Clip this… offered: anything that can be cut, or holds what can (a structure's component, a volume) */
  canClipThis(ref: string): boolean {
    const state = this.ctx.plugin.state.data;
    const cell = state.cells.get(ref);
    if (!cell) return false;
    if (RepresentationClip.isClippable(cell)) return true;
    if (cell.obj?.data instanceof Structure) return true;
    return clippableCells(state).some((c) => isBelow(state, c.transform.ref, ref));
  }

  /** The cuts whose editor is open in the section, by number (a cut just made opens alone) */
  readonly open = observableState<readonly number[]>([]);

  /** What `cut` can be turned on for: its representations and those near them (clipping.apply) */
  candidates(cut: Cut): CutTarget[] {
    return cutCandidates(this.ctx.plugin.state.data, cut);
  }

  /** The representation types drawing the selection: what a new cut can be narrowed to (the section's filter) */
  readonly selectionTypes: ReadonlyState<readonly string[]> = watchedState(
    () => typesInSelection(this.ctx.plugin),
    () => [this.ctx.plugin.managers.structure.selection.events.changed, this.ctx.plugin.state.data.events.changed],
    { compare: (a, b) => a.length === b.length && a.every((t, i) => t === b[i]) },
  );

  readonly list = this.operation({
    id: 'clipping.list',
    title: 'List the cuts',
    params: {},
    read: true,
    async run(_, ctx) {
      const cuts = findCuts(ctx.plugin);
      if (!cuts.length) return ApiResult.ok('Nothing is cut.', {}, { value: cuts });
      return ApiResult.ok(cuts.map(describe).join('\n'), {}, { value: cuts });
    },
  });

  readonly add = this.operation({
    id: 'clipping.add',
    title: 'Cut…',
    params: {
      shape: PD.Select<CutShape>(
        'plane',
        PD.arrayToOptions([...CutShapes], (s) => CutShapeLabels[s]),
        {
          description: 'plane (removes the half its normal points to), sphere, cube, cylinder (their inside), cone',
        },
      ),
      target: TargetParam('What to cut: the representations drawing it (default: the selection)'),
      node: NodeParam('Or representations by node selector (volume > representation, #ref): volumes, shapes too'),
      everything: OptionalParam(
        PD.Boolean(false, { description: 'Every representation that can be cut (nothing selected, a map alone)' }),
      ),
      types: PD.Value<readonly string[]>([], {
        description: 'Only these representation types (molecular-surface, cartoon, isosurface, …); empty: every one',
      }),
      invert: PD.Boolean(false, { description: 'Remove the other side (a plane) or the outside (a closed shape)' }),
    },
    async run({ shape, target, node, everything, types, invert }, ctx) {
      const around: { sphere?: Sphere3D } = {};
      const cells = everything ? allCells(ctx, types) : await targetCells(ctx, target, node, types, around);
      if (ApiResult.isFailure(cells)) return cells;
      const full = cells.filter(
        (c) => RepresentationClip.get(c.transform.params)!.objects.length >= Clip.MaxObjectCount,
      );
      if (full.length) {
        return ApiResult.fail(
          'blocked',
          `${representationLabel(ctx.plugin.state.data, full[0])} has ${Clip.MaxObjectCount} clip objects already.`,
          'Remove a cut first.',
        );
      }
      // centered on what was targeted (its atoms), else on what the representations draw (a volume's, a shape's)
      const drawn = representationsSphere(cells);
      if (!drawn) return ApiResult.fail('not-found', 'What it would cut draws nothing yet.');
      const sphere = around.sphere ?? drawn;
      const camera = ctx.plugin.canvas3d?.camera.state;
      // along the line of sight (not toward the camera's position): pushed toward the viewer, it stays over the target
      const view = camera ? Vec3.sub(Vec3(), camera.position, camera.target) : Vec3.create(0, 0, 1);
      const object = { ...placeCut(shape, sphere, view, drawn), invert };
      // a cut like it on these already (Add pressed again): that one, not another on top of it
      const same = findCuts(ctx.plugin).find(
        (c) =>
          nearlyEqual(c.object, object) && cells.some((cell) => c.targets.some((t) => t.ref === cell.transform.ref)),
      );
      if (same) return ApiResult.already(`Cut ${same.index} is that cut already.`, { value: same });
      await commit(ctx, cells, undefined, object);
      const cut = findCuts(ctx.plugin).find((c) => RepresentationClip.objectsEqual(c.object, object));
      return ApiResult.ok(
        `Cut ${namesOf(ctx, cells)} with a ${CutShapeLabels[shape].toLowerCase()}${cut ? ` (cut ${cut.index})` : ''}.`,
        { nodes: cells.map((c) => c.transform.ref) },
        { value: cut },
      );
    },
  });

  readonly update = this.operation({
    id: 'clipping.update',
    title: 'Change a cut',
    params: {
      cut: CutParam(),
      shape: OptionalParam(
        PD.Select<CutShape>(
          'plane',
          PD.arrayToOptions([...CutShapes], (s) => CutShapeLabels[s]),
          {
            description: 'Another shape, in its place: its position, direction and size carried over',
          },
        ),
      ),
      position: Vec3Param('Its center (Å, world): a plane through it'),
      normal: Vec3Param("A plane's normal (it removes the side this points to); a shape's axis"),
      rotation: OptionalParam(
        PD.Value<{ axis: readonly number[]; angle: number }>(
          { axis: [1, 0, 0], angle: 0 },
          { description: 'Its rotation: an axis and degrees' },
        ),
      ),
      scale: Vec3Param('Its size (Å): a sphere or cube diameter or edge; a cylinder (diameter, height, diameter)'),
      invert: OptionalParam(PD.Boolean(false, { description: 'Remove the other side, or the outside' })),
    },
    async run(v, ctx) {
      const { cut: index } = v;
      const { shape, position, normal, rotation, scale, invert } = v as {
        shape?: CutShape;
        position?: readonly number[];
        normal?: readonly number[];
        rotation?: { axis: readonly number[]; angle: number };
        scale?: readonly number[];
        invert?: boolean;
      };
      const cut = oneCut(ctx, index);
      if (ApiResult.isFailure(cut)) return cut;
      // a direction has a length (the turn gesture and the shader take a unit axis)
      const zero = (v?: readonly number[]) => !!v && Math.hypot(v[0], v[1], v[2]) < 1e-9;
      if (zero(normal) || zero(rotation?.axis)) {
        return ApiResult.fail('invalid', `A ${normal ? 'normal' : 'rotation axis'} of length 0 has no direction.`);
      }
      const o = cut.object;
      // another shape first (sized from what it cuts when it had none), the other values over it
      const shaped = shape
        ? reshape(o, shape, Math.max(representationsSphere(cellsOf(ctx, cut))?.radius ?? 0, MinPlacedRadius))
        : o;
      const next: CutObject = {
        ...shaped,
        ...(position ? { position: Vec3.fromArray(Vec3(), position, 0) } : {}),
        ...(rotation
          ? {
              rotation: {
                axis: Vec3.normalize(Vec3(), Vec3.fromArray(Vec3(), rotation.axis, 0)),
                angle: rotation.angle,
              },
            }
          : {}),
        ...(normal ? { rotation: rotationTo(shapeAxis(shaped.type), Vec3.fromArray(Vec3(), normal, 0)) } : {}),
        ...(scale ? { scale: Vec3.fromArray(Vec3(), scale, 0) } : {}),
        ...(invert !== undefined ? { invert } : {}),
      };
      if (RepresentationClip.objectsEqual(o, next)) return ApiResult.already(`Cut ${index} is already so.`);
      await commit(ctx, cellsOf(ctx, cut), o, next);
      return ApiResult.ok(`Changed cut ${index}.`, { nodes: cut.targets.map((t) => t.ref) });
    },
  });

  readonly apply = this.operation({
    id: 'clipping.apply',
    title: 'Cut representations or not',
    params: {
      cut: CutParam(),
      representation: PD.Text('', { description: "A representation's node ref (clipping.list: a cut's targets)" }),
      representations: OptionalParam(
        PD.Value<readonly string[]>([], { description: 'Several, as one step (a type of them all at once)' }),
      ),
      on: PD.Boolean(true, { description: 'Cut them (true) or leave them whole (false)' }),
    },
    async run({ cut: index, representation, representations, on }, ctx) {
      const cut = oneCut(ctx, index);
      if (ApiResult.isFailure(cut)) return cut;
      const state = ctx.plugin.state.data;
      const refs = [
        ...new Set([
          ...(representation ? [representation] : []),
          ...((representations as readonly string[] | undefined) ?? []),
        ]),
      ];
      if (!refs.length)
        return ApiResult.fail('invalid', 'No representation given.', 'representation, or representations.');
      const cells: StateObjectCell[] = [];
      for (const ref of refs) {
        const cell = state.cells.get(ref);
        if (!cell || !RepresentationClip.isClippable(cell)) {
          return ApiResult.fail('not-found', `There is no representation ${ref} to cut.`);
        }
        // only those it changes
        if (cut.targets.some((t) => t.ref === ref) !== on) cells.push(cell);
      }
      const label = cells.length === 1 ? representationLabel(state, cells[0]) : `${cells.length} representations`;
      if (!cells.length) return ApiResult.already(`Cut ${index} is ${on ? 'already on them' : 'not on them'}.`);
      if (!on && cells.length >= cut.targets.length) {
        return ApiResult.fail(
          'blocked',
          `That is all cut ${index} is on.`,
          'Remove the cut instead (clipping.remove).',
        );
      }
      const full = on
        ? cells.find((c) => RepresentationClip.get(c.transform.params)!.objects.length >= Clip.MaxObjectCount)
        : undefined;
      if (full) {
        return ApiResult.fail(
          'blocked',
          `${representationLabel(state, full)} has ${Clip.MaxObjectCount} clip objects already.`,
        );
      }
      await commit(ctx, cells, on ? undefined : cut.object, on ? cut.object : undefined);
      return ApiResult.ok(on ? `Cut ${index} now cuts ${label}.` : `Cut ${index} leaves ${label} whole.`, {
        nodes: cells.map((c) => c.transform.ref),
      });
    },
  });

  readonly remove = this.operation({
    id: 'clipping.remove',
    title: 'Remove a cut',
    params: { cut: CutParam() },
    approval: 'destructive',
    async run({ cut: index }, ctx) {
      const cut = oneCut(ctx, index);
      if (ApiResult.isFailure(cut)) return cut;
      await commit(ctx, cellsOf(ctx, cut), cut.object, undefined);
      return ApiResult.ok(`Removed cut ${index}.`, { nodes: cut.targets.map((t) => t.ref) });
    },
  });

  readonly clear = this.operation({
    id: 'clipping.clear',
    title: 'Remove every cut',
    params: {},
    approval: 'destructive',
    async run(_, ctx) {
      const cuts = findCuts(ctx.plugin);
      if (!cuts.length) return ApiResult.already('Nothing is cut.');
      const state = ctx.plugin.state.data;
      const b = state.build();
      const objects = cuts.map((c) => c.object);
      for (const cell of clippableCells(state)) {
        const clip = RepresentationClip.get(cell.transform.params)!;
        const kept = clip.objects.filter((o) => !objects.some((c) => RepresentationClip.objectsEqual(o, c)));
        if (kept.length !== clip.objects.length) {
          b.to(cell).update(RepresentationClip.set(cell.transform.params, { ...clip, objects: kept }) as never);
        }
      }
      await ctx.plugin.runTask(state.updateTree(b));
      return ApiResult.ok(`Removed ${cuts.length === 1 ? 'the cut' : `${cuts.length} cuts`}.`);
    },
  });
}

/** A cut in words: "Cut 1: plane through (12.0, 3.1, −4.0) · cartoon of Polymer, molecular-surface of Ligand" */
function describe(cut: Cut): string {
  const p = [...cut.object.position].map((x) => x.toFixed(1)).join(', ');
  const what = CutShapeLabels[cut.object.type as CutShape] ?? cut.object.type;
  return `Cut ${cut.index}: ${what.toLowerCase()} at (${p})${cut.object.invert ? ', inverted' : ''} · ${cut.targets
    .map((t) => `${t.label} (${t.ref})`)
    .join(', ')}`;
}

function oneCut(ctx: ApiContext, index: number): Cut | ApiFailure {
  const cuts = findCuts(ctx.plugin);
  const cut = cuts.find((c) => c.index === index);
  if (cut) return cut;
  return ApiResult.fail(
    'not-found',
    cuts.length ? `There is no cut ${index}.` : 'Nothing is cut.',
    cuts.length ? `The cuts: 1 to ${cuts.length} (clipping.list).` : 'Add one with clipping.add.',
  );
}

function cellsOf(ctx: ApiContext, cut: Cut): StateObjectCell[] {
  return cut.targets.flatMap((t) => {
    const cell = ctx.plugin.state.data.cells.get(t.ref);
    return cell ? [cell] : [];
  });
}

/** The object `previous` replaced by `next` (added when absent, removed when undefined) on each cell, as one update */
async function commit(
  ctx: ApiContext,
  cells: readonly StateObjectCell[],
  previous: CutObject | undefined,
  next: CutObject | undefined,
) {
  const state = ctx.plugin.state.data;
  const b = state.build();
  for (const cell of cells) b.to(cell).update(withObject(cell.transform.params, previous, next) as never);
  await ctx.plugin.runTask(state.updateTree(b));
}

function namesOf(ctx: ApiContext, cells: readonly StateObjectCell[]): string {
  if (cells.length > 2) return `${cells.length} representations`;
  return cells.map((c) => representationLabel(ctx.plugin.state.data, c).toLowerCase()).join(' and ');
}

/**
 * The representations a call reaches: by node selector (a clippable representation, or those below a node: a volume's,
 * a shape's), else by target (the representations of the components the target's atoms are in; the selection by
 * default); of the given types only.
 */
export async function targetCells(
  ctx: ApiContext,
  target: unknown,
  node: string,
  types: readonly string[],
  /** Given, it gets the bounds of the targeted atoms (none for a node selector) */
  around?: { sphere?: Sphere3D },
): Promise<StateObjectCell[] | ApiFailure> {
  const state = ctx.plugin.state.data;
  const out = new Map<string, StateObjectCell>();
  if (node) {
    const matched = selectNodes(state, node);
    if (ApiResult.isFailure(matched)) return matched;
    const all = clippableCells(state);
    for (const cell of matched) {
      if (RepresentationClip.isClippable(cell)) out.set(cell.transform.ref, cell);
      else {
        for (const c of all) if (isBelow(state, c.transform.ref, cell.transform.ref)) out.set(c.transform.ref, c);
      }
    }
  }
  const text = typeof target === 'string' ? target.trim() : target;
  if (!node || text) {
    if (!text && ctx.plugin.managers.structure.selection.stats.elementCount === 0) {
      return ApiResult.fail(
        'not-found',
        'Nothing is selected.',
        'Select what to cut, or name it (target), or give representations (node).',
      );
    }
    const resolved = await ctx.targets.resolve(text || 'selection');
    if (ApiResult.isFailure(resolved)) return resolved;
    if (resolved.empty) {
      return ApiResult.fail(
        'not-found',
        text ? `${String(text)} selects nothing.` : 'Nothing is selected.',
        'Select what to cut, or name it (target), or give representations (node).',
      );
    }
    const structures = ctx.plugin.managers.structure.hierarchy.current.structures;
    for (const part of resolved.parts) {
      if (around) {
        const s = StructureElement.Loci.getBoundary(part.loci).sphere;
        around.sphere = around.sphere ? Sphere3D.expandBySphere(Sphere3D(), around.sphere, s) : Sphere3D.clone(s);
      }
      const s = structures.find((x) => x.cell.transform.ref === part.structure.ref);
      for (const component of s?.components ?? []) {
        const data = component.cell.obj?.data;
        if (!data || StructureElement.Loci.isEmpty(StructureElement.Loci.remap(part.loci, data))) continue;
        for (const r of component.representations) {
          if (RepresentationClip.isClippable(r.cell)) out.set(r.cell.transform.ref, r.cell);
        }
      }
    }
  }
  const cells = [...out.values()].filter((c) => !types.length || types.includes(representationType(c)));
  if (!cells.length) {
    return ApiResult.fail(
      'not-found',
      out.size ? `None of it is drawn as ${types.join(' or ')}.` : 'No representation draws it.',
      out.size ? `It is drawn as ${[...new Set([...out.values()].map(representationType))].join(', ')}.` : undefined,
    );
  }
  return cells;
}

/** Every representation that can be cut, of `types` (all when empty) */
function allCells(ctx: ApiContext, types: readonly string[]): StateObjectCell[] | ApiFailure {
  const cells = clippableCells(ctx.plugin.state.data).filter(
    (c) => !types.length || types.includes(representationType(c)),
  );
  return cells.length ? cells : ApiResult.fail('not-found', 'Nothing is drawn that can be cut.');
}

function isBelow(state: ApiContext['plugin']['state']['data'], ref: string, ancestor: string): boolean {
  for (let at = state.cells.get(ref)?.transform.parent; at; ) {
    if (at === ancestor) return true;
    const next = state.cells.get(at)?.transform.parent;
    if (!next || next === at) return false;
    at = next;
  }
  return false;
}

function sameCuts(a: readonly Cut[], b: readonly Cut[]): boolean {
  return (
    a.length === b.length &&
    a.every(
      (c, i) =>
        RepresentationClip.objectsEqual(c.object, b[i].object) &&
        c.targets.length === b[i].targets.length &&
        c.targets.every((t, j) => t.ref === b[i].targets[j].ref && t.label === b[i].targets[j].label),
    )
  );
}

/** The types of the clippable representations of the components the selection is in, in the order first met */
function typesInSelection(plugin: PluginContext): string[] {
  const selection = plugin.managers.structure.selection;
  if (selection.stats.elementCount === 0) return [];
  const types: string[] = [];
  for (const s of plugin.managers.structure.hierarchy.current.structures) {
    const data = s.cell.obj?.data;
    const loci = data && selection.getLoci(data);
    if (!loci || !StructureElement.Loci.is(loci) || StructureElement.Loci.isEmpty(loci)) continue;
    for (const c of s.components) {
      const cdata = c.cell.obj?.data;
      if (!cdata || StructureElement.Loci.isEmpty(StructureElement.Loci.remap(loci, cdata))) continue;
      for (const r of c.representations) {
        const t = representationType(r.cell);
        if (RepresentationClip.isClippable(r.cell) && !types.includes(t)) types.push(t);
      }
    }
  }
  return types;
}
