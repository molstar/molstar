/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { StateObjectCell } from '@molstar/core/state';
import { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureElement } from '@molstar/model/model/structure';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import type { PluginContext } from '@molstar/plugin/context';
import { selectNodes } from '@molstar/prism-viewer/api';
import { watchedState } from '@molstar/prism-viewer/sources/watched';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from '@molstar/prism-viewer/api/area';
import { NodeParam, OptionalParam, TargetParam } from '@molstar/prism-viewer/api/params';
import {
  type Cut,
  type CutObject,
  type CutShape,
  CutShapeLabels,
  CutShapes,
  clippableCells,
  findCuts,
  placeCut,
  representationLabel,
  representationType,
  representationsSphere,
  rotationTo,
  withObject,
} from './cuts.js';

// The `clipping` area (_spec/README.md §2): cuts through the representations a target or a node selector reaches,
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
      types: PD.Value<readonly string[]>([], {
        description: 'Only these representation types (molecular-surface, cartoon, isosurface, …); empty: every one',
      }),
      invert: PD.Boolean(false, { description: 'Remove the other side (a plane) or the outside (a closed shape)' }),
    },
    async run({ shape, target, node, types, invert }, ctx) {
      const cells = await targetCells(ctx, target, node, types);
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
      const sphere = representationsSphere(cells);
      if (!sphere) return ApiResult.fail('not-found', 'What it would cut draws nothing yet.');
      const camera = ctx.plugin.canvas3d?.camera.state;
      const view = camera ? Vec3.sub(Vec3(), camera.position, sphere.center) : Vec3.create(0, 0, 1);
      const object = { ...placeCut(shape, sphere, view), invert };
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
      const { position, normal, rotation, scale, invert } = v as {
        position?: readonly number[];
        normal?: readonly number[];
        rotation?: { axis: readonly number[]; angle: number };
        scale?: readonly number[];
        invert?: boolean;
      };
      const cut = oneCut(ctx, index);
      if (ApiResult.isFailure(cut)) return cut;
      const o = cut.object;
      const next: CutObject = {
        ...o,
        ...(position ? { position: Vec3.fromArray(Vec3(), position, 0) } : {}),
        ...(rotation ? { rotation: { axis: Vec3.fromArray(Vec3(), rotation.axis, 0), angle: rotation.angle } } : {}),
        ...(normal ? { rotation: rotationTo(axisOf(o.type), Vec3.fromArray(Vec3(), normal, 0)) } : {}),
        ...(scale ? { scale: Vec3.fromArray(Vec3(), scale, 0) } : {}),
        ...(invert !== undefined ? { invert } : {}),
      };
      if (RepresentationClip.objectsEqual(o, next)) return ApiResult.already(`Cut ${index} is already so.`);
      await commit(ctx, cellsOf(ctx, cut), o, next);
      return ApiResult.ok(`Changed cut ${index}.`, { nodes: cut.targets.map((t) => t.ref) });
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
    .map((t) => t.label)
    .join(', ')}`;
}

/** A shape's local axis: a plane's normal and a cylinder's axis are +Y, a cone opens along −Z */
function axisOf(type: Clip.Props['objects'][number]['type']): Vec3 {
  return type === 'infiniteCone' ? Vec3.create(0, 0, -1) : Vec3.create(0, 1, 0);
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
