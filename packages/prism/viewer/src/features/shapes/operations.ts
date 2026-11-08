/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { StateObjectCell } from '@molstar/core/state';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import { getTransformFromParams, type TransformParam } from '@molstar/plugin/state/transforms/helpers';
import { TransformShapeProvider } from '@molstar/plugin/state/transforms/shape/transform';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { setSubtreeVisibility } from '@molstar/plugin/behavior/static/state';
import type { Loci as ModelLoci } from '@molstar/model/model/loci';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from '../../api/area.js';
import { selectNodes, updateEach } from '../../api/nodes.js';
import { NodeParam } from '../../api/params.js';
import { entryRoot, isFailure } from '../../api/scope.js';
import { type ShapeEntry, shapeEntries } from '../../sources/shapes.js';
import { currentShapeTransform, shapePivot, shapeMatrixOf, shapeTransformCell } from '../../sources/shape-transform.js';
import { IdentityTransform, rigidTransformParams, wantedTransform } from '../../api/rigid-transform.js';

// `shape` (api.md §2): the shapes loaded on their own (a mesh file, a primitives node: `sources/shapes.ts`), as the
// Structure panel's Shapes group shows them: `list`, `visibility`, `update` (the inspector's form), `remove`, `focus`.
// A `node` names a shape or a representation of one (`shape[label=tetra.obj]`, `#ref`); empty: every shape (or the
// only one, where one is meant). Their look (material, light, …) is `look.set`'s, by node.

/** A shape as `list` gives it. */
export interface ShapeInfo {
  readonly ref: string;
  readonly label: string;
  readonly hidden: boolean;
  readonly representations: readonly { readonly ref: string; readonly label: string; readonly hidden: boolean }[];
}

const ShapeNodeParam = (description = 'A shape or one of its representations by selector or ref; empty: every shape') =>
  NodeParam(description);

/** What a node selector names among the shapes: whole shapes and single representations. */
interface ShapePick {
  readonly shapes: readonly ShapeEntry[];
  readonly representations: readonly StateObjectCell[];
  /** Each representation's shape (read once a run: `shapeEntries` walks the tree). */
  readonly shapeOf: ReadonlyMap<StateObjectCell, ShapeEntry>;
}

function pickShapes(ctx: ApiContext, node: string): ShapePick | ApiFailure {
  const entries = shapeEntries(ctx.plugin.state.data);
  if (!entries.length) return ApiResult.fail('not-found', 'No shape is loaded.', 'Open a mesh file (OBJ, PLY).');
  const shapeOf = new Map<StateObjectCell, ShapeEntry>();
  for (const e of entries) for (const r of e.representations) shapeOf.set(r, e);
  if (!node) return { shapes: entries, representations: [], shapeOf };
  const cells = selectNodes(ctx.plugin.state.data, node);
  if (isFailure(cells)) return cells;
  const refs = new Set(cells.map((c) => c.transform.ref));
  const shapes = entries.filter((e) => refs.has(e.cell.transform.ref));
  const representations = entries.flatMap((e) =>
    refs.has(e.cell.transform.ref) ? [] : e.representations.filter((r) => refs.has(r.transform.ref)),
  );
  if (!shapes.length && !representations.length) {
    return ApiResult.fail(
      'not-found',
      `No shape matches '${node}'.`,
      `Shapes: ${entries.map((e) => `${e.label} (${e.cell.transform.ref})`).join(', ')}.`,
    );
  }
  return { shapes, representations, shapeOf };
}

/** The one shape a node names (a shape, or a representation's shape), or the only one loaded when empty */
function oneShape(ctx: ApiContext, node: string): ShapeEntry | ApiFailure {
  const p = pickShapes(ctx, node);
  if (isFailure(p)) return p;
  const shapes = new Set([...p.shapes, ...p.representations.map((r) => p.shapeOf.get(r)!)]);
  if (shapes.size === 1) return [...shapes][0];
  return ApiResult.fail(
    'ambiguous',
    `${shapes.size} shapes are loaded.`,
    `Pick one with node: ${[...shapes].map((e) => `${e.label} (#${e.cell.transform.ref})`).join(', ')}.`,
  );
}

/** "tetra.obj", "Mesh of tetra.obj", "3 shapes". */
function pickText(p: ShapePick): string {
  const names = [
    ...p.shapes.map((e) => e.label),
    ...p.representations.map((r) => {
      const e = p.shapeOf.get(r);
      return `${r.obj?.label ?? 'Representation'}${e ? ` of ${e.label}` : ''}`;
    }),
  ];
  return names.length <= 2 ? names.join(' and ') : `${names.length} shapes`;
}

export class ShapeApi extends ApiArea {
  readonly id = 'shape';

  readonly list = this.operation({
    id: 'shape.list',
    title: 'List the shapes',
    params: {},
    read: true,
    async run(_, ctx) {
      const value: ShapeInfo[] = shapeEntries(ctx.plugin.state.data).map((e) => ({
        ref: e.cell.transform.ref,
        label: e.label,
        hidden: !!e.cell.state.isHidden,
        representations: e.representations.map((r) => ({
          ref: r.transform.ref,
          label: r.obj?.label ?? 'Representation',
          hidden: !!r.state.isHidden,
        })),
      }));
      if (!value.length) return ApiResult.ok('No shape is loaded.', {}, { value });
      const lines = value.map(
        (s) =>
          `${s.label} (${s.ref}${s.hidden ? ', hidden' : ''})${s.representations.length ? `: ${s.representations.map((r) => `${r.label} (${r.ref})`).join(', ')}` : ''}`,
      );
      return ApiResult.ok(lines.join('\n'), {}, { value });
    },
  });

  readonly visibility = this.operation({
    id: 'shape.visibility',
    title: 'Show or hide a shape',
    params: { node: ShapeNodeParam(), visible: PD.Boolean(true) },
    async run({ node, visible }, ctx) {
      const p = pickShapes(ctx, node);
      if (isFailure(p)) return p;
      const cells = [...p.shapes.map((e) => e.cell), ...p.representations];
      const todo = cells.filter((c) => !!c.state.isHidden === visible);
      if (!todo.length) return ApiResult.already(`${pickText(p)} already ${visible ? 'shown' : 'hidden'}.`);
      for (const c of todo) setSubtreeVisibility(ctx.plugin.state.data, c.transform.ref, !visible);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${pickText(p)}.`, {
        nodes: todo.map((c) => c.transform.ref),
      });
    },
  });

  // The inspector's form: the whole params of a shape's representations (alpha, coloring, material, flat shading, …).
  readonly update = this.operation({
    id: 'shape.update',
    title: 'Update a shape representation',
    params: {
      node: ShapeNodeParam('Shape representations (or a shape: all of its) by selector or ref'),
      params: PD.Value<Record<string, unknown>>({}, { description: "The representation's params" }),
    },
    surface: { kind: 'ui-only', reason: "The inspector's form; tools use look.set" },
    async run({ node, params }, ctx) {
      const p = pickShapes(ctx, node);
      if (isFailure(p)) return p;
      const reps = [...p.shapes.flatMap((e) => e.representations), ...p.representations];
      if (!reps.length) return ApiResult.fail('not-found', `${pickText(p)} has no representation.`);
      return updateEach(
        ctx,
        reps.map((cell) => ({ cell })),
        params,
        'Shape representation',
      );
    },
  });

  readonly remove = this.operation({
    id: 'shape.remove',
    title: 'Remove a shape',
    params: { node: ShapeNodeParam('A shape or one of its representations by selector or ref; empty: the only one') },
    approval: 'destructive',
    async run({ node }, ctx) {
      const p = pickShapes(ctx, node);
      if (isFailure(p)) return p;
      if (!node && p.shapes.length > 1) {
        return ApiResult.fail(
          'ambiguous',
          `Which shape? ${p.shapes.length} are loaded.`,
          `Name one: ${p.shapes.map((e) => `${e.label} (#${e.cell.transform.ref})`).join(', ')}.`,
        );
      }
      const text = pickText(p);
      // a whole shape goes with what it was made from (its file's data), as a structure does
      const refs = [
        ...p.shapes.map((e) => entryRoot(ctx, e.cell.transform.ref)),
        ...p.representations.map((r) => r.transform.ref),
      ];
      const b = ctx.plugin.state.data.build();
      for (const ref of refs) b.delete(ref);
      await b.commit();
      return ApiResult.ok(`Removed ${text}.`, { nodes: refs });
    },
  });

  readonly focus = this.operation({
    id: 'shape.focus',
    title: 'Focus a shape',
    params: { node: ShapeNodeParam('A shape or one of its representations by selector or ref; empty: every shape') },
    history: false,
    async run({ node }, ctx) {
      const p = pickShapes(ctx, node);
      if (isFailure(p)) return p;
      const reps = [...p.shapes.flatMap((e) => e.representations), ...p.representations];
      const loci: ModelLoci[] = [];
      for (const r of reps) {
        const repr = (r.obj?.data as { repr?: { getAllLoci(): ModelLoci[] } } | undefined)?.repr;
        if (repr) loci.push(...repr.getAllLoci());
      }
      if (!loci.length) return ApiResult.fail('not-found', `${pickText(p)} shows nothing to focus.`);
      ctx.plugin.managers.camera.focusLoci(loci);
      return ApiResult.ok(`Focused ${pickText(p)}.`, { view: true });
    },
  });

  // A shape moved or turned as a whole (`TransformShapeProvider`), as `structure.transform` does a structure: its
  // values are the transform from where the shape was loaded; `compose` applies them on top of the current one.
  readonly transform = this.operation({
    id: 'shape.transform',
    title: 'Move or rotate a shape',
    params: {
      ...rigidTransformParams('shape'),
      node: ShapeNodeParam('A shape by selector or ref; empty: the only one'),
    },
    async run(v, ctx) {
      const { compose, node } = v;
      const e = oneShape(ctx, node);
      if (isFailure(e)) return e;
      const provider = e.cell.obj?.data as ShapeProvider<any, any, any> | undefined;
      if (!provider) return ApiResult.fail('blocked', `${e.label} is not loaded yet.`);
      const wanted = wantedTransform(v, 'shape.transform-reset');
      if (isFailure(wanted)) return wanted;
      const { plugin } = ctx;
      const state = plugin.state.data;
      const ref = e.cell.transform.ref;
      const existing = shapeTransformCell(state, ref);
      const now = await currentShapeTransform(plugin, ref, provider);
      // stored as a matrix: a turn about the shape's center needs the shape made, the handles read it at once; the
      // center the handles turn it about (as drawn), so a turn by the API and one by hand agree
      const pivot = await shapePivot(plugin, e, provider, now);
      let matrix = await shapeMatrixOf(plugin, wanted, provider, pivot);
      if (compose) {
        // turned about the center where it is now
        const center = Vec3.transformMat4(Vec3(), pivot, now);
        matrix = Mat4.mul(Mat4(), getTransformFromParams(wanted, center), now);
      }
      const next: TransformParam = { name: 'matrix', params: { data: matrix, transpose: false } };
      if (Mat4.areEqual(matrix, now, 1e-3)) {
        return ApiResult.already(
          existing ? `${e.label} is already placed so.` : `${e.label} is already where it was loaded.`,
        );
      }
      const build = state.build();
      if (existing) build.to(existing).update({ transform: next });
      else build.to(e.cell).insert(TransformShapeProvider, { transform: next });
      await plugin.runTask(state.updateTree(build));
      return ApiResult.ok(`Moved ${e.label}.`, { nodes: [ref] });
    },
  });

  readonly transformReset = this.operation({
    id: 'shape.transform-reset',
    title: 'Reset the shape’s transform',
    params: { node: ShapeNodeParam('A shape by selector or ref; empty: the only one') },
    async run({ node }, ctx) {
      const e = oneShape(ctx, node);
      if (isFailure(e)) return e;
      const provider = e.cell.obj?.data as ShapeProvider<any, any, any> | undefined;
      const state = ctx.plugin.state.data;
      const ref = e.cell.transform.ref;
      const existing = shapeTransformCell(state, ref);
      if (!existing || !provider || Mat4.isIdentity(await currentShapeTransform(ctx.plugin, ref, provider), 1e-6)) {
        return ApiResult.already(`${e.label} is already where it was loaded.`);
      }
      // removing the decorator would take its subtree (the representations) with it: it is set to identity instead
      await ctx.plugin.runTask(state.updateTree(state.build().to(existing).update({ transform: IdentityTransform })));
      return ApiResult.ok(`${e.label} is back where it was loaded.`, { nodes: [ref] });
    },
  });
}
