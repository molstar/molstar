/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { State, StateBuilder, StateObjectCell } from '@molstar/core/state';
import { Clip } from '@molstar/core/util/clip';
import { UUID } from '@molstar/core/util/uuid';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import {
  ClipVolumesParams,
  type ClipVolumesProps,
  ClipVolumesShape3D,
  ClipVolumesTag,
  clipVolumeTag,
} from '@molstar/plugin/state/transforms/shape/clip-volumes';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';
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
  rotated,
  rotationTo,
  shapeAxis,
  withObject,
} from '../../sources/clipping/cuts.js';
import { rotationOf } from '../../sources/clipping/gizmo/binding.js';
import {
  afterRemoval,
  layerCell,
  layerOf,
  MaxKeptIndex,
  structureOf,
  withKept,
  withoutKept,
  writeLayer,
} from '../../sources/clipping/kept.js';

// The `clipping` area (clipping.md §2): cuts through the representations a target or a node selector reaches,
// written into those representations' own params (the state tree is the record), read back as cuts. Each run is one
// step of the history.

/** A new cut's tag (`Clip` objects' `tag`): this clip manager's, unique */
export const CutTagPrefix = 'clip-';
export const newCutTag = () => `${CutTagPrefix}${UUID.create22()}`;

const CutParam = () => PD.Numeric(1, { min: 1, step: 1 }, { description: 'The cut by its number (clipping.list)' });
/** A cut by its tag instead of its number: what a view or a gizmo sends (a number shifts when a cut before it goes, and
 * runs wait their turn) */
const CutTagParam = () =>
  OptionalParam(PD.Text('', { description: 'The cut by its tag (clipping.list), instead of its number' }));
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

  /** How the cuts cut and what they leave (the section's Settings): read again when the options change */
  readonly cutSettings: ReadonlyState<ClipSettings> = watchedState(
    () => settingsOf(this.ctx.plugin),
    () => [this.ctx.plugin.managers.structure.component.events.optionsUpdated],
    { compare: (a, b) => a.variant === b.variant && a.solidSurface === b.solidSurface },
  );

  /** The cuts' gizmos in the canvas (the section's Show handles): their drags commit through the person's API */
  readonly gizmos = new ClipGizmos(this.ctx.plugin, this.cuts);

  /** Shows the gizmos, a drag committing through `api` (the person's: one `clipping.update`, one history step) */
  showGizmos(api: ViewerApi): boolean {
    return this.gizmos.show(async (from, next) => {
      // by its object: another cut added or removed during the drag shifts the numbers
      const cut = findCuts(this.ctx.plugin).find((c) => RepresentationClip.sameCut(c.object, from));
      if (!cut) return false;
      // by its tag where it has one: the run may wait its turn behind another that renumbers the cuts
      const r = await api.run('clipping.update', {
        ...(cut.object.tag ? { tag: cut.object.tag } : { cut: cut.index }),
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
   * Clip this… on a row: a plane through what the row draws, now (a component: each of its representations; a
   * representation, a volume, particles: theirs), and that cut; one like it already there is that one. The selection is
   * left as it is. The section opening, the cut active, is the caller's.
   */
  async clipThis(api: ViewerApi, ref: string): Promise<false | { readonly cut?: Cut }> {
    if (!this.ctx.plugin.state.data.cells.get(ref)) return false;
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

  /** The section's Add keeps the selection whole (clipping.add keepWhole) */
  readonly keepNew = observableState(false);

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

  // clipping.md §2 "Settings": the variant is the scene's (the component options' clip) and every cut representation's;
  // solid surface is the look's default (`look.set`: a representation with its own keeps it).
  readonly settings = this.operation({
    id: 'clipping.settings',
    title: 'Set how cuts cut',
    params: {
      variant: OptionalParam(
        PD.Select<Clip.Variant>(
          'pixel',
          [
            ['pixel', 'Pixel'],
            ['instance', 'Instance'],
          ],
          {
            description:
              'pixel: cut through, pixel by pixel; instance: leave out whole copies (instances) whose center is cut away',
          },
        ),
      ),
      solidSurface: OptionalParam(
        PD.Boolean(true, { description: 'Close a cut surface with a solid cap (also where the camera cuts it)' }),
      ),
    },
    async run({ variant, solidSurface }, ctx) {
      const manager = ctx.plugin.managers.structure.component;
      const options = manager.state.options;
      const state = ctx.plugin.state.data;
      const differ = variant
        ? clippableCells(state).filter((c) => {
            const clip = RepresentationClip.get(c.transform.params)!;
            return clip.objects.length > 0 && clip.variant !== variant;
          })
        : [];
      const setVariant = variant !== undefined && (options.clipObjects.variant !== variant || differ.length > 0);
      const setSolid = solidSurface !== undefined && options.solidSurface !== solidSurface;
      if (!setVariant && !setSolid) return ApiResult.already('The cuts are set so already.');
      if (options.clipObjects.variant !== (variant ?? options.clipObjects.variant) || setSolid) {
        await manager.setOptions({
          ...options,
          ...(setVariant ? { clipObjects: { ...options.clipObjects, variant: variant! } } : {}),
          ...(setSolid ? { solidSurface: solidSurface! } : {}),
        });
      }
      // the cut representations hold their own clip, which the options leave (RepresentationClip.replaceOwned)
      if (setVariant && differ.length) {
        const b = state.build();
        for (const cell of differ) {
          const clip = RepresentationClip.get(cell.transform.params)!;
          b.to(cell).update(RepresentationClip.set(cell.transform.params, { ...clip, variant: variant! }) as never);
        }
        await ctx.plugin.runTask(state.updateTree(b));
      }
      const said = [
        setVariant ? (variant === 'instance' ? 'cuts leave out whole copies' : 'cuts go through, pixel by pixel') : '',
        setSolid ? (solidSurface ? 'cut surfaces are solid' : 'cut surfaces are open') : '',
      ].filter(Boolean);
      return ApiResult.ok(`${said.join('; ').replace(/^./, (c) => c.toUpperCase())}.`, {
        nodes: differ.map((c) => c.transform.ref),
      });
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
      keepWhole: OptionalParam(
        PD.Boolean(false, { description: 'Keep the target atoms whole: the cut goes around them (structures only)' }),
      ),
      types: PD.Value<readonly string[]>([], {
        description: 'Only these representation types (molecular-surface, cartoon, isosurface, …); empty: every one',
      }),
      invert: PD.Boolean(false, { description: 'Remove the other side (a plane) or the outside (a closed shape)' }),
    },
    async run({ shape, target, node, everything, keepWhole, types, invert }, ctx) {
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
      // tagged as this clip manager's: the same object on its representations is one cut, and two cuts are two
      const object = { ...placeCut(shape, sphere, view, drawn), invert, tag: newCutTag() };
      // a cut like it on these already (Add pressed again): that one, not another on top of it
      const same = findCuts(ctx.plugin).find(
        (c) =>
          nearlyEqual(c.object, object) && cells.some((cell) => c.targets.some((t) => t.ref === cell.transform.ref)),
      );
      if (same) return ApiResult.already(`Cut ${same.index} is that cut already.`, { value: same });
      await commit(ctx, cells, undefined, object);
      let cut = findCuts(ctx.plugin).find((c) => RepresentationClip.sameCut(c.object, object));
      // and the target kept whole by it (its atoms: not by a node selector, not everything); not done: said
      const warnings: string[] = [];
      if (keepWhole && (everything || node)) {
        warnings.push('Nothing is kept whole: keepWhole keeps the target atoms, not a node selector or everything.');
      } else if (keepWhole && cut) {
        const r = await keepAtoms(ctx, cut.index, target, object.tag);
        if (r.ok) cut = findCuts(ctx.plugin).find((c) => RepresentationClip.sameCut(c.object, object));
        else warnings.push(`Nothing is kept whole: ${r.message}`);
      }
      return ApiResult.ok(
        `Cut ${namesOf(ctx, cells)} with a ${CutShapeLabels[shape].toLowerCase()}${cut ? ` (cut ${cut.index})` : ''}.`,
        { nodes: cells.map((c) => c.transform.ref) },
        { value: cut, ...(warnings.length ? { warnings } : {}) },
      );
    },
  });

  readonly update = this.operation({
    id: 'clipping.update',
    title: 'Change a cut',
    params: {
      cut: CutParam(),
      tag: CutTagParam(),
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
      const { tag, shape, position, normal, rotation, scale, invert } = v as {
        tag?: string;
        shape?: CutShape;
        position?: readonly number[];
        normal?: readonly number[];
        rotation?: { axis: readonly number[]; angle: number };
        scale?: readonly number[];
        invert?: boolean;
      };
      const cut = oneCut(ctx, index, tag);
      if (ApiResult.isFailure(cut)) return cut;
      // three finite numbers each (a size: positive), a finite angle
      const bad =
        badVector('position', position) ??
        badVector('normal', normal) ??
        badVector('rotation.axis', rotation?.axis) ??
        badVector('scale', scale, true);
      if (bad) return bad;
      if (rotation && !Number.isFinite(rotation.angle)) {
        return ApiResult.fail('invalid', 'rotation.angle takes a number (degrees).');
      }
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
        // turned from where its axis points now to the normal: its twist about the axis kept (a box's)
        ...(normal ? { rotation: turnedTo(shaped, Vec3.fromArray(Vec3(), normal, 0)) } : {}),
        ...(scale ? { scale: Vec3.fromArray(Vec3(), scale, 0) } : {}),
        ...(invert !== undefined ? { invert } : {}),
      };
      if (RepresentationClip.objectsEqual(o, next)) return ApiResult.already(`Cut ${cut.index} is already so.`);
      await commit(ctx, cellsOf(ctx, cut), o, next);
      return ApiResult.ok(`Changed cut ${cut.index}.`, { nodes: cut.targets.map((t) => t.ref) });
    },
  });

  readonly apply = this.operation({
    id: 'clipping.apply',
    title: 'Cut representations or not',
    params: {
      cut: CutParam(),
      tag: CutTagParam(),
      representation: PD.Text('', { description: "A representation's node ref (clipping.list: a cut's targets)" }),
      representations: OptionalParam(
        PD.Value<readonly string[]>([], { description: 'Several, as one step (a type of them all at once)' }),
      ),
      on: PD.Boolean(true, { description: 'Cut them (true) or leave them whole (false)' }),
    },
    async run({ cut: index, tag, representation, representations, on }, ctx) {
      const cut = oneCut(ctx, index, tag);
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
      if (!cells.length) return ApiResult.already(`Cut ${cut.index} is ${on ? 'already on them' : 'not on them'}.`);
      if (!on && cells.length >= cut.targets.length) {
        return ApiResult.fail(
          'blocked',
          `That is all cut ${cut.index} is on.`,
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
      return ApiResult.ok(on ? `Cut ${cut.index} now cuts ${label}.` : `Cut ${cut.index} leaves ${label} whole.`, {
        nodes: cells.map((c) => c.transform.ref),
      });
    },
  });

  readonly keep = this.operation({
    id: 'clipping.keep',
    title: 'Keep atoms whole',
    params: {
      cut: CutParam(),
      tag: CutTagParam(),
      target: TargetParam('The atoms the cut leaves whole (default: the selection)'),
    },
    run: ({ cut: index, tag, target }, ctx) => keepAtoms(ctx, index, target, tag),
  });

  readonly release = this.operation({
    id: 'clipping.release',
    title: 'Cut atoms kept whole again',
    params: { cut: CutParam(), tag: CutTagParam() },
    async run({ cut: index, tag }, ctx) {
      const cut = oneCut(ctx, index, tag);
      if (ApiResult.isFailure(cut)) return cut;
      if (!cut.kept) return ApiResult.already(`Cut ${cut.index} keeps nothing whole.`);
      const state = ctx.plugin.state.data;
      const b = state.build();
      for (const cell of cellsOf(ctx, cut)) {
        const at = objectIndex(cell, cut.object);
        if (!structureOf(cell) || at > MaxKeptIndex || !layerCell(state, cell.transform.ref)) continue;
        writeLayer(state, b, cell, withoutKept(layerOf(state, cell), 1 << at));
      }
      await ctx.plugin.runTask(state.updateTree(b));
      return ApiResult.ok(`Cut ${cut.index} cuts ${cut.kept.label} again.`, { nodes: cut.targets.map((t) => t.ref) });
    },
  });

  readonly show = this.operation({
    id: 'clipping.show',
    title: 'Draw a cut in the scene, or not',
    params: {
      cut: CutParam(),
      tag: CutTagParam(),
      on: PD.Boolean(true, { description: 'Draw it (its walls and edges, saved with the scene) or not' }),
      options: OptionalParam(
        PD.Value<Partial<ClipVolumesProps>>(
          {},
          {
            description:
              'How: walls (on, wallColor, wallOpacity), edges (on, edgeColor, edgeWidth Å), extent (× what it cuts)',
          },
        ),
      ),
    },
    async run({ cut: index, tag: byTag, on, options }, ctx) {
      let cut = oneCut(ctx, index, byTag);
      if (ApiResult.isFailure(cut)) return cut;
      // drawn by its tag: an untagged cut (set by code) is tagged first, as this clip manager's
      if (!cut.object.tag) {
        if (!on) return ApiResult.already(`Cut ${cut.index} is not drawn.`);
        const tagged = { ...cut.object, tag: newCutTag() };
        await commit(ctx, cellsOf(ctx, cut), cut.object, tagged);
        const now = oneCut(ctx, cut.index, tagged.tag);
        if (ApiResult.isFailure(now)) return now;
        cut = now;
      }
      const tag = cut.object.tag!;
      const state = ctx.plugin.state.data;
      const b = state.build();
      const drawing = cut.drawing ? state.cells.get(cut.drawing) : undefined;
      const values = (options as Partial<ClipVolumesProps> | undefined) ?? {};
      const wrong = badOptions(values);
      if (wrong) return wrong;
      if (!on) {
        if (!drawing) return ApiResult.already(`Cut ${cut.index} is not drawn.`);
        b.delete(drawing.transform.ref);
      } else if (drawing) {
        const current = drawing.transform.params as ClipVolumesProps;
        const next = { ...current, ...values, tags: [tag] };
        if (PD.areEqual(ClipVolumesParams, current, next))
          return ApiResult.already(`Cut ${cut.index} is drawn so already.`);
        b.to(drawing.transform.ref).update(next);
      } else {
        const host = cellsOf(ctx, cut)[0];
        if (!host) return ApiResult.fail('not-found', `Cut ${cut.index} is on nothing to draw it under.`);
        addDrawing(b, host, tag, values);
      }
      await ctx.plugin.runTask(state.updateTree(b));
      const said = !on ? 'is no longer drawn' : drawing ? 'is drawn differently' : 'is drawn in the scene';
      return ApiResult.ok(`Cut ${cut.index} ${said}.`, { nodes: cut.targets.map((t) => t.ref) });
    },
  });

  /** The drawing's options (without which cuts: a cut's drawing is its own) */
  readonly drawingParams = (({ tags: _, ...rest }) => rest)(ClipVolumesParams);

  /** A cut's drawing's options as they are, or the defaults */
  drawingValues(cut: Cut): ClipVolumesProps {
    const cell = cut.drawing ? this.ctx.plugin.state.data.cells.get(cut.drawing) : undefined;
    return { ...PD.getDefaultValues(ClipVolumesParams), ...((cell?.transform.params as ClipVolumesProps) ?? {}) };
  }

  readonly remove = this.operation({
    id: 'clipping.remove',
    title: 'Remove a cut',
    params: { cut: CutParam(), tag: CutTagParam() },
    approval: 'destructive',
    async run({ cut: index, tag }, ctx) {
      const cut = oneCut(ctx, index, tag);
      if (ApiResult.isFailure(cut)) return cut;
      await commit(ctx, cellsOf(ctx, cut), cut.object, undefined);
      return ApiResult.ok(`Removed cut ${cut.index}.`, { nodes: cut.targets.map((t) => t.ref) });
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
        const kept = clip.objects.filter((o) => !objects.some((c) => RepresentationClip.sameCut(o, c)));
        if (kept.length !== clip.objects.length) {
          const removed = clip.objects.flatMap((o, i) =>
            objects.some((c) => RepresentationClip.sameCut(o, c)) ? [i] : [],
          );
          remapLayer(state, b, cell, removed);
          b.to(cell).update(RepresentationClip.set(cell.transform.params, { ...clip, objects: kept }) as never);
        }
      }
      await ctx.plugin.runTask(state.updateTree(b));
      await reconcileDrawings(ctx);
      return ApiResult.ok(`Removed ${cuts.length === 1 ? 'the cut' : `${cuts.length} cuts`}.`);
    },
  });
}

/** How the cuts cut (the section's Settings, `clipping.settings`). */
export interface ClipSettings {
  readonly variant: Clip.Variant;
  readonly solidSurface: boolean;
}

/** The settings: the component options' (the scene's) */
export function settingsOf(plugin: PluginContext): ClipSettings {
  const options = plugin.managers.structure.component.state.options;
  return { variant: options.clipObjects.variant, solidSurface: options.solidSurface };
}

/** A cut in words: "Cut 1: plane through (12.0, 3.1, −4.0) · cartoon of Polymer, molecular-surface of Ligand" */
function describe(cut: Cut): string {
  const p = [...cut.object.position].map((x) => x.toFixed(1)).join(', ');
  const what = CutShapeLabels[cut.object.type as CutShape] ?? cut.object.type;
  return `Cut ${cut.index}: ${what.toLowerCase()} at (${p})${cut.object.invert ? ', inverted' : ''} · ${cut.targets
    .map((t) => `${t.label} (${t.ref})`)
    .join(', ')}`;
}

/** Drawing options that are not the drawing's (an unknown name, a value not of its kind, a number out of its range) */
function badOptions(values: Record<string, unknown>): ApiFailure | undefined {
  for (const [name, value] of Object.entries(values)) {
    const p = (ClipVolumesParams as Record<string, PD.Any>)[name];
    if (!p || name === 'tags') {
      return ApiResult.fail(
        'invalid',
        `${name} is not a drawing option.`,
        'walls, wallColor, wallOpacity, edges, edgeColor, edgeWidth, extent.',
      );
    }
    const kind = typeof p.defaultValue;
    const fits =
      typeof value === kind &&
      (kind !== 'number' ||
        (Number.isFinite(value) &&
          (p.type !== 'number' ||
            ((p.min === undefined || (value as number) >= p.min) &&
              (p.max === undefined || (value as number) <= p.max)))));
    if (!fits) return ApiResult.fail('invalid', `${name} does not take ${JSON.stringify(value)}.`);
  }
  return undefined;
}

/** A vector param that is not three finite numbers (a size: positive) */
function badVector(name: string, v: readonly number[] | undefined, positive = false): ApiFailure | undefined {
  if (v === undefined) return undefined;
  const ok = Array.isArray(v) && v.length === 3 && v.every((x) => Number.isFinite(x) && (!positive || x > 0));
  return ok ? undefined : ApiResult.fail('invalid', `${name} takes three ${positive ? 'positive ' : ''}numbers.`);
}

/** A cut's rotation turned so its axis points along `to`, by the shortest turn from where it points now */
function turnedTo(o: CutObject, to: Vec3): CutObject['rotation'] {
  const now = rotated(o.rotation, shapeAxis(o.type));
  const turn = rotationTo(now, to);
  const q = Quat.multiply(
    Quat(),
    Quat.setAxisAngle(Quat(), turn.axis, (turn.angle * Math.PI) / 180),
    Quat.setAxisAngle(Quat(), Vec3.normalize(Vec3(), o.rotation.axis), (o.rotation.angle * Math.PI) / 180),
  );
  return rotationOf(Quat.normalize(q, q));
}

function oneCut(ctx: ApiContext, index: number, tag?: string): Cut | ApiFailure {
  const cuts = findCuts(ctx.plugin);
  if (tag) {
    const tagged = cuts.find((c) => c.object.tag === tag);
    return tagged ?? ApiResult.fail('not-found', 'That cut is gone.', 'The cuts: clipping.list.');
  }
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
  for (const cell of cells) {
    // removed from it: its kept atoms' layer remapped in the same step (the objects after it moved up a place)
    if (previous && !next) {
      const at =
        RepresentationClip.get(cell.transform.params)?.objects.findIndex((o) =>
          RepresentationClip.sameCut(o, previous),
        ) ?? -1;
      if (at >= 0) remapLayer(state, b, cell, [at]);
    }
    // a cut written takes the scene's variant (clipping.settings), so every cut representation cuts the same way
    b.to(cell).update(withObject(cell.transform.params, previous, next, settingsOf(ctx.plugin).variant) as never);
  }
  await ctx.plugin.runTask(state.updateTree(b));
  await reconcileDrawings(ctx);
}

/** A structure representation's clipping layer after its clip objects at `removed` went (clipping.md §5) */
function remapLayer(state: State, b: StateBuilder.Root, cell: StateObjectCell, removed: readonly number[]) {
  if (!structureOf(cell) || !layerCell(state, cell.transform.ref)) return;
  let layer = layerOf(state, cell);
  // from the last: each removal shifts only what is after it
  for (const at of [...removed].sort((x, y) => y - x)) layer = afterRemoval(layer, at);
  writeLayer(state, b, cell, layer);
}

/**
 * The cuts' drawings where their cuts are: one whose representation no longer has its cut (it was turned off there,
 * or removed) is drawn under another the cut is on, or goes with the cut (one a cut: `clip-volume-<tag>`)
 */
async function reconcileDrawings(ctx: ApiContext) {
  const state = ctx.plugin.state.data;
  const b = state.build();
  let changed = false;
  for (const cell of state.cells.values()) {
    if (cell.transform.transformer !== ClipVolumesShape3D) continue;
    const key = cell.transform.tags?.find((t) => t !== ClipVolumesTag);
    const tag = key?.slice(clipVolumeTag('').length);
    if (!tag) continue;
    const has = (c: StateObjectCell | undefined) =>
      !!c && !!RepresentationClip.get(c.transform.params)?.objects.some((o) => o.tag === tag);
    if (has(state.cells.get(cell.transform.parent))) continue;
    changed = true;
    b.delete(cell.transform.ref);
    const other = clippableCells(state).find(has);
    if (other) addDrawing(b, other, tag, cell.transform.params as ClipVolumesProps);
  }
  if (changed) await ctx.plugin.runTask(state.updateTree(b));
}

/** A cut's drawing under a representation it is on: the clip volumes shape and its representation, tagged */
function addDrawing(b: StateBuilder.Root, host: StateObjectCell, tag: string, params: Partial<ClipVolumesProps>) {
  b.to(host)
    .apply(
      ClipVolumesShape3D,
      { ...PD.getDefaultValues(ClipVolumesParams), ...params, tags: [tag] },
      {
        tags: [ClipVolumesTag, clipVolumeTag(tag)],
      },
    )
    .apply(ShapeRepresentation3D, {}, { tags: [ClipVolumesTag] });
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

/** The target's atoms (the selection by default) kept whole by cut `index`, on each structure it is on */
async function keepAtoms(ctx: ApiContext, index: number, target: unknown, tag?: string) {
  const cut = oneCut(ctx, index, tag);
  if (ApiResult.isFailure(cut)) return cut;
  const atoms = await targetAtoms(ctx, target);
  if (ApiResult.isFailure(atoms)) return atoms;
  const state = ctx.plugin.state.data;
  const b = state.build();
  let touched = 0;
  let beyond = 0;
  for (const cell of cellsOf(ctx, cut)) {
    const structure = structureOf(cell);
    if (!structure) continue;
    const at = objectIndex(cell, cut.object);
    if (at > MaxKeptIndex) {
      beyond++;
      continue;
    }
    const loci = atomsOn(atoms, structure.root);
    if (!loci) continue;
    writeLayer(state, b, cell, withKept(layerOf(state, cell), loci, 1 << at));
    touched++;
  }
  if (!touched) {
    return beyond
      ? ApiResult.fail(
          'blocked',
          `Only the first ${MaxKeptIndex + 1} cuts on a representation can keep atoms whole.`,
          'Remove an earlier cut from it.',
        )
      : ApiResult.fail('not-found', `None of it is in what cut ${cut.index} is on (it keeps only atoms whole).`);
  }
  await ctx.plugin.runTask(state.updateTree(b));
  const kept = oneCut(ctx, cut.index, cut.object.tag);
  const label = ApiResult.isFailure(kept) ? '' : (kept.kept?.label ?? '');
  return ApiResult.ok(`Cut ${cut.index} keeps ${label || 'them'} whole.`, { nodes: cut.targets.map((t) => t.ref) });
}

/** A cut object's index on a representation (its group bit in the clipping layer) */
function objectIndex(cell: StateObjectCell, object: CutObject): number {
  return (
    RepresentationClip.get(cell.transform.params)?.objects.findIndex((o) => RepresentationClip.sameCut(o, object)) ?? -1
  );
}

/** The atoms a target names (the selection by default), each on its structure's root */
async function targetAtoms(ctx: ApiContext, target: unknown): Promise<StructureElement.Loci[] | ApiFailure> {
  const text = typeof target === 'string' ? target.trim() : target;
  if (!text && ctx.plugin.managers.structure.selection.stats.elementCount === 0) {
    return ApiResult.fail(
      'not-found',
      'Nothing is selected.',
      'Select the atoms to keep whole, or name them (target).',
    );
  }
  const resolved = await ctx.targets.resolve(text || 'selection');
  if (ApiResult.isFailure(resolved)) return resolved;
  if (resolved.empty) return ApiResult.fail('not-found', `${String(text || 'The selection')} selects nothing.`);
  return resolved.parts.map((p) => StructureElement.Loci.remap(p.loci, p.loci.structure.root));
}

/** Of `atoms`, those on `root` */
function atomsOn(atoms: readonly StructureElement.Loci[], root: Structure): StructureElement.Loci | undefined {
  let out: StructureElement.Loci | undefined;
  for (const a of atoms) {
    if (a.structure.root !== root) continue;
    out = out ? StructureElement.Loci.union(out, a) : a;
  }
  return out && !StructureElement.Loci.isEmpty(out) ? out : undefined;
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
        c.drawing === b[i].drawing &&
        c.kept?.label === b[i].kept?.label &&
        c.kept?.count === b[i].kept?.count &&
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
