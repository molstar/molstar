/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import type { PluginContext } from '@molstar/plugin/context';
import type {
  ParticleListRef,
  ParticleRepresentationRef,
} from '@molstar/plugin/state/manager/particles/hierarchy-state';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticleListUnitcell3D } from '@molstar/plugin/state/transforms/particles/unitcell';
import { ParticleListFromTrajectory, ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { AnimateParticleTrajectory } from '@molstar/plugin/state/animation/built-in/particles';
import { particleTargetSourceFromPdbIdOrUrl } from '@molstar/plugin/state/helpers/particle-targets';
import { Asset } from '@molstar/core/util/assets';
import { setSubtreeVisibility } from '@molstar/plugin/behavior/static/state';
import { Interval } from '@molstar/core/data/int';
import { type State, type StateObjectCell, StateSelection, StateTree } from '@molstar/core/state';
import type { Color } from '@molstar/core/util/color';
import { decodeColor } from '@molstar/core/util/color/utils';
import { Particle } from '@molstar/model/model/particles/particle-list';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from '../../api/area.js';
import { NodeParam, OptionalParam } from '../../api/params.js';
import { FileParam } from '../../api/files.js';
import { messageOf } from '../../api/run.js';
import { seekWhileViewing } from '../../sources/trajectory.js';
import { nodesFor, updateEach } from '../../api/nodes.js';
import { entryRoot, isFailure, namesOf, pickOne } from '../../api/scope.js';
import { commitRepresentation, setRepresentation } from '../../api/set-representation.js';

// `particles` (api.md §2): the particle lists in the scene (`get`), their representations (desired
// state, like `representation.set`), visibility and removal. It follows `volume`: a `particles` names a list by label
// or state ref, else the only one loaded; several are `ambiguous`, with their labels in the hint.

/** The feature's plugin parts: Mol*'s particles, opt-in (`ParticlesRegistry`), which its operations run on */
export const ParticlesPlugin = { registry: ParticlesRegistry, behaviors: [] } as const;

const ListParam = () =>
  PD.Text('', { description: 'A particle list by label or state ref; empty: the only one loaded' });

export interface ParticleRepresentationInfo {
  readonly ref: string;
  readonly type: string;
  readonly label: string;
  readonly hidden: boolean;
}

export interface ParticleListInfo {
  readonly ref: string;
  readonly label: string;
  readonly hidden: boolean;
  readonly count: number;
  readonly representations: readonly ParticleRepresentationInfo[];
}

export class ParticlesApi extends ApiArea {
  readonly id = 'particles';

  readonly list = this.operation({
    id: 'particles.list',
    aliases: ['particles.get'],
    title: 'Describe the particles',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = lists(ctx);
      if (!list.length) return ApiResult.ok('No particles are loaded.', {}, { value: [] as ParticleListInfo[] });
      const names = namesOf(list);
      const infos: ParticleListInfo[] = list.map((l) => ({
        ref: l.cell.transform.ref,
        label: names.get(l.cell.transform.ref)!,
        hidden: !!l.cell.state.isHidden,
        count: l.cell.obj?.data.count ?? 0,
        representations: l.representations.map((r) => ({
          ref: r.cell.transform.ref,
          type: typeOf(r),
          label: r.cell.obj?.label ?? typeOf(r),
          hidden: !!r.cell.state.isHidden,
        })),
      }));
      const lines = infos.map((l) => {
        const reps = l.representations.map((r) => `${r.type}${r.hidden ? ' (hidden)' : ''}`);
        return `${l.label} (${l.ref}): ${l.count.toLocaleString('en')} particles${reps.length ? `; ${reps.join(', ')}` : '; no representations'}`;
      });
      return ApiResult.ok(lines.join('\n'), {}, { value: infos });
    },
  });

  /** @deprecated The former name of `list` (api.md §1 "Words"), until the release. */
  readonly get = this.list;

  // api.md §2: on a particle list, the representation of a type: added if absent, its params updated if present.
  readonly representation = this.operation({
    id: 'particles.representation',
    title: 'Show particles as…',
    params: (ctx) => {
      const types = ctx.plugin.representation.particles.registry.types;
      return {
        type: PD.Select(types[0]?.[0] ?? '', types, {
          description: 'A representation type: spacefill, orientation, fibers, target, …',
        }),
        params: PD.Value<Record<string, unknown>>({}, { description: "Some of the type's params, by name" }),
        particles: ListParam(),
      };
    },
    async run({ type, params, particles }, ctx) {
      const l = oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const name = listName(ctx, l);
      const data = listData(ctx, l);
      return setRepresentation(ctx, {
        name,
        type,
        patch: params,
        registry: ctx.plugin.representation.particles.registry,
        themes: ctx.plugin.representation.particles.themes,
        data,
        existing: l.representations,
        create: (provider) =>
          commitRepresentation(
            ctx,
            ctx.plugin.state.data
              .build()
              .to(l.cell)
              .apply(ParticlesRepresentation3D, representationParams(ctx.plugin, data!, type, params ?? {})),
            `Could not show ${name} as ${provider.label}`,
          ),
      });
    },
  });

  // The inspector's form: a representation's params as a whole, merged like `representation.update`.
  readonly update = this.operation({
    id: 'particles.update',
    title: 'Update a particle representation',
    params: {
      node: RepresentationNodeParam(),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The representation's params (type, color theme, size theme)" },
      ),
    },
    surface: { kind: 'ui-only', reason: "The inspector's form; tools use representation" },
    async run({ node, params }, ctx) {
      const reps = particleRepresentationsOf(ctx, node);
      if (isFailure(reps)) return reps;
      return updateEach(ctx, reps, params, 'Particle representation');
    },
  });

  readonly visibility = this.operation({
    id: 'particles.visibility',
    title: 'Show or hide particles',
    params: {
      visible: PD.Boolean(true),
      particles: ListParam(),
      representation: PD.Text('', { description: 'One representation by type instead of the whole list' }),
      node: RepresentationNodeParam(),
    },
    async run({ visible, particles, representation, node }, ctx) {
      if (node) {
        const reps = particleRepresentationsOf(ctx, node);
        if (isFailure(reps)) return reps;
        const todo = reps.filter((r) => !!r.cell.state.isHidden === visible);
        if (!todo.length) return ApiResult.already(`They are already ${visible ? 'shown' : 'hidden'}.`);
        ctx.plugin.managers.particles.hierarchy.toggleVisibility(todo, visible ? 'show' : 'hide');
        return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${namesOfReps(ctx, todo)}.`, {
          nodes: todo.map((r) => r.cell.transform.ref),
        });
      }
      const t = representation ? findRepresentation(ctx, particles, representation) : undefined;
      if (isFailure(t)) return t;
      const l = t ? undefined : oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const target = t ? t.rep : l!;
      const what = t ? `${t.rep.cell.obj?.label ?? typeOf(t.rep)} of ${listName(ctx, t.list)}` : listName(ctx, l!);
      if (!!target.cell.state.isHidden === !visible) {
        return ApiResult.already(`${what} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      ctx.plugin.managers.particles.hierarchy.toggleVisibility([target], visible ? 'show' : 'hide');
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${what}.`, { nodes: [target.cell.transform.ref] });
    },
  });

  readonly remove = this.operation({
    id: 'particles.remove',
    title: 'Remove particles',
    params: {
      particles: ListParam(),
      representation: PD.Text('', { description: 'One representation by type instead of the whole list' }),
      node: RepresentationNodeParam(),
    },
    approval: 'destructive',
    async run({ particles, representation, node }, ctx) {
      if (node) {
        const reps = particleRepresentationsOf(ctx, node);
        if (isFailure(reps)) return reps;
        await ctx.plugin.managers.particles.hierarchy.remove(reps, false);
        return ApiResult.ok(`Removed ${namesOfReps(ctx, reps)}.`, { nodes: reps.map((r) => r.cell.transform.ref) });
      }
      if (representation) {
        const t = findRepresentation(ctx, particles, representation);
        if (isFailure(t)) return t;
        await ctx.plugin.managers.particles.hierarchy.remove([t.rep], false);
        return ApiResult.ok(`Removed ${t.rep.cell.obj?.label ?? typeOf(t.rep)} from ${listName(ctx, t.list)}.`, {
          nodes: [t.rep.cell.transform.ref],
        });
      }
      const l = oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const name = listName(ctx, l);
      await ctx.plugin.managers.particles.hierarchy.remove([entryRoot(ctx, l.cell.transform.ref)], false);
      return ApiResult.ok(`Removed ${name}.`, { nodes: [l.cell.transform.ref] });
    },
  });

  // What the particles stand for (Mol*'s `ParticleListWithTargets`): one target per call, a file or a PDB id or URL,
  // instanced at the particles of that target id. Calls add up (the same target again replaces its object); a list
  // shows its targets once one is attached (a `target` representation is added when it has none).
  readonly targets = this.operation({
    id: 'particles.targets',
    title: 'Set what the particles stand for',
    params: {
      file: FileParam({ description: 'A structure, map or shape file for the target' }),
      source: PD.Text('', { description: 'A PDB id (1crn) or a URL, instead of a file' }),
      target: PD.Text('', { description: 'The target by id or name; empty: the only one' }),
      particles: ListParam(),
    },
    approval: 'network',
    addresses: ({ source }) => (/^https?:/i.test(source.trim()) ? [source.trim()] : []),
    async run({ file, source, target, particles }, ctx) {
      const l = oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const name = listName(ctx, l);
      const data = l.cell.obj?.data;
      if (!data) return ApiResult.fail('not-found', `${name} is not loaded.`);
      if (!file && !source.trim()) {
        return ApiResult.fail('invalid', 'Which object do the particles stand for?', 'Give a file, a PDB id or a URL.');
      }
      const id = targetIdOf(data, target);
      if (isFailure(id)) return id;
      const state = ctx.plugin.state.data;
      const decorator = state.select(
        StateSelection.Generators.byRef(l.cell.transform.ref).subtree().withTransformer(ParticleListWithTargets),
      )[0];
      const src = file
        ? { name: 'file' as const, params: { file: Asset.File(file), format: 'auto' } }
        : (() => {
            const t = particleTargetSourceFromPdbIdOrUrl(ctx.plugin, source);
            const url = t.kind === 'url' ? t.url : '';
            const asset = typeof url === 'string' ? Asset.Url(url) : url;
            return {
              name: 'url' as const,
              params: {
                url: asset,
                format: (t.kind === 'url' && t.format) || 'auto',
                isBinary: t.kind === 'url' ? t.isBinary : undefined,
              },
            };
          })();
      const entry = { targetId: id, source: src };
      const old = decorator?.transform.params as
        | { files?: unknown[]; targets: { targetId: number }[]; includeParent: boolean }
        | undefined;
      const targets = [...(old?.targets ?? []).filter((t) => t.targetId !== id), entry];
      const b = state.build();
      if (decorator) b.to(decorator).update({ ...old, targets } as never);
      else b.to(l.cell).apply(ParticleListWithTargets, { files: [], targets, includeParent: true } as never);
      try {
        await ctx.plugin.runTask(state.updateTree(b));
      } catch (e) {
        return ApiResult.fail('failed', `Could not attach the target of ${name}: ${messageOf(e)}`);
      }
      const what = file ? file.name : source.trim();
      if (!l.representations.some((r) => typeOf(r) === 'target')) {
        const added = await ctx.api.run('particles.representation', {
          particles: l.cell.transform.ref,
          type: 'target',
        });
        if (ApiResult.isFailure(added)) return added;
      }
      return ApiResult.ok(
        `The particles of ${name}${data.targetInfo.size > 1 ? ` (target ${id})` : ''} stand for ${what}.`,
        {
          nodes: [l.cell.transform.ref],
        },
      );
    },
  });

  // How a list is read from its file (the format's own transform: RELION's tomograms and micrographs, the pixel size,
  // the particle radius, an mmCIF assembly, a trajectory's frame): some of its params, merged; the inspector's Source.
  readonly source = this.operation({
    id: 'particles.source',
    title: 'Set how particles are read',
    surface: {
      kind: 'agent-only',
      reason: "The list's Source form in the inspector; a format's own params have no command form",
    },
    params: {
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The format's params, by name: tomograms, micrographs, pixelSize, particleRadius, …" },
      ),
      particles: ListParam(),
    },
    async run({ params, particles }, ctx) {
      const l = oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const name = listName(ctx, l);
      const cell = l.cell;
      const parent = ctx.plugin.state.data.cells.get(cell.transform.parent);
      const def = cell.transform.transformer.definition.params?.(parent?.obj as never, ctx.plugin) as
        | PD.Params
        | undefined;
      const unknown = Object.keys(params).filter((k) => !def || !(k in def));
      if (unknown.length) {
        return ApiResult.fail(
          'invalid',
          `${name} is not read with ${unknown.join(', ')}.`,
          def ? `Its params: ${Object.keys(def).join(', ')}.` : undefined,
        );
      }
      const old = cell.transform.params as Record<string, unknown>;
      const next = { ...old, ...params };
      if (def && PD.areEqual(def, old as never, next as never)) return ApiResult.already(`${name} is already read so.`);
      const state = ctx.plugin.state.data;
      try {
        await ctx.plugin.runTask(
          state.updateTree(
            state
              .build()
              .to(cell)
              .update(next as never),
          ),
        );
      } catch (e) {
        return ApiResult.fail('failed', `Could not read ${name} so: ${messageOf(e)}`);
      }
      const count = listData(ctx, l)?.count ?? 0;
      return ApiResult.ok(`${name} now has ${count.toLocaleString('en')} particles.`, { nodes: [cell.transform.ref] });
    },
  });

  // A particle trajectory's frames (a Simularium file), as `trajectory.*` does a structure's: viewing, no history step.
  readonly frame = this.operation({
    id: 'particles.frame',
    title: 'Show a particle frame',
    params: {
      index: OptionalParam(PD.Numeric(0, { min: 0, step: 1 }, { label: 'Frame', description: 'The frame, from 0' })),
      by: OptionalParam(
        PD.Numeric(
          1,
          { step: 1 },
          { label: 'Move by', description: 'Frames to move by (negative: back), wrapping around' },
        ),
      ),
      particles: ListParam(),
    },
    history: false,
    record: false,
    async run(v, ctx) {
      const { index, by, particles } = v as { index?: number; by?: number; particles: string };
      if ((index === undefined) === (by === undefined)) {
        return ApiResult.fail('invalid', 'Give the frame to show (`index`, from 0) or how far to move (`by`).');
      }
      const l = oneFramedList(ctx, particles);
      if (isFailure(l)) return l;
      const { count, current } = framesOf(ctx, l)!;
      const name = listName(ctx, l);
      if (index !== undefined && index >= count) {
        return ApiResult.fail('invalid', `${name} has ${count} frames.`, `Give 0 to ${count - 1}.`);
      }
      const to = index !== undefined ? index : (((current + by!) % count) + count) % count;
      if (to === current) return ApiResult.already(`${name} already shows frame ${to + 1} of ${count}.`);
      await seekWhileViewing(ctx.plugin, l.cell.transform.ref, {
        ...(l.cell.transform.params as object),
        frameIndex: to,
      });
      return ApiResult.ok(`${name} shows frame ${to + 1} of ${count}.`, { nodes: [l.cell.transform.ref] });
    },
  });

  readonly play = this.operation({
    id: 'particles.play',
    title: 'Play the particle frames',
    params: {
      mode: PD.Select<'loop' | 'palindrome' | 'once'>(
        'loop',
        PD.arrayToOptions(['loop', 'palindrome', 'once'] as const),
        {
          label: 'Mode',
          description: 'Loop, back and forth (palindrome), or once through',
        },
      ),
      seconds: PD.Numeric(5, { min: 1, max: 120, step: 0.5 }, { label: 'Seconds per pass' }),
      particles: ListParam(),
    },
    history: false,
    modules: ['canvas'],
    async run({ mode, seconds, particles }, ctx) {
      const l = oneFramedList(ctx, particles);
      if (isFailure(l)) return l;
      const name = listName(ctx, l);
      const ref = l.cell.transform.ref;
      if (isPlayingParticles(ctx.plugin, ref)) return ApiResult.already(`${name}'s frames are playing.`);
      await ctx.plugin.managers.animation.play(AnimateParticleTrajectory, {
        mode:
          mode === 'palindrome'
            ? { name: 'palindrome' as const, params: {} }
            : { name: mode, params: { direction: 'forward' as const } },
        duration: { name: 'fixed', params: { durationInS: seconds } },
        refs: [ref],
      });
      return ApiResult.ok(`Playing ${name}'s frames.`, { view: true, nodes: [ref] });
    },
  });

  readonly stop = this.operation({
    id: 'particles.stop',
    title: 'Stop playing the particle frames',
    params: {},
    history: false,
    async run(_, ctx) {
      if (!isPlayingParticles(ctx.plugin)) return ApiResult.already('No particle frames are playing.');
      await ctx.plugin.managers.animation.stop();
      return ApiResult.ok('Stopped playing the particle frames.', { view: true });
    },
  });

  // A representation's color theme (the row menu's Color): the theme with its defaults for the list, or `params`.
  readonly color = this.operation({
    id: 'particles.color',
    title: 'Color particles by…',
    params: (ctx) => ({
      theme: PD.Select('particle-index', ctx.plugin.representation.particles.themes.colorThemeRegistry.types, {
        description: 'A color theme: particle-index, particle-entity, particle-attribute, uniform, …',
      }),
      params: PD.Value<Record<string, unknown>>({}, { description: "Some of the theme's params, by name" }),
      particles: ListParam(),
      representation: PD.Text('', { description: 'One representation by type; empty: every one of the list' }),
      node: RepresentationNodeParam(),
    }),
    async run({ theme, params, particles, representation, node }, ctx) {
      const reps = representationsFor(ctx, particles, representation, node);
      if (isFailure(reps)) return reps;
      const { colorThemeRegistry } = ctx.plugin.representation.particles.themes;
      if (!colorThemeRegistry.has(theme)) {
        return ApiResult.fail('invalid', `Unknown color theme '${theme}'.`, 'See particles.color for the themes.');
      }
      const provider = colorThemeRegistry.get(theme);
      const todo = reps.filter((r) => {
        const data = listData(ctx, r.list);
        return !data || provider.isApplicable({ particles: data });
      });
      if (!todo.length) return ApiResult.fail('invalid', `${provider.label} does not apply to these particles.`);
      const update = ctx.plugin.state.data.build();
      let changed = 0;
      for (const r of todo) {
        const data = listData(ctx, r.list);
        const current = (r.cell.transform.params as { colorTheme?: { name: string; params?: object } } | undefined)
          ?.colorTheme;
        const next = {
          name: theme,
          params: {
            ...PD.getDefaultValues(provider.getParams(data ? { particles: data } : {})),
            ...(current?.name === theme ? current.params : {}),
            ...params,
          },
        };
        if (current && JSON.stringify(current) === JSON.stringify(next)) continue;
        update.to(r.cell).update((old: any) => ({ ...old, colorTheme: next }));
        changed++;
      }
      if (!changed)
        return ApiResult.already(
          `${namesOfReps(ctx, todo)} already ${todo.length === 1 ? 'is' : 'are'} colored by ${provider.label}.`,
        );
      await ctx.plugin.runTask(ctx.plugin.state.data.updateTree(update));
      return ApiResult.ok(`Colored ${namesOfReps(ctx, todo)} by ${provider.label}.`, {
        nodes: todo.map((r) => r.cell.transform.ref),
      });
    },
  });

  // The camera to a particle list (the row menu's Focus): a view change, no history step.
  readonly focus = this.operation({
    id: 'particles.focus',
    title: 'Focus the particles',
    params: { particles: ListParam() },
    history: false,
    async run({ particles }, ctx) {
      const l = oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const data = listData(ctx, l);
      if (!data?.count) return ApiResult.fail('not-found', `${listName(ctx, l)} has no particles.`);
      ctx.plugin.managers.camera.focusLoci(Particle.Loci(data, Interval.ofLength(data.count)));
      return ApiResult.ok(`Focused ${listName(ctx, l)}.`, { nodes: [l.cell.transform.ref] });
    },
  });

  // A list's unit cell (its periodic box), as `structure.unit-cell` does a crystal's: shown, hidden or styled; made when
  // shown and absent.
  readonly unitCell = this.operation({
    id: 'particles.unit-cell',
    title: 'Show the particle unit cell',
    params: {
      visible: OptionalParam(
        PD.Boolean(true, {
          description: "Show or hide the list's periodic box; absent: shown when nothing else is given",
        }),
      ),
      color: OptionalParam(PD.Text('', { description: "The cell's color, a CSS color (orange, #ff8800)" })),
      opacity: OptionalParam(
        PD.Numeric(1, { min: 0.05, max: 1, step: 0.05 }, { description: 'Opacity from 0.05 to 1' }),
      ),
      particles: ListParam(),
    },
    async run(v, ctx) {
      const { color, opacity, particles } = v as { color?: string; opacity?: number; particles: string };
      const given = v as { visible?: boolean };
      const l = oneList(ctx, particles || undefined);
      if (isFailure(l)) return l;
      const name = listName(ctx, l);
      if (!l.cell.obj?.data.cell) return ApiResult.fail('blocked', `${name} has no unit cell.`);
      let cellColor: Color | undefined;
      if (color?.trim()) {
        cellColor = decodeColor(color.trim());
        if (cellColor === undefined) {
          return ApiResult.fail('invalid', `'${color}' is not a color.`, 'Use a CSS color: orange, #ff8800.');
        }
      }
      const state = ctx.plugin.state.data;
      const existing = particleUnitCellOf(state, l.cell.transform.ref);
      const current = existing?.transform.params as { cellColor: Color; alpha: number } | undefined;
      const style = {
        ...(cellColor !== undefined && cellColor !== current?.cellColor ? { cellColor } : {}),
        ...(opacity !== undefined && opacity !== current?.alpha ? { alpha: opacity } : {}),
      };
      const restyle = Object.keys(style).length > 0;
      const shown = !!existing && !existing.state.isHidden;
      const visible = given.visible ?? (color === undefined && opacity === undefined ? true : shown || !existing);
      if (!restyle && shown === visible) {
        return ApiResult.already(`The unit cell of ${name} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      if (!existing) {
        const b = state
          .build()
          .to(l.cell)
          .apply(ParticleListUnitcell3D, style as never, { state: { isHidden: !visible } });
        await ctx.plugin.runTask(state.updateTree(b));
      } else {
        if (restyle) {
          await ctx.plugin.runTask(
            state.updateTree(
              state
                .build()
                .to(existing)
                .update({ ...current, ...style } as never),
            ),
          );
        }
        if (shown !== visible) setSubtreeVisibility(state, existing.transform.ref, !visible);
      }
      const now = particleUnitCellOf(state, l.cell.transform.ref);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} the unit cell of ${name}.`, {
        nodes: now ? [now.transform.ref] : [],
      });
    },
  });
}

/**
 * A list as its representations see it: the data at the end of its decorators (`ParticleListWithTargets` adds the
 * targets; the hierarchy's ref holds the list below them)
 */
function listData(ctx: ApiContext, l: ParticleListRef): ParticleList | undefined {
  const tree = ctx.plugin.state.data.tree;
  const root = StateTree.getDecoratorRoot(tree, l.cell.transform.ref);
  return (ctx.plugin.state.data.cells.get(root)?.obj?.data as ParticleList | undefined) ?? l.cell.obj?.data;
}

function lists(ctx: ApiContext): readonly ParticleListRef[] {
  return ctx.plugin.managers.particles.hierarchy.current.lists;
}

function listName(ctx: ApiContext, l: ParticleListRef): string {
  return namesOf(lists(ctx)).get(l.cell.transform.ref) ?? l.cell.transform.ref;
}

function oneList(ctx: ApiContext, name: string | undefined): ParticleListRef | ApiFailure {
  return pickOne(
    lists(ctx),
    name,
    'particle list',
    'Open a particle file (RELION STAR, Dynamo TBL, CryoET NDJSON, EM, …) first.',
  );
}

const RepresentationNodeParam = () =>
  NodeParam('Particle representations by selector (representation[type=spheres]) or a state ref');

/** The particle representations a node selector matches (api.md §1 "Nodes"). */
function particleRepresentationsOf(ctx: ApiContext, node: string): ParticleRepresentationRef[] | ApiFailure {
  const byRef = new Map(lists(ctx).flatMap((l) => l.representations.map((r) => [r.cell.transform.ref, r] as const)));
  return nodesFor(ctx.plugin.state.data, node, 'particle representation', (cell) => byRef.get(cell.transform.ref));
}

function namesOfReps(ctx: ApiContext, reps: readonly ParticleRepresentationRef[]): string {
  return reps.length === 1
    ? `${reps[0].cell.obj?.label ?? typeOf(reps[0])} of ${listName(ctx, reps[0].list)}`
    : `${reps.length} representations`;
}

/** A list's unit cell cell (`ParticleListUnitcell3D` below it), if any */
export function particleUnitCellOf(state: State, listRef: string): StateObjectCell | undefined {
  return state.select(StateSelection.Generators.byRef(listRef).subtree().withTransformer(ParticleListUnitcell3D))[0];
}

/** The representations a call names: by node, by type in a list, else every one of the list */
function representationsFor(
  ctx: ApiContext,
  particles: string,
  representation: string,
  node: string,
): ParticleRepresentationRef[] | ApiFailure {
  if (node) return particleRepresentationsOf(ctx, node);
  if (representation) {
    const t = findRepresentation(ctx, particles, representation);
    return isFailure(t) ? t : [t.rep];
  }
  const l = oneList(ctx, particles || undefined);
  if (isFailure(l)) return l;
  if (!l.representations.length) {
    return ApiResult.fail(
      'not-found',
      `${listName(ctx, l)} has no representations.`,
      'Add one with particles.representation.',
    );
  }
  return l.representations;
}

/** A target by id or name (an entity's), or the only one */
function targetIdOf(data: ParticleList, target: string): number | ApiFailure {
  const ids = [...data.targetInfo.keys()];
  const t = target.trim();
  if (!t) {
    if (ids.length === 1) return ids[0];
    return ApiResult.fail(
      'ambiguous',
      `The particles have ${ids.length} targets.`,
      `Name one: ${ids.slice(0, 12).join(', ')}.`,
    );
  }
  if (/^\d+$/.test(t) && data.targetInfo.has(+t)) return +t;
  for (const [id, info] of data.targetInfo) {
    const entity = info.entity !== undefined ? data.entityInfo?.get(info.entity)?.name : undefined;
    if ((info.name ?? entity)?.toLowerCase() === t.toLowerCase()) return id;
  }
  return ApiResult.fail('not-found', `There is no target '${t}'.`, `Targets: ${ids.slice(0, 12).join(', ')}.`);
}

/** A list read from a trajectory's frame (`ParticleListFromTrajectory`): its frame count and the frame shown */
export function framesOf(
  ctx: { plugin: PluginContext },
  l: ParticleListRef,
): { count: number; current: number } | undefined {
  if (l.cell.transform.transformer !== ParticleListFromTrajectory) return undefined;
  const parent = ctx.plugin.state.data.cells.get(l.cell.transform.parent)?.obj;
  const count = (parent?.data as { frameCount?: number } | undefined)?.frameCount ?? 0;
  if (count < 2) return undefined;
  return { count, current: (l.cell.transform.params as { frameIndex?: number }).frameIndex ?? 0 };
}

/** Whether particle frames play: a list's (`ref`), or any */
export function isPlayingParticles(plugin: PluginContext, ref?: string): boolean {
  const animation = plugin.managers.animation;
  if (!animation.isAnimating || animation.current?.anim.name !== AnimateParticleTrajectory.name) return false;
  if (ref === undefined) return true;
  const refs = (animation.current.paramValues as { refs?: readonly string[] }).refs;
  return !refs?.length || refs.includes(ref);
}

function oneFramedList(ctx: ApiContext, name: string): ParticleListRef | ApiFailure {
  const framed = lists(ctx).filter((l) => framesOf(ctx, l));
  if (!name) {
    if (framed.length === 1) return framed[0];
    if (!framed.length)
      return ApiResult.fail('not-found', 'No particles with frames are loaded.', 'Open a Simularium file.');
  }
  const l = oneList(ctx, name || undefined);
  if (isFailure(l)) return l;
  if (!framesOf(ctx, l)) return ApiResult.fail('blocked', `${listName(ctx, l)} has one frame.`);
  return l;
}

function typeOf(r: ParticleRepresentationRef): string {
  return (r.cell.transform.params as { type?: { name?: string } } | undefined)?.type?.name ?? '';
}

/** A representation by type, among a list's (or every list's when only the representation is named). */
function findRepresentation(
  ctx: ApiContext,
  listName_: string,
  name: string,
): { rep: ParticleRepresentationRef; list: ParticleListRef } | ApiFailure {
  let pool: ParticleRepresentationRef[];
  if (listName_) {
    const l = oneList(ctx, listName_);
    if (isFailure(l)) return l;
    pool = l.representations;
  } else {
    pool = lists(ctx).flatMap((l) => l.representations);
  }
  const byType = pool.filter((r) => typeOf(r) === name);
  if (byType.length === 1) return { rep: byType[0], list: byType[0].list };
  if (byType.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `There are ${byType.length} ${name} representations.`,
      `Pick one with node: ${byType.map((r) => `#${r.cell.transform.ref}`).join(', ')}.`,
    );
  }
  return ApiResult.fail(
    'not-found',
    `There is no ${name} representation.`,
    pool.length
      ? `Representations: ${pool.map((r) => `${typeOf(r)} (${r.cell.transform.ref})`).join(', ')}.`
      : 'Add one with particles.representation.',
  );
}

/** The transform's params for a type with some of its params set: the type's default themes, as the registry has them. */
function representationParams(
  plugin: PluginContext,
  data: ParticleList,
  type: string,
  patch: Record<string, unknown>,
): never {
  const { registry, themes } = plugin.representation.particles;
  const provider = registry.get(type);
  const dataCtx = { particles: data };
  const typeParams = Object.assign(PD.getDefaultValues(provider.getParams(themes, data)), patch);
  const color = themes.colorThemeRegistry.get(provider.defaultColorTheme.name);
  const size = themes.sizeThemeRegistry.get(provider.defaultSizeTheme.name);
  return {
    type: { name: provider.name, params: typeParams },
    colorTheme: {
      name: color.name,
      params: Object.assign(PD.getDefaultValues(color.getParams(dataCtx)), provider.defaultColorTheme.props),
    },
    sizeTheme: {
      name: size.name,
      params: Object.assign(PD.getDefaultValues(size.getParams(dataCtx)), provider.defaultSizeTheme.props),
    },
  } as never;
}
