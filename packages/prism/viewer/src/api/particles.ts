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
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from './area.js';
import { mergeParams } from './params.js';
import { entryRoot, isFailure, namesOf, pickOne } from './scope.js';
import { commitRepresentation, setRepresentation } from './set-representation.js';

// `particles` (operations.md §2, plan 2.5): the particle lists in the scene (`get`), their representations (desired
// state, like `representation.set`), visibility and removal. It follows `volume`: a `particles` names a list by label
// or state ref, else the only one loaded; several are `ambiguous`, with their labels in the hint.

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

  readonly get = this.operation({
    id: 'particles.get',
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

  // operations.md §2: on a particle list, the representation of a type: added if absent, its params updated if present.
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
      const data = l.cell.obj?.data;
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
      representation: PD.Text('', { description: 'A particle representation by state ref' }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The representation's params (type, color theme, size theme)" },
      ),
    },
    surface: { kind: 'ui-only', reason: "The inspector's form; tools use representation" },
    async run({ representation, params }, ctx) {
      const rep = lists(ctx)
        .flatMap((l) => l.representations)
        .find((r) => r.cell.transform.ref === representation);
      if (!rep) return ApiResult.fail('not-found', `There is no particle representation '${representation}'.`);
      const current = rep.cell.transform.params as Record<string, unknown>;
      const next = mergeParams(current, params);
      if (JSON.stringify(next) === JSON.stringify(current)) {
        return ApiResult.already('The representation is already so.');
      }
      await ctx.plugin.state.updateTransform(ctx.plugin.state.data, representation, next, 'Particle representation');
      return ApiResult.ok('Updated the representation.', { nodes: [representation] });
    },
  });

  readonly visibility = this.operation({
    id: 'particles.visibility',
    title: 'Show or hide particles',
    params: {
      visible: PD.Boolean(true),
      particles: ListParam(),
      representation: PD.Text('', {
        description: 'One representation (by state ref or type) instead of the whole list',
      }),
    },
    async run({ visible, particles, representation }, ctx) {
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
      representation: PD.Text('', {
        description: 'One representation (by state ref or type) instead of the whole list',
      }),
    },
    approval: 'destructive',
    async run({ particles, representation }, ctx) {
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

function typeOf(r: ParticleRepresentationRef): string {
  return (r.cell.transform.params as { type?: { name?: string } } | undefined)?.type?.name ?? '';
}

/** A representation by state ref, else by type, among a list's (or every list's when only the representation is named). */
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
  const byRef = pool.find((r) => r.cell.transform.ref === name);
  if (byRef) return { rep: byRef, list: byRef.list };
  const byType = pool.filter((r) => typeOf(r) === name);
  if (byType.length === 1) return { rep: byType[0], list: byType[0].list };
  if (byType.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `There are ${byType.length} ${name} representations.`,
      `Name one by state ref: ${byType.map((r) => r.cell.transform.ref).join(', ')}.`,
    );
  }
  return ApiResult.fail(
    'not-found',
    `There is no particle representation '${name}'.`,
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
