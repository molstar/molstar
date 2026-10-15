/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection, StateTransform } from '@molstar/core/state';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Grid, Volume } from '@molstar/model/model/volume';
import { VolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/behavior';
import { InitVolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import type { VolumeRef, VolumeRepresentationRef } from '@molstar/plugin/state/manager/volume/hierarchy-state';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from './area.js';
import { checkValues } from './run.js';
import { isFailure, labelOf, oneStructure } from './scope.js';

// `volume` (operations.md §2, plan 2.5): the maps in the scene (`get`), a representation's iso value, the
// representations of a volume (desired state, like `representation.set`), visibility, removal, and volume streaming on
// a structure. A `volume` names one by label or state ref, else the only one loaded; several are `ambiguous`, with
// their labels in the hint. The data types share this shape: `particles` follows it.

const VolumeParam = () => PD.Text('', { description: 'A volume by label or state ref; empty: the only one loaded' });
const RepresentationParam = () =>
  PD.Text('', {
    description: "A volume representation by state ref or type (isosurface, …); empty: the volume's first isosurface",
  });

/** A volume's grid facts, as the classic UI and the inspector show them. */
export interface VolumeStats {
  readonly grid: readonly [number, number, number];
  /** Å per voxel along each grid axis. */
  readonly voxelSize: readonly [number, number, number];
  readonly min: number;
  readonly max: number;
  readonly mean: number;
  readonly sigma: number;
}

export interface VolumeRepresentationInfo {
  readonly ref: string;
  readonly type: string;
  readonly label: string;
  readonly hidden: boolean;
  /** For types that have one (isosurface, slice, dot). */
  readonly isoValue?: { readonly kind: 'relative' | 'absolute'; readonly relative: number; readonly absolute: number };
}

export interface VolumeInfo {
  readonly ref: string;
  readonly label: string;
  readonly hidden: boolean;
  /** What it was read from ('ccp4', 'dsn6', …); empty when unknown. */
  readonly source: string;
  readonly stats: VolumeStats;
  readonly representations: readonly VolumeRepresentationInfo[];
}

export class VolumeApi extends ApiArea {
  readonly id = 'volume';

  readonly get = this.operation({
    id: 'volume.get',
    title: 'Describe the volumes',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = volumes(ctx);
      // Streamed density hangs off its structure, not the volume hierarchy: said, so it is not "no volumes".
      const tree = ctx.plugin.state.data.tree;
      const streamed = ctx.plugin.managers.structure.hierarchy.current.structures
        .filter((s) => StateSelection.findTagInSubtree(tree, s.cell.transform.ref, VolumeStreaming.RootTag))
        .map((s) => labelOf(ctx, s));
      const streaming = streamed.length ? `Density is streamed for ${streamed.join(', ')} (volume.streaming).` : '';
      if (!list.length) {
        return ApiResult.ok(streaming || 'No volumes are loaded.', {}, { value: [] as VolumeInfo[] });
      }
      const names = volumeNames(list);
      const infos = list.map((v) => describeVolume(ctx, v, names.get(v.cell.transform.ref)!));
      const lines = infos.map((v) => {
        const reps = v.representations.map(
          (r) => `${r.type}${r.isoValue ? ` ${formatIso(r.isoValue)}` : ''}${r.hidden ? ' (hidden)' : ''}`,
        );
        return `${v.label} (${v.ref}): ${volumeFacts(v.stats)}${reps.length ? `; ${reps.join(', ')}` : '; no representations'}`;
      });
      return ApiResult.ok([...lines, ...(streaming ? [streaming] : [])].join('\n'), {}, { value: infos });
    },
  });

  // operations.md §2: the iso value of a representation (σ above the mean when relative).
  readonly isovalue = this.operation({
    id: 'volume.isovalue',
    title: 'Set the iso value',
    params: {
      value: PD.Numeric(1.5, { step: 0.05 }, { description: 'The iso value: in σ when relative, else in map units' }),
      kind: PD.Select<'relative' | 'absolute'>(
        'relative',
        [
          ['relative', 'Relative (σ)'],
          ['absolute', 'Absolute'],
        ],
        { description: 'relative: σ above the mean; absolute: the map value' },
      ),
      volume: VolumeParam(),
      representation: RepresentationParam(),
    },
    async run({ value, kind, volume, representation }, ctx) {
      if (!Number.isFinite(value)) return ApiResult.fail('invalid', "'value' takes a number.");
      const r = resolveRepresentation(ctx, volume, representation, true);
      if (isFailure(r)) return r;
      const stats = statsOf(r.volume);
      const rangeProblem = checkRange(stats, kind, value);
      if (rangeProblem) return rangeProblem;

      const params = r.rep.cell.transform.params as RepresentationValues;
      const want = kind === 'relative' ? Volume.IsoValue.relative(value) : Volume.IsoValue.absolute(value);
      const old = params.type.params.isoValue as Volume.IsoValue | undefined;
      const name = volumeLabel(ctx, r.volume);
      if (old && sameIso(old, want)) {
        return ApiResult.already(`${name} is already at ${Volume.IsoValue.toString(want)}.`);
      }
      await ctx.plugin.state.updateTransform(
        ctx.plugin.state.data,
        r.rep.cell.transform.ref,
        { ...params, type: { name: params.type.name, params: { ...params.type.params, isoValue: want } } },
        'Iso value',
      );
      return ApiResult.ok(`Set the iso value of ${name} to ${Volume.IsoValue.toString(want)}.`, {
        nodes: [r.rep.cell.transform.ref],
      });
    },
  });

  // operations.md §2: on a volume, the representation of a type: added if absent, its params updated if present.
  readonly representation = this.operation({
    id: 'volume.representation',
    title: 'Show a volume as…',
    params: (ctx) => {
      const types = ctx.plugin.representation.volume.registry.types;
      return {
        type: PD.Select(types.some((t) => t[0] === 'isosurface') ? 'isosurface' : (types[0]?.[0] ?? ''), types, {
          description: 'A representation type: isosurface, direct-volume, slice, …',
        }),
        params: PD.Value<Record<string, unknown>>({}, { description: "Some of the type's params, by name" }),
        volume: VolumeParam(),
      };
    },
    async run({ type, params, volume }, ctx) {
      const v = oneVolume(ctx, volume || undefined);
      if (isFailure(v)) return v;
      const data = v.cell.obj?.data;
      if (!data) return ApiResult.fail('blocked', `${volumeLabel(ctx, v)} has no data to show.`);
      const registry = ctx.plugin.representation.volume.registry;
      if (!registry.has(type)) {
        return ApiResult.fail(
          'invalid',
          `There is no representation type '${type}'.`,
          `Types: ${registry.types.map((t) => t[0]).join(', ')}.`,
        );
      }
      const provider = registry.get(type);
      if (!provider.isApplicable(data)) {
        return ApiResult.fail(
          'blocked',
          `${provider.label} does not apply to ${volumeLabel(ctx, v)}.`,
          `Applicable: ${registry
            .getApplicableTypes(data)
            .map((t) => t[0])
            .join(', ')}.`,
        );
      }
      const typeParams = provider.getParams(ctx.plugin.representation.volume.themes, data);
      const patch = params ?? {};
      const checked = checkValues(typeParams, patch);
      if (isFailure(checked)) {
        return ApiResult.fail(checked.code, `${provider.label}: ${checked.message}`, checked.hint);
      }

      const name = volumeLabel(ctx, v);
      const existing = v.representations.filter((r) => typeOf(r) === type);
      if (existing.length) {
        const changing = existing.filter((r) => {
          const old = (r.cell.transform.params as RepresentationValues).type.params;
          return JSON.stringify(mergeParams(old, patch)) !== JSON.stringify(old);
        });
        if (!changing.length) return ApiResult.already(`${name} already shows ${provider.label} that way.`);
        await ctx.plugin.dataTransaction(async () => {
          for (const r of changing) {
            const old = r.cell.transform.params as RepresentationValues;
            await ctx.plugin.state.updateTransform(ctx.plugin.state.data, r.cell.transform.ref, {
              ...old,
              type: { name: old.type.name, params: mergeParams(old.type.params, patch) },
            });
          }
        });
        return ApiResult.ok(`Updated ${provider.label} on ${name}.`, {
          nodes: changing.map((r) => r.cell.transform.ref),
        });
      }

      const update = ctx.plugin.state.data
        .build()
        .to(v.cell)
        .apply(
          VolumeRepresentation3D,
          createVolumeRepresentationParams(ctx.plugin, data, { type, typeParams: patch as never }),
        );
      await update.commit();
      const made = ctx.plugin.state.data.cells.get(update.selector.ref);
      if (!made || made.status === 'error') {
        // The state keeps the failed cell; the run's checkpoint restores the scene.
        return ApiResult.fail(
          'failed',
          `Could not show ${name} as ${provider.label}: ${made?.errorText ?? 'it failed'}.`,
        );
      }
      return ApiResult.ok(`Showed ${name} as ${provider.label}.`, { nodes: [update.selector.ref] });
    },
  });

  // The inspector's form: a representation's params as a whole (type, themes), merged like `representation.update`.
  readonly update = this.operation({
    id: 'volume.update',
    title: 'Update a volume representation',
    params: {
      representation: PD.Text('', { description: 'A volume representation by state ref' }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The representation's params (type, color theme, size theme)" },
      ),
    },
    surface: { kind: 'ui-only', reason: "The inspector's form; tools use representation and isovalue" },
    async run({ representation, params }, ctx) {
      const rep = allRepresentations(ctx).find((r) => r.cell.transform.ref === representation);
      if (!rep) return ApiResult.fail('not-found', `There is no volume representation '${representation}'.`);
      const current = rep.cell.transform.params as Record<string, unknown>;
      const next = mergeParams(current, params);
      if (JSON.stringify(next) === JSON.stringify(current)) {
        return ApiResult.already('The representation is already so.');
      }
      await ctx.plugin.state.updateTransform(ctx.plugin.state.data, representation, next, 'Volume representation');
      return ApiResult.ok('Updated the representation.', { nodes: [representation] });
    },
  });

  readonly visibility = this.operation({
    id: 'volume.visibility',
    title: 'Show or hide a volume',
    params: {
      visible: PD.Boolean(true),
      volume: VolumeParam(),
      representation: PD.Text('', {
        description: 'One representation (by state ref or type) instead of the whole volume',
      }),
    },
    async run({ visible, volume, representation }, ctx) {
      const t = representation ? resolveRepresentation(ctx, volume, representation, false) : undefined;
      if (isFailure(t)) return t;
      const v = t ? undefined : oneVolume(ctx, volume || undefined);
      if (isFailure(v)) return v;
      const target = t ? t.rep : v!;
      const what = t ? `${repLabel(t.rep)} of ${volumeLabel(ctx, t.volume)}` : volumeLabel(ctx, v!);
      if (!!target.cell.state.isHidden === !visible) {
        return ApiResult.already(`${what} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      ctx.plugin.managers.volume.hierarchy.toggleVisibility([target], visible ? 'show' : 'hide');
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${what}.`, { nodes: [target.cell.transform.ref] });
    },
  });

  readonly remove = this.operation({
    id: 'volume.remove',
    title: 'Remove a volume',
    params: {
      volume: VolumeParam(),
      representation: PD.Text('', {
        description: 'One representation (by state ref or type) instead of the whole volume',
      }),
    },
    approval: 'destructive',
    async run({ volume, representation }, ctx) {
      if (representation) {
        const t = resolveRepresentation(ctx, volume, representation, false);
        if (isFailure(t)) return t;
        await ctx.plugin.managers.volume.hierarchy.remove([t.rep], false);
        return ApiResult.ok(`Removed ${repLabel(t.rep)} from ${volumeLabel(ctx, t.volume)}.`, {
          nodes: [t.rep.cell.transform.ref],
        });
      }
      const v = oneVolume(ctx, volume || undefined);
      if (isFailure(v)) return v;
      const name = volumeLabel(ctx, v);
      // The whole entry goes: the format and data it was read from too, unless something else hangs off them.
      await ctx.plugin.managers.volume.hierarchy.remove([entryRoot(ctx, v.cell.transform.ref)], false);
      return ApiResult.ok(`Removed ${name}.`, { nodes: [v.cell.transform.ref] });
    },
  });

  // operations.md §2: density streamed from a server for a structure (v5's "Volume Streaming" panel, `Enable`).
  readonly streaming = this.operation({
    id: 'volume.streaming',
    title: 'Stream density for a structure',
    params: {
      enabled: PD.Boolean(true, { description: 'Turn volume streaming on or off for the structure' }),
      structure: PD.Text('', { description: 'A structure by label or state ref; empty: the one in view' }),
    },
    modules: ['network'],
    approval: 'network',
    async run({ enabled, structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const name = labelOf(ctx, s);
      const state = ctx.plugin.state.data;
      const root = StateSelection.findTagInSubtree(state.tree, s.cell.transform.ref, VolumeStreaming.RootTag);

      if (!enabled) {
        if (!root) return ApiResult.already(`Volume streaming is already off for ${name}.`);
        await ctx.plugin.managers.structure.hierarchy.remove([root], false);
        return ApiResult.ok(`Turned volume streaming off for ${name}.`, { nodes: [s.cell.transform.ref] });
      }

      if (root && state.cells.get(root)?.status === 'ok') {
        return ApiResult.already(`Volume streaming is already on for ${name}.`);
      }
      const cell = s.cell;
      if (!cell.obj || !InitVolumeStreaming.definition.isApplicable?.(cell.obj, cell.transform, ctx.plugin)) {
        return ApiResult.fail(
          'blocked',
          `${name} has no density to stream.`,
          'Streaming needs an X-ray or EM structure with a map on the server.',
        );
      }
      if (root) await ctx.plugin.managers.structure.hierarchy.remove([root], false);
      const defs = InitVolumeStreaming.definition.params?.(cell.obj, ctx.plugin) as PD.Params | undefined;
      await ctx.plugin.runTask(
        state.applyAction(InitVolumeStreaming, (defs ? PD.getDefaultValues(defs) : {}) as never, cell.transform.ref),
      );
      const made = StateSelection.findTagInSubtree(state.tree, cell.transform.ref, VolumeStreaming.RootTag);
      const madeCell = made ? state.cells.get(made) : undefined;
      if (
        !made ||
        !madeCell ||
        madeCell.status === 'error' ||
        !madeCell.obj ||
        madeCell.obj.data.entries.length === 0
      ) {
        const why = madeCell?.status === 'error' ? madeCell.errorText : 'no density entry for it on the server';
        return ApiResult.fail(
          'failed',
          `Could not turn volume streaming on for ${name}: ${why ?? 'it failed'}.`,
          'Check the network and that the entry has a map.',
        );
      }
      return ApiResult.ok(`Turned volume streaming on for ${name}.`, { nodes: [s.cell.transform.ref] });
    },
  });
}

// ---- Volumes of the hierarchy ----

type RepresentationValues = { type: { name: string; params: Record<string, any> } } & Record<string, any>;

function volumes(ctx: ApiContext): readonly VolumeRef[] {
  return ctx.plugin.managers.volume.hierarchy.current.volumes;
}

function allRepresentations(ctx: ApiContext): VolumeRepresentationRef[] {
  return volumes(ctx).flatMap((v) => v.representations);
}

/** A volume's name for people: its label, with its state ref added when another loaded one has the same. */
export function volumeNames(list: readonly VolumeRef[]): ReadonlyMap<string, string> {
  return namesOf(list);
}

/** Names for what is in a hierarchy list (volumes, particle lists): the label, told apart by the ref when it repeats. */
export function namesOf(list: readonly { cell: { obj?: { label?: string }; transform: { ref: string } } }[]) {
  const counts = new Map<string, number>();
  for (const x of list) counts.set(x.cell.obj?.label ?? '', (counts.get(x.cell.obj?.label ?? '') ?? 0) + 1);
  const names = new Map<string, string>();
  for (const x of list) {
    const label = x.cell.obj?.label || x.cell.transform.ref;
    names.set(
      x.cell.transform.ref,
      (counts.get(x.cell.obj?.label ?? '') ?? 0) > 1 ? `${label} · ${x.cell.transform.ref}` : label,
    );
  }
  return names;
}

function volumeLabel(ctx: ApiContext, v: VolumeRef): string {
  return volumeNames(volumes(ctx)).get(v.cell.transform.ref) ?? v.cell.transform.ref;
}

/**
 * One of a hierarchy list by a state ref or a label (case-insensitive), else the only one; none or several is a failure
 * that names what there is. `what` is "volume" or "particle list".
 */
export function pickOne<T extends { cell: { obj?: { label?: string }; transform: { ref: string } } }>(
  list: readonly T[],
  name: string | undefined,
  what: string,
  empty: string,
): T | ApiFailure {
  const names = namesOf(list);
  const all = list.map((x) => names.get(x.cell.transform.ref)).join(', ');
  if (!list.length) return ApiResult.fail('not-found', `No ${what}s are loaded.`, empty);
  if (name) {
    const lower = name.trim().toLowerCase();
    const found =
      list.find((x) => x.cell.transform.ref === name) ??
      list.find((x) => names.get(x.cell.transform.ref)?.toLowerCase() === lower) ??
      list.find((x) => (x.cell.obj?.label ?? '').toLowerCase() === lower);
    return found ?? ApiResult.fail('not-found', `There is no ${what} '${name}'.`, `Loaded: ${all}.`);
  }
  if (list.length === 1) return list[0];
  return ApiResult.fail('ambiguous', `Which ${what}? ${list.length} are loaded.`, `Name one: ${all}.`);
}

function oneVolume(ctx: ApiContext, name: string | undefined): VolumeRef | ApiFailure {
  return pickOne(volumes(ctx), name, 'volume', 'Open a map file (ccp4, mrc, map, cube, dx, …) first.');
}

/** A representation by state ref, else by type, among a volume's (or every volume's when `volume` is empty and a representation is named). */
function resolveRepresentation(
  ctx: ApiContext,
  volumeName: string,
  name: string,
  defaultToIso: boolean,
): { rep: VolumeRepresentationRef; volume: VolumeRef } | ApiFailure {
  let pool: VolumeRepresentationRef[];
  if (volumeName || !name) {
    const v = oneVolume(ctx, volumeName || undefined);
    if (isFailure(v)) return v;
    pool = v.representations;
  } else {
    pool = allRepresentations(ctx);
  }
  if (name) {
    const byRef = pool.find((r) => r.cell.transform.ref === name);
    const byType = pool.filter((r) => typeOf(r) === name);
    if (byRef) return { rep: byRef, volume: byRef.volume };
    if (byType.length === 1) return { rep: byType[0], volume: byType[0].volume };
    if (byType.length > 1) {
      return ApiResult.fail(
        'ambiguous',
        `There are ${byType.length} ${name} representations.`,
        `Name one by state ref: ${byType.map((r) => r.cell.transform.ref).join(', ')}.`,
      );
    }
    return ApiResult.fail(
      'not-found',
      `There is no volume representation '${name}'.`,
      pool.length
        ? `Representations: ${pool.map((r) => `${typeOf(r)} (${r.cell.transform.ref})`).join(', ')}.`
        : 'Add one with volume.representation.',
    );
  }
  // The volume's first isosurface; else any with an iso value.
  const v = pool[0]?.volume ?? oneVolume(ctx, undefined);
  if (isFailure(v)) return v;
  const rep =
    (defaultToIso ? (pool.find((r) => typeOf(r) === 'isosurface') ?? pool.find(hasIso)) : pool[0]) ?? undefined;
  if (!rep) {
    return ApiResult.fail(
      'not-found',
      `${volumeLabel(ctx, v)} has no ${defaultToIso ? 'isosurface' : 'representation'}.`,
      'Add one with volume.representation.',
    );
  }
  return { rep, volume: rep.volume };
}

function typeOf(r: VolumeRepresentationRef): string {
  return (r.cell.transform.params as { type?: { name?: string } } | undefined)?.type?.name ?? '';
}

function isoOf(r: VolumeRepresentationRef): Volume.IsoValue | undefined {
  const iso = (r.cell.transform.params as RepresentationValues | undefined)?.type?.params?.isoValue;
  return iso && typeof iso === 'object' && (iso.kind === 'relative' || iso.kind === 'absolute') ? iso : undefined;
}

const hasIso = (r: VolumeRepresentationRef) => !!isoOf(r);

/** A representation's label for people: its type's label. */
function repLabel(r: VolumeRepresentationRef): string {
  return r.cell.obj?.label ?? typeOf(r);
}

// ---- Facts ----

export function statsOf(v: VolumeRef): VolumeStats {
  const data = v.cell.obj?.data;
  return data ? volumeStats(data) : { grid: [0, 0, 0], voxelSize: [0, 0, 0], min: 0, max: 0, mean: 0, sigma: 0 };
}

/** A volume's grid size, voxel size and value statistics. */
export function volumeStats(volume: Volume): VolumeStats {
  const { grid } = volume;
  const dims = grid.cells.space.dimensions;
  const m = Grid.getGridToCartesianTransform(grid);
  const size = (i: number) => Vec3.magnitude(Vec3.create(m[i], m[i + 1], m[i + 2]));
  return {
    grid: [dims[0], dims[1], dims[2]],
    voxelSize: [size(0), size(4), size(8)],
    ...grid.stats,
  };
}

/** A value to three significant digits. */
export const shortNumber = (x: number) => String(Number(x.toPrecision(3)));
const num = shortNumber;

/** "12 × 12 × 12 grid, 1.5 Å voxel, mean 0.1, σ 0.2, range 0 … 1". */
export function volumeFacts(s: VolumeStats): string {
  return `${gridText(s)} grid, ${voxelText(s)} voxel, mean ${num(s.mean)}, σ ${num(s.sigma)}, range ${num(s.min)} … ${num(s.max)}`;
}

export function gridText(s: VolumeStats): string {
  return s.grid.join(' × ');
}

export function voxelText(s: VolumeStats): string {
  const [x, y, z] = s.voxelSize;
  const same = Math.abs(x - y) < 1e-3 * x && Math.abs(x - z) < 1e-3 * x;
  return same ? `${num(x)} Å` : `${s.voxelSize.map(num).join(' × ')} Å`;
}

function formatIso(i: NonNullable<VolumeRepresentationInfo['isoValue']>): string {
  return i.kind === 'relative' ? `${i.relative.toFixed(2)} σ` : num(i.absolute);
}

function describeVolume(ctx: ApiContext, v: VolumeRef, label: string): VolumeInfo {
  const stats = statsOf(v);
  return {
    ref: v.cell.transform.ref,
    label,
    hidden: !!v.cell.state.isHidden,
    source: v.cell.obj?.data.sourceData.kind ?? '',
    stats,
    representations: v.representations.map((r) => {
      const iso = isoOf(r);
      const abs = iso && Volume.IsoValue.toAbsolute(iso, stats).absoluteValue;
      const rel = iso && Volume.IsoValue.toRelative(iso, stats).relativeValue;
      return {
        ref: r.cell.transform.ref,
        type: typeOf(r),
        label: repLabel(r),
        hidden: !!r.cell.state.isHidden,
        ...(iso ? { isoValue: { kind: iso.kind, relative: rel!, absolute: abs! } } : {}),
      };
    }),
  };
}

function sameIso(a: Volume.IsoValue, b: Volume.IsoValue): boolean {
  if (a.kind !== b.kind) return false;
  const x = a.kind === 'relative' ? a.relativeValue : a.absoluteValue;
  const y =
    b.kind === 'relative'
      ? (b as Volume.IsoValue.Relative).relativeValue
      : (b as Volume.IsoValue.Absolute).absoluteValue;
  return Math.abs(x - y) <= 1e-9 * Math.max(1, Math.abs(x));
}

/** The iso value must lie within the map's values. */
function checkRange(s: VolumeStats, kind: 'relative' | 'absolute', value: number): ApiFailure | undefined {
  if (kind === 'absolute') {
    if (value >= s.min && value <= s.max) return undefined;
    return ApiResult.fail(
      'invalid',
      `${num(value)} is outside the map's values.`,
      `They range from ${num(s.min)} to ${num(s.max)}.`,
    );
  }
  if (!(s.sigma > 0)) return undefined;
  const lo = (s.min - s.mean) / s.sigma;
  const hi = (s.max - s.mean) / s.sigma;
  if (value >= lo && value <= hi) return undefined;
  return ApiResult.fail(
    'invalid',
    `${num(value)} σ is outside the map's values.`,
    `They range from ${num(lo)} to ${num(hi)} σ.`,
  );
}

// ---- State ----

/**
 * Params patched: nested groups merge one level, values with a `kind` (an iso value) replace as a whole, so a relative
 * one does not keep an absolute one's field.
 */
export function mergeParams(base: Record<string, any>, patch: Record<string, any>): Record<string, any> {
  const out: Record<string, any> = { ...base };
  for (const [k, v] of Object.entries(patch)) {
    out[k] = isPlain(v) && isPlain(out[k]) && !('kind' in v) ? { ...out[k], ...v } : v;
  }
  return out;
}

const isPlain = (x: unknown): x is Record<string, any> =>
  !!x && typeof x === 'object' && !Array.isArray(x) && Object.getPrototypeOf(x) === Object.prototype;

/** The highest ancestor that holds only this node (the format and data it was read from), below the root. */
export function entryRoot(ctx: ApiContext, ref: string): string {
  const tree = ctx.plugin.state.data.tree;
  let current = ref;
  for (;;) {
    const parent = tree.transforms.get(current)?.parent;
    if (!parent || parent === StateTransform.RootRef || tree.children.get(parent).size !== 1) return current;
    current = parent;
  }
}
