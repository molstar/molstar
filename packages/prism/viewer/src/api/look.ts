/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureComponentManager } from '@molstar/plugin/state/manager/structure/component';
import type { StructureRepresentationRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { type ApiContext, type ApiFailure, ApiArea, ApiResult } from './area.js';
import { selectNodes } from './nodes.js';
import { NodeParam, OptionalParam } from './params.js';
import { componentName, isFailure, labelOf, loadedStructures } from './scope.js';

// `look` (api.md §2, "Look"): how the representations are drawn (hydrogens, quality, light, material, solid surface,
// interior, animation). The scene defaults are the component manager's options: every representation follows them
// until it is given a value of its own, an override, which a later change of the default keeps
// (`StructureComponentManager.setOptions`). `set` sets the defaults (no node) or the representations under a node;
// `reset` takes overrides back to the defaults; `get` reads the defaults and the overrides.

/** The look's values by their names here, and the component manager's names for them. */
export const LookKeys = {
  hydrogens: 'hydrogens',
  quality: 'visualQuality',
  ignoreLight: 'ignoreLight',
  material: 'materialStyle',
  solidSurface: 'solidSurface',
  interior: 'interior',
  animation: 'animation',
} as const satisfies Record<string, StructureComponentManager.ReprOptionKey>;
export type LookKey = keyof typeof LookKeys;
export const LookKeyNames = Object.keys(LookKeys) as LookKey[];

/** The defaults only the scene has (not a representation's look): the clip objects, the interactions' settings. */
const SceneOnlyKeys = { clip: 'clipObjects', interactions: 'interactions' } as const;
type SceneOnlyKey = keyof typeof SceneOnlyKeys;

/** What each look value is called in a result line. */
export const LookLabels: Readonly<Record<LookKey | SceneOnlyKey, string>> = {
  hydrogens: 'hydrogens',
  quality: 'quality',
  ignoreLight: 'ignore light',
  material: 'material',
  solidSurface: 'solid surface',
  interior: 'interior',
  animation: 'animation',
  clip: 'clip objects',
  interactions: 'interactions',
};

/** The look's values, by their names here. */
export type LookValues = { [K in LookKey]: StructureComponentManager.Options[(typeof LookKeys)[K]] };

/** A representation with values of its own, as `get` reads it. */
export interface LookOverride {
  readonly ref: string;
  /** "Cartoon of Polymer (1CBS)". */
  readonly label: string;
  /** Its own values, those other than the defaults. */
  readonly values: Partial<LookValues>;
}

export interface LookState {
  readonly defaults: LookValues;
  readonly overrides: readonly LookOverride[];
}

const P = StructureComponentManager.OptionsParams;

/** The params of `set`: each value optional (absent: kept). */
function lookParams() {
  return {
    node: NodeParam(
      'The representations under these nodes (structure, component[label=Ligand], representation[type=cartoon], #ref); empty: the scene defaults',
    ),
    hydrogens: OptionalParam(P.hydrogens),
    quality: OptionalParam(P.visualQuality),
    ignoreLight: OptionalParam(P.ignoreLight),
    material: OptionalParam(P.materialStyle),
    solidSurface: OptionalParam(P.solidSurface),
    interior: OptionalParam(P.interior),
    animation: OptionalParam(P.animation),
    clip: OptionalParam({ ...P.clipObjects, description: 'The clip objects of every representation (scene only)' }),
    interactions: OptionalParam({ ...P.interactions, description: 'How interactions are computed (scene only)' }),
  };
}

export class LookApi extends ApiArea {
  readonly id = 'look';

  readonly set = this.operation({
    id: 'look.set',
    // the Structure panel's options before it, which set only the defaults
    aliases: ['component.options'],
    title: 'Set how representations are drawn',
    params: lookParams(),
    async run(values, ctx) {
      const { node, ...given } = values as Record<string, unknown> & { node: string };
      const manager = ctx.plugin.managers.structure.component;
      if (!node) {
        const current = manager.state.options;
        const next: StructureComponentManager.Options = { ...current };
        const changed: (LookKey | SceneOnlyKey)[] = [];
        for (const [name, key] of Object.entries({ ...LookKeys, ...SceneOnlyKeys }) as [
          LookKey | SceneOnlyKey,
          StructureComponentManager.OptionKey,
        ][]) {
          const value = given[name];
          if (value === undefined || PD.isParamEqual(P[key] as PD.Any, current[key], value)) continue;
          (next as Record<string, unknown>)[key] = value;
          changed.push(name);
        }
        if (!changed.length) return ApiResult.already('Those defaults are already set.');
        await manager.setOptions(next);
        const kept = changed.filter((k): k is LookKey => k in LookKeys);
        const own = overridesOf(ctx, everyRepresentation(ctx), kept);
        return ApiResult.ok(
          `Set the default ${listText(changed.map((k) => LookLabels[k]))}${own.length ? `; ${keptText(own)}` : ''}.`,
          { nodes: own.map((o) => o.ref) },
          { value: changed },
        );
      }
      const scene = (Object.keys(SceneOnlyKeys) as SceneOnlyKey[]).filter((k) => given[k] !== undefined);
      if (scene.length) {
        return ApiResult.fail(
          'invalid',
          `${listText(scene.map((k) => LookLabels[k]))} ${scene.length > 1 ? 'are' : 'is'} the scene's, not a node's.`,
          'Leave node empty to set it for the scene; clipping.cut cuts a node.',
        );
      }
      const reprs = representationsUnder(ctx, node);
      if (isFailure(reprs)) return reprs;
      const keys = LookKeyNames.filter((k) => given[k] !== undefined);
      if (!keys.length) return ApiResult.already('Nothing to set: give a value (hydrogens, quality, material, …).');
      const options = { ...manager.state.options } as StructureComponentManager.Options;
      for (const k of keys) (options as Record<string, unknown>)[LookKeys[k]] = given[k];
      const done = await writeLook(ctx, reprs, keys, options);
      if (!done.length)
        return ApiResult.already(`${reprsText(ctx, reprs)} already ${reprs.length > 1 ? 'have' : 'has'} that.`);
      return ApiResult.ok(`Set ${listText(keys.map((k) => LookLabels[k]))} of ${reprsText(ctx, done)}.`, {
        nodes: done.map(refOf),
      });
    },
  });

  readonly reset = this.operation({
    id: 'look.reset',
    title: 'Use the scene defaults again',
    params: {
      node: NodeParam('The representations under these nodes; empty: every representation'),
      keys: PD.MultiSelect<LookKey>(
        [],
        LookKeyNames.map((k) => [k, LookLabels[k]] as [LookKey, string]),
        { description: 'Which values; empty: all of them' },
      ),
    },
    async run({ node, keys }, ctx) {
      const which = keys.length ? keys : LookKeyNames;
      const reprs = node ? representationsUnder(ctx, node) : everyRepresentation(ctx);
      if (isFailure(reprs)) return reprs;
      const own = overridesOf(ctx, reprs, which);
      if (!own.length) {
        return ApiResult.already(
          `${node ? reprsText(ctx, reprs) : 'Every representation'} already ${reprs.length > 1 || !node ? 'use' : 'uses'} the scene defaults.`,
        );
      }
      const manager = ctx.plugin.managers.structure.component;
      const byRef = new Set(own.map((o) => o.ref));
      const done = await writeLook(
        ctx,
        reprs.filter((r) => byRef.has(refOf(r))),
        which,
        manager.state.options,
      );
      return ApiResult.ok(
        `${reprsText(ctx, done)} ${done.length > 1 ? 'use' : 'uses'} the scene's ${listText(which.filter((k) => own.some((o) => k in o.values)).map((k) => LookLabels[k]))} again.`,
        { nodes: done.map(refOf) },
      );
    },
  });

  readonly get = this.operation({
    id: 'look.get',
    title: 'Read how representations are drawn',
    params: { node: NodeParam('Only the overrides under these nodes; empty: all of them') },
    read: true,
    async run({ node }, ctx) {
      const reprs = node ? representationsUnder(ctx, node) : everyRepresentation(ctx);
      if (isFailure(reprs)) return reprs;
      const state: LookState = { defaults: lookDefaults(ctx), overrides: overridesOf(ctx, reprs, LookKeyNames) };
      const d = state.defaults;
      const lines = [
        `Defaults: hydrogens ${d.hydrogens}, quality ${d.quality}, ignore light ${d.ignoreLight ? 'on' : 'off'}, material ${materialText(d.material)}, solid surface ${d.solidSurface ? 'on' : 'off'}.`,
        ...state.overrides.map(
          (o) =>
            `${o.label} (${o.ref}) has its own ${(Object.keys(o.values) as LookKey[]).map((k) => `${LookLabels[k]} ${valueText(k, o.values[k])}`).join(', ')}.`,
        ),
      ];
      if (!state.overrides.length) lines.push('No representation has values of its own.');
      return ApiResult.ok(lines.join('\n'), {}, { value: state });
    },
  });
}

/** The scene defaults, by the look's names. */
export function lookDefaults(ctx: Pick<ApiContext, 'plugin'>): LookValues {
  const options = ctx.plugin.managers.structure.component.state.options;
  const out = {} as Record<LookKey, unknown>;
  for (const k of LookKeyNames) out[k] = options[LookKeys[k]];
  return out as LookValues;
}

/** A representation's values other than the defaults (of the keys asked), by the look's names. */
export function ownLookValues(
  ctx: Pick<ApiContext, 'plugin'>,
  r: StructureRepresentationRef,
  keys: readonly LookKey[] = LookKeyNames,
): Partial<LookValues> {
  const options = ctx.plugin.managers.structure.component.state.options;
  const { type } = r.cell.transform.params as { type: { name: string; params: Record<string, unknown> } };
  const out: Partial<Record<LookKey, unknown>> = {};
  for (const k of keys) {
    if (StructureComponentManager.reprHasOption(type.name, type.params, LookKeys[k], options) !== false) continue;
    out[k] = StructureComponentManager.reprOptionValue(type.name, type.params, LookKeys[k]);
  }
  return out as Partial<LookValues>;
}

/** The scene's structure representations (those the look applies to). */
export function everyRepresentation(ctx: Pick<ApiContext, 'plugin'>): StructureRepresentationRef[] {
  const out: StructureRepresentationRef[] = [];
  for (const s of loadedStructures(ctx)) {
    for (const c of s.components) {
      for (const r of c.representations) if (r.cell.transform.transformer === StructureRepresentation3D) out.push(r);
    }
  }
  return out;
}

/** The structure representations at or under the nodes a selector matches. */
export function representationsUnder(
  ctx: Pick<ApiContext, 'plugin'>,
  node: string,
): StructureRepresentationRef[] | ApiFailure {
  const state = ctx.plugin.state.data;
  const cells = selectNodes(state, node);
  if (isFailure(cells)) return cells;
  const refs = new Set(cells.map((c: StateObjectCell) => c.transform.ref));
  const under = (ref: string) => {
    for (let t = state.tree.transforms.get(ref); t; t = state.tree.transforms.get(t.parent)) {
      if (refs.has(t.ref)) return true;
      if (t.ref === t.parent) break;
    }
    return false;
  };
  const out = everyRepresentation(ctx).filter((r) => under(refOf(r)));
  if (!out.length) {
    return ApiResult.fail(
      'not-found',
      `No structure representation is at or under '${node}'.`,
      'A structure, a component or a representation (state.get lists them).',
    );
  }
  return out;
}

/** The representations among these with values of their own (of the keys asked). */
function overridesOf(
  ctx: Pick<ApiContext, 'plugin'>,
  reprs: readonly StructureRepresentationRef[],
  keys: readonly LookKey[],
): LookOverride[] {
  const out: LookOverride[] = [];
  for (const r of reprs) {
    const values = ownLookValues(ctx, r, keys);
    if (Object.keys(values).length) out.push({ ref: refOf(r), label: reprLabel(ctx, r), values });
  }
  return out;
}

/** Writes the keys' values of `options` to the representations that differ; returns those written. */
async function writeLook(
  ctx: ApiContext,
  reprs: readonly StructureRepresentationRef[],
  keys: readonly LookKey[],
  options: StructureComponentManager.Options,
): Promise<StructureRepresentationRef[]> {
  const update = ctx.plugin.state.data.build();
  const done: StructureRepresentationRef[] = [];
  for (const r of reprs) {
    const { type } = r.cell.transform.params as { type: { name: string; params: Record<string, unknown> } };
    const todo = keys.filter(
      (k) => StructureComponentManager.reprHasOption(type.name, type.params, LookKeys[k], options) === false,
    );
    if (!todo.length) continue;
    done.push(r);
    update.to(r.cell).update(StructureRepresentation3D, (old) => {
      for (const k of todo) {
        StructureComponentManager.applyReprOption(old.type.name, old.type.params, LookKeys[k], options);
      }
    });
  }
  if (done.length) await update.commit();
  return done;
}

const refOf = (r: StructureRepresentationRef) => r.cell.transform.ref;

function reprLabel(ctx: Pick<ApiContext, 'plugin'>, r: StructureRepresentationRef): string {
  const type = (r.cell.transform.params as { type: { name: string } }).type.name;
  const name = ctx.plugin.representation.structure.registry.get(type)?.label ?? type;
  return `${name} of ${componentName(r.component)} (${labelOf(ctx, r.component.structure)})`;
}

/** "cartoon of Polymer (1CBS)", "3 representations". */
function reprsText(ctx: Pick<ApiContext, 'plugin'>, reprs: readonly StructureRepresentationRef[]): string {
  return reprs.length === 1 ? reprLabel(ctx, reprs[0]) : `${reprs.length} representations`;
}

/** "Ball & Stick of Ligand (1CBS) keeps its own", "2 representations keep their own". */
function keptText(own: readonly LookOverride[]): string {
  return own.length === 1 ? `${own[0].label} keeps its own` : `${own.length} representations keep their own`;
}

function listText(items: readonly string[]): string {
  if (items.length <= 1) return items.join('');
  return `${items.slice(0, -1).join(', ')} and ${items[items.length - 1]}`;
}

/** A material as its preset's name (Matte, Glossy), else its numbers. */
export function materialText(m: { metalness: number; roughness: number; bumpiness: number }): string {
  const presets = (P.materialStyle.presets ?? []) as readonly (readonly [typeof m, string])[];
  const preset = presets.find(
    ([p]) => p.metalness === m.metalness && p.roughness === m.roughness && p.bumpiness === m.bumpiness,
  );
  return preset ? preset[1] : `metalness ${m.metalness}, roughness ${m.roughness}, bumpiness ${m.bumpiness}`;
}

/** A look value for a line. */
export function valueText(key: LookKey, value: unknown): string {
  if (key === 'material') return materialText(value as LookValues['material']);
  if (typeof value === 'boolean') return value ? 'on' : 'off';
  if (key === 'hydrogens' || key === 'quality') {
    const options = (key === 'hydrogens' ? P.hydrogens : P.visualQuality).options as readonly (readonly [
      string,
      string,
    ])[];
    return options.find(([v]) => v === value)?.[1] ?? String(value);
  }
  return 'custom';
}
