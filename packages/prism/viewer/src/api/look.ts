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
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { ClipVolumesTag } from '@molstar/plugin/state/transforms/shape/clip-volumes';
import { StructureFocusRepresentationTags } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
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
const DefaultOptions = PD.getDefaultValues(P);
const FocusTags: ReadonlySet<string> = new Set(Object.values(StructureFocusRepresentationTags));

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
      const cells = lookCellsUnder(ctx, node);
      if (isFailure(cells)) return cells;
      const keys = LookKeyNames.filter((k) => given[k] !== undefined);
      if (!keys.length) return ApiResult.already('Nothing to set: give a value (hydrogens, quality, material, …).');
      const taken = keys.filter((k) => cells.some((c) => cellTakes(c, k)));
      if (!taken.length) {
        return ApiResult.fail(
          'invalid',
          `${cellsText(ctx, cells)} ${cells.length > 1 ? 'take' : 'takes'} no ${listText(keys.map((k) => LookLabels[k]))}.`,
          'Hydrogens are a structure representation’s; solid surface a surface’s.',
        );
      }
      const options = { ...manager.state.options } as StructureComponentManager.Options;
      for (const k of keys) (options as Record<string, unknown>)[LookKeys[k]] = given[k];
      const done = await writeLook(ctx, cells, taken, options);
      if (!done.length)
        return ApiResult.already(`${cellsText(ctx, cells)} already ${cells.length > 1 ? 'have' : 'has'} that.`);
      return ApiResult.ok(`Set ${listText(taken.map((k) => LookLabels[k]))} of ${cellsText(ctx, done)}.`, {
        nodes: done.map((c) => c.transform.ref),
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
      const reprs = node ? structureReprsUnder(ctx, node) : everyRepresentation(ctx);
      if (isFailure(reprs)) return reprs;
      // a map's, particles' or a shape's look is its own: no scene default to go back to
      if (node && !reprs.length) {
        const cells = lookCellsUnder(ctx, node);
        return ApiResult.fail(
          'invalid',
          `${isFailure(cells) ? node : cellsText(ctx, cells)} ${!isFailure(cells) && cells.length > 1 ? 'follow' : 'follows'} no scene defaults: only structure representations do.`,
          'Set its values with look.set by node.',
        );
      }
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
        reprs.filter((r) => byRef.has(refOf(r))).map((r) => r.cell),
        which,
        manager.state.options,
      );
      return ApiResult.ok(
        `${cellsText(ctx, done)} ${done.length > 1 ? 'use' : 'uses'} the scene's ${listText(which.filter((k) => own.some((o) => k in o.values)).map((k) => LookLabels[k]))} again.`,
        { nodes: done.map((c) => c.transform.ref) },
      );
    },
  });

  readonly get = this.operation({
    id: 'look.get',
    title: 'Read how representations are drawn',
    params: { node: NodeParam('Only the overrides under these nodes; empty: all of them') },
    read: true,
    async run({ node }, ctx) {
      const reprs = node ? structureReprsUnder(ctx, node) : everyRepresentation(ctx);
      if (isFailure(reprs)) return reprs;
      const state: LookState = { defaults: lookDefaults(ctx), overrides: overridesOf(ctx, reprs, LookKeyNames) };
      const d = state.defaults;
      const lines = [
        `Defaults: ${LookKeyNames.map((k) => `${LookLabels[k]} ${valueText(k, d[k])}`).join(', ')}.`,
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
      for (const r of c.representations) {
        if (r.cell.transform.transformer !== StructureRepresentation3D) continue;
        // the focus behavior's are its own (their look comes through its params)
        if (r.cell.transform.tags?.some((t) => FocusTags.has(t))) continue;
        out.push(r);
      }
    }
  }
  return out;
}

/**
 * The representations at or under the nodes a selector matches that take a look: a structure's, a volume's, particles',
 * a shape's (not a cut's drawing).
 */
export function lookCellsUnder(ctx: Pick<ApiContext, 'plugin'>, node: string): StateObjectCell[] | ApiFailure {
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
  const out: StateObjectCell[] = [];
  state.cells.forEach((cell) => {
    if (cell.state.isGhost || !cell.obj || !PluginStateObject.isRepresentation3D(cell.obj)) return;
    // a cut's drawing is the cut's; the focus behavior's representations are the behavior's (it rewrites them, and
    // takes the scene defaults through its params: `setOptions`)
    const tags = cell.transform.tags;
    if (tags?.some((t) => t === ClipVolumesTag || FocusTags.has(t)) || !lookParamsOf(cell.transform.params)) return;
    if (under(cell.transform.ref)) out.push(cell);
  });
  if (!out.length) {
    return ApiResult.fail(
      'not-found',
      `No representation is at or under '${node}'.`,
      'A structure, a component, a volume, a shape or a representation (state.get lists them).',
    );
  }
  return out;
}

/** The structure representations at or under the nodes a selector matches (those that follow the scene defaults). */
export function structureReprsUnder(
  ctx: Pick<ApiContext, 'plugin'>,
  node: string,
): StructureRepresentationRef[] | ApiFailure {
  const cells = lookCellsUnder(ctx, node);
  if (isFailure(cells)) return cells;
  return structureReprsOf(ctx, cells);
}

/** The structure representations among these cells, as the hierarchy has them. */
export function structureReprsOf(
  ctx: Pick<ApiContext, 'plugin'>,
  cells: readonly StateObjectCell[],
): StructureRepresentationRef[] {
  const refs = ctx.plugin.managers.structure.hierarchy.current.refs;
  const out: StructureRepresentationRef[] = [];
  for (const c of cells) {
    const r = refs.get(c.transform.ref);
    if (r?.kind === 'structure-representation' && c.transform.transformer === StructureRepresentation3D) out.push(r);
  }
  return out;
}

/** Where a representation's look params are: its type's params (structure, volume, particles), or its own (a shape). */
export function lookParamsOf(params: unknown): Record<string, unknown> | undefined {
  const p = params as { type?: { params?: Record<string, unknown> } } & Record<string, unknown>;
  const at = p?.type?.params ?? p;
  return at && typeof at === 'object' && 'material' in at && 'ignoreLight' in at ? at : undefined;
}

/** The params with the look params `at` in place of theirs (the same path). */
function withLookParams(params: unknown, at: Record<string, unknown>): unknown {
  const p = params as { type?: { name: string; params?: Record<string, unknown> } } & Record<string, unknown>;
  return p.type?.params ? { ...p, type: { ...p.type, params: at } } : at;
}

/** The look params each value is written to, outside a structure representation (which the component manager reads). */
const LookParams: Readonly<Record<Exclude<LookKey, 'hydrogens'>, string>> = {
  quality: 'quality',
  ignoreLight: 'ignoreLight',
  material: 'material',
  solidSurface: 'solidInterior',
  interior: 'interior',
  animation: 'animation',
};

const isStructureRepr = (cell: StateObjectCell) => cell.transform.transformer === StructureRepresentation3D;

/** Whether a representation takes a look value (hydrogens: a structure's; solid surface: a surface's). */
export function cellTakes(cell: StateObjectCell, key: LookKey): boolean {
  if (isStructureRepr(cell)) {
    const { type } = cell.transform.params as { type: { name: string; params: Record<string, unknown> } };
    // undefined: the type does not take it (whichever value the options have)
    return StructureComponentManager.reprHasOption(type.name, type.params, LookKeys[key], DefaultOptions) !== undefined;
  }
  if (key === 'hydrogens') return false;
  const at = lookParamsOf(cell.transform.params);
  return !!at && LookParams[key] in at;
}

/** A representation's look value, by the look's names; undefined where it takes none. */
export function cellLookValue(cell: StateObjectCell, key: LookKey): unknown {
  if (isStructureRepr(cell)) {
    const { type } = cell.transform.params as { type: { name: string; params: Record<string, unknown> } };
    return StructureComponentManager.reprOptionValue(type.name, type.params, LookKeys[key]);
  }
  if (key === 'hydrogens') return undefined;
  return lookParamsOf(cell.transform.params)?.[LookParams[key]];
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

/**
 * Writes the keys' values of `options` to the representations that take them and differ, as one update; returns those
 * written. A structure representation through the component manager's rules (`applyReprOption`), the others' look
 * params as they are named there.
 */
async function writeLook(
  ctx: ApiContext,
  cells: readonly StateObjectCell[],
  keys: readonly LookKey[],
  options: StructureComponentManager.Options,
): Promise<StateObjectCell[]> {
  const update = ctx.plugin.state.data.build();
  const done: StateObjectCell[] = [];
  for (const cell of cells) {
    if (isStructureRepr(cell)) {
      const { type } = cell.transform.params as { type: { name: string; params: Record<string, unknown> } };
      const todo = keys.filter(
        (k) => StructureComponentManager.reprHasOption(type.name, type.params, LookKeys[k], options) === false,
      );
      if (!todo.length) continue;
      done.push(cell);
      update.to(cell).update(StructureRepresentation3D, (old) => {
        for (const k of todo) {
          StructureComponentManager.applyReprOption(old.type.name, old.type.params, LookKeys[k], options);
        }
      });
      continue;
    }
    const at = lookParamsOf(cell.transform.params);
    if (!at) continue;
    const next = { ...at };
    let changed = false;
    for (const k of keys) {
      if (k === 'hydrogens' || !(LookParams[k] in at)) continue;
      const key = LookKeys[k];
      if (PD.isParamEqual(P[key] as PD.Any, at[LookParams[k]], options[key])) continue;
      next[LookParams[k]] = options[key];
      changed = true;
    }
    if (!changed) continue;
    done.push(cell);
    update.to(cell).update(withLookParams(cell.transform.params, next) as never);
  }
  if (done.length) await update.commit();
  return done;
}

const refOf = (r: StructureRepresentationRef) => r.cell.transform.ref;

/** "Cartoon of Polymer (1CBS)": a structure representation for a line. */
export function reprLabel(ctx: Pick<ApiContext, 'plugin'>, r: StructureRepresentationRef): string {
  const type = (r.cell.transform.params as { type: { name: string } }).type.name;
  const name = ctx.plugin.representation.structure.registry.get(type)?.label ?? type;
  return `${name} of ${componentName(r.component)} (${labelOf(ctx, r.component.structure)})`;
}

/** A representation for a line: a structure representation as `reprLabel`, another by its label and its parent's. */
function cellLabel(ctx: Pick<ApiContext, 'plugin'>, cell: StateObjectCell): string {
  const [r] = structureReprsOf(ctx, [cell]);
  if (r) return reprLabel(ctx, r);
  const parent = ctx.plugin.state.data.cells.get(cell.transform.parent);
  const label = cell.obj?.label ?? cell.transform.transformer.definition.display.name;
  return parent?.obj?.label ? `${label} of ${parent.obj.label}` : label;
}

/** "Mesh of Shape", "3 representations". */
function cellsText(ctx: Pick<ApiContext, 'plugin'>, cells: readonly StateObjectCell[]): string {
  return cells.length === 1 ? cellLabel(ctx, cells[0]) : `${cells.length} representations`;
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
  if (key === 'interior') {
    const i = value as LookValues['interior'];
    return `color #${i.color.toString(16).padStart(6, '0')} at ${i.colorStrength}`;
  }
  if (key === 'animation') {
    const a = value as LookValues['animation'];
    return `wiggle ${a.wiggleAmplitude}, tumble ${a.tumbleAmplitude}`;
  }
  return String(value);
}
