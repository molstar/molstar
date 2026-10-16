/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell } from '@molstar/core/state';
import { ColorLists } from '@molstar/core/util/color/lists';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { PluginContext } from '@molstar/plugin/context';

// The data palette (README §3.7, decision 31): the categorical color themes (chains, entities, polymers, operators,
// models) take the chosen list as their default palette, so every new representation and preset uses it; the
// representations that use a theme's default palette now are recolored with the new one. A theme's own non-default
// palette (chosen by the person) is left alone. The registry keeps the themes' names: saved sessions read the same.
// Prism's own list is the default choice; Mol*'s classic one is the themes' own; the lists are the spec's
// (`spec.palettes`), else Many-Distinct (Prism) and Okabe–Ito.

/** Prism's own list (many-distinct-prism), Mol*'s (the themes' own: many-distinct), or the colorblind-safe one. */
export type DataPalette = 'prism' | 'classic' | 'colorblind';

/** The lists the palettes take (0xRRGGBB), in the order categories take them; a choice not given keeps the default. */
export interface DataPaletteLists {
  readonly prism?: readonly number[];
  readonly colorblind?: readonly number[];
}

/** The themes whose default palette follows the choice. */
export const CategoryThemes: readonly string[] = [
  'chain-id',
  'entity-id',
  'entity-source',
  'polymer-id',
  'polymer-index',
  'operator-name',
  'operator-hkl',
  'model-index',
  'structure-index',
  'unit-index',
  'trajectory-index',
];

/** A spec's lists (`spec.palettes`, §2.8): `'#rrggbb'` strings read as colors. */
export function paletteLists(
  spec: { readonly prism?: readonly (number | string)[]; readonly colorblind?: readonly (number | string)[] } = {},
): DataPaletteLists {
  const read = (list?: readonly (number | string)[]) =>
    list?.map((c) => (typeof c === 'number' ? c : Number.parseInt(c.replace(/^#/, ''), 16)));
  return { prism: read(spec.prism), colorblind: read(spec.colorblind) };
}

/** A list as a palette param's value (a set). */
const paletteOf = (colors: readonly number[]) => ({
  name: 'colors',
  params: { list: { kind: 'set' as const, colors: [...colors] } },
});

interface PaletteState {
  palette: DataPalette;
  /** The palette param values the lists give. */
  readonly values: { readonly prism: unknown; readonly colorblind: unknown };
  /** The registry's own providers, by name, while ours stand in. */
  readonly originals: Map<string, ColorTheme.Provider<any, any>>;
}

const states = new WeakMap<PluginContext, PaletteState>();

export function dataPaletteOf(plugin: PluginContext): DataPalette {
  return states.get(plugin)?.palette ?? 'classic';
}

/** A palette a session names: an older one's `default` is the classic one; none (a session from before) too. */
function readPalette(value: unknown): DataPalette {
  return value === 'prism' || value === 'colorblind' ? value : 'classic';
}

/** A provider like `p` whose `palette` param defaults to `value`. */
function withPalette(p: ColorTheme.Provider<any, any>, value: unknown): ColorTheme.Provider<any, any> {
  return {
    ...p,
    getParams: (ctx) => {
      const params = p.getParams(ctx);
      const palette = params.palette as PD.Any | undefined;
      if (!palette) return params;
      return { ...params, palette: { ...palette, defaultValue: value } };
    },
  };
}

/** Whether a representation's color theme uses its theme's default palette (so it follows the choice). */
function usesDefault(params: Record<string, any> | undefined, palette: unknown): boolean {
  if (!params?.palette) return true;
  return JSON.stringify(params.palette) === JSON.stringify(palette);
}

/** The registry's categorical themes under a palette (no recolor); their defaults before, by name. */
function swapThemes(plugin: PluginContext, state: PaletteState, palette: DataPalette): Map<string, unknown> {
  const registry = plugin.representation.structure.themes.colorThemeRegistry;
  const before = new Map<string, unknown>();
  for (const name of CategoryThemes) {
    const p = registry.get(name);
    if (!p || p.name !== name) continue;
    before.set(name, (p.getParams(NoData).palette as PD.Any | undefined)?.defaultValue);
    const original = state.originals.get(name) ?? p;
    state.originals.set(name, original);
    registry.remove(p);
    registry.add(palette === 'classic' ? original : withPalette(original, state.values[palette]));
  }
  state.palette = palette;
  return before;
}

const NoData = { structure: undefined } as never;

/**
 * Registers the palette's custom state when a viewer is made, with the lists the palettes take, whatever its starting
 * palette: a session or a checkpoint from a viewer on another palette then brings that palette back.
 */
export function initDataPalette(plugin: PluginContext, lists: DataPaletteLists = {}): void {
  stateOf(plugin, lists);
}

/** The palette's state for a plugin, registered once as its custom state: checkpoints and sessions carry it. */
function stateOf(plugin: PluginContext, lists: DataPaletteLists = {}): PaletteState {
  let state = states.get(plugin);
  if (state) return state;
  const made: PaletteState = {
    palette: 'classic',
    values: {
      prism: paletteOf(lists.prism?.length ? lists.prism : (ColorLists['many-distinct-prism'].list as number[])),
      colorblind: paletteOf(lists.colorblind?.length ? lists.colorblind : (ColorLists['okabe-ito'].list as number[])),
    },
    originals: new Map(),
  };
  state = made;
  states.set(plugin, made);
  plugin.state.registerCustomState('prism-viewer.palette', {
    get: () => ({ version: 1, palette: made.palette }),
    // An undo or a session: the registry as it was then (the representations come back with their own params).
    set: (value) => {
      const palette = readPalette((value as { palette?: unknown } | undefined)?.palette);
      if (palette !== made.palette) swapThemes(plugin, made, palette);
    },
  });
  return state;
}

/**
 * Sets the data palette: the categorical themes' default, and the representations that use the old default recolored
 * (one state update: one history step through the API). The refs recolored; undefined when it was so already.
 */
export async function setDataPalette(plugin: PluginContext, palette: DataPalette): Promise<string[] | undefined> {
  const state = stateOf(plugin);
  if (state.palette === palette) return undefined;
  const before = swapThemes(plugin, state, palette);
  const registry = plugin.representation.structure.themes.colorThemeRegistry;
  // Recolor what used the old default (the person's own palettes stay).
  const update = plugin.build();
  const refs: string[] = [];
  for (const cell of plugin.state.data.cells.values() as Iterable<StateObjectCell>) {
    const params = cell.transform.params as { colorTheme?: { name: string; params?: Record<string, any> } } | undefined;
    const theme = params?.colorTheme;
    if (!theme || !before.has(theme.name) || !usesDefault(theme.params, before.get(theme.name))) continue;
    const next = (registry.get(theme.name).getParams(NoData).palette as PD.Any | undefined)?.defaultValue;
    if (!next) continue;
    update
      .to(cell.transform.ref)
      .update({ ...params, colorTheme: { ...theme, params: { ...theme.params, palette: next } } });
    refs.push(cell.transform.ref);
  }
  if (refs.length) await update.commit();
  return refs;
}
