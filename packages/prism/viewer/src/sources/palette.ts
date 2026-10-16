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
// models) take a colorblind-safe list (Okabe–Ito) as their default palette, so every new representation and preset
// uses it; the representations that use a theme's default palette now are recolored with the new one. A theme's own
// non-default palette (chosen by the person) is left alone. The registry keeps the themes' names: saved sessions read
// the same.

export type DataPalette = 'default' | 'colorblind';

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

/** Okabe–Ito's colors as a palette param's value (a set from the core's list). */
const ColorblindPalette = {
  name: 'colors',
  params: { list: { kind: 'set' as const, colors: [...ColorLists['okabe-ito'].list] } },
};

interface PaletteState {
  palette: DataPalette;
  /** The registry's own providers, by name, while ours stand in. */
  readonly originals: Map<string, ColorTheme.Provider<any, any>>;
}

const states = new WeakMap<PluginContext, PaletteState>();

export function dataPaletteOf(plugin: PluginContext): DataPalette {
  return states.get(plugin)?.palette ?? 'default';
}

/** A provider like `p` whose `palette` param defaults to the colorblind-safe list. */
function withPalette(p: ColorTheme.Provider<any, any>): ColorTheme.Provider<any, any> {
  return {
    ...p,
    getParams: (ctx) => {
      const params = p.getParams(ctx);
      const palette = params.palette as PD.Any | undefined;
      if (!palette) return params;
      return { ...params, palette: { ...palette, defaultValue: ColorblindPalette } };
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
    registry.add(palette === 'colorblind' ? withPalette(original) : original);
  }
  state.palette = palette;
  return before;
}

const NoData = { structure: undefined } as never;

/**
 * Registers the palette's custom state when a viewer is made, whatever its starting palette: a session or a checkpoint
 * from a colorblind viewer then brings the palette back in one that started with the default.
 */
export function initDataPalette(plugin: PluginContext): void {
  stateOf(plugin);
}

/** The palette's state for a plugin, registered once as its custom state: checkpoints and sessions carry it. */
function stateOf(plugin: PluginContext): PaletteState {
  let state = states.get(plugin);
  if (state) return state;
  const made: PaletteState = { palette: 'default', originals: new Map() };
  state = made;
  states.set(plugin, made);
  plugin.state.registerCustomState('prism-viewer.palette', {
    get: () => ({ version: 1, palette: made.palette }),
    // An undo or a session: the registry as it was then (the representations come back with their own params).
    set: (value) => {
      const palette = (value as { palette?: DataPalette } | undefined)?.palette ?? 'default';
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
