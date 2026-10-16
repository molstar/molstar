/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { MmcifFormat } from '@molstar/model/formats/structure/mmcif-format';
import type { Model, Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { observableState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type {
  QuickChoice,
  RepresentationQuickChange,
  RepresentationQuickModel,
} from '../model/representation-quick.js';
import type { ScaleLegend, TableLegend } from '@molstar/core/util/legend';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import { gestureSteps } from './gesture-steps.js';

// A structure representation's quick settings over the API (model/representation-quick.ts): its params rebuilt from the
// state and the change, one `representation.update` per change, a gesture one undo step (`gestureSteps`).

export interface RepresentationQuickSource extends RepresentationQuickModel {
  /** Reads the representation again (the state tree changed). */
  refresh(): void;
}

type Params = {
  type: { name: string; params: { alpha?: number; sizeFactor?: number } & Record<string, unknown> };
  colorTheme: { name: string; params: Record<string, unknown> };
  sizeTheme: { name: string; params: Record<string, unknown> };
};

/** The types offered at once, in order (those that apply to the structure). */
const CommonTypes: readonly [string, string][] = [
  ['cartoon', 'Cartoon'],
  ['ball-and-stick', 'Ball & stick'],
  ['spacefill', 'Spacefill'],
  ['molecular-surface', 'Surface'],
  ['gaussian-surface', 'Gaussian surface'],
];

const TypeIcons: Readonly<Record<string, IconName>> = {
  cartoon: 'cartoon',
  'ball-and-stick': 'ball-and-stick',
  spacefill: 'spacefill',
  'molecular-surface': 'surface',
  'gaussian-surface': 'volume',
};

/** The color themes offered at once, in order (those that apply). */
const CommonThemes: readonly [string, string][] = [
  ['chain-id', 'Chain'],
  ['element-symbol', 'Element'],
  ['secondary-structure', 'Secondary'],
  ['residue-name', 'Residue'],
  ['sequence-id', 'Rainbow'],
  ['plddt-confidence', 'pLDDT'],
  ['uniform', 'One color'],
];

const round = (x: number) => Math.round(x * 1000) / 1000;

/**
 * A predicted model, by positive evidence: ModelCIF's local quality metrics (AlphaFold DB, ModelArchive), or no
 * experimental method at all (a predictor's PDB file, whose B-factors are pLDDT).
 */
function isPredicted(model: Model): boolean {
  const source = model.sourceData;
  if (!MmcifFormat.is(source)) return false;
  if (source.data.frame.categories.ma_qa_metric_local) return true;
  return source.data.db.exptl.method.rowCount === 0;
}
const sameChoices = (a: readonly QuickChoice[], b: readonly QuickChoice[]) =>
  a.length === b.length && a.every((c, i) => c.id === b[i].id && c.swatch === b[i].swatch);

/** A legend as a CSS background: a scale's gradient, a table's first colors as bands. */
function legendSwatch(legend: ScaleLegend | TableLegend | undefined): string | undefined {
  const hex = (c: number) => `#${c.toString(16).padStart(6, '0')}`;
  if (legend?.kind === 'scale-legend') {
    const colors = legend.colors.map((c) => (Array.isArray(c) ? c[0] : c) as number);
    return colors.length ? `linear-gradient(90deg, ${colors.map(hex).join(', ')})` : undefined;
  }
  if (legend?.kind === 'table-legend') {
    const colors = legend.table.slice(0, 6).map(([, c]) => c as number);
    if (!colors.length) return undefined;
    const w = 100 / colors.length;
    const bands = colors.map((c, i) => `${hex(c)} ${(i * w).toFixed(1)}% ${((i + 1) * w).toFixed(1)}%`);
    return `linear-gradient(90deg, ${bands.join(', ')})`;
  }
  return undefined;
}

export function createRepresentationQuick(
  plugin: PluginContext,
  api: ViewerApi,
  ref: string,
  options: { quietMs?: number; throttleMs?: number } = {},
): RepresentationQuickSource | undefined {
  const cell = () => plugin.state.data.cells.get(ref);
  const structure = () =>
    plugin.state.data.cells.get(cell()?.transform.parent ?? '')?.obj?.data as Structure | undefined;
  const params = () => cell()?.transform.params as Params | undefined;
  if (!params()?.type || !structure()) return undefined;
  const registry = plugin.representation.structure.registry;
  const themes = plugin.representation.structure.themes;

  // A type's defaults on a structure, once each: a size or opacity drag reads them at pointer rate.
  const defaultsOf = new WeakMap<Structure, Map<string, Params['type']['params'] | undefined>>();
  const typeDefaults = (name: string) => {
    const s = structure();
    if (!s) return undefined;
    let known = defaultsOf.get(s);
    if (!known) defaultsOf.set(s, (known = new Map()));
    if (!known.has(name)) {
      const provider = registry.get(name);
      known.set(
        name,
        provider ? (PD.getDefaultValues(provider.getParams(themes, s)) as Params['type']['params']) : undefined,
      );
    }
    return known.get(name);
  };
  const themeDefaults = (name: string) => {
    const s = structure();
    const provider = themes.colorThemeRegistry.get(name);
    return s && provider ? (PD.getDefaultValues(provider.getParams({ structure: s })) as Record<string, unknown>) : {};
  };
  /** The common choices that apply, and the current one when it is not among them. */
  const choices = (
    common: readonly [string, string][],
    current: string,
    ok: (id: string) => boolean,
    label: (id: string) => string,
    extra: (id: string) => Partial<QuickChoice>,
  ) => {
    const list: QuickChoice[] = common.filter(([id]) => ok(id)).map(([id, l]) => ({ id, label: l, ...extra(id) }));
    if (current && !list.some((c) => c.id === current)) {
      list.push({ id: current, label: label(current), ...extra(current) });
    }
    return list;
  };
  /** A theme's colors on this structure (its legend, with its defaults), once each. */
  const swatches = new Map<string, string | undefined>();
  const swatchOf = (id: string): Partial<QuickChoice> => {
    if (id === 'uniform') return {};
    if (!swatches.has(id)) {
      const s = structure();
      let swatch: string | undefined;
      try {
        const provider = themes.colorThemeRegistry.get(id);
        if (s && provider) swatch = legendSwatch(provider.factory({ structure: s }, themeDefaults(id) as never).legend);
      } catch {
        swatch = undefined;
      }
      swatches.set(id, swatch);
    }
    const swatch = swatches.get(id);
    return swatch ? { swatch } : {};
  };
  // Another type (the current one, not among the common: interactions, putty, …): a generic glyph; its name is the tooltip
  // and the setting's head.
  const iconOf = (id: string): Partial<QuickChoice> => ({ icon: TypeIcons[id] ?? 'shapes' });
  const typeOk = (id: string) => {
    const s = structure();
    return registry.has(id) && !!s && registry.get(id).isApplicable(s);
  };
  const themeOk = (id: string) => {
    const s = structure();
    if (!themes.colorThemeRegistry.has(id) || !s || !themes.colorThemeRegistry.get(id).isApplicable({ structure: s })) {
      return false;
    }
    // pLDDT only for a prediction: the theme takes any B-factors of a model not marked experimental, and an older file
    // (no determination methodology) of an X-ray structure passes; its B-factors are not confidences.
    return id !== 'plddt-confidence' || s.models.some(isPredicted);
  };

  const types = observableState<readonly QuickChoice[]>([]);
  const themeChoices = observableState<readonly QuickChoice[]>([]);
  const type = observableState('');
  const theme = observableState('');
  const color = observableState<number | undefined>(undefined);
  const opacity = observableState<number | undefined>(undefined);
  const size = observableState<number | undefined>(undefined);

  // The choices are listed for these (and the structure), again when they change or on a refresh (the hierarchy changed),
  // not for each step of a drag, which changes neither.
  let listed: { type: string; theme: string; structure: Structure | undefined } | undefined;
  const publish = (refreshed = false) => {
    const p = params();
    if (!p) return;
    const o = steps.overlay('rep');
    const t = o?.type ?? p.type.name;
    const th = o?.color !== undefined ? 'uniform' : (o?.theme ?? p.colorTheme.name);
    type.set(t);
    theme.set(th);
    const s = structure();
    if (refreshed || listed?.type !== t || listed.theme !== th || listed.structure !== s) {
      listed = { type: t, theme: th, structure: s };
      const typeList = choices(CommonTypes, t, typeOk, (id) => registry.get(id)?.label ?? id, iconOf);
      if (!sameChoices(types.value, typeList)) types.set(typeList);
      const themeList = choices(
        CommonThemes,
        th,
        themeOk,
        (id) => themes.colorThemeRegistry.get(id)?.label ?? id,
        swatchOf,
      );
      if (!sameChoices(themeChoices.value, themeList)) themeChoices.set(themeList);
    }
    color.set(o?.color ?? (p.colorTheme.name === 'uniform' ? (p.colorTheme.params.value as number) : undefined));
    opacity.set(o?.opacity ?? p.type.params.alpha);
    const base =
      o?.size === undefined && p.type.params.sizeFactor !== undefined
        ? typeDefaults(p.type.name)?.sizeFactor
        : undefined;
    size.set(
      o?.size ?? (p.type.params.sizeFactor !== undefined && base ? round(p.type.params.sizeFactor / base) : undefined),
    );
  };

  const steps = gestureSteps<RepresentationQuickChange>(api, {
    label: 'Representation',
    apply: async (g, _key, change) => {
      const p = params();
      if (!p) return;
      let type = p.type;
      if (change.type && change.type !== p.type.name) {
        // Another type from its defaults, keeping the opacity chosen.
        const defaults = typeDefaults(change.type);
        if (!defaults) return;
        type = {
          name: change.type,
          params: { ...defaults, ...(p.type.params.alpha !== undefined ? { alpha: p.type.params.alpha } : {}) },
        };
      }
      if (change.opacity !== undefined) type = { ...type, params: { ...type.params, alpha: change.opacity } };
      if (change.size !== undefined) {
        const base = typeDefaults(type.name)?.sizeFactor;
        if (base !== undefined) type = { ...type, params: { ...type.params, sizeFactor: round(base * change.size) } };
      }
      let colorTheme = p.colorTheme;
      if (change.theme && change.theme !== p.colorTheme.name) {
        colorTheme = { name: change.theme, params: themeDefaults(change.theme) };
      }
      if (change.color !== undefined) {
        const base = colorTheme.name === 'uniform' ? colorTheme.params : themeDefaults('uniform');
        colorTheme = { name: 'uniform', params: { ...base, value: Color(change.color) } };
      }
      await g.run('representation.update', { representation: ref, params: { type, colorTheme } });
    },
    changed: () => publish(),
    ...(options.quietMs !== undefined ? { quietMs: options.quietMs } : {}),
    ...(options.throttleMs !== undefined ? { throttleMs: options.throttleMs } : {}),
  });

  publish(true);
  return {
    types: types.asReadonly,
    themes: themeChoices.asReadonly,
    type: type.asReadonly,
    theme: theme.asReadonly,
    color: color.asReadonly,
    opacity: opacity.asReadonly,
    size: size.asReadonly,
    change: (c) => steps.change('rep', c),
    settle: () => steps.settle(),
    refresh: () => publish(true),
  };
}
