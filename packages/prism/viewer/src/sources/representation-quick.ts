/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { observableState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type {
  QuickChoice,
  RepresentationQuickChange,
  RepresentationQuickModel,
} from '../model/representation-quick.js';
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
  ['gaussian-surface', 'Gaussian'],
];

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

  const typeDefaults = (name: string) => {
    const s = structure();
    const provider = registry.get(name);
    return s && provider ? (PD.getDefaultValues(provider.getParams(themes, s)) as Params['type']['params']) : undefined;
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
  ) => {
    const list: QuickChoice[] = common.filter(([id]) => ok(id)).map(([id, l]) => ({ id, label: l }));
    if (current && !list.some((c) => c.id === current)) list.push({ id: current, label: label(current) });
    return list;
  };
  const typeOk = (id: string) => {
    const s = structure();
    return registry.has(id) && !!s && registry.get(id).isApplicable(s);
  };
  const themeOk = (id: string) => {
    const s = structure();
    return themes.colorThemeRegistry.has(id) && !!s && themes.colorThemeRegistry.get(id).isApplicable({ structure: s });
  };

  const types = observableState<readonly QuickChoice[]>([]);
  const themeChoices = observableState<readonly QuickChoice[]>([]);
  const type = observableState('');
  const theme = observableState('');
  const color = observableState<number | undefined>(undefined);
  const opacity = observableState<number | undefined>(undefined);
  const size = observableState<number | undefined>(undefined);

  const publish = () => {
    const p = params();
    if (!p) return;
    const o = steps.overlay('rep');
    const t = o?.type ?? p.type.name;
    const th = o?.color !== undefined ? 'uniform' : (o?.theme ?? p.colorTheme.name);
    type.set(t);
    theme.set(th);
    types.set(choices(CommonTypes, t, typeOk, (id) => registry.get(id)?.label ?? id));
    themeChoices.set(choices(CommonThemes, th, themeOk, (id) => themes.colorThemeRegistry.get(id)?.label ?? id));
    color.set(o?.color ?? (p.colorTheme.name === 'uniform' ? (p.colorTheme.params.value as number) : undefined));
    opacity.set(o?.opacity ?? p.type.params.alpha);
    const base = typeDefaults(p.type.name)?.sizeFactor;
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

  publish();
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
    refresh: publish,
  };
}
