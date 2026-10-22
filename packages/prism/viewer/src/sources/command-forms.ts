/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { CommandValues } from '@molstar/prism-ui/core';
import { inStructures, isFailure, someStructures } from '../api/scope.js';
import type { StyleChoice, StyleFormsModel, TypeOption } from '../model/command-forms.js';
import { CommonTypes, legendSwatch, themeApplies } from './representation-quick.js';
import { ReprIcons } from './scene.js';

// The command menu's own forms for a component's look over the plugin (model/command-forms.ts): the component and the
// structure resolved as the operations resolve them (`api/scope.ts`), the types and themes from the registries.

/** A type's main options, in order, as the form names them; a type shows those it has. */
const MainOptions: readonly (readonly [name: string, label: string])[] = [
  ['sizeFactor', 'Size'],
  ['quality', 'Quality'],
  ['alpha', 'Opacity'],
];

const str = (v: CommandValues[string]) => (typeof v === 'string' ? v : '');

export function createStyleForms(plugin: PluginContext): StyleFormsModel {
  const ctx = { plugin };
  const { registry, themes } = plugin.representation.structure;
  // The component as `representation.set` finds it; none while the values name none (or several).
  const component = (values: CommandValues) => {
    const found = inStructures(ctx, str(values.component), str(values.structure));
    return isFailure(found) ? undefined : found.component;
  };
  // A theme's colors on a structure, once each (a legend is a theme made on the whole structure).
  // The structure the values name (or the only one); none when they name several.
  const structure = (values: CommandValues) => {
    const found = someStructures(ctx, str(values.structure) || undefined);
    return isFailure(found) || found.length !== 1 ? undefined : found[0].cell.obj?.data;
  };
  // What the types are read on: the component, else (none named yet, from the command menu) its structure.
  const dataOf = (values: CommandValues) => component(values)?.cell.obj?.data ?? structure(values);
  const swatches = new WeakMap<Structure, Map<string, string | undefined>>();
  const swatchOf = (s: Structure, id: string) => {
    let known = swatches.get(s);
    if (!known) swatches.set(s, (known = new Map()));
    if (!known.has(id)) {
      let swatch: string | undefined;
      try {
        const provider = themes.colorThemeRegistry.get(id);
        const params = PD.getDefaultValues(provider.getParams({ structure: s }));
        swatch = legendSwatch(provider.factory({ structure: s }, params as never).legend);
      } catch {
        swatch = undefined;
      }
      known.set(id, swatch);
    }
    return known.get(id);
  };

  return {
    types(values) {
      const c = component(values);
      const data = dataOf(values);
      if (!data) return [];
      const shown = new Set(
        c?.representations.map((r) => (r.cell.transform.params as { type?: { name: string } })?.type?.name),
      );
      const applicable = registry.getApplicableTypes(data);
      const rank = (id: string) => {
        const i = CommonTypes.findIndex(([t]) => t === id);
        return i < 0 ? CommonTypes.length : i;
      };
      return [...applicable]
        .sort((a, b) => rank(a[0]) - rank(b[0]))
        .map(
          ([id, label]): StyleChoice => ({
            id,
            label,
            icon: ReprIcons[id] ?? 'shapes',
            ...(shown.has(id) ? { shown: true } : {}),
          }),
        );
    },
    typeOptions(type, values) {
      const data = dataOf(values);
      if (!data || !registry.has(type)) return [];
      const params = registry.get(type).getParams(themes, data);
      // Added as the component manager adds one (`representation.set`): its quality is the default.
      const quality = plugin.managers.structure.component.state.options.visualQuality;
      const out: TypeOption[] = [];
      for (const [name, label] of MainOptions) {
        const p = params[name];
        if (p?.type === 'number' && p.min !== undefined && p.max !== undefined) {
          out.push({
            kind: 'range',
            name,
            label,
            min: p.min,
            max: p.max,
            step: p.step ?? 0.01,
            default: p.defaultValue,
          });
        } else if (p?.type === 'select' && p.options.every(([v]) => typeof v === 'string')) {
          const options = p.options as readonly (readonly [string, string])[];
          const fallback = name === 'quality' && options.some(([v]) => v === quality) ? quality : p.defaultValue;
          out.push({ kind: 'choice', name, label, options, default: fallback as string });
        }
      }
      return out;
    },
    themes(values) {
      const s = structure(values);
      if (!s) return [];
      // In the registry's order, by its categories (as the inspector's theme select lists them).
      return themes.colorThemeRegistry.types
        .filter(([id]) => themeApplies(plugin, s, id))
        .map(([id, label, category]): StyleChoice => {
          const swatch = id === 'uniform' ? undefined : swatchOf(s, id);
          return { id, label, ...(category ? { category } : {}), ...(swatch ? { swatch } : {}) };
        });
    },
  };
}
