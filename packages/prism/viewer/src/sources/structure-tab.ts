/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { observableState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import { Looks, RepresentationStyles } from '../api/styles.js';
import type { PresetsModel, QuickStylesModel, StructurePreset } from '../model/scene.js';
import type { SourceScope } from './state.js';

// The Structure tab's header models (state.md §1, rows 6 and 26): the representation presets from the plugin's registry
// (grouped as it groups them; api.md §2, "registries are arguments": the preset is `structure.preset`'s id), and
// the quick styles. The component manager's options are the look's scene defaults (`sources/look.ts`).

export function createPresetsModel(scope: SourceScope, plugin: PluginContext, api: ViewerApi): PresetsModel {
  const groups = new Map<string, StructurePreset[]>();
  for (const p of plugin.builders.structure.representation.providers) {
    const group = p.display.group ?? '';
    if (!groups.has(group)) groups.set(group, []);
    groups.get(group)!.push({
      id: p.id,
      label: p.display.name,
      ...(p.display.description ? { description: p.display.description } : {}),
    });
  }
  // The one in effect until the scene changes some other way (Undo, a load, Styles): the state tree it left is the mark.
  const current = observableState<string | undefined>(undefined);
  let applying = 0;
  let applied: unknown;
  scope.add(
    plugin.state.data.events.changed.subscribe(() => {
      if (!applying && current.value !== undefined && plugin.state.data.tree !== applied) current.set(undefined);
    }),
  );
  return {
    groups: [...groups].map(([label, presets]) => ({ ...(label ? { label } : {}), presets })),
    current: current.asReadonly,
    apply: (id) => {
      applying++;
      void api
        .run('structure.preset', { preset: id })
        .then((r) => {
          if (!r.ok) return;
          applied = plugin.state.data.tree;
          current.set(id);
        })
        .finally(() => void applying--);
    },
  };
}

/** Quick styles (the Structure tab's Styles): representation styles and looks, both `styles.apply`. */
export function createQuickStylesModel(api: ViewerApi): QuickStylesModel {
  return {
    representations: Object.entries(RepresentationStyles).map(([id, s]) => ({ id, label: s.label })),
    styles: Object.entries(Looks).map(([id, s]) => ({ id, label: s.label })),
    applyRepresentation: (id) => void api.run('styles.apply', { style: id }),
    applyStyle: (id) => void api.run('styles.apply', { style: id }),
  };
}
