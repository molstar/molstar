/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { StructureComponentManager } from '@molstar/plugin/state/manager/structure/component';
import { observableState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerApi } from '../api/area.js';
import { Looks, RepresentationStyles } from '../api/styles.js';
import type { PresetsModel, QuickStylesModel, StructurePreset } from '../model/scene.js';
import type { SourceScope } from './state.js';

// The Structure tab's header models (state.md §1, rows 6 and 8): the representation presets from the plugin's registry
// (grouped as it groups them; operations.md §2, "registries are arguments": the preset is `structure.preset`'s id), and
// the component manager's options as a form (`component.options`).

export function createPresetsModel(plugin: PluginContext, api: ViewerApi): PresetsModel {
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
  const current = observableState<string | undefined>(undefined);
  return {
    groups: [...groups].map(([label, presets]) => ({ ...(label ? { label } : {}), presets })),
    current: current.asReadonly,
    apply: (id) => {
      void api.run('structure.preset', { preset: id }).then((r) => {
        if (r.ok) current.set(id);
      });
    },
  };
}

/** The component manager's options, applied to every component through `component.options`. */
export function createComponentOptionsForm(scope: SourceScope, plugin: PluginContext, api: ViewerApi) {
  const params = StructureComponentManager.OptionsParams;
  let timer: ReturnType<typeof setTimeout> | undefined;
  const form = new ParamsFormModel(params, plugin.managers.structure.component.state.options, {
    onChange: (values) => {
      clearTimeout(timer);
      timer = setTimeout(() => void api.run('component.options', values as Record<string, unknown>), 150);
    },
  });
  // Changed elsewhere (a session, Undo): the form follows.
  scope.add(
    plugin.managers.structure.component.events.optionsUpdated.subscribe(() => {
      const options = plugin.managers.structure.component.state.options;
      if (!PD.areEqual(params, form.values.value, options)) form.load(options);
    }),
  );
  scope.add(() => clearTimeout(timer));
  return form;
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
