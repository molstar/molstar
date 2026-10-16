/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerApi } from '../api/area.js';
import { actionName } from '../api/load.js';
import type { HomeAction } from '../model/home.js';
import type { SourceScope } from './state.js';

// Home's actions (state.md §1, the registry row): the state actions applicable to the root, as forms, each applied by
// `load.action` as the person. Followed as actions are registered (an extension loaded later).

const Icons: Readonly<Record<string, IconName>> = {
  'download-structure': 'download',
  'open-files': 'folder-open',
  'download-file': 'link',
  'download-density': 'volume',
};

export function createHomeActions(
  scope: SourceScope,
  plugin: PluginContext,
  api: ViewerApi,
): ReadonlyState<readonly HomeAction[]> {
  const state = observableState<readonly HomeAction[]>([]);
  // A form keeps what was typed in it and its busy flag: an action registered or removed rebuilds only its own.
  let made = new Map<string, { source: unknown; action: HomeAction }>();
  const build = () => {
    const root = plugin.state.data.cells.get(StateTransform.RootRef);
    if (!root?.obj) return;
    const list = plugin.state.data.actions.fromCell(root, plugin);
    const next = new Map<string, { source: unknown; action: HomeAction }>();
    state.set(
      list
        // Opening files is the drop area's and the picker's (load.files); sessions are their own.
        .filter((a) => actionName(a) !== 'open-files')
        .map((a): HomeAction => {
          const key = actionName(a);
          const kept = made.get(key);
          if (kept && kept.source === a) {
            next.set(key, kept);
            return kept.action;
          }
          const params = a.definition.params?.(root.obj!, plugin) as PD.Params | undefined;
          const form = params ? new ParamsFormModel(params) : undefined;
          const busy = observableState(false);
          const action: HomeAction = {
            id: key,
            title: a.definition.display.name,
            icon: Icons[key] ?? 'download',
            ...(a.definition.display.description ? { description: a.definition.display.description } : {}),
            ...(form ? { form } : {}),
            applyLabel: a.definition.display.triggerLabel ?? 'Apply',
            busy: busy.asReadonly,
            apply: async () => {
              busy.set(true);
              try {
                await api.load.action({ action: key, params: (form?.values.value ?? {}) as Record<string, unknown> });
              } finally {
                busy.set(false);
              }
            },
          };
          next.set(key, { source: a, action });
          return action;
        }),
    );
    made = next;
  };
  build();
  scope.add(plugin.state.data.actions.events.added.subscribe(build));
  scope.add(plugin.state.data.actions.events.removed.subscribe(build));
  return state.asReadonly;
}
