/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, domView, onCleanup, SectionList, SettingsTab, type ShellTab, Show } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { SettingsGroup, SettingsModel } from '@molstar/prism-ui/model/settings';
import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import type { ViewerManagers } from '../api/managers.js';
import { type SettingsGroupState, withCheckpointing } from '../sources/settings.js';
import type { ViewerSources } from '../sources/viewer-sources.js';

// The Settings tab (parity.md §5): the groups of `sources.settings` as forms that follow the live values and apply
// through `settings.set` as the person, debounced like the Canvas card. A group of the tab is a step in the History.
// The groups the host can show change with the canvas and the behaviors: the tab is made again when they do. Other
// extensions' settings sections go after the groups (the `settings` slot).

const Debounce = 150;

export function settingsTab(sources: ViewerSources, managers: ViewerManagers, api: ViewerApi): ShellTab {
  const current = () => withCheckpointing(sources.settings.groups.value, managers.history.checkpointing.value);
  const live = combineLatest([sources.settings.groups, managers.history.checkpointing]).pipe(
    map(([groups, on]) => withCheckpointing(groups, on)),
  );
  // What decides the forms' shape (their params); values are loaded into them.
  const shape = live.pipe(
    map((groups) => groups.map((g) => `${g.id}:${Object.keys(g.params).join(',')}`).join('|')),
    distinctUntilChanged(),
  );
  return {
    id: 'settings',
    title: 'Settings',
    icon: 'settings',
    order: 2.7,
    place: 'foot',
    render: domView(({ ctx }) => {
      const more = () => <SectionList shell={ctx.shell} slot="settings" />;
      return (
        <Show when={shape} fallback={() => <SettingsTab model={{ groups: [] }} more={more} />}>
          {() => <SettingsTab model={modelOf(current(), live, api)} more={more} />}
        </Show>
      );
    }),
  };
}

/** The forms for the groups, each following its group's values and applying its edits. */
function modelOf(
  groups: readonly SettingsGroupState[],
  live: Observable<SettingsGroupState[]>,
  api: ViewerApi,
): SettingsModel {
  const forms = new Map<string, ParamsFormModel>();
  const timers = new Map<string, ReturnType<typeof setTimeout>>();
  // Edits sent and not yet seen back: a form keeps what the person is editing until they are.
  const busy = new Map<string, number>();
  let latest = groups;
  onCleanup(() => {
    for (const t of timers.values()) clearTimeout(t);
  });

  const send = (id: string, changed: Set<string>) => {
    const form = forms.get(id)!;
    const values = form.values.value as Record<string, unknown>;
    busy.set(id, (busy.get(id) ?? 0) + 1);
    void api.settings
      .set({ group: id as never, values: Object.fromEntries([...changed].map((k) => [k, values[k]])) })
      .finally(() => {
        busy.set(id, (busy.get(id) ?? 1) - 1);
        // Whatever the result, the form shows what the viewer has.
        if (!busy.get(id) && !timers.has(id)) load(id);
      });
  };
  const load = (id: string) => {
    const g = latest.find((x) => x.id === id);
    if (g) forms.get(id)?.load(g.values as never);
  };

  const model = groups.map((g): SettingsGroup => {
    const pending = new Set<string>();
    const form = new ParamsFormModel(g.params, g.values as never, {
      onChange: (_values, changed) => {
        pending.add(changed);
        clearTimeout(timers.get(g.id));
        timers.set(
          g.id,
          setTimeout(() => {
            timers.delete(g.id);
            const keys = new Set(pending);
            pending.clear();
            send(g.id, keys);
          }, Debounce),
        );
      },
    });
    forms.set(g.id, form);
    return { id: g.id, title: g.title, icon: g.icon, form, wideLabels: g.wideLabels };
  });

  bind(live, (gs) => {
    latest = gs;
    for (const g of gs) if (!busy.get(g.id) && !timers.has(g.id)) load(g.id);
  });
  return { groups: model };
}
