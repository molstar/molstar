/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, isState, type ReadonlyState, State } from '@molstar/core/reactive';
import { ExpandableList, Expandable } from '../components/expandable/expandable.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { ParamsForm, type ParamsFormLabels, searchParams } from '../components/params/params-form.js';
import { TextInput } from '../components/text-input/text-input.js';
import { type Child, effect, For, type MaybeReactive, onCleanup, Show } from '../dom/index.js';
import { ViewEmpty, ViewStack, ViewToolbar } from '../shell/view-frame.js';
import { classes } from './_generated/settings-view.module.css.js';
import type { ParamsFormModel } from './params-form.js';

/**
 * A group of settings by intent (Appearance, Camera, Interaction): a form the app fills (a `ParamsFormModel` over its
 * params, changes applied as they are made). The view knows no app.
 */
export interface SettingsGroup {
  readonly id: string;
  readonly title: string;
  readonly icon?: IconName;
  readonly form: ParamsFormModel<any>;
}

/** Settings' words (styles §8): strings, or functions for what names the search. */
export interface SettingsLabels {
  /** The toolbar's accessible name. */
  readonly toolbar: string;
  readonly search: string;
  readonly searchPlaceholder: string;
  readonly noMatch: (query: string) => string;
}

export const DefaultSettingsLabels: SettingsLabels = {
  toolbar: 'Settings',
  search: 'Search settings',
  searchPlaceholder: 'Search settings',
  noMatch: (q) => `No setting matches “${q}”.`,
};

/**
 * Settings (styles §8): a search over every group's settings, then the groups as `Expandable` panels, the first open.
 * While searching, only matching rows show, their groups open and the rest hide; the words may span a group's title and
 * its settings, and a group whose title takes them all shows whole. With no match, "No setting matches". Each group's
 * search is one walk per query and values, read by the group and its rows alike. The groups may be a state (an app
 * adding an extension's). Built with the view primitives, so it sits in any view.
 */
export function SettingsView(props: {
  groups: MaybeReactive<readonly SettingsGroup[]>;
  /** The search text: the host's (a command can search), else the view's own. */
  query?: State<string>;
  /** After the groups (an app's extensions' sections). */
  more?: () => Child;
  labels?: Partial<SettingsLabels>;
  /** The groups' forms' words. */
  formLabels?: Partial<ParamsFormLabels>;
}) {
  const l = { ...DefaultSettingsLabels, ...props.labels };
  const query = props.query ?? new State('');
  const groups: ReadonlyState<readonly SettingsGroup[]> = isState(props.groups)
    ? (props.groups as ReadonlyState<readonly SettingsGroup[]>)
    : new State(props.groups as readonly SettingsGroup[]);
  const searching = derive(query, (t) => t.trim().length > 0);
  // Each group says whether it shows (by id): "No setting matches" when none does, read from their one walk each.
  const shown = new State<ReadonlyMap<string, boolean>>(new Map());
  const report = (id: string, on: boolean | undefined) => {
    if (shown.value.get(id) === on) return;
    const next = new Map(shown.value);
    if (on === undefined) next.delete(id);
    else next.set(id, on);
    shown.set(next);
  };
  const none = derive([searching, shown, groups], (s, m, list) => s && list.every((g) => m.get(g.id) === false));
  return (
    <ViewStack class={classes.settings}>
      <ViewToolbar label={l.toolbar}>
        <TextInput
          search
          class={classes.search}
          aria-label={l.search}
          placeholder={l.searchPlaceholder}
          value={query}
          onInput={(v) => query.set(v)}
        />
      </ViewToolbar>
      <ExpandableList>
        <For each={groups} by={(g) => g.id}>
          {(group, index) => (
            <Group
              group={group}
              query={query}
              searching={searching}
              open={index.value === 0}
              report={report}
              labels={props.formLabels}
            />
          )}
        </For>
      </ExpandableList>
      <Show when={none}>{() => <ViewEmpty>{derive(query, (q) => l.noMatch(q.trim()))}</ViewEmpty>}</Show>
      {props.more?.()}
    </ViewStack>
  );
}

/** A group: an `Expandable` holding its form, filtered by the search; open while searching, hidden with no match. */
function Group(props: {
  group: SettingsGroup;
  query: ReadonlyState<string>;
  searching: ReadonlyState<boolean>;
  open: boolean;
  report(id: string, shown: boolean | undefined): void;
  labels?: Partial<ParamsFormLabels>;
}) {
  const { group, query } = props;
  // The person's own open or closed; the search opens it over that, and lets it go back after.
  const own = new State(props.open);
  const open = derive([own, props.searching], (o, s) => o || s);
  // One walk per query and values: whether the group shows, and which of its rows do (the form reads it). Its title
  // takes the words it contains, its settings the rest ('': every setting in it shows).
  const search = derive([query, group.form.values], (t, values) =>
    searchParams(group.form.params, t, { values, scope: group.title }),
  );
  const hidden = derive(search, (s) => !!s && !s.any);
  effect(hidden, (h) => props.report(group.id, !h));
  onCleanup(() => props.report(group.id, undefined));
  return (
    <div hidden={hidden} data-settings-group={group.id}>
      <Expandable
        title={
          group.icon ? (
            <span class={classes.title}>
              <Icon name={group.icon} />
              {group.title}
            </span>
          ) : (
            group.title
          )
        }
        open={open}
        onToggle={(o) => own.set(o)}
      >
        {() => <ParamsForm form={group.form} search={search} labels={props.labels} />}
      </Expandable>
    </div>
  );
}
