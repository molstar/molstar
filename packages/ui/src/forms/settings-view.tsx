/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { ExpandableList, Expandable } from '../components/expandable/expandable.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { ParamsForm, paramsMatch, searchRest } from '../components/params/params-form.js';
import { TextInput } from '../components/text-input/text-input.js';
import { type Child, Show } from '../dom/index.js';
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

/** The search's words a group's title leaves for its settings (lower case; '': the title matched them all). */
function groupRest(group: SettingsGroup, text: string): string {
  return searchRest(group.title, text.trim().toLowerCase().replace(/\s+/g, ' '));
}

/** Whether a group shows for the text: its title and settings match every word between them (any group for none). */
export function groupMatches(group: SettingsGroup, text: string): boolean {
  const rest = groupRest(group, text);
  return !rest || paramsMatch(group.form.params, rest);
}

/**
 * The groups that match the text (all of them for none): a headless search over the settings. Each word may match the
 * group's title or a setting's path (`searchRest`): "perf" finds Performance whole, "perf qua" its Quality.
 */
export function matchingGroups(groups: readonly SettingsGroup[], text: string): SettingsGroup[] {
  return groups.filter((g) => groupMatches(g, text));
}

/**
 * Settings (styles §8): a search over every group's settings, then the groups as `Expandable` panels, the first open.
 * While searching, only matching rows show, their groups open and the rest hide; the words may span a group's title and
 * its settings, and a group whose title takes them all shows whole. With no match, "No setting matches".
 * Built with the view primitives, so it sits in any view.
 */
export function SettingsView(props: {
  groups: readonly SettingsGroup[];
  /** The search text: the host's (a command can search), else the view's own. */
  query?: State<string>;
  /** After the groups (an app's extensions' sections). */
  more?: () => Child;
}) {
  const query = props.query ?? new State('');
  const searching = derive(query, (t) => t.trim().length > 0);
  const none = derive(query, (t) => !!t.trim() && !props.groups.some((g) => groupMatches(g, t)));
  return (
    <ViewStack class={classes.settings}>
      <ViewToolbar>
        <TextInput
          search
          class={classes.search}
          aria-label="Search settings"
          placeholder="Search settings"
          value={query}
          onInput={(v) => query.set(v)}
        />
      </ViewToolbar>
      <ExpandableList>
        {props.groups.map((group, i) => (
          <Group group={group} query={query} searching={searching} open={i === 0} />
        ))}
      </ExpandableList>
      <Show when={none}>
        {() => <ViewEmpty>{derive(query, (q) => `No setting matches “${q.trim()}”.`)}</ViewEmpty>}
      </Show>
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
}) {
  const { group, query } = props;
  // The person's own open or closed; the search opens it over that, and lets it go back after.
  const own = new State(props.open);
  const open = derive([own, props.searching], (o, s) => o || s);
  const hidden = derive(query, (t) => !groupMatches(group, t));
  // The form gets the words its title leaves ('': every setting in it shows).
  const filter = derive(query, (t) => groupRest(group, t));
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
        {() => <ParamsForm form={group.form} filter={filter} />}
      </Expandable>
    </div>
  );
}
