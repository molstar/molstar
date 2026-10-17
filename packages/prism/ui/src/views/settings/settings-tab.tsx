/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm, paramsMatch } from '../../components/params/params-form.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import type { SettingsGroup, SettingsModel } from '../../model/settings.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { classes } from './_generated/settings-tab.module.css.js';

import { derive } from '../../reactive/state.js';
/**
 * The Settings tab: a search over every setting, then the groups by intent (essentials first), each a form that folds
 * open, then `more`. While searching, only matching rows show and their groups open; a group with no match hides.
 */
export function SettingsTab(props: {
  model: SettingsModel;
  /** More after the groups (the app's extension sections). */
  more?: () => Child;
}) {
  const { model } = props;
  const query = observableState('');
  const q = query.asReadonly;
  const none = derive(
    q,
    (text) => text.trim().length > 0 && !model.groups.some((g) => paramsMatch(g.form.params, text)),
  );
  return (
    <div class={classes.tab}>
      <div class={classes.search}>
        <Icon name="search" size={14} class={classes.searchIcon} />
        <TextInput
          search
          size="sm"
          class={classes.searchInput}
          aria-label="Search settings"
          placeholder="Search settings"
          value={q}
          onInput={(v) => query.set(v)}
        />
      </div>
      <ScrollArea class={classes.groups} contentClass={classes.groupsContent}>
        {model.groups.map((group, i) => (
          <Group group={group} query={q} open={i === 0} />
        ))}
        <Show when={none}>{() => <p class={classes.none}>No setting matches “{q}”.</p>}</Show>
        {props.more?.()}
      </ScrollArea>
    </div>
  );
}

/** One group: a header that folds, and its form, filtered by the search. */
function Group(props: { group: SettingsGroup; query: ReadonlyState<string>; open: boolean }) {
  const { group, query } = props;
  const own = observableState(props.open);
  const searching = derive(query, (t) => t.trim().length > 0);
  const open = combineLatest([own.asReadonly, searching]).pipe(
    map(([o, s]) => o || s),
    distinctUntilChanged(),
  );
  const hidden = derive(query, (t) => !paramsMatch(group.form.params, t));
  return (
    <section class={classes.group} hidden={hidden} aria-label={group.title}>
      <div class={classes.header}>
        <button type="button" class={classes.toggle} aria-expanded={open} onClick={() => own.set(!own.value)}>
          <span class={classes.title}>{group.title}</span>
          <Icon name="chevron-down" size={14} class={classes.chevron} />
        </button>
      </div>
      <Show when={open} reveal>
        {() => (
          <div class={classes.body} data-wide-labels={group.wideLabels ? '' : undefined}>
            <ParamsForm form={group.form} filter={query} />
          </div>
        )}
      </Show>
    </section>
  );
}
