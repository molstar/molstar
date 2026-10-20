/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Choices } from '../../components/choices/choices.js';
import { Icon } from '../../components/icon/icon.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { Dynamic, For } from '../../dom/control-flow.js';
import type { StatusMessage, NotificationsManager } from '../../managers/notifications.js';
import { observableState } from '../../reactive/observable-state.js';
import { derive } from '../../reactive/state.js';
import { classes } from './_generated/log-view.module.css.js';
import { KindIcons } from './canvas-chrome.js';

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

/** What the Log shows: everything, the messages the status line showed, one log-only source, or the problems. */
type LogFilter = 'all' | 'messages' | 'problems' | `source:${string}`;

const passes = (m: StatusMessage, filter: LogFilter, query: string) => {
  if (filter === 'messages' && m.source !== undefined) return false;
  if (filter === 'problems' && m.kind !== 'warning' && m.kind !== 'error') return false;
  if (filter.startsWith('source:') && m.source !== filter.slice(7)) return false;
  return !query || m.text.toLowerCase().includes(query) || !!m.detail?.toLowerCase().includes(query);
};

/**
 * The Log tab: every message the status line showed and what was written to the log only (the plugin's transforms,
 * each with how long it took), newest first, with its kind, its detail and its time. Filtered by what wrote it or the
 * problems, and by text.
 */
export function LogView(props: { status: NotificationsManager }) {
  const filter = observableState<LogFilter>('all');
  const query = observableState('');
  // The log-only sources present, in the order they first appear: a filter each.
  const sources = derive(
    props.status.logged,
    (list) => {
      const seen = new Set<string>();
      for (const m of list) if (m.source !== undefined) seen.add(m.source);
      return [...seen].sort();
    },
    (a, b) => a.length === b.length && a.every((s, i) => s === b[i]),
  );
  const choices = derive(sources, (list) => [
    { id: 'all', label: 'All' },
    { id: 'messages', label: 'Messages' },
    ...list.map((s) => ({ id: `source:${s}`, label: s })),
    { id: 'problems', label: 'Problems' },
  ]);
  const shown = derive(
    [props.status.logged, filter.asReadonly, query.asReadonly],
    (list, f, q): readonly StatusMessage[] => {
      const text = q.trim().toLowerCase();
      return f === 'all' && !text ? list : list.filter((m) => passes(m, f, text));
    },
  );
  // Not live (a log is, by default): each message is said once, by the status line's region (accessibility.md).
  return (
    <div class={classes.log}>
      <div class={classes.filters}>
        {/* Made again only when the sources change (a first transform logged). */}
        <Dynamic
          component={derive(choices, (list) => () => (
            <Choices
              label="Show"
              choices={list}
              active={filter.asReadonly}
              onPick={(id) => filter.set(id as LogFilter)}
            />
          ))}
        />
        <TextInput
          search
          size="sm"
          placeholder="Filter"
          aria-label="Filter the log"
          value={query.asReadonly}
          onInput={(v) => query.set(v)}
        />
      </div>
      <div class={classes.entries} role="log" aria-live="off" aria-label="Log">
        <For
          each={shown}
          key={(m) => m.id}
          fallback={() => (
            <p class={classes.empty}>
              {derive(props.status.logged, (h) => (h.length ? 'Nothing matches.' : 'No messages yet.'))}
            </p>
          )}
        >
          {(m) => (
            <div class={classes.entry} data-kind={m.kind} data-source={m.source}>
              <Icon name={KindIcons[m.kind]} size={14} class={classes.kind} />
              <span class={classes.text}>{m.text}</span>
              {m.detail ? <span class={classes.detail}>{m.detail}</span> : <span />}
              <time class={classes.time} dateTime={m.time.toISOString()}>
                {time.format(m.time)}
              </time>
            </div>
          )}
        </For>
      </div>
    </div>
  );
}
