/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import { Button } from '../components/button/button.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { TextInput } from '../components/text-input/text-input.js';
import { RecycledList, type VirtualSlot } from '../components/virtual-list/recycled-list.js';
import type { VirtualListHandle } from '../components/virtual-list/virtual-list.js';
import { effect, For, frames, on, setAttr, Show } from '../dom/index.js';
import { ViewEmpty, ViewRow, ViewToolbar } from '../shell/view-frame.js';
import { classes } from './_generated/log-view.module.css.js';
import { isLogAll, LogAll, type LogEntry, type LogKind, type LogModel } from './log-model.js';

/** An entry's kind as an icon (the log's rows, the status line's message). */
export const LogKindIcons: Readonly<Record<LogKind, IconName>> = {
  info: 'info',
  success: 'circle-check',
  warning: 'triangle-alert',
  error: 'circle-alert',
  debug: 'terminal',
};

/** An entry's time as the log says it (hours, minutes, seconds). */
export const logTime = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

/** An entry's description: its detail, its source, its time. */
export function logDescription(e: LogEntry): string {
  const at = logTime.format(e.time);
  if (e.detail && e.source) return `${e.detail} · ${e.source} · ${at}`;
  if (e.detail || e.source) return `${e.detail ?? e.source} · ${at}`;
  return at;
}

/** The kinds with a toggle (debug entries show under All). */
const ToggledKinds = ['info', 'success', 'warning', 'error'] as const;

/** The log's words (log §3): strings, or functions for what names a source. */
export interface LogViewLabels {
  /** The toolbar's accessible name. */
  readonly filters: string;
  readonly kinds: string;
  readonly sources: string;
  readonly all: string;
  readonly allTip: string;
  /** A kind's toggle: its name and its tooltip. */
  readonly kind: Readonly<Record<(typeof ToggledKinds)[number], { readonly label: string; readonly tip: string }>>;
  readonly messages: string;
  readonly messagesTip: string;
  readonly sourceTip: (source: string) => string;
  readonly filter: string;
  /** The text filter's accessible name. */
  readonly filterLabel: string;
  /** The list's accessible name (the `label` prop wins). */
  readonly list: string;
  readonly empty: string;
  readonly noMatch: string;
}

export const DefaultLogViewLabels: LogViewLabels = {
  filters: 'Log filters',
  kinds: 'Kinds',
  sources: 'Sources',
  all: 'All',
  allTip: 'Show everything',
  kind: {
    info: { label: 'Info', tip: 'Show the info messages' },
    success: { label: 'Success', tip: 'Show what succeeded' },
    warning: { label: 'Warnings', tip: 'Show the warnings' },
    error: { label: 'Errors', tip: 'Show the errors' },
  },
  messages: 'Messages',
  messagesTip: 'Show the messages (no source)',
  sourceTip: (s) => `Show ${s}`,
  filter: 'Filter',
  filterLabel: 'Filter the log',
  list: 'Log',
  empty: 'Nothing is logged yet.',
  noMatch: 'Nothing matches the filter.',
};

interface SourceChoice {
  /** `null`: the messages (no source). */
  readonly source: string | null;
  readonly label: string;
  readonly tip: string;
  readonly icon: IconName;
}

export interface LogViewProps {
  log: LogModel;
  /** A source's icon in the filters (a tag when absent). */
  sourceIcons?: Readonly<Record<string, IconName>>;
  /** The list's accessible name (the labels' `list`). */
  label?: string;
  labels?: Partial<LogViewLabels>;
  class?: string;
}

/**
 * A log (log §3): the filters (All, a toggle per kind, then Messages and one per source; each named in its tooltip; a
 * status line opens it at `LogProblems`, the warnings and errors toggled) and a text filter
 * (Escape empties it) over a recycled list of one-line rows (a kind icon, the text fading where cut, the detail, source
 * and time as its description, the whole in the row's side tooltip). Newest last, kept at the bottom while entries
 * come unless scrolled up. `role="log"` but not live: a status line speaks. It fills its parent (a flex column).
 */
export function LogView(props: LogViewProps) {
  const { log } = props;
  const l = { ...DefaultLogViewLabels, ...props.labels };
  // The sources' filters only once something wrote with a source (the messages alone need none).
  const sources = derive(log.sources, (list): readonly SourceChoice[] =>
    list.length
      ? [
          { source: null, label: l.messages, tip: l.messagesTip, icon: 'message-square' },
          ...list.map(
            (s): SourceChoice => ({ source: s, label: s, tip: l.sourceTip(s), icon: props.sourceIcons?.[s] ?? 'tag' }),
          ),
        ]
      : [],
  );
  const empty = derive([log.filtered, log.entries], (shown, all) =>
    shown.length ? undefined : all.length ? l.noMatch : l.empty,
  );

  // Kept at the bottom (log §3): while the list ends in view, a change scrolls to its new end, once a frame; scrolled
  // up, it stays where it is.
  let handle: VirtualListHandle | undefined;
  let atEnd = true;
  const toEnd = frames(() => {
    const n = log.filtered.value.length;
    if (atEnd && n) handle?.scrollToIndex(n - 1);
  });
  effect(log.filtered, () => {
    if (atEnd) toEnd.request();
  });
  // The list's viewport scrolls (a scroll does not bubble: heard in the capture phase, on the view's root).
  const follow = (root: HTMLElement) =>
    on(
      root,
      'scroll',
      (e: Event) => {
        const viewport = e.target as HTMLElement;
        if (!viewport.querySelector?.(':scope > [role="log"]')) return;
        // Read only (no write in this handler): within half a row of the end is at the end.
        atEnd = viewport.scrollTop + viewport.clientHeight >= viewport.scrollHeight - 12;
      },
      { capture: true, passive: true },
    );

  return (
    <div class={props.class ? `${classes.log} ${props.class}` : classes.log} ref={follow}>
      <ViewToolbar label={l.filters} class={classes.filters}>
        <div class={classes.kinds} role="group" aria-label={l.kinds}>
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="list" />}
            title={l.allTip}
            aria-label={l.all}
            pressed={derive(log.filter, isLogAll)}
            onClick={() => log.setFilter(LogAll)}
          />
          {ToggledKinds.map((kind) => (
            <Button
              variant="ghost"
              size="sm"
              icon={<Icon name={LogKindIcons[kind]} />}
              title={l.kind[kind].tip}
              aria-label={l.kind[kind].label}
              data-log-kind={kind}
              pressed={derive(log.filter, (f) => !!f.kinds?.includes(kind))}
              onClick={() => log.toggleKind(kind)}
            />
          ))}
        </div>
        <Show when={derive(sources, (list) => list.length > 0)}>
          {() => (
            <div class={classes.kinds} role="group" aria-label={l.sources}>
              <For each={sources} by={(c) => c.source ?? ''}>
                {(c) => (
                  <Button
                    variant="ghost"
                    size="sm"
                    icon={<Icon name={c.icon} />}
                    title={c.tip}
                    aria-label={c.label}
                    pressed={derive(log.filter, (f) => f.source === c.source)}
                    onClick={() => log.setSource(log.filter.value.source === c.source ? undefined : c.source)}
                  />
                )}
              </For>
            </div>
          )}
        </Show>
        <TextInput
          search
          size="sm"
          class={classes.search}
          placeholder={l.filter}
          aria-label={l.filterLabel}
          value={log.query}
          onInput={(v) => log.setQuery(v)}
        />
      </ViewToolbar>
      <Show when={empty}>{(text) => <ViewEmpty>{text}</ViewEmpty>}</Show>
      <RecycledList
        items={log.filtered}
        itemKey={(e) => e.id}
        class={classes.entries}
        edge
        role="log"
        aria-label={props.label ?? l.list}
        ref={(h) => {
          handle = h;
        }}
        contentRef={(el) => {
          el.setAttribute('aria-live', 'off');
          // At the end once laid out: the newest entry shows first.
          toEnd.request();
        }}
        slot={logRow}
      />
    </div>
  );
}

/** A recycled row of the log: made once, shown one entry, then another as the list scrolls or grows. */
function logRow(): VirtualSlot<LogEntry> {
  const label = new State('');
  const description = new State('');
  // An icon per kind, made the first time this row shows that kind.
  const icons = new Map<LogKind, SVGSVGElement>();
  let icon = Icon({ name: LogKindIcons.info, class: classes.kind });
  icons.set('info', icon);
  const el = (
    <ViewRow class={classes.row} icon={icon} label={label} description={description} data-kind="info" />
  ) as HTMLElement;
  let shown: LogEntry | undefined;
  return {
    el,
    show(e) {
      if (!e || e === shown) return;
      shown = e;
      setAttr(el, 'data-kind', e.kind);
      let next = icons.get(e.kind);
      if (!next) icons.set(e.kind, (next = Icon({ name: LogKindIcons[e.kind], class: classes.kind })));
      if (next !== icon) {
        icon.replaceWith(next);
        icon = next;
      }
      label.set(e.text);
      description.set(logDescription(e));
    },
  };
}
