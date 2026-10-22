/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { RecycledList, type VirtualSlot } from '../../components/virtual-list/recycled-list.js';
import { For, Show } from '../../dom/control-flow.js';
import { setAttr, setText } from '../../dom/writes.js';
import type { LogFilter, MessageKind, StatusMessage, NotificationsManager } from '../../managers/notifications.js';
import { observableState } from '../../reactive/observable-state.js';
import { derive } from '../../reactive/state.js';
import { classes } from './_generated/log-view.module.css.js';
import { KindIcons } from './canvas-chrome.js';

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

interface LogChoice {
  readonly id: LogFilter;
  readonly label: string;
  readonly icon: IconName;
}

/** The icons of the log-only sources the viewer writes; another source gets a tag. */
const SourceIcons: Readonly<Record<string, IconName>> = { Transforms: 'layers', Plugin: 'wrench' };

const passes = (m: StatusMessage, filter: LogFilter, query: string) => {
  if (filter === 'messages' && m.source !== undefined) return false;
  if (filter === 'problems' && m.kind !== 'warning' && m.kind !== 'error') return false;
  if (filter.startsWith('source:') && m.source !== filter.slice(7)) return false;
  return !query || m.text.toLowerCase().includes(query) || !!m.detail?.toLowerCase().includes(query);
};

/**
 * The Log tab: every message the status line showed and what was written to the log only (the plugin's transforms,
 * each with how long it took), newest first, with its kind, its detail and its time. Filtered by what wrote it or the
 * problems, and by text. The entries are a virtual list of one-line rows (perf audit 2, 11/P11): a long line fades at
 * its end and its tooltip has the rest. It fills its parent (a flex column: the panel's body, a half of History).
 */
export function LogView(props: { status: NotificationsManager }) {
  const filter = props.status.logFilter;
  const query = observableState('');
  // The log-only sources written, a filter each (the manager keeps them: no pass over the entries a flush).
  const choices = derive(props.status.logSources, (list): readonly LogChoice[] => [
    { id: 'all', label: 'All', icon: 'list' },
    { id: 'messages', label: 'Messages', icon: 'info' },
    ...list.map((s): LogChoice => ({ id: `source:${s}`, label: s, icon: SourceIcons[s] ?? 'tag' })),
    { id: 'problems', label: 'Problems', icon: 'triangle-alert' },
  ]);
  const shown = derive([props.status.logged, filter, query.asReadonly], (list, f, q): readonly StatusMessage[] => {
    const text = q.trim().toLowerCase();
    return f === 'all' && !text ? list : list.filter((m) => passes(m, f, text));
  });
  const empty = derive([shown, props.status.logged], (list, all) =>
    list.length ? undefined : all.length ? 'Nothing matches.' : 'No messages yet.',
  );
  // Not live (a log is, by default): each message is said once, by the status line's region (accessibility.md).
  return (
    <div class={classes.log}>
      {/* One line above the entries: the filters as icons, then the text filter. */}
      <div class={classes.filters}>
        <div class={classes.kinds} role="group" aria-label="Show">
          <For each={choices} by={(c) => c.id}>
            {(c) => (
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name={c.icon} size={15} />}
                title={c.label}
                pressed={derive(filter, (f) => f === c.id)}
                onClick={() => props.status.setLogFilter(c.id)}
              />
            )}
          </For>
        </div>
        <TextInput
          search
          size="sm"
          class={classes.search}
          placeholder="Filter"
          aria-label="Filter the log"
          value={query.asReadonly}
          onInput={(v) => query.set(v)}
        />
      </div>
      <Show when={empty}>{(text) => <p class={classes.empty}>{text}</p>}</Show>
      <RecycledList
        items={shown}
        itemKey={(m) => m.id}
        rowHeight={30}
        class={classes.entries}
        edge
        role="log"
        aria-label="Log"
        contentRef={(el) => el.setAttribute('aria-live', 'off')}
        slot={logRow}
      />
    </div>
  );
}

/** A recycled row of the Log: made once, shown one entry, then another as the list scrolls or grows. */
function logRow(): VirtualSlot<StatusMessage> {
  // An icon per kind, made the first time this row shows that kind.
  const icons = new Map<MessageKind, SVGSVGElement>();
  let icon: SVGSVGElement = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  const text = document.createTextNode('');
  const detail = document.createTextNode('');
  const at = document.createTextNode('');
  const stamp = (<time class={classes.time}>{at}</time>) as HTMLElement;
  const el = (
    <div class={classes.entry}>
      {icon}
      <span class={classes.text} data-tip-overflow="">
        {text}
      </span>
      <span class={classes.detail}>{detail}</span>
      {stamp}
    </div>
  ) as HTMLElement;
  let shown: StatusMessage | undefined;
  return {
    el,
    show(m) {
      if (!m || m === shown) return;
      shown = m;
      setAttr(el, 'data-kind', m.kind);
      setAttr(el, 'data-source', m.source);
      let next = icons.get(m.kind);
      if (!next) icons.set(m.kind, (next = Icon({ name: KindIcons[m.kind], size: 14, class: classes.kind })));
      if (next !== icon) {
        icon.replaceWith(next);
        icon = next;
      }
      setText(text, m.text);
      setText(detail, m.detail ?? '');
      setAttr(stamp, 'dateTime', m.time.toISOString());
      setText(at, time.format(m.time));
    },
  };
}
