/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The menu's list (commands §5.1–5.2): browsing, searching, prefixes and help, the Ask row.

import { derive, KeyedSelector, type ReadonlyState } from '@molstar/core/reactive';
import { formatChord } from '../../commands/keys.js';
import type { CommandManager } from '../../commands/manager.js';
import { fuzzyScore } from '../../commands/match.js';
import type { CommandMatch } from '../../commands/query.js';
import type { CommandGroup, CommandMenuItem, CommandMenuPrefix } from '../../commands/types.js';
import {
  effect,
  For,
  flushFrame,
  frameState,
  on,
  onCleanup,
  restartable,
  Show,
  setAttr,
  setStyle,
  setText as writeText,
  uniqueId,
  watch,
} from '../../dom/index.js';
import { state } from '../../state.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { keepRows } from '../virtual-list/virtual-list.js';
import { classes } from './_generated/command-menu.module.css.js';
import { Footer } from './footer.js';
import { CommandRow, GroupRow, ItemRow, type Pickable, type Row, RowPx, rowKey, sameRow } from './rows.js';
import type { CommandMenuAssistant, CommandMenuProps } from './types.js';

/** Rows rendered at first, and added at a time as the list scrolls or the highlight goes past them. */
const Screenful = 24;

export function ListMode(props: CommandMenuProps & { onPropose(words: string): void }) {
  const { manager } = props;
  const text = state(props.initial ?? '');
  // What the list is filtered by: the typed text, once a frame (dom §3); the input shows each key.
  const query = frameState(text.value);
  const settle = () => {
    if (query.latest !== query.value) flushFrame();
  };
  const group = state<string | undefined>(props.group);
  const highlight = state(0);
  const message = state<string | undefined>(undefined);
  const listId = uniqueId('msui-command-list');
  const optionId = (i: number) => `${listId}-${i}`;
  const prefixes = [...(props.prefixes ?? []), helpPrefix(props, (t) => setText(t))];
  const prefixOf = (t: string) => prefixes.find((p) => t.startsWith(p.char));

  // A row the same as before (the same command at the same place) keeps its object: the list keeps its node. Keyed by
  // the group browsed and by browsing or searching too: a row's path and meta are made with it. Rows are kept by what they were made from, so every reader of
  // the derived list gets the same array (`keepSame` is stateful).
  const keepSame = keepRows<Row>((r) => `${group.value ?? ''}|${query.value.trim() ? 's' : 'b'}|${rowKey(r)}`, sameRow);
  let made: { t: string; g: string | undefined; c: unknown; gs: unknown; rows: Row[] } | undefined;
  const rows = derive([query, group, manager.commands, manager.groups], (t, g, c, gs) => {
    if (made && made.t === t && made.g === g && made.c === c && made.gs === gs) return made.rows;
    made = { t, g, c, gs, rows: keepSame(listRows(manager, prefixOf, t, g)) };
    return made.rows;
  });
  let current: readonly Row[] = [];
  const picks = () => current.filter((r): r is Pickable => r.kind !== 'caption');
  effect(rows, (r) => {
    current = r;
    if (highlight.value >= picks().length) highlight.set(0);
  });
  const unmatched = derive([query, rows], (t, r) => !!t.trim() && !prefixOf(t) && !r.some((x) => x.kind !== 'caption'));

  // A screenful first, more as the list scrolls or the highlight moves past them; a spacer keeps the scroll height of
  // the rest. The highlighted option is always rendered (`aria-activedescendant` names it).
  const limit = state<{ of?: readonly Row[]; n: number }>({ n: Screenful });
  const shown = derive(
    [rows, limit],
    (r, l) => {
      const n = l.of === r ? l.n : Screenful;
      return r.length <= n ? r : r.slice(0, n);
    },
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  const rest = derive([rows, shown], (r, s) => (r.length - s.length) * RowPx);
  const grow = (n: number) => {
    const l = limit.value;
    const have = l.of === current ? l.n : Screenful;
    if (n > have && have < current.length) limit.set({ of: current, n: Math.max(n, have + Screenful) });
  };
  effect(highlight, (h) => {
    const at = current.findIndex((r) => r.kind !== 'caption' && r.index === h);
    if (at >= 0) grow(at + 1);
  });
  // A move of the highlight tells the two rows it concerns, not every row.
  const highlighted = new KeyedSelector(highlight.asReadonly);
  // The highlighted command is previewed and counted (commands §4.5); none once the menu closes.
  effect(
    derive([highlight, rows], (h, r) => {
      const row = r.find((x) => x.kind === 'match' && x.index === h);
      return row?.kind === 'match' && !row.match.preset && !row.tail && !row.match.blocked
        ? row.match.command.id
        : undefined;
    }),
    (id) => manager.browse(id),
  );
  onCleanup(() => manager.browse(undefined));
  const shared = derive([group, manager.groups], (g) => [...manager.groupPath(g)].reverse().find((x) => x.shared));
  const active = derive([highlight, rows], (h, r) => (r.some((x) => x.kind !== 'caption') ? optionId(h) : undefined));

  // What a search found, said once the typing pauses (the list itself is not live).
  const found = state('');
  const say = () => {
    const t = query.value;
    if (!t.trim()) return found.set('');
    let n = 0;
    for (const x of rows.value) if (x.kind !== 'caption') n++;
    if (n) return found.set(n === 1 ? '1 result' : `${n} results`);
    found.set(props.assistant ? 'No command matches: ask the assistant.' : 'Nothing matches.');
  };
  const foundTimer = restartable(say);
  effect(rows, () => foundTimer.start(400));

  const setText = (t: string) => {
    text.set(t);
    query.offer(t);
    settle();
    input.value = t;
    highlight.set(0);
    message.set(undefined);
  };
  const enter = (id: string | undefined) => {
    group.set(id);
    setText('');
  };
  const ask = () => {
    const assistant = props.assistant;
    if (!assistant) return;
    if (!assistant.ready()) {
      assistant.setup();
      props.onClose();
      return;
    }
    const words = text.value.trim();
    if (!words || prefixOf(text.value)) return;
    const proposal = assistant.propose(words);
    if (proposal) {
      props.onPropose(words);
      // Proposed: the sentence opens filled in, and nothing runs until Enter.
      const r = manager.open(proposal.command, proposal.tail);
      if (r.kind === 'blocked') message.set(r.reason);
      else if (r.kind === 'ran') props.onClose();
      return;
    }
    assistant.ask(words);
    props.onClose();
  };
  const pick = (row: Pickable) => {
    settle();
    if (row.kind === 'group') return enter(row.group.id);
    if (row.kind === 'item') {
      const keep = row.item.run();
      // A completion: the menu's text becomes it, and the list follows.
      if (keep && typeof keep === 'object') return setText(keep.text);
      if (!keep && !manager.active.value) props.onClose();
      return;
    }
    if (row.match.blocked) return message.set(row.match.blocked);
    const r = manager.open(row.match.command.id, row.tail, row.match.preset?.id);
    if (r.kind === 'blocked') message.set(r.reason);
    else if (r.kind === 'ran') props.onClose();
  };
  const up = () => {
    const path = manager.groupPath(group.value);
    if (!path.length) return false;
    enter(path[path.length - 2]?.id);
    return true;
  };

  const placeholder = derive([text, group], (t, g) => {
    const prefix = prefixOf(t);
    if (prefix) return prefix.placeholder ?? prefix.title;
    return g === undefined ? 'Search commands, or type ? for help' : `Search in ${manager.groupPath(g).pop()?.title}`;
  });
  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-expanded="true"
      aria-controls={listId}
      aria-autocomplete="list"
      aria-label="Command"
      autocomplete="off"
      spellcheck={false}
      onInput={(e) => {
        text.set(e.currentTarget.value);
        query.offer(e.currentTarget.value);
        highlight.set(0);
        message.set(undefined);
      }}
      onKeyDown={(e) => {
        settle();
        const n = picks().length;
        const empty = !e.currentTarget.value;
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          if (n) highlight.set((highlight.value + (e.key === 'ArrowDown' ? 1 : n - 1)) % n);
        } else if (e.key === 'Enter') {
          e.preventDefault();
          const row = picks()[highlight.value];
          if (e.metaKey || e.ctrlKey) ask();
          else if (row) pick(row);
          else ask();
        } else if (e.key === 'ArrowRight' && empty) {
          const row = picks()[highlight.value];
          if (row?.kind === 'group') {
            e.preventDefault();
            enter(row.group.id);
          }
        } else if ((e.key === 'Backspace' || e.key === 'ArrowLeft') && empty) {
          if (up()) e.preventDefault();
        } else if (e.key === 'Escape') {
          e.preventDefault();
          e.stopPropagation();
          props.onClose();
        }
      }}
    />
  ) as HTMLInputElement;
  queueMicrotask(() => input.focus());
  // The input as a function of the highlight, the placeholder and the text: one watch.
  watch([active, placeholder, text], () => {
    setAttr(input, 'aria-activedescendant', active.value);
    setAttr(input, 'placeholder', placeholder.value);
    if (input.value !== text.value) input.value = text.value;
  });
  const foundText = document.createTextNode('');
  watch([found], () => writeText(foundText, found.value));

  return (
    <>
      <div class={classes.head}>
        <span class={classes.headSlot}>
          <Icon name="search" />
        </span>
        <span class={classes.headLine}>
          <For each={derive(group, (g) => manager.groupPath(g))}>
            {(g) => (
              <button type="button" class={classes.crumb} onClick={() => enter(g.id)}>
                {g.title}
                <Icon name="chevron-right" />
              </button>
            )}
          </For>
          {input}
        </span>
      </div>
      <Show when={shared}>{(owner) => <SharedChoice manager={manager} group={owner} />}</Show>
      <Show when={message}>
        {(m) => (
          <p class={classes.message} role="alert">
            {m}
          </p>
        )}
      </Show>
      <span class={classes.srOnly} role="status" aria-live="polite">
        {foundText}
      </span>
      <ScrollArea
        class={classes.list}
        contentClass={classes.listContent}
        id={listId}
        role="listbox"
        aria-label="Commands"
        viewportRef={(el) => {
          // Reads only (a scroll event): the rows down to two screens below the top.
          on(el, 'scroll', () => grow(Math.ceil((el.scrollTop + 2 * el.clientHeight) / RowPx)), { passive: true });
        }}
      >
        <For
          each={shown}
          fallback={() => (
            <p class={classes.empty}>
              {props.assistant ? 'No command matches: ask the assistant.' : 'Nothing matches.'}
            </p>
          )}
        >
          {(row) => {
            // A caption is not an option (a listbox owns options only): text, skipped by the highlight.
            if (row.kind === 'caption')
              return (
                <div class={classes.caption} role="presentation">
                  {row.label}
                </div>
              );
            const i = row.index;
            const option = {
              id: optionId(i),
              selected: highlighted.key(i),
              onPick: () => pick(row),
              onHover: () => highlight.set(i),
            };
            if (row.kind === 'group') return <GroupRow {...option} manager={manager} group={row.group} />;
            if (row.kind === 'item') return <ItemRow {...option} item={row.item} />;
            return (
              <CommandRow
                {...option}
                manager={manager}
                row={row}
                path={relativePath(manager, row.match.command.group, group.value)}
                inGroup={group.value !== undefined || !query.value.trim()}
                query={query}
              />
            );
          }}
        </For>
        <Show when={derive(rest, (h) => h > 0)}>
          {() => {
            const spacer = (<div class={classes.rest} role="presentation" />) as HTMLElement;
            watch([rest], () => setStyle(spacer, 'height', `${rest.value}px`));
            return spacer;
          }}
        </Show>
      </ScrollArea>
      {props.assistant ? (
        <AskRow
          assistant={props.assistant}
          text={derive(text, (t) => (prefixOf(t) ? '' : t.trim()))}
          unmatched={unmatched}
          onAsk={() => {
            if (!props.assistant?.ready() || text.value.trim()) ask();
            else input.focus();
          }}
        />
      ) : null}
      <Footer
        hints={[
          ['↑↓', 'move'],
          ['↵', 'run'],
          ...(props.assistant?.ready() ? [[formatChord({ key: 'Enter', mod: true }), 'ask'] as const] : []),
          ['→', 'open group'],
          ['?', 'help'],
          ['esc', 'close'],
        ]}
      />
    </>
  );
}

/** The rows for a text in a group (commands §5.1): a prefix's items, browsing, or a search. */
function listRows(
  manager: CommandManager<any>,
  prefixOf: (t: string) => CommandMenuPrefix | undefined,
  t: string,
  g: string | undefined,
): Row[] {
  const out: Row[] = [];
  let n = 0;
  const add = (row: Row) => {
    if (row.kind !== 'caption') row.index = n++;
    out.push(row);
  };
  const prefix = prefixOf(t);
  if (prefix) {
    for (const item of prefix.items(t.slice(prefix.char.length).trimStart())) add({ kind: 'item', item, index: 0 });
    return out;
  }
  if (!t.trim()) {
    // Browsing: a group's subgroups and commands; at the top, recent ones first, then the groups, then the ungrouped
    // commands by category.
    const { groups, commands } = manager.children(g);
    const recent = g === undefined ? manager.query('').recent : [];
    if (recent.length) {
      add({ kind: 'caption', label: 'Recent' });
      for (const m of recent) add({ kind: 'match', match: m, index: 0 });
    }
    if (g === undefined && groups.length) add({ kind: 'caption', label: 'Browse' });
    for (const gr of groups) add({ kind: 'group', group: gr, index: 0 });
    if (g !== undefined) {
      // Inside a group: its commands under their categories' captions when there is more than one.
      const categories = new Set(commands.map((m) => m.command.category));
      let last: string | undefined;
      for (const m of commands) {
        if (categories.size > 1 && m.command.category !== last) {
          last = m.command.category;
          add({ kind: 'caption', label: last });
        }
        add({ kind: 'match', match: m, index: 0 });
      }
      return out;
    }
    const seen = new Set(recent.map((m) => m.command));
    const byCategory = new Map<string, CommandMatch<any>[]>();
    for (const m of commands) {
      if (seen.has(m.command)) continue;
      const list = byCategory.get(m.command.category);
      if (list) list.push(m);
      else byCategory.set(m.command.category, [m]);
    }
    for (const [category, matches] of byCategory) {
      add({ kind: 'caption', label: category });
      for (const m of matches) add({ kind: 'match', match: m, index: 0 });
    }
    return out;
  }
  // Searching: everything at any depth (inside a group, its subtree), a direct alias first, then matching groups.
  const q = manager.query(t);
  const within = (m: CommandMatch<any>) => g === undefined || manager.inGroup(m.command, g);
  if (q.direct && within(q.direct)) add({ kind: 'match', match: q.direct, tail: q.direct.tail, index: 0 });
  for (const gr of manager.groups.value) {
    if (gr.id === g || (fuzzyScore(gr.title, t) ?? -1) < 0) continue;
    if (g !== undefined && !manager.groupPath(gr.id).some((x) => x.id === g)) continue;
    add({ kind: 'group', group: gr, index: 0 });
  }
  for (const m of q.matches)
    if (m.command !== q.direct?.command && within(m)) add({ kind: 'match', match: m, index: 0 });
  return out;
}

/** A command's group path below the group being browsed (all of it at the top level). */
function relativePath(manager: CommandManager<any>, of: string | undefined, browsing: string | undefined): string[] {
  const path = manager.groupPath(of);
  const at = browsing === undefined ? -1 : path.findIndex((g) => g.id === browsing);
  return path.slice(at + 1).map((g) => g.title);
}

/**
 * The assistant under the list: it takes what is typed (⌘↵, or ↵ when no command matches: then it is the default);
 * without a model, a quiet way to the setup.
 */
function AskRow(props: {
  assistant: CommandMenuAssistant;
  text: ReadonlyState<string>;
  unmatched: ReadonlyState<boolean>;
  onAsk(): void;
}) {
  if (!props.assistant.ready())
    return (
      <button type="button" class={classes.ask} data-setup="" onClick={() => props.onAsk()}>
        <span class={classes.askSlot}>
          <Icon name="sparkles" />
        </span>
        <span class={classes.askText}>Set up an assistant to ask in your words</span>
      </button>
    );
  const words = document.createTextNode('');
  const keys = document.createTextNode('');
  const el = (
    <button type="button" class={classes.ask} onClick={() => props.onAsk()}>
      <span class={classes.askSlot}>
        <Icon name="sparkles" />
      </span>
      <span class={classes.askLabel}>Ask the assistant</span>
      <span class={classes.askText}>{words}</span>
      <kbd class={classes.askKeys}>{keys}</kbd>
    </button>
  ) as HTMLElement;
  watch([props.text, props.unmatched], () => {
    const t = props.text.value;
    const u = props.unmatched.value;
    setAttr(el, 'data-default', u);
    writeText(words, t ? `“${t}”` : 'describe what you want, in your words');
    writeText(keys, u ? '↵' : formatChord({ key: 'Enter', mod: true }));
  });
  return el;
}

/**
 * A group's shared argument (commands §2.3) as toggle buttons under the input: set once, it goes to whichever command
 * in the group runs. Tab reaches them; the input keeps the focus otherwise.
 */
function SharedChoice(props: { manager: CommandManager<any>; group: CommandGroup }) {
  const { manager, group } = props;
  const shared = group.shared!;
  const value = derive(manager.shared, () => manager.sharedValue(group.id));
  const buttons: HTMLElement[] = [];
  const el = (
    <div class={classes.shared} role="group" aria-label={`${group.title}: ${shared.title}`}>
      <span class={classes.sharedTitle}>{shared.title}</span>
      {shared.options.map(([v, label]) => (
        <button
          type="button"
          class={classes.sharedOption}
          ref={(b: HTMLElement) => buttons.push(b)}
          onClick={() => manager.setShared(group.id, v)}
        >
          {label}
        </button>
      ))}
    </div>
  );
  watch([value], () => {
    for (let i = 0; i < buttons.length; i++)
      setAttr(buttons[i], 'aria-pressed', value.value === shared.options[i][0] ? 'true' : 'false');
  });
  return el;
}

/** `?`: the prefixes (picking one types it), then every command with a key. */
function helpPrefix(props: CommandMenuProps, type: (text: string) => void): CommandMenuPrefix {
  return {
    char: '?',
    title: 'Help',
    placeholder: 'Prefixes and keys',
    items: (text) => {
      const q = text.toLowerCase();
      const prefixes: CommandMenuItem[] = (props.prefixes ?? []).map((p) => ({
        id: `prefix:${p.char}`,
        label: `${p.char.trim()}  ${p.title}`,
        ...(p.placeholder ? { detail: p.placeholder } : {}),
        icon: 'command',
        run: () => {
          type(p.char);
          return true;
        },
      }));
      const keys: CommandMenuItem[] = [];
      for (const c of props.manager.commands.value) {
        const chord = props.manager.shortcut(c);
        if (!chord) continue;
        keys.push({
          id: `keys:${c.id}`,
          label: c.title,
          detail: formatChord(chord),
          icon: c.icon ?? 'keyboard',
          run: () => {
            const r = props.manager.open(c.id);
            return r.kind === 'form';
          },
        });
      }
      return [...prefixes, ...keys].filter((i) => !q || i.label.toLowerCase().includes(q));
    },
  };
}
