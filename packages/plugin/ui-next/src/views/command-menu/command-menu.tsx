/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import type { CommandInvocation, InvocationState } from '../../commands/invocation.js';
import { formatChord } from '../../commands/keys.js';
import type { CommandManager, CommandMatch } from '../../commands/manager.js';
import { fuzzyScore, matchPositions } from '../../commands/match.js';
import type { CommandGroup } from '../../commands/types.js';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import type { FormFieldView } from '../../forms/flat-fields.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/command-menu.module.css.js';

export interface CommandMenuProps {
  manager: CommandManager<any>;
  /** Closes the command menu (Esc on an empty list, a click outside, a command that ran). */
  onClose(): void;
  /** Text to start with ("select " from a selection button). */
  initial?: string;
  /** Prefixes that turn the menu into another list (`@` scene, `#` snapshots, `:` select); `?` help is built in. */
  prefixes?: readonly CommandMenuPrefix[];
  /** Opens inside this group (the Structure tab's Preset and Styles buttons open theirs). */
  group?: string;
}

/** One thing a prefix offers (a scene node for `@`, a snapshot for `#`, …). */
export interface CommandMenuItem {
  id: string;
  label: string;
  detail?: string;
  icon?: IconName;
  /** Does it; return true to keep the menu open (it opened a command's form, or typed another prefix). */
  run(): boolean | void;
}

/**
 * A prefix that turns the menu into another list: `@` goes to a scene node, `#` to a snapshot, `:` selects by residues
 * or an expression. The host provides them (it knows the scene); `?` (help) is built in.
 */
export interface CommandMenuPrefix {
  char: string;
  title: string;
  /** What to type after it. */
  placeholder?: string;
  items(text: string): readonly CommandMenuItem[];
}

/** A row of the list: a caption, a group to open, a command (or one of its ways in), or a prefix's item. */
type Row =
  | { kind: 'caption'; label: string }
  | { kind: 'group'; group: CommandGroup; index: number }
  | { kind: 'match'; match: CommandMatch<any>; tail?: string /** Its place among the pickable rows. */; index: number }
  | { kind: 'item'; item: CommandMenuItem; index: number };
type Pickable = Exclude<Row, { kind: 'caption' }>;

export function CommandMenu(props: CommandMenuProps) {
  const { manager } = props;
  // A command that ran (or failed past its form) closes the command menu; its outcome goes to the status line.
  const sub = manager.outcomes.subscribe((o) => {
    if (o.kind === 'done') props.onClose();
  });
  onCleanup(() => {
    sub.unsubscribe();
    manager.close();
  });
  return (
    <div
      class={classes.scrim}
      onPointerDown={(e) => {
        if (e.target === e.currentTarget) props.onClose();
      }}
    >
      <div class={classes.menu} role="dialog" aria-label="Commands">
        <Show when={manager.active} fallback={() => <ListMode {...props} />}>
          {(invocation) => <FormMode manager={manager} invocation={invocation} onClose={props.onClose} />}
        </Show>
      </div>
    </div>
  );
}

function ListMode(props: CommandMenuProps) {
  const { manager } = props;
  const text = observableState(props.initial ?? '');
  const group = observableState<string | undefined>(props.group);
  const highlight = observableState(0);
  const message = observableState<string | undefined>(undefined);
  const prefixes = [...(props.prefixes ?? []), helpPrefix(props, (t) => setText(t))];
  const prefixOf = (t: string) => prefixes.find((p) => t.startsWith(p.char));

  const rows = combineLatest([text.asReadonly, group.asReadonly, manager.commands, manager.groups]).pipe(
    map(([t, g]) => {
      const out: Row[] = [];
      let n = 0;
      const add = (row: Row) => {
        if (row.kind !== 'caption') (row as Pickable).index = n++;
        out.push(row);
      };
      const prefix = prefixOf(t);
      if (prefix) {
        for (const item of prefix.items(t.slice(prefix.char.length).trimStart())) add({ kind: 'item', item, index: 0 });
        return out;
      }
      if (!t.trim()) {
        // Browsing: a group's subgroups and commands; at the top, recent ones first, then the groups, then the
        // ungrouped commands by category.
        const { groups, commands } = manager.children(g);
        if (g === undefined) {
          const q = manager.query('');
          if (q.recent.length) {
            add({ kind: 'caption', label: 'Recent' });
            for (const m of q.recent) add({ kind: 'match', match: m, index: 0 });
          }
          if (groups.length) add({ kind: 'caption', label: 'Browse' });
        }
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
        const recent = new Set(manager.query('').recent.map((m) => m.command));
        const byCategory = new Map<string, CommandMatch<any>[]>();
        for (const m of commands) {
          if (recent.has(m.command)) continue;
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
      // Searching: everything at any depth (inside a group, its subtree), groups that match first.
      const q = manager.query(t);
      const within = (m: CommandMatch<any>) => g === undefined || manager.inGroup(m.command, g);
      if (q.direct && within(q.direct)) add({ kind: 'match', match: q.direct, tail: q.direct.tail, index: 0 });
      const groups = manager.groups.value.filter(
        (gr) =>
          (g === undefined || manager.groupPath(gr.id).some((x) => x.id === g)) &&
          gr.id !== g &&
          (fuzzyScore(gr.title, t) ?? -1) >= 0,
      );
      for (const gr of groups) add({ kind: 'group', group: gr, index: 0 });
      for (const m of q.matches)
        if (m.command !== q.direct?.command && within(m)) add({ kind: 'match', match: m, index: 0 });
      return out;
    }),
  );
  let current: Row[] = [];
  const picks = () => current.filter((r): r is Pickable => r.kind !== 'caption');
  rows.subscribe((r) => {
    current = r;
    if (highlight.value >= picks().length) highlight.set(0);
  });

  const setText = (t: string) => {
    text.set(t);
    input.value = t;
    highlight.set(0);
    message.set(undefined);
  };
  const enter = (id: string | undefined) => {
    group.set(id);
    setText('');
  };
  const pick = (row: Pickable) => {
    if (row.kind === 'group') return enter(row.group.id);
    if (row.kind === 'item') {
      const keep = row.item.run();
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

  const placeholder = combineLatest([text.asReadonly, group.asReadonly]).pipe(
    map(([t, g]) => {
      const prefix = prefixOf(t);
      if (prefix) return prefix.placeholder ?? prefix.title;
      return g === undefined ? 'Search commands, or type ? for help' : `Search in ${manager.groupPath(g).pop()?.title}`;
    }),
  );
  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-expanded="true"
      aria-controls="msp-command-menu-list"
      aria-label="Command"
      placeholder={placeholder}
      spellcheck={false}
      value={text}
      onInput={(e) => {
        text.set(e.currentTarget.value);
        highlight.set(0);
        message.set(undefined);
      }}
      onKeyDown={(e) => {
        const n = picks().length;
        const empty = !e.currentTarget.value;
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          if (n) highlight.set((highlight.value + (e.key === 'ArrowDown' ? 1 : n - 1)) % n);
        } else if (e.key === 'Enter') {
          e.preventDefault();
          const row = picks()[highlight.value];
          if (row) pick(row);
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

  return (
    <>
      <div class={classes.head}>
        <Icon name="search" size={15} class={classes.headIcon} />
        <For each={group.asReadonly.pipe(map((g) => manager.groupPath(g)))}>
          {(g) => (
            <button type="button" class={classes.crumb} onClick={() => enter(g.id)}>
              {g.title}
              <Icon name="chevron-right" size={12} />
            </button>
          )}
        </For>
        {input}
      </div>
      <Show when={message}>{(m) => <p class={classes.message}>{m}</p>}</Show>
      <ScrollArea
        class={classes.list}
        contentClass={classes.listContent}
        id="msp-command-menu-list"
        role="listbox"
        aria-label="Commands"
      >
        <For each={rows} fallback={() => <p class={classes.empty}>Nothing matches.</p>}>
          {(row) => {
            if (row.kind === 'caption') return <div class={classes.caption}>{row.label}</div>;
            const i = row.index;
            const selected = highlight.asReadonly.pipe(map((h) => h === i));
            const hover = () => highlight.set(i);
            if (row.kind === 'group')
              return (
                <GroupRow
                  manager={manager}
                  group={row.group}
                  selected={selected}
                  onPick={() => pick(row)}
                  onHover={hover}
                />
              );
            if (row.kind === 'item')
              return <ItemRow item={row.item} selected={selected} onPick={() => pick(row)} onHover={hover} />;
            return (
              <CommandRow
                row={row}
                path={relativePath(manager, row.match.command.group, group.value)}
                inGroup={group.value !== undefined}
                query={text.asReadonly}
                selected={selected}
                onPick={() => pick(row)}
                onHover={hover}
              />
            );
          }}
        </For>
      </ScrollArea>
      <Footer
        hints={[
          ['↑↓', 'move'],
          ['↵', 'run'],
          ['→', 'open group'],
          ['?', 'help'],
          ['esc', 'close'],
        ]}
      />
    </>
  );
}

/** A command's group path below the group being browsed (all of it at the top level). */
function relativePath(manager: CommandManager<any>, of: string | undefined, browsing: string | undefined): string[] {
  const path = manager.groupPath(of);
  const at = browsing === undefined ? -1 : path.findIndex((g) => g.id === browsing);
  return path.slice(at + 1).map((g) => g.title);
}

/** `?`: the prefixes (picking one types it), then every command with a shortcut. */
function helpPrefix(props: CommandMenuProps, type: (text: string) => void): CommandMenuPrefix {
  return {
    char: '?',
    title: 'Help',
    placeholder: 'Prefixes and shortcuts',
    items: (text) => {
      const q = text.toLowerCase();
      const prefixes: CommandMenuItem[] = (props.prefixes ?? []).map((p) => ({
        id: `prefix:${p.char}`,
        label: `${p.char}  ${p.title}`,
        detail: p.placeholder,
        icon: 'command',
        run: () => {
          type(p.char);
          return true;
        },
      }));
      const shortcuts: CommandMenuItem[] = props.manager.commands.value
        .filter((c) => c.shortcut)
        .map((c) => ({
          id: `keys:${c.id}`,
          label: c.title,
          detail: formatChord(c.shortcut!),
          icon: c.icon ?? 'keyboard',
          run: () => {
            props.manager.open(c.id);
          },
        }));
      return [...prefixes, ...shortcuts].filter((i) => !q || i.label.toLowerCase().includes(q));
    },
  };
}

function GroupRow(props: {
  manager: CommandManager<any>;
  group: CommandGroup;
  selected: Observable<boolean>;
  onPick(): void;
  onHover(): void;
}) {
  const path = props.manager.groupPath(props.group.parent).map((g) => g.title);
  const el = (
    <div
      class={classes.row}
      role="option"
      aria-selected={props.selected}
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowIcon}>
        <Icon name={props.group.icon ?? 'list-tree'} size={15} />
      </span>
      <span class={classes.rowTitle}>{props.group.title}</span>
      <span class={classes.rowMeta}>{path.join(' › ')}</span>
      <Icon name="chevron-right" size={14} class={classes.rowChevron} />
    </div>
  ) as HTMLElement;
  props.selected.subscribe((s) => s && queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' })));
  return el;
}

function ItemRow(props: { item: CommandMenuItem; selected: Observable<boolean>; onPick(): void; onHover(): void }) {
  const { item } = props;
  const el = (
    <div
      class={classes.row}
      role="option"
      aria-selected={props.selected}
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowIcon}>{item.icon ? <Icon name={item.icon} size={15} /> : null}</span>
      <span class={classes.rowTitle}>{item.label}</span>
      <span class={classes.rowMeta}>{item.detail}</span>
    </div>
  ) as HTMLElement;
  props.selected.subscribe((s) => s && queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' })));
  return el;
}

function CommandRow(props: {
  row: Extract<Row, { kind: 'match' }>;
  /** Its group path below the group being browsed ("Expression"), shown instead of the category when it has one. */
  path: readonly string[];
  /** Browsing a group: its captions say the categories, so a row without a deeper path shows nothing more. */
  inGroup?: boolean;
  query: import('rxjs').Observable<string>;
  selected: import('rxjs').Observable<boolean>;
  onPick(): void;
  onHover(): void;
}) {
  const { match, tail } = props.row;
  const c = match.command;
  const title = match.preset ? `${c.title.replace(/…$/, '')} ${match.preset.title}` : c.title;
  const el = (
    <div
      class={classes.row}
      role="option"
      aria-selected={props.selected}
      aria-disabled={match.blocked ? 'true' : undefined}
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowIcon}>{c.icon ? <Icon name={c.icon} size={15} /> : null}</span>
      <span class={classes.rowTitle}>
        {props.query.pipe(map((q) => highlight(title, tail ? '' : q)))}
        {tail ? <span class={classes.tail}>{tail}</span> : null}
      </span>
      <span class={classes.rowMeta}>
        {match.blocked ?? (props.path.length ? props.path.join(' › ') : props.inGroup ? '' : c.category)}
      </span>
      {c.shortcut ? <kbd class={classes.keys}>{formatChord(c.shortcut)}</kbd> : null}
    </div>
  ) as HTMLElement;
  props.selected.subscribe((s) => {
    if (s) queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' }));
  });
  return el;
}

/** The title with the query's letters marked. */
function highlight(text: string, query: string): Child {
  const at = new Set(matchPositions(text, query));
  if (!at.size) return text;
  const parts: Child[] = [];
  let run = '';
  let marked = false;
  const flush = () => {
    if (!run) return;
    parts.push(marked ? <mark>{run}</mark> : run);
    run = '';
  };
  for (let i = 0; i < text.length; i++) {
    const m = at.has(i);
    if (m !== marked) {
      flush();
      marked = m;
    }
    run += text[i];
  }
  flush();
  return <>{parts}</>;
}

/** A command with arguments: its sentence, the active field typed into the input, completions under it. */
function FormMode(props: { manager: CommandManager<any>; invocation: CommandInvocation<any>; onClose(): void }) {
  const { manager, invocation } = props;
  const state = invocation.state;
  const command = invocation.command;
  const field = state.pipe(
    map((s) => s.fields.find((f) => f.id === s.field)),
    distinctUntilChanged(),
  );
  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-expanded="true"
      aria-controls="msp-command-menu-completions"
      aria-label={field.pipe(map((f) => f?.label ?? command.title))}
      placeholder={field.pipe(map((f) => (f ? f.label : 'Enter to run')))}
      spellcheck={false}
      value={state.pipe(map((s) => s.text))}
      onInput={(e) => invocation.type(e.currentTarget.value)}
      onKeyDown={(e) => {
        if (e.key === 'Enter') {
          e.preventDefault();
          void invocation.enter();
        } else if (e.key === 'Tab') {
          e.preventDefault();
          void invocation.tab(e.shiftKey ? -1 : 1);
        } else if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          invocation.move(e.key === 'ArrowDown' ? 1 : -1);
        } else if (e.key === 'Backspace' && !e.currentTarget.value) {
          e.preventDefault();
          if (!invocation.back()) manager.close();
        } else if (e.key === 'Escape') {
          e.preventDefault();
          e.stopPropagation();
          manager.close();
        }
      }}
    />
  ) as HTMLInputElement;
  queueMicrotask(() => input.focus());
  return (
    <>
      <div class={classes.head}>
        <span class={classes.command}>
          {command.icon ? <Icon name={command.icon} size={15} /> : null}
          {command.title.replace(/…$/, '')}
        </span>
        {input}
      </div>
      <div class={classes.sentence} aria-label="Arguments">
        <For each={state.pipe(map((s) => s.fields))} key={(f) => f.id}>
          {(f) => (
            <FieldToken
              field={f}
              state={state}
              onFocus={() => {
                invocation.focus(f.id);
                input.focus();
              }}
            />
          )}
        </For>
        <span class={classes.spacer} />
        <button type="button" class={classes.run} onClick={() => void invocation.enter()}>
          {state.pipe(map((s) => s.action))}
          <kbd>↵</kbd>
        </button>
      </div>
      <Show when={state.pipe(map((s) => s.error))}>
        {(err) => (
          <p class={classes.error} role="alert">
            {err.message}
            {err.hint ? <span class={classes.hint}> {err.hint}</span> : null}
          </p>
        )}
      </Show>
      <ScrollArea
        class={classes.list}
        contentClass={classes.listContent}
        id="msp-command-menu-completions"
        role="listbox"
        aria-label="Options"
      >
        <For each={state.pipe(map((s) => s.completions))}>
          {(c, i) => (
            <div
              class={classes.row}
              role="option"
              aria-selected={combineLatest([state, i]).pipe(map(([s, j]) => s.pick === j))}
              onClick={() => invocation.accept(i.value)}
              onPointerMove={() => invocation.point(i.value)}
            >
              <span class={classes.rowIcon} />
              <span class={classes.rowTitle}>{c.label}</span>
              <span class={classes.rowMeta}>{c.detail ?? (c.value !== c.label ? c.value : '')}</span>
            </div>
          )}
        </For>
      </ScrollArea>
      <Show when={state.pipe(map((s) => s.formOnly.length > 0))}>
        {() => (
          <p class={classes.more}>
            {state.pipe(
              map((s) => `+ ${s.formOnly.length} more option${s.formOnly.length > 1 ? 's' : ''} at their defaults`),
            )}
          </p>
        )}
      </Show>
      <Footer
        hints={[
          ['tab', 'next field'],
          ['↵', 'run'],
          ['⌫', 'back'],
          ['esc', 'commands'],
        ]}
      />
    </>
  );
}

/** One argument in the sentence: its quiet before-word, then its value (or its name, until set). */
function FieldToken(props: {
  field: FormFieldView;
  state: import('rxjs').Observable<InvocationState>;
  onFocus(): void;
}) {
  const { field } = props;
  const value = props.state.pipe(
    map((s) => {
      const v = s.values[field.id];
      if (v === undefined || v === '') return undefined;
      return field.options?.find((o) => o.value === String(v))?.label ?? String(v);
    }),
  );
  return (
    <span class={classes.token}>
      {field.before ? <span class={classes.before}>{field.before}</span> : null}
      <button
        type="button"
        class={classes.value}
        data-active={props.state.pipe(map((s) => (s.field === field.id ? '' : undefined)))}
        data-missing={props.state.pipe(map((s) => (s.missing.includes(field.id) ? '' : undefined)))}
        data-empty={value.pipe(map((v) => (v === undefined ? '' : undefined)))}
        onClick={() => props.onFocus()}
      >
        {value.pipe(map((v) => v ?? field.label))}
      </button>
    </span>
  );
}

function Footer(props: { hints: readonly (readonly [string, string])[] }) {
  return (
    <div class={classes.footer}>
      {props.hints.map(([k, label]) => (
        <span>
          <kbd>{k}</kbd> {label}
        </span>
      ))}
    </div>
  );
}
