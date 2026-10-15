/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, debounceTime, distinctUntilChanged, map, type Observable, shareReplay } from 'rxjs';
import type { CommandInvocation, InvocationState } from '../../commands/invocation.js';
import { formatChord, isMac } from '../../commands/keys.js';
import type { CommandManager, CommandMatch } from '../../commands/manager.js';
import { fuzzyScore, matchPositions } from '../../commands/match.js';
import type { CommandGroup } from '../../commands/types.js';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { For, Show } from '../../dom/control-flow.js';
import { rememberFocus, trapTab, uniqueId } from '../../dom/focus.js';
import type { Child } from '../../dom/jsx-types.js';
import { bind, onCleanup } from '../../dom/scope.js';
import type { FormFieldView } from '../../forms/flat-fields.js';
import { keyedSelector } from '../../reactive/keyed-selector.js';
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
  /** The assistant (when the spec turns AI on): typed words get an Ask row (design/agent.md §5). */
  assistant?: CommandMenuAssistant;
  /**
   * Opened for one command from elsewhere (a component's action): its form is the menu, so Esc (and ⌫ past its first
   * field) closes the menu rather than going back to the list.
   */
  direct?: boolean;
}

/**
 * What the command menu asks of the assistant: words read into one command (any-ui's `assist`: it proposes, nothing
 * runs until Enter), or the request handed to the agent; without a model, the setup.
 */
export interface CommandMenuAssistant {
  ready(): boolean;
  /** A command these words fit, with them read as its arguments (the filled sentence opens), if one does. */
  propose(text: string): { command: string; tail?: string } | undefined;
  /** Hands the request to the agent: a turn of its own. */
  ask(text: string): void;
  /** Opens the model setup. */
  setup(): void;
}

/** The Ask row: the assistant proposes a command (its sentence opens, filled) or takes the request as a turn. */
function askItem(
  props: CommandMenuProps,
  manager: CommandManager<any>,
  text: string,
  onPropose: (words: string) => void,
): CommandMenuItem {
  const assistant = props.assistant!;
  const words = text.trim();
  if (!assistant.ready())
    return {
      id: 'assistant-setup',
      label: 'Set up the assistant to ask',
      detail: `“${words}”`,
      icon: 'sparkles',
      run: () => assistant.setup(),
    };
  return {
    id: 'assistant-ask',
    label: 'Ask AI assistant',
    detail: `“${words}”`,
    icon: 'sparkles',
    run: () => {
      const proposal = assistant.propose(words);
      if (proposal) {
        onPropose(words);
        const r = manager.open(proposal.command, proposal.tail);
        // Proposed: the sentence opens filled in, and nothing runs until Enter.
        return r.kind === 'form';
      }
      assistant.ask(words);
    },
  };
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

/** Rows rendered at first (and added at a time): more than the list shows. */
const Screenful = 24;
/** An option's height (`.row`), for the scroll height of the rows not rendered yet. */
const RowPx = 28;

export function CommandMenu(props: CommandMenuProps) {
  const { manager } = props;
  // The words a sentence was proposed from by the assistant (shown over it), until it closes.
  const proposed = observableState<string | undefined>(undefined);
  bind(manager.active, (a) => {
    if (!a) proposed.set(undefined);
  });
  // A command that ran (or failed past its form) closes the command menu; its outcome goes to the status line.
  bind(manager.outcomes, (o) => {
    if (o.kind === 'done') props.onClose();
  });
  onCleanup(() => manager.close());
  // A modal dialog (README §3.7): Tab stays inside it, and closing it gives the focus back to what opened it.
  const restore = rememberFocus();
  let menu: HTMLElement | undefined;
  onCleanup(() => restore(menu));
  return (
    <div
      class={classes.scrim}
      data-keeps-card=""
      onPointerDown={(e) => {
        if (e.target === e.currentTarget) props.onClose();
      }}
    >
      <div
        class={classes.menu}
        role="dialog"
        aria-modal="true"
        aria-label="Commands"
        ref={(el) => {
          menu = el;
        }}
        onKeyDown={(e) => {
          trapTab(e, e.currentTarget);
          // Escape anywhere in the menu (a run button, a picker) closes it, and nothing behind it.
          if (e.key !== 'Escape' || e.defaultPrevented) return;
          e.preventDefault();
          e.stopPropagation();
          props.onClose();
        }}
      >
        <Show when={manager.active} fallback={() => <ListMode {...props} onPropose={(t) => proposed.set(t)} />}>
          {(invocation) => (
            <FormMode
              manager={manager}
              invocation={invocation}
              proposed={proposed.asReadonly}
              onClose={props.onClose}
              {...(props.direct ? { direct: true } : {})}
            />
          )}
        </Show>
      </div>
    </div>
  );
}

function ListMode(props: CommandMenuProps & { onPropose(words: string): void }) {
  const { manager } = props;
  const text = observableState(props.initial ?? '');
  const group = observableState<string | undefined>(props.group);
  const highlight = observableState(0);
  const message = observableState<string | undefined>(undefined);
  // The combobox's list and its options' ids: the input names the highlighted one (`aria-activedescendant`).
  const listId = uniqueId('msp-command-list');
  const optionId = (i: number) => `${listId}-${i}`;
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
    // One computation for its subscribers (the list, `unmatched`, `current`), not one each.
    shareReplay({ bufferSize: 1, refCount: true }),
  );
  let current: Row[] = [];
  const picks = () => current.filter((r): r is Pickable => r.kind !== 'caption');
  // The assistant takes the typed words (⌘↵ anywhere; ↵ when no command matches).
  const ask = () => {
    if (!props.assistant) return;
    if (!props.assistant.ready()) {
      props.assistant.setup();
      props.onClose();
      return;
    }
    if (!text.value.trim() || prefixOf(text.value)) return;
    const keep = askItem(props, manager, text.value, props.onPropose).run();
    if (!keep && !manager.active.value) props.onClose();
  };
  const unmatched = combineLatest([text.asReadonly, rows]).pipe(
    map(([t, r]) => !!t.trim() && !prefixOf(t) && !r.some((x) => x.kind !== 'caption')),
    distinctUntilChanged(),
  );
  bind(rows, (r) => {
    current = r;
    if (highlight.value >= picks().length) highlight.set(0);
  });
  // The rows rendered: a screenful first, more as the list scrolls or the highlight moves past them, so opening the menu
  // builds what shows rather than every command (styling.md §5, G8 perf #3). A spacer keeps the scroll height of the
  // rest. The highlighted option is always rendered (`aria-activedescendant` names it): its bind below runs before the
  // input's and the list's, so the list has grown to it by the time they read the highlight.
  const limit = observableState<{ of?: readonly Row[]; n: number }>({ n: Screenful });
  const shown = combineLatest([rows, limit.asReadonly]).pipe(
    map(([r, l]) => (r.length <= (l.of === r ? l.n : Screenful) ? r : r.slice(0, l.of === r ? l.n : Screenful))),
    distinctUntilChanged((a, b) => a.length === b.length && a[0] === b[0] && a[a.length - 1] === b[b.length - 1]),
  );
  const rest = combineLatest([rows, shown]).pipe(
    map(([r, s]) => (r.length - s.length) * RowPx),
    distinctUntilChanged(),
  );
  const grow = (n: number) => {
    const l = limit.value;
    const have = l.of === current ? l.n : Screenful;
    if (n > have && have < current.length) limit.set({ of: current, n: Math.max(n, have + Screenful) });
  };
  bind(highlight.asReadonly, (h) => {
    const at = current.findIndex((r) => r.kind !== 'caption' && r.index === h);
    if (at >= 0) grow(at + 1);
  });
  // A move of the highlight notifies the two rows it concerns, not every row.
  const highlighted = keyedSelector(highlight.asReadonly);
  const active = combineLatest([highlight.asReadonly, rows]).pipe(
    map(([h, r]) => (r.some((x) => x.kind !== 'caption') ? optionId(h) : undefined)),
    distinctUntilChanged(),
  );
  // What a search found, said once the typing pauses (the list itself is not live).
  const found = combineLatest([text.asReadonly, rows]).pipe(
    debounceTime(400),
    map(([t, r]) => {
      if (!t.trim()) return '';
      const n = r.filter((x) => x.kind !== 'caption').length;
      if (n) return n === 1 ? '1 result' : `${n} results`;
      return props.assistant ? 'No command matches: ask AI assistant.' : 'Nothing matches.';
    }),
    distinctUntilChanged(),
  );

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
      aria-controls={listId}
      aria-autocomplete="list"
      aria-activedescendant={active}
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
      <Show when={message}>
        {(m) => (
          <p class={classes.message} role="alert">
            {m}
          </p>
        )}
      </Show>
      <span class={classes.srOnly} role="status" aria-live="polite">
        {found}
      </span>
      <ScrollArea
        class={classes.list}
        contentClass={classes.listContent}
        id={listId}
        role="listbox"
        aria-label="Commands"
        viewportRef={(el) => {
          // Reads only (a scroll event, before anything is written): the rows down to two screens below the top.
          const more = () => grow(Math.ceil((el.scrollTop + 2 * el.clientHeight) / RowPx));
          el.addEventListener('scroll', more, { passive: true });
          onCleanup(() => el.removeEventListener('scroll', more));
        }}
      >
        <For
          each={shown}
          fallback={() => (
            <p class={classes.empty}>
              {props.assistant ? 'No command matches: ask AI assistant.' : 'Nothing matches.'}
            </p>
          )}
        >
          {(row) => {
            // A caption is not an option (a listbox owns options only): read as text, skipped by the highlight.
            if (row.kind === 'caption')
              return (
                <div class={classes.caption} role="presentation">
                  {row.label}
                </div>
              );
            const i = row.index;
            const id = optionId(i);
            const selected = highlighted.is(i);
            const hover = () => highlight.set(i);
            if (row.kind === 'group')
              return (
                <GroupRow
                  id={id}
                  manager={manager}
                  group={row.group}
                  selected={selected}
                  onPick={() => pick(row)}
                  onHover={hover}
                />
              );
            if (row.kind === 'item')
              return <ItemRow id={id} item={row.item} selected={selected} onPick={() => pick(row)} onHover={hover} />;
            return (
              <CommandRow
                id={id}
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
        <Show when={rest.pipe(map((h) => h > 0))}>
          {() => <div class={classes.rest} role="presentation" style={{ height: rest.pipe(map((h) => `${h}px`)) }} />}
        </Show>
      </ScrollArea>
      {props.assistant ? (
        <AskBar
          ready={props.assistant.ready()}
          text={text.asReadonly.pipe(map((t) => (prefixOf(t) ? '' : t.trim())))}
          unmatched={unmatched}
          onAsk={() => {
            // Not set up: the setup, whatever is typed. Set up: the words, or the input to type them in.
            if (!props.assistant?.ready()) {
              props.assistant?.setup();
              props.onClose();
            } else if (text.value.trim()) ask();
            else input.focus();
          }}
        />
      ) : null}
      <Footer
        hints={[
          ['↑↓', 'move'],
          ['↵', 'run'],
          ...(props.assistant?.ready() ? [[isMac() ? '⌘↵' : 'Ctrl+↵', 'ask'] as const] : []),
          ['→', 'open group'],
          ['?', 'help'],
          ['esc', 'close'],
        ]}
      />
    </>
  );
}

/**
 * The assistant, always at hand under the list (design/agent.md §5): it takes what is typed (⌘↵, or ↵ when no command
 * matches, which it then shows as the default); without a model it opens the setup.
 */
function AskBar(props: { ready: boolean; text: Observable<string>; unmatched: Observable<boolean>; onAsk(): void }) {
  if (!props.ready)
    return (
      // Not set up: a quiet way to the setup, not an offer that competes with the commands.
      <button type="button" class={classes.ask} data-setup="" onClick={() => props.onAsk()}>
        <Icon name="sparkles" size={15} class={classes.askIcon} />
        <span class={classes.askText}>Set up an AI assistant</span>
      </button>
    );
  return (
    <button
      type="button"
      class={classes.ask}
      data-default={props.unmatched.pipe(map((u) => (u ? '' : undefined)))}
      onClick={() => props.onAsk()}
    >
      <Icon name="sparkles" size={15} class={classes.askIcon} />
      <span class={classes.askLabel}>Ask AI assistant</span>
      <span class={classes.askText}>
        {props.text.pipe(map((t) => (t ? `“${t}”` : 'describe what you want, in your words')))}
      </span>
      <kbd class={classes.askKeys}>{props.unmatched.pipe(map((u) => (u ? '↵' : isMac() ? '⌘↵' : 'Ctrl+↵')))}</kbd>
    </button>
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
      // The Help card first, when the app has one: the canvas's bindings and every key, in one place.
      const card: CommandMenuItem[] = props.manager.get('help.show')
        ? [
            {
              id: 'help:card',
              label: 'Keyboard and mouse',
              detail: 'The Help card',
              icon: 'circle-help',
              run: () => {
                props.manager.open('help.show');
              },
            },
          ]
        : [];
      return [...card, ...prefixes, ...shortcuts].filter((i) => !q || i.label.toLowerCase().includes(q));
    },
  };
}

function GroupRow(props: {
  id: string;
  manager: CommandManager<any>;
  group: CommandGroup;
  selected: Observable<boolean>;
  onPick(): void;
  onHover(): void;
}) {
  const path = props.manager.groupPath(props.group.parent).map((g) => g.title);
  const el = (
    <div
      id={props.id}
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
  bind(props.selected, (s) => s && queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' })));
  return el;
}

function ItemRow(props: {
  id: string;
  item: CommandMenuItem;
  selected: Observable<boolean>;
  onPick(): void;
  onHover(): void;
}) {
  const { item } = props;
  const el = (
    <div
      id={props.id}
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
  bind(props.selected, (s) => s && queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' })));
  return el;
}

function CommandRow(props: {
  id: string;
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
      id={props.id}
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
  bind(props.selected, (s) => {
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
function FormMode(props: {
  manager: CommandManager<any>;
  invocation: CommandInvocation<any>;
  /** Proposed by the assistant from these words. */
  proposed?: Observable<string | undefined>;
  onClose(): void;
  /** Opened for this command alone: leaving the form closes the menu. */
  direct?: boolean;
}) {
  const { manager, invocation } = props;
  const state = invocation.state;
  const command = invocation.command;
  const field = state.pipe(
    map((s) => s.fields.find((f) => f.id === s.field)),
    distinctUntilChanged(),
  );
  const listId = uniqueId('msp-command-options');
  const optionId = (i: number) => `${listId}-${i}`;
  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-expanded={state.pipe(map((s) => s.completions.length > 0))}
      aria-controls={listId}
      aria-autocomplete="list"
      aria-activedescendant={state.pipe(
        map((s) => (s.pick >= 0 && s.pick < s.completions.length ? optionId(s.pick) : undefined)),
        distinctUntilChanged(),
      )}
      aria-label={field.pipe(map((f) => f?.label ?? command.title))}
      placeholder={field.pipe(map((f) => (f ? f.label : 'Enter to run')))}
      spellcheck={false}
      value={state.pipe(map((s) => s.text))}
      onInput={(e) => invocation.type(e.currentTarget.value)}
      onKeyDown={(e) => {
        if (e.key === 'Enter' && e.altKey && command.view) {
          e.preventDefault();
          goToView();
        } else if (e.key === 'Enter') {
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
          if (!invocation.back()) leave();
        } else if (e.key === 'Escape') {
          e.preventDefault();
          e.stopPropagation();
          leave();
        }
      }}
    />
  ) as HTMLInputElement;
  queueMicrotask(() => input.focus());
  // Leaving the form: back to the list, or, opened for this command alone, the menu closes.
  const leave = () => {
    manager.close();
    if (props.direct) props.onClose();
  };
  // To where the command lives in the UI: the form and the menu close, the place opens.
  const goToView = () => {
    const view = command.view;
    if (!view) return;
    manager.close();
    props.onClose();
    view.open(manager.context.value);
  };
  return (
    <>
      <Show when={props.proposed}>
        {(words) => (
          <p class={classes.proposed}>
            <Icon name="sparkles" size={13} class={classes.askIcon} />
            Proposed by the assistant from “{words}”: check it, then Enter to run
          </p>
        )}
      </Show>
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
        {command.view ? (
          <button type="button" class={classes.view} title={`Open ${command.view.title} (⌥↵)`} onClick={goToView}>
            {command.view.title}
            <Icon name="chevron-right" size={13} />
          </button>
        ) : null}
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
        id={listId}
        role="listbox"
        aria-label="Options"
      >
        <For each={state.pipe(map((s) => s.completions))}>
          {(c, i) => (
            <div
              id={i.pipe(map(optionId))}
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
          ...(command.view ? [['⌥↵', 'open in the UI'] as [string, string]] : []),
          ['⌫', 'back'],
          props.direct ? ['esc', 'close'] : ['esc', 'commands'],
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
