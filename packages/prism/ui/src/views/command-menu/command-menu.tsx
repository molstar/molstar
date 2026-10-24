/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandInvocation, InvocationState } from '../../commands/invocation.js';
import { formatChord, isMac } from '../../commands/keys.js';
import type { CommandManager, CommandMatch } from '../../commands/manager.js';
import { fuzzyScore, matchPositions } from '../../commands/match.js';
import type { CommandGroup, CommandMenuItem, CommandMenuPrefix, FormatCompletion } from '../../commands/types.js';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { keepRows } from '../../components/virtual-list/virtual-list.js';
import { For, Show } from '../../dom/control-flow.js';
import { deepActiveElement, focusables, rememberFocus, trapTab, uniqueId } from '../../dom/focus.js';
import type { Child } from '../../dom/jsx-types.js';
import { toNode } from '../../dom/render.js';
import { effect, flushFrame, frameState, onCleanup, type WatchSource, watch } from '../../dom/scope.js';
import { on, restartable } from '../../dom/scoped.js';
import { setAttr, setStyle, setText as writeText } from '../../dom/writes.js';
import type { FormFieldView } from '../../forms/flat-fields.js';
import { keyedSelector } from '../../reactive/keyed-selector.js';
import type { Listenable } from '../../reactive/listen.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
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
  /** The assistant (when the spec turns AI on): typed words get an Ask row (prism-assistant's _spec, "The Ask row"). */
  assistant?: CommandMenuAssistant;
  /**
   * Opened for one command from elsewhere (a component's action): its form is the menu, so Esc (and ⌫ past its first
   * field) closes the menu rather than going back to the list.
   */
  direct?: boolean;
  /** Commands' own forms, by command id (commands.md §5.1): drawn under the sentence in place of its option list. */
  forms?: Readonly<Record<string, CommandFormView>>;
}

/**
 * What a command's own form view gets (commands.md §5.1). It edits the open invocation (`set`, `setParam`, `focus`) and
 * runs it (`run`); the menu's input keeps the focus and its keys (Tab, Enter, ⌫, Esc), and `keys` sees each key first.
 * Tab past the sentence's last field goes into the view's own controls (those that take the focus), and back to the
 * input past them.
 */
export interface CommandFormViewProps {
  manager: CommandManager<any>;
  invocation: CommandInvocation<any>;
  /** Runs the command as Enter does. */
  run(): void;
  /** A key in the menu's input, before the menu acts on it: true takes it. */
  keys(handler: (e: KeyboardEvent) => boolean): void;
}

/** A command's own form in the command menu (commands.md §5.1), under the sentence. */
export interface CommandFormView {
  /** The fields it draws itself: while one of them is active, its options are not listed under the sentence. */
  fields?: readonly string[];
  render(props: CommandFormViewProps): Child;
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

// The prefixes' types are DOM-free (`commands/types.ts`): a host and its extensions make prefixes without the views.
export type { CommandMenuItem, CommandMenuPrefix };

/** A row of the list: a caption, a group to open, a command (or one of its ways in), or a prefix's item. */
/** A row's identity across lists: what it shows and where. */
const rowKey = (r: Row): string =>
  r.kind === 'caption'
    ? `caption:${r.label}`
    : r.kind === 'group'
      ? `group:${r.group.id}:${r.index}`
      : r.kind === 'item'
        ? `item:${r.item.label}:${r.index}`
        : `match:${r.match.command.id}:${r.match.preset?.id ?? ''}:${r.index}:${r.tail ?? ''}`;

/** The same row as before: what its view shows is the same (a query's marks are its own `watch`'s). */
const sameRow = (a: Row, b: Row): boolean => {
  if (a.kind !== b.kind) return false;
  if (a.kind === 'caption') return a.label === (b as typeof a).label;
  if (a.kind === 'group') return a.group === (b as typeof a).group && a.index === (b as typeof a).index;
  if (a.kind === 'item') return a.item === (b as typeof a).item && a.index === (b as typeof a).index;
  const m = b as typeof a;
  return (
    a.match.command === m.match.command &&
    a.match.preset === m.match.preset &&
    a.match.blocked === m.match.blocked &&
    a.index === m.index &&
    a.tail === m.tail
  );
};

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
  // Open, the menu wants every command: lazy ones (the Select group's per structure) are made now and follow changes
  // while it is open (commands.md §4, `registerLazy`).
  onCleanup(manager.want());
  // The words a sentence was proposed from by the assistant (shown over it), until it closes.
  const proposed = observableState<string | undefined>(undefined);
  effect(manager.active, (a) => {
    if (!a) proposed.set(undefined);
  });
  // A command that ran (or failed past its form) closes the command menu; its outcome goes to the status line.
  effect(manager.outcomes, (o) => {
    if (o.kind === 'done') props.onClose();
  });
  onCleanup(() => manager.close());
  // A modal dialog (accessibility.md): Tab stays inside it, and closing it gives the focus back to what opened it.
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
              {...(props.forms?.[invocation.command.id] ? { form: props.forms[invocation.command.id] } : {})}
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
  // What the list is filtered by (rendering.md §2.3, the frame lane): the typed text, once a frame. Keys typed faster
  // than frames (a held key, a fast typist, a paste in pieces) filter once, for the latest text; the input shows each.
  const query = frameState(text.value);
  // Before reading the rows (a key acting on them, a pick): the latest text filtered now.
  const settle = () => {
    if (query.latest !== query.value) flushFrame();
  };
  const group = observableState<string | undefined>(props.group);
  const highlight = observableState(0);
  const message = observableState<string | undefined>(undefined);
  // The combobox's list and its options' ids: the input names the highlighted one (`aria-activedescendant`).
  const listId = uniqueId('msp-command-list');
  const optionId = (i: number) => `${listId}-${i}`;
  const prefixes = [...(props.prefixes ?? []), helpPrefix(props, (t) => setText(t))];
  const prefixOf = (t: string) => prefixes.find((p) => t.startsWith(p.char));

  // Keyed by the group browsed too: a row shows its path below it (made with the row).
  const keepSame = keepRows<Row>((r) => `${group.value ?? ''}|${rowKey(r)}`, sameRow);
  // One computation for its listeners (the list, `unmatched`, `current`), not one each. A row the same as one before
  // (the same command at the same place) keeps its object: the list keeps its node, and a keystroke renders only the
  // marks in its title (the quick tier read 1.9 MB a keystroke before; hot.quick.ts). `keepSame` is stateful and a
  // derived value computes again when read with nothing listening (a subscription reads, then listens): the rows are
  // kept by what they were made from, so every reader gets the same array.
  let made: { t: string; g: string | undefined; c: unknown; gs: unknown; rows: Row[] } | undefined;
  const rows = derive([query, group.asReadonly, manager.commands, manager.groups], (t, g, c, gs) => {
    if (made && made.t === t && made.g === g && made.c === c && made.gs === gs) return made.rows;
    const r = keepSame(listRows(t, g));
    if (made) {
      made.t = t;
      made.g = g;
      made.c = c;
      made.gs = gs;
      made.rows = r;
    } else made = { t, g, c, gs, rows: r };
    return r;
  });
  function listRows(t: string, g: string | undefined): Row[] {
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
  }
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
  const unmatched = derive([query, rows], (t, r) => !!t.trim() && !prefixOf(t) && !r.some((x) => x.kind !== 'caption'));
  effect(rows, (r) => {
    current = r;
    if (highlight.value >= picks().length) highlight.set(0);
  });
  // The rows rendered: a screenful first, more as the list scrolls or the highlight moves past them, so opening the menu
  // builds what shows rather than every command (styling.md §5). A spacer keeps the scroll height of the
  // rest. The highlighted option is always rendered (`aria-activedescendant` names it): its effect below runs before the
  // input's and the list's watches, so the list has grown to it by the time they read the highlight.
  const limit = observableState<{ of?: readonly Row[]; n: number }>({ n: Screenful });
  const shown = derive(
    [rows, limit.asReadonly],
    (r, l) => (r.length <= (l.of === r ? l.n : Screenful) ? r : r.slice(0, l.of === r ? l.n : Screenful)),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  const rest = derive([rows, shown], (r, s) => (r.length - s.length) * RowPx);
  const grow = (n: number) => {
    const l = limit.value;
    const have = l.of === current ? l.n : Screenful;
    if (n > have && have < current.length) limit.set({ of: current, n: Math.max(n, have + Screenful) });
  };
  effect(highlight.asReadonly, (h) => {
    const at = current.findIndex((r) => r.kind !== 'caption' && r.index === h);
    if (at >= 0) grow(at + 1);
  });
  // A move of the highlight tells the two rows it concerns, not every row.
  const highlighted = keyedSelector(highlight.asReadonly);
  // The highlighted command is previewed and counted (commands.md §5): that row only; none once the menu closes.
  effect(
    derive([highlight.asReadonly, rows], (h, r) => {
      const row = r.find((x) => x.kind === 'match' && x.index === h);
      return row?.kind === 'match' && !row.match.preset && !row.tail && !row.match.blocked
        ? row.match.command.id
        : undefined;
    }),
    (id) => manager.browse(id),
  );
  onCleanup(() => manager.browse(undefined));
  // The shared argument of the group being browsed (Select's mode), shown once under the input.
  const shared = derive([group.asReadonly, manager.groups], (g) =>
    [...manager.groupPath(g)].reverse().find((x) => x.shared),
  );
  const active = derive([highlight.asReadonly, rows], (h, r) =>
    r.some((x) => x.kind !== 'caption') ? optionId(h) : undefined,
  );
  // What a search found, said once the typing pauses (the list itself is not live): read 400 ms after the last change
  // of the rows (they change with the query).
  const found = observableState('');
  const say = () => {
    const t = query.value;
    if (!t.trim()) return found.set('');
    let n = 0;
    for (const x of rows.value) if (x.kind !== 'caption') n++;
    if (n) return found.set(n === 1 ? '1 result' : `${n} results`);
    found.set(props.assistant ? 'No command matches: ask AI assistant.' : 'Nothing matches.');
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
  const pick = (row: Pickable) => {
    settle();
    if (row.kind === 'group') return enter(row.group.id);
    if (row.kind === 'item') {
      const keep = row.item.run();
      // A completion: the menu's text becomes it, and the list follows (commands.md §5).
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

  const placeholder = derive([text.asReadonly, group.asReadonly], (t, g) => {
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
  // The input as a function of the highlight, the group and the text: one watch.
  watch([active, placeholder, text.asReadonly], () => {
    setAttr(input, 'aria-activedescendant', active.value);
    setAttr(input, 'placeholder', placeholder.value);
    setValue(input, text.value);
  });
  const foundText = document.createTextNode('');
  watch([found.asReadonly], () => writeText(foundText, found.value));

  return (
    <>
      <div class={classes.head}>
        <Icon name="search" size={15} class={classes.headIcon} />
        <For each={derive(group.asReadonly, (g) => manager.groupPath(g))}>
          {(g) => (
            <button type="button" class={classes.crumb} onClick={() => enter(g.id)}>
              {g.title}
              <Icon name="chevron-right" size={12} />
            </button>
          )}
        </For>
        {input}
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
          // Reads only (a scroll event, before anything is written): the rows down to two screens below the top.
          const more = () => grow(Math.ceil((el.scrollTop + 2 * el.clientHeight) / RowPx));
          on(el, 'scroll', more, { passive: true });
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
            const option = {
              id: optionId(i),
              index: i,
              highlight: highlight.asReadonly,
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
                inGroup={group.value !== undefined}
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
        <AskBar
          ready={props.assistant.ready()}
          text={derive(text.asReadonly, (t) => (prefixOf(t) ? '' : t.trim()))}
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
 * The assistant, always at hand under the list (prism-assistant's _spec/README.md, "The Ask row"): it takes what is
 * typed (⌘↵, or ↵ when no command matches, which it then shows as the default); without a model it opens the setup.
 */
function AskBar(props: {
  ready: boolean;
  text: ReadonlyState<string>;
  unmatched: ReadonlyState<boolean>;
  onAsk(): void;
}) {
  if (!props.ready)
    return (
      // Not set up: a quiet way to the setup, not an offer that competes with the commands.
      <button type="button" class={classes.ask} data-setup="" onClick={() => props.onAsk()}>
        <Icon name="sparkles" size={15} class={classes.askIcon} />
        <span class={classes.askText}>Set up an AI assistant</span>
      </button>
    );
  const text = document.createTextNode('');
  const keys = document.createTextNode('');
  const el = (
    <button type="button" class={classes.ask} onClick={() => props.onAsk()}>
      <Icon name="sparkles" size={15} class={classes.askIcon} />
      <span class={classes.askLabel}>Ask AI assistant</span>
      <span class={classes.askText}>{text}</span>
      <kbd class={classes.askKeys}>{keys}</kbd>
    </button>
  ) as HTMLElement;
  watch([props.text, props.unmatched], () => {
    const t = props.text.value;
    const u = props.unmatched.value;
    setAttr(el, 'data-default', u);
    writeText(text, t ? `“${t}”` : 'describe what you want, in your words');
    writeText(keys, u ? '↵' : isMac() ? '⌘↵' : 'Ctrl+↵');
  });
  return el;
}

/**
 * A group's shared argument (commands.md §5: Select's mode) as a row of toggle buttons: set once, it goes to whichever
 * command in the group runs. Reachable with Tab (the menu keeps Tab inside); the input keeps the focus otherwise.
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
        .flatMap((c) => {
          const chord = props.manager.shortcut(c);
          return chord ? [{ c, chord }] : [];
        })
        .map(({ c, chord }) => ({
          id: `keys:${c.id}`,
          label: c.title,
          detail: formatChord(chord),
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

/** An input's value, written only when it differs (a write moves the caret). */
function setValue(el: HTMLInputElement, value: string) {
  if (el.value !== value) el.value = value;
}

/**
 * An option's highlight as its render writes it (rendering.md §2): `aria-selected`, and brought into view when it
 * becomes the highlighted one (after the flush's writes, a microtask later). Returns the next `was`.
 */
function writeSelected(el: HTMLElement, selected: boolean, was: boolean): boolean {
  setAttr(el, 'aria-selected', String(selected));
  if (selected && !was) queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' }));
  return selected;
}

/** What every option row is given: its id and index, the list's highlight, what it does. */
interface OptionProps {
  id: string;
  index: number;
  highlight: ReadonlyState<number>;
  /** Tells the row when it becomes or stops being the highlighted one (a move renders two rows, not all). */
  selected: Listenable;
  onPick(): void;
  onHover(): void;
}

function GroupRow(props: OptionProps & { manager: CommandManager<any>; group: CommandGroup }) {
  const path = props.manager.groupPath(props.group.parent).map((g) => g.title);
  const el = (
    <div
      id={props.id}
      class={classes.row}
      role="option"
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
  let was = false;
  watch([props.selected], () => {
    was = writeSelected(el, props.highlight.value === props.index, was);
  });
  return el;
}

function ItemRow(props: OptionProps & { item: CommandMenuItem }) {
  const { item } = props;
  const el = (
    <div
      id={props.id}
      class={classes.row}
      role="option"
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowIcon}>{item.icon ? <Icon name={item.icon} size={15} /> : null}</span>
      <span class={classes.rowTitle}>{item.label}</span>
      <span class={classes.rowMeta}>{item.detail}</span>
    </div>
  ) as HTMLElement;
  let was = false;
  watch([props.selected], () => {
    was = writeSelected(el, props.highlight.value === props.index, was);
  });
  return el;
}

/**
 * A command's row: its DOM made once, a `watch` render for what changes while it shows: the highlight and
 * (highlighted) how much it would do in place of its path; another for the query's marks in its title.
 */
function CommandRow(
  props: OptionProps & {
    manager: CommandManager<any>;
    row: Extract<Row, { kind: 'match' }>;
    /** Its group path below the group being browsed ("Expression"), shown instead of the category when it has one. */
    path: readonly string[];
    /** Browsing a group: its captions say the categories, so a row without a deeper path shows nothing more. */
    inGroup?: boolean;
    query: ReadonlyState<string>;
  },
) {
  const { manager, query } = props;
  const { match, tail } = props.row;
  const c = match.command;
  const title = match.preset ? `${c.title.replace(/…$/, '')} ${match.preset.title}` : c.title;
  const chord = manager.shortcut(c);
  const meta = match.blocked ?? (props.path.length ? props.path.join(' › ') : props.inGroup ? '' : c.category);
  const titleText = document.createElement('span');
  const metaText = document.createTextNode(meta);
  const el = (
    <div
      id={props.id}
      class={classes.row}
      role="option"
      aria-disabled={match.blocked ? 'true' : undefined}
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowIcon}>{c.icon ? <Icon name={c.icon} size={15} /> : null}</span>
      <span class={classes.rowTitle}>
        {titleText}
        {tail ? <span class={classes.tail}>{tail}</span> : null}
      </span>
      <span class={classes.rowMeta}>{metaText}</span>
      {chord ? <kbd class={classes.keys}>{formatChord(chord)}</kbd> : null}
    </div>
  ) as HTMLElement;
  let was = false;
  const sources: WatchSource[] = [props.selected];
  // Highlighted, how much it would do (`count`, commands.md §5) in place of its path or category.
  if (c.count && !match.preset) sources.push(manager.browsing);
  watch(sources, () => {
    was = writeSelected(el, props.highlight.value === props.index, was);
    if (c.count && !match.preset) {
      const b = manager.browsing.value;
      writeText(metaText, (b?.id === c.id ? b.count : undefined) ?? meta);
    }
  });
  // The query's marks in its own watch: the query is a frame state (the frame lane), the highlight a key's (at once).
  let marked: string | undefined;
  const mark = (q: string) => {
    if (q === marked) return;
    marked = q;
    // The marks are a few text nodes and `mark`s, no bindings: written whole (no marks: the text alone).
    const marks = highlight(title, q);
    if (typeof marks === 'string') titleText.textContent = marks;
    else titleText.replaceChildren(toNode(marks));
  };
  if (tail) mark('');
  else watch([query], () => mark(query.value));
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
  proposed?: ReadonlyState<string | undefined>;
  onClose(): void;
  /** Opened for this command alone: leaving the form closes the menu. */
  direct?: boolean;
  /** The command's own form, in place of the option list. */
  form?: CommandFormView;
}) {
  const { manager, invocation } = props;
  // The own form's key handler: it sees each key before the input acts on it.
  let formKeys: ((e: KeyboardEvent) => boolean) | undefined;
  // The active field is one the own form draws: its options are there, not in the list.
  const drawn = (s: InvocationState) => s.field !== undefined && !!props.form?.fields?.includes(s.field);
  const state = invocation.state;
  const command = invocation.command;
  const listId = uniqueId('msp-command-options');
  // The picked completion, a gate: the options listen to it, not to every change of the state.
  const pick = derive(state, (s) => s.pick);
  const optionId = (i: number) => `${listId}-${i}`;
  // The own form's element: Tab past the sentence's last field goes into its controls (commands.md §5.1).
  let formEl: HTMLElement | undefined;
  // Its Tab stops: not its tiles (`tabIndex` -1: the input's arrows move between them), which `focusables` keeps.
  const formControls = () => (formEl ? focusables(formEl).filter((el) => el.tabIndex >= 0) : []);
  /**
   * Tab past the last field (Shift+Tab before the first) with nothing typed to take: the focus goes to the own form's
   * first (last) control, its options the sentence does not show (a type's size, the uniform color). False: the
   * sentence's Tab.
   */
  const intoForm = (back: boolean): boolean => {
    if (!formEl) return false;
    const s = state.value;
    if (s.canAssist || (s.text.trim() && s.completions.length) || s.pick >= 0) return false;
    const i = s.fields.findIndex((f) => f.id === s.field);
    if (back ? i > 0 : i < s.fields.length - 1) return false;
    const controls = formControls();
    const to = back ? controls[controls.length - 1] : controls[0];
    if (!to) return false;
    to.focus();
    return true;
  };
  /** Tab past the own form's last control (Shift+Tab before its first): back to the input, the sentence's first (last) field. */
  const outOfForm = (e: KeyboardEvent) => {
    if (e.key !== 'Tab' || e.defaultPrevented || !formEl || e.altKey || e.metaKey || e.ctrlKey) return;
    const controls = formControls();
    const at = deepActiveElement(formEl.ownerDocument);
    if (at !== (e.shiftKey ? controls[0] : controls[controls.length - 1])) return;
    e.preventDefault();
    const fields = state.value.fields;
    const field = e.shiftKey ? fields[fields.length - 1] : fields[0];
    if (field) invocation.focus(field.id);
    input.focus();
  };
  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-controls={listId}
      aria-autocomplete="list"
      spellcheck={false}
      onInput={(e) => invocation.type(e.currentTarget.value)}
      onKeyDown={(e) => {
        if (formKeys?.(e)) {
          e.preventDefault();
          return;
        }
        if (e.key === 'Enter' && e.altKey && command.view) {
          e.preventDefault();
          goToView();
        } else if (e.key === 'Enter') {
          e.preventDefault();
          void invocation.enter();
        } else if (e.key === 'Tab') {
          e.preventDefault();
          if (!intoForm(e.shiftKey)) void invocation.tab(e.shiftKey ? -1 : 1);
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
  // The input and the run button's action as a function of the invocation's state: one watch.
  const action = document.createTextNode('');
  watch([state], () => {
    const s = state.value;
    writeText(action, s.action);
    const f = s.fields.find((x) => x.id === s.field);
    const listed = !drawn(s) && s.completions.length > 0;
    setAttr(input, 'aria-expanded', String(listed));
    setAttr(
      input,
      'aria-activedescendant',
      listed && s.pick >= 0 && s.pick < s.completions.length ? optionId(s.pick) : undefined,
    );
    setAttr(input, 'aria-label', f?.label ?? command.title);
    setAttr(input, 'placeholder', f ? f.label : 'Enter to run');
    setValue(input, s.text);
  });
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
        <For each={derive(state, (s) => s.fields)} by={(f) => f.id}>
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
          {action}
          <kbd>↵</kbd>
        </button>
      </div>
      <Show when={derive(state, (s) => s.error)}>
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
        <For each={derive(state, (s) => (drawn(s) ? NoCompletions : s.completions))}>
          {(c, i) => {
            const option = (
              <div
                class={classes.row}
                role="option"
                onClick={() => invocation.accept(i.value)}
                onPointerMove={() => invocation.point(i.value)}
              >
                <span class={classes.rowIcon} />
                <span class={classes.rowTitle}>{c.label}</span>
                <span class={classes.rowMeta}>{c.detail ?? (c.value !== c.label ? c.value : '')}</span>
              </div>
            ) as HTMLElement;
            // Its place in the list (the entry's index, a state) and whether it is the pick: only a move of the pick
            // or of the row wakes it, not a keystroke.
            watch([pick, i], () => {
              setAttr(option, 'id', optionId(i.value));
              setAttr(option, 'aria-selected', String(pick.value === i.value));
            });
            return option;
          }}
        </For>
      </ScrollArea>
      {props.form ? (
        <div
          class={classes.form}
          ref={(el: HTMLElement) => {
            formEl = el;
          }}
          onKeyDown={outOfForm}
        >
          {props.form.render({
            manager,
            invocation,
            run: () => void invocation.enter(),
            keys: (handler) => {
              formKeys = handler;
            },
          })}
        </div>
      ) : null}
      {/* An own form shows the options it offers; the rest stay at their defaults unsaid. */}
      <Show when={derive(state, (s) => !props.form && s.formOnly.length > 0)}>
        {() => {
          const text = document.createTextNode('');
          watch([state], () => {
            const n = state.value.formOnly.length;
            writeText(text, `+ ${n} more option${n > 1 ? 's' : ''} at their defaults`);
          });
          return <p class={classes.more}>{text}</p>;
        }}
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

const NoCompletions: readonly FormatCompletion[] = [];

/** One argument in the sentence: its quiet before-word, then its value (or its name, until set). */
function FieldToken(props: { field: FormFieldView; state: ReadonlyState<InvocationState>; onFocus(): void }) {
  const { field } = props;
  // What decides its look, each a gate over the invocation's state: the token renders when one of them changes.
  const value = derive(props.state, (s) => {
    const v = s.values[field.id];
    if (v === undefined || v === '') return undefined;
    return field.options?.find((o) => o.value === String(v))?.label ?? String(v);
  });
  const active = derive(props.state, (s) => s.field === field.id);
  const missing = derive(props.state, (s) => s.missing.includes(field.id));
  const text = document.createTextNode('');
  const button = (
    <button type="button" class={classes.value} onClick={() => props.onFocus()}>
      {text}
    </button>
  ) as HTMLElement;
  watch([value, active, missing], () => {
    const v = value.value;
    setAttr(button, 'data-active', active.value);
    setAttr(button, 'data-missing', missing.value);
    setAttr(button, 'data-empty', v === undefined);
    writeText(text, v ?? field.label);
  });
  return (
    <span class={classes.token}>
      {field.before ? <span class={classes.before}>{field.before}</span> : null}
      {button}
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
