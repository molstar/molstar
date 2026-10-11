/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import type { CommandInvocation, InvocationState } from '../../commands/invocation.js';
import { formatChord } from '../../commands/keys.js';
import type { CommandManager, CommandMatch } from '../../commands/manager.js';
import { matchPositions } from '../../commands/match.js';
import { Icon } from '../../components/icon/icon.js';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import type { FormFieldView } from '../../forms/flat-fields.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/palette.module.css.js';

export interface CommandPaletteProps {
  manager: CommandManager<any>;
  /** Closes the palette (Esc on an empty list, a click outside, a command that ran). */
  onClose(): void;
  /** Text to start with ("select " from a selection button). */
  initial?: string;
}

/** A row of the list: a section caption, or a command (or one of its ways in). */
type Row =
  | { kind: 'caption'; label: string }
  | { kind: 'match'; match: CommandMatch<any>; tail?: string /** Its place among the pickable rows. */; index: number };

/**
 * The command palette: everything the plugin can do, by name. A list while searching (recent commands first on an
 * empty input; matched letters marked; category, keys, and why a command cannot run now); a command with arguments
 * reads as a sentence under the input, its active field typed into the same input, completions below. Enter runs,
 * Tab moves between fields, Backspace on an empty field goes back, Esc steps out (form, then palette).
 */
export function CommandPalette(props: CommandPaletteProps) {
  const { manager } = props;
  // A command that ran (or failed past its form) closes the palette; its outcome goes to the status line.
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
      <div class={classes.palette} role="dialog" aria-label="Commands">
        <Show when={manager.active} fallback={() => <ListMode {...props} />}>
          {(invocation) => <FormMode manager={manager} invocation={invocation} onClose={props.onClose} />}
        </Show>
      </div>
    </div>
  );
}

function ListMode(props: CommandPaletteProps) {
  const { manager } = props;
  const text = observableState(props.initial ?? '');
  const highlight = observableState(0);
  const message = observableState<string | undefined>(undefined);
  const rows = combineLatest([text.asReadonly, manager.commands]).pipe(
    map(([t]) => {
      const q = manager.query(t);
      const out: Row[] = [];
      let n = 0;
      const add = (match: CommandMatch<any>, tail?: string) => out.push({ kind: 'match', match, tail, index: n++ });
      if (q.direct) add(q.direct, q.direct.tail);
      if (q.recent.length) {
        out.push({ kind: 'caption', label: 'Recent' });
        for (const m of q.recent) add(m);
      }
      if (!t.trim()) {
        // Nothing typed: everything, by category (in the order categories came), under captions.
        const byCategory = new Map<string, CommandMatch<any>[]>();
        for (const m of q.matches) {
          const list = byCategory.get(m.command.category);
          if (list) list.push(m);
          else byCategory.set(m.command.category, [m]);
        }
        for (const [category, matches] of byCategory) {
          out.push({ kind: 'caption', label: category });
          for (const m of matches) add(m);
        }
        return out;
      }
      for (const m of q.matches) if (m.command !== q.direct?.command) add(m);
      return out;
    }),
  );
  let current: Row[] = [];
  const picks = () => current.filter((r): r is Extract<Row, { kind: 'match' }> => r.kind === 'match');
  rows.subscribe((r) => {
    current = r;
    if (highlight.value >= picks().length) highlight.set(0);
  });

  const pick = (row: Extract<Row, { kind: 'match' }>) => {
    if (row.match.blocked) return message.set(row.match.blocked);
    const r = manager.open(row.match.command.id, row.tail, row.match.preset?.id);
    if (r.kind === 'blocked') message.set(r.reason);
    else if (r.kind === 'ran') props.onClose();
  };

  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-expanded="true"
      aria-controls="msp-palette-list"
      aria-label="Command"
      placeholder="Search commands"
      spellcheck={false}
      value={text}
      onInput={(e) => {
        text.set(e.currentTarget.value);
        highlight.set(0);
        message.set(undefined);
      }}
      onKeyDown={(e) => {
        const n = picks().length;
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          if (n) highlight.set((highlight.value + (e.key === 'ArrowDown' ? 1 : n - 1)) % n);
        } else if (e.key === 'Enter') {
          e.preventDefault();
          const row = picks()[highlight.value];
          if (row) pick(row);
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
        <Icon name="search" size={16} class={classes.headIcon} />
        {input}
      </div>
      <Show when={message}>{(m) => <p class={classes.message}>{m}</p>}</Show>
      <div class={classes.list} id="msp-palette-list" role="listbox" aria-label="Commands">
        <For each={rows} fallback={() => <p class={classes.empty}>No commands match.</p>}>
          {(row) => {
            if (row.kind === 'caption') return <div class={classes.caption}>{row.label}</div>;
            const i = row.index;
            return (
              <CommandRow
                row={row}
                query={text.asReadonly}
                selected={highlight.asReadonly.pipe(map((h) => h === i))}
                onPick={() => pick(row)}
                onHover={() => highlight.set(i)}
              />
            );
          }}
        </For>
      </div>
      <Footer
        hints={[
          ['↑↓', 'move'],
          ['↵', 'run'],
          ['esc', 'close'],
        ]}
      />
    </>
  );
}

function CommandRow(props: {
  row: Extract<Row, { kind: 'match' }>;
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
      <span class={classes.rowIcon}>{c.icon ? <Icon name={c.icon} size={16} /> : null}</span>
      <span class={classes.rowTitle}>
        {props.query.pipe(map((q) => highlight(title, tail ? '' : q)))}
        {tail ? <span class={classes.tail}>{tail}</span> : null}
      </span>
      <span class={classes.rowMeta}>{match.blocked ?? c.category}</span>
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
      aria-controls="msp-palette-completions"
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
      <div class={classes.list} id="msp-palette-completions" role="listbox" aria-label="Options">
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
      </div>
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
