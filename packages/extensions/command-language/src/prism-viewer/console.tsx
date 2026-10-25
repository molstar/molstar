/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, type Child, domView, effect, For, Icon, Show } from '@molstar/prism-ui';
import { chordFromEvent, chordId, formatChord } from '@molstar/prism-ui/commands/keys';
import { frames, observeResize } from '@molstar/prism-ui/dom/scoped';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { applyCompletion, complete, type Signature } from '../language/completion.js';
import { type Dialect, DialectNames, Dialects } from '../language/ir.js';
import { dialectPrefix } from '../language/script.js';
import { classes } from './_generated/console.module.css.js';
import {
  ConsoleExamples,
  type ConsoleCompletion,
  type ConsoleEntry,
  type ConsoleInspection,
  type ConsoleLine,
  ConsoleModel,
  ScriptExtensions,
} from './console-model.js';
import { ConsoleTab } from './parts.js';

// The console tab (_spec §9): a transcript of runs (each its status mark, its dialect, the command and what it said),
// and a composer whose line says what the input will do before Enter, with the completions of its last word. Enter
// runs, ⇧Enter a new line, ↑/↓ the history at the input's first and last line, Tab the first completion. Script files
// opened or dropped with the files they load. It acts through the person's API only.

const Marks: Readonly<Record<ConsoleEntry['state'], 'loader-circle' | 'check' | 'circle-alert'>> = {
  running: 'loader-circle',
  done: 'check',
  failed: 'circle-alert',
};

const timeOf = (ms: number | undefined) =>
  ms === undefined ? '' : ms < 1000 ? `${ms} ms` : `${(ms / 1000).toFixed(ms < 10000 ? 1 : 0)} s`;

function LineView(props: { line: ConsoleLine }) {
  const { line } = props;
  if (line.kind === 'value') return <div class={classes.value}>{line.text}</div>;
  if (line.kind === 'did') return <div class={classes.line}>{line.text}</div>;
  const icon = line.kind === 'error' ? 'circle-alert' : line.kind === 'warning' ? 'triangle-alert' : 'info';
  return (
    <div class={classes.notice} data-kind={line.kind}>
      <span class={classes.noticeIcon}>
        <Icon name={icon} size={13} />
      </span>
      <span>{line.text}</span>
      {line.hint ? <span class={classes.hint}>{line.hint}</span> : null}
    </div>
  );
}

function EntryView(props: { entry: ConsoleEntry; model: ConsoleModel; onSize: () => void }) {
  const { entry, model } = props;
  const label = entry.state === 'running' ? 'Running' : entry.state === 'failed' ? 'Failed' : 'Done';
  const el = (
    <div class={classes.entry} data-state={entry.state}>
      <span class={classes.mark} role="img" aria-label={label}>
        <Icon name={Marks[entry.state]} size={14} />
      </span>
      <div class={classes.command}>
        <span class={classes.chip}>{DialectNames[entry.dialect]}</span>
        <span class={classes.code}>{entry.input}</span>
      </div>
      <div class={classes.side}>
        <span class={classes.time}>{timeOf(entry.ms)}</span>
        {model.editable(entry) ? (
          <button
            type="button"
            class={classes.action}
            aria-label="Edit and run again"
            data-tip="Edit and run again"
            onClick={() => model.edit(entry)}
          >
            <Icon name="pencil" size={13} />
          </button>
        ) : null}
      </div>
      {entry.lines.length ? (
        <div class={classes.lines}>
          {entry.lines.map((l) => (
            <LineView line={l} />
          ))}
        </div>
      ) : null}
    </div>
  ) as HTMLElement;
  // Its size, from its first layout and as what it said comes in: the transcript follows (ConsoleView).
  observeResize(el, props.onSize);
  return el;
}

function EmptyView(props: { model: ConsoleModel; onPick: (text: string) => void }) {
  return (
    <div class={classes.empty}>
      <div class={classes.emptyTitle}>PyMOL and VMD commands</div>
      <div class={classes.emptyText}>
        Type the commands you know. Enter runs them as one step, which Undo takes back; the line under the input says
        what they will do first. Open or drop a .pml, .tcl or .vmd script with the files it loads.
      </div>
      <div class={classes.examples}>
        <For each={derive(props.model.dialect, (d) => ConsoleExamples[d])}>
          {(example) => (
            <button type="button" class={classes.example} onClick={() => props.onPick(example)}>
              {example}
            </button>
          )}
        </For>
      </div>
    </div>
  );
}

/** The composer's line, one row whatever it says: what Enter will do, or why it will not (the whole in its tip). */
function StatusView(props: { inspection: ConsoleInspection; usageShown?: boolean }) {
  const i = props.inspection;
  if (i.kind === 'empty') return <span class={classes.does} />;
  if (i.kind === 'error') {
    const full = i.hint ? `${i.message} ${i.hint}` : i.message;
    return (
      <span class={classes.does} data-kind="error" data-tip={full}>
        {i.message}
        {/* The hint is mostly the usage: left to the card when it shows it (still in the tip). */}
        {i.hint && !props.usageShown ? <span class={classes.hint}> {i.hint}</span> : null}
      </span>
    );
  }
  const said = [i.does, ...i.warnings].filter(Boolean).join(' · ');
  return (
    <span class={classes.does} data-tip={said}>
      {DialectNames[i.dialect]} · {said}
    </span>
  );
}

/** Signature help, the card's header: the command's usage with the argument being typed in bold, then what it does. */
function SignatureView(props: { signature: Signature; dialect: Dialect }) {
  const { signature: sig, dialect } = props;
  const sep = dialect === 'pymol' ? ', ' : ' ';
  const parts: Child[] = [<span class={classes.sigName}>{sig.name}</span>, ' '];
  sig.params.forEach((p, k) => {
    if (k) parts.push(sep);
    const text = p.optional ? `[${p.name}]` : p.name;
    parts.push(k === sig.active ? <strong class={classes.sigActive}>{text}</strong> : text);
  });
  return (
    <div class={classes.signature}>
      <div class={classes.usage}>{parts}</div>
      <div class={classes.summary}>{sig.summary}</div>
    </div>
  );
}

/** What a completion is, in a few letters (a color shows its swatch instead). */
const KindTags: Partial<Record<ConsoleCompletion['kind'], string>> = {
  command: 'cmd',
  keyword: 'sel',
  selection: 'sel',
  rep: 'rep',
  color: 'color',
  setting: 'set',
  value: 'val',
  object: 'obj',
  name: 'name',
  molecule: 'mol',
  file: 'file',
  enum: 'opt',
  bool: 'opt',
};

/** A completion's label, the part typed so far marked. */
function labelWithTyped(label: string, typed: string): Child {
  const at = typed ? label.toLowerCase().indexOf(typed.toLowerCase()) : -1;
  if (at < 0) return label;
  return [
    label.slice(0, at),
    <mark class={classes.typed}>{label.slice(at, at + typed.length)}</mark>,
    label.slice(at + typed.length),
  ];
}

/** What asks for the completion list (⌘I, Ctrl+I). */
const ListChord = { key: 'i', mod: true } as const;

/** Each console's list its own id (several viewers on a page). */
let nextList = 0;

export function ConsoleView(props: { model: ConsoleModel }) {
  const { model } = props;
  const dropping = observableState(false);

  // A log, not live: the viewer says each run once, on its own region (accessibility.md).
  const list = (
    <div class={classes.output} role="log" aria-live="off" aria-label="Console output">
      <For each={model.entries} fallback={() => <EmptyView model={model} onPick={(t) => fill(t)} />}>
        {(e) => <EntryView entry={e} model={model} onSize={() => stick()} />}
      </For>
    </div>
  ) as HTMLElement;
  // The newest at the bottom, always in view: a run appended, what it said once it ends, the box resized. Each entry
  // says when its size changes (it renders its outcome on its own); scrolled once a frame, however many say so.
  const stick = frames(() => {
    list.scrollTop = list.scrollHeight;
  }).request;
  observeResize(list, stick);

  const input = (
    <textarea
      class={classes.input}
      rows={2}
      spellcheck={false}
      aria-label="Command"
      // A textarea takes no combobox role (ARIA in HTML): a textbox that says it completes from a list, which one, and
      // the row picked in it.
      aria-autocomplete="list"
      onInput={(e) => model.setDraft((e.currentTarget as HTMLTextAreaElement).value)}
      onKeyDown={(e) => {
        // Composing text (an input method) owns its keys; Safari confirms with an Enter that says so by its code only.
        if (e.isComposing || e.keyCode === 229) return;
        const el = e.currentTarget as HTMLTextAreaElement;
        // ⌘I (Ctrl+I) asks for the list: every command on the empty input, the list again after Esc.
        const chord = chordFromEvent(e);
        if (chord && chordId(chord) === chordId(ListChord)) {
          e.preventDefault();
          dismissed.set(false);
          asked.set(true);
          return;
        }
        // The completion list, while open and the caret at the end (what it completes): ↑ ↓ move in it, Tab or Enter
        // take the one picked (Tab: the first when none is), Esc closes it; Enter with none picked runs, as always.
        const items = shown.value;
        if (card.value && e.key === 'Escape') {
          e.preventDefault();
          e.stopPropagation();
          dismissed.set(true);
          return;
        }
        if (items.length && el.selectionStart === el.value.length) {
          const at = active.value;
          if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
            e.preventDefault();
            const n = items.length;
            active.set(e.key === 'ArrowDown' ? (at + 1) % n : at <= 0 ? n - 1 : at - 1);
            return;
          }
          if (e.key === 'Tab' && !e.shiftKey) {
            e.preventDefault();
            accept(items[Math.max(0, at)]);
            return;
          }
          if (e.key === 'Enter' && !e.shiftKey && at >= 0) {
            e.preventDefault();
            accept(items[at]);
            return;
          }
        }
        if (e.key === 'Enter' && !e.shiftKey) {
          e.preventDefault();
          run();
        } else if (e.key === 'ArrowUp' && !el.value.slice(0, el.selectionStart).includes('\n')) {
          e.preventDefault();
          recall(() => model.previous());
        } else if (e.key === 'ArrowDown' && !el.value.slice(el.selectionEnd).includes('\n')) {
          e.preventDefault();
          recall(() => model.next());
        } else if (e.key === 'Tab' && !e.shiftKey && el.value.trim() && el.selectionStart === el.value.length) {
          const filled = completeLine(model, el.value);
          if (filled !== undefined) {
            e.preventDefault();
            model.setDraft(filled);
          }
        }
      }}
    />
  ) as HTMLTextAreaElement;
  const fill = (text: string) => {
    model.setDraft(text);
    input.focus();
  };
  const accept = (item: ConsoleCompletion) => fill(item.text);
  // A line from the history comes without its list, so ↑ ↓ go on through the history; typing brings the list back.
  const recall = (step: () => void) => {
    step();
    dismissed.set(true);
  };
  const run = () => {
    if (input.value.trim()) void model.submit(input.value);
  };
  effect(model.draft, (text) => {
    if (input.value !== text) {
      input.value = text;
      input.setSelectionRange(text.length, text.length);
    }
  });
  effect(model.dialect, (d) => {
    input.placeholder = model.placeholder(d);
  });

  // What helps the input, read again as the draft, the dialect or the scene change (a run that ends, a structure
  // loaded: the entries change with it): the completions of the word at its end, the command it is in, and what Enter
  // will do.
  // ⌘I on the empty input asks for the list: every command of the dialect, until the next keystroke. (Not Tab: there it
  // moves the focus on, as anywhere.)
  const asked = observableState(false);
  const help = derive([model.draft, model.dialect, model.entries, asked.asReadonly], (t, _d, _e, a) => ({
    assist: a && !t.trim() ? model.assist(t, Number.POSITIVE_INFINITY, true) : model.assist(t),
    inspection: model.inspect(t),
    dialect: dialectPrefix(t)?.dialect ?? model.dialect.value,
  }));
  // Esc closes the card (the command's help and the list) until the next keystroke; after every keystroke none is
  // picked (-1).
  const dismissed = observableState(false);
  const active = observableState(-1);
  effect(model.draft, () => {
    asked.set(false);
    dismissed.set(false);
  });
  // The card over the composer: the command's help from its name typed until the line runs (so it does not come and go
  // between arguments), and the completions of the word at the end when there are any.
  const card = derive([help, dismissed.asReadonly], (h, d) =>
    d || (!h.assist.signature && !h.assist.items.length)
      ? undefined
      : { items: h.assist.items, signature: h.assist.signature, dialect: h.dialect, typed: h.assist.typed },
  );
  const shown = derive(card, (c) => c?.items ?? []);
  // None picked whenever the rows change: a keystroke, the dialect, a run that ends (the index would name another row,
  // or none).
  effect(shown, () => active.set(-1));
  const listId = `msp-console-completions-${(nextList += 1)}`;
  effect(
    derive([shown, active.asReadonly], (l, a) => (l.length ? a : undefined)),
    (at) => {
      if (at !== undefined) input.setAttribute('aria-controls', listId);
      else input.removeAttribute('aria-controls');
      if (at !== undefined && at >= 0) {
        input.setAttribute('aria-activedescendant', `${listId}-${at}`);
        root.querySelector(`#${listId}-${at}`)?.scrollIntoView?.({ block: 'nearest' });
      } else input.removeAttribute('aria-activedescendant');
    },
  );
  const empty = derive(model.draft, (t) => !t.trim());

  const picker = (
    <input
      class={classes.picker}
      type="file"
      multiple
      onChange={(e) => {
        const el = e.currentTarget as HTMLInputElement;
        const files = [...(el.files ?? [])];
        el.value = '';
        if (files.length) void model.openFiles(files);
      }}
    />
  ) as HTMLInputElement;
  picker.accept = [...ScriptExtensions.map((x) => `.${x}`), '*/*'].join(',');

  const root = (
    <div
      class={classes.console}
      onDragEnter={(e) => {
        if ((e as DragEvent).dataTransfer?.types.includes('Files')) dropping.set(true);
      }}
      onDragOver={(e) => e.preventDefault()}
      onDragLeave={(e) => {
        if (!root.contains((e as DragEvent).relatedTarget as Node | null)) dropping.set(false);
      }}
      onDrop={(e) => {
        dropping.set(false);
        const files = [...((e as DragEvent).dataTransfer?.files ?? [])];
        if (!files.length) return;
        e.preventDefault();
        e.stopPropagation();
        void model.openFiles(files);
      }}
    >
      <div class={classes.bar}>
        <div class={classes.segments} role="group" aria-label="Dialect">
          {Dialects.map((d: Dialect) => (
            <Button
              size="sm"
              variant="ghost"
              pressed={derive(model.dialect, (x) => x === d)}
              onClick={() => model.setDialect(d)}
            >
              {DialectNames[d]}
            </Button>
          ))}
        </div>
        <span class={classes.spacer} />
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="upload" size={15} />}
          title="Import a .pml, .tcl or .vmd script with the files it loads, or structure files (or drop them here)"
          onClick={() => picker.click()}
        >
          Import
        </Button>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="trash-2" size={15} />}
          title="Clear the output"
          onClick={() => model.clear()}
        />
        {picker}
      </div>
      {list}
      <div class={classes.composer}>
        <Show when={card}>
          {(c) => {
            const { items } = c;
            // A detail every row shares says nothing the signature does not (`representation` on each): left out.
            const details = items.some((i) => i.detail !== items[0].detail);
            return (
              <div class={classes.card}>
                {c.signature ? <SignatureView signature={c.signature} dialect={c.dialect} /> : null}
                {items.length ? (
                  <div class={classes.completions} id={listId} role="listbox" aria-label="Completions">
                    {items.map((item, k) => (
                      <div
                        id={`${listId}-${k}`}
                        class={classes.completion}
                        role="option"
                        aria-selected={derive(active.asReadonly, (a) => String(a === k))}
                        onPointerMove={() => active.set(k)}
                        onMouseDown={(e) => {
                          // Keeps the input's focus.
                          e.preventDefault();
                          accept(item);
                        }}
                      >
                        <span class={classes.kind}>
                          {item.rgb !== undefined ? (
                            <span
                              class={classes.swatch}
                              style={`background: #${item.rgb.toString(16).padStart(6, '0')}`}
                            />
                          ) : (
                            (KindTags[item.kind] ?? '')
                          )}
                        </span>
                        <span class={classes.label}>{labelWithTyped(item.label, c.typed)}</span>
                        <span class={classes.detail}>{details ? (item.detail ?? '') : ''}</span>
                      </div>
                    ))}
                  </div>
                ) : null}
              </div>
            );
          }}
        </Show>
        <div class={classes.field}>
          {input}
          <button
            type="button"
            class={classes.send}
            aria-label="Run"
            data-tip="Run (Enter)"
            disabled={empty}
            onClick={run}
          >
            <Icon name="arrow-up" size={14} />
          </button>
        </div>
        <div class={classes.status}>
          <Show when={derive([help, card], (h, c) => ({ inspection: h.inspection, usageShown: !!c?.signature }))}>
            {(st) => <StatusView inspection={st.inspection} usageShown={st.usageShown} />}
          </Show>
          <span
            class={classes.keys}
            data-tip={`Enter runs · Shift+Enter a new line · Tab completes · ${formatChord(ListChord)} lists the completions · ↑ ↓ the list or the history · Esc closes the help and the list`}
          >
            ↵ run · ⇥ complete
          </span>
        </div>
      </div>
      <div class={classes.drop}>Drop a .pml, .tcl or .vmd script with the files it loads</div>
    </div>
  ) as HTMLElement;
  effect(dropping.asReadonly, (on) => {
    if (on) root.setAttribute('data-dropping', '');
    else root.removeAttribute('data-dropping');
  });
  return root;
}

/** The line with its first completion applied, or undefined when there is none. */
export function completeLine(model: ConsoleModel, text: string): string | undefined {
  const area = model.api.script;
  if (!area) return undefined;
  const prefixed = dialectPrefix(text);
  const dialect = prefixed?.dialect ?? model.dialect.value;
  const body = prefixed ? prefixed.rest : text;
  const c = complete(body, dialect, area.env());
  const first = c.items.find((i) => i.insert !== body.slice(c.from, c.to));
  if (!first) return undefined;
  return `${prefixed ? text.slice(0, prefixed.offset) : ''}${applyCompletion(body, c, first)}`;
}

export function consoleUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const model = new ConsoleModel(ctx.api, ctx.preferences);
  return {
    tabs: [
      {
        id: ConsoleTab,
        // The rail's name says what it is; the header has the dialect switch under it.
        title: 'PyMOL / VMD console',
        heading: 'Console',
        icon: 'terminal',
        place: 'foot',
        // Beside the Assistant (2.65): the other panel one drives the viewer by typing in.
        order: 2.6,
        render: domView(() => <ConsoleView model={model} />),
      },
    ],
  };
}
