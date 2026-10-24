/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, effect, For, Icon, Show } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { applyCompletion, complete } from '../language/completion.js';
import { type Dialect, DialectNames, Dialects } from '../language/ir.js';
import { dialectPrefix } from '../language/script.js';
import { classes } from './_generated/console.module.css.js';
import {
  ConsoleExamples,
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

function EntryView(props: { entry: ConsoleEntry; model: ConsoleModel }) {
  const { entry, model } = props;
  const label = entry.state === 'running' ? 'Running' : entry.state === 'failed' ? 'Failed' : 'Done';
  return (
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
  );
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

/** The composer's line: what Enter will do, or why it will not. */
function StatusView(props: { inspection: ConsoleInspection }) {
  const i = props.inspection;
  if (i.kind === 'empty') return <span class={classes.does} />;
  if (i.kind === 'error') {
    return (
      <span class={classes.does} data-kind="error">
        {i.message}
        {i.hint ? ` ${i.hint}` : ''}
      </span>
    );
  }
  const warning = i.warnings.length ? ` · ${i.warnings[0]}` : '';
  return (
    <span class={classes.does} data-tip={i.does + warning}>
      {DialectNames[i.dialect]} · {i.does}
      {warning}
    </span>
  );
}

export function ConsoleView(props: { model: ConsoleModel }) {
  const { model } = props;
  const dropping = observableState(false);

  // A log, not live: the viewer says each run once, on its own region (accessibility.md).
  const list = (
    <div class={classes.output} role="log" aria-live="off" aria-label="Console output">
      <For each={model.entries} fallback={() => <EmptyView model={model} onPick={(t) => fill(t)} />}>
        {(e) => <EntryView entry={e} model={model} />}
      </For>
    </div>
  ) as HTMLElement;
  // The newest at the bottom, in view, unless the person scrolled up to read (then the transcript stays).
  let following = true;
  list.addEventListener('scroll', () => {
    following = list.scrollHeight - list.scrollTop - list.clientHeight < 24;
  });
  effect(model.entries, () =>
    queueMicrotask(() => {
      if (following) list.scrollTop = list.scrollHeight;
    }),
  );

  const input = (
    <textarea
      class={classes.input}
      rows={2}
      spellcheck={false}
      aria-label="Command"
      onInput={(e) => model.setDraft((e.currentTarget as HTMLTextAreaElement).value)}
      onKeyDown={(e) => {
        // Composing text (an input method) owns its keys.
        if (e.isComposing) return;
        const el = e.currentTarget as HTMLTextAreaElement;
        if (e.key === 'Enter' && !e.shiftKey) {
          e.preventDefault();
          run();
        } else if (e.key === 'ArrowUp' && !el.value.slice(0, el.selectionStart).includes('\n')) {
          e.preventDefault();
          model.previous();
        } else if (e.key === 'ArrowDown' && !el.value.slice(el.selectionEnd).includes('\n')) {
          e.preventDefault();
          model.next();
        } else if (e.key === 'Tab' && !e.shiftKey && el.selectionStart === el.value.length) {
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

  // What the input will do, and its completions: read again as the draft, the dialect or the scene change (a run that
  // ends, a structure loaded: the entries change with it).
  const inspection = derive([model.draft, model.dialect, model.entries], (t) => model.inspect(t));
  const suggestions = derive([model.draft, model.dialect, model.entries], (t) =>
    t.trim() ? model.completions(t) : [],
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
            <button
              type="button"
              class={classes.segment}
              aria-pressed={derive(model.dialect, (x) => String(x === d))}
              onClick={() => model.setDialect(d)}
            >
              {DialectNames[d]}
            </button>
          ))}
        </div>
        <span class={classes.spacer} />
        <Button size="sm" variant="ghost" icon={<Icon name="file-code" size={15} />} onClick={() => picker.click()}>
          Open script…
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
        <div class={classes.suggestions}>
          <For each={suggestions} by={(s) => s.text}>
            {(s) => (
              <button type="button" class={classes.suggestion} data-tip={s.detail} onClick={() => fill(s.text)}>
                {s.label}
              </button>
            )}
          </For>
        </div>
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
          <Show when={inspection}>{(i) => <StatusView inspection={i} />}</Show>
          <span class={classes.keys} data-tip="Enter runs · Shift+Enter a new line · Tab completes · ↑ ↓ history">
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
        title: 'Console',
        icon: 'file-code',
        place: 'foot',
        order: 2.4,
        render: domView(() => <ConsoleView model={model} />),
      },
    ],
  };
}
