/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, effect, For, Icon } from '@molstar/prism-ui';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { applyCompletion, complete } from '../language/completion.js';
import { type Dialect, DialectNames, Dialects } from '../language/ir.js';
import { dialectPrefix } from '../language/script.js';
import { classes } from './_generated/console.module.css.js';
import { ConsoleTab } from './parts.js';
import { type ConsoleEntry, type ConsoleLineKind, ConsoleModel, ScriptExtensions } from './console-model.js';

// The console tab (_spec §9): the submissions and what each said, an input that takes several lines (Enter runs,
// ⇧Enter a new line, ↑/↓ the history at the input's first and last line, Tab the first completion), the dialect, and
// script files opened or dropped with the files they load. It acts through the person's API only.

const LineClass: Readonly<Record<ConsoleLineKind, string>> = {
  did: classes.line,
  value: classes.line,
  warning: classes.warning,
  note: classes.quiet,
  error: classes.error,
};

function EntryView(props: { entry: ConsoleEntry }) {
  const { entry } = props;
  return (
    <div class={classes.entry}>
      <div class={`${classes.mono} ${classes.input}`}>
        <span class={classes.prompt}>{entry.dialect === 'pymol' ? 'PyMOL>' : 'vmd>'}</span>
        <span>{entry.input}</span>
      </div>
      {entry.state === 'running' ? (
        <div class={`${classes.mono} ${classes.quiet}`}>…</div>
      ) : (
        entry.lines.map((l) => (
          <div class={`${classes.mono} ${LineClass[l.kind]}`}>
            {l.kind === 'warning' || l.kind === 'note' ? `${l.kind}: ${l.text}` : l.text}
          </div>
        ))
      )}
    </div>
  );
}

export function ConsoleView(props: { model: ConsoleModel }) {
  const { model } = props;
  // A log, not live: the viewer says each run once, on its own region (accessibility.md).
  const list = (
    <div class={classes.output} role="log" aria-live="off" aria-label="Console output">
      <For each={model.entries} fallback={() => <EmptyHint />}>
        {(e) => <EntryView entry={e} />}
      </For>
    </div>
  ) as HTMLElement;
  // The newest at the bottom, in view.
  effect(model.entries, () =>
    queueMicrotask(() => {
      list.scrollTop = list.scrollHeight;
    }),
  );

  const input = (
    <textarea
      class={`${classes.field} ${classes.mono}`}
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
          void model.submit(el.value);
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
  effect(model.draft, (text) => {
    if (input.value !== text) {
      input.value = text;
      input.setSelectionRange(text.length, text.length);
    }
  });
  effect(model.dialect, (d) => {
    input.placeholder = model.placeholder(d);
  });

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

  return (
    <div
      class={classes.console}
      onDragOver={(e) => e.preventDefault()}
      onDrop={(e) => {
        const files = [...((e as DragEvent).dataTransfer?.files ?? [])];
        if (!files.length) return;
        e.preventDefault();
        e.stopPropagation();
        void model.openFiles(files);
      }}
    >
      <div class={classes.toolbar}>
        <div class={classes.dialects} role="group" aria-label="Dialect">
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
      <div class={classes.foot}>{input}</div>
    </div>
  );
}

function EmptyHint() {
  return (
    <div class={classes.hint}>
      PyMOL or VMD commands: one a line, Enter runs them as one step (Undo takes it back). Type help for the commands;
      open or drop a .pml, .tcl or .vmd script with the files it loads.
    </div>
  );
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
