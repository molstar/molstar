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
import { ConsoleTab } from './parts.js';
import { type ConsoleEntry, type ConsoleLineKind, ConsoleModel, ScriptExtensions } from './console-model.js';

// The console tab (_spec §9): the submissions and what each said, an input that takes several lines (Enter runs,
// ⇧Enter a new line, ↑/↓ the history at the input's first and last line, Tab the first completion), the dialect, and
// script files opened or dropped with the files they load. It acts through the person's API only.

const LineColors: Readonly<Record<ConsoleLineKind, string>> = {
  did: 'var(--msp-color-text)',
  value: 'var(--msp-color-text)',
  warning: 'var(--msp-color-warning, var(--msp-color-text-muted))',
  note: 'var(--msp-color-text-muted)',
  error: 'var(--msp-color-danger, var(--msp-color-text))',
};

const Mono = 'font-family: var(--msp-font-mono, ui-monospace, monospace); font-size: var(--msp-font-size-sm)';

function EntryView(props: { entry: ConsoleEntry }) {
  const { entry } = props;
  return (
    <div style="display: grid; gap: 2px; padding: var(--msp-space-1) 0; border-bottom: 1px solid var(--msp-color-border)">
      <div style={`${Mono}; display: flex; gap: var(--msp-space-2); white-space: pre-wrap; overflow-wrap: anywhere`}>
        <span style="color: var(--msp-color-text-muted)">{entry.dialect === 'pymol' ? 'PyMOL>' : 'vmd>'}</span>
        <span>{entry.input}</span>
      </div>
      {entry.state === 'running' ? (
        <div style={`${Mono}; color: var(--msp-color-text-muted)`}>…</div>
      ) : (
        entry.lines.map((l) => (
          <div style={`${Mono}; white-space: pre-wrap; overflow-wrap: anywhere; color: ${LineColors[l.kind]}`}>
            {l.kind === 'warning' || l.kind === 'note' ? `${l.kind}: ${l.text}` : l.text}
          </div>
        ))
      )}
    </div>
  );
}

export function ConsoleView(props: { model: ConsoleModel }) {
  const { model } = props;
  const list = (
    <div
      role="log"
      aria-label="Console output"
      style="flex: 1; min-height: 0; overflow: auto; padding: 0 var(--msp-space-2)"
    >
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
      rows={2}
      spellcheck={false}
      aria-label="Command"
      style={`${Mono}; width: 100%; resize: vertical; box-sizing: border-box; padding: var(--msp-space-1) var(--msp-space-2); background: var(--msp-color-surface, transparent); color: var(--msp-color-text); border: 1px solid var(--msp-color-border); border-radius: var(--msp-radius-sm, 4px)`}
      onInput={(e) => model.setDraft((e.currentTarget as HTMLTextAreaElement).value)}
      onKeyDown={(e) => {
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
      type="file"
      multiple
      style="display: none"
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
      style="display: flex; flex-direction: column; height: 100%; min-height: 0; gap: var(--msp-space-2)"
      onDragOver={(e) => e.preventDefault()}
      onDrop={(e) => {
        const files = [...((e as DragEvent).dataTransfer?.files ?? [])];
        if (!files.length) return;
        e.preventDefault();
        e.stopPropagation();
        void model.openFiles(files);
      }}
    >
      <div style="display: flex; align-items: center; gap: var(--msp-space-1); padding: var(--msp-space-1) var(--msp-space-2) 0">
        <div role="group" aria-label="Dialect" style="display: flex; gap: 2px">
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
        <span style="flex: 1" />
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
      <div style="padding: 0 var(--msp-space-2) var(--msp-space-2)">{input}</div>
    </div>
  );
}

function EmptyHint() {
  return (
    <div style="color: var(--msp-color-text-muted); padding: var(--msp-space-2) 0; font-size: var(--msp-font-size-sm)">
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
