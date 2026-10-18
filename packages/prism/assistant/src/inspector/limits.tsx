/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LimitRowView } from '@molstar/ai/trace';
import { k } from '@molstar/ai/trace';
import { For, type ReadonlyState, ScrollArea, Select } from '@molstar/prism-ui';
import { watch } from '@molstar/prism-ui/dom/scope';
import { setAttr, setShown, setText } from '@molstar/prism-ui/dom/writes';
import { derive } from '@molstar/prism-ui/reactive/state';
import { held, present } from '../held.js';
import { classes } from './_generated/limits.module.css.js';
import { findLimit, type InspectorModel } from './model.js';
import { Meter, parts } from './parts.js';

const show = (n: number) => n.toLocaleString('en-US', { maximumFractionDigits: 2 });

const sameKeys = (a: readonly string[], b: readonly string[]) => a.length === b.length && a.every((x, i) => x === b[i]);

/** Reads what was typed: `500000`, `500,000`, `500k`, `1.5M`; undefined when it is not a number. */
export function parseLimit(text: string): number | undefined {
  const m = /^([\d.]+)\s*([kKmM]?)$/.exec(text.replace(/[,\s_]/g, ''));
  if (!m) return undefined;
  const n = Number(m[1]) * (m[2] === 'k' || m[2] === 'K' ? 1e3 : m[2] ? 1e6 : 1);
  return Number.isFinite(n) ? n : undefined;
}

/**
 * What bounds the agent (packages/ai/_spec/trace.md § The Inspector): each limit with what it bounds, its value
 * (typed, committed on Enter or leaving the field, Esc puts it back), its unit and its default, with a reset when it
 * differs; the turn budget's rows with what the latest turn used. Credits are shown, never a limit.
 */
export function LimitsView(props: { model: InspectorModel }) {
  const { model } = props;
  // Groups and rows followed by id: a change rewrites the records, and a field made anew would lose its focus.
  const groupIds = derive(model.limits, (l) => l.groups.map((g) => g.id), sameKeys);
  return (
    <ScrollArea class={classes.limits} contentClass={classes.content}>
      <For each={groupIds}>{(id) => <LimitGroup model={model} id={id} />}</For>
    </ScrollArea>
  );
}

/** A group: its title and note follow its record, its rows are followed by key. */
function LimitGroup(props: { model: InspectorModel; id: string }) {
  const { model, id } = props;
  const group = held(
    derive(model.limits, (l) => l.groups.find((g) => g.id === id)),
    present,
  );
  const keys = derive(group, (g) => g.rows.map((r) => r.key), sameKeys);
  const title = document.createTextNode('');
  const noteText = document.createTextNode('');
  const note = (<p class={parts.note}>{noteText}</p>) as HTMLElement;
  const root = (
    <section class={classes.group}>
      <h3 class={classes.title}>{title}</h3>
      {note}
      <div class={classes.rows}>
        <For each={keys}>{(key) => <LimitRow model={model} rowKey={key} />}</For>
      </div>
    </section>
  ) as HTMLElement;
  watch([group], () => {
    const g = group.value;
    setAttr(root, 'aria-label', g.title);
    setText(title, g.title);
    setText(noteText, g.note ?? '');
    setShown(note, !!g.note);
  });
  return root;
}

function LimitRow(props: { model: InspectorModel; rowKey: string }) {
  const { model } = props;
  const row = held(
    derive(model.limits, (l) => findLimit(l, props.rowKey)),
    present,
  );
  // The row's shape (whether it is a field, a choice or a fact; its label and unit) does not change; its numbers do.
  const shape = findLimit(model.limits.value, props.rowKey)!;
  const editable = model.editable && shape.editable;
  const labelOf = (r: LimitRowView, value: number) => r.choices?.find((c) => c.value === value)?.label ?? show(value);
  const usedBar =
    shape.used !== undefined ? derive(row, (r) => (r.value > 0 ? (r.used ?? 0) / r.value : 0)) : undefined;
  const usedText = document.createTextNode('');
  const valueText = !editable ? document.createTextNode('') : undefined;
  const defaultText = document.createTextNode('');
  const reset = editable
    ? ((
        <button
          type="button"
          class={parts.link}
          data-tip="Back to the default"
          onClick={() => model.setLimit(shape.key, undefined)}
        >
          Reset
        </button>
      ) as HTMLElement)
    : undefined;
  const root = (
    <div class={classes.row}>
      <div class={classes.label}>
        <span>{shape.label}</span>
        {shape.note ? <span class={classes.what}>{shape.note}</span> : null}
        {usedBar ? (
          <span class={classes.used}>
            <Meter class={classes.usedBar} share={usedBar} />
            <span>
              {usedText}
              {shape.ceiling !== undefined ? ` · up to ${k(shape.ceiling)}` : ''}
            </span>
          </span>
        ) : null}
      </div>
      {valueText ? (
        <span class={classes.value}>{valueText}</span>
      ) : shape.choices ? (
        <Select
          over
          size="sm"
          class={classes.choice}
          aria-label={shape.label}
          options={shape.choices.map((c) => [String(c.value), c.label] as const)}
          value={derive(row, (r) => String(r.value))}
          onChange={(v) => model.setLimit(shape.key, Number(v))}
        />
      ) : (
        <LimitInput model={model} row={row} shape={shape} />
      )}
      <span class={classes.unit}>{shape.unit}</span>
      <span class={classes.default}>
        {defaultText}
        {reset}
      </span>
    </div>
  ) as HTMLElement;
  // One render for the numbers of its record: what changed against what the DOM has is written.
  watch([row], () => {
    const r = row.value;
    const changed = r.value !== r.default;
    setAttr(root, 'data-changed', changed);
    setText(usedText, `${k(r.used ?? 0)} of ${k(r.value)} used`);
    if (valueText) setText(valueText, labelOf(r, r.value));
    // Changed: the default and a way back; else the word, for a field or a value that was set away from it.
    setText(
      defaultText,
      changed ? `default ${labelOf(r, r.default)}` : editable || shape.value !== shape.default ? 'default' : '',
    );
    if (reset) setShown(reset, changed);
  });
  return root;
}

/** A limit's value as typed: committed on Enter or leaving the field; Esc puts it back (a second Esc closes the page). */
function LimitInput(props: { model: InspectorModel; row: ReadonlyState<LimitRowView>; shape: LimitRowView }) {
  const { model, row, shape } = props;
  let current = shape.value;
  const input = (
    <input
      class={classes.input}
      inputMode="decimal"
      aria-label={shape.label}
      onKeyDown={(e) => {
        if (e.key === 'Enter') {
          commit();
          e.preventDefault();
        } else if (e.key === 'Escape' && input.value !== show(current)) {
          input.value = show(current);
          e.preventDefault();
          e.stopPropagation();
        }
      }}
      onBlur={() => commit()}
    />
  ) as HTMLInputElement;
  const commit = () => {
    const n = parseLimit(input.value);
    if (n !== undefined && n !== current) model.setLimit(shape.key, n);
    // Refused, or the same: the field shows the value as it is.
    input.value = show(current);
  };
  // The value it holds follows the limit's; a field being typed in is written only when the limit's value changes.
  let written: number | undefined;
  watch([row], () => {
    const v = row.value.value;
    current = v;
    if (v === written) return;
    written = v;
    if (input.value !== show(v)) input.value = show(v);
  });
  return input;
}
