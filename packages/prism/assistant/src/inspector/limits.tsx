/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LimitRowView, LimitsGroupView } from '@molstar/ai-trace';
import { k } from '@molstar/ai-trace/format';
import { For, ScrollArea, Select, Show } from '@molstar/prism-ui';
import { distinctUntilChanged, filter, map, type Observable } from 'rxjs';
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
 * What bounds the agent (analysis/viewer-mode.md §5): each limit with what it bounds, its value (typed, committed on
 * Enter or leaving the field, Esc puts it back), its unit and its default, with a reset when it differs; the turn
 * budget's rows with what the latest turn used. Credits are shown, never a limit.
 */
export function LimitsView(props: { model: InspectorModel }) {
  const { model } = props;
  // Groups and rows followed by id: a change rewrites the records, and a field made anew would lose its focus.
  const groupIds = model.limits.pipe(
    map((l) => l.groups.map((g) => g.id)),
    distinctUntilChanged(sameKeys),
  );
  const group = (id: string) =>
    model.limits.pipe(
      map((l) => l.groups.find((g) => g.id === id)),
      filter((g): g is LimitsGroupView => !!g),
    );
  return (
    <ScrollArea class={classes.limits} contentClass={classes.content}>
      <For each={groupIds}>
        {(id) => {
          const g = group(id);
          const keys = g.pipe(
            map((x) => x.rows.map((r) => r.key)),
            distinctUntilChanged(sameKeys),
          );
          return (
            <section class={classes.group} aria-label={g.pipe(map((x) => x.title))}>
              <h3 class={classes.title}>{g.pipe(map((x) => x.title))}</h3>
              <Show when={g.pipe(map((x) => x.note ?? ''))}>{(note) => <p class={parts.note}>{note}</p>}</Show>
              <div class={classes.rows}>
                <For each={keys}>{(key) => <LimitRow model={model} rowKey={key} />}</For>
              </div>
            </section>
          );
        }}
      </For>
    </ScrollArea>
  );
}

function LimitRow(props: { model: InspectorModel; rowKey: string }) {
  const { model } = props;
  const row = model.limits.pipe(
    map((l) => findLimit(l, props.rowKey)),
    filter((r): r is LimitRowView => !!r),
    distinctUntilChanged(),
  );
  // The row's shape (whether it is a field, a choice or a fact; its label and unit) does not change; its numbers do.
  const shape = findLimit(model.limits.value, props.rowKey)!;
  const editable = model.editable && shape.editable;
  const field = <K extends keyof LimitRowView>(key: K) =>
    row.pipe(
      map((r) => r[key]),
      distinctUntilChanged(),
    );
  const changed = row.pipe(
    map((r) => r.value !== r.default),
    distinctUntilChanged(),
  );
  return (
    <div class={classes.row} data-changed={changed.pipe(map((c) => (c ? '' : undefined)))}>
      <div class={classes.label}>
        <span>{shape.label}</span>
        {shape.note ? <span class={classes.what}>{shape.note}</span> : null}
        {shape.used !== undefined ? (
          <span class={classes.used}>
            <Meter class={classes.usedBar} share={row.pipe(map((r) => (r.value > 0 ? (r.used ?? 0) / r.value : 0)))} />
            <span>
              {row.pipe(map((r) => `${k(r.used ?? 0)} of ${k(r.value)} used`))}
              {shape.ceiling !== undefined ? ` · up to ${k(shape.ceiling)}` : ''}
            </span>
          </span>
        ) : null}
      </div>
      {!editable ? (
        <span class={classes.value}>
          {row.pipe(map((r) => r.choices?.find((c) => c.value === r.value)?.label ?? show(r.value)))}
        </span>
      ) : shape.choices ? (
        <Select
          over
          size="sm"
          class={classes.choice}
          aria-label={shape.label}
          options={shape.choices.map((c) => [String(c.value), c.label] as const)}
          value={field('value').pipe(map(String))}
          onChange={(v) => model.setLimit(shape.key, Number(v))}
        />
      ) : (
        <LimitInput model={model} row={row} shape={shape} />
      )}
      <span class={classes.unit}>{shape.unit}</span>
      <span class={classes.default}>
        <Show when={changed} fallback={() => (editable || shape.value !== shape.default ? <span>default</span> : null)}>
          {() => (
            <>
              <span>
                default {row.pipe(map((r) => r.choices?.find((c) => c.value === r.default)?.label ?? show(r.default)))}
              </span>
              {editable ? (
                <button
                  type="button"
                  class={parts.link}
                  data-tip="Back to the default"
                  onClick={() => model.setLimit(shape.key, undefined)}
                >
                  Reset
                </button>
              ) : null}
            </>
          )}
        </Show>
      </span>
    </div>
  );
}

/** A limit's value as typed: committed on Enter or leaving the field; Esc puts it back (a second Esc closes the page). */
function LimitInput(props: { model: InspectorModel; row: Observable<LimitRowView>; shape: LimitRowView }) {
  const { model, row, shape } = props;
  let current = shape.value;
  const input = (
    <input
      class={classes.input}
      inputMode="decimal"
      aria-label={shape.label}
      value={row.pipe(
        map((r) => {
          current = r.value;
          return show(r.value);
        }),
      )}
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
  return input;
}
