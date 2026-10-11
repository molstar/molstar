/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { ColorLists, type ColorListName } from '@molstar/core/util/color/lists';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import { colorToHex } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/extra-controls.module.css.js';
import { classes as form } from './_generated/params-form.module.css.js';
import type { FieldBinding } from './field-binding.js';

/** A row with a field that opens something under it (a list, a grid, an editor), like a dropdown. */
function ExpandingRow(props: {
  head: Child;
  field: (open: ObservableState<boolean>) => Child;
  body: (close: () => void) => Child;
}) {
  const open = observableState(false);
  return (
    <>
      <div class={form.row}>
        <span class={form.label}>{props.head}</span>
        <div class={form.control}>{props.field(open)}</div>
      </div>
      <Show when={open} reveal>
        {() => <div class={classes.under}>{props.body(() => open.set(false))}</div>}
      </Show>
    </>
  );
}

/** A field like a dropdown's: its content, a chevron; toggles `open`. */
function Field(props: { open: ObservableState<boolean>; label: string; children: Child }) {
  return (
    <button
      type="button"
      class={classes.field}
      aria-expanded={props.open}
      aria-label={props.label}
      onClick={() => props.open.set(!props.open.value)}
    >
      <span class={classes.fieldContent}>{props.children}</span>
      <Icon name="chevron-down" size={14} class={classes.chevron} />
    </button>
  );
}

/**
 * Multi-select: the chosen options' names in a dropdown-like field ("None" when empty); opened, the options under the
 * row, each toggling (a check on the chosen), the list staying open.
 */
export function MultiSelectRow(props: {
  head: Child;
  label: string;
  param: PD.MultiSelect<any>;
  binding: FieldBinding;
}) {
  const { param, binding } = props;
  const names = binding.value.pipe(
    map((v: readonly unknown[]) => {
      const chosen = param.options.filter(([o]) => v?.includes(o)).map(([, l]) => l);
      return chosen.length ? chosen.join(', ') : (param.emptyValue ?? 'None');
    }),
  );
  return (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          {names}
        </Field>
      )}
      body={() => (
        <ScrollArea
          class={classes.options}
          contentClass={classes.optionsContent}
          role="listbox"
          aria-multiselectable="true"
          aria-label={props.label}
        >
          {param.options.map(([value, label]) => (
            <button
              type="button"
              role="option"
              class={classes.option}
              aria-selected={binding.value.pipe(
                map((v: readonly unknown[]) => !!v?.includes(value)),
                distinctUntilChanged(),
              )}
              onClick={() => {
                const v: unknown[] = binding.get() ?? [];
                binding.set(v.includes(value) ? v.filter((x) => x !== value) : [...v, value]);
              }}
            >
              <span>{label}</span>
              <Icon name="check" size={14} class={classes.check} />
            </button>
          ))}
        </ScrollArea>
      )}
    />
  );
}

/** Three numbers in a row (x, y, z), each its own small field. */
export function Vec3Control(props: { param: PD.Vec3; binding: FieldBinding }) {
  const { binding } = props;
  return (
    <div class={classes.vec}>
      {(['x', 'y', 'z'] as const).map((axis, i) => (
        <input
          class={form.number}
          aria-label={axis}
          data-tip={axis}
          value={binding.value.pipe(map((v: number[]) => trim(v?.[i])))}
          onChange={(e) => {
            const n = Number(e.currentTarget.value);
            if (!Number.isFinite(n)) return;
            const v = [...(binding.get() ?? [0, 0, 0])];
            v[i] = n;
            binding.set(v);
          }}
        />
      ))}
    </div>
  );
}

function trim(v: number | undefined) {
  return v === undefined || !Number.isFinite(v) ? '' : String(Math.round(v * 1e4) / 1e4);
}

/** A 4 × 4 matrix: "Identity" or "Matrix" in a field; opened, its 16 numbers as a grid under the row (row-major). */
export function Mat4Row(props: { head: Child; label: string; param: PD.Mat4; binding: FieldBinding }) {
  const { binding } = props;
  const isIdentity = (m: ArrayLike<number>) => {
    for (let i = 0; i < 16; i++) if (m[i] !== (i % 5 === 0 ? 1 : 0)) return false;
    return true;
  };
  return (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          {binding.value.pipe(map((m: ArrayLike<number>) => (m && isIdentity(m) ? 'Identity' : '4 × 4 matrix')))}
        </Field>
      )}
      body={() => (
        <div class={classes.matrix}>
          {Array.from({ length: 16 }, (_, k) => {
            // Shown by rows; stored by columns (Mol*'s Mat4).
            const index = (k % 4) * 4 + Math.floor(k / 4);
            return (
              <input
                class={form.number}
                value={binding.value.pipe(map((m: ArrayLike<number>) => trim(m?.[index])))}
                onChange={(e) => {
                  const n = Number(e.currentTarget.value);
                  if (!Number.isFinite(n)) return;
                  const m = Array.from(binding.get() as ArrayLike<number>);
                  m[index] = n;
                  binding.set(m);
                }}
              />
            );
          })}
        </div>
      )}
    />
  );
}

/** An interval: one track, two knobs (the low and the high end), the range on the right. */
export function IntervalControl(props: { param: PD.Interval; binding: FieldBinding }) {
  const { param, binding } = props;
  const min = param.min ?? 0;
  const max = param.max ?? 1;
  const step = param.step ?? (max - min) / 100;
  const pct = (v: number) => `${(100 * (v - min)) / (max - min || 1)}%`;
  const digits = step < 1 ? Math.min(4, Math.ceil(-Math.log10(step))) : 0;
  const end = (i: 0 | 1) => (
    <input
      class={classes.intervalRange}
      type="range"
      min={String(min)}
      max={String(max)}
      step={String(step)}
      aria-label={i === 0 ? 'From' : 'To'}
      value={binding.value.pipe(map((v: [number, number]) => String(v?.[i] ?? (i ? max : min))))}
      onInput={(e) => {
        const n = Number(e.currentTarget.value);
        const v: [number, number] = [...(binding.get() ?? [min, max])] as [number, number];
        // The ends cannot cross.
        v[i] = i === 0 ? Math.min(n, v[1]) : Math.max(n, v[0]);
        binding.set(v);
        e.currentTarget.value = String(v[i]);
      }}
    />
  );
  return (
    <div class={classes.interval}>
      <div
        class={classes.intervalTrack}
        style={{
          '--msp-from': binding.value.pipe(map((v: [number, number]) => pct(v?.[0] ?? min))),
          '--msp-to': binding.value.pipe(map((v: [number, number]) => pct(v?.[1] ?? max))),
        }}
      >
        {end(0)}
        {end(1)}
      </div>
      <span class={classes.intervalValue}>
        {binding.value.pipe(
          map((v: [number, number]) => `${(v?.[0] ?? min).toFixed(digits)}–${(v?.[1] ?? max).toFixed(digits)}`),
        )}
      </span>
    </div>
  );
}

/** A color list as a CSS gradient: interpolated, or a set of hard bands. */
export function colorListGradient(value: { kind: 'interpolate' | 'set'; colors: readonly unknown[] } | undefined) {
  const colors = value?.colors ?? [];
  if (!colors.length) return 'none';
  const entries = colors.map((c, i) =>
    Array.isArray(c)
      ? { color: colorToHex(c[0]), at: c[1] as number }
      : { color: colorToHex(c as number), at: colors.length === 1 ? 0 : i / (colors.length - 1) },
  );
  if (value?.kind === 'set') {
    const n = entries.length;
    return `linear-gradient(90deg, ${entries.map((e, i) => `${e.color} ${(100 * i) / n}% ${(100 * (i + 1)) / n}%`).join(', ')})`;
  }
  return `linear-gradient(90deg, ${entries.map((e) => `${e.color} ${100 * e.at}%`).join(', ')})`;
}

const Presets = Object.entries(ColorLists) as [ColorListName, (typeof ColorLists)[ColorListName]][];

/**
 * A color list: its colors as a strip in a field (and the preset's name when it is one); opened, the presets under the
 * row, each its strip and name, the current one checked.
 */
export function ColorListRow(props: { head: Child; label: string; param: PD.ColorList; binding: FieldBinding }) {
  const { param, binding } = props;
  const kind = param.presetKind ?? 'all';
  const presets = Presets.filter(([, l]) =>
    kind === 'all' ? true : kind === 'set' ? l.type === 'qualitative' : l.type !== 'qualitative',
  );
  const nameOf = (v: { colors: readonly unknown[] } | undefined) =>
    presets.find(([, l]) => JSON.stringify(l.list) === JSON.stringify(v?.colors))?.[1].label ?? 'Custom';
  return (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          <span class={classes.strip} style={{ background: binding.value.pipe(map(colorListGradient)) }} />
          <span class={classes.stripName}>{binding.value.pipe(map(nameOf))}</span>
        </Field>
      )}
      body={(close) => (
        <ScrollArea
          class={classes.options}
          contentClass={classes.optionsContent}
          role="listbox"
          aria-label={props.label}
        >
          {presets.map(([, preset]) => {
            const value = { kind: preset.type === 'qualitative' ? 'set' : 'interpolate', colors: preset.list } as const;
            return (
              <button
                type="button"
                role="option"
                class={classes.option}
                aria-selected={binding.value.pipe(map((v) => nameOf(v) === preset.label))}
                onClick={() => {
                  binding.set({ kind: value.kind, colors: [...value.colors] });
                  close();
                }}
              >
                <span class={classes.strip} style={{ background: colorListGradient(value) }} />
                <span class={classes.stripName}>{preset.label}</span>
                <Icon name="check" size={14} class={classes.check} />
              </button>
            );
          })}
        </ScrollArea>
      )}
    />
  );
}

/** A URL: its address as text. */
export function UrlControl(props: { binding: FieldBinding }) {
  const { binding } = props;
  return (
    <input
      class={form.input}
      placeholder="https://…"
      value={binding.value.pipe(map((v: Asset.Url | undefined) => v?.url ?? ''))}
      onChange={(e) => binding.set(Asset.Url(e.currentTarget.value))}
    />
  );
}

/**
 * A line graph (a volume's transfer function): a small preview of the curve in a field, the number of points; opened,
 * an editor under the row: drag a point (between its neighbours), click to add one, double-click one to remove it (two
 * stay).
 */
export function LineGraphRow(props: { head: Child; label: string; param: PD.LineGraph; binding: FieldBinding }) {
  const { binding } = props;
  const pathOf = (pts: readonly (readonly [number, number])[] | undefined, w: number, h: number) =>
    (pts ?? []).map(([x, y], i) => `${i ? 'L' : 'M'}${(x * w).toFixed(1)} ${((1 - y) * h).toFixed(1)}`).join(' ');
  return (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          <svg class={classes.spark} viewBox="0 0 56 14" preserveAspectRatio="none" aria-hidden="true">
            <path d={binding.value.pipe(map((p) => pathOf(p, 56, 14)))} />
          </svg>
          <span class={classes.stripName}>{binding.value.pipe(map((p: unknown[]) => `${p?.length ?? 0} points`))}</span>
        </Field>
      )}
      body={() => <LineGraphEditor binding={binding} pathOf={pathOf} />}
    />
  );
}

function LineGraphEditor(props: {
  binding: FieldBinding;
  pathOf: (p: readonly (readonly [number, number])[] | undefined, w: number, h: number) => string;
}) {
  const { binding } = props;
  const W = 240;
  const H = 96;
  const hover = observableState('');
  const points = binding.value as Observable<[number, number][]>;
  const svg = (
    <svg class={classes.graph} viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="none">
      {[0.25, 0.5, 0.75].map((t) => (
        <line class={classes.grid} x1={String(t * W)} x2={String(t * W)} y1="0" y2={String(H)} />
      ))}
      {[0.25, 0.5, 0.75].map((t) => (
        <line class={classes.grid} y1={String(t * H)} y2={String(t * H)} x1="0" x2={String(W)} />
      ))}
      <path class={classes.curve} d={binding.value.pipe(map((p) => props.pathOf(p, W, H)))} />
      <For each={points.pipe(map((p) => (p ?? []).map((_, i) => i)))}>
        {(i) => (
          <circle
            class={classes.point}
            r="4"
            cx={points.pipe(map((p) => String((p?.[i]?.[0] ?? 0) * W)))}
            cy={points.pipe(map((p) => String((1 - (p?.[i]?.[1] ?? 0)) * H)))}
            data-index={String(i)}
          />
        )}
      </For>
    </svg>
  ) as SVGSVGElement;
  // Pointer: drag a point, click empty space to add, double-click a point to remove (not the ends).
  const at = (e: PointerEvent | MouseEvent): [number, number] => {
    const r = svg.getBoundingClientRect();
    const x = Math.min(1, Math.max(0, (e.clientX - r.left) / r.width));
    const y = Math.min(1, Math.max(0, 1 - (e.clientY - r.top) / r.height));
    return [Math.round(x * 1000) / 1000, Math.round(y * 1000) / 1000];
  };
  let dragging = -1;
  const onDown = (e: PointerEvent) => {
    const t = e.target as Element;
    const pts: [number, number][] = [...(binding.get() ?? [])];
    if (t.tagName === 'circle') {
      dragging = Number((t as SVGElement).dataset.index);
      svg.setPointerCapture(e.pointerId);
      return;
    }
    const p = at(e);
    const i = pts.findIndex(([x]) => x > p[0]);
    pts.splice(i < 0 ? pts.length : i, 0, p);
    binding.set(pts);
  };
  const onMove = (e: PointerEvent) => {
    const p = at(e);
    hover.set(`${p[0].toFixed(2)}, ${p[1].toFixed(2)}`);
    if (dragging < 0) return;
    const pts: [number, number][] = [...(binding.get() ?? [])];
    // A point stays between its neighbours (the curve keeps its order).
    const lo = dragging > 0 ? pts[dragging - 1][0] : 0;
    const hi = dragging < pts.length - 1 ? pts[dragging + 1][0] : 1;
    pts[dragging] = [Math.min(hi, Math.max(lo, p[0])), p[1]];
    binding.set(pts);
  };
  const onUp = () => {
    dragging = -1;
  };
  const onDouble = (e: MouseEvent) => {
    const t = e.target as Element;
    if (t.tagName !== 'circle') return;
    const i = Number((t as SVGElement).dataset.index);
    const pts: [number, number][] = [...(binding.get() ?? [])];
    // A curve keeps two points at least.
    if (pts.length <= 2) return;
    pts.splice(i, 1);
    binding.set(pts);
  };
  svg.addEventListener('pointerdown', onDown);
  svg.addEventListener('pointermove', onMove);
  svg.addEventListener('pointerup', onUp);
  svg.addEventListener('pointercancel', onUp);
  svg.addEventListener('dblclick', onDouble);
  svg.addEventListener('pointerleave', () => hover.set(''));
  onCleanup(() => {
    dragging = -1;
  });
  return (
    <div class={classes.graphEditor}>
      {/* The position under the pointer in the graph's corner: the hint below never changes the editor's height. */}
      <div class={classes.graphArea}>
        {svg}
        <span class={classes.graphAt}>{hover.asReadonly}</span>
      </div>
      <p class={classes.graphHint}>Drag a point · click to add · double-click to remove</p>
    </div>
  );
}
