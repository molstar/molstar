/** @jsxImportSource @molstar/prism-ui/dom */

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
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { frameState, watch } from '../../dom/scope.js';
import { drag, on } from '../../dom/scoped.js';
import { setAttr, setStyle, setText } from '../../dom/writes.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
import { colorToHex } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/extra-controls.module.css.js';
import { classes as form } from './_generated/params-form.module.css.js';
import type { FieldBinding } from './field-binding.js';
import { PresetsButton, PresetsList, presetIndex } from './presets.js';

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
  const el = (
    <button
      type="button"
      class={classes.field}
      aria-label={props.label}
      onClick={() => props.open.set(!props.open.value)}
    >
      <span class={classes.fieldContent}>{props.children}</span>
      <Icon name="chevron-down" size={14} class={classes.chevron} />
    </button>
  ) as HTMLElement;
  watch([props.open], () => setAttr(el, 'aria-expanded', String(props.open.value)));
  return el;
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
  const names = document.createTextNode('');
  watch([binding.value], () => {
    const v = binding.get() as readonly unknown[] | undefined;
    const chosen = param.options.filter(([o]) => v?.includes(o)).map(([, l]) => l);
    setText(names, chosen.length ? chosen.join(', ') : (param.emptyValue ?? 'None'));
  });
  return (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          {names}
        </Field>
      )}
      body={() => {
        const buttons: HTMLElement[] = [];
        const list = (
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
                ref={(el: HTMLElement) => buttons.push(el)}
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
        );
        // One watch for the options' marks.
        watch([binding.value], () => {
          const v = binding.get() as readonly unknown[] | undefined;
          for (let i = 0; i < buttons.length; i++)
            setAttr(buttons[i], 'aria-selected', String(!!v?.includes(param.options[i][0])));
        });
        return list;
      }}
    />
  );
}

/** Three numbers in a row (x, y, z), each its own small field. */
export function Vec3Control(props: { param: PD.Vec3; binding: FieldBinding; label?: string }) {
  const { binding } = props;
  const inputs: HTMLInputElement[] = [];
  const content = (
    <div class={classes.vec}>
      {(['x', 'y', 'z'] as const).map((axis, i) => (
        <input
          class={form.number}
          aria-label={props.label ? `${props.label} ${axis}` : axis}
          data-tip={axis}
          ref={(el: HTMLInputElement) => inputs.push(el)}
          onChange={(e) => {
            const input = e.currentTarget;
            const n = Number(input.value);
            const v = [...(binding.get() ?? [0, 0, 0])];
            if (Number.isFinite(n)) {
              v[i] = n;
              binding.set(v);
            }
            // The field says the value in its form, though it formats to the text that was shown ("1.50" → "1.5").
            input.value = trim(v[i]);
          }}
        />
      ))}
    </div>
  );
  watch([binding.value], () => {
    const v = binding.get() as number[] | undefined;
    for (let i = 0; i < inputs.length; i++) setValue(inputs[i], trim(v?.[i]));
  });
  return content;
}

/** An input's value, written only when it differs (a write moves the caret). */
function setValue(el: HTMLInputElement, value: string) {
  if (el.value !== value) el.value = value;
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
  const summary = document.createTextNode('');
  watch([binding.value], () => {
    const m = binding.get() as ArrayLike<number> | undefined;
    setText(summary, m && isIdentity(m) ? 'Identity' : '4 × 4 matrix');
  });
  return (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          {summary}
        </Field>
      )}
      body={() => {
        const inputs: HTMLInputElement[] = [];
        const grid = (
          <div class={classes.matrix}>
            {Array.from({ length: 16 }, (_, k) => {
              // Shown by rows; stored by columns (Mol*'s Mat4).
              const index = (k % 4) * 4 + Math.floor(k / 4);
              return (
                <input
                  class={form.number}
                  aria-label={`${props.label}, row ${Math.floor(k / 4) + 1}, column ${(k % 4) + 1}`}
                  ref={(el: HTMLInputElement) => inputs.push(el)}
                  onChange={(e) => {
                    const input = e.currentTarget;
                    const n = Number(input.value);
                    const m = Array.from(binding.get() as ArrayLike<number>);
                    if (Number.isFinite(n)) {
                      m[index] = n;
                      binding.set(m);
                    }
                    input.value = trim(m[index]);
                  }}
                />
              );
            })}
          </div>
        );
        // One watch for the 16 numbers.
        watch([binding.value], () => {
          const m = binding.get() as ArrayLike<number> | undefined;
          for (let k = 0; k < 16; k++) setValue(inputs[k], trim(m?.[(k % 4) * 4 + Math.floor(k / 4)]));
        });
        return grid;
      }}
    />
  );
}

/** An interval: one track, two knobs (the low and the high end), the range on the right. */
export function IntervalControl(props: { param: PD.Interval; binding: FieldBinding; label?: string }) {
  const { param, binding } = props;
  const min = param.min ?? 0;
  const max = param.max ?? 1;
  const step = param.step ?? (max - min) / 100;
  const pct = (v: number) => `${(100 * (v - min)) / (max - min || 1)}%`;
  const digits = step < 1 ? Math.min(4, Math.ceil(-Math.log10(step))) : 0;
  const ends: HTMLInputElement[] = [];
  let track!: HTMLElement;
  const range = document.createTextNode('');
  const end = (i: 0 | 1) => (
    <input
      class={classes.intervalRange}
      type="range"
      min={String(min)}
      max={String(max)}
      step={String(step)}
      aria-label={`${props.label ? `${props.label} ` : ''}${i === 0 ? 'from' : 'to'}`}
      ref={(el: HTMLInputElement) => (ends[i] = el)}
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
  const content = (
    <div class={classes.interval}>
      <div class={classes.intervalTrack} ref={(el: HTMLElement) => (track = el)}>
        {end(0)}
        {end(1)}
      </div>
      <span class={classes.intervalValue}>{range}</span>
    </div>
  );
  watch([binding.value], () => {
    const v = binding.get() as [number, number] | undefined;
    const from = v?.[0] ?? min;
    const to = v?.[1] ?? max;
    setValue(ends[0], String(v?.[0] ?? min));
    setValue(ends[1], String(v?.[1] ?? max));
    setStyle(track, '--msp-from', pct(from));
    setStyle(track, '--msp-to', pct(to));
    setText(range, `${from.toFixed(digits)}–${to.toFixed(digits)}`);
  });
  return content;
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
  let strip!: HTMLElement;
  const name = document.createTextNode('');
  // The current preset's name: a gate, so the options' marks are written only when the preset changes.
  const current = derive(binding.value, nameOf);
  const content = (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          <span class={classes.strip} ref={(el: HTMLElement) => (strip = el)} />
          <span class={classes.stripName}>{name}</span>
        </Field>
      )}
      body={(close) => {
        const buttons: HTMLElement[] = [];
        const list = (
          <ScrollArea
            class={classes.options}
            contentClass={classes.optionsContent}
            role="listbox"
            aria-label={props.label}
          >
            {presets.map(([, preset]) => {
              const value = {
                kind: preset.type === 'qualitative' ? 'set' : 'interpolate',
                colors: preset.list,
              } as const;
              return (
                <button
                  type="button"
                  role="option"
                  class={classes.option}
                  ref={(el: HTMLElement) => buttons.push(el)}
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
        );
        watch([current], () => {
          for (let i = 0; i < buttons.length; i++)
            setAttr(buttons[i], 'aria-selected', String(current.value === presets[i][1].label));
        });
        return list;
      }}
    />
  );
  watch([binding.value], () => {
    setStyle(strip, 'background', colorListGradient(binding.get()));
    setText(name, nameOf(binding.get()));
  });
  return content;
}

/** A URL: its address as text. */
export function UrlControl(props: { binding: FieldBinding; labelId?: string }) {
  const { binding } = props;
  const el = (
    <input
      class={form.input}
      aria-labelledby={props.labelId}
      placeholder="https://…"
      onChange={(e) => binding.set(Asset.Url(e.currentTarget.value))}
    />
  ) as HTMLInputElement;
  watch([binding.value], () => setValue(el, (binding.get() as Asset.Url | undefined)?.url ?? ''));
  return el;
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
  let spark!: SVGPathElement;
  const count = document.createTextNode('');
  const content = (
    <ExpandingRow
      head={props.head}
      field={(open) => (
        <Field open={open} label={props.label}>
          <svg class={classes.spark} viewBox="0 0 56 14" preserveAspectRatio="none" aria-hidden="true">
            <path ref={(el: SVGElement) => (spark = el as SVGPathElement)} />
          </svg>
          <span class={classes.stripName}>{count}</span>
        </Field>
      )}
      body={() => <LineGraphEditor param={props.param} binding={binding} pathOf={pathOf} />}
    />
  );
  watch([binding.value], () => {
    const p = binding.get() as [number, number][] | undefined;
    setAttr(spark, 'd', pathOf(p, 56, 14));
    setText(count, `${p?.length ?? 0} points`);
  });
  return content;
}

function LineGraphEditor(props: {
  param: PD.LineGraph;
  binding: FieldBinding;
  pathOf: (p: readonly (readonly [number, number])[] | undefined, w: number, h: number) => string;
}) {
  const { binding } = props;
  const W = 240;
  const H = 96;
  // The pointer's position is input: a frame state, one write a frame (rendering.md §2.3).
  const hover = frameState('');
  const hoverText = document.createTextNode('');
  const points = binding.value as ReadonlyState<[number, number][]>;
  // Presets come from the data at hand (a volume's statistics): asked for when opened.
  const getPresets = props.param.getPresets;
  const presetsOpen = observableState(false);
  let curve!: SVGPathElement;
  const svg = (
    <svg class={classes.graph} viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="none">
      {[0.25, 0.5, 0.75].map((t) => (
        <line class={classes.grid} x1={String(t * W)} x2={String(t * W)} y1="0" y2={String(H)} />
      ))}
      {[0.25, 0.5, 0.75].map((t) => (
        <line class={classes.grid} y1={String(t * H)} y2={String(t * H)} x1="0" x2={String(W)} />
      ))}
      <path class={classes.curve} ref={(el: SVGElement) => (curve = el as SVGPathElement)} />
      <For each={derive(points, (p) => (p ?? []).map((_, i) => i))}>
        {(i) => {
          // The point's own position: the list's other points change nothing here (a drag moves one tuple).
          const point = derive(points, (p) => p?.[i]);
          const circle = (<circle class={classes.point} r="4" data-index={String(i)} />) as SVGElement;
          watch([point], () => {
            setAttr(circle, 'cx', String((point.value?.[0] ?? 0) * W));
            setAttr(circle, 'cy', String((1 - (point.value?.[1] ?? 0)) * H));
          });
          return circle;
        }}
      </For>
    </svg>
  ) as SVGSVGElement;
  watch([binding.value], () => setAttr(curve, 'd', props.pathOf(binding.get(), W, H)));
  watch([hover], () => setText(hoverText, hover.value));
  // Pointer: drag a point, click empty space to add, double-click a point to remove (not the ends). The drag is `drag`'s
  // (dom/scoped.ts: attached once, the pointer captured, every way it ends ends it); the hover's position is its own.
  const at = (e: PointerEvent | MouseEvent): [number, number] => {
    const r = svg.getBoundingClientRect();
    const x = Math.min(1, Math.max(0, (e.clientX - r.left) / r.width));
    const y = Math.min(1, Math.max(0, 1 - (e.clientY - r.top) / r.height));
    return [Math.round(x * 1000) / 1000, Math.round(y * 1000) / 1000];
  };
  drag<number>(svg, {
    start(e) {
      const t = e.target as Element;
      // A point: dragged.
      if (t.tagName === 'circle') return Number((t as SVGElement).dataset.index);
      // Empty space: a point added there (no drag).
      const pts: [number, number][] = [...(binding.get() ?? [])];
      const p = at(e);
      const i = pts.findIndex(([x]) => x > p[0]);
      pts.splice(i < 0 ? pts.length : i, 0, p);
      binding.set(pts);
      return undefined;
    },
    move(e, dragging) {
      const p = at(e);
      const pts: [number, number][] = [...(binding.get() ?? [])];
      // A point stays between its neighbours (the curve keeps its order).
      const lo = dragging > 0 ? pts[dragging - 1][0] : 0;
      const hi = dragging < pts.length - 1 ? pts[dragging + 1][0] : 1;
      pts[dragging] = [Math.min(hi, Math.max(lo, p[0])), p[1]];
      binding.set(pts);
    },
  });
  on<PointerEvent>(svg, 'pointermove', (e) => {
    const p = at(e);
    hover.offer(`${p[0].toFixed(2)}, ${p[1].toFixed(2)}`);
  });
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
  on<MouseEvent>(svg, 'dblclick', onDouble);
  on(svg, 'pointerleave', () => hover.offer(''));
  return (
    <div class={classes.graphEditor}>
      {/* The position under the pointer in the graph's corner: the hint below never changes the editor's height. */}
      <div class={classes.graphArea}>
        {svg}
        <span class={classes.graphAt}>{hoverText}</span>
      </div>
      <div class={classes.graphFoot}>
        <p class={classes.graphHint}>Drag a point · click to add · double-click to remove</p>
        {getPresets ? <PresetsButton open={presetsOpen} label="Curve" /> : null}
      </div>
      {getPresets ? (
        <PresetsList
          open={presetsOpen}
          presets={getPresets}
          current={derive(binding.value, (v) => presetIndex(getPresets(), v))}
          onPick={(v) => binding.set(v)}
        />
      ) : null}
    </div>
  );
}
