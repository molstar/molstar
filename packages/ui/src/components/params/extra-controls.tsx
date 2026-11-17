/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, KeyedSelector, type State } from '@molstar/core/reactive';
import { Asset } from '@molstar/core/util/assets';
import { type ColorListName, ColorLists } from '@molstar/core/util/color/lists';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  activeElement,
  type Child,
  drag,
  eventTarget,
  For,
  frameState,
  on,
  Show,
  setAttr,
  setStyle,
  setText,
  watch,
} from '../../dom/index.js';
import { state } from '../../state.js';
import { colorToHex } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { dismissable } from '../popover/dismiss.js';
import { Popover } from '../popover/popover.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes as select } from '../select/_generated/select.module.css.js';
import { Floating } from '../select/select.js';
import { classes } from './_generated/extra-controls.module.css.js';
import { classes as form } from './_generated/params-form.module.css.js';
import type { FieldBinding } from './field-binding.js';
import { formatNumber, frameInput, holdNumber, NumberInput, setValue } from './number-input.js';
import { PresetsButton, presetIndex } from './presets.js';

/** A form row: the label (its text and help) on the label column, the control on the value column (styles §5). */
export function Row(props: { head: Child; children: Child }) {
  return (
    <div class={form.row} data-field="">
      <span class={form.label}>{props.head}</span>
      <div class={form.control}>{props.children}</div>
    </div>
  );
}

/**
 * A field like a select's (its content at the control inset, a chevron at its end) that opens something over the form:
 * a list over the field (`listbox`), or a panel under it (`dialog`). Toggles `open`.
 */
function Field(props: {
  open: State<boolean>;
  label: string;
  /** Labelled by these ids instead (a row's label and the field's value). */
  labelledBy?: string;
  popup: 'listbox' | 'dialog';
  ref?: (el: HTMLElement) => void;
  children: Child;
}) {
  const el = (
    <button
      type="button"
      class={classes.field}
      aria-label={props.labelledBy ? undefined : props.label}
      aria-labelledby={props.labelledBy}
      aria-haspopup={props.popup}
      onClick={() => props.open.set(!props.open.value)}
    >
      <span class={classes.fieldContent}>{props.children}</span>
      <Icon name="chevron-down" class={classes.chevron} />
    </button>
  ) as HTMLElement;
  props.ref?.(el);
  watch([props.open], () => setAttr(el, 'aria-expanded', String(props.open.value)));
  return el;
}

/** A field named by its row's label and its value (`<labelId> <labelId>-value`, so the value is heard with the name). */
const labelledBy = (labelId?: string) => (labelId ? `${labelId} ${labelId}-value` : undefined);
/** The id of the field's value, the second part of its name. */
const valueId = (labelId?: string) => (labelId ? `${labelId}-value` : undefined);

/** A row whose field opens a panel under it (a matrix, a curve's editor, a color picker), over the form. */
export function PanelRow(props: {
  head: Child;
  label: string;
  labelledBy?: string;
  field: Child;
  panel: (close: () => void) => Child;
}) {
  const open = state(false);
  let anchor!: HTMLElement;
  return (
    <>
      <Row head={props.head}>
        <Field open={open} label={props.label} labelledBy={props.labelledBy} popup="dialog" ref={(el) => (anchor = el)}>
          {props.field}
        </Field>
      </Row>
      <Show when={open}>
        {() => (
          <Popover
            anchor={anchor}
            align="start"
            aria-label={props.label}
            onClose={(refocus) => {
              open.set(false);
              if (refocus) anchor.focus();
            }}
          >
            <div class={classes.panel}>{props.panel(() => open.set(false))}</div>
          </Popover>
        )}
      </Show>
    </>
  );
}

/**
 * Options over their field, as a select's (the current one on the field's text; the select's own look, so the two read
 * as one): arrows move, Escape closes, a press outside (not on the field) closes. `multiple`: picks toggle, the list
 * stays.
 */
function FieldOptions(props: {
  anchor: HTMLElement;
  label: string;
  multiple?: boolean;
  onClose(): void;
  children: Child;
}) {
  let list!: HTMLElement;
  const options = () => [...list.querySelectorAll<HTMLElement>('[role="option"]')];
  const focus = (i: number) => {
    const all = options();
    all[Math.max(0, Math.min(all.length - 1, i))]?.focus();
  };
  // Closes as every overlay does (popover/dismiss.ts): a press outside, Escape, Tab; the focus back on the field.
  dismissable({
    layer: () => list,
    anchor: props.anchor,
    close: () => props.onClose(),
    refocus: () => props.anchor,
    tab: true,
  });
  return (
    <Floating anchor={props.anchor} placement="over">
      <ScrollArea
        class={select.list}
        contentClass={select.listContent}
        role="listbox"
        aria-label={props.label}
        aria-multiselectable={props.multiple ? 'true' : undefined}
        ref={(el: HTMLElement) => (list = el)}
        onKeyDown={(e) => {
          const all = options();
          const i = all.indexOf(activeElement(list) as HTMLElement);
          if (e.key === 'ArrowDown') focus(i + 1);
          else if (e.key === 'ArrowUp') focus(i - 1);
          else if (e.key === 'Home') focus(0);
          else if (e.key === 'End') focus(all.length - 1);
          else return;
          e.preventDefault();
          e.stopPropagation();
        }}
      >
        {props.children}
      </ScrollArea>
    </Floating>
  );
}

/**
 * Multi-select: the chosen options' names in a select-like field ("None" when empty, the names fading when long: they
 * are data); opened, the options over the field, each toggling (checked when chosen), the list staying open.
 */
export function MultiSelectRow(props: {
  head: Child;
  label: string;
  /** The row label's id: the field is named by it and its value. */
  labelId?: string;
  param: PD.MultiSelect<any>;
  binding: FieldBinding;
}) {
  const { param, binding } = props;
  const open = state(false);
  const names = document.createTextNode('');
  let anchor!: HTMLElement;
  watch([binding.value], () => {
    const v = binding.value.value as readonly unknown[] | undefined;
    const chosen = param.options.filter(([o]) => v?.includes(o)).map(([, l]) => l);
    setText(names, chosen.length ? chosen.join(', ') : (param.emptyValue ?? 'None'));
  });
  return (
    <>
      <Row head={props.head}>
        <Field
          open={open}
          label={props.label}
          labelledBy={labelledBy(props.labelId)}
          popup="listbox"
          ref={(el) => (anchor = el)}
        >
          <span class={classes.value} id={valueId(props.labelId)}>
            {names}
          </span>
        </Field>
      </Row>
      <Show when={open}>
        {() => {
          const buttons: HTMLElement[] = [];
          const list = (
            <FieldOptions anchor={anchor} label={props.label} multiple onClose={() => open.set(false)}>
              {param.options.map(([value, label]) => (
                <button
                  type="button"
                  role="option"
                  class={select.option}
                  ref={(el: HTMLElement) => buttons.push(el)}
                  onClick={() => {
                    const v: unknown[] = binding.value.value ?? [];
                    binding.set(v.includes(value) ? v.filter((x) => x !== value) : [...v, value]);
                  }}
                >
                  <span class={select.optionLabel}>{label}</span>
                  <Icon name="check" class={select.check} />
                </button>
              ))}
            </FieldOptions>
          );
          // The marks: one watch over the value (a set, so no single key to select by).
          watch([binding.value], () => {
            const v = binding.value.value as readonly unknown[] | undefined;
            for (let i = 0; i < buttons.length; i++)
              setAttr(buttons[i], 'aria-selected', String(!!v?.includes(param.options[i][0])));
          });
          return list;
        }}
      </Show>
    </>
  );
}

/** Three numbers in a row (x, y, z), each its own small field, sharing the value column; each keeps to the range. */
export function Vec3Control(props: { param: PD.Vec3; binding: FieldBinding; label?: string }) {
  const { param, binding } = props;
  const inputs = (['x', 'y', 'z'] as const).map((axis, i) => (
    <NumberInput
      class={form.number}
      aria-label={props.label ? `${props.label} ${axis}` : axis}
      data-tip={axis}
      range={param}
      current={() => binding.value.value?.[i]}
      commit={(n) => {
        const v = [...(binding.value.value ?? [0, 0, 0])];
        v[i] = n;
        binding.set(v);
      }}
    />
  )) as HTMLInputElement[];
  watch([binding.value], () => {
    const v = binding.value.value as number[] | undefined;
    for (let i = 0; i < inputs.length; i++) setValue(inputs[i], formatNumber(v?.[i], param.step));
  });
  return <div class={classes.vec}>{inputs}</div>;
}

function isIdentity(m: ArrayLike<number>) {
  for (let i = 0; i < 16; i++) if (m[i] !== (i % 5 === 0 ? 1 : 0)) return false;
  return true;
}

/** A 4 × 4 matrix: "Identity" or "4 × 4 matrix" in a field; opened, its 16 numbers as a grid in a panel (row-major). */
export function Mat4Row(props: {
  head: Child;
  label: string;
  /** The row label's id: the field is named by it and its value (`<labelId> <labelId>-value`). */
  labelId?: string;
  param: PD.Mat4;
  binding: FieldBinding;
}) {
  const { binding, labelId } = props;
  const summary = document.createTextNode('');
  watch([binding.value], () => {
    const m = binding.value.value as ArrayLike<number> | undefined;
    setText(summary, m && isIdentity(m) ? 'Identity' : '4 × 4 matrix');
  });
  return (
    <PanelRow
      head={props.head}
      label={props.label}
      labelledBy={labelledBy(labelId)}
      field={
        <span class={classes.value} id={valueId(labelId)}>
          {summary}
        </span>
      }
      panel={() => {
        // Shown by rows; stored by columns (Mol*'s Mat4).
        const inputs = Array.from({ length: 16 }, (_, k) => {
          const index = (k % 4) * 4 + Math.floor(k / 4);
          return (
            <NumberInput
              class={form.number}
              aria-label={`${props.label}, row ${Math.floor(k / 4) + 1}, column ${(k % 4) + 1}`}
              current={() => (binding.value.value as ArrayLike<number> | undefined)?.[index]}
              commit={(n) => {
                const m = Array.from(binding.value.value as ArrayLike<number>);
                m[index] = n;
                binding.set(m);
              }}
            />
          ) as HTMLInputElement;
        });
        const grid = <div class={classes.matrix}>{inputs}</div>;
        // One watch for the 16 numbers.
        watch([binding.value], () => {
          const m = binding.value.value as ArrayLike<number> | undefined;
          for (let k = 0; k < 16; k++) setValue(inputs[k], formatNumber(m?.[(k % 4) * 4 + Math.floor(k / 4)]));
        });
        return grid;
      }}
    />
  );
}

/**
 * An interval: one track, two knobs (the low and the high end), the range after it (as a slider's number). Without both
 * bounds there is no track to place the knobs on: two number fields (from, to), each keeping to the bound it has and to
 * the other end.
 */
export function IntervalControl(props: { param: PD.Interval; binding: FieldBinding; label?: string }) {
  const { param, binding } = props;
  const name = (i: 0 | 1) => `${props.label ? `${props.label} ` : ''}${i === 0 ? 'from' : 'to'}`;
  const value = () => binding.value.value as [number, number] | undefined;
  if (param.min === undefined || param.max === undefined) {
    // The ends cannot cross: each is bounded by the other.
    const ranges = [
      () => ({ min: param.min, max: value()?.[1] ?? param.max, step: param.step }),
      () => ({ min: value()?.[0] ?? param.min, max: param.max, step: param.step }),
    ];
    const fields = ([0, 1] as const).map(
      (i) =>
        (
          <NumberInput
            class={form.number}
            aria-label={name(i)}
            range={ranges[i]}
            current={() => value()?.[i]}
            commit={(n) => {
              const v = [...(value() ?? param.defaultValue)] as [number, number];
              v[i] = n;
              binding.set(v);
            }}
          />
        ) as HTMLInputElement,
    );
    watch([binding.value], () => {
      const v = value();
      for (const i of [0, 1] as const) setValue(fields[i], formatNumber(v?.[i], param.step));
    });
    return <div class={classes.pair}>{fields}</div>;
  }
  const min = param.min;
  const max = param.max;
  const step = param.step ?? (max - min) / 100;
  const pct = (v: number) => `${(100 * (v - min)) / (max - min || 1)}%`;
  const ends: HTMLInputElement[] = [];
  let track!: HTMLElement;
  const range = document.createTextNode('');
  // A drag in the frame lane (dom §3), as a slider's: the field is set once a frame.
  const drag = frameInput<[number, number]>(
    (v) => binding.set(v),
    (v) => {
      const c = value();
      return !!c && c[0] === v[0] && c[1] === v[1];
    },
  );
  const end = (i: 0 | 1) => (
    <input
      class={classes.intervalRange}
      type="range"
      min={String(min)}
      max={String(max)}
      step={String(step)}
      aria-label={name(i)}
      ref={(el: HTMLInputElement) => (ends[i] = el)}
      onInput={(e) => {
        const n = holdNumber(Number(e.currentTarget.value), param);
        // The other end as it shows (the field may not have it yet: it is set in the frame).
        const other = Number(ends[1 - i].value);
        // The ends cannot cross.
        const v: [number, number] = i === 0 ? [Math.min(n, other), other] : [other, Math.max(n, other)];
        drag(v);
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
    const v = value();
    const from = v?.[0] ?? min;
    const to = v?.[1] ?? max;
    setValue(ends[0], String(from));
    setValue(ends[1], String(to));
    setStyle(track, '--msui-params-from', pct(from));
    setStyle(track, '--msui-params-to', pct(to));
    setText(range, `${formatNumber(from, step)}–${formatNumber(to, step)}`);
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

const ColorListPresets = Object.entries(ColorLists) as [ColorListName, (typeof ColorLists)[ColorListName]][];

/**
 * A color list: its colors as a strip in a select-like field, and the preset's name when it is one; opened, the presets
 * over the field (a select's options: each its strip and name, the current one checked).
 */
export function ColorListRow(props: {
  head: Child;
  label: string;
  /** The row label's id: the field is named by it and its value (the preset's name). */
  labelId?: string;
  param: PD.ColorList;
  binding: FieldBinding;
}) {
  const { param, binding } = props;
  const kind = param.presetKind ?? 'all';
  const presets = ColorListPresets.filter(([, l]) =>
    kind === 'all' ? true : kind === 'set' ? l.type === 'qualitative' : l.type !== 'qualitative',
  );
  const nameOf = (v: { colors: readonly unknown[] } | undefined) =>
    presets.find(([, l]) => JSON.stringify(l.list) === JSON.stringify(v?.colors))?.[1].label ?? 'Custom';
  const open = state(false);
  let strip!: HTMLElement;
  let anchor!: HTMLElement;
  const name = document.createTextNode('');
  // The current preset's name: the options are marked by one selector over it (a change rewrites two options).
  const current = derive(binding.value, nameOf);
  const content = (
    <>
      <Row head={props.head}>
        <Field
          open={open}
          label={props.label}
          labelledBy={labelledBy(props.labelId)}
          popup="listbox"
          ref={(el) => (anchor = el)}
        >
          <span class={classes.strip} ref={(el: HTMLElement) => (strip = el)} />
          <span class={classes.value} id={valueId(props.labelId)}>
            {name}
          </span>
        </Field>
      </Row>
      <Show when={open}>
        {() => {
          const selected = new KeyedSelector(current);
          return (
            <FieldOptions anchor={anchor} label={props.label} onClose={() => open.set(false)}>
              {presets.map(([, preset]) => {
                const value = {
                  kind: preset.type === 'qualitative' ? 'set' : 'interpolate',
                  colors: preset.list,
                } as const;
                const button = (
                  <button
                    type="button"
                    role="option"
                    class={select.option}
                    onClick={() => {
                      binding.set({ kind: value.kind, colors: [...value.colors] });
                      open.set(false);
                    }}
                  >
                    <span class={classes.optionContent}>
                      <span class={classes.strip} style={{ background: colorListGradient(value) }} />
                      <span class={select.optionLabel}>{preset.label}</span>
                    </span>
                    <Icon name="check" class={select.check} />
                  </button>
                ) as HTMLElement;
                const on = selected.key(preset.label);
                watch([on], () => setAttr(button, 'aria-selected', String(on.value)));
                return button;
              })}
            </FieldOptions>
          );
        }}
      </Show>
    </>
  );
  watch([binding.value], () => {
    setStyle(strip, 'background', colorListGradient(binding.value.value));
    setText(name, current.value);
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
  watch([binding.value], () => setValue(el, (binding.value.value as Asset.Url | undefined)?.url ?? ''));
  return el;
}

type Points = readonly (readonly [number, number])[];

function pathOf(pts: Points | undefined, w: number, h: number) {
  return (pts ?? []).map(([x, y], i) => `${i ? 'L' : 'M'}${(x * w).toFixed(1)} ${((1 - y) * h).toFixed(1)}`).join(' ');
}

/**
 * A line graph (a volume's transfer function): a small preview of the curve in a field and the number of points;
 * opened, an editor in a panel under it: drag a point (between its neighbours), click to add one, double-click one to
 * remove it (two stay).
 */
export function LineGraphRow(props: {
  head: Child;
  label: string;
  /** The row label's id: the field is named by it and its value (the number of points). */
  labelId?: string;
  param: PD.LineGraph;
  binding: FieldBinding;
}) {
  const { binding } = props;
  let spark!: SVGPathElement;
  const count = document.createTextNode('');
  const content = (
    <PanelRow
      head={props.head}
      label={props.label}
      labelledBy={labelledBy(props.labelId)}
      field={
        <>
          <svg class={classes.spark} viewBox="0 0 48 12" preserveAspectRatio="none" aria-hidden="true">
            <path ref={(el: SVGElement) => (spark = el as SVGPathElement)} />
          </svg>
          <span class={classes.value} id={valueId(props.labelId)}>
            {count}
          </span>
        </>
      }
      panel={() => <LineGraphEditor param={props.param} binding={binding} label={props.label} />}
    />
  );
  watch([binding.value], () => {
    const p = binding.value.value as [number, number][] | undefined;
    setAttr(spark, 'd', pathOf(p, 48, 12));
    setText(count, `${p?.length ?? 0} points`);
  });
  return content;
}

function LineGraphEditor(props: { param: PD.LineGraph; binding: FieldBinding; label: string }) {
  const { binding } = props;
  const W = 240;
  const H = 96;
  // The pointer's position is input: a frame state, one write a frame (dom §3).
  const hover = frameState('');
  const hoverText = document.createTextNode('');
  const points = binding.value;
  // Presets come from the data at hand (a volume's statistics): asked for when opened.
  const getPresets = props.param.getPresets;
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
      <For each={derive(points, (p: Points | undefined) => (p ?? []).map((_, i) => i))}>
        {(i) => {
          // The point's own position: the list's other points change nothing here (a drag moves one tuple).
          const point = derive(points, (p: Points | undefined) => p?.[i]);
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
  watch([points], () => setAttr(curve, 'd', pathOf(points.value, W, H)));
  watch([hover], () => setText(hoverText, hover.value));
  // Pointer: drag a point, click empty space to add, double-click a point to remove (not below two). The drag is
  // `drag`'s (dom §6: attached once, the pointer captured, every way it ends ends it); the hover's position is its own.
  const at = (e: PointerEvent | MouseEvent): [number, number] => {
    const r = svg.getBoundingClientRect();
    const x = Math.min(1, Math.max(0, (e.clientX - r.left) / r.width));
    const y = Math.min(1, Math.max(0, 1 - (e.clientY - r.top) / r.height));
    return [Math.round(x * 1000) / 1000, Math.round(y * 1000) / 1000];
  };
  drag<number>(svg, {
    start(e) {
      const t = eventTarget(e) as Element;
      // A point: dragged.
      if (t.tagName === 'circle') return Number((t as SVGElement).dataset.index);
      // Empty space: a point added there (no drag).
      const pts: [number, number][] = [...(points.value ?? [])];
      const p = at(e);
      const i = pts.findIndex(([x]) => x > p[0]);
      pts.splice(i < 0 ? pts.length : i, 0, p);
      binding.set(pts);
      return undefined;
    },
    move(e, dragging) {
      const p = at(e);
      const pts: [number, number][] = [...(points.value ?? [])];
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
  on<MouseEvent>(svg, 'dblclick', (e) => {
    const t = eventTarget(e) as Element;
    if (t.tagName !== 'circle') return;
    const i = Number((t as SVGElement).dataset.index);
    const pts: [number, number][] = [...(points.value ?? [])];
    // A curve keeps two points at least.
    if (pts.length <= 2) return;
    pts.splice(i, 1);
    binding.set(pts);
  });
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
        {getPresets ? (
          <PresetsButton
            label={props.label}
            presets={getPresets}
            current={derive(points, (v) => presetIndex(getPresets(), v))}
            onPick={(v) => binding.set(v)}
          />
        ) : null}
      </div>
    </div>
  );
}
