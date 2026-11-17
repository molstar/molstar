/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, type State } from '@molstar/core/reactive';
import { clamp } from '@molstar/core/math/interpolate';
import { Color } from '@molstar/core/util/color/color';
import { For, type MaybeReactive, Show, setAttr, setClass, setStyle, uniqueId, watch } from '../../dom/index.js';
import { held, state } from '../../state.js';
import { Button } from '../button/button.js';
import { ColorPicker, colorToHex } from '../color/color-picker.js';
import {
  AxisLabels,
  controlClass,
  type Histogram,
  HistogramPlot,
  keySteps,
  PlotHeader,
  plotGesture,
  UnitSwitch,
  unitState,
  type ValueAxis,
  type ValueStats,
  type ValueUnit,
  valueAxis,
  valueNumber,
  valueText,
} from '../histogram/histogram.js';
import { Icon } from '../icon/icon.js';
import { PresetsButton } from '../params/presets.js';
import { Popover } from '../popover/popover.js';
import { classes } from './_generated/transfer-function.module.css.js';

/** A control point: at a value (the values' own units), an opacity (0 to 1) and a color. */
export interface TransferPoint {
  readonly x: number;
  readonly alpha: number;
  readonly color: number;
}

/** A named curve to start from; a point without a color takes the color the current curve has there. */
export interface TransferPreset {
  readonly label: string;
  readonly points: readonly { readonly x: number; readonly alpha: number; readonly color?: number }[];
}

/**
 * Colors for the points: categorical (one color a peak, in order: a peak is a run of points above zero) or sequential
 * (spread along the axis, low to high).
 */
export interface TransferPalette {
  readonly label: string;
  readonly colors: readonly number[];
  readonly sequential?: boolean;
}

export interface TransferFunctionProps {
  /** The header's title (default "Transfer function"); also names the section. */
  title?: string;
  /**
   * The values' statistics: the axis spans [min, max]; σ is counted from the mean. A state or a stream when they change
   * (another map): the control is drawn again for them.
   */
  stats: MaybeReactive<ValueStats>;
  /** Counts of the values; none (yet): an empty well. */
  histogram?: MaybeReactive<Histogram | undefined>;
  /** The control points in order of `x`, with any change made through `onChange` already in (a state the caller sets). */
  points: MaybeReactive<readonly TransferPoint[]>;
  /** The unit values are written in: shared between controls when given; else the control's own (σ). */
  unit?: State<ValueUnit>;
  /**
   * The opacity axis: logarithmic from this up to 1 (opacities per step of a volume are small), its bottom band zero;
   * none: linear from 0 to 1.
   */
  logMin?: number;
  /** Curves to start from (the Presets menu; none: no button). */
  presets?: readonly TransferPreset[];
  /** Palettes to color the points with (the palette button; none: no button). */
  palettes?: readonly TransferPalette[];
  /** The points changed (pointer rate: once a frame while dragged; a key; a color picked). */
  onChange(points: readonly TransferPoint[]): void;
  /** The end of a gesture (a drag let go, a key released, the color closed, a preset): the changes are one step. */
  onCommit?(): void;
}

/** The bottom of a log plot (this fraction of its height) is opacity zero. */
const ZeroBand = 0.04;
/** Two presses on a point within this long are a double click (the pointer is captured on the first). */
const DoublePressMs = 350;
/** A point keeps this much of the axis (a fraction) from its neighbours: the curve keeps its order. */
const Gap = 0.001;

/** The color the points give at a value: interpolated between the points around it, the ends' beyond them. */
export function transferColorAt(points: readonly TransferPoint[], x: number): number {
  if (!points.length) return 0x808080;
  if (x <= points[0].x) return points[0].color;
  for (let i = 1; i < points.length; i++) {
    const a = points[i - 1];
    const b = points[i];
    if (x > b.x) continue;
    const t = b.x > a.x ? (x - a.x) / (b.x - a.x) : 1;
    return Color.interpolate(a.color as Color, b.color as Color, t);
  }
  return points[points.length - 1].color;
}

/** A palette's color at t (0 to 1): sequential, interpolated between its colors. */
function sampleColors(colors: readonly number[], t: number): number {
  if (colors.length < 2) return colors[0] ?? 0x808080;
  const f = clamp(t, 0, 1) * (colors.length - 1);
  const i = Math.min(colors.length - 2, Math.floor(f));
  return Color.interpolate(colors[i] as Color, colors[i + 1] as Color, f - i);
}

/**
 * The points colored by a palette: a sequential one along the axis; a categorical one a color a peak (a run of points
 * above zero), in order, a point at zero taking the next peak's color (the last's after the last peak).
 */
export function applyPalette(
  points: readonly TransferPoint[],
  palette: TransferPalette,
  axis: ValueAxis,
): TransferPoint[] {
  const { colors } = palette;
  if (!colors.length) return [...points];
  if (palette.sequential) return points.map((p) => ({ ...p, color: sampleColors(colors, axis.at(p.x)) }));
  const peakOf: number[] = [];
  let peak = -1;
  let inPeak = false;
  for (const p of points) {
    if (p.alpha > 0) {
      if (!inPeak) peak++;
      inPeak = true;
      peakOf.push(peak);
    } else {
      inPeak = false;
      peakOf.push(-1);
    }
  }
  // A point at zero: the next peak's color, else the last peak's.
  let next = Math.max(0, peak);
  for (let i = points.length - 1; i >= 0; i--) {
    if (peakOf[i] >= 0) next = peakOf[i];
    else peakOf[i] = next;
  }
  return points.map((p, i) => ({ ...p, color: colors[peakOf[i] % colors.length] }));
}

/** A preset's points as control points: their colors, or the current curve's where they have none. */
function presetPoints(preset: TransferPreset, current: readonly TransferPoint[]): TransferPoint[] {
  return preset.points.map((p) => ({ x: p.x, alpha: p.alpha, color: p.color ?? transferColorAt(current, p.x) }));
}

/** The curve is a preset's when every point's place and opacity is. */
function samePoints(
  a: readonly { x: number; alpha: number }[],
  b: readonly { x: number; alpha: number }[],
  span: number,
) {
  return (
    a.length === b.length &&
    a.every((p, i) => Math.abs(p.x - b[i].x) < span * 1e-3 && Math.abs(p.alpha - b[i].alpha) < 1e-4)
  );
}

const alphaText = (a: number) => (a === 0 ? '0' : a < 0.01 ? a.toPrecision(2) : a.toFixed(3));

/**
 * A transfer function over a histogram of the values: a "Transfer function" header with the presets, the palette
 * button and the σ / abs switch; over the histogram, the opacity as a polyline through its points (x the value, y the
 * opacity, a log axis with `logMin`) with the area under it in the points' colors; under the plot the color bar with a
 * handle per point (a click opens its color); the axis' ends in the unit.
 *
 * Points: drag one (x and y, kept between its neighbours), once a frame (`drag`, the moves a frame state); a press on
 * the plot adds one there (in the color the curve has there) and drags it; a double click removes one, as do Delete
 * and the color popover's Remove (two stay). A point is focusable (a slider: its value in the unit, its opacity): the
 * arrows move it (left and right along the axis, up and down in opacity; Shift: further), Enter opens its color.
 * Model-free: the points come in as data (a state the caller sets from `onChange`), the end of a gesture `onCommit`.
 */
export function TransferFunction(props: TransferFunctionProps) {
  const stats = held<ValueStats | undefined>(props.stats, undefined);
  return <Show when={stats}>{(s) => <TransferCurve {...props} stats={s} />}</Show>;
}

/** The control for one set of statistics. */
function TransferCurve(props: TransferFunctionProps & { stats: ValueStats }) {
  const { stats } = props;
  const axis = valueAxis(stats);
  const unit = unitState(props.unit);
  const points = held<readonly TransferPoint[]>(props.points, []);
  const histogram = held<Histogram | undefined>(props.histogram, undefined);
  const title = props.title ?? 'Transfer function';
  const commit = () => props.onCommit?.();
  const change = (p: readonly TransferPoint[]) => props.onChange(p);

  // The opacity axis.
  const logMin = props.logMin;
  const logLo = logMin !== undefined ? Math.log10(logMin) : 0;
  /** An opacity's height on the plot, 0 to 1. */
  const yOf = (alpha: number) =>
    logMin === undefined
      ? clamp(alpha, 0, 1)
      : alpha <= logMin
        ? 0
        : ZeroBand + (1 - ZeroBand) * clamp((Math.log10(alpha) - logLo) / -logLo, 0, 1);
  /** The opacity at a height on the plot (on a log axis, the bottom band is zero). */
  const alphaAt = (y: number) =>
    logMin === undefined
      ? clamp(y, 0, 1)
      : y < ZeroBand
        ? 0
        : 10 ** (logLo - logLo * clamp((y - ZeroBand) / (1 - ZeroBand), 0, 1));
  const pointText = (p: TransferPoint, u: ValueUnit) => `${valueText(u, p.x, stats)}, opacity ${alphaText(p.alpha)}`;

  /** The point whose color is open. */
  const editing = state<number | undefined>(undefined);
  const palettesOpen = state(false);
  let paletteAnchor!: HTMLElement;
  const handles: HTMLElement[] = [];
  const markers: HTMLElement[] = [];

  /** The positions of the points, told only when their number changes. */
  const indices = derive(
    points,
    (l) => Array.from({ length: l.length }, (_, i) => i),
    (a, b) => a.length === b.length,
  );
  /** The point at `i`: the last one there (one removed keeps its last until its element goes). */
  const pointAt = (i: number): ReadonlyState<TransferPoint> => {
    let last = points.value[i];
    return derive(points, (l) => (last = l[i] ?? last));
  };

  /**
   * A point moved to (x, alpha), kept between its neighbours, in `base` (a gesture's own list: what it last gave out,
   * whether or not the caller has applied it yet) or the points shown. Returns the list given out.
   */
  const moveTo = (i: number, x: number, alpha: number, base: readonly TransferPoint[] = points.value) => {
    const list = [...base];
    if (!list[i]) return base;
    const lo = i > 0 ? list[i - 1].x + Gap * axis.span : stats.min;
    const hi = i < list.length - 1 ? list[i + 1].x - Gap * axis.span : stats.max;
    list[i] = { ...list[i], x: clamp(x, lo, hi), alpha: clamp(alpha, 0, 1) };
    change(list);
    return list;
  };
  const remove = (i: number) => {
    const list = [...points.value];
    // A curve keeps two points at least.
    if (list.length <= 2 || !list[i]) return;
    list.splice(i, 1);
    if (editing.value !== undefined) editing.set(undefined);
    change(list);
    commit();
  };

  // ---- Points dragged on the plot (the shared plot gesture: once a frame, the plot's box read at each step) ----
  // The gesture keeps its own list (what it last gave out): a point it added is dragged even when the caller applies
  // `onChange` later (its own state commit), and never a neighbour that is still at that index in the points shown.
  let working: readonly TransferPoint[] | undefined;
  let lastPress: { index: number; time: number } | undefined;
  const pressedTwice = (index: number) => {
    const now = performance.now();
    const twice = lastPress?.index === index && now - lastPress.time < DoublePressMs;
    lastPress = twice ? undefined : { index, time: now };
    return twice;
  };
  const plotDrag = (el: HTMLElement) =>
    plotGesture<number>(el, {
      start: (e, at) => {
        const hit = (e.target as Element).closest?.<HTMLElement>('[data-point]')?.dataset.point;
        let index: number;
        if (hit !== undefined) {
          index = Number(hit);
          if (pressedTwice(index)) {
            remove(index);
            return undefined;
          }
          working = points.value;
        } else {
          // A new point where it was pressed, in the color the curve has there, then dragged from there.
          const x = axis.valueAt(at.x);
          const list = [...points.value];
          index = list.findIndex((q) => q.x > x);
          if (index < 0) index = list.length;
          list.splice(index, 0, { x, alpha: alphaAt(clamp(at.y, 0, 1)), color: transferColorAt(list, x) });
          change(list);
          working = list;
        }
        markers[index]?.focus({ preventScroll: true });
        return index;
      },
      move: (i, at) => {
        working = moveTo(i, axis.valueAt(at.x), alphaAt(clamp(at.y, 0, 1)), working);
      },
      end: () => {
        working = undefined;
        commit();
      },
    });

  const keys = keySteps(commit);
  const onPointKey = (i: number) => (e: KeyboardEvent) => {
    const p = points.value[i];
    if (!p || e.altKey || e.ctrlKey || e.metaKey) return;
    if (e.key === 'Delete' || e.key === 'Backspace') {
      e.preventDefault();
      // The last point removed: the focus to the one before (the others keep their elements).
      if (i === points.value.length - 1 && points.value.length > 2) markers[i - 1]?.focus();
      remove(i);
      return;
    }
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      editing.set(i);
      return;
    }
    const dx = (e.shiftKey ? 0.02 : 0.005) * axis.span;
    const up = (a: number) =>
      logMin === undefined ? a + (e.shiftKey ? 0.1 : 0.02) : a <= 0 ? logMin : a * (e.shiftKey ? 2 : 1.25);
    const down = (a: number) =>
      logMin === undefined ? a - (e.shiftKey ? 0.1 : 0.02) : a <= logMin * 1.01 ? 0 : a / (e.shiftKey ? 2 : 1.25);
    const next =
      e.key === 'ArrowLeft'
        ? [p.x - dx, p.alpha]
        : e.key === 'ArrowRight'
          ? [p.x + dx, p.alpha]
          : e.key === 'ArrowUp'
            ? [p.x, up(p.alpha)]
            : e.key === 'ArrowDown'
              ? [p.x, down(p.alpha)]
              : undefined;
    if (!next) return;
    e.preventDefault();
    keys.stepped();
    moveTo(i, next[0], next[1]);
  };
  // What a point's keys do (a slider in two directions): said once, for every point.
  const keysHelp = uniqueId('msui-tf-keys');

  // ---- What all the points draw together: the curve, the area, its gradient and the color bar (one watch) ----
  const gradientId = uniqueId('msui-tf');
  const svgNs = 'http://www.w3.org/2000/svg';
  const gradient = (
    <linearGradient id={gradientId} gradientUnits="userSpaceOnUse" x1="0" x2="1" y1="0" y2="0" />
  ) as SVGElement;
  const area = (<path class={classes.area} fill={`url(#${gradientId})`} />) as SVGElement;
  const curve = (<path class={classes.curve} />) as SVGElement;
  const bar = (<div class={classes.bar} />) as HTMLElement;
  watch([points], () => {
    const list = points.value;
    const xy = list.map((q) => `${axis.at(q.x).toFixed(4)} ${(1 - yOf(q.alpha)).toFixed(4)}`);
    setAttr(curve, 'd', xy.length ? `M${xy.join('L')}` : '');
    setAttr(
      area,
      'd',
      xy.length ? `M${axis.at(list[0].x).toFixed(4)} 1L${xy.join('L')}L${axis.at(list.at(-1)!.x).toFixed(4)} 1Z` : '',
    );
    // The gradient's stops: one a point (the colors are the points'); made only when their number changes.
    while (gradient.childNodes.length > list.length) gradient.lastChild!.remove();
    while (gradient.childNodes.length < list.length) gradient.append(document.createElementNS(svgNs, 'stop'));
    list.forEach((q, i) => {
      const stop = gradient.childNodes[i] as SVGElement;
      setAttr(stop, 'offset', axis.at(q.x).toFixed(4));
      setAttr(stop, 'stop-color', colorToHex(q.color));
    });
    setStyle(
      bar,
      'background',
      list.length
        ? `linear-gradient(90deg, ${list.map((q) => `${colorToHex(q.color)} ${(axis.at(q.x) * 100).toFixed(2)}%`).join(', ')})`
        : null,
    );
  });

  const marker = (i: number) => {
    const point = pointAt(i);
    const el = (
      <span
        class={classes.point}
        data-point={String(i)}
        role="slider"
        tabIndex={0}
        aria-label={`Point ${i + 1}`}
        aria-describedby={keysHelp}
        onKeyDown={onPointKey(i)}
        onKeyUp={keys.up}
      />
    ) as HTMLElement;
    markers[i] = el;
    const track = (<span class={classes.track}>{el}</span>) as HTMLElement;
    // The point's own properties: one watch for it, the unit and whether its color is open.
    const open = derive(editing, (e) => e === i);
    watch([point, unit, open], () => {
      const p = point.value;
      const u = unit.value;
      setStyle(track, 'transform', `translate(${axis.at(p.x) * 100}%, ${(1 - yOf(p.alpha)) * 100}%)`);
      setClass(el, classes.on, open.value);
      setAttr(el, 'aria-valuemin', valueNumber(u, stats.min, stats));
      setAttr(el, 'aria-valuemax', valueNumber(u, stats.max, stats));
      setAttr(el, 'aria-valuenow', valueNumber(u, p.x, stats));
      setAttr(el, 'aria-valuetext', pointText(p, u));
    });
    return track;
  };

  const colorHandle = (i: number) => {
    const point = pointAt(i);
    const el = (
      <button
        type="button"
        class={classes.handle}
        aria-haspopup="dialog"
        onClick={() => editing.set(editing.value === i ? undefined : i)}
      />
    ) as HTMLElement;
    handles[i] = el;
    const track = (<span class={classes.handleTrack}>{el}</span>) as HTMLElement;
    const open = derive(editing, (e) => e === i);
    watch([point, unit, open], () => {
      const p = point.value;
      setStyle(track, 'transform', `translateX(${axis.at(p.x) * 100}%)`);
      setStyle(el, '--msui-stop', colorToHex(p.color));
      setAttr(el, 'aria-label', `Color at ${valueText(unit.value, p.x, stats)}`);
      setAttr(el, 'aria-expanded', String(open.value));
    });
    return track;
  };

  const presets = props.presets;
  const palettes = props.palettes;
  return (
    <section class={controlClass} aria-label={title}>
      <PlotHeader title={title}>
        {presets?.length ? (
          <PresetsButton
            label="Curve"
            presets={() => presets.map((p) => [p, p.label] as const)}
            current={derive(points, (list) => {
              const i = presets.findIndex((p) => samePoints(p.points, list, axis.span));
              return i < 0 ? undefined : i;
            })}
            onPick={(p) => {
              change(presetPoints(p, points.value));
              commit();
            }}
          />
        ) : null}
        {palettes?.length ? (
          <span ref={(el: HTMLElement) => (paletteAnchor = el)} style="display: inline-flex;">
            <Button
              variant="ghost"
              size="xs"
              icon={<Icon name="palette" />}
              title="Colors"
              pressed={palettesOpen}
              onClick={() => palettesOpen.set(!palettesOpen.value)}
            />
          </span>
        ) : null}
        <UnitSwitch unit={unit} />
      </PlotHeader>
      <HistogramPlot
        stats={stats}
        axis={axis}
        histogram={histogram}
        height={64}
        class={classes.plot}
        data-tip="Drag a point · press to add one · double-click to remove it"
        ref={plotDrag}
      >
        <svg class={classes.svg} viewBox="0 0 1 1" preserveAspectRatio="none" aria-hidden="true">
          <defs>{gradient}</defs>
          {area}
          {curve}
        </svg>
        <For each={indices}>{(i) => marker(i)}</For>
      </HistogramPlot>
      <div class={classes.colors}>
        {bar}
        <For each={indices}>{(i) => colorHandle(i)}</For>
      </div>
      <AxisLabels stats={stats} unit={unit} />
      <span id={keysHelp} class={classes.help}>
        Left and right move the point along the values, up and down change its opacity; Shift moves further; Enter opens
        its color; Delete removes it.
      </span>
      <Show when={derive(editing, (i) => (i === undefined ? undefined : { i }))}>
        {({ i }) => (
          <PointColor
            anchor={handles[i]}
            point={pointAt(i)}
            label={derive([pointAt(i), unit], (p, u) => `Color at ${valueText(u, p.x, stats)}`)}
            canRemove={derive(points, (l) => l.length > 2)}
            onColor={(color) => {
              const list = [...points.value];
              if (!list[i]) return;
              list[i] = { ...list[i], color };
              change(list);
            }}
            onRemove={() => {
              // The focus to a neighbour's handle: the one before; for the first, its own (handles are by place: it
              // shows the next point once this one goes).
              (handles[i - 1] ?? handles[i])?.focus();
              remove(i);
            }}
            onClose={(refocus) => {
              editing.set(undefined);
              commit();
              if (refocus) handles[i]?.focus();
            }}
          />
        )}
      </Show>
      {palettes?.length ? (
        <Show when={palettesOpen}>
          {() => (
            <Popover
              anchor={paletteAnchor}
              title="Colors"
              onClose={(refocus) => {
                palettesOpen.set(false);
                if (refocus) (paletteAnchor.firstElementChild as HTMLElement)?.focus();
              }}
            >
              <div class={classes.palettes}>
                {palettes.map((p) => (
                  <button
                    type="button"
                    class={classes.palette}
                    onClick={() => {
                      palettesOpen.set(false);
                      change(applyPalette(points.value, p, axis));
                      commit();
                    }}
                  >
                    <span class={classes.swatches} aria-hidden="true">
                      {p.colors.map((c) => (
                        <span class={classes.swatch} style={{ background: colorToHex(c) }} />
                      ))}
                    </span>
                    {p.label}
                  </button>
                ))}
              </div>
            </Popover>
          )}
        </Show>
      ) : null}
    </section>
  );
}

/** A point's color in a popover from its handle: the picker, and removing the point. */
function PointColor(props: {
  anchor: HTMLElement;
  point: ReadonlyState<TransferPoint>;
  label: ReadonlyState<string>;
  canRemove: ReadonlyState<boolean>;
  onColor(color: number): void;
  onRemove(): void;
  onClose(refocus?: boolean): void;
}) {
  return (
    <Popover anchor={props.anchor} align="start" title={props.label} onClose={props.onClose}>
      <div class={classes.panel}>
        <ColorPicker value={derive(props.point, (p) => p.color)} onChange={props.onColor} />
        <div class={classes.foot}>
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="trash-2" />}
            disabled={derive(props.canRemove, (c) => !c)}
            onClick={props.onRemove}
          >
            Remove point
          </Button>
        </div>
      </div>
    </Popover>
  );
}
