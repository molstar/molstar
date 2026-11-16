/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, type State } from '@molstar/core/reactive';
import { clamp } from '@molstar/core/math/interpolate';
import { type Child, setAttr, watch } from '../../dom/index.js';
import { state } from '../../state.js';
import { classes } from './_generated/histogram.module.css.js';

// What the controls over a map's values share (the iso level, the transfer function): the values' statistics and
// histogram as data, the plot's value ↔ x mapping, the σ / map-units switch and how a value is written in either, the
// plot's well with the histogram and the mean, and the axis' ends under it.

/** The values' statistics: what σ is measured against and the axis spans. */
export interface ValueStats {
  readonly min: number;
  readonly max: number;
  readonly mean: number;
  /** The standard deviation (0: σ is not defined, steps fall back to a part of the range). */
  readonly sigma: number;
}

/** Counts per bin over [min, max], equal bins. Its range need not be the axis': the bins are placed by value. */
export interface Histogram {
  readonly min: number;
  readonly max: number;
  readonly counts: ArrayLike<number>;
}

/** How values are shown and typed: σ from the mean, or the values' own units. */
export type ValueUnit = 'sigma' | 'absolute';

/** The plot's axis: values in [min, max] to x in [0, 1] and back. */
export interface ValueAxis {
  readonly min: number;
  readonly max: number;
  /** max − min (1 for an empty range, so nothing divides by zero). */
  readonly span: number;
  /** A value's place, 0 to 1 (clamped). */
  at(value: number): number;
  /** The value at a place (0 to 1, clamped). */
  valueAt(x: number): number;
  /** A value kept in [min, max]. */
  clamp(value: number): number;
}

export function valueAxis(stats: ValueStats): ValueAxis {
  const { min, max } = stats;
  const span = max - min || 1;
  return {
    min,
    max,
    span,
    at: (v) => clamp((v - min) / span, 0, 1),
    valueAt: (x) => min + span * clamp(x, 0, 1),
    clamp: (v) => clamp(v, min, max),
  };
}

/** A value in σ from the mean (σ 0: counted in units). */
export function sigmaOf(value: number, stats: ValueStats): number {
  return (value - stats.mean) / (stats.sigma || 1);
}

/** A value from σ. */
export function fromSigma(sigmas: number, stats: ValueStats): number {
  return stats.mean + sigmas * (stats.sigma || 1);
}

const minus = (s: string) => s.replace('-', '−');

/** A value as the number written in a unit: "1.50" (σ, two decimals), "0.0734" (three significant digits). */
export function formatValue(unit: ValueUnit, value: number, stats: ValueStats): string {
  if (unit === 'sigma') return minus(sigmaOf(value, stats).toFixed(2));
  return minus(String(Number(value.toPrecision(3))));
}

/** A value with its unit: "1.50 σ", "0.0734". */
export function valueText(unit: ValueUnit, value: number, stats: ValueStats): string {
  return unit === 'sigma' ? `${formatValue(unit, value, stats)} σ` : formatValue(unit, value, stats);
}

/** A value in a unit as the plain number ARIA wants (an ASCII minus; the typographic one is for text). */
export function valueNumber(unit: ValueUnit, value: number, stats: ValueStats): number {
  return unit === 'sigma' ? Number(sigmaOf(value, stats).toFixed(2)) : Number(value.toPrecision(3));
}

/**
 * The histogram as a filled step path in a 1 × 1 box over the axis (each bin placed by its values, the parts outside
 * the axis cut off); log counts, as a map's tails are thin.
 */
export function histogramPath(h: Histogram | undefined, axis: ValueAxis): string {
  if (!h || !h.counts.length) return '';
  const n = h.counts.length;
  let top = 0;
  for (let i = 0; i < n; i++) top = Math.max(top, h.counts[i]);
  const scale = top > 0 ? 1 / Math.log1p(top) : 0;
  const width = (h.max - h.min) / n;
  const parts = [`M${axis.at(h.min).toFixed(4)} 1`];
  for (let i = 0; i < n; i++) {
    const y = (1 - Math.log1p(h.counts[i]) * scale).toFixed(4);
    const x0 = axis.at(h.min + i * width).toFixed(4);
    const x1 = axis.at(h.min + (i + 1) * width).toFixed(4);
    parts.push(`L${x0} ${y}L${x1} ${y}`);
  }
  parts.push(`L${axis.at(h.max).toFixed(4)} 1Z`);
  return parts.join('');
}

/** The unit shown: the caller's (shared between controls), or the control's own, σ. */
export function unitState(unit: State<ValueUnit> | undefined): State<ValueUnit> {
  return unit ?? state<ValueUnit>('sigma');
}

const UnitLabels: Record<ValueUnit, readonly [label: string, title: string]> = {
  sigma: ['σ', 'σ from the mean'],
  absolute: ['abs', 'Absolute values'],
};

/**
 * σ or the values' units: a radio group of two xs controls (20px), the one in use pressed (WAI-ARIA radio group: one
 * Tab stop, the arrows pick the other).
 */
export function UnitSwitch(props: { unit: State<ValueUnit> }) {
  const { unit } = props;
  const units = ['sigma', 'absolute'] as const;
  const buttons = units.map(
    (u) =>
      (
        <button type="button" role="radio" class={classes.unit} data-tip={UnitLabels[u][1]} onClick={() => unit.set(u)}>
          {UnitLabels[u][0]}
        </button>
      ) as HTMLButtonElement,
  );
  watch([unit], () => {
    units.forEach((u, i) => {
      const on = unit.value === u;
      setAttr(buttons[i], 'aria-checked', String(on));
      buttons[i].tabIndex = on ? 0 : -1;
    });
  });
  const onKey = (e: KeyboardEvent) => {
    if (!['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown', 'Home', 'End'].includes(e.key)) return;
    e.preventDefault();
    const i =
      e.key === 'Home' ? 0 : e.key === 'End' ? 1 : (units.indexOf(unit.value as (typeof units)[number]) + 1) % 2;
    unit.set(units[i]);
    buttons[i].focus();
  };
  return (
    <span class={classes.units} role="radiogroup" aria-label="Unit" onKeyDown={onKey}>
      {buttons}
    </span>
  );
}

/** A control's header: its title on the inset, its tools at the plot's right edge. */
export function PlotHeader(props: { title: Child; children?: Child }) {
  return (
    <div class={classes.head}>
      <span class={classes.title}>{props.title}</span>
      <span class={classes.tools}>{props.children}</span>
    </div>
  );
}

/**
 * The values' well: the histogram (a muted area) and the mean (a dashed line), then `children` over them (each placed
 * absolutely over the plot). `height` in px (on the 4px grid; default 48).
 */
export function HistogramPlot(props: {
  stats: ValueStats;
  axis: ValueAxis;
  histogram: ReadonlyState<Histogram | undefined>;
  height?: number;
  class?: string;
  /** The plot's element (a drag is attached to it). */
  ref?: (el: HTMLElement) => void;
  'data-tip'?: string;
  children?: Child;
}) {
  const { axis } = props;
  return (
    <div
      class={props.class ? `${classes.plot} ${props.class}` : classes.plot}
      style={props.height ? { '--msui-plot-height': `${props.height}px` } : undefined}
      data-tip={props['data-tip']}
      ref={props.ref}
    >
      <svg class={classes.layer} viewBox="0 0 1 1" preserveAspectRatio="none" aria-hidden="true">
        <path class={classes.bars} d={derive(props.histogram, (h) => histogramPath(h, axis))} />
      </svg>
      <span
        class={classes.mean}
        style={{ transform: `translateX(${axis.at(props.stats.mean) * 100}%)` }}
        aria-hidden="true"
      />
      {props.children}
    </div>
  );
}

/** The axis' ends under the plot, in the unit. */
export function AxisLabels(props: { stats: ValueStats; unit: ReadonlyState<ValueUnit> }) {
  const { stats, unit } = props;
  return (
    <div class={classes.axis} aria-hidden="true">
      <span>{derive(unit, (u) => valueText(u, stats.min, stats))}</span>
      <span>{derive(unit, (u) => valueText(u, stats.max, stats))}</span>
    </div>
  );
}

/** The control's root class (it pads the rest of the inset in a surface, styles §5). */
export const controlClass: string = classes.control;
