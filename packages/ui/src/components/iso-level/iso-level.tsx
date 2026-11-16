/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, type State } from '@molstar/core/reactive';
import {
  drag,
  effect,
  For,
  frameState,
  type MaybeReactive,
  setAttr,
  setClass,
  setStyle,
  watch,
} from '../../dom/index.js';
import { held } from '../../state.js';
import { colorToHex } from '../color/color-picker.js';
import {
  AxisLabels,
  controlClass,
  type Histogram,
  HistogramPlot,
  PlotHeader,
  UnitSwitch,
  unitState,
  type ValueStats,
  type ValueUnit,
  valueAxis,
  valueNumber,
  valueText,
} from '../histogram/histogram.js';
import { classes } from './_generated/iso-level.module.css.js';

/** A level: an iso value on the values' axis, in a color. */
export interface IsoLevel {
  /** Its key (stable while it lives: a marker keeps its element). */
  readonly id: string;
  /** Its name, for its slider's accessible name ("<label> iso value"). */
  readonly label: string;
  /** The iso value, in the values' own units. */
  readonly value: number;
  /** Its marker's color; none: the text's. */
  readonly color?: number;
  /** Not shown: its marker and what it encloses faded. */
  readonly hidden?: boolean;
}

export interface IsoLevelProps {
  /** The header's title (default "Levels"); also names the section. */
  title?: string;
  /** The values' statistics: the axis spans [min, max]; σ is counted from the mean. */
  stats: ValueStats;
  /** Counts of the values; none (yet): an empty well. */
  histogram?: MaybeReactive<Histogram | undefined>;
  /** The levels, with any change made through `onChange` already in (the control shows what it is given). */
  levels: MaybeReactive<readonly IsoLevel[]>;
  /** The unit values are written in: shared between controls when given; else the control's own (σ). */
  unit?: State<ValueUnit>;
  /** A level moved (pointer rate: once a frame while dragged; a key). The value is within [min, max]. */
  onChange(id: string, value: number): void;
  /** The end of a gesture (a drag let go, a key released): the changes since are one step. */
  onCommit?(): void;
}

/** Arrow-key steps in σ (Shift, PageUp/PageDown: the big one); without σ, a hundredth and a tenth of the range. */
const Step = 0.1;
const BigStep = 1;

/**
 * Iso levels over a histogram of the values: a "Levels" header with the σ / abs switch; the histogram in a well, the
 * mean dashed, a marker per level (a line in its color, a handle on top) with what the level encloses tinted (above it;
 * below one under the mean); the axis' ends under it in the unit. Drag anywhere on the plot: the marker under the
 * pointer, else the nearest, follows it, once a frame (`drag`, the moves a frame state). A marker is a slider (WAI-ARIA):
 * the arrows step it by 0.1 σ (Shift: 1 σ), PageUp/PageDown by 1 σ, Home and End to the ends; its value is read in the
 * unit. Model-free: the levels come in as data, a move goes out through `onChange`, the end of a gesture `onCommit`.
 */
export function IsoLevelControl(props: IsoLevelProps) {
  const { stats } = props;
  const axis = valueAxis(stats);
  const unit = unitState(props.unit);
  const levels = held<readonly IsoLevel[]>(props.levels, []);
  const histogram = held<Histogram | undefined>(props.histogram, undefined);
  const title = props.title ?? 'Levels';
  const ids = derive(
    levels,
    (ls) => ls.map((l) => l.id),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  /** A level by id, as a state: the last one found (one removed keeps its last until its marker goes). */
  const levelOf = (id: string): ReadonlyState<IsoLevel> => {
    let last = levels.value.find((l) => l.id === id) as IsoLevel;
    return derive(levels, (ls) => (last = ls.find((l) => l.id === id) ?? last));
  };
  const commit = () => props.onCommit?.();

  // Dragging: the pointer's x is input, published once a frame (a frame state; its value only a tick), and the level
  // follows it then. The drag's listeners are attached once (`drag`).
  let rect: DOMRect | undefined;
  let dragging: string | undefined;
  let pointerX = 0;
  let moves = 0;
  const moved = frameState(0);
  const step = () => {
    if (dragging === undefined || !rect) return;
    props.onChange(dragging, axis.valueAt((pointerX - rect.left) / (rect.width || 1)));
  };
  effect(moved, step);
  const plotDrag = (el: HTMLElement) =>
    drag(el, {
      start: (e) => {
        if (e.button !== 0) return undefined;
        rect = el.getBoundingClientRect();
        const hit = (e.target as Element).closest?.<HTMLElement>('[data-level]')?.dataset.level;
        const x = (e.clientX - rect.left) / (rect.width || 1);
        let nearest: IsoLevel | undefined;
        for (const l of levels.value) {
          if (!nearest || Math.abs(axis.at(l.value) - x) < Math.abs(axis.at(nearest.value) - x)) nearest = l;
        }
        dragging = hit ?? nearest?.id;
        if (dragging === undefined) return undefined;
        e.preventDefault();
        for (const m of el.querySelectorAll<HTMLElement>('[data-level]')) {
          if (m.dataset.level === dragging) m.focus({ preventScroll: true });
        }
        pointerX = e.clientX;
        step();
        return dragging;
      },
      move: (e) => {
        pointerX = e.clientX;
        moved.offer(++moves);
      },
      // Let go: where it was let go, now. Ended otherwise (cancelled, the scope paused or gone): no more steps. Either way the
      // gesture is over.
      end: (_id, how, e) => {
        if (how === 'up' && e) {
          pointerX = e.clientX;
          step();
        }
        dragging = undefined;
        rect = undefined;
        commit();
      },
    });

  let keyed = false;
  const onKey = (id: string) => (e: KeyboardEvent) => {
    const level = levels.value.find((l) => l.id === id);
    if (!level || e.altKey || e.ctrlKey || e.metaKey) return;
    const small = (e.shiftKey ? BigStep : Step) * (stats.sigma || axis.span / 100);
    const big = BigStep * (stats.sigma || axis.span / 10);
    const to =
      e.key === 'ArrowRight' || e.key === 'ArrowUp'
        ? level.value + small
        : e.key === 'ArrowLeft' || e.key === 'ArrowDown'
          ? level.value - small
          : e.key === 'PageUp'
            ? level.value + big
            : e.key === 'PageDown'
              ? level.value - big
              : e.key === 'Home'
                ? stats.min
                : e.key === 'End'
                  ? stats.max
                  : undefined;
    if (to === undefined) return;
    e.preventDefault();
    keyed = true;
    props.onChange(id, axis.clamp(to));
  };
  // A key held or pressed again before it is released is one gesture: it ends when the key is up.
  const onKeyUp = () => {
    if (!keyed) return;
    keyed = false;
    commit();
  };

  const marker = (id: string) => {
    const level = levelOf(id);
    const inside = (<span class={classes.inside} aria-hidden="true" />) as HTMLElement;
    const handle = (
      <span class={classes.marker} data-level={id} role="slider" tabIndex={0} onKeyDown={onKey(id)} onKeyUp={onKeyUp} />
    ) as HTMLElement;
    const track = (<span class={classes.track}>{handle}</span>) as HTMLElement;
    // The marker's own properties: one watch for the level and the unit.
    watch([level, unit], () => {
      const l = level.value;
      const u = unit.value;
      const x = axis.at(l.value);
      const color = l.color === undefined ? null : colorToHex(l.color);
      setStyle(track, 'transform', `translateX(${x * 100}%)`);
      // What the surface encloses: above the level, or below one under the mean (a difference map's negative level).
      setStyle(inside, 'transform', `translateX(${l.value >= stats.mean ? x * 100 : x * 100 - 100}%)`);
      setStyle(track, '--msui-level', color);
      setStyle(inside, '--msui-level', color);
      setClass(handle, classes.off, !!l.hidden);
      setClass(inside, classes.off, !!l.hidden);
      setAttr(handle, 'aria-label', `${l.label} iso value`);
      setAttr(handle, 'aria-valuemin', valueNumber(u, stats.min, stats));
      setAttr(handle, 'aria-valuemax', valueNumber(u, stats.max, stats));
      setAttr(handle, 'aria-valuenow', valueNumber(u, l.value, stats));
      setAttr(handle, 'aria-valuetext', valueText(u, l.value, stats));
    });
    return (
      <>
        {inside}
        {track}
      </>
    );
  };

  return (
    <section class={controlClass} aria-label={title}>
      <PlotHeader title={title}>
        <UnitSwitch unit={unit} />
      </PlotHeader>
      <HistogramPlot
        stats={stats}
        axis={axis}
        histogram={histogram}
        class={classes.plot}
        data-tip="Drag to move a level"
        ref={plotDrag}
      >
        <For each={ids}>{(id) => marker(id)}</For>
      </HistogramPlot>
      <AxisLabels stats={stats} unit={unit} />
    </section>
  );
}
