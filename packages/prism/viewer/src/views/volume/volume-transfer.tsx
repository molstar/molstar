/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, ColorPicker, colorToHex, Floating, PresetsButton, PresetsList } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { combineLatest, distinctUntilChanged, filter, map, type Observable, shareReplay } from 'rxjs';
import type { VolumeLevelUnit } from '../../model/volume-levels.js';
import {
  type TransferColors,
  type TransferPoint,
  type TransferStop,
  transferColorAt,
  transferPeaks,
  type VolumeTransferModel,
} from '../../model/volume-transfer.js';
import { classes } from './_generated/volume-transfer.module.css.js';
import { formatIn, histogramPath } from './volume-levels.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// A direct volume's transfer function (model/volume-transfer.ts): on one axis, the map's values as a histogram, the
// opacity as a curve through its points (a log scale: opacities per step are small) with the area under it in the
// colors the renderer paints there, and under it the color strip with the theme's stops. Points: drag, click to add,
// double-click to remove, arrow keys and Delete when focused. Stops: the same on the strip; a picked stop opens its
// color; dragged near a point, a stop snaps to it, and a stop on a point moves with it. Colors at peaks (a popover of
// palettes, or the theme's own colors) puts stops at the ends of each of the curve's peaks; a preset moves them with it.

/** Samples of the colors for the strip and the area's gradient. */
const Samples = 64;
/** A stop dragged this close to a point (px) snaps to it. */
const SnapPx = 6;
/** A stop this close to a point (on the axis, 0 to 1) is on it, and moves with it. */
const LinkEps = 0.0015;
/** The bottom of the plot (this fraction of its height) is opacity zero. */
const ZeroBand = 0.04;

type Picked = { kind: 'point' | 'stop'; index: number } | { kind: 'uniform' } | undefined;

const clamp01 = (x: number) => Math.min(1, Math.max(0, x));
/** Two presses on a point or a stop within this long are a double click. */
const DoublePressMs = 350;

/** A curve is a preset's when every point is (as stored: x to three decimals). */
function samePoints(a: readonly TransferPoint[], b: readonly TransferPoint[]) {
  return (
    a.length === b.length && a.every((p, i) => Math.abs(p.x - b[i].x) < 1e-3 && Math.abs(p.alpha - b[i].alpha) < 1e-4)
  );
}

/** The colors along the axis, as [x, css color] samples; none: not on the axis (another theme). */
function samplesOf(colors: TransferColors): [number, string][] | undefined {
  if (colors.kind === 'theme') return undefined;
  const out: [number, string][] = [];
  for (let i = 0; i <= Samples; i++) {
    const x = i / Samples;
    const c = transferColorAt(colors, x);
    out.push([x, c === undefined ? 'transparent' : colorToHex(c)]);
  }
  return out;
}

export function VolumeTransfer(props: { model: VolumeTransferModel }) {
  const { model } = props;
  const { stats } = model;
  const span = stats.max - stats.min || 1;
  const logLo = Math.log10(model.logMin);
  /** An opacity's height on the plot, 0 to 1. */
  const yOf = (alpha: number) =>
    alpha <= model.logMin ? 0 : ZeroBand + (1 - ZeroBand) * clamp01((Math.log10(alpha) - logLo) / -logLo);
  /** The opacity at a height on the plot (the bottom band is zero). */
  const alphaAt = (y: number) => (y < ZeroBand ? 0 : 10 ** (logLo - logLo * clamp01((y - ZeroBand) / (1 - ZeroBand))));
  const valueText = (x: number, u: VolumeLevelUnit) =>
    `${formatIn(u, stats.min + x * span, stats)}${u === 'sigma' ? ' σ' : ''}`;
  const alphaText = (a: number) => (a === 0 ? '0' : a < 0.01 ? a.toPrecision(2) : a.toFixed(3));

  const picked = observableState<Picked>(undefined);
  const presetsOpen = observableState(false);
  const palettesOpen = observableState(false);
  let paletteAnchor: HTMLElement | undefined;
  const gradientId = uniqueId('msp-tf');
  // Computed once per change of the colors for the many bindings below (the gradient's stops, the strip, the classes).
  const samples = model.colors.pipe(map(samplesOf), shareReplay({ bufferSize: 1, refCount: true }));
  const points = model.points;
  const stops = derive(model.colors, (c) => (c.kind === 'value' ? c.stops : []));
  const indices = (o: Observable<readonly unknown[]>) =>
    o.pipe(
      map((l) => l.length),
      distinctUntilChanged(),
      map((n) => Array.from({ length: n }, (_, i) => i)),
    );
  const pointAt = (i: number) =>
    points.pipe(
      map((p) => p[i]),
      filter((p): p is TransferPoint => !!p),
    );
  const stopAt = (i: number) =>
    stops.pipe(
      map((s) => s[i]),
      filter((s): s is TransferStop => !!s),
    );
  const curve = derive(points, (p) =>
    p.length ? p.map((q, i) => `${i ? 'L' : 'M'}${q.x.toFixed(4)} ${(1 - yOf(q.alpha)).toFixed(4)}`).join('') : '',
  );
  const area = derive(points, (p) =>
    p.length
      ? `M${p[0].x.toFixed(4)} 1${p.map((q) => `L${q.x.toFixed(4)} ${(1 - yOf(q.alpha)).toFixed(4)}`).join('')}L${p[p.length - 1].x.toFixed(4)} 1Z`
      : '',
  );

  // ---- Points: dragged on the plot, once a frame ----
  let rect: DOMRect | undefined;
  let dragging: { kind: 'point' | 'stop'; index: number } | undefined;
  let frame = 0;
  let pointer = { x: 0, y: 0 };
  const at = (clientX: number, clientY: number) => ({
    x: clamp01((clientX - rect!.left) / (rect!.width || 1)),
    y: clamp01(1 - (clientY - rect!.top) / (rect!.height || 1)),
  });
  const movePoint = (i: number, x: number, y: number) => {
    const p = [...points.value];
    if (!p[i]) return;
    // Between its neighbours: the curve keeps its order.
    const lo = i > 0 ? p[i - 1].x + 0.001 : 0;
    const hi = i < p.length - 1 ? p[i + 1].x - 0.001 : 1;
    const from = p[i].x;
    p[i] = { x: Math.min(hi, Math.max(lo, x)), alpha: alphaAt(y) };
    model.change({ points: p, ...linkedStops(from, p[i].x) });
  };
  /** Stops on a point move with it (Colors at peaks puts them there; a stop dragged near a point snaps onto it). */
  const linkedStops = (from: number, to: number): { stops?: TransferStop[] } => {
    const colors = model.colors.value;
    if (colors.kind !== 'value' || from === to) return {};
    const s = [...colors.stops];
    let moved = false;
    for (let k = 0; k < s.length; k++) {
      if (Math.abs(s[k].x - from) > LinkEps) continue;
      const lo = Math.max(colors.domain[0], k > 0 ? s[k - 1].x : colors.domain[0]);
      const hi = Math.min(colors.domain[1], k < s.length - 1 ? s[k + 1].x : colors.domain[1]);
      s[k] = { ...s[k], x: Math.min(hi, Math.max(lo, to)) };
      moved = true;
    }
    return moved ? { stops: s } : {};
  };
  const moveStop = (i: number, x: number) => {
    const colors = model.colors.value;
    if (colors.kind !== 'value') return;
    const s = [...colors.stops];
    if (!s[i]) return;
    // Within the theme's domain (outside it the renderer clamps) and between its neighbours; near a point, on it.
    const [d0, d1] = colors.domain;
    const lo = Math.max(d0, i > 0 ? s[i - 1].x + 0.001 : d0);
    const hi = Math.min(d1, i < s.length - 1 ? s[i + 1].x - 0.001 : d1);
    let to = Math.min(hi, Math.max(lo, x));
    const snap = points.value.find((p) => Math.abs(p.x - to) * (rect?.width ?? 0) < SnapPx);
    if (snap && snap.x >= lo && snap.x <= hi) to = snap.x;
    s[i] = { ...s[i], x: to };
    model.change({ stops: s });
  };
  const flush = () => {
    if (frame) cancelAnimationFrame(frame);
    frame = 0;
    if (!dragging || !rect) return;
    const p = at(pointer.x, pointer.y);
    if (dragging.kind === 'point') movePoint(dragging.index, p.x, p.y);
    else moveStop(dragging.index, p.x);
  };
  const end = () => {
    if (!dragging) return;
    flush();
    dragging = undefined;
    rect = undefined;
    model.settle();
  };
  // Gone mid-drag (another row is picked, the node removed): no more frames, and the step closes.
  onCleanup(() => {
    if (frame) cancelAnimationFrame(frame);
    frame = 0;
    if (dragging) {
      dragging = undefined;
      rect = undefined;
      model.settle();
    }
  });
  const onMove = (e: PointerEvent) => {
    if (!dragging) return;
    pointer = { x: e.clientX, y: e.clientY };
    if (!frame) frame = requestAnimationFrame(flush);
  };
  const capture = (e: PointerEvent) => {
    e.preventDefault();
    (e.currentTarget as HTMLElement).setPointerCapture(e.pointerId);
    pointer = { x: e.clientX, y: e.clientY };
  };

  const onPlotDown = (e: PointerEvent) => {
    if (e.button !== 0) return;
    rect = (e.currentTarget as HTMLElement).getBoundingClientRect();
    const hit = (e.target as HTMLElement).closest<HTMLElement>('[data-point]')?.dataset.point;
    let index: number;
    if (hit !== undefined) {
      index = Number(hit);
      if (pressedTwice('point', index)) {
        removePoint(index);
        return;
      }
    } else {
      // A new point where it was clicked, then dragged from there.
      const p = at(e.clientX, e.clientY);
      const list = [...points.value];
      index = list.findIndex((q) => q.x > p.x);
      if (index < 0) index = list.length;
      list.splice(index, 0, { x: p.x, alpha: alphaAt(p.y) });
      model.change({ points: list });
    }
    picked.set({ kind: 'point', index });
    dragging = { kind: 'point', index };
    capture(e);
  };
  // A double click removes a point or a stop: told on the second press (the pointer is captured on the first, so the
  // dblclick event lands on the plot, not on what was pressed).
  let lastPress: { kind: 'point' | 'stop'; index: number; time: number } | undefined;
  const pressedTwice = (kind: 'point' | 'stop', index: number) => {
    const now = performance.now();
    const twice = lastPress?.kind === kind && lastPress.index === index && now - lastPress.time < DoublePressMs;
    lastPress = twice ? undefined : { kind, index, time: now };
    return twice;
  };
  const removePoint = (i: number) => {
    const list = [...points.value];
    // A curve keeps two points at least.
    if (list.length <= 2 || !list[i]) return;
    list.splice(i, 1);
    model.change({ points: list });
    model.settle();
    picked.set(undefined);
  };
  const onPointKey = (i: number) => (e: KeyboardEvent) => {
    const p = points.value[i];
    if (!p) return;
    const dx = e.shiftKey ? 0.02 : 0.005;
    const up = (a: number) => (a <= 0 ? model.logMin : Math.min(1, a * (e.shiftKey ? 2 : 1.25)));
    const down = (a: number) => (a <= model.logMin * 1.01 ? 0 : a / (e.shiftKey ? 2 : 1.25));
    const next =
      e.key === 'ArrowLeft'
        ? { x: p.x - dx, alpha: p.alpha }
        : e.key === 'ArrowRight'
          ? { x: p.x + dx, alpha: p.alpha }
          : e.key === 'ArrowUp'
            ? { x: p.x, alpha: up(p.alpha) }
            : e.key === 'ArrowDown'
              ? { x: p.x, alpha: down(p.alpha) }
              : undefined;
    if (e.key === 'Delete' || e.key === 'Backspace') {
      e.preventDefault();
      removePoint(i);
      return;
    }
    if (!next) return;
    e.preventDefault();
    const list = [...points.value];
    const lo = i > 0 ? list[i - 1].x + 0.001 : 0;
    const hi = i < list.length - 1 ? list[i + 1].x - 0.001 : 1;
    list[i] = { x: Math.min(hi, Math.max(lo, next.x)), alpha: next.alpha };
    model.change({ points: list, ...linkedStops(p.x, list[i].x) });
  };

  // ---- Stops: dragged on the strip ----
  const onStripDown = (e: PointerEvent) => {
    if (e.button !== 0) return;
    const colors = model.colors.value;
    if (colors.kind === 'uniform') {
      picked.set({ kind: 'uniform' });
      return;
    }
    if (colors.kind !== 'value') return;
    rect = (e.currentTarget as HTMLElement).getBoundingClientRect();
    const hit = (e.target as HTMLElement).closest<HTMLElement>('[data-stop]')?.dataset.stop;
    let index: number;
    if (hit !== undefined) {
      index = Number(hit);
      if (pressedTwice('stop', index)) {
        removeStop(index);
        return;
      }
    } else {
      const x = Math.min(colors.domain[1], Math.max(colors.domain[0], at(e.clientX, e.clientY).x));
      const list = [...colors.stops];
      index = list.findIndex((s) => s.x > x);
      if (index < 0) index = list.length;
      list.splice(index, 0, { x, color: transferColorAt(colors, x) ?? 0x808080 });
      model.change({ stops: list });
    }
    picked.set({ kind: 'stop', index });
    dragging = { kind: 'stop', index };
    capture(e);
  };
  const removeStop = (i: number) => {
    const colors = model.colors.value;
    if (colors.kind !== 'value' || colors.stops.length <= 1) return;
    const list = [...colors.stops];
    list.splice(i, 1);
    model.change({ stops: list });
    model.settle();
    picked.set(undefined);
  };
  const setStopColor = (i: number, color: number) => {
    const colors = model.colors.value;
    if (colors.kind !== 'value' || !colors.stops[i]) return;
    const list = [...colors.stops];
    // A peak's color: the stops of the same color within the peak the stop is on (Colors at peaks puts one at each end).
    const at = list[i];
    const peak = transferPeaks(points.value).find((p) => at.x >= p.from - LinkEps && at.x <= p.to + LinkEps);
    for (let k = 0; k < list.length; k++) {
      const same =
        k === i ||
        (peak && list[k].color === at.color && list[k].x >= peak.from - LinkEps && list[k].x <= peak.to + LinkEps);
      if (same) list[k] = { ...list[k], color };
    }
    model.change({ stops: list });
  };

  const strip = samples.pipe(
    map((s) => (s ? `linear-gradient(90deg, ${s.map(([x, c]) => `${c} ${(x * 100).toFixed(2)}%`).join(', ')})` : '')),
  );
  const pickedPoint = derive([picked.asReadonly, points], (p, list) =>
    p?.kind === 'point' ? list[p.index] : undefined,
  );
  const pickedStop = combineLatest([picked.asReadonly, stops]).pipe(
    map(([p, list]) => (p?.kind === 'stop' ? list[p.index] : undefined)),
  );
  const hasPeaks = derive(points, (p) => transferPeaks(p).length > 0);
  /** A picked stop's (or the one) color; its picker opens from the chevron. */
  const colorOpen = observableState(false);
  const pickedColor = combineLatest([pickedStop, model.colors]).pipe(
    map(([s, c]) => s?.color ?? (c.kind === 'uniform' ? c.color : 0x808080)),
  );

  return (
    <section
      class={classes.transfer}
      aria-label={model.title ? `${model.title} transfer function` : 'Transfer function'}
    >
      <header class={classes.head}>
        <span class={classes.title}>{model.title ?? 'Transfer function'}</span>
        <span class={classes.tools}>
          <PresetsButton open={presetsOpen} label="Curve" />
          <button
            type="button"
            class={classes.iconTool}
            disabled={hasPeaks.pipe(map((h) => !h))}
            aria-label="Colors at peaks"
            aria-haspopup="listbox"
            aria-expanded={palettesOpen}
            title="Colors at peaks: a color for each of the curve's peaks"
            ref={(el: HTMLElement) => (paletteAnchor = el)}
            onClick={() => palettesOpen.set(!palettesOpen.value)}
          >
            <Icon name="palette" size={13} />
          </button>
          <span class={classes.units} role="radiogroup" aria-label="Unit">
            {(['sigma', 'absolute'] as const).map((u) => (
              <button
                type="button"
                role="radio"
                class={classes.unit}
                aria-checked={derive(model.unit, (x) => x === u)}
                title={u === 'sigma' ? 'σ above the mean' : 'Map units'}
                onClick={() => model.setUnit(u)}
              >
                {u === 'sigma' ? 'σ' : 'abs'}
              </button>
            ))}
          </span>
        </span>
      </header>
      {/* Mol*'s curves from the map's statistics, as a form's presets show (the one in effect marked). */}
      <PresetsList
        open={presetsOpen}
        presets={() => model.presets().map((p) => [p.points, p.label] as [readonly TransferPoint[], string])}
        current={derive(points, (p) => {
          const i = model.presets().findIndex((q) => samePoints(q.points, p));
          return i < 0 ? undefined : i;
        })}
        // The color stops move onto the new curve's peaks with their colors (one step).
        onPick={(p) => model.applyPreset(p)}
      />
      <Show when={palettesOpen}>
        {() =>
          paletteAnchor ? (
            <Floating anchor={paletteAnchor}>
              <PalettePicker model={model} onClose={() => palettesOpen.set(false)} />
            </Floating>
          ) : null
        }
      </Show>
      <div
        class={classes.plot}
        title="Drag a point · click to add one · double-click to remove it"
        onPointerDown={onPlotDown}
        onPointerMove={onMove}
        onPointerUp={end}
        onPointerCancel={end}
      >
        <svg class={classes.svg} viewBox="0 0 1 1" preserveAspectRatio="none" aria-hidden="true">
          <defs>
            <linearGradient id={gradientId} x1="0" x2="1" y1="0" y2="0">
              {Array.from({ length: Samples + 1 }, (_, i) => (
                <stop
                  offset={String(i / Samples)}
                  stop-color={samples.pipe(map((all) => all?.[i]?.[1] ?? 'transparent'))}
                />
              ))}
            </linearGradient>
          </defs>
          <path class={classes.bars} d={derive(model.histogram, histogramPath)} />
          <path
            class={classes.area}
            classList={{ [classes.uncolored]: samples.pipe(map((s) => !s)) }}
            d={area}
            fill={`url(#${gradientId})`}
          />
          <path class={classes.curve} d={curve} />
        </svg>
        <For each={indices(points)} key={(i) => i}>
          {(i) => {
            const point = pointAt(i);
            return (
              <span
                class={classes.track}
                style={{ transform: point.pipe(map((p) => `translate(${p.x * 100}%, ${(1 - yOf(p.alpha)) * 100}%)`)) }}
              >
                <span
                  class={classes.point}
                  classList={{ [classes.on]: derive(picked.asReadonly, (p) => p?.kind === 'point' && p.index === i) }}
                  data-point={String(i)}
                  role="slider"
                  tabIndex={0}
                  aria-label={combineLatest([point, model.unit]).pipe(
                    map(([p, u]) => `Point at ${valueText(p.x, u)}, opacity ${alphaText(p.alpha)}`),
                  )}
                  onFocus={() => picked.set({ kind: 'point', index: i })}
                  onKeyDown={onPointKey(i)}
                />
              </span>
            );
          }}
        </For>
      </div>
      <div
        class={classes.strip}
        classList={{ [classes.themed]: samples.pipe(map((s) => !s)) }}
        style={{ background: strip }}
        title={derive(model.colors, (c) =>
          c.kind === 'value'
            ? 'Click to add a color · drag one · double-click to remove'
            : c.kind === 'uniform'
              ? 'One color: click to change it'
              : `Colored by ${c.label}`,
        )}
        onPointerDown={onStripDown}
        onPointerMove={onMove}
        onPointerUp={end}
        onPointerCancel={end}
      >
        <For each={indices(stops)} key={(i) => i}>
          {(i) => {
            const stop = stopAt(i);
            return (
              <span class={classes.stopTrack} style={{ transform: stop.pipe(map((s) => `translateX(${s.x * 100}%)`)) }}>
                <span
                  class={classes.stop}
                  classList={{ [classes.on]: derive(picked.asReadonly, (p) => p?.kind === 'stop' && p.index === i) }}
                  style={{ '--msp-stop': stop.pipe(map((s) => colorToHex(s.color))) }}
                  data-stop={String(i)}
                  role="button"
                  tabIndex={0}
                  aria-label={combineLatest([stop, model.unit]).pipe(map(([s, u]) => `Color at ${valueText(s.x, u)}`))}
                  onKeyDown={(e) => {
                    if (e.key === 'Delete' || e.key === 'Backspace') removeStop(i);
                    if (e.key === 'Enter' || e.key === ' ') picked.set({ kind: 'stop', index: i });
                  }}
                />
              </span>
            );
          }}
        </For>
      </div>
      <div class={classes.axis} aria-hidden="true">
        <span>{derive(model.unit, (u) => valueText(0, u))}</span>
        <span>{derive(model.unit, (u) => valueText(1, u))}</span>
      </div>
      {/* What is picked, in one quiet row: a point's value and opacity; a stop's (or the one) color, whose picker opens
          from the chevron only, as a form's color does. */}
      <Show when={derive(picked.asReadonly, (p) => p?.kind)}>
        {(kind) => (
          <div class={classes.picked}>
            {kind === 'point' ? (
              <span class={classes.pickedText}>
                {combineLatest([pickedPoint, model.unit]).pipe(
                  map(([p, u]) => (p ? `${valueText(p.x, u)} · opacity ${alphaText(p.alpha)}` : '')),
                )}
              </span>
            ) : (
              <button
                type="button"
                class={classes.colorRow}
                aria-expanded={colorOpen}
                onClick={() => colorOpen.set(!colorOpen.value)}
              >
                <span class={classes.swatch} style={{ background: pickedColor.pipe(map(colorToHex)) }} />
                <span class={classes.pickedText}>
                  {combineLatest([pickedStop, model.unit]).pipe(
                    map(([s, u]) => (s ? `Color at ${valueText(s.x, u)}` : 'Color')),
                  )}
                </span>
                <Icon name="chevron-down" size={14} class={classes.chevron} />
              </button>
            )}
            {kind !== 'uniform' ? (
              <button
                type="button"
                class={classes.tool}
                onClick={() => {
                  const s = picked.value;
                  if (s?.kind === 'point') removePoint(s.index);
                  else if (s?.kind === 'stop') removeStop(s.index);
                }}
              >
                Remove
              </button>
            ) : null}
            <Button
              variant="ghost"
              size="sm"
              icon={<Icon name="x" />}
              title="Done"
              onClick={() => picked.set(undefined)}
            />
          </div>
        )}
      </Show>
      <Show
        when={derive(
          [picked.asReadonly, colorOpen.asReadonly],
          (p, o) => o && (p?.kind === 'stop' || p?.kind === 'uniform'),
        )}
        reveal
      >
        {() => (
          <ColorPicker
            value={pickedColor}
            onChange={(color) => {
              const s = picked.value;
              if (s?.kind === 'stop') setStopColor(s.index, color);
              else model.change({ uniform: color });
            }}
          />
        )}
      </Show>
    </section>
  );
}

/**
 * Colors at peaks: the theme's own colors, or a palette (each its swatches), in a popover under the palette button; a
 * pick colors the peaks and closes it; Esc or a press outside closes it.
 */
function PalettePicker(props: { model: VolumeTransferModel; onClose(): void }) {
  const { model } = props;
  const pick = (palette?: string) => {
    model.colorsAtPeaks(palette);
    props.onClose();
  };
  const list = (
    <div
      class={classes.palettes}
      role="listbox"
      aria-label="Colors at peaks"
      onKeyDown={(e) => {
        if (e.key === 'Escape') {
          e.stopPropagation();
          props.onClose();
        }
      }}
    >
      <button type="button" role="option" class={classes.palette} onClick={() => pick()}>
        <span class={classes.paletteLabel}>The theme's colors</span>
      </button>
      {model.palettes.map((p) => (
        <button type="button" role="option" class={classes.palette} onClick={() => pick(p.id)}>
          <span class={classes.paletteSwatches} aria-hidden="true">
            {p.colors.map((c) => (
              <span class={classes.paletteSwatch} style={{ background: colorToHex(c) }} />
            ))}
          </span>
          <span class={classes.paletteLabel}>{p.label}</span>
        </button>
      ))}
    </div>
  ) as HTMLElement;
  // A press outside (not on the palette button, which toggles) closes it.
  const onDown = (e: PointerEvent) => {
    const path = e.composedPath();
    if (path.includes(list) || path.some((t) => (t as Element).getAttribute?.('aria-label') === 'Colors at peaks')) {
      return;
    }
    props.onClose();
  };
  document.addEventListener('pointerdown', onDown, true);
  onCleanup(() => document.removeEventListener('pointerdown', onDown, true));
  queueMicrotask(() => list.querySelector<HTMLElement>('button')?.focus());
  return list;
}
