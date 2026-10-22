/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { colorToHex } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { effect, frameState } from '@molstar/prism-ui/dom/scope';
import { drag } from '@molstar/prism-ui/dom/scoped';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type {
  VolumeLevel,
  VolumeLevelHistogram,
  VolumeLevelsModel,
  VolumeLevelUnit,
} from '../../model/volume-levels.js';
import { QuickColorRow, QuickRange } from '../common/quick-controls.js';
import { classes } from './_generated/volume-levels.module.css.js';

// A map's levels (model/volume-levels.ts; v5's density controls, improved): the map's values as a histogram with a
// marker per level, dragged (or stepped with the arrow keys) to set the iso value, the part of the map the surface
// encloses tinted in its color; under it a row per level: its color (a click opens the picker and the opacity), its
// value typed in σ or map units, the other unit before it. A drag is one undo step.

/** σ steps for the arrow keys (Shift: the big one). */
const Step = 0.1;
const BigStep = 1;

const minus = (s: string) => s.replace('-', '−');

/** "1.50 σ", "0.0734": a value in a unit. */
export function formatIn(unit: VolumeLevelUnit, value: number, stats: VolumeLevelsModel['stats']): string {
  if (unit === 'sigma') return minus(((value - stats.mean) / (stats.sigma || 1)).toFixed(2));
  return minus(String(Number(value.toPrecision(3))));
}

/** A value in σ as the number ARIA wants (a plain minus; the typographic one is for the text, `formatIn`). */
const sigmas = (value: number, stats: VolumeLevelsModel['stats']) =>
  ((value - stats.mean) / (stats.sigma || 1)).toFixed(2);

/** A typed value in a unit, as map units; undefined: not a number. */
function parseIn(unit: VolumeLevelUnit, text: string, stats: VolumeLevelsModel['stats']): number | undefined {
  const v = Number(text.replace('−', '-').replace('σ', '').trim());
  if (!text.trim() || !Number.isFinite(v)) return undefined;
  return unit === 'sigma' ? stats.mean + v * stats.sigma : v;
}

/** The histogram as a filled step path in a 1 × 1 box (log counts: a map's tails are thin). */
export function histogramPath(h: VolumeLevelHistogram | undefined): string {
  if (!h || !h.counts.length) return '';
  const n = h.counts.length;
  let top = 0;
  for (let i = 0; i < n; i++) top = Math.max(top, h.counts[i]);
  const scale = top > 0 ? 1 / Math.log1p(top) : 0;
  const parts = ['M0 1'];
  for (let i = 0; i < n; i++) {
    const y = (1 - Math.log1p(h.counts[i]) * scale).toFixed(4);
    parts.push(`L${(i / n).toFixed(4)} ${y}L${((i + 1) / n).toFixed(4)} ${y}`);
  }
  parts.push('L1 1Z');
  return parts.join('');
}

export function VolumeLevels(props: { model: VolumeLevelsModel }) {
  const { model } = props;
  const { stats } = model;
  const span = stats.max - stats.min || 1;
  const clamp = (v: number) => Math.min(stats.max, Math.max(stats.min, v));
  /** A value's place on the plot, 0 to 1. */
  const at = (v: number) => (clamp(v) - stats.min) / span;
  /** A level by id, as a state: the last one found (one removed keeps its last until its rows go). */
  const levelOf = (id: string): ReadonlyState<VolumeLevel> => {
    let last = model.levels.value.find((l) => l.id === id) as VolumeLevel;
    return derive(model.levels, (ls) => (last = ls.find((l) => l.id === id) ?? last));
  };
  const ids = derive(
    model.levels,
    (ls) => ls.map((l) => l.id),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );

  // Dragging on the plot: the marker under the pointer, else the nearest, follows it; once a frame (the moves, counted,
  // are a frame state: its value is only a tick). The drag's listeners are attached once (`drag`, dom/scoped.ts).
  let rect: DOMRect | undefined;
  let dragging: string | undefined;
  let pointerX = 0;
  let moves = 0;
  const moved = frameState(0);
  const valueAt = (clientX: number) =>
    clamp(stats.min + span * Math.min(1, Math.max(0, (clientX - rect!.left) / (rect!.width || 1))));
  const step = () => {
    if (dragging) model.change(dragging, { value: valueAt(pointerX) });
  };
  effect(moved, step);
  const plotDrag = (el: HTMLElement) =>
    drag(el, {
      start: (e) => {
        if (e.button !== 0) return undefined;
        rect = el.getBoundingClientRect();
        const levels = model.levels.value;
        const hit = (e.target as HTMLElement).closest<HTMLElement>('[data-level]')?.dataset.level;
        const x = (e.clientX - rect.left) / (rect.width || 1);
        const nearest = levels.reduce<VolumeLevel | undefined>(
          (best, l) => (!best || Math.abs(at(l.value) - x) < Math.abs(at(best.value) - x) ? l : best),
          undefined,
        );
        dragging = hit ?? nearest?.id;
        if (!dragging) return undefined;
        e.preventDefault();
        el.querySelector<HTMLElement>(`[data-level="${CSS.escape(dragging)}"]`)?.focus({ preventScroll: true });
        pointerX = e.clientX;
        step();
        return dragging;
      },
      move: (e) => {
        pointerX = e.clientX;
        moved.offer(++moves);
      },
      // Ended by the pointer: the last place now. Gone mid-drag (another row is picked, the node removed, the panel
      // hidden): no more steps. Either way the step closes.
      end: (_id, _how, e) => {
        if (e) step();
        dragging = undefined;
        rect = undefined;
        model.settle();
      },
    });
  const onKey = (id: string) => (e: KeyboardEvent) => {
    const level = model.levels.value.find((l) => l.id === id);
    if (!level) return;
    const step = (e.shiftKey ? BigStep : Step) * (stats.sigma || span / 100);
    const to =
      e.key === 'ArrowRight' || e.key === 'ArrowUp'
        ? level.value + step
        : e.key === 'ArrowLeft' || e.key === 'ArrowDown'
          ? level.value - step
          : e.key === 'PageUp'
            ? level.value + BigStep * (stats.sigma || span / 10)
            : e.key === 'PageDown'
              ? level.value - BigStep * (stats.sigma || span / 10)
              : e.key === 'Home'
                ? stats.min
                : e.key === 'End'
                  ? stats.max
                  : undefined;
    if (to === undefined) return;
    e.preventDefault();
    // Keys held or pressed in a row are one step (the model's quiet ends it).
    model.change(id, { value: clamp(to) });
  };

  const unitLabel = (u: VolumeLevelUnit) => (u === 'sigma' ? 'σ' : 'abs');
  return (
    <section class={classes.levels} aria-label={model.title ? `${model.title} levels` : 'Levels'}>
      <header class={classes.head}>
        <span class={classes.title}>{model.title ?? 'Levels'}</span>
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
              {unitLabel(u)}
            </button>
          ))}
        </span>
      </header>
      <div class={classes.plot} ref={plotDrag}>
        <svg class={classes.bars} viewBox="0 0 1 1" preserveAspectRatio="none" aria-hidden="true">
          <path d={derive(model.histogram, histogramPath)} />
        </svg>
        <span class={classes.mean} style={{ transform: `translateX(${at(stats.mean) * 100}%)` }} aria-hidden="true" />
        <For each={ids} by={(id) => id}>
          {(id) => {
            const level = levelOf(id);
            const color = derive(level, (l) => (l.color === undefined ? undefined : colorToHex(l.color)));
            const x = derive(level, (l) => at(l.value));
            const hidden = derive(level, (l) => l.hidden);
            // What the surface encloses: above the level, or below a negative one (Fo-Fc's −).
            const inside = derive(level, (l) => (l.value >= stats.mean ? at(l.value) * 100 : at(l.value) * 100 - 100));
            return (
              <>
                <span
                  class={classes.inside}
                  classList={{ [classes.off]: hidden }}
                  style={{ transform: derive(inside, (p) => `translateX(${p}%)`), '--msp-level': color }}
                  aria-hidden="true"
                />
                <span
                  class={classes.track}
                  style={{ transform: derive(x, (p) => `translateX(${p * 100}%)`), '--msp-level': color }}
                >
                  <span
                    class={classes.marker}
                    classList={{ [classes.off]: hidden }}
                    data-level={id}
                    role="slider"
                    tabIndex={0}
                    aria-label={derive(level, (l) => `${l.label} iso value`)}
                    aria-valuemin={sigmas(stats.min, stats)}
                    aria-valuemax={sigmas(stats.max, stats)}
                    aria-valuenow={derive(level, (l) => sigmas(l.value, stats))}
                    aria-valuetext={derive(level, (l) => `${formatIn('sigma', l.value, stats)} σ`)}
                    onKeyDown={onKey(id)}
                  />
                </span>
              </>
            );
          }}
        </For>
      </div>
      <div class={classes.axis} aria-hidden="true">
        <span>{derive(model.unit, (u) => formatIn(u, stats.min, stats))}</span>
        <span>{derive(model.unit, (u) => formatIn(u, stats.max, stats))}</span>
      </div>
      <For each={ids} by={(id) => id}>
        {(id) => <LevelRow model={model} id={id} level={levelOf(id)} />}
      </For>
    </section>
  );
}

/** A level's row: its color (opening the picker and the opacity), its name, its value in the unit and the other. */
function LevelRow(props: { model: VolumeLevelsModel; id: string; level: ReadonlyState<VolumeLevel> }) {
  const { model, id, level } = props;
  const { stats } = model;
  const open = observableState(false);
  const color = derive(level, (l) => (l.color === undefined ? undefined : colorToHex(l.color)));
  return (
    <div class={classes.level}>
      <div class={classes.row}>
        <button
          type="button"
          class={classes.swatch}
          classList={{ [classes.themed]: derive(color, (c) => c === undefined) }}
          style={{ '--msp-level': color }}
          aria-expanded={open}
          title={derive(color, (c) => (c === undefined ? 'Colored by a theme: pick a color' : 'Color and opacity'))}
          onClick={() => open.set(!open.value)}
        />
        <span class={classes.name} classList={{ [classes.dim]: derive(level, (l) => l.hidden) }}>
          {derive(level, (l) => l.label)}
        </span>
        {/* The other unit before the field: its width changes as the value does, and the field stays put. */}
        <span class={classes.other}>
          {derive([level, model.unit], (l, u) =>
            u === 'sigma' ? formatIn('absolute', l.value, stats) : `${formatIn('sigma', l.value, stats)} σ`,
          )}
        </span>
        <input
          class={classes.value}
          inputMode="decimal"
          spellcheck={false}
          aria-label={derive(level, (l) => `${l.label} iso value`)}
          value={derive([level, model.unit], (l, u) => formatIn(u, l.value, stats))}
          onChange={(e) => {
            const v = parseIn(model.unit.value, e.currentTarget.value, stats);
            if (v === undefined) return;
            model.change(id, { value: Math.min(stats.max, Math.max(stats.min, v)) });
            model.settle();
          }}
          onKeyDown={(e) => {
            if (e.key === 'Enter') e.currentTarget.blur();
          }}
        />
        <span class={classes.suffix}>{derive(model.unit, (u) => (u === 'sigma' ? 'σ' : ''))}</span>
      </div>
      <Show when={open} reveal>
        {() => (
          <div class={classes.editor}>
            <QuickColorRow
              label={`${model.levels.value.find((l) => l.id === id)?.label ?? 'Level'} color`}
              value={derive(level, (l) => l.color)}
              onPick={(c) => {
                model.change(id, { color: c });
                model.settle();
              }}
              onInput={(c) => model.change(id, { color: c })}
              onClose={() => open.set(false)}
            />
            <QuickRange
              label="Opacity"
              min={0}
              max={1}
              step={0.01}
              value={derive(level, (l) => l.opacity)}
              format={(v) => `${Math.round(v * 100)}%`}
              onInput={(v) => model.change(id, { opacity: v })}
              onCommit={() => model.settle()}
            />
          </div>
        )}
      </Show>
    </div>
  );
}
