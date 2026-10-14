/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { For, Show } from '../../dom/control-flow.js';
import type { MaybeReactive } from '../../dom/scope.js';
import { classes } from './_generated/legend.module.css.js';

/** One entry of a legend: what a color means. */
export interface LegendItem {
  label: string;
  /** A CSS color; `null` for "no color" (missing data), drawn as a quiet empty mark. */
  color: string | null;
  /** The mark: a square (default), a dot or a line. */
  glyph?: 'swatch' | 'dot' | 'line';
}

/** A continuous scale: the colors along it and what its two ends stand for. */
export interface LegendRamp {
  /** Colors at positions `t` in 0..1, in order. */
  stops: readonly { t: number; color: string }[];
  /** What the low and the high end mean (`0`, `100`, `Hydrophobic`, …). */
  lo: string;
  hi: string;
  caption?: string;
}

export interface LegendProps {
  /** A heading over the legend. */
  title?: MaybeReactive<string | undefined>;
  items?: MaybeReactive<readonly LegendItem[]>;
  ramps?: MaybeReactive<readonly LegendRamp[]>;
  class?: string;
}

function gradient(ramp: LegendRamp): string {
  const stops = ramp.stops.length > 0 ? ramp.stops : [{ t: 0, color: 'transparent' }];
  const parts = stops.map((s) => `${s.color} ${(Math.min(1, Math.max(0, s.t)) * 100).toFixed(1)}%`);
  return `linear-gradient(to right, ${parts.join(', ')})`;
}

/**
 * What the colors over the canvas mean: a swatch (or dot, or line) and a label per item, and for a continuous scale
 * a gradient with its low and high ends. Presentational and transparent: it takes its text color from where it sits, so
 * it reads on the canvas layer in its light and dark ink. Reactive: items and ramps may be streams.
 */
export function Legend(props: LegendProps) {
  const items = props.items ?? [];
  const ramps = props.ramps ?? [];
  return (
    <div class={props.class ? `${classes.legend} ${props.class}` : classes.legend} role="group" aria-label="Legend">
      <Show when={props.title}>{(title) => <span class={classes.title}>{title}</span>}</Show>
      <For each={items}>
        {(item) => (
          <div class={classes.item}>
            <span
              class={`${classes.mark} ${classes[item.glyph ?? 'swatch']}`}
              classList={{ [classes.none]: item.color === null }}
              style={{ '--legend-color': item.color ?? undefined }}
              aria-hidden="true"
            />
            <span class={classes.label}>{item.label}</span>
          </div>
        )}
      </For>
      <For each={ramps}>
        {(ramp) => (
          <div class={classes.ramp}>
            {ramp.caption ? <span class={classes.caption}>{ramp.caption}</span> : null}
            <span
              class={classes.bar}
              style={{ background: gradient(ramp) }}
              role="img"
              aria-label={`${ramp.lo} to ${ramp.hi}`}
            />
            <span class={classes.ends}>
              <span>{ramp.lo}</span>
              <span>{ramp.hi}</span>
            </span>
          </div>
        )}
      </For>
    </div>
  );
}
