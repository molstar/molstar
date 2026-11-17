/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { type Child, type MaybeReactive, Show } from '../../dom/index.js';
import { held } from '../../state.js';
import { classes } from './_generated/skeleton.module.css.js';

/** Text sizes a placeholder line takes (its line box: the size's leading). */
export type SkeletonTextSize = 'xs' | 'sm' | 'md';

/**
 * Where content is loading: holds placeholder shapes (hidden from assistive technology) and its words for it, `label`
 * ("Loading" by default), in a polite status region: said once as it appears. It shows only after a moment
 * (`--msui-skeleton-delay`, 300ms), so a quick load never flashes it. Put it where the content goes and swap it for the
 * content (`Loaded` does both, and marks the area busy meanwhile).
 */
export function Loading(props: { label?: string; children: Child; class?: string }) {
  return (
    <div class={props.class ? `${classes.loading} ${props.class}` : classes.loading} role="status">
      <span class={classes.label}>{props.label ?? 'Loading'}</span>
      {props.children}
    </div>
  );
}

/**
 * Content that loads, with its placeholder in place until it does: while `when` is falsy, `placeholder` (shapes, in a
 * `Loading` region); then `children` with the value. The two share one place (a grid cell), so the change is a
 * crossfade in place (styles §7.1): the placeholder fades out (`--msui-exit`) as the content fades in, and nothing
 * moves when the shapes match the content. Content there from the first render shows at once.
 */
export function Loaded<T>(props: {
  when: MaybeReactive<T>;
  placeholder: () => Child;
  children: (value: NonNullable<T>) => Child;
  /** What the loading region is called for assistive technology (default "Loading"). */
  label?: string;
  class?: string;
}) {
  const value = held<T | undefined>(props.when, undefined);
  return (
    // The area is busy while its content loads (`aria-busy` on what loads, not on the status region: busy there would
    // hold its words back).
    <div
      class={props.class ? `${classes.swap} ${props.class}` : classes.swap}
      aria-busy={derive(value, (v) => (v ? undefined : 'true'))}
    >
      <Show when={derive(value, (v) => !v)} reveal="fade">
        {() => <Loading label={props.label}>{props.placeholder()}</Loading>}
      </Show>
      <Show when={value} reveal="fade">
        {(v) => <div class={classes.layer}>{props.children(v as NonNullable<T>)}</div>}
      </Show>
    </div>
  );
}

/**
 * One placeholder shape: a block of the given size (a plot, an image, a control), `md` corners by default; `round`
 * for an avatar or a dot. Sizes are CSS lengths (on the grid: tokens or multiples of 4px).
 */
export function Skeleton(props: { width?: string; height?: string; round?: boolean; radius?: 'sm' | 'md' | 'lg' }) {
  const radius = props.round ? '50%' : `var(--msui-radius-${props.radius ?? 'md'})`;
  return (
    <span
      class={classes.shape}
      aria-hidden="true"
      style={{ width: props.width ?? '100%', height: props.height ?? 'var(--msui-row)', borderRadius: radius }}
    />
  );
}

/** Widths of successive lines and rows: varied, so placeholders read as text rather than as a grid of bars. */
const Widths = ['72%', '56%', '84%', '64%', '48%', '76%'];

/**
 * Lines of text in a size's line boxes, a bar on each line; the last of several is shorter, as a paragraph's last line.
 */
export function SkeletonText(props: { lines?: number; size?: SkeletonTextSize; width?: string }) {
  const lines = Math.max(1, props.lines ?? 1);
  const out: Child[] = [];
  for (let i = 0; i < lines; i++) {
    const width = lines > 1 && i === lines - 1 ? '40%' : (props.width ?? (lines > 1 ? '100%' : Widths[0]));
    out.push(
      <span class={classes.line} data-size={props.size ?? 'sm'} aria-hidden="true">
        <span class={classes.shape} style={{ width }} />
      </span>,
    );
  }
  return out;
}

/**
 * List rows (styles §5, the row primitive): a row's height, its icon slot at the inset (a square of the icon's size;
 * `icon: false` for a list without icons) and a bar on its text column, the widths varied.
 */
export function SkeletonRows(props: { count?: number; icon?: boolean; round?: boolean }) {
  const count = Math.max(1, props.count ?? 3);
  const icon = props.icon ?? true;
  const out: Child[] = [];
  for (let i = 0; i < count; i++) {
    out.push(
      <div class={icon ? 'msui-row' : 'msui-row msui-row-plain'} aria-hidden="true">
        {icon ? (
          <span class="msui-slot">
            <span class={`${classes.shape} ${classes.icon}${props.round ? ` ${classes.round}` : ''}`} />
          </span>
        ) : null}
        <span class={`${classes.shape} ${classes.text}`} style={{ width: Widths[i % Widths.length] }} />
        <span />
      </div>,
    );
  }
  return out;
}
