/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, ColorPalette, ColorPicker, colorToHex } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import { classes } from './_generated/quick-controls.module.css.js';

// The quick controls above a full form (a map's levels, a structure representation's): the settings people change
// most, at hand. A color is a row of swatches (one click applies one), the whole picker behind its chevron; a number is
// a slider with its value. Pointer-rate changes go to `onInput`, the end of a gesture to `onCommit` (one undo step).

/** The colors offered at once: Mol*'s density blue, the palette's hues, a light grey. */
export const QuickColors: readonly number[] = [0x3362b2, ...ColorPalette, 0xcccccc];

export function QuickColorRow(props: {
  label: string;
  value: ReadonlyState<number | undefined>;
  /** A swatch picked (one click: a whole change). */
  onPick(color: number): void;
  /** The picker dragged across (pointer rate). */
  onInput(color: number): void;
  /** A close at the row's end (an editor opened from a swatch); none: the row stays. */
  onClose?(): void;
}) {
  const more = observableState(false);
  const shown = derive(props.value, (v) => v ?? 0x808080);
  return (
    <div class={classes.colors}>
      <div class={classes.quick} role="listbox" aria-label={props.label}>
        {QuickColors.map((c) => (
          <button
            type="button"
            role="option"
            class={classes.swatch}
            style={{ background: colorToHex(c) }}
            aria-label={colorToHex(c)}
            aria-selected={derive(props.value, (v) => v === c)}
            onClick={() => props.onPick(c)}
          />
        ))}
        <button
          type="button"
          class={classes.more}
          aria-expanded={more}
          title="More colors"
          onClick={() => more.set(!more.value)}
        >
          <Icon name="chevron-down" size={14} class={classes.chevron} />
        </button>
        {props.onClose ? (
          <Button variant="ghost" size="sm" icon={<Icon name="x" />} title="Done" onClick={() => props.onClose?.()} />
        ) : null}
      </div>
      <Show when={more} reveal>
        {() => <ColorPicker value={shown} onChange={(c) => props.onInput(c)} />}
      </Show>
    </div>
  );
}

/** A labelled slider with its value; `format` writes the value. */
export function QuickRange(props: {
  label: string;
  min: number;
  max: number;
  step: number;
  value: ReadonlyState<number>;
  format(value: number): string;
  onInput(value: number): void;
  onCommit(): void;
}) {
  const span = props.max - props.min || 1;
  const value = props.value;
  return (
    <label class={classes.range}>
      <span class={classes.rangeLabel}>{props.label}</span>
      <input
        class={classes.slider}
        type="range"
        min={String(props.min)}
        max={String(props.max)}
        step={String(props.step)}
        value={derive(value, String)}
        style={{ '--msp-progress': derive(value, (v) => `${(100 * (v - props.min)) / span}%`) }}
        onInput={(e) => props.onInput(Number(e.currentTarget.value))}
        onChange={() => props.onCommit()}
      />
      <span class={classes.rangeValue}>{derive(value, (v) => props.format(v))}</span>
    </label>
  );
}
