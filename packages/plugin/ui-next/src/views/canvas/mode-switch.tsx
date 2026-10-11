/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { Show } from '../../dom/control-flow.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { classes } from './_generated/mode-switch.module.css.js';

/** What a click on the canvas does: interact (focus, rotate, measure) or select. */
export type CanvasMode = 'interact' | 'select';

export interface ModeSwitchProps {
  mode: ReadonlyState<CanvasMode>;
  setMode(mode: CanvasMode): void;
  /** What a click selects (residue, chain, …), shown next to Select while it is on. */
  granularity?: ReadonlyState<string>;
  granularities?: readonly (readonly [value: string, label: string])[];
  setGranularity?(value: string): void;
}

/**
 * The canvas mode, always visible at the top right of the canvas: two labelled segments, so the difference between
 * interacting with the scene and selecting in it is visible before the first click. While selecting, the switch also
 * shows what a click selects.
 */
export function ModeSwitch(props: ModeSwitchProps) {
  const selecting = props.mode.pipe(
    map((m) => m === 'select'),
    distinctUntilChanged(),
  );
  return (
    <div class={classes.switch} role="group" aria-label="Canvas mode" data-canvas-mode={props.mode}>
      <button
        type="button"
        class={classes.segment}
        aria-pressed={selecting.pipe(map((s) => !s))}
        title="Interact: drag to rotate, click to focus"
        onClick={() => props.setMode('interact')}
      >
        <Icon name="rotate-3d" size={15} />
        <span class={classes.label}>Interact</span>
      </button>
      <button
        type="button"
        class={classes.segment}
        aria-pressed={selecting}
        title="Select: click to select, Shift+click to extend; drag still rotates"
        onClick={() => props.setMode('select')}
      >
        <Icon name="mouse-pointer-click" size={15} />
        <span class={classes.label}>Select</span>
      </button>
      {props.granularities && props.granularity ? (
        <Show when={selecting}>{() => <Granularity {...props} />}</Show>
      ) : null}
    </div>
  );
}

function Granularity(props: ModeSwitchProps) {
  const options = props.granularities!;
  const index = props.granularity!.pipe(map((g) => String(options.findIndex(([v]) => v === g))));
  return (
    <select
      class={classes.granularity}
      aria-label="Select by"
      title="What a click selects"
      value={index}
      onChange={(e) => props.setGranularity?.(options[+(e.currentTarget as HTMLSelectElement).value][0])}
    >
      {options.map(([, label], i) => (
        <option value={String(i)}>{label}</option>
      ))}
    </select>
  );
}
