/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { Icon, type IconName } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { classes } from './_generated/mode-switch.module.css.js';

/** What a click on the canvas does: interact (focus, rotate, measure) or select. */
export type CanvasMode = 'interact' | 'select';

/** A picking level (Mol*'s granularity): its value, label, icon, and an optional group (multi-instance levels). */
export interface Granularity {
  value: string;
  label: string;
  icon: IconName;
  group?: string;
}

export interface ModeSwitchProps {
  mode: ReadonlyState<CanvasMode>;
  setMode(mode: CanvasMode): void;
  /** What a click selects (residue, chain, …), offered below the switch while selecting. */
  granularity?: ReadonlyState<string>;
  granularities?: readonly Granularity[];
  setGranularity?(value: string): void;
}

/**
 * The canvas mode, at the top of the canvas stack: two toggle buttons, Interact and Select (the titles explain the
 * difference; Select takes the accent and carries the picking level as a small badge). `PickLevel` changes the level:
 * it shows beside Select while the pointer is on Select or on it.
 */
export function ModeSwitch(props: ModeSwitchProps) {
  const selecting = props.mode.pipe(
    map((m) => m === 'select'),
    distinctUntilChanged(),
  );
  // The picking level as a badge on Select: what a click picks, at a glance.
  const levels = props.granularities;
  const current =
    levels && props.granularity
      ? props.granularity.pipe(map((g) => levels.find((l) => l.value === g) ?? levels[0]))
      : undefined;
  return (
    <div class={classes.switch} role="group" aria-label="Canvas mode" data-canvas-mode={props.mode}>
      <button
        type="button"
        class={classes.segment}
        aria-pressed={selecting.pipe(map((s) => !s))}
        data-tip="Interact: drag to rotate, click to focus"
        data-tip-side="top"
        aria-label="Interact"
        onClick={() => props.setMode('interact')}
      >
        <Icon name="rotate-3d" size={20} />
      </button>
      <button
        type="button"
        class={`${classes.segment} ${classes.select}`}
        aria-pressed={selecting}
        data-tip={current ? undefined : 'Select: click to select, Shift+click to extend'}
        data-tip-side="top"
        aria-label="Select"
        onClick={() => props.setMode('select')}
      >
        <Icon name="mouse-pointer-click" size={20} />
        {current ? <Show when={current}>{(l) => <Icon name={l.icon} size={10} class={classes.badge} />}</Show> : null}
      </button>
    </div>
  );
}

/**
 * The picking levels, above the canvas toolbar over Select (outside the toolbar's frame, so it is not clipped): the list
 * itself, one click from the toolbar: the current level marked; picking one sets it and turns select mode on. Shown while
 * the pointer is on Select or on it, or Select or a row has the keyboard focus (Tab goes from Select into the list);
 * on touch screens (no hover) a tap on Select opens it. It opens upward.
 */
export function PickLevel(props: ModeSwitchProps) {
  const levels = props.granularities;
  const granularity = props.granularity;
  if (!levels || !granularity) return null;
  return (
    <div class={classes.level}>
      <div class={classes.levels} role="menu" aria-label="What a click selects">
        {levels.map((l, i) => (
          <>
            {l.group && l.group !== levels[i - 1]?.group ? <div class={classes.group}>{l.group}</div> : null}
            <button
              type="button"
              role="menuitemradio"
              class={classes.levelRow}
              aria-checked={granularity.pipe(map((g) => g === l.value))}
              onClick={() => {
                props.setGranularity?.(l.value);
                props.setMode('select');
              }}
            >
              <Icon name={l.icon} size={15} />
              {l.label}
            </button>
          </>
        ))}
      </div>
    </div>
  );
}
