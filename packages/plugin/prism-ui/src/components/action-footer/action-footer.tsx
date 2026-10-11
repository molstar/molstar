/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, map, type Observable, of } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import type { MaybeReactive } from '../../dom/scope.js';
import { Icon, type IconName } from '../icon/icon.js';
import { classes } from './_generated/action-footer.module.css.js';

export interface ActionFooterProps {
  /** The action's icon and what applying does ("Download", "Create"). */
  icon: IconName;
  label: MaybeReactive<string>;
  /** Running: the action shows progress and is disabled. */
  busy?: Observable<boolean>;
  /** Not possible yet (not enough selected): the action is disabled, no progress. */
  disabled?: Observable<boolean>;
  onApply(): void;
  /** Set defaults (the form's reset); without it, no reset button. */
  onReset?(): void;
  /** A quiet line on the left (the output's size). */
  info?: MaybeReactive<string>;
  /** A second, lesser action before the main one (Copy, Clear): in the text color. */
  secondary?: { icon: IconName; label: string; onClick(): void };
}

/**
 * The last row of a block with a form (Mol*'s action footer, flat): set defaults on the left, the action at the right
 * (its icon and label in the accent color, only they are the button); optionally a quiet info line and a lesser action
 * before it. Edge to edge: it undoes its container's 8px padding at the
 * sides and bottom. No lines or fills; hovering darkens the text.
 */
export function ActionFooter(props: ActionFooterProps) {
  const busy = props.busy ?? of(false);
  const disabled = combineLatest([busy, props.disabled ?? of(false)]).pipe(map(([b, d]) => b || d));
  return (
    <div class={classes.footer}>
      {props.onReset ? (
        <button
          type="button"
          class={classes.reset}
          data-tip="Set defaults"
          aria-label="Set defaults"
          onClick={() => props.onReset?.()}
        >
          {/* Not the refresh arrow: back to the defaults (a list, reset). */}
          <Icon name="list-restart" size={15} />
        </button>
      ) : null}
      {props.info !== undefined ? <span class={classes.info}>{props.info}</span> : null}
      {props.secondary ? (
        <button type="button" class={classes.secondary} disabled={busy} onClick={() => props.secondary?.onClick()}>
          <span>{props.secondary.label}</span>
        </button>
      ) : null}
      <button type="button" class={classes.apply} disabled={disabled} onClick={() => props.onApply()}>
        {/* No icon: the label says it; a spinner only while it runs. */}
        <Show when={busy}>{() => <Icon name="loader-circle" size={15} class={classes.spin} />}</Show>
        <span>{props.label}</span>
      </button>
    </div>
  );
}
