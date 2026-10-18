/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Show } from '../../dom/control-flow.js';
import { effect, type MaybeReactive, requireScope, runInScope } from '../../dom/scope.js';
import { on } from '../../dom/scoped.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive, isState } from '../../reactive/state.js';
import { Icon, type IconName } from '../icon/icon.js';
import { classes } from './_generated/action-footer.module.css.js';

export interface ActionFooterProps {
  /** The action's icon and what applying does ("Download", "Create"). */
  icon: IconName;
  label: MaybeReactive<string>;
  /** Running: the action shows progress and is disabled. */
  busy?: MaybeReactive<boolean>;
  /** Not possible yet (not enough selected): the action is disabled, no progress. */
  disabled?: MaybeReactive<boolean>;
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
  const busy = held(props.busy);
  const disabled = derive([busy, held(props.disabled)], (b, d) => b || d);
  const scope = requireScope('ActionFooter');
  // Enter in a field of the block above it applies, as in Mol* v5's actions (a PDB id typed, Enter: it loads).
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Enter' || e.defaultPrevented || e.isComposing || e.altKey || e.ctrlKey || e.metaKey) return;
    const field = e.target as HTMLElement;
    if (!submitsOnEnter(field) || disabled.value) return;
    e.preventDefault();
    // The field's value first (forms commit text and numbers on `change`, which a browser may send after the key).
    field.dispatchEvent(new Event('change', { bubbles: true }));
    queueMicrotask(() => {
      if (!disabled.value) props.onApply();
    });
  };
  const listen = (footer: HTMLElement) =>
    // The block is its parent, there once the footer is placed in it.
    queueMicrotask(() => {
      const block = footer.parentElement;
      if (block && !scope.disposed) runInScope(scope, () => on(block, 'keydown', onKey));
    });
  return (
    <div class={classes.footer} ref={(el: HTMLElement) => listen(el)}>
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

/** A one-line field whose Enter applies: text, number and the like; not a text area, a list box, a check or a slider. */
function submitsOnEnter(el: HTMLElement): boolean {
  if (!(el instanceof HTMLInputElement)) return false;
  if (el.getAttribute('role') === 'combobox' || el.getAttribute('aria-expanded') === 'true') return false;
  return ['text', 'number', 'search', 'url', 'email', 'tel', ''].includes(el.type);
}

/** A flag as a state: a state as it is; a stream (or a plain value) held in one, with the current scope. */
function held(source: MaybeReactive<boolean> | undefined): ReadonlyState<boolean> {
  if (isState(source)) return source as ReadonlyState<boolean>;
  const state = observableState(false);
  effect(source ?? false, (v) => state.set(v));
  return state.asReadonly;
}
