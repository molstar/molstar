/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { type MaybeReactive, on, requireScope, runInScope, Show } from '../../dom/index.js';
import { held } from '../../state.js';
import { Icon } from '../icon/icon.js';
import { classes } from './_generated/action-footer.module.css.js';

export interface ActionFooterProps {
  /** What applying does ("Download", "Create"). */
  label: MaybeReactive<string>;
  /** Running: the action shows progress and is disabled. */
  busy?: MaybeReactive<boolean>;
  /** Not possible yet (not enough selected): the action is disabled, no progress. */
  disabled?: MaybeReactive<boolean>;
  onApply(): void;
  /** Set defaults (the form's reset); without it, no reset button. */
  onReset?(): void;
  /** A quiet line after the reset (the output's size). */
  info?: MaybeReactive<string>;
  /** A second, lesser action before the main one (Copy, Clear): in the muted text color. */
  secondary?: { label: string; onClick(): void };
}

/**
 * The last row of a block with a form: set defaults at the start, the action at the end (its label in the accent
 * color, only it is the button); optionally a quiet info line and a lesser action before the main one. A row (styles
 * §4, §5): as tall as one and padded as one, so the reset sits in the icon slot and the info line starts at the text
 * column of the rows above it, and the action's label ends at the inset. No lines or fills; hovering underlines.
 *
 * Enter in a one-line field of the block it ends applies it (a code typed, Enter: it loads).
 */
export function ActionFooter(props: ActionFooterProps) {
  const busy = held(props.busy, false);
  const disabled = derive([busy, held(props.disabled, false)], (b, d) => b || d);
  const scope = requireScope('ActionFooter');
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
          <Icon name="list-restart" />
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
        <Show when={busy}>{() => <Icon name="loader-circle" class={classes.spin} />}</Show>
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
