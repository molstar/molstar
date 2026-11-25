/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/core/reactive';
import { observeResize, setStyle, uniqueId } from '../../dom/index.js';
import { Icon, type IconName } from '../icon/icon.js';
import { classes } from './_generated/choices.module.css.js';
import { radioGroup } from './radio-group.js';

export interface Choice {
  id: string;
  label: string;
  icon?: IconName;
  /** Its tooltip (`data-tip`). An icon-only choice shows its label there without one. */
  title?: string;
  /** Not available now: greyed, not picked, skipped by the arrows; a string says why (its tooltip). */
  disabled?: boolean | string;
}

export interface ChoicesProps {
  /** What the choices are: in the label column (`--msui-label`), so a `Choices` lines up with a form's rows. */
  label?: string;
  choices: readonly Choice[];
  /**
   * The choice in effect: given, the choices are one single-select group (a radio group: one Tab stop, the arrows,
   * Home and End move and pick); without it, a row of one-click actions (plain buttons).
   */
  active?: ReadonlyState<string | undefined>;
  onPick(id: string): void;
  /** Only the icons: each choice's label is its accessible name and its tooltip. */
  iconOnly?: boolean;
  /**
   * This many equal columns, wrapping into more rows; by default as many as fit, each at least as wide as the widest
   * choice (labels are never cut: the row wraps instead).
   */
  columns?: number;
  /** The group's accessible name when it has no `label`. */
  'aria-label'?: string;
  class?: string;
}

/**
 * A labelled row of equal-width choices (quick styles: Default, Cartoon, Spacefill, Surface; wiggle: Dynamics,
 * Uncertainty, Clear), at control height, each row of them a form row (styles §4): equal columns that share the width,
 * each at least as wide as the widest choice, wrapping into more rows when they do not fit; a label is never cut. With `active`, a single-select group (WAI-ARIA radio group: `role="radio"`, `aria-checked`, a
 * roving tab stop on the checked choice, else the first enabled; `radioGroup`): a change rewrites the two choices it
 * concerns. Without it, plain buttons.
 */
export function Choices(props: ChoicesProps) {
  const { active, choices, iconOnly } = props;
  const single = active !== undefined;
  const labelId = props.label ? uniqueId('msui-choices') : undefined;
  const buttons: HTMLButtonElement[] = [];

  const pick = (i: number) => {
    const c = choices[i];
    if (!c || c.disabled) return;
    props.onPick(c.id);
  };

  const choice = (c: Choice, i: number) => {
    const reason = typeof c.disabled === 'string' ? c.disabled : undefined;
    const tip = reason ?? c.title ?? (iconOnly ? c.label : undefined);
    const el = (
      <button
        type="button"
        class={classes.choice}
        role={single ? 'radio' : undefined}
        aria-label={iconOnly ? c.label : undefined}
        aria-disabled={c.disabled ? 'true' : undefined}
        data-tip={tip}
        onClick={() => pick(i)}
      >
        {/* Its content at its own width: what the row measures for its columns. */}
        <span class={classes.content}>
          {c.icon ? <Icon name={c.icon} /> : null}
          {iconOnly ? null : c.label}
        </span>
      </button>
    ) as HTMLButtonElement;
    buttons.push(el);
    return el;
  };

  const row = (
    <div
      ref={(el: HTMLElement) => {
        // Equal columns at least as wide as the widest choice's content (a UI label is never cut, styles §2): the row's
        // minimum column, measured from the contents (their intrinsic width: the row's own width does not change it),
        // when it is first laid out and on each resize after (a font that loads late).
        observeResize(el, () => {
          let widest = 0;
          for (const b of buttons) {
            const content = b.firstElementChild as HTMLElement;
            const pad = Number.parseFloat(getComputedStyle(b).paddingInlineStart) * 2;
            widest = Math.max(widest, Math.ceil(content.getBoundingClientRect().width + pad));
          }
          if (widest > 0) setStyle(el, '--msui-choices-min', `${widest}px`);
        });
      }}
      class={classes.row}
      role={single ? 'radiogroup' : 'group'}
      aria-labelledby={labelId}
      aria-label={labelId ? undefined : props['aria-label']}
      data-icons={iconOnly ? '' : undefined}
      data-columns={props.columns ? '' : undefined}
      style={props.columns ? { '--msui-choices-columns': String(props.columns) } : undefined}
      onKeyDown={single ? (e: KeyboardEvent) => onKey?.(e) : undefined}
    >
      {choices.map(choice)}
    </div>
  );

  // The radio group's marks and keys (radio-group.ts), shared with the other single-select rows.
  const onKey = single
    ? radioGroup({
        keys: choices.map((c) => c.id),
        buttons,
        disabled: (i) => !!choices[i].disabled,
        active,
        pick,
      })
    : undefined;

  return (
    <div
      class={props.class ? `${classes.choices} ${props.class}` : classes.choices}
      data-labelled={props.label ? '' : undefined}
    >
      {props.label ? (
        <span class={classes.label} id={labelId}>
          {props.label}
        </span>
      ) : null}
      {row}
    </div>
  );
}
