/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, isState, State } from '@molstar/core/reactive';
import { type MaybeReactive, toStream, uniqueId } from '../../dom/index.js';
import { held } from '../../state.js';
import { radioGroup } from '../choices/radio-group.js';
import { Icon, type IconName } from '../icon/icon.js';
import { classes } from './_generated/segmented.module.css.js';

export interface SegmentedOption<T> {
  value: T;
  /** Its name: the segment's text, or with `iconOnly` its accessible name and tooltip; the heading shows it. */
  label: string;
  icon?: IconName;
  /** Not available now: greyed, not picked, skipped by the arrows; a string says why (its tooltip). */
  disabled?: boolean | string;
}

export interface SegmentedProps<T> {
  options: readonly SegmentedOption<T>[];
  /**
   * The option in effect. A `State` is set when a segment is picked; a plain value is where the picker starts (it
   * keeps the pick itself); a read-only state or a stream is followed (`onChange` updates its source).
   */
  value: State<T> | MaybeReactive<T | undefined>;
  /** Told the picked value (after a `State` or a plain value took it). */
  onChange?(value: T): void;
  /** What is picked: a heading over the segments, the chosen option's name after it; the group's accessible name. */
  label?: string;
  /** Only the icons: each option's label is its accessible name and its tooltip. */
  iconOnly?: boolean;
  /**
   * `md` (default): segments a row high (24px) in a well a header high (32px), 16px icons (an icon names an object, a
   * representation type); `sm`: segments a control high (20px) in a well a row high, 12px icons, for a form row.
   */
  size?: 'md' | 'sm';
  /** The group's accessible name when it has no `label`. */
  'aria-label'?: string;
  class?: string;
}

/**
 * A segmented picker: one option of a few (a representation's type: cartoon, ball & stick, spacefill, surface), as a
 * row of equal segments in a quiet well (`--msui-surface-sunken`), the chosen one raised in it (a raised segment with a
 * hairline, its icon and text in the accent); with `label`, a heading over it, the label muted, then the chosen
 * option's name. Segments share the width equally, each at least as wide as its content (a label is never cut, styles
 * §2): a picker with long text labels needs their room (`Choices` wraps instead). A single-select WAI-ARIA radio group
 * (`radioGroup`: one tab stop, the arrows, Home and End move and pick), so a change rewrites the two segments it
 * concerns. Use it where the options are few and seen at once, their icons the way they are told apart; `Choices` for a
 * form's row of labelled choices or of actions.
 */
export function Segmented<T>(props: SegmentedProps<T>) {
  const { options, iconOnly } = props;
  // A writable state or a plain value takes the pick here; a read-only source is the caller's to update.
  const own =
    props.value instanceof State
      ? (props.value as State<T | undefined>)
      : !isState(props.value) && !toStream(props.value)
        ? new State(props.value as T | undefined)
        : undefined;
  const value = own ?? held<T | undefined>(props.value as MaybeReactive<T | undefined>, undefined);
  const labelId = props.label ? uniqueId('msui-segmented') : undefined;
  const buttons: HTMLButtonElement[] = [];

  const pick = (i: number) => {
    const o = options[i];
    if (!o || o.disabled) return;
    own?.set(o.value);
    props.onChange?.(o.value);
  };

  const segment = (o: SegmentedOption<T>, i: number) => {
    const reason = typeof o.disabled === 'string' ? o.disabled : undefined;
    const el = (
      <button
        type="button"
        role="radio"
        class={classes.segment}
        aria-label={iconOnly ? o.label : undefined}
        aria-disabled={o.disabled ? 'true' : undefined}
        data-tip={reason ?? (iconOnly ? o.label : undefined)}
        onClick={() => pick(i)}
      >
        {o.icon ? <Icon name={o.icon} /> : null}
        {iconOnly ? null : o.label}
      </button>
    ) as HTMLButtonElement;
    buttons.push(el);
    return el;
  };

  const well = (
    <div
      class={classes.well}
      role="radiogroup"
      aria-labelledby={labelId}
      aria-label={labelId ? undefined : props['aria-label']}
      data-icons={iconOnly ? '' : undefined}
      onKeyDown={(e: KeyboardEvent) => onKey(e)}
    >
      {options.map(segment)}
    </div>
  );

  const onKey = radioGroup({
    keys: options.map((o) => o.value),
    buttons,
    disabled: (i) => !!options[i].disabled,
    active: value,
    pick,
  });

  return (
    <div class={props.class ? `${classes.segmented} ${props.class}` : classes.segmented} data-size={props.size ?? 'md'}>
      {props.label ? (
        <div class={classes.head}>
          <span class={classes.label} id={labelId}>
            {props.label}
          </span>
          <span class={classes.value}>
            {derive(value, (v) => options.find((o) => Object.is(o.value, v))?.label ?? '')}
          </span>
        </div>
      ) : null}
      {well}
    </div>
  );
}
