/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from '../../dom/index.js';
import type { MaybeReactive } from '../../dom/index.js';
import { classes } from './_generated/button.module.css.js';

export interface ButtonProps {
  variant?: 'default' | 'primary' | 'ghost' | 'danger';
  /**
   * `md` (default): a row's height (24px), a 16px icon; `sm`: a control's height (20px), a 12px icon, for a form row;
   * `xs`: 20px with small text, a 12px icon and a tighter padding, next to captions and counts.
   */
  size?: 'md' | 'sm' | 'xs';
  /** An icon before the label; without a label the button is square and needs a `title`. */
  icon?: Child;
  children?: Child;
  /** Tooltip (the app's, `data-tip`: no native `title`) and accessible name (required for icon-only buttons). */
  title?: MaybeReactive<string | undefined>;
  disabled?: MaybeReactive<boolean | undefined>;
  /** A toggle button's state (`aria-pressed`). */
  pressed?: MaybeReactive<boolean | undefined>;
  /** Stretches across its container (the apply button ending a form). */
  block?: boolean;
  class?: string;
  onClick?: (e: MouseEvent) => void;
}

export function Button(props: ButtonProps) {
  const variant = props.variant ?? 'default';
  const iconOnly = props.icon !== undefined && props.children === undefined;
  const className = [
    classes.button,
    variant !== 'default' ? classes[variant] : '',
    props.size === 'sm' ? classes.small : props.size === 'xs' ? classes.tiny : '',
    iconOnly ? classes.iconOnly : '',
    props.block ? classes.block : '',
    props.class ?? '',
  ]
    .filter(Boolean)
    .join(' ');
  return (
    <button
      type="button"
      class={className}
      data-tip={props.title}
      aria-label={iconOnly ? props.title : undefined}
      disabled={props.disabled}
      aria-pressed={props.pressed}
      onClick={props.onClick}
    >
      {props.icon !== undefined ? <span class={classes.icon}>{props.icon}</span> : null}
      {props.children}
    </button>
  );
}
