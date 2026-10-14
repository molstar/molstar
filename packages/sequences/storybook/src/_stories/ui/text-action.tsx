/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { clsx } from 'clsx';
import type { ButtonHTMLAttributes, MouseEvent, ReactNode } from 'react';
import styles from './text-action.module.css';
import { Tooltip } from './tooltip.js';

export interface TextActionProps extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, 'className' | 'title'> {
  /** `default` for most actions, `strong` for the one that matters, `quiet` for meta actions like "log". */
  tone?: 'default' | 'strong' | 'quiet' | 'accent';
  className?: string;
  /** Its tooltip (never the DOM's `title`). Disabled, it stays hoverable, so it still says why. */
  title?: ReactNode;
}

/** The workhorse action: underlined text, not a box. */
export function TextAction({
  tone = 'default',
  className,
  type = 'button',
  title,
  disabled,
  onClick,
  ...props
}: TextActionProps) {
  const button = (
    <button
      type={type}
      className={clsx(styles.root, styles[tone], className)}
      // Disabled but still focusable and hoverable, so its tooltip can say why.
      aria-disabled={disabled || undefined}
      data-disabled={disabled || undefined}
      onClick={disabled ? (e: MouseEvent) => e.preventDefault() : onClick}
      {...props}
    />
  );
  return title ? <Tooltip content={title}>{button}</Tooltip> : button;
}
