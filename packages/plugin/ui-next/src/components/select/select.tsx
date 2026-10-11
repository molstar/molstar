/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MaybeReactive } from '../../dom/scope.js';
import { Icon } from '../icon/icon.js';
import { classes } from './_generated/select.module.css.js';

export interface SelectProps {
  /** The options as [value, label]. */
  options: readonly (readonly [value: string, label: string])[];
  value?: MaybeReactive<string>;
  onChange?(value: string): void;
  size?: 'sm' | 'md';
  class?: string;
  title?: string;
  'aria-label'?: string;
}

/**
 * A dropdown that looks the same on every platform: the native `<select>` (its keyboard handling, accessibility and
 * option list) with the platform's own look removed, drawn with the theme's control color, radius, and a chevron.
 */
export function Select(props: SelectProps) {
  return (
    <span class={`${classes.select} ${props.class ?? ''}`} data-size={props.size ?? 'md'}>
      <select
        class={classes.native}
        value={props.value}
        title={props.title}
        aria-label={props['aria-label']}
        onChange={(e) => props.onChange?.((e.currentTarget as HTMLSelectElement).value)}
      >
        {props.options.map(([value, label]) => (
          <option value={value}>{label}</option>
        ))}
      </select>
      <Icon name="chevron-down" size={14} class={classes.chevron} />
    </span>
  );
}
