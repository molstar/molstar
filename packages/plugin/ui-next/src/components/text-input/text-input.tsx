/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MaybeReactive } from '../../dom/scope.js';
import { classes } from './_generated/text-input.module.css.js';

export interface TextInputProps {
  value?: MaybeReactive<string>;
  onInput?(value: string): void;
  /** Enter pressed. */
  onSubmit?(value: string): void;
  placeholder?: string;
  size?: 'sm' | 'md';
  class?: string;
  'aria-label'?: string;
}

/** A one-line text field in the theme's control style (no border: the control color sets it apart), like `Select`. */
export function TextInput(props: TextInputProps) {
  return (
    <input
      type="text"
      class={`${classes.input} ${props.class ?? ''}`}
      data-size={props.size ?? 'md'}
      value={props.value}
      placeholder={props.placeholder}
      aria-label={props['aria-label']}
      onInput={(e) => props.onInput?.((e.currentTarget as HTMLInputElement).value)}
      onKeyDown={(e) => {
        if (e.key === 'Enter') props.onSubmit?.((e.currentTarget as HTMLInputElement).value);
      }}
    />
  );
}
