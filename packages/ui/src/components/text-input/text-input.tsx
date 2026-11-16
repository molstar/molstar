/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type MaybeReactive, setAttr, watch } from '../../dom/index.js';
import { held } from '../../state.js';
import { classes } from './_generated/text-input.module.css.js';

export interface TextInputProps {
  value?: MaybeReactive<string>;
  onInput?(value: string): void;
  /** Enter pressed. */
  onSubmit?(value: string): void;
  placeholder?: MaybeReactive<string>;
  /** A secret (an API key): its characters hidden. */
  type?: 'text' | 'password';
  /** A search or filter field: Escape empties it (then, empty, lets Escape through to close what it is in). */
  search?: boolean;
  /**
   * As `Button`'s sizes: `md` (default) a row's height (24px), for a field on its own (a toolbar's filter); `sm` a
   * control (`--msui-control-height`, 20px), centered in a form row (styles §4).
   */
  size?: 'sm' | 'md';
  disabled?: MaybeReactive<boolean | undefined>;
  /** Tooltip (`data-tip`, the app's tooltips). */
  title?: string;
  class?: string;
  'aria-label'?: string;
}

/**
 * A one-line text field in the theme's control style (no border: the control color sets it apart), like `Select`; its
 * text starts at the control inset (styles §5). The reactive properties are written by one watch.
 */
export function TextInput(props: TextInputProps) {
  const value = held(props.value, '');
  const placeholder = held(props.placeholder, '');
  const disabled = held(props.disabled, false);
  const input = (
    <input
      type={props.type ?? 'text'}
      class={props.class ? `${classes.input} ${props.class}` : classes.input}
      data-size={props.size ?? 'md'}
      data-tip={props.title}
      aria-label={props['aria-label']}
      spellcheck={false}
      onInput={(e) => props.onInput?.((e.currentTarget as HTMLInputElement).value)}
      onKeyDown={(e) => {
        const el = e.currentTarget as HTMLInputElement;
        if (e.key === 'Enter') props.onSubmit?.(el.value);
        else if (e.key === 'Escape' && props.search && el.value) {
          el.value = '';
          props.onInput?.('');
          e.preventDefault();
          e.stopPropagation();
        }
      }}
    />
  ) as HTMLInputElement;
  watch([value, placeholder, disabled], () => {
    // Written only when it differs: an echo of what was typed leaves the caret where it is.
    const v = value.value ?? '';
    if (input.value !== v) input.value = v;
    setAttr(input, 'placeholder', placeholder.value || undefined);
    input.disabled = !!disabled.value;
  });
  return input;
}
