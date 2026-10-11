/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { bind } from '../../dom/scope.js';
import { Button } from '../button/button.js';
import { classes } from './_generated/color-picker.module.css.js';

/**
 * The palette, row by row (eight per row): neutrals, vivid hues, their light tints, their deep shades. Hues in the same
 * column across rows are one family, so a lighter or darker variant is a step up or down.
 */
export const ColorPalette: readonly (readonly number[])[] = [
  [0x000000, 0x262626, 0x4d4d4d, 0x808080, 0xb3b3b3, 0xd9d9d9, 0xf2f2f2, 0xffffff],
  [0xe5484d, 0xf76b15, 0xffc53d, 0x46a758, 0x12a594, 0x0090ff, 0x6e56cf, 0xd6409f],
  [0xf9c6c6, 0xffcca7, 0xffe3a2, 0xb2ddb5, 0xa5dcd3, 0xacd8fc, 0xd0c8f5, 0xf3c6e2],
  [0x8c1d18, 0xa33c00, 0x9e6c00, 0x2a7e3b, 0x067a6f, 0x0b5cad, 0x4c3a9e, 0xa1266f],
];

export function colorToHex(c: number): string {
  return `#${((c >>> 0) & 0xffffff).toString(16).padStart(6, '0')}`;
}

/** `#rrggbb`, `rrggbb` or `#rgb`; undefined when it is not a color. */
export function parseHex(text: string): number | undefined {
  let t = text.trim().replace(/^#/, '');
  if (/^[0-9a-f]{3}$/i.test(t)) t = t.replace(/./g, (d) => d + d);
  return /^[0-9a-f]{6}$/i.test(t) ? Number.parseInt(t, 16) : undefined;
}

/** Mixes a color toward white (`amount` > 0) or black (< 0), per channel. */
export function shade(c: number, amount: number): number {
  const target = amount > 0 ? 255 : 0;
  const t = Math.abs(amount);
  const ch = (shift: number) => {
    const v = (c >> shift) & 0xff;
    return Math.round(v + (target - v) * t) << shift;
  };
  return ch(16) | ch(8) | ch(0);
}

export interface ColorPickerProps {
  value: Observable<number>;
  onChange(color: number): void;
}

/**
 * A color picker that looks and works the same everywhere (no system picker): the palette as swatches, the current
 * color ringed; a hex field for any other color; lighten and darken in steps.
 */
export function ColorPicker(props: ColorPickerProps) {
  const value = props.value.pipe(distinctUntilChanged());
  let current = 0;
  bind(value, (v) => {
    current = v;
  });
  return (
    <div class={classes.picker}>
      <div class={classes.palette} role="listbox" aria-label="Palette">
        {ColorPalette.flat().map((c) => (
          <button
            type="button"
            role="option"
            class={classes.swatch}
            style={{ background: colorToHex(c) }}
            aria-label={colorToHex(c)}
            aria-selected={value.pipe(map((v) => v === c))}
            onClick={() => props.onChange(c)}
          />
        ))}
      </div>
      <div class={classes.custom}>
        <span class={classes.preview} style={{ background: value.pipe(map(colorToHex)) }} />
        <input
          type="text"
          class={classes.hex}
          aria-label="Hex color"
          spellcheck={false}
          value={value.pipe(map(colorToHex))}
          onInput={(e) => {
            const c = parseHex(e.currentTarget.value);
            if (c !== undefined) props.onChange(c);
          }}
          onBlur={(e) => {
            e.currentTarget.value = colorToHex(current);
          }}
        />
        <Button size="sm" variant="ghost" onClick={() => props.onChange(shade(current, 0.15))}>
          Lighten
        </Button>
        <Button size="sm" variant="ghost" onClick={() => props.onChange(shade(current, -0.15))}>
          Darken
        </Button>
      </div>
    </div>
  );
}
