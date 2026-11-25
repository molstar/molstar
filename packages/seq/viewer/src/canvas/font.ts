/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultFonts, type FontStacks } from '@molstar/seqx-core/render';

// A draw list's font as a CSS font (_spec/README.md "Canvas2D backend"): the shorthand names a family role (`12px sans`,
// `bold 10px mono`), written with the token table's stack. Shared by the drawing and the text metrics, so a width
// measured is a width drawn.

/** The shorthand written as a CSS `font` value with the stack of its role (a family that names none is kept). */
export function cssFont(font: string, stacks: FontStacks = DefaultFonts): string {
  const size = /(\d+(?:\.\d+)?)px/.exec(font);
  const weight = /\b(bold|[1-9]00)\b/.exec(font);
  const family = /\bmono(space)?\b/i.test(font) ? stacks.mono : /\bsans\b/i.test(font) ? stacks.sans : undefined;
  if (!family) return font;
  return `${weight ? `${weight[1]} ` : ''}${size ? size[1] : '12'}px ${family}`;
}

/** `cssFont` kept per shorthand: drawing allocates none once a font has been seen. */
export class FontCache {
  private readonly _css = new Map<string, string>();

  constructor(public stacks: FontStacks = DefaultFonts) {}

  css(font: string): string {
    let c = this._css.get(font);
    if (c === undefined) this._css.set(font, (c = cssFont(font, this.stacks)));
    return c;
  }

  /** New stacks (a theme changed): the kept strings are written again. */
  setStacks(stacks: FontStacks): void {
    if (stacks === this.stacks) return;
    this.stacks = stacks;
    this._css.clear();
  }
}
