/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/icon.module.css.js';
import { LucideIcons } from './_generated/lucide.js';
import { MolecularIcons } from './molecular.js';

export type IconName = keyof typeof LucideIcons | keyof typeof MolecularIcons;

export interface IconProps {
  name: IconName;
  /**
   * Width and height in pixels, overriding the context's size. Without it the icon takes `--msui-icon-size` from where
   * it sits (a button, a slot, a header), else `--msui-icon` (16px).
   */
  size?: number;
  /** An accessible name; without it the icon is decorative (`aria-hidden`). */
  label?: string;
  class?: string;
}

const SVG_NS = 'http://www.w3.org/2000/svg';

const markups: Readonly<Record<string, string>> = { ...LucideIcons, ...MolecularIcons };
const templates = new Map<IconName, SVGSVGElement>();

/** The icon's parsed `<svg>`, created once per name; `Icon` clones it. */
function template(name: IconName): SVGSVGElement {
  let svg = templates.get(name);
  if (svg) return svg;
  const markup = markups[name];
  if (markup === undefined) throw new Error(`Unknown icon "${name}"`);
  svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('viewBox', '0 0 24 24');
  svg.setAttribute('fill', 'none');
  svg.setAttribute('stroke', 'currentColor');
  svg.setAttribute('stroke-width', '2');
  svg.setAttribute('stroke-linecap', 'round');
  svg.setAttribute('stroke-linejoin', 'round');
  svg.setAttribute('focusable', 'false');
  svg.innerHTML = markup;
  templates.set(name, svg);
  return svg;
}

/** A stroke icon (Lucide or a molecular glyph) in `currentColor`, sized by its context (styles §4.1) or `size`. */
export function Icon(props: IconProps): SVGSVGElement {
  const svg = template(props.name).cloneNode(true) as SVGSVGElement;
  // Sized by where it sits (styles §4.1: `--msui-icon-size`, set by a button, a slot, a header); `size` overrides it.
  if (props.size !== undefined) {
    svg.style.width = `${props.size}px`;
    svg.style.height = `${props.size}px`;
  }
  svg.setAttribute('class', props.class ? `${classes.icon} ${props.class}` : classes.icon);
  if (props.label !== undefined) {
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', props.label);
  } else {
    svg.setAttribute('aria-hidden', 'true');
  }
  return svg;
}
