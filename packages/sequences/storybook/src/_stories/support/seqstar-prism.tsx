/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every Seq* story page in Mol* Prism's tokens and rules (`@molstar/prism-ui`'s), themed by the
// Storybook's toolbar (the host sets `data-msp-theme` on the document, as Mol*'s Storybook does), and the canvas the
// viewer draws on, whose color and theme may differ from the UI's (Prism's canvas layer: `ShellView`'s `canvasColor`).

import type { CSSProperties, ReactNode } from 'react';
import p from './seqstar-prism.module.css';
import { useSeqStyles } from './seqstar-styles.js';

/** A story page: its title, what it shows, then its blocks. Carries `.msp-root` (Prism's reset). */
export function PrismPage({
  title,
  caption,
  children,
}: {
  title: ReactNode;
  caption?: ReactNode;
  children: ReactNode;
}) {
  useSeqStyles('tokens', 'page');
  return (
    <div className={`msp-root ${p.page}`}>
      <header className={p.header}>
        <h1 className={p.title}>{title}</h1>
        {caption ? <div className={p.caption}>{caption}</div> : null}
      </header>
      {children}
    </div>
  );
}

/** Prism's themes; a canvas takes `light` or `dark` by its color. */
export type PrismTheme = 'light' | 'dark' | 'blue';

/** The relative luminance (WCAG, 0..1) of `#rrggbb` (Prism's model/canvas-theme.ts). */
function luminance(hex: string): number {
  const n = Number.parseInt(hex.slice(1), 16);
  const channel = (c: number) => {
    const v = c / 255;
    return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
  };
  return 0.2126 * channel((n >> 16) & 0xff) + 0.7152 * channel((n >> 8) & 0xff) + 0.0722 * channel(n & 0xff);
}

/** The chrome's theme over a canvas of this color: dark (light text) below 0.179, where both texts contrast equally. */
export function canvasThemeOf(hex: string): 'light' | 'dark' {
  return luminance(hex) < 0.179 ? 'dark' : 'light';
}

/** Canvas colors a story offers: the UI theme's base (none), or one of its own as Mol*'s canvas background is. */
export const CANVAS_COLORS = { 'UI theme': null, white: '#ffffff', black: '#000000', slate: '#1f2733' } as const;
export type CanvasColor = keyof typeof CANVAS_COLORS;

/**
 * The viewer's canvas: `color` sets its base and its theme by luminance; `theme` pins a theme with that theme's base;
 * neither follows the UI theme. The viewer's skin is derived inside it (seqstar-prism.module.css), so its drawing and
 * its chrome follow the canvas, not the page. The canvas reads the skin again when a theme attribute or class changes
 * on it or above it; a new `color` is an inline style, which it does not watch: a story that changes `color` remounts
 * the viewer (a `key`), as a host would call `refreshSkin`.
 */
export function PrismCanvas({ color, theme, children }: { color?: string; theme?: PrismTheme; children: ReactNode }) {
  useSeqStyles('tokens');
  const t = color ? canvasThemeOf(color) : theme;
  const style = color ? ({ '--msp-color-bg': color } as CSSProperties) : undefined;
  return (
    <div className={p.canvas} data-msp-theme={t} data-msp-box="" style={style}>
      {children}
    </div>
  );
}

/** The page's class names, for a story's blocks (`stack`, `pair`) and a card's controls (`controls`, `wide`, `actions`). */
export const prism = p;
