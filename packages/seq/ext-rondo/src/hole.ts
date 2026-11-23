/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LayoutContext, ResolvedSection } from '@molstar/seqx-core/extension';
import type { RondoLay, RondoPrep } from './lay.js';

// What an empty hole says (_spec/README.md, Hits, emphasis, words: The hole's words).

/**
 * A source's failure as the viewer whose ring it is said it (its store's issue's message), when a viewer is composed
 * with `rondoViewer` (it sets `of`), found by the resolved section laid out (the viewer's view's, its export's too);
 * headless: none, the hole says "Links failed".
 */
export const HoleIssue: { of?: (section: ResolvedSection, source: string) => string | undefined } = {};

/**
 * The hole's words: pending "Loading links…"; failed the issue its source's failure raised as the section's viewer holds
 * it (`HoleIssue`), else "Links failed". On a ring only (a row's band has no hole); measured once at layout.
 */
export function holeWords(lay: RondoLay, prep: RondoPrep, ctx: LayoutContext): void {
  if (!lay.frame.polar || (lay.empty !== 'pending' && lay.empty !== 'failed')) return;
  const source = prep.source.region?.src.source;
  lay.holeText =
    lay.empty === 'pending'
      ? 'Loading links…'
      : ((source !== undefined ? HoleIssue.of?.(ctx.section, source) : undefined) ?? 'Links failed');
  lay.holeW = ctx.text.measure(lay.font, lay.holeText);
}
