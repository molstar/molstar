/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GlassTiming } from './glass.js';
import { GlassMark } from './glass-mark.js';
import type { WordmarkLabels } from './wordmark.js';

/** How far the pieces travel, as a factor of a splash's: a heading has little room. */
const Travel = 0.35;
/** A heading's assembly: quicker than a splash's. */
const Quick: GlassTiming = { duration: 0.9, stagger: 0.3 };
/** Room around the word for the pieces (units). */
const Margin = 22;
/** The cut: the heading's own (the splash's is another). */
const Seed = 23;

/**
 * The Mol* wordmark as glass, small: a view's heading (the shell's Home). Each time it shows (the view opened, its group
 * unfolded) it assembles quickly, its pieces travelling little; the rest is `GlassMark`'s.
 */
export function GlassWordmark(props: {
  /** Whole pixels a unit (default 1: 57 × 22). */
  scale?: number;
  class?: string;
  /** Whether it assembles as it shows (asked each time it shows); otherwise it shows whole. Default: always. */
  assemble?: () => boolean;
  /** Its accessible name (`mark`: the word it draws). */
  labels?: Partial<WordmarkLabels>;
}): Node {
  return (
    <GlassMark
      scale={props.scale ?? 1}
      margin={Margin}
      seed={Seed}
      travel={Travel}
      timing={Quick}
      assemble={props.assemble}
      class={props.class}
      labels={props.labels}
    />
  ) as Node;
}
