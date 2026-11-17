/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GlassMark } from '../logo/glass-mark.js';
import { classes } from './_generated/splash.module.css.js';

// The splash (stage §4): the Mol* wordmark as glass, large, over an empty stage. The frame loop is `GlassMark`'s.

/** Room around the word for the pieces (units): a splash's pieces fly far (travel 1). */
const Margin = 64;
/** The cut: the splash's own (the heading's is another). */
const Seed = 7;

export interface SplashProps {
  /** Whole pixels a unit (default 3: 171 × 66 px). */
  scale?: number;
  class?: string;
}

/**
 * The glass wordmark, large (stage §4): each time it shows it assembles, its pieces drifting in from around it; under
 * the pointer they lean away, a press shatters it for a moment before it gathers again. No frames at rest, hidden or
 * off screen; one still frame under reduced motion.
 */
export function Splash(props: SplashProps = {}): Node {
  return (
    <GlassMark
      scale={props.scale ?? 3}
      margin={Margin}
      seed={Seed}
      travel={1}
      class={props.class ? `${classes.splash} ${props.class}` : classes.splash}
      data="data-msui-splash"
    />
  ) as Node;
}
