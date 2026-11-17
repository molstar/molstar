/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from '../dom/index.js';
import { GlassSplash } from '../components/splash/splash.js';
import type { EmptyMarkContext } from './empty-stage.js';

// The Mol* glass as an empty stage's mark (stage §4), opt in: `Stage`'s `emptyMark={GlassEmptyMark}`. A module of its
// own, which the stage does not import: an app without it does not bundle the glass.

/** The Mol* wordmark as glass (`GlassSplash`), smaller on a compact stage. */
export const GlassEmptyMark = (ctx: EmptyMarkContext): Child => <GlassSplash scale={ctx.compact ? 2 : 3} />;
