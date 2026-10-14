/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stream manager as the viewer builds it, for tests, smoke runs and demos that drive it by hand (a pass, the
// dynamic nodes held): the root entry gives hosts its contract only (`StreamManager`, VT-10).

import type { SeqViewerContext } from '../context.js';
import type { Streams } from '../streams.js';

export type { Streams } from '../streams.js';
/** The context's stream manager, with its pass and held nodes. */
export const streamsOf = (ctx: SeqViewerContext): Streams => ctx.streams as Streams;
