/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Chain, Labels } from './numbering.js';

/** What the views hold of the recorded files: each chain's letters (1N8Z) and its IMGT labels (SAbDab's). */
export type Derived = Readonly<Record<Chain, { readonly letters: string; readonly imgt: Labels }>>;
