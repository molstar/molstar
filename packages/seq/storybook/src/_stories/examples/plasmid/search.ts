/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RotatedOrigin } from '@molstar/seqx-examples/plasmid/derive';
import { PublishedLetters } from '@molstar/seqx-examples/plasmid/record';

// The Plasmid's search story's patterns (its story and its headless twin): letters across the rotated origin.

/**
 * Twelve letters of the rotated pUC19c across its origin (inside the primer site, plus strand): the published record's
 * letters on either side of the rotated origin.
 */
export const AcrossTheOrigin = PublishedLetters.slice(RotatedOrigin - 6, RotatedOrigin + 6);

const Complement: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A' };

/** The same site's reverse complement: found on the minus strand. */
export const AcrossTheOriginMinus = [...AcrossTheOrigin]
  .reverse()
  .map((c) => Complement[c] ?? c)
  .join('');
