/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Pyramid, PyramidLevel } from './pyramid.js';

// The level pair the rows draw their image tiles from (_spec/README.md "Render", The rows layer; performance.md,
// Levels): while the view moves, a level not built yet is built `rows.buildBudget` cells a frame and the level drawn last
// stays drawn meanwhile (its tiles placed at the new zoom), none when none of the base was drawn; a settled frame builds
// the rest.

/** The level pair each pyramid drew last: its base's version, row level and column level (kept, reused). */
const Drawn = new WeakMap<Pyramid, Int32Array>();

/**
 * The level to draw at `rowLevel` × `columnLevel`: that one when built (or settled: built now), else while the view moves
 * the one drawn last, or none (no level of this base drawn yet: a zoom out of the letters, a fold that changed the base)
 * until the wanted one is whole. Why none, not a whole build: a 10,000 × 300 level from its base is two frames of about
 * 100 ms; the rows' band shows no cells for the few frames its build takes, never a long frame.
 */
export function rowsLevel(
  pyr: Pyramid,
  rowLevel: number,
  columnLevel: number,
  settled: boolean,
  cells: number,
): PyramidLevel | undefined {
  let drawn = Drawn.get(pyr);
  if (!drawn) Drawn.set(pyr, (drawn = new Int32Array([-1, 0, 0])));
  let level = settled ? pyr.level(rowLevel, columnLevel) : pyr.levelWithin(rowLevel, columnLevel, cells);
  // the level drawn last, while it is the base's (its tiles are held), until the wanted one is whole
  if (!level && drawn[0] === pyr.version && pyr.isBuilt(drawn[1], drawn[2])) level = pyr.level(drawn[1], drawn[2]);
  if (!level) return undefined;
  drawn[0] = pyr.version;
  drawn[1] = level.rowLevel;
  drawn[2] = level.columnLevel;
  return level;
}
