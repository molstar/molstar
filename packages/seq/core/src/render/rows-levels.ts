/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Pyramid, PyramidLevel } from './pyramid.js';

// The level pair the rows draw their image tiles from (_spec/README.md "Render", The rows layer; performance.md,
// Levels): while the view moves, a level not built yet is built `rows.buildBudget` cells a frame and the level drawn last
// stays drawn meanwhile (its tiles placed at the new zoom); a settled frame builds the rest.

/** The level pair each pyramid drew last: its base's version, row level and column level (kept, reused). */
const Drawn = new WeakMap<Pyramid, Int32Array>();

/** The level to draw at `rowLevel` × `columnLevel`: that one when built (or settled), else the one drawn last. */
export function rowsLevel(
  pyr: Pyramid,
  rowLevel: number,
  columnLevel: number,
  settled: boolean,
  cells: number,
): PyramidLevel {
  let drawn = Drawn.get(pyr);
  if (!drawn) Drawn.set(pyr, (drawn = new Int32Array([-1, 0, 0])));
  let level = settled ? undefined : pyr.levelWithin(rowLevel, columnLevel, cells);
  // the level drawn last, while it is the base's (its tiles are held), until the wanted one is whole; none drawn yet
  // (a first frame): built now
  if (!level && !settled && drawn[0] === pyr.version && pyr.isBuilt(drawn[1], drawn[2]))
    level = pyr.level(drawn[1], drawn[2]);
  level ??= pyr.level(rowLevel, columnLevel);
  drawn[0] = pyr.version;
  drawn[1] = level.rowLevel;
  drawn[2] = level.columnLevel;
  return level;
}
