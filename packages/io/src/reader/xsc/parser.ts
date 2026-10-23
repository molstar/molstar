/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import type { StringLike } from '@molstar/core/util/string-like';
import { ReaderResult as Result } from '../result.js';

// NAMD's extended system files: `.xsc` (one line, the last step's box) and `.xst` (a line per saved step). A header
// `#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z ...` names the columns; `o` is the box's center (NAMD's
// `cellOrigin`), Å. A file without the labels has NAMD's column order.

export interface XscEntry {
  /** The simulation step */
  readonly step: number;
  /** The box vectors, Å */
  readonly a: [number, number, number];
  readonly b: [number, number, number];
  readonly c: [number, number, number];
  /** The box's center, Å */
  readonly origin: [number, number, number];
}

export interface XscFile {
  /** In file order (by step) */
  readonly entries: readonly XscEntry[];
}

const DefaultLabels = ['step', 'a_x', 'a_y', 'a_z', 'b_x', 'b_y', 'b_z', 'c_x', 'c_y', 'c_z', 'o_x', 'o_y', 'o_z'];

function parseInternal(data: string): Result<XscFile> {
  let labels = DefaultLabels;
  const entries: XscEntry[] = [];
  for (const raw of data.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith('#$LABELS')) {
      labels = line.slice('#$LABELS'.length).trim().split(/\s+/);
      continue;
    }
    if (line.startsWith('#')) continue;
    const values = line.split(/\s+/).map(Number);
    const at = (name: string) => {
      const i = labels.indexOf(name);
      return i >= 0 && i < values.length ? values[i] : NaN;
    };
    const v = (p: string): [number, number, number] => [at(`${p}_x`), at(`${p}_y`), at(`${p}_z`)];
    const entry: XscEntry = { step: at('step'), a: v('a'), b: v('b'), c: v('c'), origin: v('o') };
    const numbers = [entry.step, ...entry.a, ...entry.b, ...entry.c];
    if (numbers.some((n) => !Number.isFinite(n))) {
      return Result.error(`not a NAMD extended system line: '${line.slice(0, 80)}'`);
    }
    // a file without an origin column: NAMD's default, the center at 0
    if (entry.origin.some((n) => !Number.isFinite(n))) entry.origin.fill(0);
    entries.push(entry);
  }
  if (!entries.length) return Result.error('no NAMD extended system lines');
  return Result.success({ entries });
}

export function parseXsc(data: StringLike) {
  // a small text file: as one string
  return Task.create<Result<XscFile>>('Parse NAMD XSC', async () => parseInternal(data.toString()));
}
