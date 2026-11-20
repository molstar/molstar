/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect } from 'vitest';
import type { FormatRuntime, OpenContext, SourceEntry } from '../protocol.js';
import { checkCoverage } from '../seq-format.js';
import { checkAnswer } from '../tables.js';

// The format runtimes' shared contract (sources.md, Consequences: one suite every runtime passes): a coverage that
// checks; answers that check over the whole of each system, an empty region, its last position and a region past its
// end, at every level the coverage lists; a system it does not have not served; an aborted read rejected as aborted.
// The records suite's (`records.test.ts`) with the levels: a runtime with levels of its own (bigWig) is read at each.

const signal = new AbortController().signal;

export async function sourceSuite(runtime: FormatRuntime, entry: SourceEntry, context: OpenContext) {
  const s = await runtime.open(entry, context);
  expect(checkCoverage(s.coverage)).toBeUndefined();
  for (const { system, shape } of s.coverage.shapes) {
    const length = s.coverage.systems.find((x) => x.name === system)?.length ?? 1024;
    for (const level of s.coverage.levels) {
      // the regions at a level start on its multiples; a tile of 1,024 bins at most
      const whole = Math.min(length, 1024 * level);
      const last = Math.max(0, Math.floor((length - 1) / level) * level);
      for (const [start, end] of [
        [0, whole],
        [0, 0],
        [last, length],
        [Math.ceil(length / level) * level, Math.ceil(length / level) * level + 7 * level],
      ]) {
        const request = { system, start, end, level, shapes: [shape] };
        const answer = await s.read(request, signal);
        expect(checkAnswer(answer, request), `${system} ${shape} ${level} ${start}-${end}`).toBeUndefined();
      }
    }
  }
  const unknown = await s.read({ system: '\u0000none', start: 0, end: 1, level: 1, shapes: ['ranges'] }, signal);
  expect(unknown).toMatchObject({ ok: false, kind: 'not-served' });
  const stop = new AbortController();
  stop.abort();
  const first = s.coverage.shapes[0];
  await expect(
    s.read({ system: first.system, start: 0, end: 1, level: 1, shapes: [first.shape] }, stop.signal),
  ).rejects.toMatchObject({ name: 'AbortError' });
  s.dispose();
}
