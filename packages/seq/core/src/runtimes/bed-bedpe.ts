/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { abortError, type Coverage } from '../data/protocol.js';
import { wholeLevel } from '../data/seq-format.js';
import { binPairs } from '../data/summaries-pairs.js';
import { indexPairs, PairsBuilder, pairsIn, pairsList, pairsSchema, rowBytes } from './bed-pairs.js';
import type { BedLine, DeclaredType } from './bed-parse.js';
import { type BedpeSource, dataCells, dataLines, type Opened, refuseField } from './bed-read.js';
import { failure, refuse, wholeSource } from './internal/answers.js';

// A BEDPE read whole and indexed by both ends (_spec/README.md "Data", BED family): served as the protocol's `pairs`
// shape by region, every pair with an end in it (EC-208), the items at any level when asked, and coarser levels as
// bin pairs made per tile from the held rows (the pairs level). Loaded when the first BEDPE source opens, so a host
// reading BED alone never loads the pairs index.

export function openPairs(o: Opened, text: BedLine[]): BedpeSource {
  const { source, ctx, causes, options } = o;
  const { columns, count } = dataLines(text);
  const schema = pairsSchema(
    columns || 6,
    source.fields,
    source.types as Readonly<Record<string, DeclaredType>> | undefined,
  );
  if ('error' in schema) refuseField(source.address, schema.error);
  const check = (bytes: number) => {
    if (bytes > options.maxPairsBytes)
      refuse(
        source.address,
        'seq.source.too_large',
        `The BEDPE source holds ${bytes} bytes with its two orders, over ${options.maxPairsBytes}`,
        'A format indexed by both ends (pairix), or a levelled pairs format',
      );
  };
  // the floor of what the lines will hold, refused before any is parsed
  check(count * rowBytes(schema, 1));
  const builder = new PairsBuilder(schema);
  for (const { cells, place } of dataCells(text, causes)) builder.add(cells, place, causes);
  const bytes = builder.rows * rowBytes(schema, builder.names.length);
  check(bytes);
  const rows = builder.finish();
  const index = indexPairs(rows);
  // levels of factor 16 up to the one whose tile holds the longest end's system (`wholeLevel`), bins made per tile
  let extent = 1;
  for (const c of [rows.aEnd, rows.bEnd]) for (let i = 0; i < c.length; i++) if (c[i] > extent) extent = c[i];
  const levels = [1];
  for (let l = 16, top = wholeLevel(extent); l <= top; l *= 16) levels.push(l);
  const coverage: Coverage = {
    systems: rows.names.map((name) => ({ name })),
    shapes: rows.names.map((system) => ({ system, shape: 'pairs' as const })),
    levels,
  };
  const found = (system: string, start: number, end: number): Int32Array | undefined => {
    const s = rows.names.indexOf(system);
    return s < 0 ? undefined : pairsIn(rows, index, s, start, end);
  };
  const maxRows = ctx.limits.maxRows;
  return {
    ...wholeSource(coverage, causes.issues(source.address), (r) => {
      if (!r.shapes.includes('pairs')) return failure('not-served', 'A BEDPE serves pairs only');
      const got = found(r.system, r.start, r.end)!;
      if (got.length > maxRows) return failure('too-large', `${got.length} pairs in one region`);
      if (r.level === 1 || r.items) return { ok: true, tables: { pairs: pairsList(rows, got) } };
      // the pairs level: bins of the pairs with an end's first position in the tile (each bin pair's every pair)
      const s = rows.names.indexOf(r.system);
      const anchored = got.filter(
        (k) =>
          (rows.aSys[k] === s && rows.aStart[k] >= r.start && rows.aStart[k] < r.end) ||
          (rows.bSys[k] === s && rows.bStart[k] >= r.start && rows.bStart[k] < r.end),
      );
      return { ok: true, tables: { pairs: binPairs(pairsList(rows, anchored), r.level) } };
    }),
    summaryBytes: bytes,
    pairs(system, start, end, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      const got = found(system, start, end);
      if (!got) return Promise.resolve(failure('not-served', `No sequence ${system}`));
      return Promise.resolve(
        got.length > maxRows
          ? failure('too-large', `${got.length} pairs in one region`)
          : { ok: true, pairs: pairsList(rows, got) },
      );
    },
  };
}
