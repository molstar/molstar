/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { abortError, type Coverage } from '../data/protocol.js';
import { indexPairs, PairsBuilder, pairsIn, pairsList, pairsSchema, rowBytes } from './bed-pairs.js';
import type { BedLine, DeclaredType } from './bed-parse.js';
import {
  type BedpeSource,
  dataCells,
  dataLines,
  failure,
  type Opened,
  refuse,
  refuseField,
  wholeSource,
} from './bed-read.js';

// A BEDPE read whole and indexed by both ends (_spec/README.md "Data", BED family): loaded when the first BEDPE source
// opens, so a host reading BED alone never loads the pairs index.

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
  const coverage: Coverage = { systems: rows.names.map((name) => ({ name })), shapes: [], levels: [1] };
  return {
    ...wholeSource(coverage, causes.issues(source.address), () => failure('not-served', 'Pairs only (`pairs`)')),
    summaryBytes: bytes,
    pairs(system, start, end, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      const s = rows.names.indexOf(system);
      if (s < 0) return Promise.resolve(failure('not-served', `No sequence ${system}`));
      const found = pairsIn(rows, index, s, start, end);
      return Promise.resolve(
        found.length > ctx.limits.maxRows
          ? failure('too-large', `${found.length} pairs in one region`)
          : { ok: true, pairs: pairsList(rows, found) },
      );
    },
  };
}
