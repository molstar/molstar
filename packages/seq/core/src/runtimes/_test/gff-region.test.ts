/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readTabixIndex, regionBytes } from '@molstar/io/common/tabix';
import { describe, expect, it } from 'vitest';
import { sourceSuite } from '../../data/_test/source-suite.js';
import {
  DefaultOpenLimits,
  type OpenContext,
  type RangesList,
  type RegionAnswer,
  type SourceEntry,
} from '../../data/protocol.js';
import { checkAnswer } from '../../data/tables.js';
import { GffRuntime } from '../gff.js';
import { densityOf } from '../internal/levels.js';
import { memoryAssets, plainRows } from './bed-fixtures.js';
import { generateGff3, generateGtf, indexedGff } from './gff-fixtures.js';

// GFF3 and GTF by region (_spec/README.md "Data", GFF3 and GTF: by region): each tile's lines assembled alone, by type,
// from the hull of the lines meeting it; every tiling of a generated file equal to the whole read clipped, at every
// level and in every `features` mode (EC-227); a tile inside an intron gives the transcript's span row (EC-229); the
// levels and a tile's read cost from the index.

const signal = new AbortController().signal;
const names = ['chr1', 'chr2'];

const context = (files: Record<string, Uint8Array | string>): OpenContext => ({
  assets: memoryAssets(files),
  signal,
  limits: DefaultOpenLimits,
});
const entry = (format: string, more: Partial<SourceEntry> = {}): SourceEntry => ({
  address: '/data/g',
  format,
  uri: 'f',
  ...more,
});

/**
 * An answer's rows without the attributes they lack: a tile's fields are the fixed ones and the attributes its items
 * have, the whole read's every attribute of the file (null where an item has none).
 */
function rowsOf(a: RegionAnswer) {
  return (plainRows(a) as Record<string, unknown>[])
    .map((r) => Object.fromEntries(Object.entries(r).filter(([, v]) => v !== null && v !== undefined)))
    .sort(byStart);
}

/** Rows in order of start, then end and id: the order among rows of one start is not the protocol's. */
function byStart(a: Record<string, unknown>, b: Record<string, unknown>) {
  return (
    (a.start as number) - (b.start as number) ||
    (a.end as number) - (b.end as number) ||
    String(a.id).localeCompare(String(b.id))
  );
}

/** A tiling's rows deduped as a reader does: by id and start. */
function dedupe(rows: Record<string, unknown>[]) {
  const seen = new Set<string>();
  return rows.filter((r) => {
    const k = `${r.id}|${r.start}|${r.end}`;
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
}

describe('GFF3 and GTF by region', () => {
  const gff3 = generateGff3(24, names, 3, 200_000);
  const gtf = generateGtf(24, names, 4, true, 200_000);

  it('passes the source suite at every level (GFF3, GTF), and lists the levels of its index', async () => {
    for (const [format, g] of [
      ['gff3', gff3],
      ['gtf', gtf],
    ] as const) {
      const b = indexedGff(g.text, names);
      await sourceSuite(GffRuntime, entry(format, { index: 'i' }), context({ f: b.file, i: b.index }));
      const s = await GffRuntime.open(entry(format, { index: 'i' }), context({ f: b.file, i: b.index }));
      expect(s.coverage.levels).toEqual([1, 16, 256, 4096]);
      // `##sequence-region` lengths (GFF3's), the index's names (GTF's)
      expect(s.coverage.systems).toEqual(
        format === 'gff3' ? names.map((name) => ({ name, length: 250_000 })) : names.map((name) => ({ name })),
      );
      s.dispose();
    }
  });

  for (const [format, g] of [
    ['gff3', generateGff3(12, names, 5, 200_000)],
    ['gtf', generateGtf(12, names, 6, true, 200_000)],
  ] as const)
    for (const features of ['transcripts', 'genes', 'rows'] as const)
      it(`${format}, features ${features}: every tiling equals the whole read clipped, at every level (EC-227)`, async () => {
        const b = indexedGff(g.text, names);
        const more = { features } as Partial<SourceEntry>;
        const whole = await GffRuntime.open(entry(format, more), context({ f: g.text }));
        const region = await GffRuntime.open(
          entry(format, { index: 'i', ...more }),
          context({ f: b.file, i: b.index }),
        );
        for (const level of [1, 16, 256]) {
          const width = level * 1024 * (level === 1 ? 4 : 1);
          for (const system of names) {
            const tiles: Record<string, unknown>[] = [];
            for (let s = 0; s < 210_000; s += width) {
              const t = {
                system,
                start: s,
                end: s + width,
                level: level === 1 ? 1 : level,
                shapes: ['ranges' as const],
              };
              const items = await region.read({ ...t, items: true }, signal);
              expect(checkAnswer(items, { ...t, items: true })).toBeUndefined();
              const w = await whole.read({ ...t, level: 1 }, signal);
              expect(rowsOf(items), `${level} ${system}:${s}`).toEqual(rowsOf(w));
              tiles.push(...rowsOf(items));
              if (level > 1) {
                const d = await region.read(t, signal);
                expect(checkAnswer(d, t)).toBeUndefined();
                if (!d.ok || !w.ok) throw new Error('not read');
                expect(d.tables.ranges).toEqual(densityOf(w.tables.ranges as RangesList, t));
              }
            }
            // the tiles together, deduped by id and start: the whole read's pieces (a span row is a tile's own)
            const all = rowsOf(
              await whole.read({ system, start: 0, end: 210_000, level: 1, shapes: ['ranges'] }, signal),
            );
            const pieces = (rows: Record<string, unknown>[]) => dedupe(rows.filter((r) => r.start !== r.end));
            const order = (a: Record<string, unknown>, b: Record<string, unknown>) =>
              (a.start as number) - (b.start as number) || String(a.id).localeCompare(String(b.id));
            expect(pieces(tiles).sort(order)).toEqual(pieces(all).sort(order));
          }
        }
        whole.dispose();
        region.dispose();
      });

  it('a tile inside an intron gives the transcript’s span row, its thick part and phase the whole’s (EC-229)', async () => {
    const text = [
      '##gff-version 3',
      '##sequence-region chr1 1 100000',
      'chr1\tt\tgene\t1001\t40000\t.\t+\t.\tID=g1;Name=G1',
      'chr1\tt\tmRNA\t1001\t40000\t.\t+\t.\tID=t1;Parent=g1;Name=T1',
      'chr1\tt\texon\t1001\t1200\t.\t+\t.\tParent=t1',
      'chr1\tt\tCDS\t1101\t1200\t.\t+\t0\tParent=t1',
      'chr1\tt\texon\t39001\t40000\t.\t+\t.\tParent=t1',
      'chr1\tt\tCDS\t39001\t39500\t.\t+\t2\tParent=t1',
    ].join('\n');
    const b = indexedGff(`${text}\n`, ['chr1']);
    const s = await GffRuntime.open(entry('gff3', { index: 'i' }), context({ f: b.file, i: b.index }));
    const a = await s.read({ system: 'chr1', start: 20_480, end: 21_504, level: 1, shapes: ['ranges'] }, signal);
    expect(rowsOf(a)).toMatchObject([
      {
        id: 't1',
        start: 20_480,
        end: 20_480,
        itemStart: 1000,
        itemEnd: 40_000,
        label: 'T1',
        thick_start: 1100,
        thick_end: 39_500,
        cds_phase: 0,
        gene: 'g1',
      },
    ]);
    // the hull read: the tile read the transcript's lines (its CDS hull) though none of them is in it
    const ix = await readTabixIndex(b.index);
    expect(s.readCost!({ system: 'chr1', start: 20_480, end: 21_504, level: 1, shapes: ['ranges'] })).toBe(
      regionBytes(ix, 'chr1', 20_480, 21_504),
    );
  });
});
