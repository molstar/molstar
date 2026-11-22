/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import type { BlockPlacement } from '../layout/blocks.js';
import type { Frame } from '../renderer.js';
import { hitView } from '../tracks-hit.js';
import { settle } from './helpers.js';

// Data over an assembly read by sequence (_spec/README.md "Resolve", Assemblies: data over an assembly): each block of
// a section of the assembly reads its own sequence's rows, inline (split by `target`) or from a source by region.

beforeAll(async () => {
  await loadGenome();
});

const Toy = [
  { name: 'chr1', length: 3000 },
  { name: 'chr2', length: 2000 },
  { name: 'chr3', length: 1000 },
];

const Size = { width: 900, height: 200 };

/** The keys hit along each strip's band `band`, by strip (its sequence's name). */
function keysByStrip(f: Frame, band: number): Record<string, string[]> {
  const sec = f.layout.sections[0];
  const p = sec.placement as BlockPlacement;
  const s = p.strips;
  const out: Record<string, string[]> = {};
  for (let i = 0; i < s.n; i++) {
    const a = s.axes[i];
    const keys = new Set<string>();
    const y0 = s.bandY[i * s.bands + band];
    for (let x = a.x0; x < a.x0 + a.width; x += 1)
      for (let y = y0; y < y0 + s.bandH[i * s.bands + band]; y += 3) {
        const h = hitView(f.layout, x, y);
        if (h?.element && h.layer.endsWith('/bed')) keys.add(String(h.element.key));
      }
    out[a.system.replace(/^.*:name /, '').replace(/\)$/, '')] = [...keys].sort();
  }
  return out;
}

/** A toy assembly's section with a `bed` ranges layer over `rows` (inline) or `source` (a seq source's value). */
function bedSpec(over: { rows?: Record<string, unknown>; source?: SeqData; include?: string[] }): Node {
  const b = createSeqSpec({ title: 'bed over an assembly' });
  const g = b.data.assembly('toy', assembly('Toy', Toy, { order: 'source' } as never));
  const src = over.source && b.data.source('src', { format: 'seq', value: seqToValue(over.source) as never });
  const bed = b.data.ranges('bed', { target: g, ...(over.rows ? { rows: over.rows } : { source: src }) } as never);
  b.views.view('v', (v) =>
    genome(v, 'genome', g, over.include ? { include: over.include } : {}, (s) =>
      s.track('bed').layer('bed', { type: 'ranges', data: bed }),
    ),
  );
  return b.build();
}

describe('ranges over an assembly', () => {
  it('inline rows are split by their `target` through the names rule, each drawn on its block (EC-189)', async () => {
    const s = await settle(
      bedSpec({
        rows: {
          start: [100, 500, 200, 10, 0],
          end: [400, 900, 600, 300, 10],
          id: ['a', 'b', 'c', 'd', 'x'],
          target: ['chr1', '1', '2', 'chr3', 'chr99'],
        },
      }),
    );
    const f = await s.frame({ size: Size });
    expect(keysByStrip(f, 1)).toEqual({ chr1: ['a', 'b'], chr2: ['c'], chr3: ['d'] });
    // a row naming no sequence is counted, once
    const said = s.model.issues.filter((i) => i.code === 'seq.resolve.assembly.name');
    expect(said.map((i) => i.message)).toEqual(['Rows name 1 sequences the assembly lacks (chr99)']);
    expect(s.model.issues.filter((i) => i.code === 'seq.resolve.param.unread')).toEqual([]);
  });

  it('a source is read by region on each block, its names mapped both ways (`1` is `chr1`, EC-187)', async () => {
    const source: SeqData = {
      coverage: {
        systems: [
          { name: '1', length: 3000 },
          { name: '3', length: 1000 },
          { name: 'chrUn', length: 50 },
        ],
        shapes: ['1', '3', 'chrUn'].map((system) => ({ system, shape: 'ranges' as const })),
        levels: [1],
      },
      pairs: [
        { system: '1', end: 3000, id: ['p', 'q'], start: [10, 1000], stop: [200, 2000] },
        { system: '3', end: 1000, id: ['r'], start: [500], stop: [900] },
        { system: 'chrUn', end: 50, id: ['u'], start: [0], stop: [10] },
      ].map((x) => ({
        request: { system: x.system, start: 0, end: x.end, level: 1, shapes: ['ranges' as const] },
        answer: {
          ok: true as const,
          tables: {
            ranges: {
              form: 'list' as const,
              start: Int32Array.from(x.start),
              end: Int32Array.from(x.stop),
              id: categoryColumn(x.id),
            },
          },
        },
      })),
    };
    const s = await settle(bedSpec({ source }), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    expect(keysByStrip(f, 1)).toEqual({ chr1: ['p', 'q'], chr2: [], chr3: ['r'] });
    const said = s.model.issues.filter((i) => i.code === 'seq.resolve.assembly.name');
    expect(said.map((i) => i.message)).toEqual(['The source names 1 sequences the assembly lacks (chrUn)']);
  });

  it('a sequence the section’s `include` leaves out reads nothing', async () => {
    const s = await settle(
      bedSpec({
        rows: { start: [100, 200], end: [400, 600], id: ['a', 'c'], target: ['chr1', 'chr2'] },
        include: ['chr2', 'chr3'],
      }),
    );
    const f = await s.frame({ size: Size });
    expect(keysByStrip(f, 1)).toEqual({ chr2: ['c'], chr3: [] });
  });
});
