/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import {
  type Coverage,
  DefaultOpenLimits,
  type FormatRuntime,
  type GridDense,
  type MatrixItems,
  type OpenContext,
  type OpenedSource,
  type PairsList,
  type RegionRequest,
  type SourceEntry,
  withSummaries,
} from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { checkAdapter, checkRegionReads } from '../../data/_test/conformance.js';
import type { GridLay } from '../grid.js';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';

// A host's adapter serving the region shapes (_spec/README.md "Data", Protocol; "Seq io"): a host that holds its data in
// its own forms (links as objects, contacts as a function of two bins, a genotype table as arrays) and answers a
// region by cutting them, with no `seq` file and no asset read. The two consumers of one adapter in one test: the
// conformance suite (the protocol's checks and the region reads equal to the whole read clipped, at every level) and the
// renderer drawing a pairs layer, a grid and a matrix from the same sources (EC-208, EC-220, EC-223).

const Length = 2048;
const Systems = [{ name: 'chr', length: Length }];

/** The host's own data. */
interface Link {
  readonly a: number;
  readonly b: number;
  readonly id: string;
}
const links: Link[] = Array.from({ length: 300 }, (_, i) => ({
  a: (i * 37) % (Length - 20),
  b: (i * 91 + 500) % (Length - 20),
  id: `link${i}`,
}));
/** A symmetric contact count: falling with the distance. */
const contact = (i: number, j: number) => Math.floor(500 / (1 + Math.abs(i - j)));
const variants = Array.from({ length: 120 }, (_, i) => ({ start: 10 + i * 16, end: 11 + i * 16, id: `v${i}` }));
variants.push({ start: 600, end: 900, id: 'del' });
variants.sort((x, y) => x.start - y.start);
const States = ['s1', 's2', 's3'];
const genotype = (variant: string, state: number) =>
  variant === 'del' ? (state === 0 ? 1 : 0) : (variant.length + state) % 3;

/** The host's adapter: it answers a region from the data it holds, whichever shape the source `ref` names. */
const HostRuntime: FormatRuntime = {
  name: 'host-regions',
  formats: ['host-regions'],
  async open(source: SourceEntry, _ctx: OpenContext): Promise<OpenedSource> {
    const kind = source.ref;
    const coverage: Coverage =
      kind === 'pairs'
        ? { systems: Systems, shapes: [{ system: 'chr', shape: 'pairs' }], levels: [1] }
        : kind === 'grid'
          ? {
              systems: Systems,
              shapes: [{ system: 'chr', shape: 'grid', y: ['chr'], symmetric: true, aggregate: 'sum' }],
              levels: [1],
            }
          : { systems: Systems, shapes: [{ system: 'chr', shape: 'matrix', states: States }], levels: [1] };
    return {
      coverage,
      async read(request: RegionRequest, signal) {
        if (signal.aborted) throw Object.assign(new Error('aborted'), { name: 'AbortError' });
        const { start, end } = request;
        if (kind === 'pairs') {
          const rows = links.filter((l) => (l.a < end && l.a + 10 > start) || (l.b < end && l.b + 10 > start));
          const pairs: PairsList = {
            form: 'list',
            aStart: Int32Array.from(rows, (l) => l.a),
            aEnd: Int32Array.from(rows, (l) => l.a + 10),
            bStart: Int32Array.from(rows, (l) => l.b),
            bEnd: Int32Array.from(rows, (l) => l.b + 10),
            id: categoryColumn(rows.map((l) => l.id)),
          };
          return { ok: true, tables: { pairs } };
        }
        if (kind === 'grid') {
          const y = request.y ?? { system: 'chr', start, end };
          const nx = end - start;
          const ny = y.end - y.start;
          const value = new Float32Array(nx * ny);
          for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) value[i * ny + j] = contact(start + i, y.start + j);
          const grid: GridDense = { form: 'dense', x: start, y: y.start, nx, ny, value };
          return { ok: true, tables: { grid } };
        }
        // items with an end in the region, uncut
        const rows = variants.filter((v) => v.start < end && v.end > start);
        const codes = new Uint8Array(rows.length * States.length);
        rows.forEach((v, k) => States.forEach((_, s) => (codes[k * States.length + s] = genotype(v.id, s))));
        const matrix: MatrixItems = {
          form: 'items',
          start: Int32Array.from(rows, (v) => v.start),
          end: Int32Array.from(rows, (v) => v.end),
          states: States.length,
          codes,
          id: categoryColumn(rows.map((v) => v.id)),
        };
        return { ok: true, tables: { matrix } };
      },
      dispose() {},
    };
  },
};

const entry = (ref: string): SourceEntry => ({ address: `/data/${ref}`, format: 'host-regions', ref });
/** No asset is ever read: the host holds the data. */
const reads: string[] = [];
const assets = {
  readText: async (uri: string) => (reads.push(uri), ''),
  readBytes: async (uri: string) => (reads.push(uri), new Uint8Array()),
};

describe('a host’s adapter serving region shapes: the conformance suite and the viewer’s renderer', () => {
  it('pairs, grid tiles and a matrix pass the protocol’s checks and the region reads’ (EC-208, EC-220, EC-223)', async () => {
    for (const ref of ['pairs', 'grid', 'matrix']) {
      // `checkAdapter` reads every shape over a whole system with no second axis, which a grid is refused (`checkAnswer`:
      // "a grid without a second axis asked"): a grid's own checks are the region reads' below, which ask `y`
      if (ref !== 'grid') {
        const opened = await checkAdapter({ runtime: HostRuntime, entry: entry(ref), assets });
        expect(opened.coverage.shapes[0].shape).toBe(ref);
      }
      // the region reads equal the whole read clipped: at the host's level, and at the store's summaries' levels
      const live = await HostRuntime.open(entry(ref), {
        assets,
        signal: new AbortController().signal,
        limits: DefaultOpenLimits,
      });
      await checkRegionReads(live, 256);
      if (ref !== 'matrix') await checkRegionReads(withSummaries(live), 256);
    }
    expect(reads).toEqual([]);
  });

  it('draw: the same sources give arcs, a grid’s tiles and a matrix’s cells', async () => {
    const b = createSeqSpec({ title: 'a host’s region reads' });
    const src = (id: string) => b.data.source(id, { format: 'host-regions', ref: id } as never);
    const [pairsSource, gridSource, matrixSource] = ['pairs', 'grid', 'matrix'].map(src);
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: matrixSource, name: 'chr' });
    const pairs = b.data.pairs('links', { a: chr, source: pairsSource } as never);
    const contacts = b.data.grid('contacts', { x: chr, source: gridSource } as never);
    const gt = b.data.matrix('gt', { target: chr, source: matrixSource } as never);
    b.views.view('main', (v) => {
      const top = v.section('chr', { of: chr });
      top.track('links', { height: 80 }).layer('links', { type: 'pairs', data: pairs, legend: 'none' } as never);
      top.track('gt', { height: 60 }).layer('gt', { type: 'matrix', data: gt });
      v.section('plot', { layout: 'grid', of: chr, follow: top.address } as never).track('plot', {}, (t) => {
        t.layer('contacts', { type: 'grid', data: contacts });
      });
    });
    const s = await settle(b.build(), [HostRuntime]);
    const f = await s.frame({ size: { width: 720, height: 900 } });
    expect(f.layout.issues.filter((i) => i.severity === 'error')).toEqual([]);
    // the pairs: every link in view an arc
    const arcs = f.layout.sections[0].sectionLays[0].lay as unknown as { n: number };
    expect(arcs.n).toBeGreaterThan(100);
    const matrix = f.layout.sections[0].tracks.find((t) => t.track.label === 'gt')!.layers[0]
      .lay as unknown as MatrixLay;
    expect(matrix.states).toBe(States.length);
    expect(matrix.n).toBeGreaterThan(States.length * 50);
    const grid = f.layout.sections[1].tracks[0].layers[0].lay as unknown as GridLay;
    expect(grid.images.length).toBeGreaterThan(0);
    expect(grid.pend.length).toBe(0);
    // nothing came from a file
    expect(s.store.issues.value).toEqual([]);
  });
});
