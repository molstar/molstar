/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `renderSvg`: waits for the data its view needs and renders it; says when it did not settle; reads no assets unless
// given a reader; reports a missing view; the render codes' catalog.

import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { p53Spec } from '../../lang/_test/p53.js';
import { RenderIssues, renderSvg } from '../index.js';
import { blocks } from './helpers.js';
import { BigRuntime, bigSpec, p53TracksSpec } from './fixtures.js';

const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

describe('renderSvg', () => {
  it('waits for the p53 example’s data and renders it', async () => {
    const r = await renderSvg(p53TracksSpec(), { size: { width: 640, height: 120 }, viewport: { start: 0, end: 30 } });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    const [letters, items, values] = blocks(r.lists);
    expect(letters.letters).toBe(30);
    expect(items.rects).toBeGreaterThan(0);
    expect(values.rects).toBe(30);
    expect(r.svg).toContain('data-seq-settled="true"');
    expect(r.svg).toMatch(
      /data-seq-domains="\/data\/conservation|data-seq-domains="\/views\/main\/protein\/conservation:0=0,/,
    );
  });

  it('waits for a source read by region (a length its coverage gives)', async () => {
    const r = await renderSvg(bigSpec(), {
      size: { width: 640, height: 120 },
      runtimes: [BigRuntime],
      viewport: { start: 1000, end: 1030 },
    });
    expect(r.settled).toBe(true);
    expect(blocks(r.lists)[0].letters).toBe(30);
  });

  it('says it did not settle when a source never answers, drawing what is pending as pending', async () => {
    const silent: FormatRuntime = {
      name: 'silent',
      formats: ['synthetic'],
      open: async () => ({
        coverage: {
          systems: [{ name: 'chr', length: 1000 }],
          shapes: [{ system: 'chr', shape: 'ranges' }],
          levels: [1],
        },
        read: () => new Promise(() => {}),
        dispose() {},
      }),
    };
    const r = await renderSvg(bigSpec(), { size: { width: 640, height: 120 }, runtimes: [silent], timeoutMs: 0 });
    expect(r.settled).toBe(false);
    expect(codes(r.issues)).toContain('seq.render.unsettled');
    expect(r.svg).toContain('data-seq-settled="false"');
    // the items track is one placeholder over the view
    expect(blocks(r.lists)[1].rects).toBe(1);
  });

  it('reads no assets by default: a source by location is not opened (the gene still follows the protein)', async () => {
    const reads: string[] = [];
    const fasta: FormatRuntime = {
      name: 'fasta',
      formats: ['fasta'],
      open: async (_source, ctx) => {
        await ctx.assets.readText('data/chr17-region.fa', ctx.signal).then((t) => reads.push(t));
        throw new Error('unreachable');
      },
    };
    const r = await renderSvg(p53Spec(), { size: { width: 640, height: 120 }, runtimes: [fasta] });
    expect(reads).toEqual([]);
    expect(codes(r.issues)).toEqual(['seq.render.source']);
    expect(r.settled).toBe(true);
  });

  it('reports a view the spec does not have', async () => {
    const r = await renderSvg(p53TracksSpec(), { size: { width: 100, height: 50 }, view: '/views/none' });
    expect(codes(r.issues)).toEqual(['seq.render.view']);
    expect(r.svg).toContain('width="100" height="50"');
  });
});

describe('seq.render: issue codes', () => {
  it('every code a source uses is cataloged, and every cataloged code is used', () => {
    const known = new Set(Object.keys(RenderIssues.codes));
    const dir = fileURLToPath(new URL('../', import.meta.url));
    const used = new Set<string>();
    for (const name of readdirSync(dir)) {
      if (!name.endsWith('.ts') || name === 'issues.ts') continue;
      for (const m of readFileSync(join(dir, name), 'utf8').matchAll(/'(seq\.render\.[a-z_.]*[a-z_])'/g))
        used.add(m[1]);
    }
    expect([...used].filter((c) => !known.has(c))).toEqual([]);
    expect([...known].filter((c) => !used.has(c))).toEqual([]);
  });
});
