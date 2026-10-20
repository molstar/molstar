/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 14 and the genome story rendered headless to SVG (SV1b), in Node and with the committed tables' widths: well-formed,
// every number finite, the text where the layout put it, in a bounded time and size at two scales. The PNG of the same
// two, rasterized by Chrome as `seqv-render --png` does, with the pixels checked in _test/canvas/export.test.ts. Chrome,
// the times and the 3,000,000-nucleotide genome are the heavy lane's; the quick lane draws a tenth of the genome.

import { deflateSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { case14 } from '../../../core/_test/format/cases.js';
import { genomeDoc, synthGenome } from '../../../core/_test/format/genome-fixture.js';
import { itWithin } from '../../../core/_test/model/timing.js';
import { rasterizeSvg } from '../../src/cli/rasterize.js';
import { DefaultFontMetrics, metricsMeasurer, renderSvg } from '../../src/svg/index.js';
import { decodePng } from './png.js';
import { all, parseXml } from './xml.js';

const ports = { measurer: metricsMeasurer(DefaultFontMetrics), deflate: (b: Uint8Array) => deflateSync(b) };

describe.each([1, 2])('case 14 at scale %i', (scale) => {
  itWithin("a protein over its coding DNA: letters, the gene's tracks, the codons", async (within) => {
    // About 30 ms.
    const xml = await within(5000, () => renderSvg(case14(), { width: 1200, scale, overlay: true, nav: true, ports }));
    expect(xml.length).toBeLessThan(100_000); // about 15 KB
    expect(xml).not.toMatch(/NaN|Infinity|undefined/);
    const svg = parseXml(xml);
    expect(svg.attrs.width).toBe(String(1200 * scale));
    expect(svg.attrs.viewBox?.startsWith('0 0 1200 ')).toBe(true);
    const letters = all(svg, 'text').map((e) => e.text);
    for (const want of ['DNA', 'Protein', 'Exons', 'Domain', 'Phosphosite']) expect(letters, want).toContain(want);
    // Baselines are explicit offsets (the tables carry an em box), and nothing is left to a reader's dominant-baseline.
    expect(xml).not.toContain('dominant-baseline');
  });
});

describe('as PNG, through Chrome', { timeout: 120_000, tags: ['heavy'] }, () => {
  /** A PNG of the SVG, and of the same size with nothing drawn: Chrome's PNG filters, so only sizes are compared. */
  async function drawnAndBlank(svg: string) {
    const [w, h] = [/ width="([\d.]+)"/.exec(svg)?.[1], / height="([\d.]+)"/.exec(svg)?.[1]];
    const blank = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"></svg>`;
    const t = performance.now();
    const drawn = await rasterizeSvg(svg);
    return { drawn, blank: await rasterizeSvg(blank), ms: performance.now() - t, w: Number(w), h: Number(h) };
  }
  it("case 14 at scale 2: the SVG's size, ink on it, in seconds", async () => {
    const svg = await renderSvg(case14(), { width: 900, scale: 2, ports });
    const r = await drawnAndBlank(svg);
    const d = decodePng(r.drawn);
    expect([d.w, d.h]).toEqual([1800, Math.ceil(r.h)]);
    expect(r.drawn.length).toBeGreaterThan(r.blank.length * 3);
    expect(r.drawn.length).toBeLessThan(1_000_000);
    expect(r.ms).toBeLessThan(60_000);
  });
  it('the genome story', async () => {
    const svg = await renderSvg(genomeDoc(synthGenome()), { width: 1200, ports });
    const r = await drawnAndBlank(svg);
    expect(decodePng(r.drawn).w).toBe(1200);
    expect(r.drawn.length).toBeGreaterThan(r.blank.length * 3);
    expect(r.drawn.length).toBeLessThan(3_000_000);
    expect(r.ms).toBeLessThan(60_000);
  });
});

/** The genome story at two scales: every number finite, its tracks named, the larger under four times the smaller. */
async function genomeAtTwoScales(doc: ReturnType<typeof genomeDoc>) {
  const run = async (scale: number) => {
    const t = performance.now();
    const xml = await renderSvg(doc, { width: 1200, scale, overlay: true, nav: true, ports });
    return { xml, ms: performance.now() - t };
  };
  const [a, b] = [await run(1), await run(2)];
  for (const r of [a, b]) {
    expect(r.xml).not.toMatch(/NaN|Infinity|undefined/);
    expect(all(parseXml(r.xml), 'text').map((e) => e.text)).toContain('chr1');
    for (const want of ['Cytobands', 'Genes', 'GC content', 'Repeats'])
      expect(
        all(parseXml(r.xml), 'text').map((e) => e.text),
        want,
      ).toContain(want);
  }
  expect(b.xml.length).toBeLessThan(a.xml.length * 4);
  return [a, b];
}

describe('the genome story', () => {
  it('a tenth of it, 300,000 nucleotides: two scales', async () => {
    await genomeAtTwoScales(
      genomeDoc(synthGenome({ total: 300_000, genes: 300, proteins: 30, variants: 2000, repeats: 1000 })),
    );
  });

  it('24 chromosomes, 3,000,000 nucleotides: two scales in a bounded time and size', {
    tags: ['heavy'],
    timeout: 60_000,
  }, async () => {
    for (const r of await genomeAtTwoScales(genomeDoc(synthGenome()))) {
      expect(r.ms).toBeLessThan(30_000); // about 1.3 s
      expect(r.xml.length).toBeLessThan(2_000_000); // 0.2 to 0.45 MB
    }
  });
});
