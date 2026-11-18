/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, symlinkSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createSeqContext } from '@molstar/seqx-core/context';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import { GeneFile, geneSeqText } from '../p53/data.js';
import { RegionFile, regionSeqText, Transcripts } from '../p53/region.js';
import { type ExampleView, ExampleViews } from '../views.js';

// Every view: built, validated, resolved, laid out with no unexpected issue, rendered to SVG; two keep a draw-list
// record. A fixture file is compared with what generates it: `SEQ_RECORDS=update` writes records and fixtures again.

const update = process.env.SEQ_RECORDS === 'update';

/** The context a view needs: the default runtime and the view's own, answering at once (no timer in tests). */
function context(view?: ExampleView) {
  return createSeqContext({
    registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view?.registry?.({ delayMs: 0 }) ?? [])],
    ports: { assets: () => nodeFixtureAssets() },
  });
}

function compare(file: string, text: string) {
  if (update) writeFileSync(file, text);
  // a missing file is a failure: it is written on purpose, never by accident
  expect(existsSync(file), `no ${file}: run with SEQ_RECORDS=update`).toBe(true);
  expect(text).toBe(readFileSync(file, 'utf8'));
}

const record = (name: string, text: string) =>
  compare(fileURLToPath(new URL(`./records/${name}.txt`, import.meta.url)), text);

describe('the fixtures', () => {
  it('the committed gene file is what its generator writes', () => {
    compare(`${PackageFolder}${GeneFile}`, geneSeqText());
  });

  it('the committed region file is what its generator writes', () => {
    compare(`${PackageFolder}${RegionFile}`, regionSeqText());
  });
});

describe.each(ExampleViews.map((v) => [`${v.dataset}/${v.id}`, v] as const))('%s', (_, view) => {
  it('builds a spec the language accepts, which resolves and lays out with no issue', async () => {
    const ctx = context(view);
    const spec = view.build();
    expect(ctx.check(spec)).toEqual([]);
    const model = ctx.resolve(spec);
    expect(model.issues).toEqual([]);
    const resolved = ctx.view(model);
    expect(resolved?.issues).toEqual([]);
    expect(resolved?.sections.length).toBeGreaterThan(0);
    ctx.retire();
  });

  it('renders to SVG settled, with no issue', async () => {
    const ctx = context(view);
    const r = await ctx.renderSvg(view.build(), { size: { width: 720, height: 320 } });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    expect(r.svg.startsWith('<svg ')).toBe(true);
    expect(r.svg.endsWith('</svg>')).toBe(true);
    // the same spec renders the same: no clock, no random
    const again = await ctx.renderSvg(view.build(), { size: { width: 720, height: 320 } });
    expect(again.svg).toBe(r.svg);
    ctx.retire();
  });
});

describe('draw-list records', () => {
  it('p53 tracks: letters, domains in lanes, conservation with missing values', async () => {
    const ctx = context();
    const spec = ExampleViews.find((v) => v.id === 'tracks')!.build();
    const r = await ctx.renderSvg(spec, {
      size: { width: 640, height: 116 },
      viewport: { start: 84, end: 112 },
      emphasis: {
        hover: { layer: '/views/main/protein/letters:0', section: 'protein', position: 90, system: '/data/p53' },
      },
    });
    expect(r.issues).toEqual([]);
    record('p53-tracks', recordDrawLists(r.lists));
    ctx.retire();
  });

  it('p53 protein on its gene: the gene read from its file, the domains cut into exons', async () => {
    const ctx = context();
    const spec = ExampleViews.find((v) => v.id === 'protein-on-gene')!.build();
    const r = await ctx.renderSvg(spec, { size: { width: 640, height: 144 }, viewport: { start: 150, end: 160 } });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    record('p53-gene', recordDrawLists(r.lists));
    ctx.retire();
  });
});

describe('p53 gene by region', () => {
  const view = ExampleViews.find((v) => v.id === 'gene-by-region')!;

  it('shows a track for each transcript the regions brought, and TP53’s domains through their codon maps', async () => {
    const ctx = context(view);
    const r = await ctx.renderSvg(view.build(), { size: { width: 720, height: 320 } });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    for (const id of Object.keys(Transcripts)) expect(r.svg).toContain(`Transcript ${id}<`);
    for (const id of ['T1', 'T2', 'T3']) expect(r.svg).toContain(`Domains ${id}-domains<`);
    expect(r.svg).not.toContain('Domains W1');
    ctx.retire();
  });

  it('zoomed in on a stretch, only the transcripts its regions bring', async () => {
    const ctx = context(view);
    const r = await ctx.renderSvg(view.build(), {
      size: { width: 720, height: 320 },
      viewport: { start: 4200, end: 4800 },
    });
    expect(r.issues).toEqual([]);
    expect(r.svg).toContain('Transcript W1<');
    expect(r.svg).not.toContain('Transcript T1<');
    ctx.retire();
  });
});

describe('the Node assets adapter', () => {
  const read = (uri: string) => nodeFixtureAssets().reader.readText(uri, new AbortController().signal);

  it('reads a fixture by its location relative to the package folder', async () => {
    const text = await read(GeneFile);
    expect(text).toContain('"seq"');
    expect(text).toBe(geneSeqText());
  });

  it('refuses the network, an absolute path and a path above the folder', async () => {
    await expect(read('https://example.org/a.seq')).rejects.toThrow(/not read/);
    await expect(read('http://localhost/a.seq')).rejects.toThrow(/not read/);
    await expect(read('file:///etc/passwd')).rejects.toThrow(/not read/);
    await expect(read('//example.org/a.seq')).rejects.toThrow(/not read/);
    await expect(read('../package.json')).rejects.toThrow(/not read/);
    await expect(read('fixtures/../../package.json')).rejects.toThrow(/not read/);
  });

  it('reads the fixtures folder alone: not the package’s other files, its modules or an absolute path', async () => {
    await expect(read('package.json')).rejects.toThrow(/not read/);
    await expect(read('src/index.ts')).rejects.toThrow(/not read/);
    await expect(read('node_modules/typescript/package.json')).rejects.toThrow(/not read/);
    await expect(read('/etc/passwd')).rejects.toThrow(/not read/);
    await expect(read(`${PackageFolder}package.json`)).rejects.toThrow(/not read/);
    await expect(read('fixtures/p53/none.seq')).rejects.toThrow(/not read/);
    await expect(read('fixtures/p53')).rejects.toThrow();
  });

  it('refuses a link that leads out of the folder, by the real path', async () => {
    const dir = mkdtempSync(join(tmpdir(), 'seqx-assets-'));
    try {
      mkdirSync(join(dir, 'fixtures'));
      writeFileSync(join(dir, 'secret.txt'), 'secret');
      writeFileSync(join(dir, 'fixtures', 'ok.txt'), 'fine');
      symlinkSync(join(dir, 'secret.txt'), join(dir, 'fixtures', 'link.txt'));
      symlinkSync(dir, join(dir, 'fixtures', 'up'));
      const reader = nodeFixtureAssets(dir).reader;
      const signal = new AbortController().signal;
      expect(await reader.readText('fixtures/ok.txt', signal)).toBe('fine');
      await expect(reader.readText('fixtures/link.txt', signal)).rejects.toThrow(/not read/);
      await expect(reader.readText('fixtures/up/secret.txt', signal)).rejects.toThrow(/not read/);
    } finally {
      rmSync(dir, { recursive: true, force: true });
    }
  });

  it('a render without the port reads nothing: the source is an issue, never a read', async () => {
    const spec = ExampleViews.find((v) => v.id === 'protein-on-gene')!.build();
    const none = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    const r = await none.renderSvg(spec, { size: { width: 640, height: 144 } });
    expect(r.issues.map((i) => [i.code, i.address])).toContainEqual(['seq.source.unreadable', '/data/tp53-gene']);
    none.retire();
  });
});
