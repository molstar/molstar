/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { GenomeFormatsViews } from '@molstar/seqx-examples';
import { MisreadGffSection, Tp53Section } from '@molstar/seqx-examples/genomes/gff';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import type { RangesLay } from '@molstar/seqx-core/render';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { withViewRegistry } from '../../../examples/view-registry.js';
import { statusIssues } from '../../../examples/status-issues.js';
import { findAndFocus } from '../names.js';

// Features/Genome formats/GFF3 and GTF, headless on a fake clock, composed as the stories compose them
// (`example-root.ts`): the default viewer spec with the view's registry (the `gff` runtime, its files from their
// recordings, each answer after the twin's delay); the canvas and the SVG of each view; a name found and focused as the
// FindByName stories do; the misread file's issues as its status line says them. The fixtures are read up front.

const DelayMs = 150;
const raf = DomHost.requestAnimationFrame;
const viewers: SeqViewer[] = [];

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  for (const v of viewers.splice(0)) v.retire();
  vi.useRealTimers();
});

const decoder = new TextDecoder();
const files = new Map<string, Uint8Array>();
for (const f of [
  'fixtures/genomes/Homo_sapiens.GRCh38.116.chr.TP53.gff3.gz.pages.json',
  'fixtures/genomes/gencode.v50.annotation.TP53.gtf.gz.pages.json',
  'fixtures/genomes/misread.gff3',
])
  files.set(f, new Uint8Array(readFileSync(`${PackageFolder}${f}`)));
const bytesOf = async (uri: string) => {
  const b = files.get(uri);
  if (!b) throw new Error(`no fixture ${uri}`);
  return b;
};
const reader = { readText: async (uri: string) => decoder.decode(await bytesOf(uri)), readBytes: bytesOf };

const view = (id: string) => GenomeFormatsViews.find((v) => v.id === id)!;

/** The view's registry with each runtime's opens kept: a test waits on them, not on a count of turns. */
function watched(entries: readonly SeqRegistryEntry[], opens: Promise<unknown>[]): SeqRegistryEntry[] {
  return entries.map((e) =>
    e.kind === 'format-runtime'
      ? {
          ...e,
          runtime: {
            ...e.runtime,
            open: (source, ctx) => {
              const p = e.runtime.open(source, ctx);
              opens.push(p);
              return p;
            },
          } satisfies FormatRuntime,
        }
      : e,
  );
}

async function open(id: string) {
  const base = DefaultSeqViewerSpec();
  const opens: Promise<unknown>[] = [];
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: withViewRegistry(base.registry ?? [], watched(view(id).registry?.({ delayMs: DelayMs }) ?? [], opens)),
    ports: { assets: () => ({ reader }) },
  });
  viewers.push(v);
  v.resize(900);
  void v.load(view(id).build()).then(() => v.settle());
  // the frames that start the opens, then the opens themselves (a gzip file inflates on the platform's thread pool, by
  // no clock: awaited), then the twin's delays and the rest on the fake clock
  await vi.advanceTimersByTimeAsync(200);
  expect(opens.length).toBeGreaterThan(0);
  for (let seen = -1; seen !== opens.length; ) {
    seen = opens.length;
    await Promise.allSettled(opens);
    await vi.advanceTimersByTimeAsync(200);
  }
  await vi.advanceTimersByTimeAsync(5000);
  v.draw();
  return v;
}

const layOf = (v: SeqViewer, label: string) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === label)?.layers[0].lay as RangesLay | undefined;

/** The item keys a gene-line layer drew, with their lanes (one each). */
function drawn(lay: RangesLay) {
  const lanes = new Map<string, Set<number>>();
  for (let j = 0; j < lay.m; j++) {
    const k = String(lay.keys[lay.mitem[j]]);
    if (!lanes.has(k)) lanes.set(k, new Set());
    lanes.get(k)!.add(lay.mlane[j]);
  }
  return lanes;
}

describe('Features/Genome formats/GFF3 and GTF, headless', () => {
  for (const [id, canonical] of [
    ['gff-tp53', 'transcript:ENST00000269305'],
    ['gtf-tp53', 'ENST00000269305.9'],
  ] as const)
    it(`${id}: TP53’s transcripts as gene models (one lane each, thick and thin boxes, chevrons), its genes; the SVG settled with no issue`, async () => {
      const v = await open(id);
      expect(v.layoutIssues.value).toEqual([]);
      expect(v.view.get(Tp53Section)).toBeDefined();
      const tx = layOf(v, 'Transcripts')!;
      expect(tx.pending).toEqual([]);
      expect(tx.failed).toEqual([]);
      const lanes = drawn(tx);
      expect(lanes.size).toBeGreaterThan(15);
      expect(lanes.has(canonical)).toBe(true);
      for (const l of lanes.values()) expect(l.size).toBe(1);
      expect(Array.from(tx.mthin.subarray(0, tx.m)).some((t) => t === 1)).toBe(true);
      expect(Array.from(tx.mthin.subarray(0, tx.m)).some((t) => t === 0)).toBe(true);
      // TP53 on the minus strand: its hairline's chevrons point left
      const k = Array.from({ length: tx.cn }, (_, j) => j).find((j) => String(tx.keys[tx.citem[j]]) === canonical)!;
      expect(tx.cdir[k]).toBe(2);
      const genes = drawn(layOf(v, 'Genes')!);
      expect(genes.size).toBe(2);
      const svg = await v.renderSvg({ idPrefix: `${id}-` });
      expect(svg.issues).toEqual([]);
      expect(svg.settled).toBe(true);
      expect(svg.svg).toContain('TP53');
    });

  it('a transcript is one item: hovered, every exon lights', async () => {
    const v = await open('gff-tp53');
    const tx = layOf(v, 'Transcripts')!;
    const j = Array.from({ length: tx.m }, (_, i) => i).find(
      (i) => String(tx.keys[tx.mitem[i]]) === 'transcript:ENST00000269305',
    )!;
    const x = (tx.mx0[j] + tx.mx1[j]) / 2;
    const y = tx.laneY + tx.mlane[j] * tx.laneHeight + tx.laneHeight / 2;
    expect(v.hit(x, y)?.element?.key).toBe('transcript:ENST00000269305');
    v.input.move(x, y);
    v.draw();
    const marks = Array.from({ length: tx.m }, (_, i) => i).filter(
      (i) => String(tx.keys[tx.mitem[i]]) === 'transcript:ENST00000269305',
    ).length;
    const hovered = recordDrawLists(v.frame!.lists, ['hover'])
      .split('\n')
      .filter((l) => l.startsWith('rect ')).length;
    expect(hovered).toBeGreaterThanOrEqual(marks);
  });

  it('a name found by the file’s index and focused, as FindByName does (EC-228)', async () => {
    const v = await open('gff-tp53');
    const before = v.view.get(Tp53Section)!;
    const text = findAndFocus(v, 'gff3', 'TP53-201', reader);
    await vi.advanceTimersByTimeAsync(2000);
    const said = await text;
    expect(said.split('\n')[1]).toMatch(/^TP53-201 {2}\(Name\) {2}17:/);
    await vi.advanceTimersByTimeAsync(2000);
    const after = v.view.get(Tp53Section)!;
    expect([after.start, after.end]).not.toEqual([before.start, before.end]);
    // TP53-201's extent (7,668,420 to 7,687,490 in release 116, 0-based) on the section over 7,660,000 to 7,690,000
    expect(after.start).toBeLessThanOrEqual(8420);
    expect(after.start).toBeGreaterThan(7000);
    expect(after.end).toBeGreaterThanOrEqual(27_490);
    expect(after.end).toBeLessThan(29_000);
    expect(said).toMatch(/The first focused: /);
  });

  it('the misread GFF3: its good transcript drawn, every fault said once with its line numbers', async () => {
    const v = await open('gff-misread');
    expect(v.view.get(MisreadGffSection)).toBeDefined();
    const rows = drawn(layOf(v, 'Read')!);
    expect(rows.has('tx:A1')).toBe(true);
    expect(rows.has('gene:B')).toBe(true);
    const svg = await v.renderSvg({ idPrefix: 'gff-' });
    expect(svg.settled).toBe(true);
    const issues = statusIssues(svg.issues, [...v.modelIssues.value, ...v.store.issues.value]);
    const said = issues.map((i) => `${i.code}: ${i.message}`);
    expect(said).toEqual(
      expect.arrayContaining([
        expect.stringMatching(/^seq\.source\.gff\.line: An end before its start: 1 \(line 9\)$/),
        expect.stringMatching(/^seq\.source\.gff\.line: A position is not a whole number.*\(line 10\)$/),
        expect.stringMatching(/^seq\.source\.gff\.line: Not nine tab-separated columns: 1 \(line 11\)$/),
        expect.stringMatching(/^seq\.source\.gff\.parent: A `Parent` no line has: 1 \(line 12\)$/),
        expect.stringMatching(/^seq\.source\.gff\.parent: A cycle of `Parent`: 2 \(line 13, line 14\)$/),
        expect.stringMatching(/^seq\.source\.gff\.type: `gene` read as a transcript: 1$/),
      ]),
    );
  });
});
