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
import type { HoverCardBehavior } from '@molstar/seqx-viewer/behaviors/hover-card';
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

/** Where item `i`'s label would start: its model's start on the screen (its hairline's, else its first box's) less
 * the gap and the label's width. */
function labelRoom(lay: RangesLay, i: number): number {
  for (let c = 0; c < lay.cn; c++) if (lay.citem[c] === i) return lay.canchor[c] - 4 - lay.gene.labelW[i];
  let x = Number.POSITIVE_INFINITY;
  for (let j = 0; j < lay.m; j++) if (lay.keys[lay.mitem[j]] === lay.keys[i]) x = Math.min(x, lay.mx0[j]);
  return x - 4 - lay.gene.labelW[i];
}

/** Scrolls the Transcripts band (its lanes scroll inside the capped track) until the item `key` is drawn. */
function reveal(v: SeqViewer, key: string): RangesLay {
  const track = '/views/main/chr17/transcripts';
  const first = layOf(v, 'Transcripts')!;
  for (let top = 0; top <= first.scrollMax + first.room; top += Math.max(1, first.room - 2)) {
    v.view.setBand(track, { bandTop: top });
    v.draw();
    const lay = layOf(v, 'Transcripts')!;
    if (drawn(lay).has(key)) return lay;
  }
  throw new Error(`${key} is in no lane`);
}

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
      // the track holds its cap's lanes (16) and scrolls to the rest
      expect(lanes.size).toBeGreaterThan(15);
      expect(tx.scrollMax).toBeGreaterThan(0);
      expect(drawn(reveal(v, canonical)).has(canonical)).toBe(true);
      for (const l of lanes.values()) expect(l.size).toBe(1);
      expect(Array.from(tx.mthin.subarray(0, tx.m)).some((t) => t === 1)).toBe(true);
      expect(Array.from(tx.mthin.subarray(0, tx.m)).some((t) => t === 0)).toBe(true);
      // TP53 on the minus strand: its hairline's chevrons point left
      const k = Array.from({ length: tx.cn }, (_, j) => j).find((j) => String(tx.keys[tx.citem[j]]) === canonical)!;
      expect(tx.cdir[k]).toBe(2);
      // one model a gene, each in a lane of its own and no density row (no stray feature, no overflow)
      const gl = layOf(v, 'Genes')!;
      expect([...drawn(gl).keys()].sort()).toEqual(
        id === 'gff-tp53'
          ? ['gene:ENSG00000141499', 'gene:ENSG00000141510']
          : ['ENSG00000141499.20', 'ENSG00000141510.21'],
      );
      expect(gl.overflow).toBe(0);
      expect(tx.overflow).toBe(0);
      // each transcript whose label's room lies inside the band has its label, whole, as text left of its model
      const labels = recordDrawLists(v.frame!.lists, ['labels'])
        .split('\n')
        .filter((l) => l.startsWith('text '));
      expect(labels.some((l) => l.includes('…'))).toBe(false);
      const g = tx.gene;
      const room = Array.from({ length: tx.items }, (_, i) => i).filter(
        (i) => g.repOf[i] === i && tx.lanes[i] >= 0 && g.labelW[i] > 0 && labelRoom(tx, i) >= tx.box.x,
      );
      expect(room.length).toBeGreaterThan(10);
      // a label with no room before its model floats at the data's left edge: every shown lane is labelled
      const shown = Array.from({ length: tx.items }, (_, i) => i).filter(
        (i) => g.repOf[i] === i && tx.lanes[i] >= 0 && g.labelW[i] > 0,
      );
      expect(g.ln).toBe(shown.length);
      expect(g.ln - room.length).toBe(Array.from(g.lfloat.subarray(0, g.ln)).filter((f) => f === 1).length);
      const svg = await v.renderSvg({ idPrefix: `${id}-` });
      expect(svg.issues).toEqual([]);
      expect(svg.settled).toBe(true);
      expect(svg.svg).toContain('TP53');
    });

  it('a transcript is one item: hovered, every exon lights', async () => {
    const v = await open('gff-tp53');
    const tx = reveal(v, 'transcript:ENST00000269305');
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

  it('the Transcripts track is sized to the lanes its source read, to its cap; its lanes start at its top, scroll inside it', async () => {
    const v = await open('gff-tp53');
    const tx = layOf(v, 'Transcripts')!;
    const config = v.context.renderConfig;
    // 47 transcripts: past the cap (240 px, 16 lanes and padding)
    expect(tx.box.height).toBe(config.track.maxPx);
    expect(tx.laneCount).toBeGreaterThan(tx.room);
    expect(tx.scrollMax).toBe(tx.laneCount - tx.room);
    // the first lane at the band's padding: no blank band above the lanes
    expect(tx.laneY - tx.box.y).toBe(config.lanePadding);
    // the Genes track (WRAP53 and TP53) keeps its own height: two lanes, not the Transcripts' (tracks of one kind share)
    const genes = layOf(v, 'Genes')!;
    expect(genes.box.height).toBe(2 * config.laneHeight + 2 * config.lanePadding);
  });

  it('the hover card of a transcript deep in the tall track is anchored by its lane, not the track’s top', async () => {
    const v = await open('gff-tp53');
    const tx = layOf(v, 'Transcripts')!;
    let j = 0;
    for (let i = 0; i < tx.m; i++) if (tx.mlane[i] > tx.mlane[j]) j = i;
    const x = (tx.mx0[j] + tx.mx1[j]) / 2;
    const y = tx.laneY + tx.mlane[j] * tx.laneHeight + tx.laneHeight / 2;
    expect(tx.mlane[j]).toBeGreaterThan(10);
    v.input.move(x, y);
    v.draw();
    await vi.advanceTimersByTimeAsync(2000);
    const card = (v.behaviors.get('hoverCard') as HoverCardBehavior).shown.value;
    expect(card).toBeDefined();
    const a = card!.place.anchor;
    expect(a.y).toBeLessThanOrEqual(y);
    expect(a.y + a.height).toBeGreaterThanOrEqual(y);
    expect(a.height).toBeLessThanOrEqual(2 * tx.laneHeight);
  });

  it('a name found by the file’s index and focused, as FindByName does (EC-228)', async () => {
    const v = await open('gff-tp53');
    const before = v.view.get(Tp53Section)!;
    const text = findAndFocus(v, 'gff3', 'TP53-201', reader);
    await vi.advanceTimersByTimeAsync(2000);
    const said = await text;
    const lines = said.split('\n');
    expect(lines[0]).toBe(
      'Searched “TP53-201” in Name, ID, Alias (exact names first, then names starting with it): 1 hit',
    );
    expect(lines[1]).toBe('Found: TP53-201 at 17:7,668,421-7,687,490 (its Name)');
    expect(lines[2]).toMatch(/^The view went to 17:7,66\d,\d{3}-7,68\d,\d{3}: TP53-201’s extent$/);
    await vi.advanceTimersByTimeAsync(2000);
    const after = v.view.get(Tp53Section)!;
    expect([after.start, after.end]).not.toEqual([before.start, before.end]);
    // TP53-201's extent (7,668,420 to 7,687,490 in release 116, 0-based) on the section over 7,660,000 to 7,690,000
    expect(after.start).toBeLessThanOrEqual(8420);
    expect(after.start).toBeGreaterThan(7000);
    expect(after.end).toBeGreaterThanOrEqual(27_490);
    expect(after.end).toBeLessThan(29_000);
  });

  it('`.` and `,` step the cursor to the next and previous feature of the track, over the stretches with none', async () => {
    const v = await open('gff-tp53');
    const section = v.frame!.layout.sections[0];
    const track = section.tracks.find((t) => t.layers.some((l) => l.type.items?.(l.lay).length))!;
    const starts = [...new Set(track.layers.flatMap((l) => (l.type.items?.(l.lay) ?? []).map((i) => i.start)))].sort(
      (a, b) => a - b,
    );
    expect(starts.length).toBeGreaterThan(1);
    const system = section.section.axis.system;
    // before the first feature, then each feature in turn, then the last: no step
    v.selection.setCursor({ section: section.section.address, system, position: starts[0] - 1, keyboard: true });
    const at = () => v.selection.cursor.value?.position;
    expect(v.input.key({ key: '.' })).toBe(true);
    expect(at()).toBe(starts[0]);
    expect(v.input.key({ key: '.' })).toBe(true);
    expect(at()).toBe(starts[1]);
    expect(v.input.key({ key: ',' })).toBe(true);
    expect(at()).toBe(starts[0]);
    // between two features: the next one starts after the cursor, the previous before it
    const gap = Math.floor((starts[0] + starts[1]) / 2);
    expect(gap).toBeGreaterThan(starts[0]);
    v.selection.setCursor({ section: section.section.address, system, position: gap, keyboard: true });
    v.input.key({ key: '.' });
    expect(at()).toBe(starts[1]);
    v.selection.setCursor({ section: section.section.address, system, position: gap, keyboard: true });
    v.input.key({ key: ',' });
    expect(at()).toBe(starts[0]);
    v.selection.setCursor({ section: section.section.address, system, position: starts.at(-1)!, keyboard: true });
    expect(v.input.key({ key: '.' })).toBe(false);
    expect(at()).toBe(starts.at(-1));
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
