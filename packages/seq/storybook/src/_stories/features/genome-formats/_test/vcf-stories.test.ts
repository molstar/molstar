/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { GenomeFormatsViews } from '@molstar/seqx-examples';
import { ClinvarSection, GeneratedVcfSection, MisreadVcfSection } from '@molstar/seqx-examples/genomes/vcf';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import type { RangesLay } from '@molstar/seqx-core/render';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { statusIssues } from '../../../examples/status-issues.js';
import { withViewRegistry } from '../../../examples/view-registry.js';
import { clickAt } from '../../selection/compose.js';
import { fileText, filesOf, printSources } from '../../files/file-head.js';

// Features/Genome formats/VCF, headless on a fake clock, composed as the story composes it (`example-root.ts`): the
// default viewer spec with the view's registry (the `vcf` runtime and the `bed` one for chrom.sizes, the files from
// their recordings, each answer after the twin's delay); the canvas and the SVG of each view, and the header summary
// and source nodes the story prints. The fixtures are read up front (the fake clock drives no disk read).

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
/** The fixtures' files in memory (the recordings and the plain files), read by location as the page reads them. */
const files = new Map<string, Uint8Array>();
for (const f of [
  'fixtures/genomes/clinvar_20240107.vcf.gz.pages.json',
  'fixtures/genomes/clinvar_20240107.vcf.gz.tbi.pages.json',
  'fixtures/genomes/generated-variants.vcf',
  'fixtures/genomes/misread.vcf',
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
  // the frames that start the opens, then the opens themselves (BGZF blocks inflate on the platform's thread pool, by
  // no clock: awaited), then the twin's delays and the rest on the fake clock
  await vi.advanceTimersByTimeAsync(200);
  expect(opens.length).toBeGreaterThan(0);
  for (let seen = -1; seen !== opens.length; ) {
    seen = opens.length;
    await Promise.allSettled(opens);
    await vi.advanceTimersByTimeAsync(200);
  }
  await settled(v);
  return v;
}

/** Whether a track's layer still waits for data. */
const waiting = (v: SeqViewer) =>
  v.frame!.layout.sections.some((s) =>
    s.tracks.some((t) => ((t.layers[0]?.lay as RangesLay)?.pending?.length ?? 0) > 0),
  );

/**
 * The viewer's region reads done (each a BGZF inflate on the platform's thread pool, by no clock: a turn of the event
 * loop each step; then the twin's delay on the fake clock), then a frame: never a count of turns, a bound on the steps.
 */
async function settled(v: SeqViewer) {
  for (let k = 0; k < 200; k++) {
    await vi.advanceTimersByTimeAsync(100);
    await new Promise<void>((done) => setImmediate(done));
    v.draw();
    if (k > 10 && !waiting(v)) break;
  }
  v.draw();
}

const layOf = (v: SeqViewer, label: string) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === label)?.layers[0].lay as RangesLay | undefined;

/** The hit mark of the item `key` (its own, before `bars` merges a lane's runs): its centre on the canvas. */
function markOf(lay: RangesLay, key: string): [number, number] | undefined {
  for (let j = 0; j < lay.hn; j++) {
    if (String(lay.keys[lay.hitem[j]]) !== key) continue;
    return [(lay.hx0[j] + lay.hx1[j]) / 2, lay.laneY + lay.hlane[j] * lay.laneHeight + lay.laneHeight / 2];
  }
  return undefined;
}

/** Each item's hit mark's width, by its key. */
function widthsOf(lay: RangesLay): Map<string, number> {
  const out = new Map<string, number>();
  for (let j = 0; j < lay.hn; j++) out.set(String(lay.keys[lay.hitem[j]]), lay.hx1[j] - lay.hx0[j]);
  return out;
}

/** TP53's region, the section's range (its positions are the section's from 7,660,000). */
const From = 7_660_000;

describe('Features/Genome formats/VCF, headless', () => {
  it('TP53’s ClinVar variants: read by region, an allele a mark coloured by category; the SVG settled with no issue', async () => {
    const v = await open('vcf-clinvar');
    expect(v.layoutIssues.value).toEqual([]);
    const r = v.view.get(ClinvarSection)!;
    expect([r.start + From, r.end + From]).toEqual([7_660_000, 7_690_000]);
    const lay = layOf(v, 'ClinVar')!;
    expect(lay.pending).toEqual([]);
    expect(lay.failed).toEqual([]);
    expect(lay.items).toBeGreaterThan(2000);
    // the categories' colours: more than one fill among the marks
    expect(new Set(lay.fills.slice(0, lay.items)).size).toBeGreaterThan(2);
    const svg = await v.renderSvg({ idPrefix: 'vcf-' });
    expect(svg.issues).toEqual([]);
    expect(svg.settled).toBe(true);
  });

  it('TP53’s drawn marks carry their category’s fill at rest: a merged run takes its majority category’s, no rect is neutral', async () => {
    const v = await open('vcf-clinvar');
    const lay = layOf(v, 'ClinVar')!;
    const band = lay.box;
    const hex = (c: number) => `#${(c >>> 0).toString(16).padStart(8, '0')}`;
    const neutral = hex(lay.tokens.colors['seq-range']);
    // the data pass's rects inside the ClinVar band, by fill (the rect line's fifth number)
    const rects = recordDrawLists(v.frame!.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => l.split(' '))
      .filter((p) => Number(p[2]) >= band.y && Number(p[2]) < band.y + band.height);
    const fills = new Set(rects.map((p) => p[5]));
    // SNVs, deletions, insertions, indels, MNVs: their categories' fills
    expect(fills.size).toBeGreaterThanOrEqual(4);
    expect(rects.filter((p) => p[5] === neutral)).toEqual([]);
    // no merged run is the neutral mark (`-1`): a mixed one is its majority's
    for (let j = 0; j < lay.m; j++) expect(lay.mitem[j], `run ${j} at ${lay.mx0[j]}`).not.toBe(-1);
  });

  it('an allele: hovered it lights, clicked its REF extent is selected, double-clicked it is focused', async () => {
    const v = await open('vcf-clinvar');
    // TP53 R175H (ClinVar 12374): zoomed in about it first, so its mark is its own
    v.runCommand('address.focus', { address: '/data/chr17/(range :from 7675060 :to 7675120)' });
    await settled(v);
    const lay = layOf(v, 'ClinVar')!;
    const at = markOf(lay, '7675088:C>T');
    expect(at, 'R175H drawn').toBeDefined();
    const [x, y] = at!;
    expect(v.hit(x, y)?.element?.key).toBe('7675088:C>T');
    v.input.move(x, y);
    v.draw();
    expect(recordDrawLists(v.frame!.lists, ['hover'])).toContain('rect ');
    clickAt(v, [x, y]);
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['7675088:C>T']);
    v.input.doubleClick(x, y);
    await settled(v);
    const focused = v.view.get(ClinvarSection)!;
    // the allele's extent (7,675,087 to 7,675,088) in the view, the view narrower than before
    expect(focused.start + From).toBeLessThanOrEqual(7_675_087);
    expect(focused.end + From).toBeGreaterThanOrEqual(7_675_088);
    expect(focused.end - focused.start).toBeLessThan(60);
  });

  it('prints ClinVar’s header summary (read by its index from the recording) and its source node', async () => {
    const spec = view('vcf-clinvar').build();
    const [file] = filesOf(spec).filter((f) => f.format === 'vcf');
    expect(file).toEqual({
      uri: 'fixtures/genomes/clinvar_20240107.vcf.gz',
      format: 'vcf',
      index: 'fixtures/genomes/clinvar_20240107.vcf.gz.tbi',
    });
    vi.useRealTimers();
    const text = await fileText(reader, file, new AbortController().signal);
    const lines = text.split('\n');
    expect(lines[0]).toMatch(/^VCFv4\.1: \d+ header lines/);
    expect(text).toContain('Contigs: none');
    expect(text).toContain('CLNSIG (Number=., String)');
    expect(text).toContain('Samples: 0');
    const printed = printSources(spec);
    for (const word of ['clinvar_20240107.vcf.gz.tbi', 'chrom-sizes', 'LicenseRef-NCBI-public-data'])
      expect(printed).toContain(word);
  });

  it('ClinVar’s chr17 at the fit: the too-many-tiles state (the band hatched as waiting, never empty), nothing read', async () => {
    const v = await open('vcf-zoom-in');
    const r = v.view.get(ClinvarSection)!;
    expect([r.start, r.end]).toEqual([0, 83_257_441]);
    const lay = layOf(v, 'ClinVar')!;
    expect(lay.items).toBe(0);
    expect(lay.failed).toEqual([]);
    // the whole band waiting: 81,307 level-1 tiles, more than a want asks (`maxTilesPerWant`), so none is read
    expect(lay.pending.length).toBeGreaterThan(0);
    expect(Math.min(...lay.pending.filter((_, i) => i % 3 === 0))).toBeLessThanOrEqual(lay.box.x + 1);
    expect(Math.max(...lay.pending.filter((_, i) => i % 3 === 1))).toBeGreaterThanOrEqual(
      lay.box.x + lay.box.width - 1,
    );
  });

  it('a long deletion over SNVs: one item each, a bar across the SNVs it covers (EC-222)', async () => {
    const v = await open('vcf-long-deletion');
    expect(v.layoutIssues.value).toEqual([]);
    expect(v.view.get(GeneratedVcfSection)).toBeDefined();
    const lay = layOf(v, 'Variants')!;
    expect(lay.pending).toEqual([]);
    const widths = widthsOf(lay);
    const sv1 = widths.get('4001:N><DEL>@9000')!;
    // 4,000 to 9,000 of 20,000 at 900 px: a bar about a quarter of the axis wide, the SNVs it covers marks beside it
    expect(sv1).toBeGreaterThan(150);
    const ids = [...widths.keys()];
    expect(ids.filter((k) => /^\d+:[ACGT]>[ACGT]$/.test(k)).length).toBeGreaterThan(20);
    expect(ids.some((k) => /^12001:600:[0-9a-f]{8}>[ACGT]$/.test(k))).toBe(true);
    const svg = await v.renderSvg({ idPrefix: 'vcf-' });
    expect(svg.issues).toEqual([]);
    expect(svg.settled).toBe(true);
  });

  it('the misread file: its good rows drawn, every fault said once with its line numbers', async () => {
    const v = await open('vcf-misread');
    expect(v.view.get(MisreadVcfSection)).toBeDefined();
    const lay = layOf(v, 'Read')!;
    expect(lay.pending).toEqual([]);
    // good1, bad-depth, star's C, past-end (clipped) and good2
    expect(lay.items).toBe(5);
    const svg = await v.renderSvg({ idPrefix: 'vcf-' });
    expect(svg.settled).toBe(true);
    const issues = statusIssues(svg.issues, [...v.modelIssues.value, ...v.store.issues.value]);
    const said = issues.map((i) => `${i.code}: ${i.message}`);
    expect(said).toEqual(
      expect.arrayContaining([
        expect.stringMatching(/^seq\.source\.vcf\.header: .*: 1$/),
        expect.stringMatching(/^seq\.source\.vcf\.line: A position that is not a whole number.*\(line 7\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.line: No REF: 1 \(line 8\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.line: Fewer than eight columns: 1 \(line 9\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.cell: A QUAL that is no number.*\(line 10\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.cell: An integer field’s cell is no int32.*\(line 10\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.info: An INFO key the header does not declare: 1 \(XX\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.allele: An ALT `\*`.*\(line 11\)$/),
        expect.stringMatching(/^seq\.source\.vcf\.range: .*: 2 \(line 12, line 13\)$/),
      ]),
    );
  });
});
