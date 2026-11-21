/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { GenomeFormatsViews } from '@molstar/seqx-examples';
import { CytobandsSection, MisreadSection } from '@molstar/seqx-examples/genomes/bed';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { statusIssues } from '../../../examples/status-issues.js';
import { fileBytes, fileHead, filesOf, printSources } from '../../files/file-head.js';

// Features/Genome formats/BED, headless on a fake clock, composed as the story composes it (`example-root.ts`): the
// default viewer spec with the view's registry (the `bed` runtime, its files from their recordings, each answer after
// the twin's delay); the canvas and the SVG of each view, and the file head and source nodes the story prints. The
// fixtures are read up front (the fake clock drives no disk read).

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
/** The fixtures' files in memory (a recording and the plain files), read by location as the page reads them. */
const files = new Map<string, Uint8Array>();
for (const f of [
  'fixtures/genomes/cytoBandIdeo.txt.gz.pages.json',
  'fixtures/genomes/generated-genes.bed',
  'fixtures/genomes/misread.bed',
])
  files.set(f, new Uint8Array(readFileSync(`${PackageFolder}${f}`)));
const bytesOf = async (uri: string) => {
  const b = files.get(uri);
  if (!b) throw new Error(`no fixture ${uri}`);
  return b;
};
const reader = { readText: async (uri: string) => decoder.decode(await bytesOf(uri)), readBytes: bytesOf };

const view = (id: string) => GenomeFormatsViews.find((v) => v.id === id)!;

async function open(id: string) {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), ...(view(id).registry?.({ delayMs: DelayMs }) ?? [])],
    ports: { assets: () => ({ reader }) },
  });
  viewers.push(v);
  v.resize(720);
  void v.load(view(id).build()).then(() => v.settle());
  // the fake clock's time, and between its steps a real turn: a gzip file inflates on the platform's thread pool
  for (let k = 0; k < 50; k++) {
    await vi.advanceTimersByTimeAsync(100);
    await new Promise((resolve) => setImmediate(resolve));
  }
  v.draw();
  return v;
}

/** What these tests read of an items layer's layout: its placeholders and its rows. */
interface ItemsBase {
  readonly pending: readonly unknown[];
  readonly failed: readonly unknown[];
  readonly rows: { readonly n: number };
}

/** A track's first layer's layout (an items layer's). */
const layOf = (v: SeqViewer, label: string) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === label)?.layers[0].lay as ItemsBase | undefined;

describe('Features/Genome formats/BED, headless', () => {
  it('the cytobands: chr17’s bands at the fit, the ideogram above, the SVG settled with no issue', async () => {
    const v = await open('bed-cytobands');
    expect(v.layoutIssues.value).toEqual([]);
    const r = v.view.get(CytobandsSection)!;
    expect([r.start, r.end]).toEqual([0, 83_257_441]);
    const bands = layOf(v, 'Bands')!;
    expect(bands.pending).toEqual([]);
    expect(bands.failed).toEqual([]);
    // chr17's 24 bands in hg38's cytoBandIdeo
    expect(bands.rows.n).toBe(24);
    const svg = await v.renderSvg({ idPrefix: 'bed-', navigators: 'include' });
    expect(svg.issues).toEqual([]);
    expect(svg.settled).toBe(true);
    expect(svg.svg).toContain('p13.3');
  });

  it('prints the cytobands’ head (the gzip file inflated from its recording) and its source node', async () => {
    const spec = view('bed-cytobands').build();
    const [file] = filesOf(spec).filter((f) => f.format === 'bed');
    expect(file.uri).toBe('fixtures/genomes/cytoBandIdeo.txt.gz');
    vi.useRealTimers();
    const head = fileHead(await fileBytes(reader, file.uri, new AbortController().signal), file.format).split('\n');
    expect(head[0]).toBe('chr1\t0\t2300000\tp36.33\tgneg');
    expect(head.length).toBe(15);
    const printed = printSources(spec);
    for (const word of ['chrom-sizes', 'stain', 'cytoBandIdeo']) expect(printed).toContain(word);
  });

  it('the BED12 genes: a gene’s blocks one item each, drawn and named; no issue', async () => {
    const v = await open('bed-genes');
    expect(v.layoutIssues.value).toEqual([]);
    const genes = layOf(v, 'Genes')!;
    expect(genes.pending).toEqual([]);
    expect(genes.rows.n).toBeGreaterThan(24);
    const svg = await v.renderSvg({ idPrefix: 'bed-' });
    expect(svg.issues).toEqual([]);
    expect(svg.settled).toBe(true);
    expect(svg.svg).toContain('gene1');
  });

  it('the misread file: its good rows drawn, every fault said once with its line numbers', async () => {
    const v = await open('bed-misread');
    expect(v.view.get(MisreadSection)).toBeDefined();
    const rows = layOf(v, 'Rows read')!;
    expect(rows.pending).toEqual([]);
    // 8 lines, a comment and 2 left out: 5 rows
    expect(rows.rows.n).toBe(5);
    const svg = await v.renderSvg({ idPrefix: 'bed-' });
    expect(svg.settled).toBe(true);
    // what the source said when it opened: the story's status line (`example-root.ts`) shows the store's reports
    const issues = statusIssues(svg.issues, [...v.modelIssues.value, ...v.store.issues.value]);
    const said = issues.map((i) => `${i.code}: ${i.message}`);
    expect(said).toEqual(
      expect.arrayContaining([
        expect.stringMatching(/^seq\.source\.bed\.line: An end before its start: 1 \(line 3\)$/),
        expect.stringMatching(/^seq\.source\.bed\.line: .*: 1 \(line 4\)$/),
        expect.stringMatching(/^seq\.source\.bed\.cells: A strand not \+, - or \.: read as \.: 1 \(line 5\)$/),
        expect.stringMatching(/^seq\.source\.bed\.cells: A number field’s cell is no number.*\(line 6\)$/),
        expect.stringMatching(/^seq\.source\.bed\.whitespace: .*\(line 7\)$/),
        expect.stringMatching(/^seq\.source\.bed\.cells: More columns than the first line.*\(line 8\)$/),
      ]),
    );
    expect(new Set(issues.map((i) => i.code))).toEqual(
      new Set(['seq.source.bed.line', 'seq.source.bed.cells', 'seq.source.bed.whitespace']),
    );
  });
});
