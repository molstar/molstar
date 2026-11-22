/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { type ExampleView, FilesViews, GenomeFormatsViews } from '@molstar/seqx-examples';
import { MisreadGenbankSection, PlasmidRingSection } from '@molstar/seqx-examples/genomes/genbank';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { statusIssues } from '../../../examples/status-issues.js';
import { withViewRegistry } from '../../../examples/view-registry.js';
import { fileText, filesOf, printSources } from '../../files/file-head.js';

// Features/Genome formats/GenBank, headless on a fake clock, composed as the story composes it (`example-root.ts`): the
// default viewer spec with the view's registry (the `genbank` runtime, each answer after the twin's delay; the GFF3
// view's records runtime), the GenBank view and the GFF3 view under it in one page; their canvas draw lists, their SVGs,
// the misread file's status line and the records and head the story prints. The fixtures are read up front (the fake
// clock drives no disk read).

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
  'fixtures/plasmid/pUC19c-rotated.gb',
  'fixtures/plasmid/pUC19c-rotated.seq',
  'fixtures/plasmid/pUC19c-rotated-features.gff3',
  'fixtures/genomes/misread.gb',
])
  files.set(f, new Uint8Array(readFileSync(`${PackageFolder}${f}`)));
const bytesOf = async (uri: string) => {
  const b = files.get(uri);
  if (!b) throw new Error(`no fixture ${uri}`);
  return b;
};
const reader = { readText: async (uri: string) => decoder.decode(await bytesOf(uri)), readBytes: bytesOf };

const view = (id: string) => GenomeFormatsViews.find((v) => v.id === id)!;
const gff3 = FilesViews.find((v) => v.id === 'plasmid-gff3')!;

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

async function open(v: ExampleView) {
  const base = DefaultSeqViewerSpec();
  const opens: Promise<unknown>[] = [];
  const viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: withViewRegistry(base.registry ?? [], watched(v.registry?.({ delayMs: DelayMs }) ?? [], opens)),
    ports: { assets: () => ({ reader }) },
  });
  viewers.push(viewer);
  viewer.resize(720);
  void viewer.load(v.build()).then(() => viewer.settle());
  for (let k = 0; k < 40; k++) {
    await vi.advanceTimersByTimeAsync(200);
    await Promise.allSettled(opens);
  }
  viewer.draw();
  return viewer;
}

describe('Features/Genome formats/GenBank, headless', () => {
  it('the plasmid from GenBank and from GFF3 in one page: their canvas draw lists equal, their SVGs settled', async () => {
    const gb = await open(view('genbank-plasmid'));
    const gff = await open(gff3);
    expect(gb.layoutIssues.value).toEqual([]);
    expect(gb.view.get(PlasmidRingSection)).toBeDefined();
    const a = recordDrawLists(gb.frame!.lists);
    const b = recordDrawLists(gff.frame!.lists);
    expect(a.split('\n').length).toBeGreaterThan(100);
    expect(a).toBe(b);
    // the twin is not empty: bla's label drawn
    const svg = await gb.renderSvg({ idPrefix: 'gb-' });
    expect(svg.settled).toBe(true);
    expect(svg.issues).toEqual([]);
    expect(svg.svg).toContain('bla (AmpR)');
  });

  it('a misread file: what reads drawn, each fault said once with its lines in the status line', async () => {
    const v = await open(view('genbank-misread'));
    expect(v.view.get(MisreadGenbankSection)).toBeDefined();
    const issues = statusIssues([], [...v.modelIssues.value, ...v.store.issues.value]);
    const codes = new Set(issues.map((i) => i.code));
    for (const c of ['location', 'approximate', 'length', 'record', 'header', 'translation'])
      expect(codes, c).toContain(`seq.source.genbank.${c}`);
    const location = issues.find((i) => i.code === 'seq.source.genbank.location')!;
    expect(issues.filter((i) => i.code === 'seq.source.genbank.location').length).toBe(4);
    expect(location.message).toMatch(/line \d+/);
    const svg = await v.renderSvg({ idPrefix: 'gb-' });
    expect(svg.settled).toBe(true);
  });

  it('prints the file’s records and its head, and the source node', async () => {
    const spec = view('genbank-plasmid').build();
    const [file] = filesOf(spec);
    expect(file).toEqual({ uri: 'fixtures/plasmid/pUC19c-rotated.gb', format: 'genbank' });
    vi.useRealTimers();
    const text = await fileText(reader, file, new AbortController().signal);
    expect(text.split('\n')[0]).toBe('pUC19c: 2686 bp, circular, DNA');
    expect(text).toContain('Cloning vector pUC19c; 6 features');
    expect(text).toContain('LOCUS       pUC19c');
    expect(printSources(spec)).toContain('genbank');
  });
});
