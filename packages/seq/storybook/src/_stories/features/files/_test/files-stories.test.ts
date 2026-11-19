/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { type ExampleView, FilesViews, TableFiles } from '@molstar/seqx-examples';
import { nodeFixtureAssets, PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { fileHead, filesOf, printSources } from '../file-head.js';

// The Files stories, headless: the view of each composed as the story composes it (the default viewer spec; its files read
// through the Node adapter of the assets port, the way the page reads them through the Vite plugin), its canvas and its
// SVG, and the file's head and source node the story prints beside it.

const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
});

async function open(view: ExampleView): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), ...(view.registry?.() ?? [])],
    ports: { assets: () => nodeFixtureAssets() },
  });
  viewers.push(v);
  v.resize(720);
  await v.load(view.build());
  await v.settle();
  v.draw();
  return v;
}

const view = (id: string) => FilesViews.find((v) => v.id === id)!;
const bytesOf = (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`));

describe('Features/Files, headless', () => {
  it.each(FilesViews.map((v) => [v.id, v] as const))('%s: canvas and SVG settled, with no issue', async (_, v) => {
    const viewer = await open(v);
    expect(viewer.layoutIssues.value).toEqual([]);
    const svg = await viewer.renderSvg({ idPrefix: 'f-' });
    expect(svg.issues).toEqual([]);
    expect(svg.settled).toBe(true);
    expect(svg.svg.startsWith('<svg ')).toBe(true);
  });

  it('every form of the domains table draws the same canvas', async () => {
    const records: string[] = [];
    for (const [, id] of TableFiles) {
      const viewer = await open(view(id));
      records.push(recordDrawLists(viewer.frame!.lists, ['data']));
    }
    for (const r of records) expect(r).toBe(records[0]);
    expect(records[0]).toContain('rect');
  });

  it('prints the head of the file and the source node, as the story does', () => {
    // UniProt's entry: the JSON as the service writes it, and the transform written out in the source node
    const entry = view('uniprot').build();
    const [file] = filesOf(entry);
    expect(file).toEqual({ uri: 'fixtures/p53/P04637.json', format: 'json' });
    const head = fileHead(bytesOf(file.uri), file.format);
    expect(head).toMatch(/"primaryAccession": "P04637"/);
    expect(head.split('\n').length).toBeLessThanOrEqual(15);
    const source = printSources(entry);
    for (const word of ['location.start.value', 'sequence.value', 'features']) expect(source).toContain(word);

    const tsv = filesOf(view('table-tsv').build())[0];
    expect(fileHead(bytesOf(tsv.uri), tsv.format).split('\n')[0]).toBe('id\tlabel\tstart\tend\tcategory');
    const csv = filesOf(view('table-csv').build())[0];
    expect(fileHead(bytesOf(csv.uri), csv.format).split('\n')[0]).toBe('id,label,start,end,category');
    // a MessagePack file prints as the JSON it holds
    const packed = filesOf(view('table-msgpack-aos').build())[0];
    expect(fileHead(bytesOf(packed.uri), packed.format)).toMatch(/^MessagePack, \d+ bytes, as JSON:\n\[/);
    const gff = filesOf(view('plasmid-gff3').build())[0];
    expect(fileHead(bytesOf(gff.uri), gff.format).split('\n')[0]).toBe('##gff-version 3');
    expect(printSources(view('plasmid-gff3').build())).toContain('gff3');
  });
});
