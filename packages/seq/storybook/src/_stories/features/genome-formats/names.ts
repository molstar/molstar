/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetReader, NameHit } from '@molstar/seqx-core/data';
import { findTp53, hitAddress, hitRange, Tp53Section } from '@molstar/seqx-examples/genomes/gff';
import { Tp53Region } from '@molstar/seqx-examples/genomes/files';
import type { ExampleView } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { withParams } from '../../with-params.js';

// A name found in a GFF3 or GTF (EC-228), as the stories show it: the file's `findNames` (the GFF runtime's index, exact
// names first, then by prefix), the first hit's extent focused (`address.focus` of its range, which the hit carries: no
// tile read first), and what was searched, what was found and where the view went said. Until the locus box falls back
// to it, the story's own find box asks the index.

/** The page's own files by `fetch` (the fixtures served beside the stories). */
const fetched = (uri: string) =>
  fetch(uri).then((r) => (r.ok ? r.arrayBuffer() : Promise.reject(new Error(`${r.status}`))));
export const PageReader: AssetReader = {
  readText: (uri) => fetched(uri).then((b) => new TextDecoder().decode(b)),
  readBytes: (uri) => fetched(uri).then((b) => new Uint8Array(b)),
};

/** A position range as people read it (1-based, inclusive, thousands separated). */
const range = (system: string, start: number, end: number) =>
  `${system}:${(start + 1).toLocaleString('en-US')}-${end.toLocaleString('en-US')}`;

/** A hit as a line: its label, the field it matched, its place and its item's address (none: a range). */
export function hitLine(h: NameHit): string {
  return `${h.label}  (${h.field})  ${range(h.system, h.start, h.end)}  ${hitAddress(h) ?? '(no item: a range)'}`;
}

/** The view's region of the TP53 section once it is placed after a move (a few frames), in chr17's positions. */
async function placed(viewer: SeqViewer, before: string | undefined, system: string): Promise<string> {
  const read = () => {
    const r = viewer.view.get(Tp53Section);
    return r ? range(system, r.start + Tp53Region[0], r.end + Tp53Region[0]) : undefined;
  };
  for (let k = 0; k < 120; k++) {
    const now = read();
    if (now !== undefined && now !== before) return now;
    await new Promise<void>((done) => DomHost.requestAnimationFrame(() => done()));
  }
  return read() ?? 'nowhere (the section is not placed)';
}

/**
 * Finds `text` in the view's file (exact names first, then by prefix), focuses the first hit's extent, and says what
 * was searched, what was found and where the view went.
 */
export async function findAndFocus(
  viewer: SeqViewer,
  format: 'gff3' | 'gtf',
  text: string,
  reader: AssetReader = PageReader,
): Promise<string> {
  const fields = format === 'gff3' ? 'Name, ID, Alias' : 'gene_name, gene_id, transcript_name, transcript_id';
  const searched = `Searched “${text}” in ${fields} (exact names first, then names starting with it)`;
  const hits = await findTp53(format, text, reader, new AbortController().signal);
  if (!hits.length) return `${searched}: nothing found`;
  const first = hits[0];
  const r = viewer.view.get(Tp53Section);
  const before = r ? range(first.system, r.start + Tp53Region[0], r.end + Tp53Region[0]) : undefined;
  viewer.runCommand('address.focus', { address: hitRange(first) });
  const went = await placed(viewer, before, first.system);
  return [
    `${searched}: ${hits.length} hit${hits.length === 1 ? '' : 's'}`,
    `Found: ${first.label} at ${range(first.system, first.start, first.end)} (its ${first.field})`,
    `The view went to ${went}: ${first.label}’s extent`,
    '',
    'Every hit:',
    ...hits.map(hitLine),
  ].join('\n');
}

/**
 * A TP53 view with every transcript lane shown (its Transcripts track and its section at `height: 'content'`, as the
 * Tooltips' By data does): the canvas and its SVG twin draw the same lanes, no inner scroll on the canvas to find the
 * transcript a name leads to.
 */
export function everyLane(view: ExampleView): ExampleView {
  return {
    ...view,
    build: () =>
      withParams(withParams(view.build(), 'track', 'transcripts', { height: 'content' }), 'section', 'chr17', {
        height: 'content',
      }),
  };
}
