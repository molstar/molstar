/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetReader, NameHit } from '@molstar/seqx-core/data';
import { findTp53, hitAddress, hitRange } from '@molstar/seqx-examples/genomes/gff';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// A name found in a GFF3 or GTF (EC-228), as the stories show it: the file's `findNames` (the GFF runtime's index, exact
// names first, then by prefix), its hits printed, the first focused (`address.focus` of its extent, which the hit
// carries: no tile read first). Until the locus box falls back to it, the story asks the index itself.

/** The page's own files by `fetch` (the fixtures served beside the stories). */
const fetched = (uri: string) =>
  fetch(uri).then((r) => (r.ok ? r.arrayBuffer() : Promise.reject(new Error(`${r.status}`))));
export const PageReader: AssetReader = {
  readText: (uri) => fetched(uri).then((b) => new TextDecoder().decode(b)),
  readBytes: (uri) => fetched(uri).then((b) => new Uint8Array(b)),
};

/** A hit as a line: its label, the field it matched, its place and what it focuses. */
export function hitLine(h: NameHit): string {
  return `${h.label}  (${h.field})  ${h.system}:${h.start + 1}-${h.end}  ${hitAddress(h) ?? '(no item: a range)'}`;
}

/** Finds `text` in the view's file, focuses the first hit, and gives the hits as lines. */
export async function findAndFocus(
  viewer: SeqViewer,
  format: 'gff3' | 'gtf',
  text: string,
  reader: AssetReader = PageReader,
): Promise<string> {
  const hits = await findTp53(format, text, reader, new AbortController().signal);
  if (!hits.length) return `No name matches “${text}”`;
  const focused = viewer.runCommand('address.focus', { address: hitRange(hits[0]) });
  const head = `findNames(“${text}”, prefix): ${hits.length} hit${hits.length === 1 ? '' : 's'}`;
  return [head, ...hits.map(hitLine), `The first focused: ${focused}`].join('\n');
}
