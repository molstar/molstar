/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetReader } from '@molstar/seqx-core/data';
import { findGene, GenomeSection, geneRange } from '@molstar/seqx-examples/genomes/tracks';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { PageReader } from '../../features/genome-formats/names.js';

// The Genomes dataset's find box (EC-228): a name looked up in the GFF3's name index, the first hit's extent focused,
// and what was searched, found and focused said.

export async function findGeneAndFocus(viewer: SeqViewer, text: string, reader: AssetReader = PageReader) {
  const hits = await findGene(text, reader, new AbortController().signal);
  const searched = `Searched “${text}” in Name and ID (exact names first, then names starting with it)`;
  if (!hits.length) return `${searched}: nothing found`;
  const first = hits[0];
  viewer.runCommand('address.focus', { address: geneRange(first) });
  return [
    `${searched}: ${hits.length} hit${hits.length === 1 ? '' : 's'}`,
    `Focused ${first.label} on ${first.system} from ${first.start + 1} to ${first.end} (${GenomeSection})`,
  ].join('\n');
}
