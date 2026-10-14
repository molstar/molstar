/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `.svsx` fixture's source (SX1): a document that reads a FASTA and an annotation file by absolute URI, and the
// fetch that serves them offline. `fixtures/kinase.svsx` is `bundleUris` and `writeSvsx` over it (make-svsx-fixture.ts,
// run by hand); svsx.test.ts reads that file back; the Seq*/Cases/Format story opens it and bundles this document again.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '../../src/format/index.js';

export const FASTA_URI = 'https://example.org/kinase/abl.fa';
export const SITES_URI = 'https://example.org/kinase/sites.tsv';

/** URI → file text. */
export const SVSX_FILES: Readonly<Record<string, string>> = {
  [FASTA_URI]:
    '>abl-fragment ABL kinase, a short excerpt\nGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNS\n',
  [SITES_URI]: 'start\tstop\tlabel\n16\t20\tATP site\n47\t49\tgatekeeper\n88\t90\tSH3 contact\n',
};

/** The document with absolute URIs: a sequence from the FASTA, ranges from the TSV, one view. */
export function absoluteDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const src = b.data.source({ ref: 'fa', uri: FASTA_URI, format: 'fasta' });
  const seq = b.data.sequence({ ref: 'abl', source_ref: src, label: 'ABL kinase excerpt' });
  b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: seq,
    uri: SITES_URI,
    format: 'tsv',
  } as never);
  const v = b.view({ ref: 'v', label: 'ABL excerpt' });
  const t = v.block({ sequence: seq }).track({ label: 'Sequence' });
  t.letters();
  return b.getDocument();
}

/** A `fetch` over the fixture files: a `Response` per URI, 404 for any other; `requests` lists what it was asked. */
export function fixtureFetch(): { fetch: (uri: string) => Promise<Response>; requests: string[] } {
  const requests: string[] = [];
  return {
    requests,
    async fetch(uri) {
      requests.push(uri);
      const text = SVSX_FILES[uri];
      return text === undefined ? new Response('not found', { status: 404 }) : new Response(text);
    },
  };
}
