/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records the mmCIF runtime's fixtures, by hand, never in tests (examples.md, "Recorded by hand"): fetches each PDB
// entry's gzipped mmCIF from RCSB once, as served (not trimmed: the runtime's conformance suite reads the whole file),
// into `fixtures/structures/`, and writes `structures/sources.json` (URL, date, licence, attribution, size, SHA-256).
// Run it from the package folder:
//
//   node --experimental-strip-types scripts/record-structures.ts
//
// The entries: 1EHZ (tRNA-Phe, modified nucleotides), 1FIN (CDK2 and cyclin A, ATP and a magnesium), 4HHB (human
// deoxyhaemoglobin, hemes and phosphates) and 3SGJ (an IgG1 Fc with its N-glycans in `pdbx_branch_scheme`), all CC0-1.0.

import { createHash } from 'node:crypto';
import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';

const here = dirname(fileURLToPath(import.meta.url));
const out = join(here, '../fixtures/structures');

const Entries: readonly { readonly id: string; readonly entry: string }[] = [
  { id: '1EHZ', entry: 'PDB 1EHZ, yeast tRNA-Phe (Shi, Moore 2000), mmCIF' },
  { id: '1FIN', entry: 'PDB 1FIN, CDK2 and cyclin A with ATP (Jeffrey et al. 1995), mmCIF' },
  { id: '4HHB', entry: 'PDB 4HHB, human deoxyhaemoglobin (Fermi, Perutz, Shaanan, Fourme 1984), mmCIF' },
  { id: '3SGJ', entry: 'PDB 3SGJ, a human IgG1 Fc with its N-glycans (Ferrara et al. 2011), mmCIF' },
];

const sha256 = (b: Uint8Array) => createHash('sha256').update(b).digest('hex');
const recorded = new Date().toISOString().slice(0, 10);
mkdirSync(out, { recursive: true });
const rows = [];
for (const { id, entry } of Entries) {
  const url = `https://files.rcsb.org/download/${id}.cif.gz`;
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${url}: ${res.status}`);
  const bytes = new Uint8Array(await res.arrayBuffer());
  const text = gunzipSync(bytes);
  const file = `structures/${id}.cif.gz`;
  writeFileSync(join(out, `${id}.cif.gz`), bytes);
  rows.push({
    file,
    entry,
    url,
    recorded,
    licence: ['CC0-1.0'],
    attribution: `wwPDB, PDB ${id}, CC0`,
    trimmed: 'none: the file as served, gzipped',
    bytesFixture: bytes.length,
    bytesInflated: text.length,
    sha256Fixture: sha256(bytes),
    sha256Inflated: sha256(text),
  });
  console.log(`${file}: ${bytes.length} bytes (${text.length} inflated)`);
}
writeFileSync(join(out, 'sources.json'), `${JSON.stringify(rows, null, 2)}\n`);
