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
//   node --experimental-strip-types scripts/record-structures.ts [ID ...]
//
// With ids, only those entries are fetched and their rows replaced, the others' kept. The entries: 1EHZ (tRNA-Phe,
// modified nucleotides), 1FIN (CDK2 and cyclin A, an ATP on each CDK2), 4HHB (human deoxyhaemoglobin, hemes and
// phosphates), 3SGJ (an IgG1 Fc with its N-glycans in `pdbx_branch_scheme`), 2OCJ (p53's core domain, four chains, a
// zinc each: p53's structure) and 3AVE (a human IgG1 Fc, fucosylated, its N297 glycans: the Antibody's Fc glycan), all
// CC0-1.0.

import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
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
  { id: '2OCJ', entry: 'PDB 2OCJ, the human p53 core domain, four chains (Wang, Rosengarth, Luecke 2007), mmCIF' },
  { id: '3AVE', entry: 'PDB 3AVE, a fucosylated human IgG1 Fc with its N297 glycans (Matsumiya et al. 2007), mmCIF' },
];

const sha256 = (b: Uint8Array) => createHash('sha256').update(b).digest('hex');
const recorded = new Date().toISOString().slice(0, 10);
mkdirSync(out, { recursive: true });
const asked = process.argv.slice(2);
const unknown = asked.filter((id) => !Entries.some((e) => e.id === id));
if (unknown.length) throw new Error(`Unknown entries: ${unknown.join(', ')}`);
const sourcesFile = join(out, 'sources.json');
const kept: { file: string }[] =
  asked.length && existsSync(sourcesFile) ? JSON.parse(readFileSync(sourcesFile, 'utf8')) : [];
const rows: { file: string }[] = [];
for (const { id, entry } of Entries.filter((e) => !asked.length || asked.includes(e.id))) {
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
// in the entries' order: a re-recorded row replaces its own, the others kept as they were
const all = Entries.map(
  ({ id }) =>
    rows.find((r) => r.file === `structures/${id}.cif.gz`) ?? kept.find((r) => r.file === `structures/${id}.cif.gz`),
).filter((r) => r !== undefined);
writeFileSync(sourcesFile, `${JSON.stringify(all, null, 2)}\n`);
