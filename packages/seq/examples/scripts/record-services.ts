/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records the services showcase's fixtures (`fixtures/services/`), by hand, never in tests: PDBe's SIFTS mappings of
// PDB 2OCJ and 1FIN (as published), AlphaFold DB's per-residue confidence of P04637 (as published). UniProt's entries of P04637, P24941 and P20248 are the
// examples' own (`fixtures/p53/`, `fixtures/complex/`, recorded once). Each file's URL, date and licence go to
// `fixtures/services/sources.json`. Run from the package folder:
//
//   node --experimental-strip-types scripts/record-services.ts

import { createHash } from 'node:crypto';
import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const Folder = join(dirname(fileURLToPath(import.meta.url)), '..', 'fixtures', 'services');
const sha256 = (b: Uint8Array | string) => createHash('sha256').update(b).digest('hex');
const Recorded = new Date().toISOString().slice(0, 10);
// a generic client name: the AlphaFold DB CDN refuses the default one of command-line tools
const Headers = { 'User-Agent': 'Mozilla/5.0 (fixture recorder)', Accept: 'application/json' };

const Sifts = {
  licence: 'LicenseRef-EMBL-EBI-terms',
  attribution: 'PDBe, SIFTS (Dana et al. 2019), EMBL-EBI terms of use',
};
const Afdb = { licence: 'CC-BY-4.0', attribution: 'AlphaFold DB (Jumper et al. 2021; Varadi et al. 2024), CC BY 4.0' };

interface Job {
  readonly file: string;
  readonly entry: string;
  readonly url: string;
  readonly licence: string;
  readonly attribution: string;
  readonly trimmed: string;
}

const jobs: Job[] = [
  ...['2ocj', '1fin'].map((id) => ({
    file: `sifts-${id}.json`,
    entry: `PDBe SIFTS: PDB ${id.toUpperCase()}'s chains mapped to UniProt (mappings/uniprot/${id})`,
    url: `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/${id}`,
    ...Sifts,
    trimmed: 'none: as published, byte for byte',
  })),
  {
    file: 'AF-P04637-F1-confidence_v6.json',
    entry: 'AlphaFold DB AF-P04637-F1, model version 6, per-residue confidence (pLDDT)',
    url: 'https://alphafold.ebi.ac.uk/files/AF-P04637-F1-confidence_v6.json',
    ...Afdb,
    trimmed: 'none: as published, byte for byte',
  },
];

mkdirSync(Folder, { recursive: true });
const rows = [];
for (const j of jobs) {
  const res = await fetch(j.url, { headers: Headers });
  if (!res.ok) throw new Error(`${j.url}: ${res.status}`);
  const downloaded = new Uint8Array(await res.arrayBuffer());
  const raw = new TextDecoder().decode(downloaded);
  const text = raw;
  writeFileSync(join(Folder, j.file), text);
  const fixture = new TextEncoder().encode(text);
  rows.push({
    file: `services/${j.file}`,
    entry: j.entry,
    url: j.url,
    recorded: Recorded,
    licence: [j.licence],
    attribution: j.attribution,
    trimmed: j.trimmed,
    derived: 'none',
    bytesDownloaded: downloaded.byteLength,
    bytesFixture: fixture.byteLength,
    sha256Download: sha256(downloaded),
    sha256Fixture: sha256(fixture),
  });
}
writeFileSync(join(Folder, 'sources.json'), `${JSON.stringify(rows, null, 2)}\n`);
console.log(`recorded ${rows.length} files`);
