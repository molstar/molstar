/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records hg38's cytobands for the ideogram navigator's demos and tests (ST5): UCSC's cytoBandIdeo for the primary
// chromosomes (chr1..22, X, Y), as `_test/fixtures/hg38-cytobands.json`. Run by hand (`node scripts/record-cytobands.ts`),
// never in tests: it needs the network (curl, since Node's fetch is refused in some sandboxes) and overwrites the file.
// Rows are UCSC's, 0-based and half-open: [start, stop, band name, Giemsa stain].

import { execFileSync } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';

const URL_ = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz';
const primary = [...Array.from({ length: 22 }, (_, k) => `chr${k + 1}`), 'chrX', 'chrY'];
const text = gunzipSync(execFileSync('curl', ['-s', '-f', '--retry', '3', URL_], { maxBuffer: 1 << 26 })).toString();
const chromosomes: Record<string, [number, number, string, string][]> = Object.fromEntries(primary.map((k) => [k, []]));
for (const line of text.split('\n')) {
  const [chrom, start, stop, name, stain] = line.split('\t');
  if (chrom in chromosomes) chromosomes[chrom].push([Number(start), Number(stop), name, stain]);
}
const out = new URL('../_test/fixtures/hg38-cytobands.json', import.meta.url);
mkdirSync(new URL('./', out), { recursive: true });
// One chromosome a line: small diffs, readable.
const rows = Object.entries(chromosomes).map(([k, v]) => `${JSON.stringify(k)}:${JSON.stringify(v)}`);
writeFileSync(out, `{"source":${JSON.stringify(URL_)},"assembly":"hg38","chromosomes":{\n${rows.join(',\n')}\n}}\n`);
console.log(
  `${Object.keys(chromosomes).length} chromosomes, ${Object.values(chromosomes).flat().length} bands -> ${out.pathname}`,
);
