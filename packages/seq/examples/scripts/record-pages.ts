/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records a real file as pages (`src/recorded-pages.ts`), by hand, never in tests: the file fetched (or read from a file
// already downloaded, `--from`), the pages asked kept byte for byte in `fixtures/<dataset>/<name>.pages.json`, and its
// row (URL, date, licence, attribution, sizes and checksums of the download and of the recording) put in the dataset's
// `sources.json` in place of the file's earlier row. Run from the package folder:
//
//   node --experimental-strip-types scripts/record-pages.ts <dataset> <url> --licence <id> --attribution <text> \
//     --entry <text> [--pages all | 0,3-5] [--from <downloaded file>] [--name <file name>]
//
// `--pages` defaults to `all` (a small file whole). For a file read by region, run the runtime once over the download
// through `recordingRanges(bytes)` and pass the pages it noted. The licence is checked by hand when the file is fetched
// and must be one of `src/files/record-trim.ts`' `AllowedLicences`.

import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { basename, dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { encodePages, PagesSuffix } from '../src/recorded-pages.ts';

const sha256 = (b: Uint8Array | string) => createHash('sha256').update(b).digest('hex');

/** `0,3-5` as page numbers. */
function parsePages(spec: string): number[] | 'all' {
  if (spec === 'all') return 'all';
  return spec.split(',').flatMap((part) => {
    const [a, b] = part.split('-').map(Number);
    if (!Number.isInteger(a) || (b !== undefined && !Number.isInteger(b))) throw new Error(`--pages: ${part}`);
    return Array.from({ length: (b ?? a) - a + 1 }, (_, i) => a + i);
  });
}

export interface RecordJob {
  readonly dataset: string;
  readonly url: string;
  readonly entry: string;
  readonly licence: string;
  readonly attribution: string;
  readonly pages: number[] | 'all';
  /** A file already downloaded (else the URL is fetched). */
  readonly from?: string;
  /** The recorded file's name (default: the URL's last segment). */
  readonly name?: string;
}

/** Records one file: the recording written, the dataset's `sources.json` row in place of the file's earlier one. */
export async function recordPages(job: RecordJob, fixtures: string): Promise<string> {
  let bytes: Uint8Array;
  if (job.from) bytes = new Uint8Array(readFileSync(job.from));
  else {
    const res = await fetch(job.url, { headers: { 'User-Agent': 'fixture recorder' } });
    if (!res.ok) throw new Error(`${job.url}: ${res.status}`);
    bytes = new Uint8Array(await res.arrayBuffer());
  }
  const name = job.name ?? basename(new URL(job.url).pathname);
  const recording = encodePages(name, bytes, job.pages);
  const text = `${JSON.stringify(recording)}\n`;
  const folder = join(fixtures, job.dataset);
  mkdirSync(folder, { recursive: true });
  const file = `${name}${PagesSuffix}`;
  writeFileSync(join(folder, file), text);
  const kept = Object.keys(recording.pages).length;
  const total = Math.ceil(bytes.length / recording.pageSize);
  const row = {
    file: `${job.dataset}/${file}`,
    entry: job.entry,
    url: job.url,
    recorded: new Date().toISOString().slice(0, 10),
    licence: [job.licence],
    attribution: job.attribution,
    trimmed:
      kept === total
        ? `none: every page (${total} of ${recording.pageSize} bytes), byte for byte, as recorded pages`
        : `pages ${Object.keys(recording.pages).join(', ')} of ${total} (${recording.pageSize} bytes each), byte for byte; a read of another page fails naming it`,
    derived: 'none',
    bytesDownloaded: bytes.byteLength,
    bytesFixture: Buffer.byteLength(text),
    sha256Download: sha256(bytes),
    sha256Fixture: sha256(text),
  };
  const sources = join(folder, 'sources.json');
  const rows = existsSync(sources) ? (JSON.parse(readFileSync(sources, 'utf8')) as { file: string }[]) : [];
  const at = rows.findIndex((r) => r.file === row.file);
  if (at < 0) rows.push(row);
  else rows[at] = row;
  writeFileSync(sources, `${JSON.stringify(rows, null, 2)}\n`);
  return `${row.file}: ${kept} of ${total} pages, ${row.bytesFixture} bytes`;
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const args = process.argv.slice(2);
  const opt = (k: string) => {
    const i = args.indexOf(`--${k}`);
    return i < 0 ? undefined : args[i + 1];
  };
  const [dataset, url] = args;
  const licence = opt('licence');
  const attribution = opt('attribution');
  const entry = opt('entry');
  if (!dataset || !url || dataset.startsWith('--') || !licence || !attribution || !entry) {
    console.error(
      'usage: record-pages.ts <dataset> <url> --licence <id> --attribution <text> --entry <text> [--pages all|0,3-5] [--from <file>] [--name <file>]',
    );
    process.exit(2);
  }
  const fixtures = join(dirname(fileURLToPath(import.meta.url)), '..', 'fixtures');
  console.log(
    await recordPages(
      {
        dataset,
        url,
        entry,
        licence,
        attribution,
        pages: parsePages(opt('pages') ?? 'all'),
        from: opt('from'),
        name: opt('name'),
      },
      fixtures,
    ),
  );
}
