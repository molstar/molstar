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
//     --entry <text> [--pages all | 0,3-5] [--from <downloaded file>] [--name <file name>] [--cut <seq>:<start>-<end>]
//
// `--cut` keeps a tab-line file's lines over a region (1-based, inclusive: a GFF3's or GTF's TP53), with its header lines
// (a `##sequence-region` of that sequence only), read as a stream (a whole genome's annotation inflates past what a
// string holds), the cut gzipped and recorded whole: the row's `derived` says the cut.
//
// `--pages` defaults to `all` (a small file whole). For a file read by region, run the runtime once over the download
// through `recordingRanges(bytes)` and pass the pages it noted. The licence is checked by hand when the file is fetched
// and must be one of `src/files/record-trim.ts`' `AllowedLicences`.

import { createHash } from 'node:crypto';
import { createReadStream, existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { basename, dirname, join } from 'node:path';
import { createInterface } from 'node:readline';
import { fileURLToPath } from 'node:url';
import { createGunzip, gzipSync } from 'node:zlib';
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

/** A region a cut keeps: a sequence's name and its 1-based inclusive positions. */
export interface Cut {
  readonly seq: string;
  readonly start: number;
  readonly end: number;
}

/** `chr17:7668402-7687550` as a cut. */
export function parseCut(spec: string): Cut {
  const m = /^([^:\s]+):(\d+)-(\d+)$/.exec(spec.replaceAll(',', ''));
  if (!m || Number(m[3]) < Number(m[2])) throw new Error(`--cut: ${spec}`);
  return { seq: m[1], start: Number(m[2]), end: Number(m[3]) };
}

/**
 * The lines a cut keeps: the header lines before the first data line (`#`; a `##sequence-region` of another sequence
 * left out), and the data lines of its sequence (column 1) whose columns 4 and 5 (1-based, inclusive) meet it; nothing
 * after `##FASTA`.
 */
export async function cutLines(lines: AsyncIterable<string> | Iterable<string>, cut: Cut): Promise<string[]> {
  const out: string[] = [];
  let header = true;
  for await (const line of lines) {
    if (line.startsWith('##FASTA')) break;
    if (line.startsWith('#')) {
      if (!header) continue;
      if (line.startsWith('##sequence-region') && line.slice(17).trim().split(/\s+/)[0] !== cut.seq) continue;
      out.push(line);
      continue;
    }
    if (!line.trim()) continue;
    header = false;
    const c = line.split('\t');
    if (c[0] === cut.seq && Number(c[3]) <= cut.end && Number(c[4]) >= cut.start) out.push(line);
  }
  return out;
}

/** A downloaded file's lines over a region, gzipped (a gzip download inflated as a stream). */
async function cutFile(from: string, cut: Cut): Promise<Uint8Array> {
  const head = readFileSync(from).subarray(0, 2);
  const raw = createReadStream(from);
  const stream = head[0] === 0x1f && head[1] === 0x8b ? raw.pipe(createGunzip()) : raw;
  const lines = await cutLines(createInterface({ input: stream, crlfDelay: Number.POSITIVE_INFINITY }), cut);
  return new Uint8Array(gzipSync(`${lines.join('\n')}\n`, { level: 9 }));
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
  /** Keep only a region's lines (a downloaded file, `from`), gzipped. */
  readonly cut?: Cut;
}

/** A dataset's or file's name that stays in its folder: one path segment, no `..`, no separator, nothing hidden. */
export function isPlainName(name: string): boolean {
  return /^[A-Za-z0-9][A-Za-z0-9._-]*$/.test(name) && !name.includes('..');
}

/** Records one file: the recording written, the dataset's `sources.json` row in place of the file's earlier one. */
export async function recordPages(job: RecordJob, fixtures: string): Promise<string> {
  // refused before anything is fetched or written: a dataset or name never leaves the fixtures folder
  if (!isPlainName(job.dataset)) throw new Error(`The dataset “${job.dataset}” is not a plain folder name`);
  if (job.name !== undefined && !isPlainName(job.name))
    throw new Error(`The name “${job.name}” is not a plain file name`);
  let bytes: Uint8Array;
  let downloaded: Uint8Array | undefined;
  if (job.cut) {
    if (!job.from) throw new Error('--cut reads a downloaded file: give it with --from');
    downloaded = new Uint8Array(readFileSync(job.from));
    bytes = await cutFile(job.from, job.cut);
  } else if (job.from) bytes = new Uint8Array(readFileSync(job.from));
  else {
    const res = await fetch(job.url, { headers: { 'User-Agent': 'fixture recorder' } });
    if (!res.ok) throw new Error(`${job.url}: ${res.status}`);
    bytes = new Uint8Array(await res.arrayBuffer());
  }
  const name = job.name ?? basename(new URL(job.url).pathname);
  if (!isPlainName(name)) throw new Error(`The name “${name}” is not a plain file name: give one with --name`);
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
    derived: job.cut
      ? `the lines over ${job.cut.seq}:${job.cut.start}-${job.cut.end} (1-based, inclusive) and the header lines (a ##sequence-region of ${job.cut.seq} only), gzipped (scripts/record-pages.ts --cut)`
      : 'none',
    bytesDownloaded: (downloaded ?? bytes).byteLength,
    bytesFixture: Buffer.byteLength(text),
    sha256Download: sha256(downloaded ?? bytes),
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
      'usage: record-pages.ts <dataset> <url> --licence <id> --attribution <text> --entry <text> [--pages all|0,3-5] [--from <file>] [--name <file>] [--cut <seq>:<start>-<end>]',
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
        ...(opt('cut') ? { cut: parseCut(opt('cut')!) } : {}),
      },
      fixtures,
    ),
  );
}
