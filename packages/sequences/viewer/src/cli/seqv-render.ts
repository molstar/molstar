#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `seqv-render <document> [--view <ref>] [--width 1200] [--scale 1] [-o <file>] [--png]`
// (seqspec:viewer#the-render-backend): a SeqViewSpec document, JSON or MessagePack, rendered headless by `renderSvg`
// and written as SVG (stdout, or `-o`), or as PNG (`--png`, or `-o` ending in .png): the SVG rasterized by Chrome,
// which is how the CLI finds one (`playwright`, the Chrome the browser tests use); only this flag needs it. Sources are
// files, never the network: a source's `uri` or `hostFile` is read relative to `--root` (default the document's
// folder), and never outside it (`..`, an absolute path or a link out is refused). Text widths are the committed
// tables' (`metricsMeasurer`). Exit 0: written; 1: nothing written (the document does not load, the render or `--png`
// failed); 2: usage or an unreadable file; 3: written short of data, a source that could not be read, an issue of the
// document's that leaves something out (a warning or an error: a streamed source, which the CLI does not stream, is
// `svs.provider.unavailable`) or sources still loading after `--timeout` ms (default 120,000), each said on stderr, so
// a script can tell a figure with holes from none. Runs with Node's type stripping.

import { readFileSync, realpathSync, writeFileSync } from 'node:fs';
import { basename, dirname, isAbsolute, join, relative, resolve, sep } from 'node:path';
import { deflateSync } from 'node:zlib';
import { type Loaders, readDocument } from '@molstar/seq-core/format';
import { DefaultFontMetrics, metricsMeasurer, renderSvg } from '../svg/index.js';
import { rasterizeSvg } from './rasterize.js';

const FLAGS = new Set(['--overlay', '--nav', '--background', '--png']);
const OPTIONS = new Set([
  '--view',
  '--width',
  '--scale',
  '--sections',
  '--title',
  '--gap',
  '--root',
  '--timeout',
  '-o',
]);
const args = process.argv.slice(2);
const flags = new Set<string>();
const opt = new Map<string, string>();
const files: string[] = [];
let bad = args.length === 0;
for (let i = 0; i < args.length; i++) {
  const a = args[i];
  if (FLAGS.has(a)) flags.add(a);
  else if (OPTIONS.has(a) && i + 1 < args.length) opt.set(a, args[++i]);
  else if (a.startsWith('-')) bad = true;
  else files.push(a);
}
const num = (name: string, fallback: number) => (opt.has(name) ? Number(opt.get(name)) : fallback);
const [width, scale, gap, timeout] = [
  num('--width', 1200),
  num('--scale', 1),
  num('--gap', 8),
  num('--timeout', 120_000),
];
if (
  bad ||
  files.length !== 1 ||
  ![width, scale, gap, timeout].every((n) => Number.isFinite(n) && n >= 0) ||
  width === 0 ||
  scale === 0
) {
  console.error(
    'usage: seqv-render <document.svsj|.svsm> [--view <ref>] [--width 1200] [--scale 1] [--sections a,b] [--overlay] [--nav]\n' +
      '                   [--background] [--title <text>] [--gap 8] [--root <dir>] [--timeout 120000] [--png] [-o <file>]',
  );
  process.exit(2);
}

let bytes: Uint8Array;
try {
  bytes = readFileSync(files[0]);
} catch (e) {
  console.error(`seqv-render: ${e instanceof Error ? e.message : e}`);
  process.exit(2);
}
const root = resolve(opt.get('--root') ?? dirname(resolve(files[0])));
/** A source's file inside `--root` (links resolved, so none leads out); else it is refused, as the network is. */
const inside = (where: string) => {
  // A path that is not there: its nearest existing folder resolved, the rest appended (a missing file is "not
  // found", not "outside", under a root whose path holds a link).
  const real = (p: string): string => {
    try {
      return realpathSync(p);
    } catch {
      const up = dirname(p);
      return up === p ? p : join(real(up), basename(p));
    }
  };
  const p = real(resolve(root, where));
  const r = relative(real(root), p);
  if (r === '..' || r.startsWith(`..${sep}`) || isAbsolute(r))
    throw new Error('seqv-render reads no file outside --root');
  return p;
};
let failed = false;
const loaders: Loaders = {
  async load(request) {
    const where = request.hostFile ?? request.uri ?? '';
    try {
      if (/^[a-z][a-z0-9+.-]*:\/\//i.test(where)) throw new Error('seqv-render does not use the network');
      return readFileSync(inside(where));
    } catch (e) {
      failed = true;
      console.error(`seqv-render: ${where}: ${e instanceof Error ? e.message : e}`);
      throw e;
    }
  },
};

try {
  const { doc, issues } = readDocument(bytes);
  if (!doc) throw new Error(`not a document: ${issues.map((i) => i.message).join('; ')}`);
  const svg = await renderSvg(doc, {
    width,
    scale,
    gap,
    view: opt.get('--view'),
    sections: opt.get('--sections')?.split(','),
    overlay: flags.has('--overlay'),
    nav: flags.has('--nav'),
    background: flags.has('--background'),
    title: opt.get('--title'),
    loaders,
    ports: { measurer: metricsMeasurer(DefaultFontMetrics), deflate: (b) => deflateSync(b) },
    timeout,
    unsettled: (inFlight) => {
      failed = true;
      console.error(
        `seqv-render: sources still loading after ${timeout} ms (${inFlight} stream requests in flight); drawn as they stand`,
      );
    },
    // A warning or an error is data left out (an info is a note on how valid input was read).
    issues: (list) => {
      for (const i of list)
        if (i.severity !== 'info') {
          failed = true;
          console.error(`seqv-render: ${i.severity} ${i.code}${i.ref ? ` (${i.ref})` : ''}: ${i.message}`);
        }
    },
  });
  const out = opt.get('-o');
  const png = flags.has('--png') || out?.endsWith('.png');
  const data = png ? await rasterizeSvg(svg) : svg;
  if (out) writeFileSync(out, data);
  else process.stdout.write(data);
  if (failed) process.exitCode = 3;
} catch (e) {
  console.error(`seqv-render: ${e instanceof Error ? e.message : e}`);
  process.exit(1);
}
