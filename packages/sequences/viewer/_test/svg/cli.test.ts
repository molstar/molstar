/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// `seqv-render` (SV1b; seqspec:viewer#the-render-backend) run as a process: a document file in, the SVG the library
// writes out, the options, the exit codes, file sources read relative to the root and the network refused, the
// executable bit, and `--png` through Chrome. The documents are case 14 and one made of
// packages/sequences/core/_test/fixtures/1iep.json; nothing is fetched. A process each: the heavy lane's, but for the
// first, the quick lane's one run of the CLI.

import { spawnSync } from 'node:child_process';
import {
  existsSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  rmdirSync,
  rmSync,
  statSync,
  symlinkSync,
  unlinkSync,
  writeFileSync,
} from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { deflateSync } from 'node:zlib';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { afterAll, describe, expect, it } from 'vitest';
import { case14 } from '../../../core/_test/format/cases.js';
import { IEP } from '../../../core/_test/format/load-fixtures.js';
import { DefaultFontMetrics, metricsMeasurer, renderSvg } from '../../src/svg/index.js';
import { decodePng } from './png.js';
import { all, parseXml } from './xml.js';

// The CLI runs built (`tsc -b`: the heavy lane builds first); the source carries the executable bit and the shebang.
const source = join(import.meta.dirname, '../../src/cli/seqv-render.ts');
const cli = join(import.meta.dirname, '../../lib/cli/seqv-render.js');
const dir = mkdtempSync(join(tmpdir(), 'seqv-render-'));
const made: string[] = [];
const write = (name: string, data: string | Uint8Array) => {
  const file = join(dir, name);
  writeFileSync(file, data);
  made.push(file);
  return file;
};
afterAll(() => {
  for (const f of made) unlinkSync(f);
  rmdirSync(dir);
});
const run = (...args: string[]) =>
  spawnSync(process.execPath, [cli, ...args], { encoding: 'utf8', maxBuffer: 1 << 28 });
const ports = { measurer: metricsMeasurer(DefaultFontMetrics), deflate: (b: Uint8Array) => deflateSync(b) };

function iepDoc() {
  const b = createSeqViewSpec({ title: '1IEP' });
  const iep = b.data.sequence({
    ref: 'iep',
    label: 'ABL1 kinase',
    alphabet: 'protein',
    value: IEP.sequence.slice(0, 60),
  });
  const v = b.view({ ref: 'v', label: '1IEP' }).block({ sequence: iep });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Sequence' }).letters();
  return b.getDocument();
}

// Each case starts Node (about 0.5 s, several times that under load).
describe('seqv-render', { timeout: 30_000 }, () => {
  const c14 = write('case14.json', JSON.stringify(case14()));

  // The quick lane's one run of the CLI: it waits for a `lib` (packages/sequences/AGENTS.md, Tests).
  it.skipIf(!existsSync(cli))('writes to stdout what renderSvg writes for the document', async () => {
    const r = run(c14, '--width', '900');
    expect(r.status).toBe(0);
    expect(r.stderr).toBe('');
    expect(r.stdout).toBe(await renderSvg(case14(), { width: 900, scale: 1, gap: 8, ports }));
    const svg = parseXml(r.stdout);
    expect(svg.attrs.width).toBe('900');
    expect(all(svg, 'text').length).toBeGreaterThan(20);
    expect(r.stdout).not.toMatch(/NaN|undefined/);
  });

  it('a document made of the vendored 1IEP fixture: -o, the scale, the overlay, the navigator, a background and a title', {
    tags: ['heavy'],
  }, () => {
    const doc = write('iep.json', JSON.stringify(iepDoc()));
    const out = join(dir, 'iep.svg');
    made.push(out);
    const r = run(
      doc,
      '--width',
      '700',
      '--scale',
      '2',
      '--overlay',
      '--nav',
      '--background',
      '--title',
      '1IEP',
      '-o',
      out,
    );
    expect(r.status).toBe(0);
    expect(r.stdout).toBe('');
    const svg = parseXml(readFileSync(out, 'utf8'));
    expect([svg.attrs.width, svg.attrs.viewBox?.split(' ')[2]]).toEqual(['1400', '700']);
    expect(all(svg, 'title')[0].text).toBe('1IEP');
    expect(all(svg, 'rect').length).toBeGreaterThan(0);
    // 60 residues in 700 px: the letters are drawn, the ruler's numbers too.
    expect(all(svg, 'text').map((t) => t.text)).toContain('M');
  });

  it('--view and --sections pick what is drawn', { tags: ['heavy'] }, () => {
    const r = run(c14, '--view', 'gene-view', '--sections', 'no-such');
    expect(r.status).toBe(1);
    expect(r.stderr).toMatch(/no section "no-such"/);
    expect(run(c14, '--view', 'gene-view').status).toBe(0);
  });

  it('exit 2 for usage and an unreadable file, 1 for a document that does not load or render', {
    tags: ['heavy'],
  }, () => {
    expect(run().status).toBe(2);
    expect(run().stderr).toMatch(/^usage: seqv-render/);
    expect(run(c14, '--bogus').status).toBe(2);
    expect(run(c14, '--width').status).toBe(2);
    expect(run(c14, '--width', '0').status).toBe(2);
    expect(run(c14, '--scale', 'x').status).toBe(2);
    expect(run(c14, c14).status).toBe(2);
    const missing = run(join(dir, 'missing.json'));
    expect([missing.status, missing.stderr]).toEqual([2, expect.stringMatching(/ENOENT/)]);
    const junk = run(write('junk.json', '{ not json'));
    expect([junk.status, junk.stderr]).toEqual([1, expect.stringMatching(/seqv-render: not a document/)]);
  });

  it('sources are files relative to the document, or to --root; the network is refused', { tags: ['heavy'] }, () => {
    const b = createSeqViewSpec({ title: 'from files' });
    const file = b.data.source({ ref: 'f', host_file: 's.fasta', format: 'fasta' });
    const s = b.data.sequence({ ref: 's', label: 'S', alphabet: 'protein', source_ref: file, record: 's' });
    const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: s });
    v.track({ label: 'Sequence' }).letters();
    write('s.fasta', `>s\n${IEP.sequence.slice(0, 40)}\n`);
    const doc = write('files.json', JSON.stringify(b.getDocument()));
    const ok = run(doc, '--width', '600');
    expect(ok.status, ok.stderr).toBe(0);
    expect(all(parseXml(ok.stdout), 'text').map((t) => t.text)).toContain('M');
    // A source that cannot be read: the figure is written, short of it, and the exit says so.
    const holes = run(doc, '--root', join(dir, 'nowhere'));
    expect([holes.status, holes.stderr]).toEqual([3, expect.stringMatching(/seqv-render: s\.fasta: .*ENOENT/)]);
    expect(parseXml(holes.stdout).name).toBe('svg');
    const net = createSeqViewSpec({ title: 'from the network' });
    const remote = net.data.source({ ref: 'f', uri: 'https://example.org/s.fasta', format: 'fasta' });
    const n = net.data.sequence({ ref: 's', label: 'S', alphabet: 'protein', source_ref: remote, record: 's' });
    net.view({ ref: 'v', label: 'V' }).block({ sequence: n }).track({ label: 'Sequence' }).letters();
    const r = run(write('net.json', JSON.stringify(net.getDocument())));
    expect(r.status).toBe(3);
    expect(r.stderr).toMatch(/does not use the network/);
  });

  it('a source outside --root is refused, as the network is: `..`, an absolute path, a link out (VB-06)', {
    tags: ['heavy'],
  }, () => {
    const box = mkdtempSync(join(tmpdir(), 'seqv-render-root-'));
    try {
      mkdirSync(join(box, 'root'));
      mkdirSync(join(box, 'outside'));
      const secret = join(box, 'outside', 'secret.fasta');
      writeFileSync(secret, '>s\nWWWWWWWWWWWWWWWWWWWW\n');
      symlinkSync(secret, join(box, 'root', 'link.fasta'));
      writeFileSync(join(box, 'root', 'in.fasta'), '>s\nMKVLMKVL\n');
      const docOf = (where: string) => {
        const b = createSeqViewSpec({ title: 'confined' });
        const file = b.data.source({ ref: 'f', host_file: where, format: 'fasta' });
        const q = b.data.sequence({ ref: 's', label: 'S', alphabet: 'protein', source_ref: file, record: 's' });
        b.view({ ref: 'v', label: 'V' }).block({ sequence: q }).track({ label: 'Sequence' }).letters();
        const doc = join(box, 'root', 'doc.json');
        writeFileSync(doc, JSON.stringify(b.getDocument()));
        return doc;
      };
      const ok = run(docOf('in.fasta'));
      expect(ok.status, ok.stderr).toBe(0);
      for (const where of ['../outside/secret.fasta', secret, 'link.fasta']) {
        const r = run(docOf(where));
        expect([r.status, r.stderr]).toEqual([3, expect.stringMatching(/outside --root/)]);
        expect(all(parseXml(r.stdout), 'text').map((t) => t.text)).not.toContain('W');
      }
      // A file not there is not found, not outside, though the root's path holds a link (macOS's /var) (LB-35).
      const missing = run(docOf('missing.fasta'));
      expect([missing.status, missing.stderr]).toEqual([3, expect.stringMatching(/missing\.fasta: .*ENOENT/)]);
    } finally {
      rmSync(box, { recursive: true });
    }
  });

  it('a streamed source, which the CLI does not stream: written short of it, exit 3, the issue on stderr (LB-13)', {
    tags: ['heavy'],
  }, () => {
    const b = createSeqViewSpec({ title: 'streamed' });
    const aln = b.data.alignment({
      ref: 'lazy',
      row_count: 300,
      column_count: 1200,
      source: { uri: 'https://example.org/streamed-lazy', protocol: 'seqstar-tiles/1' },
    } as never);
    const v = b.view({ ref: 'rows', label: 'Lazy', height: 'box' });
    v.block({ alignment: aln });
    v.track({ label: 'Rows' }).rows();
    const r = run(write('streamed.json', JSON.stringify(b.getDocument())), '--timeout', '2000');
    expect([r.status, r.stderr]).toEqual([3, expect.stringMatching(/seqv-render: warning svs\.provider\.unavailable/)]);
    expect(parseXml(r.stdout).name).toBe('svg');
  });

  it('is executable, starts with a shebang, and package.json lists it as a bin', () => {
    // Windows has no execute bits (git keeps the mode in the index)
    if (process.platform !== 'win32') expect(statSync(source).mode & 0o111).toBe(0o111);
    expect(readFileSync(source, 'utf8')).toMatch(/^#!\/usr\/bin\/env node\r?\n/);
    const pkg = JSON.parse(readFileSync(join(import.meta.dirname, '../../package.json'), 'utf8'));
    expect(pkg.bin).toEqual({ 'seqv-render': './lib/cli/seqv-render.js' });
  });

  it('--png without playwright: exit 1, nothing written, and the message names what is missing', {
    tags: ['heavy'],
  }, () => {
    const out = join(dir, 'never.png');
    // A resolve hook that finds no playwright, as on a machine without it.
    const hook =
      'export async function resolve(s, c, next) { if (s === "playwright") throw Object.assign(new Error("Cannot find package \'playwright\'"), { code: "ERR_MODULE_NOT_FOUND" }); return next(s, c); }';
    const register = `import { register } from "node:module"; register(${JSON.stringify(`data:text/javascript,${encodeURIComponent(hook)}`)});`;
    const r = spawnSync(
      process.execPath,
      ['--import', `data:text/javascript,${encodeURIComponent(register)}`, cli, c14, '--png', '-o', out],
      { encoding: 'utf8' },
    );
    expect(r.status).toBe(1);
    expect(r.stderr).toMatch(
      /^seqv-render: --png needs the `playwright` package \(npm install playwright\) and Google Chrome: Cannot find package 'playwright'/,
    );
    expect(() => statSync(out)).toThrow();
  });

  it("--png rasterizes the SVG in Chrome to a PNG of the SVG's own size", { timeout: 90_000, tags: ['heavy'] }, () => {
    const out = join(dir, 'case14.png');
    made.push(out);
    const r = run(c14, '--width', '600', '--scale', '2', '--background', '--png', '-o', out);
    expect(r.status, r.stderr).toBe(0);
    const d = decodePng(readFileSync(out));
    expect([d.w, d.h]).toEqual([
      1200,
      Math.ceil(Number(parseXml(run(c14, '--width', '600', '--scale', '2').stdout).attrs.height)),
    ]);
  });
});
