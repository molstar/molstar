/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The render CLI in process: SVG and PNG outputs, the spec's encodings, the exit codes, the issues printed one a
// line, and the headless asset policy (the spec's folder only, no network without a host).

import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { afterAll, describe, expect, it } from 'vitest';
import { ExitCode, formatIssue, run, Usage } from '../cli.js';
import { decodePng, JsonText, p53FromFile, p53Letters, P53, p53Tracks, scratch, SexprText } from './fixture.js';

const dir = scratch();
afterAll(dir.dispose);

async function cli(...args: string[]) {
  let out = '';
  let err = '';
  const code = await run(args, { out: (t) => (out += t), err: (t) => (err += t) });
  return { code, out, err };
}

describe('seqx-render', () => {
  const spec = dir.file('p53.json', JsonText(p53Tracks()));

  it('writes an SVG', async () => {
    const o = join(dir.dir, 'a.svg');
    const r = await cli(spec, '-o', o, '--width', '600');
    expect([r.code, r.err]).toEqual([0, '']);
    const svg = readFileSync(o, 'utf8');
    expect(svg).toMatch(/^<svg [^>]*width="600"/);
    expect(svg).toContain('data-seq-settled="true"');
  });

  it('writes a PNG at the width, scale and dpr asked', async () => {
    const o = join(dir.dir, 'a.png');
    expect((await cli(spec, '-o', o, '--width', '300', '--height', '120', '--scale', '2', '--dpr', '1.5')).code).toBe(
      0,
    );
    const img = decodePng(readFileSync(o));
    expect([img.width, img.height]).toEqual([900, 360]);
  });

  it('draws the theme asked: the PNG’s background is the theme’s', async () => {
    const o = join(dir.dir, 'dark.png');
    expect((await cli(spec, '-o', o, '--width', '200', '--height', '60', '--theme', 'dark')).code).toBe(0);
    expect(decodePng(readFileSync(o)).at(0, 0)).toEqual([0x12, 0x14, 0x19, 255]);
  });

  it('reads the spec as an S-expression or as JSON, told by its content', async () => {
    const a = join(dir.dir, 'json.svg');
    const b = join(dir.dir, 'sexpr.svg');
    expect((await cli(dir.file('letters.json', JsonText(p53Letters())), '-o', a)).code).toBe(0);
    const sx = await cli(dir.file('p53.seq', SexprText(p53Letters())), '-o', b);
    expect([sx.code, sx.err]).toEqual([0, '']);
    expect(readFileSync(b, 'utf8')).toBe(readFileSync(a, 'utf8'));
  });

  it('is deterministic: the same arguments give the same files', async () => {
    const a = join(dir.dir, 'd1.png');
    const b = join(dir.dir, 'd2.png');
    await cli(spec, '-o', a);
    await cli(spec, '-o', b);
    expect(readFileSync(a).equals(readFileSync(b))).toBe(true);
  });

  it('prints its usage for --help and for a wrong command line (exit 2)', async () => {
    expect(await cli('--help')).toEqual({ code: 0, out: Usage, err: '' });
    for (const args of [[], [spec], [spec, '-o', 'x.jpg'], [spec, '-o', 'x.svg', '--theme', 'blue'], ['--nope']]) {
      const r = await cli(...args);
      expect([args, r.code]).toEqual([args, ExitCode.Usage]);
      expect(r.err).toContain('Usage: seqx-render');
    }
    for (const bad of [
      ['--width', '0'],
      ['--height', 'tall'],
      ['--dpr', 'x'],
      ['--scale', '9', '--dpr', '9'],
    ])
      expect((await cli(spec, '-o', join(dir.dir, 'x.svg'), ...bad)).code, bad.join(' ')).toBe(ExitCode.Usage);
  });

  it('exits 4 for a spec it cannot read or an output it cannot write', async () => {
    expect((await cli(join(dir.dir, 'missing.json'), '-o', join(dir.dir, 'm.svg'))).code).toBe(ExitCode.Io);
    const r = await cli(spec, '-o', join(dir.dir, 'no-such-folder', 'm.svg'));
    expect(r.code).toBe(ExitCode.Io);
    expect(r.err).toContain('cannot write');
  });

  it('exits 1 for a spec that is not one, printing its issues a line each', async () => {
    const r = await cli(dir.file('bad.json', '{"not": '), '-o', join(dir.dir, 'bad.svg'));
    expect(r.code).toBe(ExitCode.Spec);
    expect(r.err).toMatch(/^error repr\.codec\.json\.syntax/);
    const wrong = await cli(dir.file('wrong.json', '{"kind":"spec"}'), '-o', join(dir.dir, 'w.svg'));
    expect(wrong.code).toBe(ExitCode.Spec);
    expect(
      wrong.err
        .split('\n')
        .filter(Boolean)
        .every((l) => /^(error|warning|info) [\w.]+/.test(l)),
    ).toBe(true);
  });

  it('exits 1 for a view the spec does not have, and renders another by --view', async () => {
    const r = await cli(spec, '-o', join(dir.dir, 'v.svg'), '--view', 'nope');
    expect(r.code).toBe(ExitCode.Spec);
    expect(r.err).toContain('seq.render.view');
    expect((await cli(spec, '-o', join(dir.dir, 'v2.svg'), '--view', '/views/main')).code).toBe(0);
  });

  it('formats an issue as severity, code, place and message', () => {
    expect(formatIssue({ code: 'a.b', severity: 'error', message: 'Wrong', address: '/data/x', hint: 'Fix it' })).toBe(
      'error a.b /data/x: Wrong (Fix it)',
    );
  });

  describe('assets (the headless policy)', () => {
    it('reads a file under the spec’s folder', async () => {
      mkdirSync(join(dir.dir, 'with-data', 'data'), { recursive: true });
      writeFileSync(join(dir.dir, 'with-data', 'data', 'p53.fa'), `>p53\n${P53}\n`);
      const s = dir.file('with-data/spec.json', JsonText(p53FromFile('data/p53.fa')));
      const o = join(dir.dir, 'with-data', 'fa.svg');
      const r = await cli(s, '-o', o, '--width', '400');
      expect([r.code, r.err]).toEqual([0, '']);
      expect(readFileSync(o, 'utf8')).toContain('data-seq-settled="true"');
    });

    it('refuses a file outside the spec’s folder and every URL without --allow-host', async () => {
      writeFileSync(join(dir.dir, 'outside.fa'), `>p53\n${P53}\n`);
      mkdirSync(join(dir.dir, 'inner'), { recursive: true });
      for (const uri of [
        '../outside.fa',
        join(dir.dir, 'outside.fa'),
        'https://example.org/p53.fa',
        'file:///etc/hosts',
      ]) {
        const s = dir.file('inner/spec.json', JsonText(p53FromFile(uri)));
        const o = join(dir.dir, 'inner', 'refused.svg');
        const r = await cli(s, '-o', o, '--timeout', '2000');
        expect([uri, r.code === ExitCode.Ok]).toEqual([uri, false]);
        expect(r.err, uri).toMatch(/not read|refused|seq\./);
        expect(() => readFileSync(o)).toThrow();
      }
    });
  });
});
