/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mkdirSync, mkdtempSync, readFileSync, rmSync, symlinkSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { describe, expect, it } from 'vitest';
import { fixtureFile, seqExampleFixtures } from '../vite.js';

// The fixtures plugin: what it serves in development and emits in a build, and what it never serves.

describe('seqExampleFixtures', () => {
  it('names the fixture files under /fixtures/ and nothing else', () => {
    const file = fixtureFile('/fixtures/p53/tp53-gene.seq');
    expect(file?.endsWith('/examples/fixtures/p53/tp53-gene.seq')).toBe(true);
    expect(fixtureFile('/fixtures/p53/tp53-gene.seq?t=1')).toBe(file);
    expect(fixtureFile('/fixtures/p53/none.seq')).toBeUndefined();
    expect(fixtureFile('/fixtures/p53')).toBeUndefined();
    expect(fixtureFile('/other/p53/tp53-gene.seq')).toBeUndefined();
  });

  it('never serves a file above the folder, however the path is spelled', () => {
    expect(fixtureFile('/fixtures/../package.json')).toBeUndefined();
    expect(fixtureFile('/fixtures/%2e%2e/package.json')).toBeUndefined();
    expect(fixtureFile('/fixtures/p53/../../../package.json')).toBeUndefined();
    expect(fixtureFile('/fixtures/%')).toBeUndefined();
    expect(fixtureFile('/fixtures/p53/tp53-gene.seq%00')).toBeUndefined();
  });

  it('cuts the query and the fragment before decoding: a percent sign there is not the path’s', () => {
    const file = fixtureFile('/fixtures/p53/tp53-gene.seq');
    expect(fixtureFile('/fixtures/p53/tp53-gene.seq?x=%')).toBe(file);
    expect(fixtureFile('/fixtures/p53/tp53-gene.seq#%zz')).toBe(file);
    expect(fixtureFile('/fixtures/p53/tp53-gene.seq?/../..')).toBe(file);
  });

  it('refuses a link that leads out of the folder, by the real path', () => {
    const dir = mkdtempSync(join(tmpdir(), 'seqx-fixtures-'));
    try {
      mkdirSync(join(dir, 'root'));
      writeFileSync(join(dir, 'secret.txt'), 'secret');
      writeFileSync(join(dir, 'root', 'ok.txt'), 'fine');
      symlinkSync(join(dir, 'secret.txt'), join(dir, 'root', 'link.txt'));
      const root = join(dir, 'root');
      expect(fixtureFile('/fixtures/ok.txt', root)).toBe(join(root, 'ok.txt'));
      expect(fixtureFile('/fixtures/link.txt', root)).toBeUndefined();
    } finally {
      rmSync(dir, { recursive: true, force: true });
    }
  });

  it('serves a file in development and passes everything else on', () => {
    const plugin = seqExampleFixtures();
    let handler: (req: { url?: string; method?: string }, res: any, next: () => void) => void = () => {};
    plugin.configureServer({ middlewares: { use: (fn) => (handler = fn) } });
    const sent: unknown[] = [];
    const res = { setHeader: () => {}, end: (b: unknown) => sent.push(b) };
    let passed = 0;
    handler({ url: '/fixtures/p53/tp53-gene.seq' }, res, () => passed++);
    handler({ url: '/index.html' }, res, () => passed++);
    expect(passed).toBe(1);
    expect(sent.length).toBe(1);
    // HEAD is answered without a body; any other method is passed on
    handler({ url: '/fixtures/p53/tp53-gene.seq', method: 'HEAD' }, res, () => passed++);
    expect(sent[1]).toBe('');
    handler({ url: '/fixtures/p53/tp53-gene.seq', method: 'POST' }, res, () => passed++);
    handler({ url: '/fixtures/p53/tp53-gene.seq', method: 'PUT' }, res, () => passed++);
    expect(passed).toBe(3);
    expect(sent.length).toBe(2);
    expect(Buffer.from(sent[0] as Uint8Array).toString()).toContain('"seq"');
  });

  it('emits the files as assets of a build, under fixtures/', () => {
    const emitted: { fileName: string; source: Uint8Array }[] = [];
    seqExampleFixtures().generateBundle.call({ emitFile: (f) => emitted.push(f) });
    const names = emitted.map((e) => e.fileName);
    // every fixture under the folder, the Files dataset's and the recorded entry's too
    expect(names).toEqual(expect.arrayContaining(['fixtures/p53/gene-region.seq', 'fixtures/p53/tp53-gene.seq']));
    expect(names).toEqual(
      expect.arrayContaining([
        'fixtures/p53/P04637.json',
        'fixtures/files/domains.tsv',
        'fixtures/plasmid/L09137.2.seq',
      ]),
    );
    const gene = emitted.find((e) => e.fileName === 'fixtures/p53/tp53-gene.seq')!;
    expect(Buffer.from(gene.source).toString()).toBe(readFileSync(fixtureFile('/fixtures/p53/tp53-gene.seq')!, 'utf8'));
  });
});
