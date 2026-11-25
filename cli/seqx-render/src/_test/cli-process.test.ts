/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

// The built CLI in a process (`npx tsc -b` first): its files, its exit codes, and the cold-start target (a dataset
// view to SVG under 1 s, from the process's start to the written file; the PNG, which loads the rasterizer and the
// fonts, within the same bound).

import { execFile } from 'node:child_process';
import { existsSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { promisify } from 'node:util';
import { itWithin } from '@molstar/test-harness/timing';
import { afterAll, describe, expect, it } from 'vitest';
import { decodePng, JsonText, p53Tracks, scratch } from './fixture.js';

const Bin = fileURLToPath(new URL('../../lib/seqx-render.js', import.meta.url));
const exec = promisify(execFile);

const dir = scratch();
afterAll(dir.dispose);
const spec = dir.file('p53.json', JsonText(p53Tracks()));

/** The built CLI's run: its exit code and streams, never throwing. */
async function seqxRender(...args: string[]) {
  if (!existsSync(Bin)) throw new Error(`The CLI is not built (${Bin}): run npx tsc -b first`);
  try {
    const r = await exec(process.execPath, [Bin, ...args]);
    return { code: 0, out: r.stdout, err: r.stderr };
  } catch (e) {
    const x = e as { code: number; stdout: string; stderr: string };
    return { code: x.code, out: x.stdout, err: x.stderr };
  }
}

describe('seqx-render (built, in a process)', () => {
  it('writes an SVG and a PNG of p53’s tracks, exit 0', async () => {
    const svg = join(dir.dir, 'p53.svg');
    const png = join(dir.dir, 'p53.png');
    expect(await seqxRender(spec, '-o', svg, '--width', '640')).toEqual({ code: 0, out: '', err: '' });
    expect(await seqxRender(spec, '-o', png, '--width', '640', '--theme', 'dark', '--scale', '2')).toEqual({
      code: 0,
      out: '',
      err: '',
    });
    expect(readFileSync(svg, 'utf8')).toMatch(/^<svg /);
    const img = decodePng(readFileSync(png));
    expect([img.width, img.at(0, 0)]).toEqual([1280, [0x12, 0x14, 0x19, 255]]);
  });

  it('exits 2 with its usage, 1 with the spec’s issues, 4 for a missing file', async () => {
    const usage = await seqxRender(spec);
    expect(usage.code).toBe(2);
    expect(usage.err).toContain('Usage: seqx-render');
    const bad = await seqxRender(dir.file('bad.json', '{'), '-o', join(dir.dir, 'bad.svg'));
    expect(bad.code).toBe(1);
    expect(bad.err).toMatch(/^error repr\.codec\.json\.syntax/);
    expect((await seqxRender(join(dir.dir, 'none.json'), '-o', join(dir.dir, 'none.svg'))).code).toBe(4);
  });

  itWithin(
    'cold start: p53’s tracks to a written SVG and PNG, process start to file',
    async (within) => {
      const svg = join(dir.dir, 'cold.svg');
      const png = join(dir.dir, 'cold.png');
      // each is a fresh process: nothing is warm but the disk cache
      expect((await within(1000, () => seqxRender(spec, '-o', svg))).code).toBe(0);
      expect((await within(1000, () => seqxRender(spec, '-o', png))).code).toBe(0);
      expect(existsSync(svg) && existsSync(png)).toBe(true);
    },
    30_000,
  );
});
