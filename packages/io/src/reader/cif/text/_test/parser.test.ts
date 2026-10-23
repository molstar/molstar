/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { ChunkedBigString } from '@molstar/core/util/string-like';
import type { CifCategory, CifFile, CifFrame } from '../../data-model.js';
import { parseCifText } from '../parser.js';
import { parseCifTextFromString } from './string-parser.js';

/** The array behind a column or field: its `'array'` storage. */
function arrayOf(c: { storage?(): { kind: string; array?: ArrayLike<any> } | undefined }) {
  const s = c.storage?.();
  return s?.kind === 'array' ? s.array : void 0;
}

// The byte tokenizer against the string tokenizer it replaced (`string-parser.ts`): every accessor of every field of
// every row, the file structure, and the errors with their lines, on the repo's text CIF files and on the syntax corner
// cases. A large file is in the medium lane below.

/** U+FEFF, the byte order mark (Biome writes an escape of it as the character itself) */
const BOM = String.fromCharCode(0xfeff);

const root = path.resolve(import.meta.dirname, '../../../../../../..');
const Files = [
  'data/examples/1crn.cif',
  'data/examples/1cbs_updated.cif',
  'data/examples/1grm_updated.cif',
  'data/examples/1bna_confal_pyramids.cif',
  'data/examples/867861-core.cif',
  'data/examples/ihm-coarse-8zzh.cif',
  'data/examples/7qpd.fw2.cif',
  'packages/prism/viewer/src/sequence/molstar/_test/fixtures/structures/ins1.cif',
  'smoke/fixtures/tiny.dic',
];

/** Syntax corner cases: quotes, multi-line and triple-quoted values, save frames, `_import.get`, line ends, UTF-8. */
const Cases: Record<string, string> = {
  quotes: `data_q\n_a.x 'it''s' _a.y "a "b" c" _a.z 'end\n_a.w 'x'y'\n#\n`,
  multiline: `data_m\nloop_\n_a.x\n_a.y\n;line 1\nline 2\n;\nv2\n;\n;\n;x;y\n# comment\n;\n`,
  tripleQuotes: `data_t\n_a.x '''a 'b' ''c''' _a.y '''\n'''\n`,
  crlf: `data_c\r\nloop_\r\n_a.x\r\n_a.y\r\n1 2\r\n;multi\r\nline\r\n;\r\n3\r\n`,
  cr: `data_r\rloop_\r_a.x\r1\r2\r`,
  unicode: `data_u\n_a.x 'Jürgen Šťastný' _a.y ☃ _a.z "🧬 x"\nloop_\n_b.name\nÅngström\nμ\n.\n?\n`,
  missing: `data_k\nloop_\n_a.x\n_a.y\n. ?\n'.' '?'\n"" ''\n`,
  numbers: `data_n\nloop_\n_a.x\n-1.5e3\n+2\n1.\n.5\n12abc\n2147483648\n-0\n`,
  flat: `data_f\nloop_\n_a\n_b\n1 2\n3 4\n`,
  saveFrames: `data_d\n_d.title x\nsave_one\n_i.x 1\nloop_\n_j.a\n_j.b\n1 2\nsave_\nsave_two\n_import.get [{'save':orient_matrix  'file':templ_attr.cif}]\n_k.y 2\nsave_\n`,
  blocks: `data_1\n_a.x 1\ndata_2\n_a.x 2\nloop_\n_b.y\n1\n2\n3\n`,
  keywords: `data_k\nloop_\n_a.x\nDATA_x\nloop_\n_b.y\n1\n`,
  semicolonInLine: `data_s\n_a.x 1;2 _a.y ;x\n`,
  bom: `${BOM}data_bom\n_a.x 1\n`,
  // U+FEFF inside the file is text: kept at the start of a value
  bomInValues: `data_v\n_a.x ${BOM}abc _a.y "${BOM}"\nloop_\n_b.v\n${BOM}x\n${BOM}\n`,
  long: `data_l\nloop_\n_a.x\n_a.y\n${'x'.repeat(300)} ${'y'.repeat(254)}\n${'z'.repeat(255)} 1\n`,
  // values of 255+ bytes in a loop that grows: the ends' map is re-keyed
  longGrow: `data_g\nloop_\n_a.x\n_a.y\n${Array.from({ length: 200 }, (_, i) => `${String(i).repeat(100 + i)} ${i}`).join('\n')}\n`,
  empty: ``,
  whitespace: `  \n\t\n`,
};

// values of 0 to 9 bytes drawn from 600, repeated: the string cache's hits, misses and slot collisions
Cases.cache = (() => {
  let seed = 11;
  const next = () => (seed = (Math.imul(seed, 1103515245) + 12345) >>> 0);
  const words = Array.from({ length: 600 }, (_, i) =>
    Array.from({ length: i % 10 }, () => 'ABCabc01é'[next() % 9]).join(''),
  ).map((w) => (w ? (w.includes('é') ? `'${w}'` : w) : "''"));
  const rows = Array.from({ length: 5000 }, () => `${words[next() % 600]} ${words[next() % 37]}`);
  return `data_c\nloop_\n_a.x\n_a.y\n${rows.join('\n')}\n`;
})();

const Errors: Record<string, string> = {
  notMultiple: `data_e\nloop_\n_a.x\n_a.y\n1 2\n3\n`,
  noValue: `data_e\n_a.x\n_a.y 1\n`,
  nestedSave: `data_e\nsave_a\nsave_b\n`,
  unfinishedSave: `data_e\nsave_a\n_x.y 1\n`,
  dataInSave: `data_e\nsave_a\n_x.y 1\ndata_f\n`,
  unexpected: `data_e\n\n\nvalue\n`,
  notMultipleCr: `data_e\rloop_\r_a.x\r_a.y\r1 2\r3 4\r\r5\r_b.x 1\r`,
  notMultipleCrLf: `data_e\r\nloop_\r\n_a.x\r\n_a.y\r\n1 2\r\n;x\r\ny\r\n;\r\n5 6\r\n_b.x 1\r\n`,
  unexpectedCr: `data_e\r\r_a.x 1\r\rvalue\r`,
};

async function parse(task: ReturnType<typeof parseCifText>) {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  try {
    return await task.run();
  } finally {
    log.mockRestore();
  }
}

function expectSameCategory(a: CifCategory, b: CifCategory, where: string) {
  expect([a.name, a.rowCount, a.fieldNames], where).toEqual([b.name, b.rowCount, b.fieldNames]);
  for (const name of a.fieldNames) {
    const fa = a.getField(name)!;
    const fb = b.getField(name)!;
    expect(fa.rowCount, `${where}.${name}`).toBe(fb.rowCount);
    const n = fa.rowCount;
    for (let row = 0; row < n; row++) {
      const same =
        fa.str(row) === fb.str(row) &&
        Object.is(fa.int(row), fb.int(row)) &&
        Object.is(fa.float(row), fb.float(row)) &&
        fa.valueKind(row) === fb.valueKind(row) &&
        fa.areValuesEqual(row, 0) === fb.areValuesEqual(row, 0) &&
        fa.areValuesEqual(row, row > 0 ? row - 1 : 0) === fb.areValuesEqual(row, row > 0 ? row - 1 : 0);
      if (!same) {
        expect([row, fa.str(row), fa.int(row), fa.float(row), fa.valueKind(row)], `${where}.${name}`).toEqual([
          row,
          fb.str(row),
          fb.int(row),
          fb.float(row),
          fb.valueKind(row),
        ]);
      }
    }
    expect(fa.toStringArray({ start: 0, end: Math.min(n, 5) })).toEqual(
      fb.toStringArray({ start: 0, end: Math.min(n, 5) }),
    );
    expect(arrayOf(fa)).toBe(undefined);
    expect(fa.binaryEncoding).toBe(undefined);
  }
}

function expectSameFrame(a: CifFrame, b: CifFrame, where: string) {
  expect([a.header, a.categoryNames], where).toEqual([b.header, b.categoryNames]);
  for (const c of a.categoryNames) expectSameCategory(a.categories[c], b.categories[c], `${where}/${c}`);
}

export function expectSameFile(a: CifFile, b: CifFile, where: string) {
  expect(a.blocks.length, where).toBe(b.blocks.length);
  for (let i = 0; i < a.blocks.length; i++) {
    expectSameFrame(a.blocks[i], b.blocks[i], `${where}#${i}`);
    expect(a.blocks[i].saveFrames.length).toBe(b.blocks[i].saveFrames.length);
    for (let j = 0; j < a.blocks[i].saveFrames.length; j++) {
      expectSameFrame(a.blocks[i].saveFrames[j], b.blocks[i].saveFrames[j], `${where}#${i}/save${j}`);
    }
  }
}

/** Parsed from bytes and from the string, both equal to the string tokenizer's result. */
export async function expectSameAsStringParser(text: string, where: string) {
  const expected = await parse(parseCifTextFromString(text));
  for (const input of [new TextEncoder().encode(text), text]) {
    const actual = await parse(parseCifText(input));
    expect(actual.isError, where).toBe(expected.isError);
    if (actual.isError || expected.isError) {
      expect([actual.message, (actual as any).line], where).toEqual([expected.message, (expected as any).line]);
    } else {
      expectSameFile(actual.result, expected.result, where);
    }
  }
}

describe('text CIF byte tokenizer', () => {
  it('reads the repo text CIF files as the string tokenizer did', async () => {
    for (const f of Files) await expectSameAsStringParser(fs.readFileSync(path.join(root, f), 'utf8'), f);
  });

  it('reads the syntax corner cases as the string tokenizer did', async () => {
    for (const [name, text] of Object.entries(Cases)) {
      if (name === 'bom') continue;
      await expectSameAsStringParser(text, name);
    }
  });

  it('reports errors at the lines the string tokenizer did', async () => {
    for (const [name, text] of Object.entries(Errors)) {
      const r = await parse(parseCifTextFromString(text));
      expect(r.isError, name).toBe(true);
      await expectSameAsStringParser(text, name);
    }
  });

  it('decodes UTF-8 values', async () => {
    const r = await parse(parseCifText(new TextEncoder().encode(Cases.unicode)));
    if (r.isError) throw new Error(r.message);
    const a = r.result.blocks[0].categories['a'];
    expect([a.getField('x')!.str(0), a.getField('y')!.str(0), a.getField('z')!.str(0)]).toEqual([
      'Jürgen Šťastný',
      '☃',
      '🧬 x',
    ]);
    const b = r.result.blocks[0].categories['b'].getField('name')!;
    expect(b.toStringArray()).toEqual(['Ångström', 'μ', '', '']);
    expect(b.areValuesEqual(0, 1)).toBe(false);
  });

  it('skips a UTF-8 byte order mark', async () => {
    const r = await parse(parseCifText(new TextEncoder().encode(Cases.bom)));
    if (r.isError) throw new Error(r.message);
    expect(r.result.blocks[0].header).toBe('bom');
    expect(r.result.blocks[0].categories['a'].getField('x')!.int(0)).toBe(1);
  });

  it('reads a ChunkedBigString', async () => {
    const text = fs.readFileSync(path.join(root, 'data/examples/1cbs_updated.cif'), 'utf8');
    const expected = await parse(parseCifText(text));
    const actual = await parse(parseCifText(ChunkedBigString.fromString(text, 12)));
    if (expected.isError || actual.isError) throw new Error('parse failed');
    expectSameFile(actual.result, expected.result, '1cbs chunked');
  });
});

describe('text CIF byte tokenizer, large', { tags: ['medium'] }, () => {
  it('reads a large atom_site loop as the string tokenizer did', async () => {
    // 7qpd's atom_site rows 25 times: the loop grows past its first estimate and is trimmed
    const text = fs.readFileSync(path.join(root, 'data/examples/7qpd.fw2.cif'), 'utf8');
    const lines = text.split('\n');
    const rows = lines.filter((l) => l.startsWith('ATOM') || l.startsWith('HETATM'));
    const first = lines.indexOf(rows[0]);
    const big = [
      ...lines.slice(0, first),
      ...Array.from({ length: 25 }, () => rows).flat(),
      ...lines.slice(first + rows.length),
    ];
    await expectSameAsStringParser(big.join('\n'), '7qpd x25');
  });
});
