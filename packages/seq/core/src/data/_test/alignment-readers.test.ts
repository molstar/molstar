/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { CategoryColumn } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import {
  AlignmentFileRuntime,
  firstNumberOf,
  matchLeaves,
  readAlignedFasta,
  readNewick,
  readStockholm,
  structurePairs,
} from '../alignment-readers.js';
import { type CellsCodes, DefaultOpenLimits } from '../protocol.js';

const Stockholm = `# STOCKHOLM 1.0
#=GF ID   test
#=GS HBA_HUMAN/2-8 AC P69905
HBA_HUMAN/2-8   VLS-PAD
HBB_HUMAN/3-9   HLTPEEK
#=GC SS_cons    <<..>>.
#=GC RF         xx.xxxx

HBA_HUMAN/2-8   kt
HBB_HUMAN/3-9   sa
#=GC SS_cons    ..
#=GC RF         ..
//
`;

describe('the alignment readers (EC-38, EC-39, EC-40)', () => {
  it('Stockholm: rows joined across blocks, insert columns from RF, pairs from SS_cons, numbers from labels', () => {
    const a = readStockholm(Stockholm);
    expect(a.issues).toEqual([]);
    expect(a.columns).toBe(9);
    expect(a.rows.map((r) => [r.label, r.letters, r.first])).toEqual([
      ['HBA_HUMAN/2-8', 'VLS-PADkt', 2],
      ['HBB_HUMAN/3-9', 'HLTPEEKsa', 3],
    ]);
    expect(Array.from(a.inserts!)).toEqual([0, 0, 1, 0, 0, 0, 0, 1, 1]);
    expect(a.pairs).toEqual([
      [0, 5],
      [1, 4],
    ]);
  });

  it('Stockholm: `#=GS` lines are a row’s fields by tag (a tag twice joined), served on the cells’ tiles', async () => {
    const text = `# STOCKHOLM 1.0
#=GS A/1-4 OS Homo sapiens
#=GS A/1-4 DE Hemoglobin
#=GS A/1-4 DE subunit beta
#=GS B/2-5 AC P02088
#=GS B/2-5
A/1-4 AC-D
B/2-5 ACCD
//
`;
    const a = readStockholm(text);
    expect(a.rows.map((r) => r.fields)).toEqual([
      { OS: 'Homo sapiens', DE: 'Hemoglobin subunit beta' },
      { AC: 'P02088' },
    ]);
    expect(a.issues).toEqual([{ line: 6, message: 'a #=GS line without a name, a tag and its text' }]);
    // the runtime serves them on the rows' tiles, a column a tag, missing where a row has none
    const s = await AlignmentFileRuntime.open(
      { address: '/data/a', format: 'stockholm', value: text },
      {
        assets: { readText: () => Promise.reject(new Error('no')), readBytes: () => Promise.reject(new Error('no')) },
        signal: new AbortController().signal,
        limits: DefaultOpenLimits,
      },
    );
    const r = await s.read(
      { system: 'rows', start: 0, end: 2, level: 1, shapes: ['cells'] },
      new AbortController().signal,
    );
    const cells = r.ok ? (r.tables.cells as CellsCodes) : undefined;
    const field = (name: string) => {
      const c = cells?.fields?.find((f) => f.name === name)?.column as CategoryColumn | undefined;
      return [c?.get(0), c?.get(1)];
    };
    expect(field('OS')).toEqual(['Homo sapiens', null]);
    expect(field('AC')).toEqual([null, 'P02088']);
  });

  it('Stockholm: `#=GS` lines across blocks; one naming no row said once with its names, by the runtime', async () => {
    const text = `# STOCKHOLM 1.0
#=GS A OS Homo sapiens
A AC
B AC

#=GS B OS Mus musculus
#=GS X OS Nobody
#=GS Y OS Nobody
A -D
B CD
//
`;
    const a = readStockholm(text);
    expect(a.rows.map((r) => [r.letters, r.fields?.OS])).toEqual([
      ['AC-D', 'Homo sapiens'],
      ['ACCD', 'Mus musculus'],
    ]);
    expect(a.unknownAnnotated).toEqual(['X', 'Y']);
    const s = await AlignmentFileRuntime.open(
      { address: '/data/a', format: 'stockholm', value: text },
      {
        assets: { readText: () => Promise.reject(new Error('no')), readBytes: () => Promise.reject(new Error('no')) },
        signal: new AbortController().signal,
        limits: DefaultOpenLimits,
      },
    );
    expect(s.issues?.map((i) => [i.code, i.severity, i.message])).toEqual([
      ['seq.source.rows.annotation', 'info', '#=GS lines name 2 sequences with no row: X, Y'],
    ]);
  });

  it('aligned FASTA: lowercase-only columns are inserts; a row of another width is said and left out', () => {
    const a = readAlignedFasta('>a/10-14\nAC-gD\nE\n>b\nACcgDE\n>c\nAC\n');
    expect(a.rows.map((r) => r.letters)).toEqual(['AC-gDE', 'ACcgDE']);
    expect(a.rows[0].first).toBe(10);
    expect(Array.from(a.inserts!)).toEqual([0, 0, 1, 1, 0, 0]);
    expect(a.issues).toEqual([{ line: 6, message: expect.stringContaining('the row c is 2 columns wide') }]);
  });

  it('aligned FASTA: the most common width wins (a stray first row does not set it); a label twice said', () => {
    const a = readAlignedFasta('>x\nAC\n>a\nACDE\n>b\nAC-E\n>a\nAQDE\n');
    expect(a.columns).toBe(4);
    expect(a.rows.map((r) => r.label)).toEqual(['a', 'b', 'a']);
    expect(a.issues.map((i) => i.message)).toEqual([
      expect.stringContaining('the row x is 2 columns wide, not 4'),
      expect.stringContaining('the label a names more than one row'),
    ]);
  });

  it('aligned FASTA: a row with a character that is no letter or gap left out; `~` read as a gap', () => {
    const a = readAlignedFasta('>a\nAC~E\n>b\nA1DE\n>c\nA.*E\n');
    expect(a.rows.map((r) => [r.label, r.letters])).toEqual([
      ['a', 'AC-E'],
      ['c', 'A.*E'],
    ]);
    expect(a.issues).toEqual([{ line: 3, message: expect.stringContaining('the row b holds a character') }]);
  });

  it('dot-bracket pairs nest by bracket kind; pseudoknot letters pair with their lowercase', () => {
    expect(structurePairs('((.[[..))..]]').pairs).toEqual([
      [0, 8],
      [1, 7],
      [3, 12],
      [4, 11],
    ]);
    expect(structurePairs('A..<..a..>').pairs).toEqual([
      [0, 6],
      [3, 9],
    ]);
    expect(structurePairs('((.)').unpaired).toBe(1);
    expect(firstNumberOf('X/-3-9')).toBe(-3);
    expect(firstNumberOf('X')).toBe(1);
  });

  it('Newick: labels, quoted labels, lengths, comments; leaves in tree order; deep trees cost no stack', () => {
    const t = readNewick("((A:0.1,'B c':0.2)n1:0.3,[x]C:0.4)root;");
    expect(t.issues).toEqual([]);
    expect(t.leaves.map((l) => t.label[l])).toEqual(['A', 'B c', 'C']);
    expect(t.label[0]).toBe('root');
    expect(t.label[t.parent[t.leaves[0]]]).toBe('n1');
    expect(t.length[t.leaves[2]]).toBeCloseTo(0.4);
    const deep = `${'('.repeat(50_000)}a${')'.repeat(50_000)};`;
    expect(readNewick(deep).leaves.length).toBe(1);
    expect(readNewick('(a,b').issues[0].message).toContain('inside parentheses');
    expect(readNewick('(a,b));').issues[0].message).toContain('closes nothing');
  });

  it('a leaf matches its row by label, else by the label before its last `/` (EC-43)', () => {
    const t = readNewick('(HBA_HUMAN,(HBB_HUMAN/1-146,ZZZ));');
    const m = matchLeaves(t, [{ label: 'HBB_HUMAN/1-146' }, { label: 'HBA_HUMAN/1-141' }]);
    expect(Array.from(m)).toEqual([1, 0, -1]);
  });

  it('a row two leaves name goes to the first; the second leaf names none', () => {
    const t = readNewick('(a,(a/1-5,b));');
    expect(Array.from(matchLeaves(t, [{ label: 'a/1-5' }, { label: 'b' }]))).toEqual([0, -1, 1]);
  });

  it('fuzz: malformed text is reported, never thrown', () => {
    const rng = new PCG(11);
    const tree = "((A:0.1,'B c':0.2)n1:0.3,[x]C:0.4)root;";
    const fasta = '>a\nAC-GD\n>b\nACCGD\n';
    const chars = "()[],:;'>#=/ \n\t.-ACgx/0123456789eE+";
    const mutate = (s: string) => {
      const a = s.split('');
      const n = 1 + (rng.int() % 6);
      for (let k = 0; k < n; k++) {
        const at = rng.int() % (a.length + 1);
        const op = rng.int() % 3;
        const c = chars[rng.int() % chars.length];
        if (op === 0) a.splice(at, 1);
        else if (op === 1) a.splice(at, 0, c);
        else a[at] = c;
      }
      return a.join('');
    };
    for (let i = 0; i < 2000; i++) {
      for (const [read, base] of [
        [readNewick, tree],
        [readStockholm, Stockholm],
        [readAlignedFasta, fasta],
      ] as const) {
        const text = mutate(base);
        const out = (read as (t: string) => { issues: readonly unknown[] })(text);
        expect(Array.isArray(out.issues)).toBe(true);
      }
    }
    expect(readNewick(undefined as never).issues.length).toBe(1);
    expect(readStockholm(42 as never).rows).toEqual([]);
    expect(readAlignedFasta(null as never).rows).toEqual([]);
  });
});
