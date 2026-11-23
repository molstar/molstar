/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { readVcfHeader } from '../vcf-header.js';
import { VcfCauses } from '../vcf-issues.js';
import { parseLine, VcfRowsBuilder, vcfSchema } from '../vcf-parse.js';
import { flagsLine, tablessLines, VcfHeaderText } from './vcf-fixtures.js';

// A line costs its own length (_spec/README.md "Data", VCF): the line parser's
// searches stay in the line and in an INFO field, counted here as the bytes and characters they scan (no clock), and
// `Number=A` and `R` values are split once a line. Their wall-clock bounds are in `vcf-time.test.ts` (heavy).

afterEach(() => vi.restoreAllMocks());

const header = readVcfHeader(VcfHeaderText);
if (!header.ok) throw new Error(header.message);

/** Each data line of `text` through the one line parser. */
function parseAll(text: Uint8Array) {
  const causes = new VcfCauses();
  const schema = vcfSchema(header.ok ? header.header : (undefined as never), undefined, undefined, 256, causes);
  const builder = new VcfRowsBuilder(schema, undefined);
  const options = { maxAlleles: 256, builder: () => builder, rowsLeft: () => 1_000_000 };
  const encoding = { replaced: false };
  for (let s = 0; s < text.length; ) {
    let e = text.indexOf(10, s);
    if (e < 0) e = text.length;
    parseLine(text, s, e, 'line', options, causes, encoding);
    s = e + 1;
  }
  return { builder, causes };
}

/** The elements a search scans: from its start to what it found, or to the end of what it searched. */
const scanned = (length: number, found: number, from = 0) => (found < 0 ? length : found + 1) - Math.max(0, from);

describe('a line costs its own length', () => {
  it('lines without their tabs: each tab search stays in its line (scanned bytes at most the file’s)', () => {
    const text = tablessLines(3000);
    const proto = Object.getPrototypeOf(Uint8Array.prototype) as { indexOf(v: number, from?: number): number };
    const indexOf = proto.indexOf;
    let total = 0;
    vi.spyOn(proto, 'indexOf').mockImplementation(function (this: Uint8Array, v: number, from?: number) {
      const r = indexOf.call(this, v, from);
      if (v === 9) total += scanned(this.length, r, from);
      return r;
    });
    const { builder, causes } = parseAll(text);
    expect(builder.rows).toBe(0);
    expect(causes.issues('/data/v')[0].message).toMatch(/^Fewer than eight columns: 3000/);
    expect(total).toBeLessThanOrEqual(text.length);
  });

  it('INFO of many flags: `=` looked for in each field only (scanned characters at most the INFO’s)', () => {
    const text = flagsLine(5000);
    const indexOf = String.prototype.indexOf;
    let total = 0;
    vi.spyOn(String.prototype, 'indexOf').mockImplementation(function (this: string, v: string, from?: number) {
      const r = indexOf.call(this, v, from);
      if (v === '=') total += scanned(this.length, r, from);
      return r;
    });
    const { builder } = parseAll(text);
    expect(builder.rows).toBe(1);
    expect(total).toBeLessThanOrEqual(text.length);
  });

  it('`Number=A` and `R` values split once a line, whatever its alleles', () => {
    const alts = ['C', 'G', 'T', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA'];
    const line = `chr1\t100\t.\tA\t${alts.join(',')}\t.\t.\tAF=${alts.map((_, k) => k / 10).join(',')};RC=${[0, ...alts.map((_, k) => k)].join(',')}`;
    const split = String.prototype.split;
    let commas = 0;
    vi.spyOn(String.prototype, 'split').mockImplementation(function (this: string, sep: unknown, limit?: number) {
      if (sep === ',') commas++;
      return (split as (this: string, s: unknown, l?: number) => string[]).call(this, sep, limit);
    });
    const { builder } = parseAll(new TextEncoder().encode(line));
    expect(builder.rows).toBe(10);
    // the ALTs, AF and RC: once each
    expect(commas).toBe(3);
  });
});
