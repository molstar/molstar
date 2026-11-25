/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { checkFai, FastaRuntime, faidx, openFasta, parseFai } from '../fasta.js';
import { DefaultOpenLimits, type OpenedSource } from '../protocol.js';
import { SourceOpenError } from '../seq-runtime.js';

// FASTA with its index: the index written from a file is the one samtools writes and reads the records back; a
// malformed or stale index is refused when the source opens, never read past.

const Fasta = '>chr1 the first\nACGTACGTAC\nGTACGTAC\n>chrM\nTTTTGGGGCC\nCCAA\n>empty\n';
const signal = new AbortController().signal;
const read = (fasta: string, index?: string) => openFasta(fasta, index);

async function letters(src: OpenedSource, system: string, start: number, end: number) {
  const a = await src.read({ system, start, end, level: 1, shapes: ['letters'] }, signal);
  if (!a.ok || !a.tables.letters) throw new Error('no letters');
  return a.tables.letters.letters;
}

function refusal(fasta: string, index?: string) {
  try {
    read(fasta, index);
  } catch (e) {
    return (e as SourceOpenError).issues[0];
  }
  return undefined;
}

describe('FASTA with its index', () => {
  it('builds the index samtools writes, and reads the records back by it (a round trip)', async () => {
    const fai = faidx(Fasta);
    expect(fai).toBe('chr1\t18\t16\t10\t11\nchrM\t14\t42\t10\t11\nempty\t0\t65\t0\t0\n');
    expect(parseFai(fai).records?.map((r) => r.name)).toEqual(['chr1', 'chrM', 'empty']);
    expect(checkFai(parseFai(fai).records!, Fasta)).toBeUndefined();
    const src = read(Fasta, fai);
    expect(src.coverage.systems).toEqual([
      { name: 'chr1', length: 18 },
      { name: 'chrM', length: 14 },
      { name: 'empty', length: 0 },
    ]);
    expect(await letters(src, 'chr1', 0, 18)).toBe('ACGTACGTACGTACGTAC');
    expect(await letters(src, 'chr1', 8, 12)).toBe('ACGT');
    expect(await letters(src, 'chrM', 10, 1024)).toBe('CCAA');
    // without an index, the same records
    expect(await letters(read(Fasta), 'chrM', 0, 14)).toBe('TTTTGGGGCCCCAA');
    // CRLF lines: two bytes a break
    const crlf = Fasta.replace(/\n/g, '\r\n');
    expect(await letters(read(crlf, faidx(crlf)), 'chrM', 8, 12)).toBe('CCCC');
  });

  it('refuses a malformed index: fields, integers, duplicates', () => {
    expect(parseFai('chr1\t18\t16\t10\n').issue).toMatch(/4 fields/);
    expect(parseFai('chr1\t18\t16\tten\t11\n').issue).toMatch(/integers/);
    expect(parseFai('a\t1\t3\t1\t2\na\t1\t3\t1\t2\n').issue).toMatch(/twice/);
    expect(refusal(Fasta, 'chr1\t18\t16\t10\n')?.code).toBe('seq.source.malformed');
  });

  it('refuses an index that does not fit the file: an offset, a header, a line, past the end', () => {
    expect(refusal(Fasta, 'chr1\t18\t17\t10\t11\n')?.message).toMatch(/offset/);
    expect(refusal(Fasta, 'chrX\t18\t16\t10\t11\n')?.message).toMatch(/header/);
    expect(refusal(Fasta, 'chr1\t18\t16\t9\t10\n')?.message).toMatch(/line/);
    expect(refusal(Fasta, 'chrM\t40\t42\t10\t11\n')?.message).toMatch(/past the file/);
    expect(refusal(Fasta, 'chr1\t18\t16\t10\t11\n')).toBeUndefined();
    expect(refusal('>a\nACÉ\n')?.message).toMatch(/ASCII/);
  });

  it('refuses a record longer than it is, a wrong CRLF line width, and overlapping records', () => {
    const two = '>a\nACGT\n>b\nACGT\n';
    expect(refusal(two, faidx(two))).toBeUndefined();
    expect(refusal(two, 'a\t6\t3\t4\t5\n')?.message).toMatch(/another record/);
    // the last base followed by neither a break nor the end of the file
    expect(refusal(two, 'a\t3\t3\t4\t5\n')?.message).toMatch(/last line/);
    const crlf = '>a\r\nACGT\r\nAC\r\n';
    expect(refusal(crlf, faidx(crlf))).toBeUndefined();
    expect(refusal(crlf, 'a\t6\t4\t4\t5\n')?.message).toMatch(/line/);
    expect(refusal(two, 'a\t4\t3\t4\t5\nb\t4\t11\t4\t5\nc\t2\t3\t4\t5\n')?.message).toMatch(/header|overlap/);
    const dup = '>a\nACGT\n>a x\nACGT\n';
    expect(refusal(dup, 'a\t4\t3\t4\t5\nb\t2\t3\t4\t5\n')).toBeDefined();
  });

  it('an empty line ends a record: a later sequence line in it is refused, as samtools does', () => {
    expect(() => faidx('>a\nACGT\n\nTTGG\n')).toThrow(/empty line/);
    expect(refusal('>a\nACGT\n\nTTGG\n')?.code).toBe('seq.source.malformed');
    expect(faidx('>a\nACGT\n\n>b\nTT\n')).toBe('a\t4\t3\t4\t5\nb\t2\t12\t2\t3\n');
  });

  it('without an index, duplicate names are the parse issue, never an empty source', () => {
    expect(refusal('>a\nAC\n>a\nGT\n')?.message).toMatch(/twice/);
  });

  it('as a runtime: the file and its index read through the assets port', async () => {
    const files: Record<string, string> = { 'g.fa': Fasta, 'g.fa.fai': faidx(Fasta), 'bad.fai': 'x\n' };
    const assets = {
      readText: (uri: string) => (uri in files ? Promise.resolve(files[uri]) : Promise.reject(new Error('no'))),
      readBytes: (uri: string) =>
        uri in files ? Promise.resolve(new TextEncoder().encode(files[uri])) : Promise.reject(new Error('no')),
    };
    const entry = { address: '/data/g', format: 'fasta', uri: 'g.fa' };
    const ok = await FastaRuntime.open({ ...entry, index: 'g.fa.fai' }, { assets, signal, limits: DefaultOpenLimits });
    expect(await letters(ok, 'chr1', 0, 4)).toBe('ACGT');
    await expect(
      FastaRuntime.open({ ...entry, index: 'bad.fai' }, { assets, signal, limits: DefaultOpenLimits }),
    ).rejects.toBeInstanceOf(SourceOpenError);
    await expect(
      FastaRuntime.open({ ...entry, uri: 'none.fa' }, { assets, signal, limits: DefaultOpenLimits }),
    ).rejects.toMatchObject({
      issues: [{ code: 'seq.source.unreadable', address: '/data/g' }],
    });
    const outside = await ok.read({ system: 'chrZ', start: 0, end: 4, level: 1, shapes: ['letters'] }, signal);
    expect(outside).toMatchObject({ ok: false, kind: 'not-served' });
  });
});
