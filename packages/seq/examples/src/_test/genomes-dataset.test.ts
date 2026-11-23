/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { existsSync, readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { ungzip } from '@molstar/core/util/zip/zip';
import { describe, expect, it, vi } from 'vitest';
import { HlaName, NGap } from '../genomes/generate/assemblies.js';
import { Inversion, TranslocationLength, toB } from '../genomes/generate/alignment.js';
import { GeneCount } from '../genomes/generate/genes.js';
import { EstimatedBytes, GenomesFiles, generateGenomes } from '../genomes/generate/index.js';
import { CpgCount, MissingShare, MethylationSamples } from '../genomes/generate/methylation.js';
import { DeletionExtent, deletionSnvPositions, SampleCount, SiteCount } from '../genomes/generate/variants.js';
import { ToyBSequences } from '../genomes/toys.js';
import { PackageFolder } from '../node-assets.js';

// The Genomes dataset's generator (examples.md, Genomes): its output byte-identical to the committed fixtures (by
// checksum), each file under its estimated size, the same bytes on a second run. `SEQ_RECORDS=update` writes the
// fixtures again (the record script's `genomes` case writes them with their `sources.json` rows).

const sha = (b: Uint8Array) => createHash('sha256').update(b).digest('hex');
const generated = generateGenomes();

describe('the Genomes dataset, generated', () => {
  it('is byte-identical to the committed fixtures (checksums), each file under its estimated size', async () => {
    const { files } = await generated;
    expect([...files.keys()].sort()).toEqual(Object.values(GenomesFiles).sort());
    const sizes: string[] = [];
    for (const path of Object.values(GenomesFiles)) {
      const bytes = files.get(path)!;
      const file = `${PackageFolder}${path}`;
      if (process.env.SEQ_RECORDS === 'update') writeFileSync(file, bytes);
      expect(existsSync(file), `no ${path}: run with SEQ_RECORDS=update`).toBe(true);
      expect(sha(new Uint8Array(readFileSync(file))), path).toBe(sha(bytes));
      sizes.push(`${path}\t${bytes.length}`);
    }
    if (process.env.SEQ_SIZES) console.log(sizes.join('\n'));
    for (const [k, path] of Object.entries(GenomesFiles))
      expect(files.get(path)!.length, path).toBeLessThan(EstimatedBytes[k as keyof typeof GenomesFiles]);
  });

  it('is the same on every run and every engine: no clock, no Math.random, no transcendental (an engine’s rounding)', async () => {
    // `exact.ts`: Math's approximations are not required to agree across engines, `+ - * /` round correctly everywhere
    const Approximate = ['sin', 'cos', 'tan', 'exp', 'expm1', 'log', 'log1p', 'log2', 'log10', 'pow', 'cbrt', 'atan2'];
    for (const f of Approximate)
      vi.spyOn(Math, f as 'sin').mockImplementation(() => {
        throw new Error(`Math.${f} in the generator`);
      });
    const again = await generateGenomes();
    vi.restoreAllMocks();
    const { files } = await generated;
    for (const [path, bytes] of files) expect(sha(again.files.get(path)!), path).toBe(sha(bytes));
    // nor `**`, Math.pow's operator (a spy cannot see it)
    const folder = `${PackageFolder}src/genomes/generate/`;
    for (const f of readdirSync(folder).filter((n) => n.endsWith('.ts')))
      expect(readFileSync(`${folder}${f}`, 'utf8').replace(/\/\*[\s\S]*?\*\/|\/\/.*$/gm, ''), f).not.toMatch(/\*\*/);
  });
});

/**
 * A region's lines of a bgzipped file (its BGZF members inflated whole, each line kept whose extent meets the region):
 * the index is read by the runtimes in the views' tests.
 */
async function region(file: string, _index: string, name: string, start: number, end: number): Promise<string[]> {
  const { files } = await generated;
  const bytes = files.get(file)! as Uint8Array<ArrayBuffer>;
  const text = new TextDecoder().decode(await ungzip(SyncRuntimeContext, bytes));
  const vcf = file.endsWith('.vcf.gz');
  return text.split('\n').filter((l) => {
    const c = l.split('\t');
    if (l.startsWith('#') || c[0] !== name) return false;
    const [b, e] = vcf
      ? [Number(c[1]) - 1, Number(c[1]) - 1 + c[3].length]
      : file.includes('bedgraph')
        ? [Number(c[1]), Number(c[2])]
        : [Number(c[3]) - 1, Number(c[4])];
    return b < end && e > start;
  });
}

describe('the Genomes dataset’s facts (anchors)', () => {
  it('toy A and toy B: A prefixed with chrM and the HLA contig, B unprefixed, the letters A’s lengths and N in chr1’s gap', async () => {
    const { truth } = await generated;
    expect(truth.a.sequences.map((s) => s.name)).toEqual(['chr1', 'chr2', 'chr3', 'chrM', HlaName]);
    expect(truth.alignment.b.sequences).toEqual(ToyBSequences);
    for (const s of truth.a.sequences) expect(truth.lettersA.get(s.name)!.length).toBe(s.length);
    expect(truth.lettersA.get('chr1')!.slice(NGap[0], NGap[1])).toBe('N'.repeat(NGap[1] - NGap[0]));
    // B's 2 holds A's chr2 and A's chr3 tail: about 260 kb (the small indels move it a little)
    const b2 = truth.alignment.b.sequences.find((s) => s.name === '2')!.length;
    expect(Math.abs(b2 - 200_000 - TranslocationLength)).toBeLessThan(2_000);
  });

  it('the alignment truth: the inversion reversed, the translocation onto B’s 2, chr1’s run of N unaligned (EC-206)', async () => {
    const { truth } = await generated;
    const blocks = truth.alignment.blocks;
    for (const k of blocks) expect(k.aEnd - k.aStart).toBe(k.bEnd - k.bStart);
    const inv = blocks.filter((k) => k.a === 'chr1' && k.aStart >= Inversion[0] && k.aEnd <= Inversion[1]);
    expect(inv.length).toBeGreaterThan(5);
    for (const k of inv) expect(k.strand).toBe(-1);
    // the first block of the inversion is the last on B (read reversed)
    expect(inv[0].bStart).toBeGreaterThan(inv[inv.length - 1].bStart);
    expect(blocks.some((k) => k.a === 'chr3' && k.b === '2')).toBe(true);
    expect(toB(blocks, 'chr1', (NGap[0] + NGap[1]) / 2)).toBeUndefined();
    // an aligned base is A's letter on B (1 % substituted): at least 95 % of a block's letters agree
    const k = blocks.find((b) => b.a === 'chr2')!;
    const a = truth.lettersA.get('chr2')!.slice(k.aStart, k.aEnd);
    const b = truth.alignment.letters.get(k.b)!.slice(k.bStart, k.bEnd);
    let same = 0;
    for (let i = 0; i < a.length; i++) if (a[i] === b[i]) same++;
    expect(same / a.length).toBeGreaterThan(0.95);
  });

  it('150 genes, 1 to 4 transcripts of 3 to 12 exons, both strands; the GFF3 and GTF read by region (EC-227)', async () => {
    const { truth } = await generated;
    expect(truth.genes.length).toBe(GeneCount);
    expect(new Set(truth.genes.map((g) => g.strand))).toEqual(new Set(['+', '-']));
    for (const g of truth.genes) {
      expect(g.transcripts.length).toBeGreaterThanOrEqual(1);
      expect(g.transcripts.length).toBeLessThanOrEqual(4);
      expect(g.transcripts[0].exons.length).toBeGreaterThanOrEqual(3);
      expect(g.transcripts[0].exons.length).toBeLessThanOrEqual(12);
    }
    const g = truth.genes[40];
    const gff = await region(GenomesFiles.gff3, GenomesFiles.gff3Index, g.chr, g.start, g.end);
    expect(gff.some((l) => l.includes(`ID=gene:${g.id};`))).toBe(true);
    const gtf = await region(GenomesFiles.gtf, GenomesFiles.gtfIndex, g.chr, g.start, g.end);
    expect(gtf.some((l) => l.includes(`gene_id "${g.id}"`))).toBe(true);
  });

  it('the VCF: 2,000 sites of 40 samples, chrM haploid, the 600-base deletion over 12 SNVs (EC-220, EC-222)', async () => {
    const { truth } = await generated;
    expect([...truth.variants.sites.values()].reduce((n, v) => n + v, 0)).toBe(SiteCount);
    const [d0, d1] = DeletionExtent;
    const lines = await region(GenomesFiles.vcf, GenomesFiles.vcfIndex, 'chr1', d0, d1);
    const del = lines.find((l) => l.split('\t')[2] === 'del600')!.split('\t');
    expect(del[3].length - del[4].length).toBe(600);
    const snvs = lines.filter((l) => l.split('\t')[3].length === 1).map((l) => Number(l.split('\t')[1]));
    expect(snvs).toEqual(deletionSnvPositions());
    expect(del.length).toBe(9 + SampleCount);
    const mt = await region(GenomesFiles.vcf, GenomesFiles.vcfIndex, 'chrM', 0, 16_569);
    expect(mt.length).toBeGreaterThan(0);
    for (const l of mt) for (const c of l.split('\t').slice(9)) expect(c).toMatch(/^[.\d]$/);
  });

  it('the bedGraph: 100 bp bins from 0 to 10^4, none over chr1’s run of N (missing, not 0)', async () => {
    const lines = await region(GenomesFiles.signal, GenomesFiles.signalIndex, 'chr1', NGap[0] - 1000, NGap[1] + 1000);
    const starts = lines.map((l) => Number(l.split('\t')[1]));
    expect(starts.filter((s) => s >= NGap[0] && s < NGap[1])).toEqual([]);
    expect(starts.length).toBe(20);
    const b = await region(GenomesFiles.signal, GenomesFiles.signalIndex, 'chr1', 0, 300_000);
    const values = b.map((l) => Number(l.split('\t')[3]));
    expect(Math.min(...values)).toBeGreaterThanOrEqual(0);
    expect(Math.max(...values)).toBeLessThanOrEqual(10_000);
    expect(Math.max(...values)).toBeGreaterThan(1_000);
  });

  it('methylation: 250 CpGs at A’s own CG letters, 20 samples, about 8 % missing', async () => {
    const { truth } = await generated;
    expect(truth.cpgs.positions.length).toBe(CpgCount);
    const chr1 = truth.lettersA.get('chr1')!;
    for (const p of truth.cpgs.positions) expect(chr1.slice(p, p + 2)).toBe('CG');
    const share = truth.cpgs.missing / (CpgCount * MethylationSamples.length);
    expect(Math.abs(share - MissingShare)).toBeLessThan(0.02);
  });

  it('the liftover: a duplicated position has two images, a deleted one none, a shifted one its offset (EC-212)', async () => {
    const { truth } = await generated;
    const at = (chr: string, p: number) => truth.liftover.find((t) => t.chr === chr && t.position === p)!.images;
    expect(at('chr1', 100_000)).toEqual([100_000]);
    expect(at('chr1', 205_000)).toEqual([205_000, 215_000]);
    expect(at('chr1', 250_000)).toEqual([260_000]);
    expect(at('chr2', 60_000)).toEqual([65_000]);
    expect(at('chr3', 32_000)).toEqual([]);
    expect(at('chr3', 50_000)).toEqual([46_000]);
  });
});
