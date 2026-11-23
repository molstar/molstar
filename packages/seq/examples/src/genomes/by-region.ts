/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { HlaName } from './generate/assemblies.js';
import { type Files, fileSource, GenomesFiles, LettersA, LettersFormat } from './toys.js';

// The genome by region with its names aliased (examples.md, Genomes; EC-186 to EC-189, EC-200): toy A from the
// in-memory letters adapter (its sequences, lengths and letters by region, from the seed: no manifest, no file, EC-188)
// with its alias table (a `seq` source of systems and their aliases: `1`, `NC_900001.1` for chr1; chrM its accession,
// no `MT`), data written in other namings drawn on it (`1`, an accession, `chr3`, and `MT`, refused: never guessed
// onto chrM, EC-189), and the HLA contig, whose name carries `*` and `:`, from the small FASTA (the built reader,
// EC-200).

export const NamesSection = '/views/main/toy-a';
export const HlaSection = '/views/main/hla';

/** The marks' rows: each named as another file would name its sequence. */
export const AliasedMarks = {
  target: ['1', 'NC_900002.1', 'chr3', 'MT', HlaName],
  start: [10_000, 20_000, 30_000, 1_000, 500],
  end: [12_000, 23_000, 34_000, 2_000, 1_500],
  id: ['by number', 'by accession', 'by name', 'MT (refused)', 'on the HLA contig'],
};

export function namesAliased(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Toy A by region, its names aliased' });
  const a = b.data.assembly('toy-a', {
    name: 'Toy A',
    source: fileSource(b, 'letters', LettersFormat, LettersA),
    alias_source: fileSource(b, 'aliases', 'seq', files.aliases),
  } as never);
  const marks = b.data.ranges('marks', { target: a, rows: { ...AliasedMarks } } as never);
  const fasta = fileSource(b, 'small-fasta', 'fasta', files.fasta);
  const hla = b.data.sequence('hla', { alphabet: 'dna', source: fasta, name: HlaName } as never);
  b.views.view('main', (v) => {
    genome(v, 'toy-a', a, {}, (s) => {
      s.track('letters', { label: 'Letters', height: 18 }).layer('letters', { type: 'letters' });
      s.track('marks', { label: 'Named elsewhere' }).layer('marks', { type: 'ranges', data: marks } as never);
    });
    v.section('hla', { of: hla }, (s) => {
      s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
      s.track('letters', { label: 'HLA-A*01:01:01:01 (FASTA)', height: 18 }).layer('letters', { type: 'letters' });
    });
  });
  return b.build();
}
