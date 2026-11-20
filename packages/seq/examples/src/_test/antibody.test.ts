/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { decodeSeq } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { cdrAddress, cdrRows, Data, Files, Lengths, SystemNames, theDesign, wildType } from '../antibody/data.js';
import { DesignEdits, design } from '../antibody/design.js';
import { type Chain, ImgtCdrs, imgtCdrLabels, imgtOrder, kabatFromImgt } from '../antibody/numbering.js';
import { fabFasta, numberingSeq, segmentLabels } from '../antibody/record-trim.js';
import { PackageFolder } from '../node-assets.js';

// The Antibody's fixtures and rules against the data (examples.md, Antibody): the files are what the derived data and
// the generators say; the anchors of both schemes; IMGT's rule gives SAbDab's published labels on the wild type; the
// Kabat labels derived from them are trastuzumab's known ones; the design's numberings by the same rules (EC-7).

const read = (f: string) => readFileSync(`${PackageFolder}${f}`, 'utf8');
const Chains: readonly Chain[] = ['light', 'heavy'];

/** A numbering file's labels of a chain, a label a position (`null` where none). */
function fileLabels(file: string, chain: Chain): (string | null)[] {
  const { data } = decodeSeq(read(file));
  const pair = data!.pairs.find((p) => p.request.system === SystemNames[chain])!;
  const r = pair.answer.ok ? pair.answer.tables.ranges : undefined;
  if (r?.form !== 'list') throw new Error('no ranges');
  const out: (string | null)[] = Array.from({ length: pair.request.end }, () => null);
  const col = r.fields![0].column as unknown as { dict: string[]; codes: ArrayLike<number> };
  // a row a position, its label missing outside the variable domain
  expect(Array.from(r.start)).toEqual(out.map((_, i) => i));
  for (let k = 0; k < r.start.length; k++) out[r.start[k]] = col.codes[k] < 0 ? null : col.dict[col.codes[k]];
  return out;
}

describe('the fixtures', () => {
  it('the FASTA holds both chains as deposited (214 and 220 residues), the derived data’s letters', () => {
    expect(Data.light.letters).toHaveLength(Lengths.light);
    expect(Data.heavy.letters).toHaveLength(Lengths.heavy);
    const chains = {
      light: { letters: Data.light.letters, description: '' },
      heavy: { letters: Data.heavy.letters, description: '' },
    };
    const text = read(Files.fasta);
    expect(text.replaceAll(/^>(\S+).*$/gm, '>$1')).toBe(fabFasta(chains).replaceAll(/^>(\S+).*$/gm, '>$1'));
  });

  it('the IMGT file is the derived data’s labels; the Kabat file is what its generator writes from them', () => {
    for (const c of Chains) expect(fileLabels(Files.imgt, c)).toEqual(Data[c].imgt);
    const kabat = Chains.map((c) => ({ name: SystemNames[c], labels: kabatFromImgt(c, Data[c].imgt) }));
    expect(read(Files.kabat)).toBe(numberingSeq(kabat));
  });

  it('SAbDab’s rows say where their licence is stated, quoted, and hash the entry that chose the segments', () => {
    const rows = JSON.parse(read('fixtures/antibody/sources.json')) as {
      file: string;
      licence: string[];
      licenceEvidence?: { url: string; quoted: string };
      selectedBy?: { url: string; sha256: string };
    }[];
    const sabdab = rows.filter((r) => r.file.endsWith('.seq') && r.file.includes('1N8Z'));
    expect(sabdab).toHaveLength(2);
    for (const r of sabdab) {
      expect(r.licence).toEqual(['CC-BY-4.0']);
      expect(r.licenceEvidence?.url).toBe('https://sabdab.opig.stats.ox.ac.uk/api/openapi.json');
      expect(JSON.parse(r.licenceEvidence!.quoted).name).toBe(
        'Creative Commons Attribution 4.0 International (CC BY 4.0)',
      );
      expect(r.selectedBy?.url).toBe('https://sabdab.opig.stats.ox.ac.uk/api/pdb/pdb_00001n8z');
      expect(r.selectedBy?.sha256).toMatch(/^[0-9a-f]{64}$/);
    }
  });

  it('a segment residue that is not its chain’s is refused', () => {
    expect(() =>
      segmentLabels(
        {
          sequence: 'AC',
          numbering: [
            [[1, ' '], 'A'],
            [[2, ' '], 'W'],
          ],
        },
        'GACG',
      ),
    ).toThrow();
    expect(
      segmentLabels(
        {
          sequence: 'AC',
          numbering: [
            [[1, ' '], 'A'],
            [[2, ' '], '-'],
            [[3, 'A'], 'C'],
          ],
        },
        'GACG',
      ),
    ).toEqual([null, '1', '3A', null]);
  });
});

describe('the anchors', () => {
  const wt = wildType();
  const at = (c: Chain, scheme: 'kabat' | 'imgt', label: string) => wt[c].letters[wt[c][scheme].indexOf(label)];

  it('Kabat L23, L88, H22 and H92 are cysteines', () => {
    for (const l of ['23', '88']) expect(at('light', 'kabat', l)).toBe('C');
    for (const l of ['22', '92']) expect(at('heavy', 'kabat', l)).toBe('C');
  });

  it('IMGT 23 and 104 are cysteines, 41 a tryptophan, 118 a phenylalanine or tryptophan, in both chains', () => {
    for (const c of Chains) {
      expect(at(c, 'imgt', '23')).toBe('C');
      expect(at(c, 'imgt', '104')).toBe('C');
      expect(at(c, 'imgt', '41')).toBe('W');
      expect(['F', 'W']).toContain(at(c, 'imgt', '118'));
    }
  });

  it('the heavy chain’s Kabat insertions are 52A, 82A–C and 100A–C; the light chain has none', () => {
    const ins = (c: Chain) => wt[c].kabat.filter((l) => l !== null && /[A-Z]$/.test(l));
    expect(ins('heavy')).toEqual(['52A', '82A', '82B', '82C', '100A', '100B', '100C']);
    expect(ins('light')).toEqual([]);
  });

  it('trastuzumab’s CDR-H3 is 13 residues by IMGT (105–117, no insertion), 11 by Kabat (95–102)', () => {
    const h3 = cdrRows('heavy', 'imgt', wt.heavy.imgt);
    expect(wt.heavy.letters.slice(h3.start[2], h3.end[2])).toBe('SRWGGDGFYAMDY');
    const k3 = cdrRows('heavy', 'kabat', wt.heavy.kabat);
    expect(wt.heavy.letters.slice(k3.start[2], k3.end[2])).toBe('WGGDGFYAMDY');
    expect(wt.light.letters.slice(...spanOf(cdrRows('light', 'kabat', wt.light.kabat), 0))).toBe('RASQDVNTAVA');
  });

  it('the variable domains end where SAbDab’s numbering does: light at 107, heavy at 113 (Kabat)', () => {
    expect(wt.light.kabat.filter((l) => l !== null).at(-1)).toBe('107');
    expect(wt.heavy.kabat.filter((l) => l !== null).at(-1)).toBe('113');
  });
});

const spanOf = (rows: { start: number[]; end: number[] }, k: number) => [rows.start[k], rows.end[k]] as const;

describe('the Kabat counting refuses frameworks it does not assume', () => {
  /** The wild type's IMGT labels with one framework label moved to another position. */
  const moved = (c: Chain, from: string, to: string) => Data[c].imgt.map((l) => (l === from ? to : l));

  it('a heavy chain with IMGT 73 occupied (74 empty), a light chain with 81 occupied (a lambda chain’s)', () => {
    expect(() => kabatFromImgt('heavy', moved('heavy', '74', '73'))).toThrow(/IMGT 73 is occupied/);
    expect(() => kabatFromImgt('light', moved('light', '80', '81'))).toThrow(/IMGT 80 is empty/);
  });

  it('an FR2 residue missing, and a framework insertion, are refused', () => {
    expect(() => kabatFromImgt('heavy', moved('heavy', '45', null as unknown as string))).toThrow(/IMGT 45 is empty/);
    expect(() => kabatFromImgt('heavy', moved('heavy', '45', '44A'))).toThrow(/framework insertion/);
  });
});

describe('IMGT’s rule reproduces the published numbering', () => {
  it('each CDR of the wild type, numbered by its length, is SAbDab’s labels', () => {
    for (const c of Chains)
      for (const cdr of [1, 2, 3] as const) {
        const [first, last] = ImgtCdrs[cdr];
        const published = Data[c].imgt.filter((l) => l !== null && Number(l) >= first && Number(l) <= last);
        expect(imgtCdrLabels(cdr, published.length), `${c} CDR${cdr}`).toEqual(published);
      }
  });

  it('a CDR3 longer than 13 gains 111.1, 112.1, 111.2, …, 112.n before 112', () => {
    expect(imgtCdrLabels(3, 15).slice(6, 10)).toEqual(['111', '111.1', '112.1', '112']);
    const l17 = imgtCdrLabels(3, 17);
    expect(l17.slice(6, 12)).toEqual(['111', '111.1', '111.2', '112.2', '112.1', '112']);
    expect([...l17].sort((a, b) => imgtOrder(a) - imgtOrder(b))).toEqual(l17);
  });
});

describe('the design', () => {
  const d = theDesign();

  it('its edits: CDR-H2 and CDR-H3 substituted, two residues inserted in CDR-H3, one deleted from CDR-L1', () => {
    expect(d.heavy.letters).toHaveLength(Lengths.heavy + 2);
    expect(d.light.letters).toHaveLength(Lengths.light - 1);
    const k3 = cdrRows('heavy', 'kabat', d.heavy.kabat);
    expect(d.heavy.letters.slice(k3.start[2], k3.end[2])).toBe('WGGEGFGYYAMDY');
    const k2 = cdrRows('heavy', 'kabat', d.heavy.kabat);
    expect(d.heavy.letters.slice(k2.start[1], k2.end[1])).toBe('RIYPTQGHTRYADSVKG');
    expect(d.light.letters.slice(...spanOf(cdrRows('light', 'kabat', d.light.kabat), 0))).toBe('RASQDVTAVA');
    expect(DesignEdits).toHaveLength(5);
  });

  it('IMGT’s CDR-H3 has 111.1 and 112.1 out of index order (EC-7); Kabat’s gains 100D and 100E', () => {
    const imgt = d.heavy.imgt;
    const i = imgt.indexOf('111.1');
    expect(imgt.slice(i - 1, i + 3)).toEqual(['111', '111.1', '112.1', '112']);
    expect(d.heavy.kabat.filter((l) => l?.startsWith('100'))).toEqual(['100', '100A', '100B', '100C', '100D', '100E']);
  });

  // this generator's rule for a shorter CDR-L1 (Kabat: from its insertion point down; IMGT: the loop's top, the
  // C-terminal side first), not checked against ANARCI: no oracle is recorded
  it('CDR-L1 one shorter, by this generator’s rule (not checked against ANARCI): Kabat loses 27, IMGT loses 36', () => {
    expect(d.light.kabat).not.toContain('27');
    expect(d.light.kabat).toContain('28');
    expect(d.light.imgt.filter((l) => l !== null && Number(l) >= 27 && Number(l) <= 38)).toEqual([
      '27',
      '28',
      '29',
      '37',
      '38',
    ]);
  });

  it('the alignment to the wild type: as wide in both rows, each row its letters, the indels gaps', () => {
    for (const c of Chains) {
      const { design: a, wildType: b } = d[c].aligned;
      expect(a).toHaveLength(b.length);
      expect(a.replaceAll('-', '')).toBe(d[c].letters);
      expect(b.replaceAll('-', '')).toBe(Data[c].letters);
    }
    expect([...d.heavy.aligned.wildType].filter((x) => x === '-')).toHaveLength(2);
    expect([...d.light.aligned.design].filter((x) => x === '-')).toHaveLength(1);
  });

  it('the columns pair the residues the design kept, whose labels the schemes change around an indel', () => {
    const wt = wildType().heavy;
    const { design: a, wildType: b } = d.heavy.aligned;
    // a column's residue in each row (-1 at a gap)
    const at = (row: string, col: number) => (row[col] === '-' ? -1 : row.slice(0, col).replaceAll('-', '').length);
    const col = (pred: (wi: number) => boolean) => [...b].findIndex((_, c) => at(b, c) >= 0 && pred(at(b, c)));
    // the wild type's F (Kabat 100, IMGT 112) and the design's F in one column: IMGT 111.1 in the design
    const f = col((wi) => wt.kabat[wi] === '100');
    expect([a[f], b[f]]).toEqual(['F', 'F']);
    expect([d.heavy.imgt[at(a, f)], wt.imgt[at(b, f)]]).toEqual(['111.1', '112']);
    // the two inserted residues right after it, gaps in the wild type
    expect([b[f + 1], b[f + 2]]).toEqual(['-', '-']);
    // the wild type's Y100A is the design's 100C
    const y = col((wi) => wt.kabat[wi] === '100A');
    expect([wt.kabat[at(b, y)], d.heavy.kabat[at(a, y)]]).toEqual(['100A', '100C']);
  });

  it('an edit outside the CDRs that changes the length is refused; with no edit the design is the wild type', () => {
    expect(() => design(Data, [{ chain: 'heavy', at: '10', delete: true }])).toThrow(/not in a CDR/);
    const same = design(Data, []);
    for (const c of Chains) {
      expect(same[c].letters).toBe(Data[c].letters);
      expect(same[c].imgt).toEqual(Data[c].imgt);
    }
  });

  it('a CDR’s address spans its labels as the chain has them', () => {
    expect(cdrAddress('heavy', 3, 'kabat', d.heavy.kabat)).toBe(
      '/data/heavy/(range :from 95 :to 102 :numbering kabat)',
    );
    expect(cdrAddress('heavy', 1, 'kabat', wildType().heavy.kabat)).toBe(
      '/data/heavy/(range :from 31 :to 35 :numbering kabat)',
    );
  });
});
