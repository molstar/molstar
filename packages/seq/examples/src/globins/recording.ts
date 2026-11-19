/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AlignmentOp, globalAlign } from '@molstar/core/algorithms/sequence/needleman-wunsch';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { SeqData } from '@molstar/seqx-core/data';
import { cifLoop } from '../files/cif.js';

// What the record script does with the Globins entries (`scripts/record.ts globins`): the Pfam seed cut to the rows
// nearest human hemoglobin, 4HHB's chains, helices and heme contacts read from its mmCIF. Kept in the package's sources
// (no relative import: the script runs it under Node as it is) so a unit test holds the rules.

/** The seed's attribution (CC0: given as a courtesy). */
export const SeedAttribution = 'Pfam (EMBL-EBI), PF00042 (Globin) seed alignment, Pfam 39.0';

/** UniProt's attribution, with the changes made and the licence's link. */
export const UniprotAttribution =
  'UniProt Consortium, UniProtKB (https://www.uniprot.org), CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/); changed: the organism and length kept, the class derived from the taxonomic lineage, each record keyed by its Pfam seed row';

/** How many of the seed's rows are kept (examples.md, Globins: a seed over 60 rows is cut to 60). */
export const SeedRows = 60;

const Codes = 'ACDEFGHIKLMNPQRSTVWYX';
const code = (s: string) => Uint8Array.from(s.toUpperCase(), (c) => Math.max(0, Codes.indexOf(c)));

/**
 * The identity of `a` to `b`: the pairs of one letter in a global alignment (core’s Needleman–Wunsch, +4 a match, −1
 * a mismatch, gaps 8 and 1, the ends free), over the shorter's length.
 */
export function identity(a: string, b: string): number {
  const { ops } = globalAlign(code(a), code(b), {
    match: 4,
    mismatch: -1,
    gapOpen: 8,
    gapExtend: 1,
    leadingGaps: 'free',
    trailingGaps: 'free',
    gapTies: 'open',
  });
  let i = 0;
  let j = 0;
  let same = 0;
  for (const op of ops) {
    if (op === AlignmentOp.Pair) {
      if (a[i].toUpperCase() === b[j].toUpperCase()) same++;
      i++;
      j++;
    } else if (op === AlignmentOp.OnlyA) i++;
    else j++;
  }
  return same / Math.min(a.length, b.length);
}

/** A Stockholm row's residues: its letters less the gaps. */
export const residuesOf = (letters: string) => letters.replace(/[-.]/g, '');

/** The seed's rows (name and letters, one line each as Pfam writes them). */
export function seedRows(text: string): { name: string; letters: string }[] {
  const out: { name: string; letters: string }[] = [];
  for (const line of text.split('\n')) {
    if (!line || line.startsWith('#') || line.startsWith('//')) continue;
    const [name, letters] = line.trim().split(/\s+/);
    if (name && letters) out.push({ name, letters });
  }
  return out;
}

/**
 * The seed cut to `SeedRows`: the rows nearest human α or β by identity kept (each row's best of the two), the rest
 * removed with their `#=GS` and `#=GR` lines, `#=GF SQ` (the row count) said again; every other line as published.
 * Returns the text and the rows kept and removed. A seed of `SeedRows` rows or fewer comes back as it is.
 */
export function cutSeed(text: string, human: readonly string[]): { text: string; removed: string[]; kept: string[] } {
  const rows = seedRows(text);
  const near = rows.map((r) => Math.max(...human.map((h) => identity(residuesOf(r.letters), h))));
  const order = rows.map((_, k) => k).sort((x, y) => near[y] - near[x] || x - y);
  const keep = new Set(order.slice(0, SeedRows).map((k) => rows[k].name));
  const removed = rows.filter((r) => !keep.has(r.name)).map((r) => r.name);
  const kept = rows.filter((r) => keep.has(r.name)).map((r) => r.name);
  if (!removed.length) return { text, removed, kept };
  const lines = text
    .split('\n')
    .filter((line) => {
      const t = line.trim().split(/\s+/);
      if (line.startsWith('#=GS') || line.startsWith('#=GR')) return keep.has(t[1]);
      if (!line || line.startsWith('#') || line.startsWith('//')) return true;
      return keep.has(t[0]);
    })
    .map((line) => (line.startsWith('#=GF SQ') ? line.replace(/\d+\s*$/, String(kept.length)) : line));
  return { text: lines.join('\n'), removed, kept };
}

/** The entry and its version, from `#=GF AC` (`PF00042.29`). */
export function seedVersion(text: string): string | undefined {
  return /^#=GF AC\s+(\S+)/m.exec(text)?.[1];
}

/** The seed's `#=GS <row> AC <accession>.<version>` lines: each row's UniProt accession, unversioned. */
export function seedAccessions(text: string): Map<string, string> {
  const out = new Map<string, string>();
  for (const line of text.split('\n')) {
    const t = line.trim().split(/\s+/);
    if (t[0] === '#=GS' && t[2] === 'AC') out.set(t[1], t[3].replace(/\.\d+$/, ''));
  }
  return out;
}

/**
 * The rows' columns file from UniProt's TSV (`accession, id, organism_name, length, lineage`): a record a row of the
 * seed, in its order: `id` the row's label (derived: the row whose `AC` is the entry), `OS` the organism as UniProt
 * names it, `class` the lineage's name of rank class (derived; empty where it has none), `uniprot_length` the entry's
 * length (the whole protein, the row being its globin domain).
 */
export function columnsFromUniprot(uniprotTsv: string, seedText: string): string {
  const [head, ...lines] = uniprotTsv.trim().split('\n');
  const col = head.split('\t');
  const at = (name: string) => col.indexOf(name);
  const byAcc = new Map(lines.map((l) => l.split('\t')).map((c) => [c[at('Entry')], c] as const));
  const out = ['id\tOS\tclass\tuniprot_length'];
  for (const [row, acc] of seedAccessions(seedText)) {
    const c = byAcc.get(acc);
    if (!c) throw new Error(`UniProt has no ${acc} (${row})`);
    const cls = c[at('Taxonomic lineage')]
      .split(', ')
      .find((x) => x.endsWith(' (class)'))
      ?.slice(0, -' (class)'.length);
    out.push(`${row}\t${c[at('Organism')]}\t${cls ?? ''}\t${c[at('Length')]}`);
  }
  return `${out.join('\n')}\n`;
}

export { cifLoop };

const OneLetter: Readonly<Record<string, string>> = {
  ALA: 'A',
  ARG: 'R',
  ASN: 'N',
  ASP: 'D',
  CYS: 'C',
  GLN: 'Q',
  GLU: 'E',
  GLY: 'G',
  HIS: 'H',
  ILE: 'I',
  LEU: 'L',
  LYS: 'K',
  MET: 'M',
  PHE: 'F',
  PRO: 'P',
  SER: 'S',
  THR: 'T',
  TRP: 'W',
  TYR: 'Y',
  VAL: 'V',
};

/** A chain of 4HHB: its letters (`pdbx_poly_seq_scheme`, by `label_seq_id`) and its `auth` number of the first. */
export interface CifChain {
  readonly asym: string;
  readonly letters: string;
  readonly authFirst: number;
}

/** The polymer chains, by `asym_id`. */
export function cifChains(text: string): CifChain[] {
  const s = cifLoop(text, 'pdbx_poly_seq_scheme');
  const by = new Map<string, { letters: string; authFirst: number }>();
  s.asym_id.forEach((asym, k) => {
    const c = by.get(asym) ?? { letters: '', authFirst: Number(s.auth_seq_num[k]) };
    c.letters += OneLetter[s.mon_id[k]] ?? 'X';
    by.set(asym, c);
  });
  return [...by].map(([asym, c]) => ({ asym, ...c }));
}

/** A helix of `struct_conf`: its chain, its name (`pdbx_PDB_helix_id`, `AA`…), its `label_seq_id` range (1-based). */
export interface CifHelix {
  readonly asym: string;
  readonly name: string;
  readonly start: number;
  readonly end: number;
}

export function cifHelices(text: string): CifHelix[] {
  const s = cifLoop(text, 'struct_conf');
  return s.id.map((_, k) => ({
    asym: s.beg_label_asym_id[k],
    name: s.pdbx_PDB_helix_id[k],
    start: Number(s.beg_label_seq_id[k]),
    end: Number(s.end_label_seq_id[k]),
  }));
}

/** The heavy-atom distance a residue contacts a heme within, Å (examples.md, Globins). */
export const ContactDistance = 4;

/**
 * Each heme's contacting residues (derived): the residues of the heme's own chain (`auth_asym_id`) with a heavy atom
 * within `ContactDistance` of one of the heme's, as `label_seq_id`s ascending; by the heme's `auth_asym_id`.
 */
export function hemeContacts(text: string): Map<string, number[]> {
  const a = cifLoop(text, 'atom_site');
  const n = a.id.length;
  const out = new Map<string, number[]>();
  const hemes = new Map<string, number[]>();
  for (let k = 0; k < n; k++)
    if (a.label_comp_id[k] === 'HEM' && a.type_symbol[k] !== 'H') {
      const c = a.auth_asym_id[k];
      hemes.set(c, [...(hemes.get(c) ?? []), k]);
    }
  const xyz = (k: number) => [Number(a.Cartn_x[k]), Number(a.Cartn_y[k]), Number(a.Cartn_z[k])];
  const d2 = ContactDistance * ContactDistance;
  for (const [chain, atoms] of hemes) {
    const ps = atoms.map(xyz);
    const seen = new Set<number>();
    for (let k = 0; k < n; k++) {
      if (a.group_PDB[k] !== 'ATOM' || a.auth_asym_id[k] !== chain || a.type_symbol[k] === 'H') continue;
      const [x, y, z] = xyz(k);
      if (ps.some(([u, v, w]) => (x - u) ** 2 + (y - v) ** 2 + (z - w) ** 2 <= d2)) seen.add(Number(a.label_seq_id[k]));
    }
    out.set(
      chain,
      [...seen].sort((x, y) => x - y),
    );
  }
  return out;
}

/**
 * 4HHB as `seq` data: the α and β chains' letters (chains A and B; C and D are their copies) with their helices as
 * ranges (`id` the helix's letter, A to H; `label` `helix A`…), and the heme as a system of one position.
 */
export function hemoglobinSeq(text: string): SeqData {
  const chains = cifChains(text);
  const helices = cifHelices(text);
  const systems = [
    { name: 'hba', asym: 'A' },
    { name: 'hbb', asym: 'B' },
  ].map((s) => ({ ...s, letters: chains.find((c) => c.asym === s.asym)!.letters }));
  const pairs: SeqData['pairs'][number][] = systems.map((s) => {
    const own = helices.filter((h) => h.asym === s.asym);
    return {
      request: { system: s.name, start: 0, end: s.letters.length, level: 1, shapes: ['letters', 'ranges'] },
      answer: {
        ok: true,
        tables: {
          letters: { form: 'letters', start: 0, letters: s.letters },
          ranges: {
            form: 'list',
            start: Int32Array.from(own, (h) => h.start - 1),
            end: Int32Array.from(own, (h) => h.end),
            id: categoryColumn(own.map((h) => h.name.slice(1))),
            fields: [{ name: 'label', column: categoryColumn(own.map((h) => `helix ${h.name.slice(1)}`)) }],
          },
        },
      },
    };
  });
  pairs.push({
    request: { system: 'hem', start: 0, end: 1, level: 1, shapes: ['letters'] },
    answer: { ok: true, tables: { letters: { form: 'letters', start: 0, letters: 'X' } } },
  });
  return {
    coverage: {
      systems: [
        ...systems.map((s) => ({ name: s.name, length: s.letters.length, alphabet: 'protein' })),
        { name: 'hem', length: 1, alphabet: 'protein' },
      ],
      shapes: [
        ...systems.map((s) => ({ system: s.name, shape: 'letters' as const })),
        ...systems.map((s) => ({
          system: s.name,
          shape: 'ranges' as const,
          fields: [{ name: 'label', type: 'category' as const }],
        })),
        { system: 'hem', shape: 'letters' as const },
      ],
      levels: [1],
    },
    pairs,
  };
}
