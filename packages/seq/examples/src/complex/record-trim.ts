/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import { encodeSeq, inlineSeq, type SeqData } from '@molstar/seqx-core/data';

// What the record script keeps of PDB 1FIN's mmCIF and of the CDK2 entries its ligands come from (`scripts/record.ts
// complex`): the chains' letters and secondary structure, and what is derived from the coordinates (residue contacts
// within `ContactCutoff` between heavy atoms across the chains of a copy, and each ligand's contacting residues). Kept
// here, in the package's sources, so a unit test holds the rules; the script imports them.

/** Heavy atoms closer than this (Å) are in contact (examples' `contacts.cutoff`). */
export const ContactCutoff = 4;

/** The rows of a `loop_` (or a single-row category) of mmCIF text, each a record by item name. */
export function cifRows(text: string, category: string): Record<string, string>[] {
  const lines = text.split('\n');
  const prefix = `${category}.`;
  let i = lines.findIndex(
    (l, k) => l.startsWith(prefix) || (l.startsWith('loop_') && lines[k + 1]?.startsWith(prefix)),
  );
  if (i < 0) return [];
  const loop = lines[i].startsWith('loop_');
  if (loop) i++;
  const names: string[] = [];
  const tokens: string[] = [];
  const block = () => {
    let t = lines[i].slice(1);
    while (!lines[++i].startsWith(';')) t += `\n${lines[i]}`;
    tokens.push(t);
  };
  if (loop) {
    for (; lines[i]?.startsWith(prefix); i++) names.push(lines[i].trim().slice(prefix.length));
    for (; i < lines.length; i++) {
      const l = lines[i];
      if (l.startsWith('#') || l.startsWith('loop_') || l.startsWith('_')) break;
      if (l.startsWith(';')) block();
      else tokens.push(...tokenize(l));
    }
  } else {
    for (; i < lines.length && (lines[i].startsWith(prefix) || lines[i].startsWith(';')); i++) {
      if (lines[i].startsWith(';')) {
        block();
        continue;
      }
      const m = /^\S+?\.(\S+)\s*(.*)$/.exec(lines[i])!;
      names.push(m[1]);
      if (m[2].trim()) tokens.push(...tokenize(m[2]));
    }
  }
  const rows: Record<string, string>[] = [];
  for (let r = 0; r + names.length <= tokens.length; r += names.length) {
    const row: Record<string, string> = {};
    names.forEach((n, k) => (row[n] = tokens[r + k]));
    rows.push(row);
    if (!loop) break;
  }
  return rows;
}

function tokenize(line: string): string[] {
  const out: string[] = [];
  const re = /'((?:[^']|'(?=\S))*)'(?=\s|$)|"((?:[^"]|"(?=\S))*)"(?=\s|$)|(\S+)/g;
  for (const m of line.matchAll(re)) out.push(m[1] ?? m[2] ?? m[3]);
  return out;
}

interface Atom {
  readonly asym: string;
  readonly comp: string;
  readonly entity: string;
  /** `label_seq_id` (from 1), 0 for a non-polymer. */
  readonly seq: number;
  readonly auth: number;
  readonly name: string;
  readonly x: number;
  readonly y: number;
  readonly z: number;
}

/** The heavy atoms of the first model, the first alternate location of each. */
function atoms(cif: string): Atom[] {
  const out: Atom[] = [];
  for (const r of cifRows(cif, '_atom_site')) {
    if (r.type_symbol === 'H' || r.type_symbol === 'D' || r.label_comp_id === 'HOH') continue;
    if (r.pdbx_PDB_model_num !== '1' || (r.label_alt_id !== '.' && r.label_alt_id !== 'A')) continue;
    out.push({
      asym: r.label_asym_id,
      comp: r.label_comp_id,
      entity: r.label_entity_id,
      seq: r.label_seq_id === '.' ? 0 : Number(r.label_seq_id),
      auth: Number(r.auth_seq_id),
      name: r.label_atom_id,
      x: Number(r.Cartn_x),
      y: Number(r.Cartn_y),
      z: Number(r.Cartn_z),
    });
  }
  return out;
}

const close = (a: Atom, b: Atom) => (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2 < ContactCutoff ** 2;

/** Each atom of `a` in contact with one of `b`: the pairs of indices (brute force: a recording, not a test). */
function contactsOf(a: readonly Atom[], b: readonly Atom[]): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < a.length; i++) for (let j = 0; j < b.length; j++) if (close(a[i], b[j])) out.push([i, j]);
  return out;
}

/** A crystal chain: its label asym id, the entity's system and the copy of the complex it belongs to. */
export interface Chain {
  readonly id: string;
  readonly system: 'cdk2' | 'cyclin';
  readonly copy: number;
}

/** A secondary structure element: on its chain's system, 0-based half-open, `helix` or `strand`. */
export interface Element {
  readonly chain: string;
  readonly start: number;
  readonly end: number;
  readonly kind: 'helix' | 'strand';
}

/** A residue–residue contact across a copy's chains: CDK2's residue, cyclin A's (0-based, on their systems). */
export interface Contact {
  readonly copy: number;
  readonly cdk2: number;
  readonly cyclin: number;
}

/** A ligand's contact: the residue of CDK2 (0-based), the ligand's heavy atom (its index in `atoms`). */
export interface LigandContact {
  readonly residue: number;
  readonly atom: number;
}

/** A ligand of CDK2 in an entry: its component, chain, heavy atoms and contacting residues, the entry's resolution. */
export interface Ligand {
  readonly comp: string;
  readonly entry: string;
  /** The CDK2 chain it sits on (label asym id), its copy in 1FIN. */
  readonly chain: string;
  readonly resolution: number;
  readonly atoms: readonly string[];
  readonly contacts: readonly LigandContact[];
}

/** What is kept of 1FIN. */
export interface Trimmed1fin {
  readonly cdk2: string;
  readonly cyclin: string;
  /** The author number of the cyclin A construct's first residue (P20248's numbering). */
  readonly cyclinFirst: number;
  readonly chains: readonly Chain[];
  readonly secondary: readonly Element[];
  readonly contacts: readonly Contact[];
  /** Each copy's ATP. */
  readonly atp: readonly Ligand[];
  readonly resolution: number;
}

const resolutionOf = (cif: string) => Number(cifRows(cif, '_refine')[0]?.ls_d_res_high);

/** The entity whose sequence is CDK2's (the entity a chain of the ligand's contacts is on, its id). */
function lettersOf(cif: string): Map<string, { letters: string; chains: string[] }> {
  const out = new Map<string, { letters: string; chains: string[] }>();
  for (const r of cifRows(cif, '_entity_poly'))
    out.set(r.entity_id, {
      letters: r.pdbx_seq_one_letter_code_can.replaceAll(/\s/g, ''),
      chains: [],
    });
  for (const r of cifRows(cif, '_struct_asym')) out.get(r.entity_id)?.chains.push(r.id);
  return out;
}

/** A ligand's heavy atoms and the polymer residues of chain `chain` within the cutoff (by author number, 0-based). */
function ligandOn(all: readonly Atom[], asym: string, chain: string): { atoms: string[]; contacts: LigandContact[] } {
  const lig = all.filter((a) => a.asym === asym);
  const poly = all.filter((a) => a.asym === chain && a.seq > 0);
  const seen = new Set<string>();
  const contacts: LigandContact[] = [];
  for (const [i, j] of contactsOf(lig, poly)) {
    const residue = poly[j].auth - 1;
    const k = `${residue}|${i}`;
    if (seen.has(k)) continue;
    seen.add(k);
    contacts.push({ residue, atom: i });
  }
  contacts.sort((a, b) => a.residue - b.residue || a.atom - b.atom);
  return { atoms: lig.map((a) => a.name), contacts };
}

/** 1FIN's chains, their letters and secondary structure, the contacts across each copy and each copy's ATP's. */
export function trim1fin(cif: string): Trimmed1fin {
  const entities = lettersOf(cif);
  const [cdk2, cyclin] = [entities.get('1')!, entities.get('2')!];
  const all = atoms(cif);
  // a copy is a CDK2 chain and the cyclin A chain it touches most
  const chains: Chain[] = [];
  const contacts: Contact[] = [];
  const atp: Ligand[] = [];
  const nonpoly = cifRows(cif, '_pdbx_nonpoly_scheme');
  cdk2.chains.forEach((a, k) => {
    const copy = k + 1;
    const ca = all.filter((x) => x.asym === a);
    let best = '';
    let pairs: [number, number][] = [];
    for (const b of cyclin.chains) {
      const p = contactsOf(
        ca,
        all.filter((x) => x.asym === b),
      );
      if (p.length > pairs.length) [best, pairs] = [b, p];
    }
    const cb = all.filter((x) => x.asym === best);
    chains.push({ id: a, system: 'cdk2', copy }, { id: best, system: 'cyclin', copy });
    const seen = new Set<string>();
    for (const [i, j] of pairs) {
      const c = { copy, cdk2: ca[i].seq - 1, cyclin: cb[j].seq - 1 };
      const key = `${c.cdk2}|${c.cyclin}`;
      if (!seen.has(key)) {
        seen.add(key);
        contacts.push(c);
      }
    }
    // the copy's ATP: the one whose author chain is the CDK2 chain's
    const auth = all.find((x) => x.asym === a)!;
    const atpAsym = nonpoly.find((r) => r.mon_id === 'ATP' && r.pdb_strand_id === authChain(cif, a))?.asym_id;
    if (atpAsym && auth)
      atp.push({ comp: 'ATP', entry: '1FIN', chain: a, resolution: resolutionOf(cif), ...ligandOn(all, atpAsym, a) });
  });
  contacts.sort((x, y) => x.copy - y.copy || x.cdk2 - y.cdk2 || x.cyclin - y.cyclin);
  const system = (asym: string) => chains.find((c) => c.id === asym);
  const secondary: Element[] = [];
  const add = (kind: Element['kind'], asym: string, beg: string, end: string) => {
    if (!system(asym)) return;
    secondary.push({ chain: asym, start: Number(beg) - 1, end: Number(end), kind });
  };
  for (const r of cifRows(cif, '_struct_conf'))
    if (r.conf_type_id.startsWith('HELX')) add('helix', r.beg_label_asym_id, r.beg_label_seq_id, r.end_label_seq_id);
  for (const r of cifRows(cif, '_struct_sheet_range'))
    add('strand', r.beg_label_asym_id, r.beg_label_seq_id, r.end_label_seq_id);
  secondary.sort((x, y) => (x.chain < y.chain ? -1 : x.chain > y.chain ? 1 : x.start - y.start));
  const first = cifRows(cif, '_pdbx_poly_seq_scheme').find((r) => r.asym_id === cyclin.chains[0]);
  return {
    cdk2: cdk2.letters,
    cyclin: cyclin.letters,
    cyclinFirst: Number(first?.auth_seq_num),
    chains,
    secondary,
    contacts,
    atp,
    resolution: resolutionOf(cif),
  };
}

/** A label asym id's author chain. */
function authChain(cif: string, asym: string): string | undefined {
  return cifRows(cif, '_pdbx_poly_seq_scheme').find((r) => r.asym_id === asym)?.pdb_strand_id;
}

/**
 * A ligand of an entry of CDK2: the first instance of component `comp` and the CDK2 chain it touches most, its
 * contacts on CDK2's numbering (the entries number CDK2 as UniProt P24941 does: checked against `cdk2`'s letters).
 */
export function trimLigand(cif: string, entry: string, comp: string, cdk2: string): Ligand {
  const all = atoms(cif);
  const entities = lettersOf(cif);
  const kinase = [...entities.values()].find((e) => e.letters.includes(cdk2.slice(20, 60)));
  if (!kinase) throw new Error(`${entry}: no CDK2 entity`);
  const asyms = [...new Set(all.filter((a) => a.comp === comp).map((a) => a.asym))];
  if (!asyms.length) throw new Error(`${entry}: no ${comp}`);
  let best: { chain: string; atoms: string[]; contacts: LigandContact[] } | undefined;
  for (const chain of kinase.chains) {
    const r = ligandOn(all, asyms[0], chain);
    if (!best || r.contacts.length > best.contacts.length) best = { chain, ...r };
  }
  // the author numbering is UniProt's: each contacting residue's letter is CDK2's there
  const poly = all.filter((a) => a.asym === best!.chain && a.seq > 0);
  for (const c of best!.contacts) {
    const a = poly.find((x) => x.auth - 1 === c.residue)!;
    if (ThreeToOne[a.comp] !== cdk2[c.residue]) throw new Error(`${entry}: residue ${c.residue + 1} is not CDK2's`);
  }
  return {
    comp,
    entry,
    chain: best!.chain,
    resolution: resolutionOf(cif),
    atoms: best!.atoms,
    contacts: best!.contacts,
  };
}

const ThreeToOne: Record<string, string> = {
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

/** The systems' names in the `seq` files. */
export const Systems = { cdk2: 'cdk2', cyclin: 'cyclin-a' } as const;

/** Several systems' inline data as one `seq` file. */
function merged(parts: readonly SeqData[]): string {
  const data: SeqData = {
    coverage: {
      systems: parts.flatMap((p) => p.coverage.systems),
      shapes: parts.flatMap((p) => p.coverage.shapes),
      levels: [1],
    },
    pairs: parts.flatMap((p) => p.pairs),
  };
  return `${encodeSeq(data)}\n`;
}

/** 1FIN's `seq` file: each entity's letters and its chains' secondary structure (ranges with `chain` and `category`). */
export function moleculeSeq(t: Trimmed1fin): string {
  const part = (system: 'cdk2' | 'cyclin', letters: string) => {
    const rows = t.secondary
      .filter((e) => t.chains.find((c) => c.id === e.chain)?.system === system)
      .sort((a, b) => a.start - b.start || (a.chain < b.chain ? -1 : 1));
    return inlineSeq(
      { name: Systems[system], length: letters.length, alphabet: 'protein' },
      {
        letters: { form: 'letters', start: 0, letters },
        ranges: {
          form: 'list',
          start: Int32Array.from(rows, (e) => e.start),
          end: Int32Array.from(rows, (e) => e.end),
          id: categoryColumn(rows.map((e) => `${e.chain}${e.start + 1}`)),
          fields: [
            { name: 'chain', column: categoryColumn(rows.map((e) => e.chain)) },
            { name: 'category', column: categoryColumn(rows.map((e) => e.kind)) },
          ],
        },
      },
    );
  };
  return merged([part('cdk2', t.cdk2), part('cyclin', t.cyclin)]);
}

/**
 * The contacts' `seq` file, on CDK2's system: a range a contacting residue, the partner (`with`: cyclin A or a ligand's
 * component), its position on its own system (`partner`: a residue of cyclin A, a ligand's heavy atom), the copy
 * (`copy`) and the entry (`entry`).
 */
export function contactsSeq(t: Trimmed1fin, ligands: readonly Ligand[]): string {
  const rows: { residue: number; partner: number; with: string; copy: string; entry: string }[] = [];
  for (const c of t.contacts)
    rows.push({ residue: c.cdk2, partner: c.cyclin, with: Systems.cyclin, copy: String(c.copy), entry: '1FIN' });
  for (const l of [...t.atp, ...ligands]) {
    const copy = String(t.chains.find((c) => c.id === l.chain && l.entry === '1FIN')?.copy ?? 1);
    for (const c of l.contacts) rows.push({ residue: c.residue, partner: c.atom, with: l.comp, copy, entry: l.entry });
  }
  // a ranges list is sorted by start; each row keeps its place among those of one start
  rows.sort((a, b) => a.residue - b.residue);
  return `${encodeSeq(
    inlineSeq(
      { name: Systems.cdk2, length: t.cdk2.length, alphabet: 'protein' },
      {
        ranges: {
          form: 'list',
          start: Int32Array.from(rows, (r) => r.residue),
          end: Int32Array.from(rows, (r) => r.residue + 1),
          fields: [
            { name: 'partner', column: Int32Array.from(rows, (r) => r.partner) },
            { name: 'with', column: categoryColumn(rows.map((r) => r.with)) },
            { name: 'copy', column: categoryColumn(rows.map((r) => r.copy)) },
            { name: 'entry', column: categoryColumn(rows.map((r) => r.entry)) },
          ],
        },
      },
    ),
  )}\n`;
}
