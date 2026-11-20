/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AlignmentOp, globalAlign } from '@molstar/core/algorithms/sequence/needleman-wunsch';
import { categoryColumn } from '@molstar/core/repr/tree';
import { encodeSeq, inlineSeq } from '@molstar/seqx-core/data';
import { type Cif, cifRows } from '../files/cif.js';

// The record rules of p53's structure (examples.md, p53; structure-syncing.md, Datasets), from PDB 2OCJ's mmCIF as
// committed (`fixtures/structures/2OCJ.cif.gz`), by `scripts/record.ts p53-structure`: the entity's letters with chain
// A's secondary structure as Seq*'s `seq` file, and the fake structure viewer's own chain table (every chain instance,
// its residues with their label and author numbers, observed or not, and the model it holds beside them: entities,
// compounds, bonds, a molecule's atoms, a residue's B-factor), read from the same file and held apart from
// Seq*'s nodes. Here, in the package's sources, so a unit test holds the rules; the script imports them.

/** The entity's system in 2OCJ's `seq` file. */
export const StructureSystem = '2OCJ-1';

/** A residue of a chain instance in the fake viewer's table. */
export interface ChainResidue {
  /** `label_seq_id` (1-based), null for a non-polymer's residue. */
  readonly label: number | null;
  readonly auth: string;
  /** `pdbx_PDB_ins_code`, empty for none. */
  readonly ins: string;
  readonly comp: string;
  readonly observed: boolean;
  /** The mean B-factor of the residue's heavy atoms (the model's per-residue value); absent: no atoms. */
  readonly b?: number;
}

/** A chain instance: `label_asym_id`, `auth_asym_id`, its entity and operators, its residues in order. */
export interface ChainInstance {
  readonly label_asym_id: string;
  readonly auth_asym_id: string;
  readonly entity: string;
  readonly operators: readonly string[];
  readonly residues: readonly ChainResidue[];
}

/** An entity of the host's model (`entity`, `entity_poly`): what a host that holds a structure knows of each. */
export interface ModelEntity {
  readonly id: string;
  readonly type: string;
  readonly name?: string;
  readonly polymerType?: string;
  /** `entity_poly.pdbx_seq_one_letter_code_can`, whitespace removed. */
  readonly letters?: string;
}

/** A component the model describes (`chem_comp`): what a residue's code means. */
export interface ModelCompound {
  readonly id: string;
  readonly name?: string;
  readonly type?: string;
  readonly formula?: string;
  readonly parents?: readonly string[];
}

/** One end of a bond: a chain instance, the residue's position in its entity (a non-polymer's: 0) and the atom. */
export interface ModelBondEnd {
  readonly asym: string;
  readonly position: number;
  readonly atom?: string;
}

/** A bond between residues of the model (`struct_conn`: `covalent`, `disulfide`, `metal`, `glycosidic`). */
export interface ModelBond {
  readonly type: string;
  readonly a: ModelBondEnd;
  readonly b: ModelBondEnd;
}

/** A non-polymer entity's heavy atoms, from its first chain instance: its residue's position, name and element. */
export interface ModelAtoms {
  readonly entity: string;
  readonly atoms: readonly (readonly [position: number, name: string, element: string])[];
}

/** What the host's model holds beside its chain instances: entities, compounds, bonds, a molecule's atoms. */
export interface HostModel {
  readonly entities: readonly ModelEntity[];
  readonly compounds: readonly ModelCompound[];
  readonly bonds: readonly ModelBond[];
  readonly atoms: readonly ModelAtoms[];
}

/** The fake viewer's chain table: the entry's id and its chain instances of the deposited model, waters left out. */
export interface ChainTable {
  readonly entry: string;
  readonly chains: readonly ChainInstance[];
  /** The model the host holds beside the chains (the recorder writes it; a table of a test may have none). */
  readonly model?: HostModel;
}

const value = (s: string | undefined) => (s === undefined || s === '.' || s === '?' ? '' : s);

/** A residue's heavy atoms' B-factors summed: their mean is the model's per-residue value. */
interface BSum {
  sum: number;
  n: number;
}

const natural = (a: string, b: string) => a.localeCompare(b, 'en', { numeric: true });

/**
 * The chain table: polymer chains from `pdbx_poly_seq_scheme` (every residue of the entity, observed where `atom_site`
 * of the first model has it), non-polymers from `pdbx_nonpoly_scheme` (observed: deposited), waters left out; one
 * operator, the deposited model's (`1_555`). Each observed residue carries its mean heavy-atom B-factor, and `model`
 * what a host that holds the structure knows beside the chains (`hostModel`).
 */
export function chainTable(cif: Cif, entry: string): ChainTable {
  const water = new Set(
    cifRows(cif, '_entity')
      .filter((e) => e.type === 'water')
      .map((e) => e.id),
  );
  const modeled = new Set<string>();
  // a residue's heavy atoms' B-factors, each atom name once (the first alternate location), by chain and residue (a
  // non-polymer's one residue: its chain)
  const bs = new Map<string, BSum>();
  const seen = new Set<string>();
  const atoms = new Map<string, [position: number, name: string, element: string][]>();
  const atomsOf = new Map<string, string>();
  for (const a of cifRows(cif, '_atom_site')) {
    if (a.pdbx_PDB_model_num !== '1') continue;
    const polymer = a.label_seq_id !== '.';
    if (polymer) modeled.add(`${a.label_asym_id} ${a.label_seq_id}`);
    const residue = polymer ? `${a.label_asym_id} ${a.label_seq_id}` : a.label_asym_id;
    const atom = `${residue} ${a.label_atom_id}`;
    if (seen.has(atom)) continue;
    seen.add(atom);
    if (/^(H|D)$/i.test(a.type_symbol)) continue;
    const b = bs.get(residue) ?? { sum: 0, n: 0 };
    b.sum += Number(a.B_iso_or_equiv);
    b.n++;
    bs.set(residue, b);
    if (!polymer) {
      const rows = atoms.get(a.label_asym_id) ?? [];
      rows.push([0, a.label_atom_id, a.type_symbol]);
      atoms.set(a.label_asym_id, rows);
      atomsOf.set(a.label_asym_id, a.label_entity_id);
    }
  }
  const mean = (key: string): number | undefined => {
    const b = bs.get(key);
    return b?.n ? b.sum / b.n : undefined;
  };
  const chains = new Map<string, { auth: string; entity: string; residues: ChainResidue[] }>();
  const chainOf = (asym: string, auth: string, entity: string) => {
    let c = chains.get(asym);
    if (!c) chains.set(asym, (c = { auth, entity, residues: [] }));
    return c;
  };
  const withB = (r: ChainResidue, b: number | undefined): ChainResidue =>
    r.observed && b !== undefined ? { ...r, b } : r;
  // the author number is `pdb_seq_num`, which `atom_site.auth_seq_id` is (1FIN's ATP: 299; its `auth_seq_num` is 1)
  for (const r of cifRows(cif, '_pdbx_poly_seq_scheme'))
    chainOf(r.asym_id, r.pdb_strand_id, r.entity_id).residues.push(
      withB(
        {
          label: Number(r.seq_id),
          auth: value(r.pdb_seq_num) || value(r.auth_seq_num),
          ins: value(r.pdb_ins_code),
          comp: r.mon_id,
          observed: modeled.has(`${r.asym_id} ${r.seq_id}`),
        },
        mean(`${r.asym_id} ${r.seq_id}`),
      ),
    );
  for (const r of cifRows(cif, '_pdbx_nonpoly_scheme'))
    if (!water.has(r.entity_id))
      chainOf(r.asym_id, r.pdb_strand_id, r.entity_id).residues.push(
        withB(
          {
            label: null,
            auth: value(r.pdb_seq_num) || value(r.auth_seq_num),
            ins: value(r.pdb_ins_code),
            comp: r.mon_id,
            observed: true,
          },
          mean(r.asym_id),
        ),
      );
  const order = cifRows(cif, '_struct_asym').map((a) => a.id);
  const instances = order
    .filter((asym) => chains.has(asym))
    .map((asym): ChainInstance => {
      const c = chains.get(asym)!;
      return {
        label_asym_id: asym,
        auth_asym_id: c.auth,
        entity: c.entity,
        operators: ['1_555'],
        residues: c.residues,
      };
    });
  return { entry, chains: instances, model: hostModel(cif, instances, water, atoms, atomsOf) };
}

const BondTypes: Readonly<Record<string, string>> = { covale: 'covalent', disulf: 'disulfide', metalc: 'metal' };

/**
 * The model a host holds beside its chain instances, from the same file by rules of their own: the entities (id, type,
 * description, polymer type, canonical letters), the compounds the residues name (`chem_comp`, a modified residue's
 * parent from `pdbx_struct_mod_residue`), the bonds of `struct_conn` between the chains kept (a polymer's backbone
 * bond to its neighbor is the chain, not a bond), and each non-polymer entity's heavy atoms from its first chain.
 */
function hostModel(
  cif: Cif,
  chains: readonly ChainInstance[],
  water: ReadonlySet<string>,
  atoms: ReadonlyMap<string, readonly [number, string, string][]>,
  atomsOf: ReadonlyMap<string, string>,
): HostModel {
  const kept = new Set(chains.map((c) => c.entity));
  const poly = new Map(cifRows(cif, '_entity_poly').map((r) => [r.entity_id, r]));
  const entities = cifRows(cif, '_entity')
    .filter((e) => !water.has(e.id) && kept.has(e.id))
    .map((e): ModelEntity => {
      const p = poly.get(e.id);
      const letters = p?.pdbx_seq_one_letter_code_can?.replace(/\s+/g, '');
      return {
        id: e.id,
        type: e.type,
        ...(value(e.pdbx_description) ? { name: value(e.pdbx_description) } : {}),
        ...(p && value(p.type) ? { polymerType: p.type } : {}),
        ...(letters && /^[A-Za-z]+$/.test(letters) ? { letters } : {}),
      };
    });
  const used = new Set(chains.flatMap((c) => c.residues.map((r) => r.comp)));
  const modParent = new Map<string, string>();
  for (const r of cifRows(cif, '_pdbx_struct_mod_residue'))
    if (r.label_comp_id && r.parent_comp_id && !modParent.has(r.label_comp_id))
      modParent.set(r.label_comp_id, r.parent_comp_id);
  const compounds: ModelCompound[] = [];
  const seen = new Set<string>();
  for (const r of cifRows(cif, '_chem_comp')) {
    if (!used.has(r.id) || seen.has(r.id)) continue;
    seen.add(r.id);
    const own = value(r.mon_nstd_parent_comp_id);
    const parents = own ? own.split(',').map((x) => x.trim()) : modParent.has(r.id) ? [modParent.get(r.id)!] : [];
    compounds.push({
      id: r.id,
      ...(value(r.name) ? { name: value(r.name) } : {}),
      ...(value(r.type) ? { type: value(r.type) } : {}),
      ...(value(r.formula) ? { formula: value(r.formula) } : {}),
      ...(parents.length ? { parents } : {}),
    });
  }
  const entityOf = new Map(chains.map((c) => [c.label_asym_id, c.entity]));
  const polymer = new Set(entities.filter((e) => e.type === 'polymer').map((e) => e.id));
  const residueOf = (asym: string, auth: string, ins: string) =>
    chains.find((c) => c.label_asym_id === asym)?.residues.find((r) => r.auth === auth && r.ins === ins);
  const bonds: ModelBond[] = [];
  for (const r of cifRows(cif, '_struct_conn')) {
    const type = BondTypes[r.conn_type_id?.toLowerCase()];
    if (!type) continue;
    const end = (n: 1 | 2): ModelBondEnd | undefined => {
      const asym = r[`ptnr${n}_label_asym_id`];
      const symmetry = value(r[`ptnr${n}_symmetry`]);
      const entity = entityOf.get(asym);
      if (entity === undefined || (symmetry && symmetry !== '1_555')) return undefined;
      const atom = value(r[`ptnr${n}_label_atom_id`]);
      const seq = r[`ptnr${n}_label_seq_id`];
      if (polymer.has(entity)) {
        if (seq === '.' || seq === undefined || seq === '?') return undefined;
        return { asym, position: Number(seq) - 1, ...(atom ? { atom } : {}) };
      }
      return residueOf(asym, r[`ptnr${n}_auth_seq_id`], value(r[`pdbx_ptnr${n}_PDB_ins_code`]))
        ? { asym, position: 0, ...(atom ? { atom } : {}) }
        : undefined;
    };
    const a = end(1);
    const b = end(2);
    if (!a || !b) continue;
    const backbone =
      a.asym === b.asym &&
      polymer.has(entityOf.get(a.asym)!) &&
      type === 'covalent' &&
      Math.abs(a.position - b.position) === 1;
    if (backbone) continue;
    bonds.push({ type: type === 'covalent' && /glycosylation/i.test(value(r.pdbx_role)) ? 'glycosidic' : type, a, b });
  }
  // a non-polymer entity's atoms: its first chain instance, by natural order of `label_asym_id`
  const first = new Map<string, string>();
  for (const asym of [...atoms.keys()].sort(natural)) {
    const entity = atomsOf.get(asym)!;
    if (entityOf.has(asym) && !polymer.has(entity) && !first.has(entity)) first.set(entity, asym);
  }
  return {
    entities,
    compounds,
    bonds,
    atoms: [...first].map(([entity, asym]) => ({ entity, atoms: atoms.get(asym)! })),
  };
}

/** The order of a residue row's values in the committed table. */
export const ResidueFields = ['label', 'auth', 'ins', 'comp', 'observed', 'b'] as const;

/**
 * The chain table as committed: JSON, a chain instance a line, its residues rows of `ResidueFields` (`observed` 1 or
 * 0), so the file stays small and reads back with `readChainTable`.
 */
export function chainTableJson(t: ChainTable): string {
  const chain = (c: ChainInstance) =>
    JSON.stringify({
      ...c,
      residues: c.residues.map((r) => [r.label, r.auth, r.ins, r.comp, r.observed ? 1 : 0, r.b ?? null]),
    });
  const model = t.model ? `,"model":${JSON.stringify(t.model)}` : '';
  return `{"entry":${JSON.stringify(t.entry)},"fields":${JSON.stringify(ResidueFields)},"chains":[\n${t.chains.map(chain).join(',\n')}\n]${model}}\n`;
}

/** The committed table read back. */
export function readChainTable(text: string): ChainTable {
  const raw = JSON.parse(text) as {
    entry: string;
    chains: (Omit<ChainInstance, 'residues'> & {
      residues: [number | null, string, string, string, 0 | 1, (number | null)?][];
    })[];
    model?: HostModel;
  };
  return {
    entry: raw.entry,
    chains: raw.chains.map((c) => ({
      ...c,
      residues: c.residues.map(([label, auth, ins, comp, observed, b]) => ({
        label,
        auth,
        ins,
        comp,
        observed: observed === 1,
        ...(typeof b === 'number' ? { b } : {}),
      })),
    })),
    ...(raw.model ? { model: raw.model } : {}),
  };
}

/**
 * An entry's polymer entity as a `seq` file (2OCJ's, 3AVE's: one polymer entity each), named `system`: its canonical
 * letters (`entity_poly.pdbx_seq_one_letter_code_can`) and chain A's helices and strands (`struct_conf`,
 * `struct_sheet_range`, by `label_seq_id`) as ranges with a `category`.
 */
export function structureSeq(cif: Cif, system: string = StructureSystem): string {
  const [poly] = cifRows(cif, '_entity_poly');
  const letters = poly.pdbx_seq_one_letter_code_can.replace(/\s+/g, '');
  const rows: { start: number; end: number; kind: string }[] = [];
  for (const r of cifRows(cif, '_struct_conf'))
    if (r.beg_label_asym_id === 'A' && r.conf_type_id.startsWith('HELX'))
      rows.push({ start: Number(r.beg_label_seq_id) - 1, end: Number(r.end_label_seq_id), kind: 'helix' });
  for (const r of cifRows(cif, '_struct_sheet_range'))
    if (r.beg_label_asym_id === 'A')
      rows.push({ start: Number(r.beg_label_seq_id) - 1, end: Number(r.end_label_seq_id), kind: 'strand' });
  rows.sort((a, b) => a.start - b.start || a.end - b.end);
  const data = inlineSeq(
    { name: system, length: letters.length, alphabet: 'protein' },
    {
      letters: { form: 'letters', start: 0, letters },
      ranges: {
        form: 'list',
        start: Int32Array.from(rows, (e) => e.start),
        end: Int32Array.from(rows, (e) => e.end),
        id: categoryColumn(rows.map((e) => `${e.kind}${e.start + 1}`)),
        fields: [{ name: 'category', column: categoryColumn(rows.map((e) => e.kind)) }],
      },
    },
  );
  return `${encodeSeq(data)}\n`;
}

const Codes = 'ACDEFGHIKLMNPQRSTVWYX';
const code = (s: string) => Uint8Array.from(s.toUpperCase(), (c) => Math.max(0, Codes.indexOf(c)));

/**
 * The runs of aligned positions of `a` and `b` in their global alignment (core's Needleman–Wunsch: +4 a match, −1 a
 * mismatch, gaps 8 and 1, the ends free): what a `segments` map of the pair holds (`from` on `a`, `to` on `b`).
 */
export function alignedSegments(a: string, b: string): { from: number[]; to: number[]; length: number[] } {
  const { ops } = globalAlign(code(a), code(b), {
    match: 4,
    mismatch: -1,
    gapOpen: 8,
    gapExtend: 1,
    leadingGaps: 'free',
    trailingGaps: 'free',
    gapTies: 'open',
  });
  const out = { from: [] as number[], to: [] as number[], length: [] as number[] };
  let i = 0;
  let j = 0;
  let run = false;
  for (const op of ops) {
    if (op === AlignmentOp.Pair) {
      if (run) out.length[out.length.length - 1]++;
      else {
        out.from.push(i);
        out.to.push(j);
        out.length.push(1);
      }
      run = true;
      i++;
      j++;
    } else {
      run = false;
      if (op === AlignmentOp.OnlyA) i++;
      else j++;
    }
  }
  return out;
}
