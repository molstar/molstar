/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues, type Issue, issueOf } from '@molstar/core/repr/tree';
import type { CifCategory, CifField, CifFrame } from '@molstar/io/reader/cif/data-model';
import type { EntityType, StructureDescription } from '../io/structure.js';
import { natural } from '../io/structure-read.js';

// An mmCIF frame read into a structure description (_spec/README.md "Data", mmCIF; molecules.md, Reading): the CIF
// categories as `@molstar/io`'s parser gives them, read by field name (no schema: its size), into the columns
// `structureNodes` turns into nodes. Only what the description holds is read: entities and their components, the
// units of the deposited model or of the assembly asked, their modeled residues from `atom_site` (the first model,
// the first `label_alt_id` an atom name has), links from `struct_conn` and the branches' own, compounds from
// `chem_comp`, a molecule's atoms from its first unit and the per-residue values (`b_factor`, `occupancy`, `plddt`).

/** The mmCIF runtime's own codes (the rest are `seq.io.structure.*`, said by `structureNodes`, and `seq.source.*`). */
export const MmcifIssues = defineIssues('seq.io.mmcif', {
  'seq.io.mmcif.assembly': {
    severity: 'error',
    description: 'The assembly asked is not in the file: the source is not read',
    hint: 'Name an assembly of `pdbx_struct_assembly_gen`, or none for the deposited model',
  },
  'seq.io.mmcif.altloc': {
    severity: 'info',
    description:
      'Atoms at alternate locations: the first `label_alt_id` an atom name has is read, the others dropped (a position is one atom)',
  },
  'seq.io.mmcif.links': {
    severity: 'info',
    description:
      '`struct_conn` rows not read as links: an end on a symmetry mate, or on a residue no unit read has (a water when waters are not read)',
  },
});

/** How the file is read (molecules.md, Reading; the source params). */
export interface MmcifReadOptions {
  /** `pdbx_struct_assembly_gen.assembly_id`; absent: the deposited model (units without an operator). */
  readonly assembly?: string;
  /** Water entities and their units: left out unless true. */
  readonly water?: boolean;
  /** A molecule's name: `description` (`entity.pdbx_description`, absent) or a branched one's `composition`. */
  readonly name?: 'description' | 'composition';
  /** Units made at most (an assembly's operators times its chains): one past it lets `structureNodes` say the bound. */
  readonly maxUnits: number;
}

export interface MmcifRead {
  readonly description: StructureDescription;
  readonly issues: readonly Issue[];
}

/** The assembly asked is not in the file. */
export class MmcifAssemblyError extends Error {}

/** The link types read from `struct_conn.conn_type_id` (molecules.md, Links). */
const ConnTypes: Readonly<Record<string, string>> = { covale: 'covalent', disulf: 'disulfide', metalc: 'metal' };
const EntityTypes = new Set<string>(['polymer', 'non-polymer', 'branched', 'water']);
/** `Column.ValueKinds.Present` (a const enum of `@molstar/core/data/db`, whose module this one does not load). */
const Present = 0;

type Field = CifField | undefined;
const fieldOf = (c: CifCategory | undefined, name: string): Field => c?.getField(name);
const has = (f: Field, row: number) => f !== undefined && f.valueKind(row) === Present;
const strOf = (f: Field, row: number) => (has(f, row) ? f!.str(row) : undefined);

/** The groups of an operator expression (`1`, `1,2`, `(1-60)`, `(X0)(1-5)`): each a list of operator ids. */
export function operatorGroups(expression: string): string[][] {
  const groups = expression.includes('(') ? [...expression.matchAll(/\(([^)]*)\)/g)].map((m) => m[1]) : [expression];
  return groups.map((g) =>
    g
      .split(',')
      .map((s) => s.trim())
      .filter((s) => s)
      .flatMap((s) => {
        const m = /^(-?\d+)-(-?\d+)$/.exec(s);
        if (!m) return [s];
        const [a, b] = [Number(m[1]), Number(m[2])];
        const out: string[] = [];
        for (let i = a; i <= b && out.length <= 1_000_000; i++) out.push(String(i));
        return out;
      }),
  );
}

/** The operators of an expression as Mol\*'s instance ids (`ASM-1`, `ASM-X0-5`), at most `max`. */
export function assemblyOperators(expression: string, max: number): string[] {
  const groups = operatorGroups(expression);
  let combos: string[][] = [[]];
  for (const g of groups) {
    const next: string[][] = [];
    for (const c of combos) for (const id of g) if (next.length < max) next.push([...c, id]);
    combos = next;
  }
  return combos.filter((c) => c.length).map((c) => `ASM-${c.join('-')}`);
}

/** Residues of one asym, by position: their author numbers and the values read from their atoms. */
interface Residue {
  readonly position: number;
  readonly auth: number;
  readonly ins: string;
  /** The atom names seen, with the `label_alt_id` each read. */
  readonly alt: Map<string, string>;
  b: number;
  occ: number;
  /** Heavy atoms the values are summed over. */
  n: number;
}

/**
 * Reads a frame into a structure description and what reading it said; throws `MmcifAssemblyError` for an assembly
 * the file lacks. Pure: the frame's categories by name, nothing else.
 */
export function mmcifDescription(frame: CifFrame, options: MmcifReadOptions, address: string): MmcifRead {
  const issues: Issue[] = [];
  const cats = frame.categories;

  // entities and their polymer types and canonical letters
  const ent = cats.entity;
  const entityIds: string[] = [];
  const entityTypes: EntityType[] = [];
  const entityNames: (string | undefined)[] = [];
  const entityRow = new Map<string, number>();
  for (let i = 0; i < (ent?.rowCount ?? 0); i++) {
    const id = strOf(fieldOf(ent, 'id'), i);
    if (id === undefined || entityRow.has(id)) continue;
    const t = strOf(fieldOf(ent, 'type'), i)?.toLowerCase() ?? 'non-polymer';
    entityRow.set(id, entityIds.length);
    entityIds.push(id);
    entityTypes.push((EntityTypes.has(t) ? t : 'non-polymer') as EntityType);
    entityNames.push(strOf(fieldOf(ent, 'pdbx_description'), i));
  }
  const polymerType: (string | undefined)[] = entityIds.map(() => undefined);
  const letters: (string | undefined)[] = entityIds.map(() => undefined);
  const poly = cats.entity_poly;
  for (let i = 0; i < (poly?.rowCount ?? 0); i++) {
    const e = entityRow.get(strOf(fieldOf(poly, 'entity_id'), i) ?? '');
    if (e === undefined) continue;
    polymerType[e] = strOf(fieldOf(poly, 'type'), i);
    const can = strOf(fieldOf(poly, 'pdbx_seq_one_letter_code_can'), i)?.replace(/\s+/g, '');
    if (can && /^[A-Za-z]+$/.test(can)) letters[e] = can;
  }

  // components a position: a polymer's `entity_poly_seq`, a branch's list, a non-polymer's one component
  const compEntity: number[] = [];
  const compPosition: number[] = [];
  const compCode: string[] = [];
  /** An entity's sequence numbers (`num`) to positions, in their order. */
  const positionOf: Map<number, number>[] = entityIds.map(() => new Map());
  const listed = (catName: string, codeName: string) => {
    const c = cats[catName];
    if (!c) return;
    const [eF, nF, kF] = [fieldOf(c, 'entity_id'), fieldOf(c, 'num'), fieldOf(c, codeName)];
    const rows: { e: number; num: number; code: string }[] = [];
    for (let i = 0; i < c.rowCount; i++) {
      const e = entityRow.get(strOf(eF, i) ?? '');
      const code = strOf(kF, i);
      if (e === undefined || code === undefined || !has(nF, i)) continue;
      rows.push({ e, num: nF!.int(i), code });
    }
    for (const r of rows) positionOf[r.e].set(r.num, 0);
    for (const m of positionOf) {
      const nums = [...m.keys()].sort((a, b) => a - b);
      nums.forEach((n, k) => m.set(n, k));
    }
    for (const r of rows) {
      compEntity.push(r.e);
      compPosition.push(positionOf[r.e].get(r.num)!);
      compCode.push(r.code);
    }
  };
  listed('entity_poly_seq', 'mon_id');
  listed('pdbx_entity_branch_list', 'comp_id');
  const nonpoly = cats.pdbx_entity_nonpoly;
  const nonpolyCode = new Map<number, string>();
  for (let i = 0; i < (nonpoly?.rowCount ?? 0); i++) {
    const e = entityRow.get(strOf(fieldOf(nonpoly, 'entity_id'), i) ?? '');
    const code = strOf(fieldOf(nonpoly, 'comp_id'), i);
    if (e !== undefined && code !== undefined && !nonpolyCode.has(e)) nonpolyCode.set(e, code);
  }

  // units: the deposited model's chains, or the assembly's chains under its operators
  const asymEntity = new Map<string, number>();
  const asymOrder: string[] = [];
  const sa = cats.struct_asym;
  for (let i = 0; i < (sa?.rowCount ?? 0); i++) {
    const id = strOf(fieldOf(sa, 'id'), i);
    const e = entityRow.get(strOf(fieldOf(sa, 'entity_id'), i) ?? '');
    if (id === undefined || e === undefined || asymEntity.has(id)) continue;
    asymEntity.set(id, e);
    asymOrder.push(id);
  }
  const keep = (e: number) => entityTypes[e] !== 'water' || !!options.water;

  // atom_site: the first model, its residues a unit's (modeled) positions with their author numbers and values
  const as = cats.atom_site;
  const nA = as?.rowCount ?? 0;
  const A = {
    model: fieldOf(as, 'pdbx_PDB_model_num'),
    asym: fieldOf(as, 'label_asym_id'),
    entity: fieldOf(as, 'label_entity_id'),
    seq: fieldOf(as, 'label_seq_id'),
    comp: fieldOf(as, 'label_comp_id'),
    atom: fieldOf(as, 'label_atom_id'),
    alt: fieldOf(as, 'label_alt_id'),
    element: fieldOf(as, 'type_symbol'),
    auth: fieldOf(as, 'auth_seq_id'),
    ins: fieldOf(as, 'pdbx_PDB_ins_code'),
    authAsym: fieldOf(as, 'auth_asym_id'),
    b: fieldOf(as, 'B_iso_or_equiv'),
    occ: fieldOf(as, 'occupancy'),
  };
  // a branched asym's residues by their author number (`atom_site.auth_seq_id` is the scheme's `pdb_seq_num`)
  const branchNum = new Map<string, number>();
  const bs = cats.pdbx_branch_scheme;
  for (let i = 0; i < (bs?.rowCount ?? 0); i++) {
    const asym = strOf(fieldOf(bs, 'asym_id'), i);
    const seqNum = strOf(fieldOf(bs, 'pdb_seq_num'), i) ?? strOf(fieldOf(bs, 'auth_seq_num'), i);
    const num = fieldOf(bs, 'num');
    if (asym !== undefined && seqNum !== undefined && has(num, i)) branchNum.set(`${asym} ${seqNum}`, num!.int(i));
  }
  const firstModel = nA && has(A.model, 0) ? A.model!.int(0) : undefined;
  const residues = new Map<string, Map<number, Residue>>();
  const authAsym = new Map<string, string>();
  /** A molecule's atoms are its first unit's (natural order of `label_asym_id`). */
  const atomsAsym = new Map<number, string>();
  for (const id of [...asymOrder].sort(natural)) {
    const e = asymEntity.get(id)!;
    if (entityTypes[e] !== 'polymer' && keep(e) && !atomsAsym.has(e)) atomsAsym.set(e, id);
  }
  const atomEntity: number[] = [];
  const atomResidue: number[] = [];
  const atomName: string[] = [];
  const atomElement: string[] = [];
  /** A water asym's residues in their order (each water a position). */
  const waterIndex = new Map<string, Map<string, number>>();
  let dropped = 0;
  let current: { key: string; residue: Residue | undefined } = { key: '', residue: undefined };
  for (let i = 0; i < nA; i++) {
    if (firstModel !== undefined && has(A.model, i) && A.model!.int(i) !== firstModel) continue;
    const asym = strOf(A.asym, i);
    if (asym === undefined) continue;
    let e = asymEntity.get(asym);
    if (e === undefined) {
      // no struct_asym row: the atom's own entity
      e = entityRow.get(strOf(A.entity, i) ?? '');
      if (e === undefined) continue;
      asymEntity.set(asym, e);
      asymOrder.push(asym);
    }
    if (!keep(e)) continue;
    if (!authAsym.has(asym)) {
      const a = strOf(A.authAsym, i);
      if (a !== undefined) authAsym.set(asym, a);
    }
    const auth = has(A.auth, i) ? A.auth!.int(i) : Number.NaN;
    const ins = strOf(A.ins, i) ?? '';
    const type = entityTypes[e];
    let position: number | undefined;
    if (type === 'polymer') position = has(A.seq, i) ? positionOf[e].get(A.seq!.int(i)) : undefined;
    else if (type === 'branched') {
      const num = branchNum.get(`${asym} ${strOf(A.auth, i)}`);
      position = num === undefined ? undefined : positionOf[e].get(num);
    } else if (type === 'water') {
      let m = waterIndex.get(asym);
      if (!m) waterIndex.set(asym, (m = new Map()));
      const k = `${auth}${ins}`;
      if (!m.has(k)) m.set(k, m.size);
      position = m.get(k);
    } else {
      position = 0;
      if (!nonpolyCode.has(e)) {
        const code = strOf(A.comp, i);
        if (code !== undefined) nonpolyCode.set(e, code);
      }
    }
    if (position === undefined) continue;
    const key = `${asym} ${position}`;
    if (current.key !== key) {
      let byPosition = residues.get(asym);
      if (!byPosition) residues.set(asym, (byPosition = new Map()));
      let r = byPosition.get(position);
      if (!r) {
        r = { position, auth: Number.isNaN(auth) ? position + 1 : auth, ins, alt: new Map(), b: 0, occ: 0, n: 0 };
        byPosition.set(position, r);
      }
      current = { key, residue: r };
    }
    const r = current.residue!;
    const name = strOf(A.atom, i) ?? '';
    const alt = strOf(A.alt, i) ?? '';
    const seen = r.alt.get(name);
    if (seen === undefined) r.alt.set(name, alt);
    else {
      if (seen !== alt) dropped++;
      continue;
    }
    const element = strOf(A.element, i) ?? '';
    const heavy = !/^(H|D)$/i.test(element);
    if (heavy) {
      r.b += has(A.b, i) ? A.b!.float(i) : Number.NaN;
      r.occ += has(A.occ, i) ? A.occ!.float(i) : Number.NaN;
      r.n++;
    }
    if (atomsAsym.get(e) === asym) {
      atomEntity.push(e);
      atomResidue.push(position);
      atomName.push(name);
      atomElement.push(element);
    }
  }
  if (dropped)
    issues.push(
      issueOf(MmcifIssues, 'seq.io.mmcif.altloc', `${dropped} atoms at a second alternate location dropped`, {
        address,
      }),
    );

  // a non-polymer's component, a water's one per water of its longest unit
  const waterLength = new Map<number, number>();
  for (const [asym, m] of waterIndex) {
    const e = asymEntity.get(asym)!;
    waterLength.set(e, Math.max(waterLength.get(e) ?? 0, m.size));
  }
  for (const [e, code] of nonpolyCode) {
    const n = entityTypes[e] === 'water' ? (waterLength.get(e) ?? 0) : 1;
    for (let p = 0; p < n; p++) {
      compEntity.push(e);
      compPosition.push(p);
      compCode.push(code);
    }
  }

  // identifiers: a non-polymer's CCD code, a branch's WURCS and LINUCS, a PRD's id
  const idEntity: number[] = [];
  const idNamespace: string[] = [];
  const idValue: string[] = [];
  const identifier = (e: number, namespace: string, value: string) => {
    for (let k = 0; k < idEntity.length; k++)
      if (idEntity[k] === e && idNamespace[k] === namespace && idValue[k] === value) return;
    idEntity.push(e);
    idNamespace.push(namespace);
    idValue.push(value);
  };
  for (const [e, code] of nonpolyCode) if (entityTypes[e] === 'non-polymer') identifier(e, 'ccd', code);
  const bd = cats.pdbx_entity_branch_descriptor;
  for (let i = 0; i < (bd?.rowCount ?? 0); i++) {
    const e = entityRow.get(strOf(fieldOf(bd, 'entity_id'), i) ?? '');
    const type = strOf(fieldOf(bd, 'type'), i)?.toLowerCase();
    const value = strOf(fieldOf(bd, 'descriptor'), i)?.trim();
    if (e !== undefined && value && (type === 'wurcs' || type === 'linucs')) identifier(e, type, value);
  }
  const pm = cats.pdbx_molecule;
  for (let i = 0; i < (pm?.rowCount ?? 0); i++) {
    const e = asymEntity.get(strOf(fieldOf(pm, 'asym_id'), i) ?? '');
    const prd = strOf(fieldOf(pm, 'prd_id'), i);
    if (e !== undefined && prd) identifier(e, 'prd', prd);
  }

  // the units
  const units: { asym: string; operator?: string }[] = [];
  const max = options.maxUnits + 1;
  if (options.assembly === undefined) {
    for (const asym of asymOrder) if (keep(asymEntity.get(asym)!) && units.length < max) units.push({ asym });
  } else {
    const gen = cats.pdbx_struct_assembly_gen;
    const seen = new Set<string>();
    let found = false;
    for (let i = 0; i < (gen?.rowCount ?? 0); i++) {
      if (strOf(fieldOf(gen, 'assembly_id'), i) !== options.assembly) continue;
      found = true;
      const operators = assemblyOperators(strOf(fieldOf(gen, 'oper_expression'), i) ?? '', max);
      const asyms = (strOf(fieldOf(gen, 'asym_id_list'), i) ?? '').split(',').map((s) => s.trim());
      for (const operator of operators)
        for (const asym of asyms) {
          const e = asymEntity.get(asym);
          const k = `${asym} ${operator}`;
          if (e === undefined || !keep(e) || seen.has(k) || units.length >= max) continue;
          seen.add(k);
          units.push({ asym, operator });
        }
    }
    if (!found) throw new MmcifAssemblyError(`The assembly ${JSON.stringify(options.assembly)} is not in the file`);
  }
  const unitRow = new Map<string, number>();
  units.forEach((u, k) => unitRow.set(`${u.asym} ${u.operator ?? ''}`, k));
  const operatorsOf = new Map<string, (string | undefined)[]>();
  for (const u of units) {
    const ops = operatorsOf.get(u.asym) ?? [];
    ops.push(u.operator);
    operatorsOf.set(u.asym, ops);
  }

  // residues a unit, with their values
  const plddt = plddtOf(frame, positionOf, asymEntity);
  const resChain: number[] = [];
  const resPosition: number[] = [];
  const resAuth: number[] = [];
  const resIns: string[] = [];
  const bFactor: number[] = [];
  const occupancy: number[] = [];
  const confidence: number[] = [];
  units.forEach((u, k) => {
    const byPosition = residues.get(u.asym);
    if (!byPosition) return;
    for (const r of [...byPosition.values()].sort((a, b) => a.position - b.position)) {
      resChain.push(k);
      resPosition.push(r.position);
      resAuth.push(r.auth);
      resIns.push(r.ins);
      bFactor.push(r.n ? r.b / r.n : Number.NaN);
      occupancy.push(r.n ? r.occ / r.n : Number.NaN);
      confidence.push(plddt?.get(`${u.asym} ${r.position}`) ?? Number.NaN);
    }
  });

  // links: `struct_conn`'s within one copy, and each branch's own
  const link = {
    aChain: [] as number[],
    aPosition: [] as number[],
    bChain: [] as number[],
    bPosition: [] as number[],
    type: [] as string[],
    aAtom: [] as (string | null)[],
    bAtom: [] as (string | null)[],
  };
  const addLink = (a: string, ap: number, b: string, bp: number, type: string, aa?: string, ba?: string) => {
    for (const op of operatorsOf.get(a) ?? []) {
      const ar = unitRow.get(`${a} ${op ?? ''}`);
      const br = unitRow.get(`${b} ${op ?? ''}`);
      if (ar === undefined || br === undefined) continue;
      link.aChain.push(ar);
      link.aPosition.push(ap);
      link.bChain.push(br);
      link.bPosition.push(bp);
      link.type.push(type);
      link.aAtom.push(aa ?? null);
      link.bAtom.push(ba ?? null);
    }
  };
  const byAuth = (asym: string, auth: string | undefined, ins: string | undefined) => {
    for (const r of residues.get(asym)?.values() ?? []) if (`${r.auth}` === auth && r.ins === (ins ?? '')) return r;
    return undefined;
  };
  const sc = cats.struct_conn;
  let unread = 0;
  for (let i = 0; i < (sc?.rowCount ?? 0); i++) {
    const conn = strOf(fieldOf(sc, 'conn_type_id'), i)?.toLowerCase() ?? '';
    if (!Object.hasOwn(ConnTypes, conn)) continue;
    const end = (n: 1 | 2) => {
      const asym = strOf(fieldOf(sc, `ptnr${n}_label_asym_id`), i);
      const symmetry = strOf(fieldOf(sc, `ptnr${n}_symmetry`), i);
      if (asym === undefined || (symmetry !== undefined && symmetry !== '1_555')) return undefined;
      const e = asymEntity.get(asym);
      if (e === undefined) return undefined;
      const seq = fieldOf(sc, `ptnr${n}_label_seq_id`);
      const position =
        entityTypes[e] === 'polymer' && has(seq, i)
          ? positionOf[e].get(seq!.int(i))
          : byAuth(
              asym,
              strOf(fieldOf(sc, `ptnr${n}_auth_seq_id`), i),
              strOf(fieldOf(sc, `pdbx_ptnr${n}_PDB_ins_code`), i),
            )?.position;
      if (position === undefined) return undefined;
      return { asym, e, position, atom: strOf(fieldOf(sc, `ptnr${n}_label_atom_id`), i) };
    };
    const a = end(1);
    const b = end(2);
    if (!a || !b || !operatorsOf.has(a.asym) || !operatorsOf.has(b.asym)) {
      unread++;
      continue;
    }
    // a branch's own links are read from its entity's list, below; a polymer's backbone bond (a modified residue's
    // to its neighbor) is the chain itself, not a link
    if (a.asym === b.asym && entityTypes[a.e] === 'branched') continue;
    if (
      a.asym === b.asym &&
      entityTypes[a.e] === 'polymer' &&
      conn === 'covale' &&
      Math.abs(a.position - b.position) === 1
    )
      continue;
    const role = strOf(fieldOf(sc, 'pdbx_role'), i) ?? '';
    const type = conn === 'covale' && /glycosylation/i.test(role) ? 'glycosidic' : ConnTypes[conn];
    addLink(a.asym, a.position, b.asym, b.position, type, a.atom, b.atom);
  }
  if (unread)
    issues.push(
      issueOf(MmcifIssues, 'seq.io.mmcif.links', `${unread} struct_conn rows not read as links`, { address }),
    );
  const bl = cats.pdbx_entity_branch_link;
  for (let i = 0; i < (bl?.rowCount ?? 0); i++) {
    const e = entityRow.get(strOf(fieldOf(bl, 'entity_id'), i) ?? '');
    const [n1, n2] = [fieldOf(bl, 'entity_branch_list_num_1'), fieldOf(bl, 'entity_branch_list_num_2')];
    if (e === undefined || !has(n1, i) || !has(n2, i)) continue;
    const p1 = positionOf[e].get(n1!.int(i));
    const p2 = positionOf[e].get(n2!.int(i));
    if (p1 === undefined || p2 === undefined) continue;
    for (const asym of asymOrder) {
      if (asymEntity.get(asym) !== e) continue;
      const [a1, a2] = [strOf(fieldOf(bl, 'atom_id_1'), i), strOf(fieldOf(bl, 'atom_id_2'), i)];
      addLink(asym, p1, asym, p2, 'glycosidic', a1, a2);
    }
  }

  // compounds: those of the entities read, their parents from `chem_comp` or `pdbx_struct_mod_residue`
  const used = new Set<string>();
  for (let k = 0; k < compCode.length; k++) if (keep(compEntity[k])) used.add(compCode[k]);
  const modParent = new Map<string, string>();
  const mr = cats.pdbx_struct_mod_residue;
  for (let i = 0; i < (mr?.rowCount ?? 0); i++) {
    const code = strOf(fieldOf(mr, 'label_comp_id'), i);
    const parent = strOf(fieldOf(mr, 'parent_comp_id'), i);
    if (code && parent && !modParent.has(code)) modParent.set(code, parent);
  }
  const compounds = { id: [] as string[], name: [] as (string | undefined)[], type: [] as (string | undefined)[] };
  const parents: (string[] | undefined)[] = [];
  const cc = cats.chem_comp;
  for (let i = 0; i < (cc?.rowCount ?? 0); i++) {
    const id = strOf(fieldOf(cc, 'id'), i);
    if (id === undefined || !used.has(id) || compounds.id.includes(id)) continue;
    compounds.id.push(id);
    compounds.name.push(strOf(fieldOf(cc, 'name'), i));
    compounds.type.push(strOf(fieldOf(cc, 'type'), i));
    const own = strOf(fieldOf(cc, 'mon_nstd_parent_comp_id'), i);
    const p = own ? own.split(',').map((s) => s.trim()) : modParent.has(id) ? [modParent.get(id)!] : undefined;
    parents.push(p?.length ? p : undefined);
  }

  // names: the entity's description, or a branch's composition in list order (`NAG2 BMA MAN2`)
  const names = entityIds.map((_, e) =>
    options.name === 'composition' && entityTypes[e] === 'branched'
      ? composition(e, compEntity, compPosition, compCode)
      : entityNames[e],
  );

  const entry = strOf(fieldOf(cats.entry, 'id'), 0);
  const description: StructureDescription = {
    ...(entry !== undefined ? { name: entry } : {}),
    entities: { id: entityIds, type: entityTypes, polymerType, name: names, letters },
    components: {
      entity: Int32Array.from(compEntity),
      position: Int32Array.from(compPosition),
      component: compCode,
    },
    identifiers: { entity: Int32Array.from(idEntity), namespace: idNamespace, value: idValue },
    chains: {
      label_asym_id: units.map((u) => u.asym),
      auth_asym_id: units.map((u) => authAsym.get(u.asym)),
      entity: Int32Array.from(units, (u) => asymEntity.get(u.asym)!),
      operator: units.map((u) => u.operator),
    },
    residues: {
      chain: Int32Array.from(resChain),
      position: Int32Array.from(resPosition),
      auth_seq_id: Int32Array.from(resAuth),
      ins_code: resIns,
    },
    links: {
      a_chain: Int32Array.from(link.aChain),
      a_position: Int32Array.from(link.aPosition),
      b_chain: Int32Array.from(link.bChain),
      b_position: Int32Array.from(link.bPosition),
      type: link.type,
      a_atom: link.aAtom,
      b_atom: link.bAtom,
    },
    compounds: { ...compounds, parents },
    atoms: {
      entity: Int32Array.from(atomEntity),
      residue: Int32Array.from(atomResidue),
      name: atomName,
      element: atomElement,
    },
    values: {
      b_factor: Float32Array.from(bFactor),
      occupancy: Float32Array.from(occupancy),
      ...(plddt ? { plddt: Float32Array.from(confidence) } : {}),
    },
  };
  return { description, issues };
}

/** A branch's components in list order, a run of one code written once with its count (`NAG2 BMA MAN2`). */
function composition(e: number, entity: number[], position: number[], code: string[]): string {
  const parts: string[] = [];
  let last = '';
  let n = 0;
  const flush = () => last && parts.push(n > 1 ? `${last}${n}` : last);
  for (let k = 0; k < entity.length; k++) {
    if (entity[k] !== e || (k > 0 && entity[k - 1] === e && position[k - 1] === position[k])) continue;
    if (code[k] === last) n++;
    else {
      flush();
      last = code[k];
      n = 1;
    }
  }
  flush();
  return parts.join(' ');
}

/**
 * A model's per-residue pLDDT (ModelCIF's `ma_qa_metric_local` of the metric whose type is `pLDDT`), by asym and
 * position; undefined when the file has none.
 */
function plddtOf(
  frame: CifFrame,
  positionOf: readonly Map<number, number>[],
  asymEntity: ReadonlyMap<string, number>,
): Map<string, number> | undefined {
  const metric = frame.categories.ma_qa_metric;
  const local = frame.categories.ma_qa_metric_local;
  if (!metric || !local) return undefined;
  const ids = new Set<string>();
  for (let i = 0; i < metric.rowCount; i++)
    if (
      /^plddt$/i.test(strOf(fieldOf(metric, 'type'), i) ?? '') ||
      /plddt/i.test(strOf(fieldOf(metric, 'name'), i) ?? '')
    )
      ids.add(strOf(fieldOf(metric, 'id'), i) ?? '');
  if (!ids.size) return undefined;
  const out = new Map<string, number>();
  const [asymF, seqF, idF, valueF] = ['label_asym_id', 'label_seq_id', 'metric_id', 'metric_value'].map((n) =>
    fieldOf(local, n),
  );
  for (let i = 0; i < local.rowCount; i++) {
    if (!ids.has(strOf(idF, i) ?? '') || !has(seqF, i) || !has(valueF, i)) continue;
    const asym = strOf(asymF, i) ?? '';
    const e = asymEntity.get(asym);
    const p = e === undefined ? undefined : positionOf[e].get(seqF!.int(i));
    if (p !== undefined) out.set(`${asym} ${p}`, valueF!.float(i));
  }
  return out;
}
