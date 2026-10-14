/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `molstar:structure` complex provider (seqspec:libraries#the-mol-adapter): a structure's entities as Seq*
// molecules and its chain instances as units, with the `label`, `auth`, `atom` and SIFTS `uniprot:<acc>` numberings.

import {
  type ComplexProviderEntry,
  int,
  issueOf,
  list,
  literal,
  nullable,
  OptionalField,
  SimpleParamsSchema,
  str,
  union,
} from '@molstar/seq-core/format';
import type { MoleculeInput, NumberingInput, UnitInput } from '@molstar/seq-core/model';
import type { CifFrame } from '@molstar/io/reader/cif';
import type { Model, Structure } from '@molstar/model/model/structure';
import { type ChainInstance, chainsOf } from './chains.js';
import { siftsNumberings } from './sifts.js';
import { STRUCTURE_FORMATS, StructureSource, structureOf } from './source.js';

export type Include = 'polymer' | 'ligand' | 'ion' | 'branched';
export interface StructureOptions {
  model: number | 'all';
  assembly: string | null;
  include: Include[];
  chains: string[] | null;
  unit_name: string;
  numberings: ('label' | 'auth' | 'atom')[];
  default_numbering: string | null;
  sifts: 'best' | 'all' | 'none';
}
export const StructureOptionsSchema = SimpleParamsSchema({
  model: OptionalField(union(int, literal('all')), 0, 'the model, 0-based; all: units per model (poses, NMR)'),
  assembly: OptionalField(nullable(str), null, 'an assembly id; null: the deposited model'),
  include: OptionalField(
    list(literal('polymer', 'ligand', 'ion', 'branched')),
    ['polymer', 'ligand', 'ion', 'branched'],
    'what becomes a unit (water: LATER)',
  ),
  chains: OptionalField(nullable(list(str)), null, 'keep only the units with these names'),
  unit_name: OptionalField(str, 'label_asym_id', 'label_asym_id, auth_asym_id, or a template over their columns'),
  numberings: OptionalField(list(literal('label', 'auth', 'atom')), ['label', 'auth', 'atom'], 'which it writes'),
  default_numbering: OptionalField(nullable(str), 'auth', 'the one marked default on polymer units'),
  sifts: OptionalField(
    literal('best', 'all', 'none'),
    'best',
    "SIFTS uniprot:<acc> numberings: the best mapping's accessions, every accession (isoforms), or none",
  ),
});
const DEFAULTS: StructureOptions = {
  model: 0,
  assembly: null,
  include: ['polymer', 'ligand', 'ion', 'branched'],
  chains: null,
  unit_name: 'label_asym_id',
  numberings: ['label', 'auth', 'atom'],
  default_numbering: 'auth',
  sifts: 'best',
};

const ALPHABET = { protein: 'protein', DNA: 'dna', RNA: 'rna' } as const;

/** A unit's name by `unit_name`: a column (and `@instance` off the identity), or a `{column}` template as written. */
function nameOf(c: ChainInstance, template: string, model: number | null, compId: string): string {
  const cols: Record<string, string> = {
    label_asym_id: c.labelAsymId,
    auth_asym_id: c.authAsymId,
    entity_id: c.entityId,
    comp_id: compId,
    operator: c.instanceId,
    model: model === null ? '0' : `${model}`,
  };
  if (template === 'label_asym_id' || template === 'auth_asym_id') {
    const base = `${cols[template]}${c.identity ? '' : `@${c.instanceId}`}`;
    return model === null ? base : `${model}/${base}`;
  }
  return template.replace(/\{(\w+)\}/g, (_, k: string) => cols[k] ?? `{${k}}`);
}

/** Auth labels per entity position: atom_site's for observed residues, else `pdbx_poly_seq_scheme`'s `pdb_seq_num`. */
/** Insertion codes that are not letters, read as none: the first (chain, number, code) and how many. */
interface OddCodes {
  first: string | null;
  /** Per chain instance and position, once. */
  readonly seen: Set<string>;
}
const LETTERS = /^[A-Za-z]+$/;

function authLabels(
  c: ChainInstance,
  length: number,
  first: number,
  frame: CifFrame | null,
  odd: OddCodes,
): (number | string | null)[] {
  const out: (number | string | null)[] = new Array(length).fill(null);
  // A code that is not letters would join its number into another (52 and "1" read 521): none, and reported.
  const label = (num: string | number, ins: string, k: number) => {
    if (!ins || ins === '?' || ins === '.') return Number(num);
    if (LETTERS.test(ins)) return `${num}${ins}`;
    odd.first ??= `chain ${c.authAsymId} residue ${num} has insertion code "${ins}"`;
    odd.seen.add(`${c.labelAsymId}@${c.instanceId}:${k}`);
    return Number(num);
  };
  const scheme = frame?.categories.pdbx_poly_seq_scheme;
  if (scheme) {
    const [asym, seq, num, ins] = ['asym_id', 'seq_id', 'pdb_seq_num', 'pdb_ins_code'].map((f) => scheme.getField(f));
    for (let i = 0; asym && seq && num && i < scheme.rowCount; i++) {
      if (asym.str(i) !== c.labelAsymId || num.valueKind(i) !== 0) continue;
      const k = seq.int(i) - first;
      if (k >= 0 && k < length) out[k] = label(num.str(i), ins?.valueKind(i) === 0 ? ins.str(i) : '', k);
    }
  }
  const res = c.unit.model.atomicHierarchy.residues;
  for (let x = 0; x < c.count; x++) {
    const k = c.labels[x] - first;
    if (k >= 0 && k < length)
      out[k] = label(res.auth_seq_id.value(c.residues[x]), res.pdbx_PDB_ins_code.value(c.residues[x]), k);
  }
  return out;
}

/** Runs of present positions, [start, stop). */
function runsOf(present: Uint8Array): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < present.length; ) {
    if (!present[i]) {
      i++;
      continue;
    }
    let j = i;
    while (j < present.length && present[j]) j++;
    out.push([i, j]);
    i = j;
  }
  return out;
}

export interface BuiltComplex {
  molecules: MoleculeInput[];
  units: UnitInput[];
}

/** The complex of one or more structures (one per model with `model: "all"`), as `molstar:structure` yields it. */
export function complexInputs(
  structures: readonly { structure: Structure; model: number | null }[],
  frame: CifFrame | null,
  options: Partial<StructureOptions>,
  report: (code: string, path: string, message: string) => void,
  ligandFile = false,
): BuiltComplex {
  const o = { ...DEFAULTS, ...options };
  const molecules = new Map<string, MoleculeInput>();
  const entities = new Map<string, number>();
  const units: { input: UnitInput; order: number }[] = [];
  const names = new Map<string, number>();
  const odd: OddCodes = { first: null, seen: new Set() };
  const predicted = !!frame?.categories.ma_qa_metric_local;
  for (const { structure, model } of structures) {
    for (const c of chainsOf(structure).chains) {
      const m: Model = c.unit.model;
      const e = m.entities.getEntityIndex(c.entityId);
      if (c.kind === 'water') continue;
      const isIon = c.kind === 'ligand' && c.count === 1;
      const included = c.kind === 'polymer' ? 'polymer' : c.kind === 'branched' ? 'branched' : isIon ? 'ion' : 'ligand';
      if (!o.include.includes(included) || c.count === 0) continue;
      const compId = m.atomicHierarchy.atoms.label_comp_id.value(c.unit.elements[c.start[0]]);
      // An SDF or MOL2 record is a model with one ligand unit, `ligand`, over `entity:ligand`.
      const name = ligandFile ? (model === null ? 'ligand' : `${model}/ligand`) : nameOf(c, o.unit_name, model, compId);
      if (o.chains && !o.chains.includes(name)) continue;
      names.set(name, (names.get(name) ?? 0) + 1);
      const mkey = ligandFile ? 'entity:ligand' : `entity:${c.entityId}`;
      const description = m.entities.data.pdbx_description.value(e)?.[0];
      let mol = molecules.get(mkey);
      const seq = c.kind === 'polymer' ? m.sequence.byEntityKey[e]?.sequence : undefined;
      const first = seq && seq.length > 0 ? seq.seqId.value(0) : 1;
      if (!mol) {
        if (seq) {
          let residues = '';
          for (let i = 0; i < seq.length; i++) residues += seq.code.value(i) || 'X';
          const numberings: NumberingInput[] = o.numberings.includes('label')
            ? [{ name: 'label', method: 'offset', first }]
            : [];
          const sifts =
            o.sifts === 'none'
              ? { numberings: [], accessions: [] }
              : siftsNumberings(frame, c.entityId, first, o.sifts === 'all');
          mol = {
            key: mkey,
            alphabet: ALPHABET[seq.kind as keyof typeof ALPHABET] ?? 'generic',
            residues,
            label: description ?? mkey,
            numberings: [...numberings, ...sifts.numberings],
            identifiers: sifts.accessions.map((value) => ({ namespace: 'uniprot', value })),
          };
        } else {
          const atoms = m.atomicHierarchy.atoms;
          const el = (x: number) => c.unit.elements[c.start[x]];
          mol = {
            key: mkey,
            kind: 'ligand',
            elements: Array.from({ length: c.count }, (_, x) => atoms.type_symbol.value(el(x)) as string),
            chemCompIds: Array.from({ length: c.count }, (_, x) => atoms.label_comp_id.value(el(x))),
            label: c.kind === 'branched' ? (description ?? mkey) : compId,
          };
        }
        molecules.set(mkey, mol);
        entities.set(mkey, e);
      }
      const unit: UnitInput = {
        key: name,
        molecule: mkey,
        origin: predicted ? 'predicted' : 'experimental',
        custom: { mvs: { label_asym_id: c.labelAsymId, auth_asym_id: c.authAsymId, instance_id: c.instanceId } },
      };
      if (seq) {
        const present = new Uint8Array(seq.length);
        for (let x = 0; x < c.count; x++) {
          const k = c.labels[x] - first;
          if (k >= 0 && k < seq.length) present[k] = 1;
        }
        const observed = runsOf(present);
        if (observed.length !== 1 || observed[0][0] !== 0 || observed[0][1] !== seq.length) unit.observed = observed;
        if (o.numberings.includes('auth')) {
          const labels = authLabels(c, seq.length, first, frame, odd);
          unit.numberings = [{ name: 'auth', method: 'labels', labels, isDefault: o.default_numbering === 'auth' }];
        }
      } else {
        const len = mol.elements?.length ?? 0;
        if (c.count < len) unit.observed = [[0, c.count]];
        if (o.numberings.includes('atom')) {
          const atoms = m.atomicHierarchy.atoms;
          const names: (string | null)[] = Array.from({ length: len }, (_, x) =>
            x < c.count ? atoms.label_atom_id.value(c.unit.elements[c.start[x]]) : null,
          );
          unit.numberings = [{ name: 'atom', method: 'names', names, isDefault: true }];
        }
      }
      units.push({ input: unit, order: e });
    }
  }
  for (const [name, n] of names) {
    if (n > 1) {
      report(
        'molstar.unit_name.collision',
        '/options/unit_name',
        `"${o.unit_name}" names ${n} units "${name}": write a template that tells them apart, such as "{label_asym_id}" or "{auth_asym_id}-{entity_id}"`,
      );
      return { molecules: [], units: [] };
    }
  }
  if (odd.seen.size > 0) {
    const more = odd.seen.size > 1 ? ` (${odd.seen.size} residues in all)` : '';
    report(
      'molstar.ins_code.invalid',
      '/source',
      `${odd.first}, not letters${more}: read with no code, so its auth label is the number alone and may repeat its neighbor's; select it by label_seq_id (numbering "label") to tell them apart`,
    );
  }
  // Entity order, then chain order (a stable sort keeps the units' order within an entity).
  units.sort((a, b) => a.order - b.order);
  return { molecules: [...molecules.values()], units: units.map((u) => u.input) };
}

/** `molstar:structure`, over a `StructureSource` its formats give (or a host's `Loaders.parsed`). */
export const StructureProvider: ComplexProviderEntry = {
  name: 'molstar:structure',
  options: StructureOptionsSchema,
  formats: STRUCTURE_FORMATS,
  async build(ctx, source, options) {
    const src = source.kind === 'object' ? source.value : null;
    if (!StructureSource.is(src)) throw new Error('not a structure parsed by @molstar/seq-molstar-poc');
    const o = { ...DEFAULTS, ...(options as Partial<StructureOptions>) };
    const models = o.model === 'all' ? Array.from({ length: src.trajectory.frameCount }, (_, i) => i) : [o.model];
    const structures = await Promise.all(
      models.map(async (m) => ({
        structure: await structureOf(src, m, o.assembly),
        model: o.model === 'all' ? m : null,
      })),
    );
    const ligandFile = src.format === 'sdf' || src.format === 'mol2';
    const report = (code: string, path: string, message: string) =>
      ctx.report(issueOf(ctx.registries, code, path, message));
    return complexInputs(structures, src.frame, o, report, ligandFile);
  },
};
