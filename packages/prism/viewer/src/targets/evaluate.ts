/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import {
  StructureProperties as P,
  QueryContext,
  type Structure,
  StructureElement,
  StructureSelection,
  Unit,
} from '@molstar/model/model/structure';
import { compile } from '@molstar/model/script/runtime/query/compiler';
import type { Expression } from '@molstar/query-language/language/expression';
import { compactNumbers } from './addresses.js';
import { failure, messageOf } from './compile.js';
import { compileTarget } from './target.js';
import { selectorExpression } from './selector.js';
import type {
  AtomRow,
  ChainRow,
  CompileOptions,
  CompiledTarget,
  Failure,
  Level,
  MoleculeRow,
  ProjectedRow,
  Projection,
  ResidueRow,
  SchemaItems,
  Selected,
  SelectedChain,
  SelectionBounds,
  SelectorT,
  Target,
} from './types.js';

// Evaluating a compiled target against a structure (targets.md §2): one evaluator, one set of semantics, whoever asked.
// Ported from Stories Untold's `any-molql` (evaluate.ts).

/** The atoms a MolQL expression selects in a structure (Mol* loci, with the structure's own units). */
export function queryLoci(structure: Structure, molql: Expression): StructureElement.Loci {
  const query = compile<StructureSelection>(molql);
  const selection = query(new QueryContext(structure));
  return StructureSelection.toLociWithSourceUnits(selection);
}

/** The atoms a selector selects in a structure. */
export function selectLoci(structure: Structure, selector: SelectorT): StructureElement.Loci {
  return queryLoci(structure, selectorExpression(selector));
}

export interface Evaluated {
  ok: true;
  loci: StructureElement.Loci;
  /** What it selected (computed with the loci, so no second evaluation is needed). */
  selected: Selected;
}

/** A MolQL expression evaluated: the loci and what they hold. Failures of the query (not of the text) come back as a value. */
export function evaluateExpression(
  structure: Structure,
  molql: Expression,
  options: SummarizeOptions = {},
): Evaluated | Failure {
  try {
    const loci = queryLoci(structure, molql);
    return { ok: true, loci, selected: summarize(loci, options) };
  } catch (e) {
    return failure(`Mol* could not run the query: ${messageOf(e)}`, 'Check the target; try a simpler one.');
  }
}

export type EvaluatedTarget = (Evaluated & { target: Extract<CompiledTarget, { ok: true }> }) | Failure;

/** A target compiled and evaluated against a structure. */
export function evaluateTarget(
  structure: Structure,
  target: Target,
  options: CompileOptions & SummarizeOptions = {},
): EvaluatedTarget {
  const compiled = compileTarget(target, options);
  if (!compiled.ok) return compiled;
  const evaluated = evaluateExpression(structure, compiled.molql, options);
  return evaluated.ok ? { ...evaluated, target: compiled } : evaluated;
}

/** The bounds of what a selection selects (Mol*'s boundary of the loci); undefined when it selects nothing. */
export function selectionBounds(loci: StructureElement.Loci, margin = 0): SelectionBounds | undefined {
  const atoms = StructureElement.Loci.size(loci);
  if (atoms === 0) return undefined;
  const { box, sphere } = StructureElement.Loci.getBoundary(loci);
  const v = (a: ArrayLike<number>, d: number): [number, number, number] => [
    (a[0] as number) + d,
    (a[1] as number) + d,
    (a[2] as number) + d,
  ];
  return {
    box: { min: v(box.min, -margin), max: v(box.max, margin) },
    sphere: { center: v(sphere.center, 0), radius: sphere.radius + margin },
    atoms,
  };
}

const count = (n: number, one: string, many = `${one}s`) => `${n.toLocaleString('en-US')} ${n === 1 ? one : many}`;
const round = (v: number) => Math.round(v * 10) / 10;

/** "in chain A", "in chains A and B", "in chains A, B and C", "in 12 chains". */
function inChains(chains: readonly string[]): string {
  if (chains.length === 0) return '';
  if (chains.length === 1) return ` in chain ${chains[0]}`;
  if (chains.length > 6) return ` in ${chains.length} chains`;
  return ` in chains ${chains.slice(0, -1).join(', ')} and ${chains[chains.length - 1]}`;
}

function describe(atoms: number, residues: number, chains: readonly string[]): string {
  return atoms === 0 ? 'selects nothing' : `${count(residues, 'residue')} (${count(atoms, 'atom')})${inChains(chains)}`;
}

export interface SummarizeOptions {
  /** At most this many ranges per chain are listed in `Selected` (default 50); the rest is counted in `more`. */
  maxRangesPerChain?: number;
}

interface ChainAccumulator {
  chain: string;
  residues: number;
  atoms: number;
  /** Distinct author residues (a residue of an assembly is in several units). */
  numbers: Map<string, { seq: number; ins: string }>;
}

/**
 * What a selection selected (targets.md §2): atom count, residue count, author chains, and the residues by author
 * number per chain as compact ranges. Walks residues, not atoms, in atomic units.
 */
export function summarize(loci: StructureElement.Loci, options: SummarizeOptions = {}): Selected {
  const limit = options.maxRangesPerChain ?? 50;
  const chains = new Map<string, ChainAccumulator>();
  const instances = new Set<string>();
  let atoms = 0;
  let residues = 0;
  const location = StructureElement.Location.create(loci.structure);
  const add = (n: number) => {
    const chainId = P.chain.auth_asym_id(location);
    let chain = chains.get(chainId);
    if (!chain) {
      chain = { chain: chainId, residues: 0, atoms: 0, numbers: new Map() };
      chains.set(chainId, chain);
    }
    const seq = P.residue.auth_seq_id(location);
    const ins = P.residue.pdbx_PDB_ins_code(location) || '';
    chain.residues++;
    chain.atoms += n;
    residues++;
    const key = `${seq}|${ins}`;
    if (!chain.numbers.has(key)) chain.numbers.set(key, { seq, ins });
  };
  for (const { unit, indices } of loci.elements) {
    const size = OrderedSet.size(indices);
    if (size === 0) continue;
    atoms += size;
    location.unit = unit;
    instances.add(unit.conformation.operator.instanceId);
    if (Unit.isAtomic(unit)) {
      const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index;
      let previous = -1;
      let inResidue = 0;
      for (let i = 0; i < size; i++) {
        const element = unit.elements[OrderedSet.getAt(indices, i)];
        const residue = residueIndex[element];
        if (residue !== previous) {
          if (inResidue) add(inResidue);
          previous = residue;
          inResidue = 0;
          location.element = element;
        }
        inResidue++;
      }
      if (inResidue) add(inResidue);
    } else {
      // Coarse elements (spheres, gaussians) are their own residue-like pieces.
      for (let i = 0; i < size; i++) {
        location.element = unit.elements[OrderedSet.getAt(indices, i)];
        add(1);
      }
    }
  }
  const byChain: SelectedChain[] = [];
  const addresses: string[] = [];
  for (const acc of chains.values()) {
    const all = compactNumbers([...acc.numbers.values()]);
    const ranges = all.slice(0, limit);
    byChain.push({
      chain: acc.chain,
      residues: acc.residues,
      atoms: acc.atoms,
      ranges,
      more: all.length - ranges.length,
    });
    for (const range of ranges) addresses.push(`${acc.chain}:${range}`);
  }
  const chainIds = byChain.map((c) => c.chain);
  return {
    atoms,
    residues,
    chains: chainIds,
    byChain,
    addresses,
    instances: instances.size > 1 ? [...instances] : [],
    line: describe(atoms, residues, chainIds),
  };
}

/**
 * What a selection selects, projected to a level: atoms, residues, chains or molecules (entities), with counts at
 * every level and up to `limit` rows at the one asked for. (Ported from Stories; `summarize` is the cheaper account.)
 */
export function project(loci: StructureElement.Loci, level: Level = 'residues', limit = 80): Projection {
  const atoms: AtomRow[] = [];
  const residues = new Map<string, ResidueRow>();
  const chains = new Map<string, ChainRow & { min?: number; max?: number; residueKeys: Set<string> }>();
  const molecules = new Map<string, MoleculeRow & { residueKeys: Set<string> }>();
  let atomCount = 0;
  StructureElement.Loci.forEachLocation(loci, (l) => {
    atomCount++;
    const chain = P.chain.auth_asym_id(l);
    const seq = P.residue.auth_seq_id(l);
    const ins = P.residue.pdbx_PDB_ins_code(l);
    const comp = P.atom.label_comp_id(l);
    const residueAddress = `${chain}:${seq}${ins ?? ''}`;
    // Units of one chain under several operators (assemblies) are told apart by the unit id.
    const unitKey = `${l.unit.id}|${residueAddress}`;
    const residue = residues.get(unitKey);
    if (residue) residue.atoms++;
    else residues.set(unitKey, { kind: 'residue', address: residueAddress, comp, atoms: 1 });
    if (level === 'atoms' && atoms.length < limit) {
      atoms.push({
        kind: 'atom',
        address: `${residueAddress}@${P.atom.auth_atom_id(l)}`,
        element: String(P.atom.type_symbol(l)),
        comp,
        xyz: [round(P.atom.x(l)), round(P.atom.y(l)), round(P.atom.z(l))],
      });
    }
    const label = P.chain.label_asym_id(l);
    const entity = P.entity.id(l);
    const molecule = P.entity.pdbx_description(l).join(', ') || comp;
    const chainKey = `${l.unit.id}|${chain}|${label}`;
    let c = chains.get(chainKey);
    if (!c) {
      c = { kind: 'chain', auth: chain, label, entity, molecule, residues: 0, atoms: 0, residueKeys: new Set() };
      chains.set(chainKey, c);
    }
    c.atoms++;
    if (!c.residueKeys.has(unitKey)) {
      c.residueKeys.add(unitKey);
      c.residues++;
      c.min = c.min === undefined ? seq : Math.min(c.min, seq);
      c.max = c.max === undefined ? seq : Math.max(c.max, seq);
    }
    let m = molecules.get(entity);
    if (!m) {
      m = {
        kind: 'molecule',
        entity,
        molecule,
        type: P.entity.type(l),
        chains: [],
        residues: 0,
        atoms: 0,
        residueKeys: new Set(),
      };
      molecules.set(entity, m);
    }
    m.atoms++;
    if (!m.chains.includes(chain)) m.chains.push(chain);
    if (!m.residueKeys.has(unitKey)) {
      m.residueKeys.add(unitKey);
      m.residues++;
    }
  });
  const chainIds = [...new Set([...chains.values()].map((c) => c.auth))];
  const all: ProjectedRow[] =
    level === 'atoms'
      ? atoms
      : level === 'residues'
        ? [...residues.values()]
        : level === 'chains'
          ? [...chains.values()].map(({ residueKeys: _k, min, max, ...c }) => ({
              ...c,
              ...(min !== undefined && max !== undefined
                ? { range: `${c.auth}:${min}${max !== min ? `-${max}` : ''}` }
                : {}),
            }))
          : [...molecules.values()].map(({ residueKeys: _k, ...m }) => m);
  const total = level === 'atoms' ? atomCount : all.length;
  const rows = all.slice(0, limit);
  const counts = { atoms: atomCount, residues: residues.size, chains: chains.size };
  return {
    level,
    ...counts,
    chainIds,
    rows,
    more: Math.max(0, total - rows.length),
    line: describe(atomCount, counts.residues, chainIds),
  };
}

/**
 * What a selection holds, as MVS schema items (Struct of Arrays): one row per residue, keyed by `label_asym_id`,
 * `auth_seq_id` and `pdbx_PDB_ins_code` (the column only when the file uses insertion codes), with `instance_id` in an
 * assembly (several operators); one row per atom (`atom_id`) when the selection cuts through residues. `auth_asym_id`
 * and `label_comp_id` (and `auth_atom_id`) are there to name rows for people, not to key them.
 */
export function selectionItems(loci: StructureElement.Loci): SchemaItems {
  const whole = StructureElement.Loci.extendToWholeResidues(loci);
  const atomLevel = StructureElement.Loci.size(whole) !== StructureElement.Loci.size(loci);
  const instances = new Set(loci.structure.units.map((u) => u.conformation.operator.instanceId));
  const assembly = instances.size > 1;
  const cols: Required<Pick<SchemaItems, 'label_asym_id' | 'auth_seq_id' | 'auth_asym_id' | 'label_comp_id'>> &
    SchemaItems = { label_asym_id: [], auth_seq_id: [], pdbx_PDB_ins_code: [], auth_asym_id: [], label_comp_id: [] };
  if (assembly) cols.instance_id = [];
  if (atomLevel) {
    cols.atom_id = [];
    cols.auth_atom_id = [];
  }
  let insertions = false;
  const seen = new Set<string>();
  StructureElement.Loci.forEachLocation(loci, (l) => {
    const ins = P.residue.pdbx_PDB_ins_code(l) || null;
    const instance = P.unit.instance_id(l);
    if (!atomLevel) {
      const key = `${instance}|${P.chain.label_asym_id(l)}|${P.residue.auth_seq_id(l)}|${ins ?? ''}`;
      if (seen.has(key)) return;
      seen.add(key);
    }
    if (ins) insertions = true;
    cols.label_asym_id.push(P.chain.label_asym_id(l));
    cols.auth_seq_id.push(P.residue.auth_seq_id(l));
    cols.pdbx_PDB_ins_code?.push(ins);
    cols.auth_asym_id.push(P.chain.auth_asym_id(l));
    cols.label_comp_id.push(P.atom.label_comp_id(l));
    cols.instance_id?.push(instance);
    cols.atom_id?.push(P.atom.id(l));
    cols.auth_atom_id?.push(P.atom.auth_atom_id(l));
  });
  if (!insertions) delete cols.pdbx_PDB_ins_code;
  return cols;
}
