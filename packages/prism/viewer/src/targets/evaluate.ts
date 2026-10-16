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
import type { Expression } from '@molstar/query-language/language/expression';
import { compactNumbers } from './addresses.js';
import { compiledQuery, failure, messageOf } from './compile.js';
import { compileTarget } from './target.js';
import { selectorExpression } from './selector.js';
import type { CompileOptions, CompiledTarget, Failure, Selected, SelectedChain, SelectorT, Target } from './types.js';

// Evaluating a compiled target against a structure (targets.md §2): one evaluator, one set of semantics, whoever asked.
// Ported from Stories Untold's `any-molql` (evaluate.ts).

/** The atoms a MolQL expression selects in a structure (Mol* loci, with the structure's own units). */
export function queryLoci(structure: Structure, molql: Expression): StructureElement.Loci {
  const selection = compiledQuery(molql)(new QueryContext(structure));
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

const count = (n: number, one: string, many = `${one}s`) => `${n.toLocaleString('en-US')} ${n === 1 ? one : many}`;

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
  const piece = (chainId: string, n: number) => {
    let chain = chains.get(chainId);
    if (!chain) {
      chain = { chain: chainId, residues: 0, atoms: 0, numbers: new Map() };
      chains.set(chainId, chain);
    }
    chain.residues++;
    chain.atoms += n;
    residues++;
    return chain;
  };
  const number = (chain: ChainAccumulator, seq: number, ins: string) => {
    const key = `${seq}|${ins}`;
    if (!chain.numbers.has(key)) chain.numbers.set(key, { seq, ins });
  };
  const add = (n: number) => {
    const chain = piece(P.chain.auth_asym_id(location), n);
    number(chain, P.residue.auth_seq_id(location), P.residue.pdbx_PDB_ins_code(location) || '');
  };
  // A coarse element has no author chain or numbers: its label chain, and the residues it covers (label numbers).
  const addCoarse = () => {
    const chain = piece(P.chain.label_asym_id(location), 1);
    const end = P.coarse.seq_id_end(location);
    for (let seq = P.coarse.seq_id_begin(location); seq <= end; seq++) number(chain, seq, '');
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
        addCoarse();
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
