/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '@molstar/core/data/db';
import { SortedArray } from '@molstar/core/data/int';
import type { ChainIndex, ElementIndex, ResidueIndex } from '../../../model/indexing.js';
import type { Model } from '../../../model.js';
import type { CoarseElements } from '../../../model/properties/coarse.js';
import { filterInPlace, range, sortIfNeeded } from '@molstar/core/util/array';
import { type Mapping, MultiMap, NumberMap } from './utils.js';

export interface IndicesAndSortings {
  atomic?: AtomicIndicesAndSortings;
  spheres?: CoarseIndicesAndSortings;
  gaussians?: CoarseIndicesAndSortings;
}

/** The things the indices (and anything else derived from the hierarchy and the atom ids, but not from the coordinates) depend on.
 * Models with the same `atomicHierarchy`, `coarseHierarchy` and `atomId` column share one scope: these are the frames of a trajectory
 * (a copy of a model with other coordinates), so what is computed for one frame is valid for all of them. */
export interface HierarchyScope {
  readonly coarseHierarchy: Model['coarseHierarchy'];
  readonly atomId: Model['atomicConformation']['atomId'];
  /** Cached `IndicesAndSortings.get` */
  indices?: IndicesAndSortings;
}

const HierarchyScopes = new WeakMap<object, HierarchyScope[]>();
/** The scope of each model, so that a lookup does not have to search the scopes of its hierarchy (an ensemble of models which share a hierarchy
 * but have their own `atomId` columns has one per model) */
const ModelScopes = new WeakMap<Model, HierarchyScope>();

export const HierarchyScope = {
  /** Get the scope of a model (created if the model is the first with its hierarchies and atom ids). Does not keep the model alive. */
  get(model: Model): HierarchyScope {
    let scope = ModelScopes.get(model);
    if (!scope) ModelScopes.set(model, (scope = HierarchyScope.find(model)));
    return scope;
  },

  find(model: Model): HierarchyScope {
    const { atomicHierarchy, coarseHierarchy } = model;
    const atomId = model.atomicConformation.atomId;
    let scopes = HierarchyScopes.get(atomicHierarchy);
    if (scopes) {
      for (const s of scopes) {
        if (s.coarseHierarchy === coarseHierarchy && s.atomId === atomId) return s;
      }
    } else {
      HierarchyScopes.set(atomicHierarchy, (scopes = []));
    }
    const scope: HierarchyScope = { coarseHierarchy, atomId };
    scopes.push(scope);
    return scope;
  },
};

export const IndicesAndSortings = {
  /** Get `IndicesAndSortings` for a model (use a cached value or create if not available yet).
   * The cache is shared by the models of a trajectory (see `HierarchyScope`), the indices are created once. */
  get(model: Model): IndicesAndSortings {
    const scope = HierarchyScope.get(model);
    return (scope.indices ??= this.create(model));
  },

  /** Create `IndicesAndSortings` for a model */
  create(model: Model): IndicesAndSortings {
    return {
      atomic: createAtomicIndicesAndSortings(model),
      spheres: createCoarseIndicesAndSortings(model.coarseHierarchy.spheres),
      gaussians: createCoarseIndicesAndSortings(model.coarseHierarchy.gaussians),
    };
  },
};

/** Auxiliary data structure for efficiently finding chains/residues/atoms in an atomic model by their properties */
export interface AtomicIndicesAndSortings {
  chainsByLabelEntityId: Mapping<string, readonly ChainIndex[]>;
  chainsByLabelAsymId: Mapping<string, readonly ChainIndex[]>;
  chainsByAuthAsymId: Mapping<string, readonly ChainIndex[]>;
  residuesSortedByLabelSeqId: Mapping<ChainIndex, Sorting<ResidueIndex, number>>;
  residuesSortedByAuthSeqId: Mapping<ChainIndex, Sorting<ResidueIndex, number>>;
  residuesSortedBySourceIndex: Mapping<ChainIndex, Sorting<ResidueIndex, number>>;
  residuesByInsCode: Mapping<ChainIndex, Mapping<string, readonly ResidueIndex[]>>;
  residuesByLabelCompId: Mapping<ChainIndex, Mapping<string, readonly ResidueIndex[]>>;
  /** Indicates if each residue is listed only once in `residuesByLabelCompId` (i.e. if each residue has only one label_comp_id) */
  residuesByLabelCompIdIsPure: boolean;
  residuesByAuthCompId: Mapping<ChainIndex, Mapping<string, readonly ResidueIndex[]>>;
  /** Indicates if each residue is listed only once in `residuesByAuthCompId` (i.e. if each residue has only one auth_comp_id) */
  residuesByAuthCompIdIsPure: boolean;
  /** The first atom with each `atom_id` (ids are not necessarily unique, e.g. serial numbers wrap around in some formats) */
  atomsById: Mapping<number, ElementIndex>;
  /** All atoms (ascending) of each `atom_id` which is used by more than one atom, other ids are not listed */
  atomsWithRepeatedId: Mapping<number, readonly ElementIndex[]>;
  atomsBySourceIndex: Mapping<number, ElementIndex>;
  /** Lookups over all chains at once (for rows without chain criteria), each built on first use */
  global: GlobalResidueIndices;
}

/** Create `AtomicIndicesAndSortings` for a model */
function createAtomicIndicesAndSortings(model: Model): AtomicIndicesAndSortings | undefined {
  const h = model.atomicHierarchy;
  const nAtoms = h.atoms._rowCount;
  if (nAtoms === 0) return undefined;

  const nChains = h.chains._rowCount;
  const { label_entity_id, label_asym_id, auth_asym_id } = h.chains;
  const { label_seq_id, auth_seq_id, pdbx_PDB_ins_code } = h.residues;
  const { label_comp_id, auth_comp_id } = h.atoms;
  const { Present } = Column.ValueKind;

  const chainsByLabelEntityId = new MultiMap<string, ChainIndex>();
  const chainsByLabelAsymId = new MultiMap<string, ChainIndex>();
  const chainsByAuthAsymId = new MultiMap<string, ChainIndex>();
  const residuesSortedByLabelSeqId = new Map<ChainIndex, Sorting<ResidueIndex, number>>();
  const residuesSortedByAuthSeqId = new Map<ChainIndex, Sorting<ResidueIndex, number>>();
  const residuesSortedBySourceIndex = new Map<ChainIndex, Sorting<ResidueIndex, number>>();
  const residuesByInsCode = new Map<ChainIndex, MultiMap<string, ResidueIndex>>();
  const residuesByLabelCompId = new Map<ChainIndex, MultiMap<string, ResidueIndex>>();
  let residuesByLabelCompIdIsPure = true;
  const residuesByAuthCompId = new Map<ChainIndex, MultiMap<string, ResidueIndex>>();
  let residuesByAuthCompIdIsPure = true;
  const atomsById = new NumberMap<number, ElementIndex>(nAtoms + 1);
  const atomsBySourceIndex = new NumberMap<number, ElementIndex>(nAtoms);
  const atomsWithRepeatedId = new Map<number, ElementIndex[]>();

  const _labelCompIdSet = new Set<string>();
  const _authCompIdSet = new Set<string>();

  for (let iChain = 0 as ChainIndex; iChain < nChains; iChain++) {
    chainsByLabelEntityId.add(label_entity_id.value(iChain), iChain);
    chainsByLabelAsymId.add(label_asym_id.value(iChain), iChain);
    chainsByAuthAsymId.add(auth_asym_id.value(iChain), iChain);

    const iResFrom = h.residueAtomSegments.index[h.chainAtomSegments.offsets[iChain]];
    const iResTo = h.residueAtomSegments.index[h.chainAtomSegments.offsets[iChain + 1] - 1] + 1;

    const residuesWithLabelSeqId = filterInPlace(
      range(iResFrom, iResTo) as ResidueIndex[],
      (iRes) => label_seq_id.valueKind(iRes) === Present,
    );
    residuesSortedByLabelSeqId.set(iChain, Sorting.create(residuesWithLabelSeqId, label_seq_id.value));

    const residuesWithAuthSeqId = filterInPlace(
      range(iResFrom, iResTo) as ResidueIndex[],
      (iRes) => auth_seq_id.valueKind(iRes) === Present,
    );
    residuesSortedByAuthSeqId.set(iChain, Sorting.create(residuesWithAuthSeqId, auth_seq_id.value));

    const residuesWithSourceIndex = range(iResFrom, iResTo) as ResidueIndex[];
    residuesSortedBySourceIndex.set(iChain, Sorting.create(residuesWithSourceIndex, h.residueSourceIndex.value));

    const residuesHereByInsCode = new MultiMap<string, ResidueIndex>();
    const residuesHereByLabelCompId = new MultiMap<string, ResidueIndex>();
    const residuesHereByAuthCompId = new MultiMap<string, ResidueIndex>();
    for (let iRes = iResFrom; iRes < iResTo; iRes++) {
      if (pdbx_PDB_ins_code.valueKind(iRes) === Present) {
        residuesHereByInsCode.add(pdbx_PDB_ins_code.value(iRes), iRes);
      }
      const iAtomFrom = h.residueAtomSegments.offsets[iRes];
      const iAtomTo = h.residueAtomSegments.offsets[iRes + 1];
      for (let iAtom = iAtomFrom; iAtom < iAtomTo; iAtom++) {
        _labelCompIdSet.add(label_comp_id.value(iAtom));
        _authCompIdSet.add(auth_comp_id.value(iAtom));
      }
      if (_labelCompIdSet.size > 1) residuesByLabelCompIdIsPure = false;
      if (_authCompIdSet.size > 1) residuesByAuthCompIdIsPure = false;
      for (const labelCompId of _labelCompIdSet) residuesHereByLabelCompId.add(labelCompId, iRes);
      for (const authCompId of _authCompIdSet) residuesHereByAuthCompId.add(authCompId, iRes);
      _labelCompIdSet.clear();
      _authCompIdSet.clear();
    }
    residuesByInsCode.set(iChain, residuesHereByInsCode);
    residuesByLabelCompId.set(iChain, residuesHereByLabelCompId);
    residuesByAuthCompId.set(iChain, residuesHereByAuthCompId);
  }

  const atomId = model.atomicConformation.atomId.value;
  const atomIndex = h.atomSourceIndex.value;
  for (let iAtom = 0 as ElementIndex; iAtom < nAtoms; iAtom++) {
    const id = atomId(iAtom);
    const first = atomsById.get(id);
    if (first === undefined) {
      atomsById.set(id, iAtom);
    } else {
      let repeated = atomsWithRepeatedId.get(id);
      if (!repeated) atomsWithRepeatedId.set(id, (repeated = [first]));
      repeated.push(iAtom);
    }
    atomsBySourceIndex.set(atomIndex(iAtom), iAtom);
  }

  return {
    chainsByLabelEntityId,
    chainsByLabelAsymId,
    chainsByAuthAsymId,
    residuesSortedByLabelSeqId,
    residuesSortedByAuthSeqId,
    residuesSortedBySourceIndex,
    residuesByInsCode,
    residuesByLabelCompId,
    residuesByLabelCompIdIsPure,
    residuesByAuthCompId,
    residuesByAuthCompIdIsPure,
    atomsById,
    atomsWithRepeatedId,
    atomsBySourceIndex,
    global: new GlobalResidueIndices(h),
  };
}

/** Lookups of residues over the whole model, for rows without chain criteria: the per-chain indices would need a lookup in every chain
 * (cost proportional to the number of chains even when there is a single match). Unlike them, these are built from one pass
 * over the residues without sorting, and only when first needed. Residues in each list are ascending. */
export class GlobalResidueIndices {
  private _bySourceIndex: Map<number, ResidueIndex> | null | undefined = undefined;
  private _byInsCode?: MultiMap<string, ResidueIndex>;
  private _byLabelCompId?: MultiMap<string, ResidueIndex>;
  private _byAuthCompId?: MultiMap<string, ResidueIndex>;

  constructor(private readonly h: Model['atomicHierarchy']) {}

  /** `residue_index` -> residue, or `null` when `residue_index` is not unique in the model (then use the per-chain indices) */
  get residueBySourceIndex(): Map<number, ResidueIndex> | null {
    if (this._bySourceIndex === undefined) {
      const { residueSourceIndex } = this.h;
      const nResidues = this.h.residues._rowCount;
      const map = new Map<number, ResidueIndex>();
      let unique = true;
      for (let iRes = 0 as ResidueIndex; iRes < nResidues; iRes++) {
        const index = residueSourceIndex.value(iRes);
        if (map.has(index)) {
          unique = false;
          break;
        }
        map.set(index, iRes);
      }
      this._bySourceIndex = unique ? map : null;
    }
    return this._bySourceIndex;
  }

  get residuesByInsCode(): Mapping<string, readonly ResidueIndex[]> {
    if (!this._byInsCode) {
      const { pdbx_PDB_ins_code } = this.h.residues;
      const map = new MultiMap<string, ResidueIndex>();
      const nResidues = this.h.residues._rowCount;
      for (let iRes = 0 as ResidueIndex; iRes < nResidues; iRes++) {
        if (pdbx_PDB_ins_code.valueKind(iRes) === Column.ValueKind.Present)
          map.add(pdbx_PDB_ins_code.value(iRes), iRes);
      }
      this._byInsCode = map;
    }
    return this._byInsCode;
  }

  /** Residues having an atom with the given `label_comp_id` */
  get residuesByLabelCompId(): Mapping<string, readonly ResidueIndex[]> {
    return (this._byLabelCompId ??= this.byCompId(this.h.atoms.label_comp_id));
  }

  /** Residues having an atom with the given `auth_comp_id` */
  get residuesByAuthCompId(): Mapping<string, readonly ResidueIndex[]> {
    return (this._byAuthCompId ??= this.byCompId(this.h.atoms.auth_comp_id));
  }

  private byCompId(compId: Column<string>) {
    const map = new MultiMap<string, ResidueIndex>();
    const { offsets } = this.h.residueAtomSegments;
    const nResidues = this.h.residues._rowCount;
    const seen = new Set<string>();
    for (let iRes = 0 as ResidueIndex; iRes < nResidues; iRes++) {
      for (let iAtom = offsets[iRes], stop = offsets[iRes + 1]; iAtom < stop; iAtom++) seen.add(compId.value(iAtom));
      for (const id of seen) map.add(id, iRes);
      seen.clear();
    }
    return map;
  }
}

/** Auxiliary data structure for efficiently finding chains/elements in a coarse model by their properties */
export interface CoarseIndicesAndSortings {
  /** Coarse equivalent to `model.atomicHierarchy.chains` */
  chains: {
    /** Number of chains */
    count: number;
    /** Maps chain index to `label_entity_id` value */
    label_entity_id: string[];
    /** Maps chain index to `label_asym_id` value */
    label_asym_id: string[];
  };
  chainsByEntityId: Mapping<string, readonly ChainIndex[]>;
  chainsByAsymId: Mapping<string, readonly ChainIndex[]>;
  /** Coarse elements (per chain) sorted by `seq_id_begin`.
   * This is used to get the range of elements which may overlap with a certain seq_id interval.
   *
   * (Filtering coarse elements by seq_id range is an interval search problem, so the worst-case-efficient solution would be to use a data structure optimized for that.
   * But that would be overkill if we expect that in most cases the coarse elements cover non-overlapping seq_id ranges.
   * So the current solution should be sufficient (fast for non-overlapping elements, while still correct if there are overlaps).) */
  elementsSortedBySeqIdBegin: Mapping<
    ChainIndex,
    Sorting<ElementIndex, number> & {
      /** Non-decreasing upper bound for `seq_id_end` values of elements as listed in `keys` (`seq_id_end.value(keys[i]) <= endUpperBounds[i]`) */
      endUpperBounds: SortedArray;
    }
  >;
}

/** Create `CoarseIndicesAndSortings` for a coarse elements hierarchy */
function createCoarseIndicesAndSortings(coarseElements: CoarseElements): CoarseIndicesAndSortings | undefined {
  if (coarseElements.count === 0) return undefined;
  const { entity_id, asym_id, seq_id_begin, seq_id_end, chainElementSegments } = coarseElements;
  const { Present } = Column.ValueKind;
  const nChains = Math.max(chainElementSegments.count, 0); // chainElementSegments.count is -1 when there are no coarse elements

  const chainsByEntityId = new MultiMap<string, ChainIndex>();
  const chainsByAsymId = new MultiMap<string, ChainIndex>();
  const elementsSortedBySeqIdBegin = new Map<
    ChainIndex,
    Sorting<ElementIndex, number> & { endUpperBounds: SortedArray }
  >();
  const chains = {
    count: nChains,
    label_entity_id: new Array<string>(nChains),
    label_asym_id: new Array<string>(nChains),
  };

  for (let iChain = 0 as ChainIndex; iChain < nChains; iChain++) {
    const iElemFrom = chainElementSegments.offsets[iChain];
    const iElemTo = chainElementSegments.offsets[iChain + 1];
    const entityId = entity_id.value(iElemFrom);
    const asymId = asym_id.value(iElemFrom);
    chains.label_entity_id[iChain] = entityId;
    chains.label_asym_id[iChain] = asymId;
    chainsByEntityId.add(entityId, iChain);
    chainsByAsymId.add(asymId, iChain);

    const elementsWithSeqIds = filterInPlace(
      range(iElemFrom, iElemTo) as ElementIndex[],
      (iElem) => seq_id_begin.valueKind(iElem) === Present && seq_id_end.valueKind(iElem) === Present,
    );
    const sorting = Sorting.create(elementsWithSeqIds, seq_id_begin.value);
    const endBounds = sorting.keys.map(seq_id_end.value);
    // Ensure non-decreasing endBounds:
    for (let i = 1; i < endBounds.length; i++) {
      if (endBounds[i - 1] > endBounds[i]) {
        endBounds[i] = endBounds[i - 1];
      }
    }
    elementsSortedBySeqIdBegin.set(iChain, { ...sorting, endUpperBounds: SortedArray.ofSortedArray(endBounds) });
  }

  return { chains, chainsByEntityId, chainsByAsymId, elementsSortedBySeqIdBegin };
}

/** Represents a set of things (keys) of type `K`, sorted by some property (value) of type `V` */
export interface Sorting<K, V extends number> {
  /** Keys sorted by their corresponding values */
  keys: readonly K[];
  /** Sorted values corresponding to each key (value for `keys[i]` is `values[i]`) */
  values: SortedArray<V>;
}

export const Sorting = {
  /** Create a `Sorting` from an array of keys and a function returning their corresponding values.
   * If two keys have the same value, the smaller key will come first.
   * This function modifies `keys` - create a copy if you need the original order! */
  create<K extends number, V extends number>(keys: K[], valueFunction: (k: K) => V): Sorting<K, V> {
    sortIfNeeded(keys, (a, b) => valueFunction(a) - valueFunction(b) || a - b);
    const values: SortedArray<V> = SortedArray.ofSortedArray(keys.map(valueFunction));
    return { keys, values };
  },

  /** Return a newly allocated array of keys which have value equal to `target`.
   * The returned keys are sorted by their value. */
  getKeysWithValue<K, V extends number>(sorting: Sorting<K, V>, target: V): K[] {
    return Sorting.getKeysWithValueInRange(sorting, target, target);
  },

  /** Set `out.from`, `out.to` so that `keys.slice(out.from, out.to)` are the keys which have value within interval `[min, max]`
   * (inclusive; undefined `min` is negative infinity, undefined `max` is positive infinity). Does not allocate.
   * Bounds are compared as numbers (a numeric string works as its number, anything else selects nothing),
   * and need not be integers. */
  findIndexRange<K, V extends number>(
    sorting: Sorting<K, V>,
    min: V | undefined,
    max: V | undefined,
    out: { from: number; to: number },
  ) {
    const { keys, values } = sorting;
    const lo = min !== undefined ? Number(min) : -Infinity;
    const hi = max !== undefined ? Number(max) : Infinity;
    if (!keys || Number.isNaN(lo) || Number.isNaN(hi)) {
      out.from = 0;
      out.to = 0;
      return out;
    }
    const from = min !== undefined ? SortedArray.findPredecessorIndex(values, lo as V) : 0;
    // first index with a value greater than `hi`
    let a = from;
    let b = values.length;
    while (a < b) {
      const mid = (a + b) >> 1;
      if (values[mid] <= hi) a = mid + 1;
      else b = mid;
    }
    out.from = from;
    out.to = a;
    return out;
  },

  /** Return a newly allocated array of keys which have value within interval `[min, max]` (inclusive).
   * The returned keys are sorted by their value.
   * Undefined `min` is interpreted as negative infitity, undefined `max` is interpreted as positive infinity. */
  getKeysWithValueInRange<K, V extends number>(sorting: Sorting<K, V>, min: V | undefined, max: V | undefined): K[] {
    const r = Sorting.findIndexRange(sorting, min, max, { from: 0, to: 0 });
    return sorting.keys ? sorting.keys.slice(r.from, r.to) : [];
  },
};
