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
import { range, sortIfNeeded } from '@molstar/core/util/array';
import { ElementRanges, ElementRangesBuffer } from './element-ranges.js';
import { IntBuffer, MaxRetainedBufferCapacity, ScratchMemory } from '@molstar/core/util/int-buffer';
import { hasValue, hasWrongTypes } from './fields.js';
import { SchemaSelectionScratch } from './scratch.js';
import {
  type AtomicIndicesAndSortings,
  type CoarseIndicesAndSortings,
  IndicesAndSortings,
  Sorting,
} from './indexing.js';
import type { SchemaItem } from '../schema.js';

const EmptyArray: readonly any[] = [];

// ATOMIC SELECTIONS

/** Reusable buffers for `getAtomRangesForRowInto`, so that processing many rows does not allocate per row.
 * (JS is single-threaded and the functions are not reentrant, so one scratch can be shared by callers which never interleave.) */
export class AtomRangesScratch {
  /** Result of the last `getAtomRangesForRowInto` call which used this scratch (overwritten by the next call). */
  readonly ranges = new ElementRangesBuffer();
  readonly chains = new IntBuffer();
  readonly residues = new IntBuffer();

  /** Drop buffers which grew very large (so that one huge call does not pin memory forever). */
  trim(limit = MaxRetainedBufferCapacity) {
    this.ranges.trim(limit);
    this.chains.trim(limit);
    this.residues.trim(limit);
  }
}

const SharedScratch = new AtomRangesScratch();
const RowsBag = new ElementRangesBuffer();
ScratchMemory.register(SchemaSelectionScratch, () => {
  SharedScratch.trim(0);
  RowsBag.trim(0);
});
const IndexRange = { from: 0, to: 0 };

/** The symmetry operator of a unit, as far as rows can refer to it (`instance_id`, `operator_name`) */
export interface RowOperator {
  /** `SymmetryOperator.instanceId`, e.g. `ASM-X0-5` */
  readonly instanceId: string;
  /** `SymmetryOperator.name`, e.g. `1_555` */
  readonly name: string;
}

/** Does `row` apply to units with this operator? A row which names an instance or an operator applies only to it. */
export function rowAppliesToOperator(row: SchemaItem, operator: RowOperator): boolean {
  return (
    (!hasValue(row.instance_id) || row.instance_id === operator.instanceId) &&
    (!hasValue(row.operator_name) || row.operator_name === operator.name)
  );
}

/** `true` if some of `rows` refer to an instance or an operator (otherwise the selected elements do not depend on the unit's operator) */
export function rowsUseOperators(rows: readonly SchemaItem[]): boolean {
  return rows.some((row) => hasValue(row.instance_id) || hasValue(row.operator_name));
}
const SingleResidue: ResidueIndex[] = [0 as ResidueIndex];

/** Is the `auth_seq_id` range of `row` a range of residues in file order (it has an insertion code on a bound)? */
export function hasInsCodeRange(row: SchemaItem): boolean {
  return hasValue(row.beg_auth_ins_code) || hasValue(row.end_auth_ins_code);
}

/** An insertion code of a range bound without the `auth_seq_id` of the bound is not a bound of anything, the row selects nothing. */
function hasInsCodeWithoutSeqId(row: SchemaItem): boolean {
  return (
    (hasValue(row.beg_auth_ins_code) && !hasValue(row.beg_auth_seq_id)) ||
    (hasValue(row.end_auth_ins_code) && !hasValue(row.end_auth_seq_id))
  );
}

const ResolveRange = { from: 0, to: 0 };
/** The first (or last) residue of `iChain` with this `auth_seq_id` (and insertion code, unless there is no `insCode`) in file order, or -1. */
function findResidueByAuthId(
  model: Model,
  indices: AtomicIndicesAndSortings,
  iChain: ChainIndex,
  authSeqId: number,
  insCode: string | undefined | null,
  first: boolean,
): number {
  const sorting = indices.residuesSortedByAuthSeqId.get(iChain)!;
  Sorting.findIndexRange(sorting, authSeqId, authSeqId, ResolveRange);
  const { pdbx_PDB_ins_code } = model.atomicHierarchy.residues;
  const hasInsCode = hasValue(insCode);
  let found = -1;
  for (let k = ResolveRange.from; k < ResolveRange.to; k++) {
    const iRes = sorting.keys[k];
    if (hasInsCode && pdbx_PDB_ins_code.value(iRes) !== insCode) continue;
    if (found < 0 || (first ? iRes < found : iRes > found)) found = iRes;
  }
  return found;
}

/** Is the chain of waters (its entity is)? */
function isWaterChain(model: Model, iChain: ChainIndex): boolean {
  const { entities, atomicHierarchy } = model;
  const iEntity = entities.getEntityIndex(atomicHierarchy.chains.label_entity_id.value(iChain));
  return iEntity >= 0 && entities.data.type.value(iEntity) === 'water';
}

/** The first (or last) residue with this `auth_seq_id` (and insertion code, if there is one) in the order of the file in the chains, or -1. */
function findResidueInChains(
  model: Model,
  indices: AtomicIndicesAndSortings,
  chains: readonly ChainIndex[],
  authSeqId: number,
  insCode: string | undefined | null,
  first: boolean,
): number {
  let found = -1;
  for (const iChain of chains) {
    const iRes = findResidueByAuthId(model, indices, iChain, authSeqId, insCode, first);
    if (iRes >= 0 && (found < 0 || (first ? iRes < found : iRes > found))) found = iRes;
  }
  return found;
}

/** The residues of an auth chain (all the chains with an `auth_asym_id`, e.g. a polymer and its ligands which have label chains of their own) in the range of `row` with
 * insertion codes: from `from` to `to` (both inclusive) in the order of the file, or `null` if a bound is not in any of the chains or they are in reverse order.
 * A bound is looked for in the chains which are not waters first, and in the waters only if they are the only ones which have it (waters are numbered
 * again from 1 in many files, so their numbers are not those of the polymer which they would otherwise take the range to the end of). */
type InsCodeRange = { from: number; to: number } | null;

function resolveInsCodeRangeOfAuthChain(
  model: Model,
  row: SchemaItem,
  indices: AtomicIndicesAndSortings,
  authAsymId: string,
  iChain: ChainIndex,
): InsCodeRange {
  const { residueAtomSegments, chainAtomSegments } = model.atomicHierarchy;
  const chains: readonly ChainIndex[] = indices.chainsByAuthAsymId.get(authAsymId) ?? [iChain];
  const notWaters = chains.length > 1 ? chains.filter((c) => !isWaterChain(model, c)) : chains;
  const bound = (seqId: number, insCode: string | undefined | null, first: boolean) => {
    const iRes = notWaters.length > 0 ? findResidueInChains(model, indices, notWaters, seqId, insCode, first) : -1;
    return iRes >= 0 || notWaters.length === chains.length
      ? iRes
      : findResidueInChains(model, indices, chains, seqId, insCode, first);
  };
  // (without a bound, the range starts at the first residue of the auth chain and ends at its last)
  let from: number = Infinity;
  let to: number = -Infinity;
  for (const iC of chains) {
    from = Math.min(from, residueAtomSegments.index[chainAtomSegments.offsets[iC]]);
    to = Math.max(to, residueAtomSegments.index[chainAtomSegments.offsets[iC + 1] - 1]);
  }
  if (hasValue(row.beg_auth_seq_id)) from = bound(row.beg_auth_seq_id, row.beg_auth_ins_code, true);
  if (hasValue(row.end_auth_seq_id)) to = bound(row.end_auth_seq_id, row.end_auth_ins_code, false);
  if (from < 0 || to < 0 || from > to) return null;
  return { from, to };
}

/** The residues of `iChain` selected by the `auth_seq_id` range of `row` which has insertion codes: `out.from..out.to` (both inclusive). The range is
 * resolved in the auth chain of the chain (see `InsCodeRange`), and the chain gets the part of it which is its own. `false` if the chain has none.
 * `cache` keeps the range of each auth chain for the rest of the row. */
function resolveInsCodeRange(
  model: Model,
  row: SchemaItem,
  indices: AtomicIndicesAndSortings,
  iChain: ChainIndex,
  out: { from: number; to: number },
  cache?: Map<string, InsCodeRange>,
): boolean {
  const { residueAtomSegments, chainAtomSegments } = model.atomicHierarchy;
  const authAsymId = model.atomicHierarchy.chains.auth_asym_id.value(iChain);
  let range = cache?.get(authAsymId);
  if (range === undefined) {
    range = resolveInsCodeRangeOfAuthChain(model, row, indices, authAsymId, iChain);
    cache?.set(authAsymId, range);
  }
  if (!range) return false;
  const first = residueAtomSegments.index[chainAtomSegments.offsets[iChain]];
  const last = residueAtomSegments.index[chainAtomSegments.offsets[iChain + 1] - 1];
  out.from = Math.max(range.from, first);
  out.to = Math.min(range.to, last);
  return out.from <= out.to;
}

/** Return atom ranges in `model` which satisfy criteria given by `row` (a new object). */
export function getAtomRangesForRow(
  row: SchemaItem,
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
): ElementRanges | undefined {
  const r = getAtomRangesForRowInto(row, model, operator, indices, SharedScratch);
  const result = r?.toRanges();
  SharedScratch.trim();
  return result;
}

/** Like `getAtomRangesForRow`, but writes the result into `scratch.ranges` and returns it (no allocation per call).
 * The returned buffer is only valid until the next call with the same `scratch`; copy it if you need to keep it. */
export function getAtomRangesForRowInto(
  row: SchemaItem,
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
  scratch: AtomRangesScratch,
): ElementRangesBuffer | undefined {
  if (!indices.atomic) return undefined;
  if (!rowAppliesToOperator(row, operator)) return undefined;

  const atomicIndices = indices.atomic;
  const h = model.atomicHierarchy;
  const nAtoms = h.atoms._rowCount;
  if (nAtoms === 0) return undefined;

  const out = scratch.ranges.clear();
  // (values are matched as they are, so a row with a value of the wrong type cannot select anything; neither can an insertion code without its number)
  if (hasWrongTypes(row) || hasInsCodeWithoutSeqId(row)) return out;

  const hasAtomIds = hasValue(row.atom_id) || hasValue(row.atom_index);
  const hasLabelAtomId = hasValue(row.label_atom_id);
  const hasAuthAtomId = hasValue(row.auth_atom_id);
  const hasTypeSymbol = hasValue(row.type_symbol);
  const hasLabelCompId = hasValue(row.label_comp_id);
  const hasAuthCompId = hasValue(row.auth_comp_id);
  const labelCompIdPerAtom = hasLabelCompId && !atomicIndices.residuesByLabelCompIdIsPure;
  const authCompIdPerAtom = hasAuthCompId && !atomicIndices.residuesByAuthCompIdIsPure;
  const hasAtomFilter = hasLabelAtomId || hasAuthAtomId || hasTypeSymbol || labelCompIdPerAtom || authCompIdPerAtom;
  const hasResidueFilter =
    hasValue(row.label_seq_id) ||
    hasValue(row.auth_seq_id) ||
    hasValue(row.pdbx_PDB_ins_code) ||
    hasValue(row.beg_label_seq_id) ||
    hasValue(row.end_label_seq_id) ||
    hasValue(row.beg_auth_seq_id) ||
    hasValue(row.end_auth_seq_id) ||
    hasLabelCompId ||
    hasAuthCompId ||
    hasValue(row.residue_index);
  const hasChainFilter = hasValue(row.label_asym_id) || hasValue(row.auth_asym_id) || hasValue(row.label_entity_id);

  if (hasAtomIds) {
    addAtomsForAtomIds(out, model, row, atomicIndices);
    return out;
  }

  if (!hasChainFilter && !hasResidueFilter && !hasAtomFilter) {
    return out.add(0, nAtoms);
  }

  const qualifyingChains = getQualifyingChains(model, row, atomicIndices, scratch.chains);
  if (!hasResidueFilter && !hasAtomFilter) {
    const chainOffsets = h.chainAtomSegments.offsets;
    // (`qualifyingChains` is only `undefined` (= all chains) without any chain filter, which was handled above)
    for (let i = 0, n = qualifyingChains!.length; i < n; i++) {
      const iChain = qualifyingChains![i];
      out.add(chainOffsets[iChain], chainOffsets[iChain + 1]);
    }
    return out;
  }

  const qualifyingResidues = getQualifyingResidues(model, row, atomicIndices, qualifyingChains, scratch.residues);
  const residues = qualifyingResidues.data;
  const nResidues = qualifyingResidues.n;
  const residueOffsets = h.residueAtomSegments.offsets;
  if (!hasAtomFilter) {
    for (let i = 0; i < nResidues; i++) {
      const iRes = residues[i];
      out.add(residueOffsets[iRes], residueOffsets[iRes + 1]);
    }
    return out;
  }

  const { label_atom_id, auth_atom_id, type_symbol, label_comp_id, auth_comp_id } = h.atoms;
  const typeSymbolUpper = hasTypeSymbol ? row.type_symbol!.toUpperCase() : undefined;
  for (let i = 0; i < nResidues; i++) {
    const iRes = residues[i];
    for (let iAtom = residueOffsets[iRes], stop = residueOffsets[iRes + 1]; iAtom < stop; iAtom++) {
      if (hasLabelAtomId && label_atom_id.value(iAtom) !== row.label_atom_id) continue;
      if (hasAuthAtomId && auth_atom_id.value(iAtom) !== row.auth_atom_id) continue;
      if (hasTypeSymbol && type_symbol.value(iAtom) !== typeSymbolUpper) continue;
      if (labelCompIdPerAtom && label_comp_id.value(iAtom) !== row.label_comp_id) continue;
      if (authCompIdPerAtom && auth_comp_id.value(iAtom) !== row.auth_comp_id) continue;
      out.add(iAtom, iAtom + 1);
    }
  }
  return out;
}

/** Return atom ranges in `model` which satisfy criteria given by any of `rows` (atoms that satisfy more rows are still included only once) */
export function getAtomRangesForRows(
  rows: readonly SchemaItem[],
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
): ElementRanges {
  // Per-row results go to a shared buffer (`SharedScratch`), are collected in a shared bag (`RowsBag`) and merged once;
  // only the final result is a new object (callers may keep it).
  const bag = RowsBag.clear();
  for (const row of rows) {
    const r = getAtomRangesForRowInto(row, model, operator, indices, SharedScratch);
    if (r) bag.append(r);
  }
  const result = ElementRanges.normalize(bag);
  SharedScratch.trim();
  RowsBag.trim(MaxRetainedBufferCapacity);
  return result;
}

/** Keep elements of `src` for which `column.value(element) === value`, writing them to `out`
 * (`src` may be a view of `out` itself, then it is filtered in place). Return a view of the result. */
function filterByColumn(src: ArrayLike<number>, out: IntBuffer, column: Column<string>, value: string): Int32Array {
  out.ensureCapacity(src.length);
  const data = out.data;
  let n = 0;
  for (let i = 0; i < src.length; i++) {
    const x = src[i];
    if (column.value(x) === value) data[n++] = x;
  }
  out.n = n;
  return out.view();
}

/** Return chain indexes which satisfy criteria given by `row` (ascending order), or `undefined` if the row has no chain criteria (i.e. all chains qualify).
 * The returned array can be one of the shared indices or a view of `out` (do not modify it, it is valid until `out` is modified). */
function getQualifyingChains(
  model: Model,
  row: SchemaItem,
  indices: AtomicIndicesAndSortings,
  out: IntBuffer,
): ArrayLike<number> | undefined {
  const { auth_asym_id, label_entity_id } = model.atomicHierarchy.chains;
  let result: ArrayLike<number> | undefined = undefined;
  if (hasValue(row.label_asym_id)) {
    result = indices.chainsByLabelAsymId.get(row.label_asym_id) ?? EmptyArray;
  }
  if (hasValue(row.auth_asym_id)) {
    if (result) {
      result = filterByColumn(result, out, auth_asym_id, row.auth_asym_id);
    } else {
      result = indices.chainsByAuthAsymId.get(row.auth_asym_id) ?? EmptyArray;
    }
  }
  if (hasValue(row.label_entity_id)) {
    if (result) {
      result = filterByColumn(result, out, label_entity_id, row.label_entity_id);
    } else {
      result = indices.chainsByLabelEntityId.get(row.label_entity_id) ?? EmptyArray;
    }
  }
  return result;
}

/** Which residue index is used to get the candidate residues in `getQualifyingResidues` (the other criteria are then checked on each candidate). */
const enum ResidueSource {
  AllInChain,
  LabelSeqId,
  AuthSeqId,
  ResidueIndex,
  InsCode,
  LabelSeqIdRange,
  AuthSeqIdRange,
  /** The residues of a range of `auth_seq_id` with insertion codes (in file order): consecutive residues */
  InsCodeRange,
  LabelCompId,
  AuthCompId,
}

const InsRange = { from: 0, to: 0 };
const InsGroups = new Map<string, InsCodeRange>();

/** Write residue indexes (ascending) which satisfy criteria given by `row` into `out` (clearing it first) and return it.
 * `fromChains` are ascending chain indexes to search in, `undefined` means all chains.
 *
 * For each chain, the candidates are taken from the index of the first criterion present (so the work is proportional
 * to the number of candidates, not the size of the chain) and all remaining criteria are checked in one pass,
 * without intermediate arrays. */
function getQualifyingResidues(
  model: Model,
  row: SchemaItem,
  indices: AtomicIndicesAndSortings,
  fromChains: ArrayLike<number> | undefined,
  out: IntBuffer,
): IntBuffer {
  out.clear();
  const { label_seq_id, auth_seq_id, pdbx_PDB_ins_code } = model.atomicHierarchy.residues;
  const { label_comp_id, auth_comp_id } = model.atomicHierarchy.atoms;
  const { residueAtomSegments, chainAtomSegments, residueSourceIndex } = model.atomicHierarchy;
  const resOffsets = residueAtomSegments.offsets;
  const { Present } = Column.ValueKind;

  const hasLabelSeq = hasValue(row.label_seq_id);
  const hasAuthSeq = hasValue(row.auth_seq_id);
  const hasResIdx = hasValue(row.residue_index);
  const hasIns = hasValue(row.pdbx_PDB_ins_code);
  // An empty insertion code means "no insertion code" (as in MolScript, where it is what residues without one have). It is a
  // check on candidates and never a lookup: the indices only list residues which have an insertion code.
  const hasInsLookup = hasIns && row.pdbx_PDB_ins_code !== '';
  const hasBegLabel = hasValue(row.beg_label_seq_id);
  const hasEndLabel = hasValue(row.end_label_seq_id);
  const hasBegAuth = hasValue(row.beg_auth_seq_id);
  const hasEndAuth = hasValue(row.end_auth_seq_id);
  const hasLabelRange = hasBegLabel || hasEndLabel;
  const hasInsRange = hasInsCodeRange(row);
  if (hasInsRange) InsGroups.clear();
  // (a range with insertion codes is a range of residues in file order, not of `auth_seq_id` values)
  const hasAuthRange = (hasBegAuth || hasEndAuth) && !hasInsRange;
  const hasLabelComp = hasValue(row.label_comp_id);
  const hasAuthComp = hasValue(row.auth_comp_id);
  const labelCompPure = indices.residuesByLabelCompIdIsPure;
  const authCompPure = indices.residuesByAuthCompIdIsPure;

  const source = hasLabelSeq
    ? ResidueSource.LabelSeqId
    : hasAuthSeq
      ? ResidueSource.AuthSeqId
      : hasResIdx
        ? ResidueSource.ResidueIndex
        : hasInsLookup
          ? ResidueSource.InsCode
          : hasLabelRange
            ? ResidueSource.LabelSeqIdRange
            : hasAuthRange
              ? ResidueSource.AuthSeqIdRange
              : hasInsRange
                ? ResidueSource.InsCodeRange
                : hasLabelComp
                  ? ResidueSource.LabelCompId
                  : hasAuthComp
                    ? ResidueSource.AuthCompId
                    : ResidueSource.AllInChain;

  // Without chain criteria the candidates do not have to be looked up in every chain: some criteria select few residues in the whole
  // model, so a lookup over all chains at once is better (`covered`: the criterion which the global lookup satisfies).
  // (Not for sequence-based criteria: they select about one residue in each chain anyway, which is what looking in each chain finds.)
  let globalKeys: readonly ResidueIndex[] | undefined = undefined;
  let covered = source;
  // (not for a range with insertion codes: it is resolved in each chain)
  if (!fromChains && !hasInsRange) {
    const global = indices.global;
    const bySeq = hasLabelSeq || hasAuthSeq || hasLabelRange || hasAuthRange;
    if (hasResIdx && global.residueBySourceIndex) {
      const iRes = global.residueBySourceIndex.get(row.residue_index!);
      SingleResidue[0] = iRes as ResidueIndex;
      globalKeys = iRes === undefined ? EmptyArray : SingleResidue;
      covered = ResidueSource.ResidueIndex;
    } else if (hasInsLookup && !bySeq) {
      globalKeys = global.residuesByInsCode.get(row.pdbx_PDB_ins_code!) ?? EmptyArray;
      covered = ResidueSource.InsCode;
    } else if (hasLabelComp && !bySeq) {
      globalKeys = global.residuesByLabelCompId.get(row.label_comp_id!) ?? EmptyArray;
      covered = ResidueSource.LabelCompId;
    } else if (hasAuthComp && !bySeq) {
      globalKeys = global.residuesByAuthCompId.get(row.auth_comp_id!) ?? EmptyArray;
      covered = ResidueSource.AuthCompId;
    }
  }

  // criteria which still have to be checked on candidates (the one used as the source of candidates is already satisfied by its index)
  const checkLabelSeq = hasLabelSeq && covered !== ResidueSource.LabelSeqId;
  const checkAuthSeq = hasAuthSeq && covered !== ResidueSource.AuthSeqId;
  const checkResIdx = hasResIdx && covered !== ResidueSource.ResidueIndex;
  const checkIns = hasIns && covered !== ResidueSource.InsCode;
  const checkLabelRange = hasLabelRange && covered !== ResidueSource.LabelSeqIdRange;
  const checkAuthRange = hasAuthRange && covered !== ResidueSource.AuthSeqIdRange;
  const checkInsRange = hasInsRange && covered !== ResidueSource.InsCodeRange;
  const checkLabelComp = hasLabelComp && covered !== ResidueSource.LabelCompId;
  const checkAuthComp = hasAuthComp && covered !== ResidueSource.AuthCompId;

  // one pass over the candidates of a chain, of the whole model, or of a global lookup
  const nPasses =
    globalKeys || (!fromChains && source === ResidueSource.AllInChain)
      ? 1
      : fromChains
        ? fromChains.length
        : model.atomicHierarchy.chains._rowCount;
  for (let iFrom = 0; iFrom < nPasses; iFrom++) {
    const iChain = (fromChains ? fromChains[iFrom] : iFrom) as ChainIndex;
    if (hasInsRange && !resolveInsCodeRange(model, row, indices, iChain, InsRange, InsGroups)) continue;

    // candidates are `keys[from..to)` (or consecutive residues `from..to` when `keys` is undefined)
    let keys: readonly ResidueIndex[] | undefined = undefined;
    let from = 0;
    let to = 0;
    if (globalKeys) {
      keys = globalKeys;
      to = globalKeys.length;
    } else if (!fromChains && source === ResidueSource.AllInChain) {
      to = model.atomicHierarchy.residues._rowCount; // (all residues in all chains)
    } else
      switch (source) {
        case ResidueSource.LabelSeqId:
        case ResidueSource.LabelSeqIdRange: {
          const sorting = indices.residuesSortedByLabelSeqId.get(iChain)!;
          keys = sorting.keys;
          if (source === ResidueSource.LabelSeqId)
            Sorting.findIndexRange(sorting, row.label_seq_id, row.label_seq_id, IndexRange);
          else
            Sorting.findIndexRange(
              sorting,
              hasBegLabel ? row.beg_label_seq_id : undefined,
              hasEndLabel ? row.end_label_seq_id : undefined,
              IndexRange,
            );
          from = IndexRange.from;
          to = IndexRange.to;
          break;
        }
        case ResidueSource.AuthSeqId:
        case ResidueSource.AuthSeqIdRange: {
          const sorting = indices.residuesSortedByAuthSeqId.get(iChain)!;
          keys = sorting.keys;
          if (source === ResidueSource.AuthSeqId)
            Sorting.findIndexRange(sorting, row.auth_seq_id, row.auth_seq_id, IndexRange);
          else
            Sorting.findIndexRange(
              sorting,
              hasBegAuth ? row.beg_auth_seq_id : undefined,
              hasEndAuth ? row.end_auth_seq_id : undefined,
              IndexRange,
            );
          from = IndexRange.from;
          to = IndexRange.to;
          break;
        }
        case ResidueSource.ResidueIndex: {
          const sorting = indices.residuesSortedBySourceIndex.get(iChain)!;
          keys = sorting.keys;
          Sorting.findIndexRange(sorting, row.residue_index, row.residue_index, IndexRange);
          from = IndexRange.from;
          to = IndexRange.to;
          break;
        }
        case ResidueSource.InsCode:
          keys = indices.residuesByInsCode.get(iChain)!.get(row.pdbx_PDB_ins_code!) ?? EmptyArray;
          to = keys.length;
          break;
        case ResidueSource.LabelCompId:
          keys = indices.residuesByLabelCompId.get(iChain)!.get(row.label_comp_id!) ?? EmptyArray;
          to = keys.length;
          break;
        case ResidueSource.AuthCompId:
          keys = indices.residuesByAuthCompId.get(iChain)!.get(row.auth_comp_id!) ?? EmptyArray;
          to = keys.length;
          break;
        case ResidueSource.InsCodeRange:
          from = InsRange.from;
          to = InsRange.to + 1;
          break;
        case ResidueSource.AllInChain:
          from = residueAtomSegments.index[chainAtomSegments.offsets[iChain]];
          to = residueAtomSegments.index[chainAtomSegments.offsets[iChain + 1] - 1] + 1;
          break;
      }

    const start = out.n;
    for (let k = from; k < to; k++) {
      const iRes = (keys ? keys[k] : k) as ResidueIndex;
      if (checkLabelSeq && !(label_seq_id.valueKind(iRes) === Present && label_seq_id.value(iRes) === row.label_seq_id))
        continue;
      if (checkAuthSeq && !(auth_seq_id.valueKind(iRes) === Present && auth_seq_id.value(iRes) === row.auth_seq_id))
        continue;
      if (checkInsRange && (iRes < InsRange.from || iRes > InsRange.to)) continue;
      if (checkResIdx && residueSourceIndex.value(iRes) !== row.residue_index) continue;
      if (checkIns && pdbx_PDB_ins_code.value(iRes) !== row.pdbx_PDB_ins_code) continue;
      if (checkLabelRange) {
        if (
          hasBegLabel &&
          !(label_seq_id.valueKind(iRes) === Present && label_seq_id.value(iRes) >= row.beg_label_seq_id!)
        )
          continue;
        if (
          hasEndLabel &&
          !(label_seq_id.valueKind(iRes) === Present && label_seq_id.value(iRes) <= row.end_label_seq_id!)
        )
          continue;
      }
      if (checkAuthRange) {
        if (hasBegAuth && !(auth_seq_id.valueKind(iRes) === Present && auth_seq_id.value(iRes) >= row.beg_auth_seq_id!))
          continue;
        if (hasEndAuth && !(auth_seq_id.valueKind(iRes) === Present && auth_seq_id.value(iRes) <= row.end_auth_seq_id!))
          continue;
      }
      if (checkLabelComp && !residueHasCompId(label_comp_id, resOffsets, iRes, row.label_comp_id!, labelCompPure))
        continue;
      if (checkAuthComp && !residueHasCompId(auth_comp_id, resOffsets, iRes, row.auth_comp_id!, authCompPure)) continue;
      out.push(iRes);
    }
    // Candidates from a value range are ordered by value, not necessarily by residue index. (Candidates for a single value, from
    // the other indices and from consecutive residues are already ascending; chains are ascending, so the whole `out` is then too.)
    if (!globalKeys && (source === ResidueSource.LabelSeqIdRange || source === ResidueSource.AuthSeqIdRange)) {
      out.data.subarray(start, out.n).sort(); // (native numeric sort of the typed array)
    }
  }
  return out;
}

/** Does residue `iRes` have `compId` (as the comp_id of its first atom if `pure`, otherwise of any of its atoms)? */
function residueHasCompId(
  compIdColumn: Column<string>,
  residueOffsets: ArrayLike<number>,
  iRes: ResidueIndex,
  compId: string,
  pure: boolean,
): boolean {
  if (pure) return compIdColumn.value(residueOffsets[iRes]) === compId;
  for (let iAtom = residueOffsets[iRes], stop = residueOffsets[iRes + 1]; iAtom < stop; iAtom++) {
    if (compIdColumn.value(iAtom) === compId) return true;
  }
  return false;
}

/** Add (ascending) the atoms in `model` which satisfy criteria given by `row` to `out`.
 * Only works when `row.atom_id` and/or `row.atom_index` is defined (otherwise use `getAtomRangesForRowInto`).
 * `atom_index` is unique, `atom_id` is not necessarily (all atoms with the id are candidates). */
function addAtomsForAtomIds(
  out: ElementRangesBuffer,
  model: Model,
  row: SchemaItem,
  indices: AtomicIndicesAndSortings,
) {
  if (!hasValue(row.atom_id) && !hasValue(row.atom_index))
    throw new Error('ArgumentError: at least one of row.atom_id, row.atom_index must be defined.');
  if (hasInsCodeRange(row)) InsGroups.clear(); // (the range of each auth chain is resolved once for the atoms of the row)
  if (hasValue(row.atom_index)) {
    // (`atomQualifies` also checks `atom_id`, if given)
    const iAtom = indices.atoms.bySourceIndex.get(row.atom_index);
    if (iAtom !== undefined && atomQualifiesInInsCodeRange(model, iAtom, row, indices)) out.add(iAtom, iAtom + 1);
    return;
  }
  const id = row.atom_id!;
  const first = indices.atoms.byId.get(id);
  if (first === undefined) return;
  const candidates = indices.atoms.repeatedIds.get(id) ?? [first];
  for (let i = 0; i < candidates.length; i++) {
    const iAtom = candidates[i];
    if (atomQualifiesInInsCodeRange(model, iAtom, row, indices)) out.add(iAtom, iAtom + 1);
  }
}

/** `atomQualifies`, and the residue is in the range with insertion codes of `row` (if it has one) */
function atomQualifiesInInsCodeRange(
  model: Model,
  iAtom: ElementIndex,
  row: SchemaItem,
  indices: AtomicIndicesAndSortings,
): boolean {
  if (!atomQualifies(model, iAtom, row)) return false;
  if (!hasInsCodeRange(row)) return true;
  const h = model.atomicHierarchy;
  const iRes = h.residueAtomSegments.index[iAtom];
  return (
    resolveInsCodeRange(model, row, indices, h.chainAtomSegments.index[iAtom], InsRange, InsGroups) &&
    InsRange.from <= iRes &&
    iRes <= InsRange.to
  );
}

/** Return true if `iAtom`-th atom in `model` satisfies all selection criteria given by `row`
 * (a range with insertion codes, `beg_auth_ins_code` and `end_auth_ins_code`, needs the indices of the model and is not checked here). */
export function atomQualifies(model: Model, iAtom: ElementIndex, row: SchemaItem): boolean {
  const h = model.atomicHierarchy;

  const iChain = h.chainAtomSegments.index[iAtom];
  const label_asym_id = h.chains.label_asym_id.value(iChain);
  const auth_asym_id = h.chains.auth_asym_id.value(iChain);
  const label_entity_id = h.chains.label_entity_id.value(iChain);
  if (!matches(row.label_asym_id, label_asym_id)) return false;
  if (!matches(row.auth_asym_id, auth_asym_id)) return false;
  if (!matches(row.label_entity_id, label_entity_id)) return false;

  const iRes = h.residueAtomSegments.index[iAtom];
  const label_seq_id =
    h.residues.label_seq_id.valueKind(iRes) === Column.ValueKind.Present
      ? h.residues.label_seq_id.value(iRes)
      : undefined;
  const auth_seq_id =
    h.residues.auth_seq_id.valueKind(iRes) === Column.ValueKind.Present
      ? h.residues.auth_seq_id.value(iRes)
      : undefined;
  const pdbx_PDB_ins_code = h.residues.pdbx_PDB_ins_code.value(iRes);
  const residue_index = h.residueSourceIndex.value(iRes);
  if (!matches(row.label_seq_id, label_seq_id)) return false;
  if (!matches(row.auth_seq_id, auth_seq_id)) return false;
  if (!matches(row.pdbx_PDB_ins_code, pdbx_PDB_ins_code)) return false;
  if (!matchesRange(row.beg_label_seq_id, row.end_label_seq_id, label_seq_id)) return false;
  if (!hasInsCodeRange(row) && !matchesRange(row.beg_auth_seq_id, row.end_auth_seq_id, auth_seq_id)) return false;
  if (!matches(row.residue_index, residue_index)) return false;

  const label_comp_id = h.atoms.label_comp_id.value(iAtom);
  const auth_comp_id = h.atoms.auth_comp_id.value(iAtom);
  const label_atom_id = h.atoms.label_atom_id.value(iAtom);
  const auth_atom_id = h.atoms.auth_atom_id.value(iAtom);
  const type_symbol = h.atoms.type_symbol.value(iAtom);
  const atom_id = model.atomicConformation.atomId.value(iAtom);
  const atom_index = h.atomSourceIndex.value(iAtom);
  if (!matches(row.label_comp_id, label_comp_id)) return false;
  if (!matches(row.auth_comp_id, auth_comp_id)) return false;
  if (!matches(row.label_atom_id, label_atom_id)) return false;
  if (!matches(row.auth_atom_id, auth_atom_id)) return false;
  if (!matches(row.type_symbol?.toUpperCase(), type_symbol)) return false;
  if (!matches(row.atom_id, atom_id)) return false;
  if (!matches(row.atom_index, atom_index)) return false;

  return true;
}

/** Return true if `value` equals `requiredValue` or if `requiredValue` if not defined.  */
function matches<T>(requiredValue: T | undefined | null, value: T | undefined): boolean {
  return !hasValue(requiredValue) || value === requiredValue;
}

/** Return true if `requiredMin <= value <= requiredMax`.
 * Undefined `requiredMin` behaves like negative infinity.
 * Undefined `requiredMax` behaves like positive infinity. */
function matchesRange<T>(
  requiredMin: T | undefined | null,
  requiredMax: T | undefined | null,
  value: T | undefined,
): boolean {
  if (hasValue(requiredMin) && (!hasValue(value) || value < requiredMin)) return false;
  if (hasValue(requiredMax) && (!hasValue(value) || value > requiredMax)) return false;
  return true;
}

// COARSE SELECTIONS

/** Return sphere ranges in `model` which satisfy criteria given by `row` */
export function getSphereRangesForRow(
  row: SchemaItem,
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
): ElementRanges | undefined {
  if (!indices.spheres) return undefined;
  if (!rowAppliesToOperator(row, operator)) return undefined;
  return getCoarseElementRangesForRow(row, model.coarseHierarchy.spheres, indices.spheres);
}

/** Return sphere ranges in `model` which satisfy criteria given by any of `rows` (spheres that satisfy more rows are still included only once) */
export function getSphereRangesForRows(
  rows: readonly SchemaItem[],
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
): ElementRanges {
  return ElementRanges.union(rows.map((row) => getSphereRangesForRow(row, model, operator, indices)));
}

/** Return gaussian ranges in `model` which satisfy criteria given by `row` */
export function getGaussianRangesForRow(
  row: SchemaItem,
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
): ElementRanges | undefined {
  if (!indices.gaussians) return undefined;
  if (!rowAppliesToOperator(row, operator)) return undefined;
  return getCoarseElementRangesForRow(row, model.coarseHierarchy.gaussians, indices.gaussians);
}

/** Return gaussian ranges in `model` which satisfy criteria given by any of `rows` (gaussians that satisfy more rows are still included only once) */
export function getGaussianRangesForRows(
  rows: readonly SchemaItem[],
  model: Model,
  operator: RowOperator,
  indices: IndicesAndSortings,
): ElementRanges {
  return ElementRanges.union(rows.map((row) => getGaussianRangesForRow(row, model, operator, indices)));
}

/** Return ranges of coarse elements (spheres or gaussians) which satisfy criteria given by `row` */
export function getCoarseElementRangesForRow(
  row: SchemaItem,
  coarseElements: CoarseElements,
  indices: CoarseIndicesAndSortings,
): ElementRanges | undefined {
  const nElements = coarseElements.count;
  // (an element contains residues, which have integer numbers)
  if (nElements === 0 || hasWrongTypes(row) || (hasValue(row.label_seq_id) && !Number.isInteger(row.label_seq_id)))
    return undefined;

  const hasResidueFilter =
    hasValue(row.label_seq_id) || hasValue(row.beg_label_seq_id) || hasValue(row.end_label_seq_id);
  const hasChainFilter = hasValue(row.label_asym_id) || hasValue(row.label_entity_id);
  const hasInvalidFilter =
    hasValue(row.auth_asym_id) ||
    hasValue(row.auth_seq_id) ||
    hasValue(row.pdbx_PDB_ins_code) ||
    hasValue(row.beg_auth_seq_id) ||
    hasValue(row.end_auth_seq_id) ||
    hasValue(row.beg_auth_ins_code) ||
    hasValue(row.end_auth_ins_code) ||
    hasValue(row.label_comp_id) ||
    hasValue(row.auth_comp_id) ||
    hasValue(row.residue_index) ||
    hasValue(row.label_atom_id) ||
    hasValue(row.auth_atom_id) ||
    hasValue(row.type_symbol) ||
    hasValue(row.atom_id) ||
    hasValue(row.atom_index);

  if (hasInvalidFilter) {
    printCoarseSelectorWarning();
    return undefined;
  }

  // A range which ends before it begins is empty (as it is for atoms). MolScript would select the elements which happen to span the
  // two bounds there (or by some other rule of its overlap test for a combination with `label_seq_id`), which nobody can mean.
  if (hasValue(row.beg_label_seq_id) && hasValue(row.end_label_seq_id) && row.beg_label_seq_id > row.end_label_seq_id) {
    return undefined;
  }

  if (!hasChainFilter && !hasResidueFilter) {
    return ElementRanges.single(0 as ElementIndex, nElements as ElementIndex);
  }

  const qualifyingChains = getQualifyingCoarseChains(coarseElements, row, indices);
  if (!hasResidueFilter) {
    const chainOffsets = coarseElements.chainElementSegments.offsets;
    const ranges = ElementRanges.empty();
    for (const iChain of qualifyingChains) {
      ElementRanges.add(ranges, chainOffsets[iChain], chainOffsets[iChain + 1]);
    }
    return ranges;
  }

  const qualifyingElements = getQualifyingCoarseElements(coarseElements, row, indices, qualifyingChains);
  const ranges = ElementRanges.empty();
  for (const iElem of qualifyingElements) {
    ElementRanges.add(ranges, iElem, (iElem + 1) as ElementIndex);
  }
  return ranges;
}

/** Return an array of chain indexes which satisfy criteria given by `row` */
function getQualifyingCoarseChains(
  coarseElements: CoarseElements,
  row: SchemaItem,
  indices: CoarseIndicesAndSortings,
): readonly ChainIndex[] {
  let result: readonly ChainIndex[] | undefined = undefined;
  if (hasValue(row.label_asym_id)) {
    result = indices.chainsByAsymId.get(row.label_asym_id) ?? EmptyArray;
  }
  if (hasValue(row.label_entity_id)) {
    if (result) {
      result = result.filter(
        (iChain) =>
          coarseElements.entity_id.value(coarseElements.chainElementSegments.offsets[iChain]) === row.label_entity_id,
      );
    } else {
      result = indices.chainsByEntityId.get(row.label_entity_id) ?? EmptyArray;
    }
  }
  result ??= range(coarseElements.chainElementSegments.count) as ChainIndex[];
  return result;
}

/** Return an array of residue indexes which satisfy criteria given by `row` */
function getQualifyingCoarseElements(
  coarseElements: CoarseElements,
  row: SchemaItem,
  indices: CoarseIndicesAndSortings,
  fromChains: readonly ChainIndex[],
): ElementIndex[] {
  const result: ElementIndex[] = [];
  for (const iChain of fromChains) {
    const sorting = indices.elementsSortedBySeqIdBegin.get(iChain)!;
    const queryStart = Math.max(
      hasValue(row.label_seq_id) ? row.label_seq_id : -Infinity,
      hasValue(row.beg_label_seq_id) ? row.beg_label_seq_id : -Infinity,
    );
    const queryEnd = Math.min(
      hasValue(row.label_seq_id) ? row.label_seq_id : Infinity,
      hasValue(row.end_label_seq_id) ? row.end_label_seq_id : Infinity,
    ); // inclusive
    const iStart = SortedArray.findPredecessorIndex(sorting.endUpperBounds, queryStart); // select elements potentially ending >=queryStart (necessary condition)
    const iStop = SortedArray.findPredecessorIndex(sorting.values, queryEnd + 1); // select elements starting <=queryEnd (necessary and suffient condition) // exclusive

    for (let i = iStart; i < iStop; i++) {
      const iElem = sorting.keys[i];
      if (coarseElements.seq_id_end.value(iElem) >= queryStart) {
        // rechecking seq_id_end, as the condition was not sufficient
        result.push(iElem);
      }
    }
    // `queryStart > queryEnd` is possible when `label_seq_id` is outside of a valid range [beg_label_seq_id, end_label_seq_id]: elements which contain
    // `label_seq_id` and overlap the range qualify then, as with MolScript (an element is coarse: it can contain one and overlap the other).
  }
  sortIfNeeded(result, (a, b) => a - b);
  return result;
}

let coarseSelectorWarningPrinted = false;
function printCoarseSelectorWarning() {
  if (!coarseSelectorWarningPrinted) {
    console.warn(
      'Using unsupported selector fields (auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, beg_auth_seq_id, end_auth_seq_id, beg_auth_ins_code, end_auth_ins_code, label_comp_id, auth_comp_id, residue_index, label_atom_id, auth_atom_id, type_symbol, atom_id, atom_index) on a coarse structure. The resulting selection will be empty.',
    );
    coarseSelectorWarningPrinted = true;
  }
}
