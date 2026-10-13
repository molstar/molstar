/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '@molstar/core/data/db';
import { SortedArray } from '@molstar/core/data/int';
import {
  type ChainIndex,
  type ElementIndex,
  Model,
  type ResidueIndex,
  StructureElement,
} from '@molstar/model/model/structure';
import type { CoarseElements } from '@molstar/model/model/structure/model/properties/coarse';
import { Expression } from '@molstar/query-language/language/expression';
import { range, sortIfNeeded } from '@molstar/core/util/array';
import { ElementRanges, ElementRangesBuffer } from './element-ranges.js';
import { IntBuffer } from './int-buffer.js';
import {
  type AtomicIndicesAndSortings,
  type CoarseIndicesAndSortings,
  IndicesAndSortings,
  Sorting,
} from './indexing.js';
import type { MVSAnnotationRow } from './schemas.js';
import { isAnyDefined, isDefined } from './utils.js';

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

const MaxRetainedBufferCapacity = 1 << 20;
const SharedScratch = new AtomRangesScratch();
const RowsBag = new ElementRangesBuffer();
const IndexRange = { from: 0, to: 0 };
const SingleResidue: ResidueIndex[] = [0 as ResidueIndex];

/** Return atom ranges in `model` which satisfy criteria given by `row` (a new object). */
export function getAtomRangesForRow(
  row: MVSAnnotationRow,
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
): ElementRanges | undefined {
  const r = getAtomRangesForRowInto(row, model, instanceId, indices, SharedScratch);
  return r?.toRanges();
}

/** Like `getAtomRangesForRow`, but writes the result into `scratch.ranges` and returns it (no allocation per call).
 * The returned buffer is only valid until the next call with the same `scratch`; copy it if you need to keep it. */
export function getAtomRangesForRowInto(
  row: MVSAnnotationRow,
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
  scratch: AtomRangesScratch,
): ElementRangesBuffer | undefined {
  if (!indices.atomic) return undefined;
  if (isDefined(row.instance_id) && row.instance_id !== instanceId) return undefined;

  const atomicIndices = indices.atomic;
  const h = model.atomicHierarchy;
  const nAtoms = h.atoms._rowCount;
  if (nAtoms === 0) return undefined;

  const out = scratch.ranges.clear();

  const hasAtomIds = isDefined(row.atom_id) || isDefined(row.atom_index);
  const hasLabelAtomId = isDefined(row.label_atom_id);
  const hasAuthAtomId = isDefined(row.auth_atom_id);
  const hasTypeSymbol = isDefined(row.type_symbol);
  const hasLabelCompId = isDefined(row.label_comp_id);
  const hasAuthCompId = isDefined(row.auth_comp_id);
  const labelCompIdPerAtom = hasLabelCompId && !atomicIndices.residuesByLabelCompIdIsPure;
  const authCompIdPerAtom = hasAuthCompId && !atomicIndices.residuesByAuthCompIdIsPure;
  const hasAtomFilter = hasLabelAtomId || hasAuthAtomId || hasTypeSymbol || labelCompIdPerAtom || authCompIdPerAtom;
  const hasResidueFilter =
    isDefined(row.label_seq_id) ||
    isDefined(row.auth_seq_id) ||
    isDefined(row.pdbx_PDB_ins_code) ||
    isDefined(row.beg_label_seq_id) ||
    isDefined(row.end_label_seq_id) ||
    isDefined(row.beg_auth_seq_id) ||
    isDefined(row.end_auth_seq_id) ||
    hasLabelCompId ||
    hasAuthCompId ||
    isDefined(row.residue_index);
  const hasChainFilter = isDefined(row.label_asym_id) || isDefined(row.auth_asym_id) || isDefined(row.label_entity_id);

  if (hasAtomIds) {
    const theAtom = getTheAtomForRow(model, row, atomicIndices);
    if (theAtom === undefined) return undefined;
    return out.add(theAtom, theAtom + 1);
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
  rows: readonly MVSAnnotationRow[],
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
): ElementRanges {
  // Per-row results go to a shared buffer (`SharedScratch`), are collected in a shared bag (`RowsBag`) and merged once;
  // only the final result is a new object (callers may keep it).
  const bag = RowsBag.clear();
  for (const row of rows) {
    const r = getAtomRangesForRowInto(row, model, instanceId, indices, SharedScratch);
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
  row: MVSAnnotationRow,
  indices: AtomicIndicesAndSortings,
  out: IntBuffer,
): ArrayLike<number> | undefined {
  const { auth_asym_id, label_entity_id } = model.atomicHierarchy.chains;
  let result: ArrayLike<number> | undefined = undefined;
  if (isDefined(row.label_asym_id)) {
    result = indices.chainsByLabelAsymId.get(row.label_asym_id) ?? EmptyArray;
  }
  if (isDefined(row.auth_asym_id)) {
    if (result) {
      result = filterByColumn(result, out, auth_asym_id, row.auth_asym_id);
    } else {
      result = indices.chainsByAuthAsymId.get(row.auth_asym_id) ?? EmptyArray;
    }
  }
  if (isDefined(row.label_entity_id)) {
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
  LabelCompId,
  AuthCompId,
}

/** Write residue indexes (ascending) which satisfy criteria given by `row` into `out` (clearing it first) and return it.
 * `fromChains` are ascending chain indexes to search in, `undefined` means all chains.
 *
 * For each chain, the candidates are taken from the index of the first criterion present (so the work is proportional
 * to the number of candidates, not the size of the chain) and all remaining criteria are checked in one pass,
 * without intermediate arrays. */
function getQualifyingResidues(
  model: Model,
  row: MVSAnnotationRow,
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

  const hasLabelSeq = isDefined(row.label_seq_id);
  const hasAuthSeq = isDefined(row.auth_seq_id);
  const hasResIdx = isDefined(row.residue_index);
  const hasIns = isDefined(row.pdbx_PDB_ins_code);
  const hasBegLabel = isDefined(row.beg_label_seq_id);
  const hasEndLabel = isDefined(row.end_label_seq_id);
  const hasBegAuth = isDefined(row.beg_auth_seq_id);
  const hasEndAuth = isDefined(row.end_auth_seq_id);
  const hasLabelRange = hasBegLabel || hasEndLabel;
  const hasAuthRange = hasBegAuth || hasEndAuth;
  const hasLabelComp = isDefined(row.label_comp_id);
  const hasAuthComp = isDefined(row.auth_comp_id);
  const labelCompPure = indices.residuesByLabelCompIdIsPure;
  const authCompPure = indices.residuesByAuthCompIdIsPure;

  const source = hasLabelSeq
    ? ResidueSource.LabelSeqId
    : hasAuthSeq
      ? ResidueSource.AuthSeqId
      : hasResIdx
        ? ResidueSource.ResidueIndex
        : hasIns
          ? ResidueSource.InsCode
          : hasLabelRange
            ? ResidueSource.LabelSeqIdRange
            : hasAuthRange
              ? ResidueSource.AuthSeqIdRange
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
  if (!fromChains) {
    const global = indices.global;
    const bySeq = hasLabelSeq || hasAuthSeq || hasLabelRange || hasAuthRange;
    if (hasResIdx && global.residueBySourceIndex) {
      const iRes = global.residueBySourceIndex.get(row.residue_index!);
      SingleResidue[0] = iRes as ResidueIndex;
      globalKeys = iRes === undefined ? EmptyArray : SingleResidue;
      covered = ResidueSource.ResidueIndex;
    } else if (hasIns && !bySeq) {
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
          else Sorting.findIndexRange(sorting, row.beg_label_seq_id, row.end_label_seq_id, IndexRange);
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
          else Sorting.findIndexRange(sorting, row.beg_auth_seq_id, row.end_auth_seq_id, IndexRange);
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

/** Return index of atom in `model` which satistfies criteria given by `row`, if any.
 * Only works when `row.atom_id` and/or `row.atom_index` is defined (otherwise use `getAtomRangesForRow`). */
function getTheAtomForRow(
  model: Model,
  row: MVSAnnotationRow,
  indices: AtomicIndicesAndSortings,
): ElementIndex | undefined {
  let iAtom: ElementIndex | undefined = undefined;
  if (!isDefined(row.atom_id) && !isDefined(row.atom_index))
    throw new Error('ArgumentError: at least one of row.atom_id, row.atom_index must be defined.');
  if (isDefined(row.atom_id) && isDefined(row.atom_index)) {
    const a1 = indices.atomsById.get(row.atom_id);
    const a2 = indices.atomsBySourceIndex.get(row.atom_index);
    if (a1 !== a2) return undefined;
    iAtom = a1;
  }
  if (isDefined(row.atom_id)) {
    iAtom = indices.atomsById.get(row.atom_id);
  }
  if (isDefined(row.atom_index)) {
    iAtom = indices.atomsBySourceIndex.get(row.atom_index);
  }
  if (iAtom === undefined) return undefined;
  if (!atomQualifies(model, iAtom, row)) return undefined;
  return iAtom;
}

/** Return true if `iAtom`-th atom in `model` satisfies all selection criteria given by `row`. */
export function atomQualifies(model: Model, iAtom: ElementIndex, row: MVSAnnotationRow): boolean {
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
  if (!matchesRange(row.beg_auth_seq_id, row.end_auth_seq_id, auth_seq_id)) return false;
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
  return !isDefined(requiredValue) || value === requiredValue;
}

/** Return true if `requiredMin <= value <= requiredMax`.
 * Undefined `requiredMin` behaves like negative infinity.
 * Undefined `requiredMax` behaves like positive infinity. */
function matchesRange<T>(
  requiredMin: T | undefined | null,
  requiredMax: T | undefined | null,
  value: T | undefined,
): boolean {
  if (isDefined(requiredMin) && (!isDefined(value) || value < requiredMin)) return false;
  if (isDefined(requiredMax) && (!isDefined(value) || value > requiredMax)) return false;
  return true;
}

// COARSE SELECTIONS

/** Return sphere ranges in `model` which satisfy criteria given by `row` */
export function getSphereRangesForRow(
  row: MVSAnnotationRow,
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
): ElementRanges | undefined {
  if (!indices.spheres) return undefined;
  if (isDefined(row.instance_id) && row.instance_id !== instanceId) return undefined;
  return getCoarseElementRangesForRow(row, model.coarseHierarchy.spheres, indices.spheres);
}

/** Return sphere ranges in `model` which satisfy criteria given by any of `rows` (spheres that satisfy more rows are still included only once) */
export function getSphereRangesForRows(
  rows: readonly MVSAnnotationRow[],
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
): ElementRanges {
  return ElementRanges.union(rows.map((row) => getSphereRangesForRow(row, model, instanceId, indices)));
}

/** Return gaussian ranges in `model` which satisfy criteria given by `row` */
export function getGaussianRangesForRow(
  row: MVSAnnotationRow,
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
): ElementRanges | undefined {
  if (!indices.gaussians) return undefined;
  if (isDefined(row.instance_id) && row.instance_id !== instanceId) return undefined;
  return getCoarseElementRangesForRow(row, model.coarseHierarchy.gaussians, indices.gaussians);
}

/** Return gaussian ranges in `model` which satisfy criteria given by any of `rows` (gaussians that satisfy more rows are still included only once) */
export function getGaussianRangesForRows(
  rows: readonly MVSAnnotationRow[],
  model: Model,
  instanceId: string,
  indices: IndicesAndSortings,
): ElementRanges {
  return ElementRanges.union(rows.map((row) => getGaussianRangesForRow(row, model, instanceId, indices)));
}

/** Return ranges of coarse elements (spheres or gaussians) which satisfy criteria given by `row` */
export function getCoarseElementRangesForRow(
  row: MVSAnnotationRow,
  coarseElements: CoarseElements,
  indices: CoarseIndicesAndSortings,
): ElementRanges | undefined {
  const nElements = coarseElements.count;
  if (nElements === 0) return undefined;

  const hasResidueFilter = isAnyDefined(row.label_seq_id, row.beg_label_seq_id, row.end_label_seq_id);
  const hasChainFilter = isAnyDefined(row.label_asym_id, row.label_entity_id);
  const hasInvalidFilter = isAnyDefined(
    row.auth_asym_id,
    row.auth_seq_id,
    row.pdbx_PDB_ins_code,
    row.beg_auth_seq_id,
    row.end_auth_seq_id,
    row.label_comp_id,
    row.auth_comp_id,
    row.residue_index,
    row.label_atom_id,
    row.auth_atom_id,
    row.type_symbol,
    row.atom_id,
    row.atom_index,
  );

  if (hasInvalidFilter) {
    printCoarseSelectorWarning();
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
  row: MVSAnnotationRow,
  indices: CoarseIndicesAndSortings,
): readonly ChainIndex[] {
  let result: readonly ChainIndex[] | undefined = undefined;
  if (isDefined(row.label_asym_id)) {
    result = indices.chainsByAsymId.get(row.label_asym_id) ?? EmptyArray;
  }
  if (isDefined(row.label_entity_id)) {
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
  row: MVSAnnotationRow,
  indices: CoarseIndicesAndSortings,
  fromChains: readonly ChainIndex[],
): ElementIndex[] {
  const result: ElementIndex[] = [];
  for (const iChain of fromChains) {
    const sorting = indices.elementsSortedBySeqIdBegin.get(iChain)!;
    const queryStart = Math.max(row.label_seq_id ?? -Infinity, row.beg_label_seq_id ?? -Infinity);
    const queryEnd = Math.min(row.label_seq_id ?? Infinity, row.end_label_seq_id ?? Infinity); // inclusive
    const iStart = SortedArray.findPredecessorIndex(sorting.endUpperBounds, queryStart); // select elements potentially ending >=queryStart (necessary condition)
    const iStop = SortedArray.findPredecessorIndex(sorting.values, queryEnd + 1); // select elements starting <=queryEnd (necessary and suffient condition) // exclusive

    for (let i = iStart; i < iStop; i++) {
      const iElem = sorting.keys[i];
      if (coarseElements.seq_id_end.value(iElem) >= queryStart) {
        // rechecking seq_id_end, as the condition was not sufficient
        result.push(iElem);
      }
    }
    // This implementation can yield some elements even when queryStart>queryEnd (e.g. { beg_label_seq_id: 70, end_label_seq_id: 58, label_seq_id: 60 } -> sphere 51-100 qualifies ).
    // This is on purpose, to have the same behavior as MolScript.
  }
  sortIfNeeded(result, (a, b) => a - b);
  return result;
}

let coarseSelectorWarningPrinted = false;
function printCoarseSelectorWarning() {
  if (!coarseSelectorWarningPrinted) {
    console.warn(
      'Using unsupported selector fields (auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, beg_auth_seq_id, end_auth_seq_id, label_comp_id, auth_comp_id, residue_index, label_atom_id, auth_atom_id, type_symbol, atom_id, atom_index) on a coarse structure. The resulting selection will be empty.',
    );
    coarseSelectorWarningPrinted = true;
  }
}

// GENERAL

/** Convert an annotation row into a MolScript expression */
export function rowToExpression(row: MVSAnnotationRow): Expression {
  return StructureElement.Schema.toExpression(row);
}

/** Convert multiple annotation rows into a MolScript expression.
 * (with union semantics, i.e. an atom qualifies if it qualifies for at least one of the rows) */
export function rowsToExpression(rows: readonly MVSAnnotationRow[]): Expression {
  return StructureElement.Schema.toExpression({
    items: rows as StructureElement.SchemaItem[],
  });
}
