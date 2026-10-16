/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MmcifFormat } from '@molstar/model/formats/structure/mmcif';
import type { ChainIndex, Model, ResidueIndex } from '@molstar/model/model/structure';
import { AtomicHierarchy } from '@molstar/model/model/structure/model/properties/atomic';
import { ApiArea, ApiResult } from '@molstar/prism-viewer/api/area';
import { StructureParam } from '@molstar/prism-viewer/api/params';
import { isFailure, labelOf, oneStructure } from '@molstar/prism-viewer/api/scope';
import { QualityAssessment } from '../quality-assessment/prop.js';
import { type PaeMatrix, type PaeResidue, paeChains } from './data.js';

// The `model-archive` area (operations.md §2, an extension's own area): the predicted aligned error of a model, read
// from the ModelCIF's `ma_qa_metric_local_pairwise` (Model Archive's models carry it; AlphaFold DB ships its PAE apart
// from the model). `model-archive.pae` says whether a structure has it and how big it is; the PAE card reads the matrix
// through `paeMatrix` (a read for views, as the animation area's `state`), and acts through `select.set`; the Structure
// panel's row says its size through `paeSize`, which does not read the matrix.

/** What `model-archive.pae` gives. */
export interface PaeInfo {
  readonly structure: string;
  readonly hasPae: boolean;
  readonly metric?: string;
  /** Residues on each side. */
  readonly size?: number;
  /** The highest error, Å. */
  readonly max?: number;
  readonly chains?: readonly string[];
}

/** What the Structure panel's PAE row says of a matrix: its residues on each side and their chains. */
export interface PaeSize {
  readonly residues: number;
  readonly chains: number;
}

/**
 * The model's polymer residues and their chains (by author name): what a ModelCIF's pairwise metric covers (its rows
 * name residues by `label_seq_id`), counted from the hierarchy without reading the N² pairwise rows.
 */
function polymerSizeOf(model: Model): PaeSize {
  const h = model.atomicHierarchy;
  const types = model.entities.data.type;
  const chains = new Set<string>();
  let residues = 0;
  for (let c = 0 as ChainIndex, n = h.chains._rowCount; c < n; c++) {
    const e = h.index.getEntityFromChain(c);
    if (e < 0 || types.value(e) !== 'polymer') continue;
    const count = AtomicHierarchy.chainEndResidueIndexExcl(h, c) - AtomicHierarchy.chainStartResidueIndex(h, c);
    if (count <= 0) continue;
    residues += count;
    chains.add(h.chains.auth_asym_id.value(c));
  }
  return { residues, chains: chains.size };
}

/** The model's PAE metric as a dense matrix; undefined when it has none. */
function matrixOf(model: Model): PaeMatrix | undefined {
  if (!MmcifFormat.is(model.sourceData)) return undefined;
  const frame = model.sourceData.data.frame;
  const found = QualityAssessment.findModelArchiveCIFPAEMetrics(frame)[0];
  if (!found) return undefined;
  const metric = QualityAssessment.pairwiseMetricFromModelArchiveCIF(model, frame, found.id);
  if (!metric) return undefined;
  const [lo, hi] = metric.residueRange;
  if (hi < lo) return undefined;
  const size = hi - lo + 1;
  const values = new Float32Array(size * size).fill(Number.NaN);
  for (let i = 0; i < size; i++) {
    const row = metric.values[(lo + i) as ResidueIndex];
    if (!row) continue;
    for (let j = 0; j < size; j++) {
      const v = row[(lo + j) as ResidueIndex];
      if (typeof v === 'number') values[i * size + j] = v;
    }
  }
  const h = model.atomicHierarchy;
  const residues: PaeResidue[] = [];
  for (let k = 0; k < size; k++) {
    const r = (lo + k) as ResidueIndex;
    const c = AtomicHierarchy.residueChainIndex(h, r);
    residues.push({
      chain: h.chains.auth_asym_id.value(c),
      seq: h.residues.auth_seq_id.value(r),
      insCode: h.residues.pdbx_PDB_ins_code.value(r) ?? '',
      name: h.atoms.label_comp_id.value(h.residueAtomSegments.offsets[r]),
    });
  }
  return { name: metric.name, size, values, max: metric.valueRange[1], residues, plddt: plddtOf(model, lo, size) };
}

/** The pLDDT of residues `lo`… (0–100), from the model's local metric; undefined when it has none. */
function plddtOf(model: Model, lo: number, size: number): Float32Array | undefined {
  const metric = QualityAssessment.fromModel(model).local.find((m) => m.kind === 'pLDDT');
  if (!metric?.values.size) return undefined;
  // A metric in [0, 1] is scaled to the bands' 0–100.
  const scale = metric.domain?.[1] === 1 ? 100 : 1;
  const out = new Float32Array(size).fill(Number.NaN);
  for (let k = 0; k < size; k++) {
    const v = metric.values.get((lo + k) as ResidueIndex);
    if (typeof v === 'number') out[k] = v * scale;
  }
  return out;
}

export class ModelArchiveApi extends ApiArea {
  readonly id = 'model-archive';

  /** By model: a matrix is read once (null: the model has none). */
  readonly #matrices = new WeakMap<Model, PaeMatrix | null>();

  #matrix(model: Model): PaeMatrix | undefined {
    let m = this.#matrices.get(model);
    if (m === undefined) {
      m = matrixOf(model) ?? null;
      this.#matrices.set(model, m);
    }
    return m ?? undefined;
  }

  #modelOf(ref: string): Model | undefined {
    const s = this.ctx.plugin.managers.structure.hierarchy.current.structures.find((x) => x.cell.transform.ref === ref);
    return s?.cell.obj?.data.model;
  }

  #hasPae(model: Model | undefined): model is Model {
    if (!model || !MmcifFormat.is(model.sourceData)) return false;
    const known = this.#matrices.get(model);
    if (known !== undefined) return !!known;
    return QualityAssessment.findModelArchiveCIFPAEMetrics(model.sourceData.data.frame).length > 0;
  }

  /** Whether a structure (by state ref) has PAE, without reading the matrix: for the card's toolbar button. */
  hasPae(ref: string): boolean {
    return this.#hasPae(this.#modelOf(ref));
  }

  /**
   * A structure's PAE size (by state ref), for its row in the Structure panel; undefined when it has none. The matrix's
   * once it has been read; before, the model's polymer residues, so a row does not read the N² pairwise rows of a
   * metric nobody has opened.
   */
  paeSize(ref: string): PaeSize | undefined {
    const model = this.#modelOf(ref);
    if (!this.#hasPae(model)) return undefined;
    const m = this.#matrices.get(model);
    if (!m) return polymerSizeOf(model);
    return { residues: m.size, chains: new Set(m.residues.map((r) => r.chain)).size };
  }

  /** A structure's PAE matrix (by state ref), for the card; undefined when it has none. */
  paeMatrix(ref: string): PaeMatrix | undefined {
    const model = this.#modelOf(ref);
    return model ? this.#matrix(model) : undefined;
  }

  readonly pae = this.operation({
    id: 'model-archive.pae',
    title: 'Read the PAE of a model',
    params: { structure: StructureParam() },
    read: true,
    run: async ({ structure }, ctx): Promise<ApiResult<PaeInfo>> => {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const label = labelOf(ctx, s);
      const model = s.cell.obj?.data.model;
      const m = model && this.#matrix(model);
      if (!m) {
        return ApiResult.ok(
          `${label} has no PAE (a ModelCIF's ma_qa_metric_local_pairwise; AlphaFold DB ships it apart from the model).`,
          {},
          { value: { structure: label, hasPae: false } },
        );
      }
      const chains = [...new Set(paeChains(m).map((c) => c.label))];
      return ApiResult.ok(
        `${label}: ${m.name}, ${m.size} × ${m.size} residues (chains ${chains.join(', ')}), up to ${m.max.toFixed(1)} Å.`,
        {},
        { value: { structure: label, hasPae: true, metric: m.name, size: m.size, max: m.max, chains } },
      );
    },
  });
}
