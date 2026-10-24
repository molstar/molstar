/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type StateBuilder, StateTransform } from '@molstar/core/state';
import type { CifFrame } from '@molstar/io/reader/cif';
import { MmcifFormat } from '@molstar/model/formats/structure/mmcif';
import type { ChainIndex, Model, ResidueIndex } from '@molstar/model/model/structure';
import { AtomicHierarchy } from '@molstar/model/model/structure/model/properties/atomic';
import { type ApiAreaContext, type ApiContext, ApiArea, ApiResult } from '@molstar/prism-viewer/api/area';
import { StructureParam } from '@molstar/prism-viewer/api/params';
import { isFailure, labelOf, oneStructure } from '@molstar/prism-viewer/api/scope';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { type Observable, type Subscription, Subject } from 'rxjs';
import { QualityAssessment } from '../quality-assessment/prop.js';
import { type AfdbEntry, AfdbError, afdbEntryOf, paeValuesFromAfdb, resolveAfdbPaeUrl } from './afdb.js';
import { type AssociatedPae, associatedPaeOf } from './associated.js';
import { type PaeMatrix, type PaeResidue, paeChains } from './data.js';
import { FetchAssociatedPae, FetchPaeJson } from './fetch.js';

// The `model-archive` area (api.md §2, an extension's own area): the predicted aligned error of a model, read
// from the ModelCIF's `ma_qa_metric_local_pairwise`. A model may carry it; Model Archive often ships it apart, in an
// associated file the model names (`associated.ts`), and AlphaFold DB as a JSON beside the model. `model-archive.fetch-pae`
// gets it from either: a download and a parse in one cell under the model's cell (`fetch.ts`), tagged with the model:
// it goes with the structure, is undoable and saved in a session. `model-archive.pae` says whether a structure has it and how big it is; the PAE card reads the matrix
// through `paeMatrix` (a read for views, as the animation area's `state`), and acts through `select.set`; the Structure
// panel's row says its size through `paeSize`, which does not read the matrix.

/** What `model-archive.pae` gives. */
export interface PaeInfo {
  readonly structure: string;
  readonly hasPae: boolean;
  /** Set when a fetched PAE is of another version than the model file's. */
  readonly note?: string;
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

/**
 * Whether the ModelCIF's pairwise rows cover every pair of the model's polymer residues, once per metric: then the
 * polymer's size is the metric's (counted without reading a row); a partial metric (a domain, some chains) is read.
 */
function coversPolymer(model: Model, size: PaeSize): boolean {
  if (!MmcifFormat.is(model.sourceData)) return false;
  const { rows, metrics } = QualityAssessment.pairwiseSize(model.sourceData.data.frame);
  return metrics > 0 && rows === size.residues * size.residues * metrics;
}

/** The PAE metric of a ModelCIF frame (the model's own, or its associated file's) as a dense matrix over the model. */
function matrixOf(model: Model, frame: CifFrame): PaeMatrix | undefined {
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
  return {
    name: metric.name,
    size,
    values,
    max: metric.valueRange[1],
    residues,
    plddt: plddtOf(model, size, (k) => (lo + k) as ResidueIndex),
  };
}

/** The pLDDT of the matrix's residues (0–100), from the model's local metric; undefined when it has none. */
function plddtOf(model: Model, size: number, residueAt: (k: number) => ResidueIndex): Float32Array | undefined {
  const metric = QualityAssessment.fromModel(model).local.find((m) => m.kind === 'pLDDT');
  if (!metric?.values.size) return undefined;
  // A metric in [0, 1] is scaled to the bands' 0–100.
  const scale = metric.domain?.[1] === 1 ? 100 : 1;
  const out = new Float32Array(size).fill(Number.NaN);
  for (let k = 0; k < size; k++) {
    const v = metric.values.get(residueAt(k));
    if (typeof v === 'number') out[k] = v * scale;
  }
  return out;
}

/** The tag of the cells that hold an AlphaFold DB PAE (a download and its parsed JSON). */
export const AfdbPaeTag = 'model-archive.afdb-pae';
/** The tag of the cells that hold a model's associated PAE file (Model Archive's, parsed CIF). */
export const AssociatedPaeTag = 'model-archive.associated-pae';
/** A fetched cell's own tag: which model's it is (`<tag>:<model cell ref>`). */
const tagOf = (tag: string, modelRef: string) => `${tag}:${modelRef}`;

/**
 * AlphaFold DB's PAE for a model: its polymer residues, in order, are the matrix's tokens. A message when the JSON is
 * not a PAE or does not fit the model (another entry's, another length).
 */
function attachedMatrix(model: Model, json: unknown): PaeMatrix | string {
  const p = paeValuesFromAfdb(json);
  if (typeof p === 'string') return p;
  const h = model.atomicHierarchy;
  const types = model.entities.data.type;
  const residueIndices: ResidueIndex[] = [];
  for (let c = 0 as ChainIndex, n = h.chains._rowCount; c < n; c++) {
    const e = h.index.getEntityFromChain(c);
    if (e < 0 || types.value(e) !== 'polymer') continue;
    for (
      let r = AtomicHierarchy.chainStartResidueIndex(h, c), end = AtomicHierarchy.chainEndResidueIndexExcl(h, c);
      r < end;
      r++
    )
      residueIndices.push(r as ResidueIndex);
  }
  if (residueIndices.length !== p.size)
    return `the PAE has ${p.size} residues, the model ${residueIndices.length} polymer residues`;
  const residues = residueIndices.map((r): PaeResidue => {
    const c = AtomicHierarchy.residueChainIndex(h, r);
    return {
      chain: h.chains.auth_asym_id.value(c),
      seq: h.residues.auth_seq_id.value(r),
      insCode: h.residues.pdbx_PDB_ins_code.value(r) ?? '',
      name: h.atoms.label_comp_id.value(h.residueAtomSegments.offsets[r]),
    };
  });
  return {
    name: 'PAE',
    size: p.size,
    values: p.values,
    max: p.max,
    residues,
    plddt: plddtOf(model, p.size, (k) => residueIndices[k]),
  };
}

/** A fetched PAE in the state tree: AlphaFold DB's parsed JSON, or an associated file's CIF frame. */
type Fetched =
  | { readonly kind: 'afdb'; readonly data: object }
  | { readonly kind: 'associated'; readonly data: CifFrame };

/** A transform's address parameter (a text or an asset's `{ url }`), if any. */
function urlOf(params: unknown): string | undefined {
  const url = (params as { url?: unknown } | undefined)?.url;
  if (typeof url === 'string') return url || undefined;
  const inner = (url as { url?: unknown } | undefined)?.url;
  return typeof inner === 'string' && inner ? inner : undefined;
}

/** What one pass over the state tree finds: the structures' models, the fetched PAE cells, where a PAE can come from. */
interface TreeIndex {
  readonly structures: Map<string, { readonly model: Model; readonly modelRef: string }>;
  /** A model's fetched PAE, by the model's ref. */
  readonly fetched: Map<string, Fetched>;
  /** By structure ref, read when asked. */
  readonly entries: Map<string, AfdbEntry | undefined>;
  readonly associated: Map<string, AssociatedPae | undefined>;
}

/** Where a fetch gets a PAE from, and how it reads it. */
interface PaeSource {
  /** Named in messages ("ma-bak-cepc-0003", "AF-P04637-F1"). */
  readonly id: string;
  /** "Model Archive", "AlphaFold DB". */
  readonly service: string;
  readonly url: string;
  readonly note?: string;
  readonly tag: string;
  /** Adds the cell under the model; its ref. */
  apply(b: StateBuilder.Root, modelRef: string): string;
  /** The matrix from the cell's object, or why it does not fit the model; `key` keys the area's cache. */
  read(model: Model, obj: unknown): { key: object; matrix: PaeMatrix } | string;
}

export class ModelArchiveApi extends ApiArea {
  readonly id = 'model-archive';

  /** By model: a matrix is read once (null: the model has none). */
  readonly #matrices = new WeakMap<Model, PaeMatrix | null>();
  /** By model: whether its file carries a PAE metric (the metric list is read once). */
  readonly #carries = new WeakMap<Model, boolean>();
  /** By parsed JSON or CIF frame (a fetched PAE in the state tree): its matrix, read once (null: it does not fit). */
  #attachedMatrices = new WeakMap<object, PaeMatrix | null>();
  readonly #changes = new Subject<void>();
  readonly #subscription: Subscription;
  /** The tree read in one pass, once per change (dropped by the next). */
  #index: TreeIndex | undefined;
  /** By model ref: the fetch in progress, so a second ask shares it instead of adding the cells twice. */
  readonly #fetching = new Map<string, Promise<ApiResult<PaeInfo>>>();

  /** Emits when the state changes, so a view re-asks what has a PAE (a fetched one appears without a structure change). */
  get changes(): Observable<void> {
    return this.#changes;
  }

  constructor(ctx: ApiAreaContext) {
    super(ctx);
    // The tree changes in one step for a fetch (download, parse); a view asks again after it.
    this.#subscription = ctx.plugin.state.data.events.changed.subscribe(() => {
      this.#index = undefined;
      this.#changes.next();
    });
  }

  dispose(): void {
    this.#subscription.unsubscribe();
    this.#index = undefined;
    this.#attachedMatrices = new WeakMap();
    this.#fetching.clear();
    this.#changes.complete();
  }

  /** The tree in one pass (structures with their models, fetched PAE cells), kept until the state changes. */
  #tree(): TreeIndex {
    if (this.#index) return this.#index;
    const state = this.ctx.plugin.state.data;
    const structures = new Map<string, { model: Model; modelRef: string }>();
    const fetched = new Map<string, Fetched>();
    for (const cell of state.cells.values()) {
      const obj = cell.obj;
      if (PluginStateObject.Molecule.Structure.is(obj)) {
        const modelRef = this.#modelRefOf(cell.transform.ref);
        if (modelRef) structures.set(cell.transform.ref, { model: obj.data.model, modelRef });
      } else if (PluginStateObject.Format.Json.is(obj) && cell.transform.tags?.includes(AfdbPaeTag)) {
        const data: unknown = obj.data;
        if (data && typeof data === 'object') fetched.set(cell.transform.parent, { kind: 'afdb', data });
      } else if (PluginStateObject.Format.Cif.is(obj) && cell.transform.tags?.includes(AssociatedPaeTag)) {
        const frame = obj.data.blocks[0];
        if (frame) fetched.set(cell.transform.parent, { kind: 'associated', data: frame });
      }
    }
    return (this.#index = { structures, fetched, entries: new Map(), associated: new Map() });
  }

  #modelOf(ref: string): Model | undefined {
    return this.#tree().structures.get(ref)?.model;
  }

  /** The state ref of a structure's model (its ancestor). */
  #modelRefOf(ref: string): string | undefined {
    const state = this.ctx.plugin.state.data;
    for (
      let cell = state.cells.get(ref);
      cell && cell.transform.ref !== StateTransform.RootRef;
      cell = state.cells.get(cell.transform.parent)
    ) {
      if (PluginStateObject.Molecule.Model.is(cell.obj)) return cell.transform.ref;
    }
    return undefined;
  }

  /** The PAE fetched for a structure's model (its parsed cell), if any. */
  #attached(ref: string): { fetched: Fetched; model: Model } | undefined {
    const t = this.#tree();
    const s = t.structures.get(ref);
    const fetched = s && t.fetched.get(s.modelRef);
    return s && fetched ? { fetched, model: s.model } : undefined;
  }

  #matrix(model: Model): PaeMatrix | undefined {
    let m = this.#matrices.get(model);
    if (m === undefined) {
      m = (MmcifFormat.is(model.sourceData) && matrixOf(model, model.sourceData.data.frame)) || null;
      this.#matrices.set(model, m);
    }
    return m ?? undefined;
  }

  #embedded(model: Model | undefined): model is Model {
    if (!model || !MmcifFormat.is(model.sourceData)) return false;
    const known = this.#matrices.get(model);
    if (known !== undefined) return !!known;
    let has = this.#carries.get(model);
    if (has === undefined) {
      has = QualityAssessment.findModelArchiveCIFPAEMetrics(model.sourceData.data.frame).length > 0;
      this.#carries.set(model, has);
    }
    return has;
  }

  #fetched(ref: string): PaeMatrix | undefined {
    const a = this.#attached(ref);
    if (!a) return undefined;
    const { fetched, model } = a;
    let m = this.#attachedMatrices.get(fetched.data);
    if (m === undefined) {
      const r =
        fetched.kind === 'afdb'
          ? attachedMatrix(model, fetched.data)
          : (matrixOf(model, fetched.data) ?? 'it names none of its residues');
      if (typeof r === 'string') this.ctx.plugin.log.warn(`The fetched PAE does not fit the model: ${r}.`);
      m = typeof r === 'string' ? null : r;
      this.#attachedMatrices.set(fetched.data, m);
    }
    return m ?? undefined;
  }

  /** Whether a structure (by state ref) has PAE, without reading the matrix: for the card's toolbar button. */
  hasPae(ref: string): boolean {
    return this.#embedded(this.#modelOf(ref)) || !!this.#attached(ref);
  }

  /**
   * The associated file a structure's model names for its PAE (Model Archive's, `associated.ts`), by state ref;
   * undefined when it names none. A relative address is resolved against the one the model was loaded from.
   */
  associatedPae(ref: string): AssociatedPae | undefined {
    const known = this.#tree().associated;
    if (!known.has(ref)) {
      const model = this.#modelOf(ref);
      known.set(
        ref,
        model && MmcifFormat.is(model.sourceData)
          ? associatedPaeOf(model.sourceData.data.frame, this.#loadedFrom(ref))
          : undefined,
      );
    }
    return known.get(ref);
  }

  /** The nearest address among a structure's ancestors' parameters (what it was downloaded from), if any. */
  #loadedFrom(ref: string): string | undefined {
    const state = this.ctx.plugin.state.data;
    for (
      let cell = state.cells.get(ref);
      cell && cell.transform.ref !== StateTransform.RootRef;
      cell = state.cells.get(cell.transform.parent)
    ) {
      const url = urlOf(cell.transform.params);
      if (url) return url;
    }
    return undefined;
  }

  /**
   * The AlphaFold DB entry a structure is (by state ref), from its entry id, labels, and the address or file it came
   * from; undefined for any other structure.
   */
  afdbEntry(ref: string): AfdbEntry | undefined {
    const known = this.#tree().entries;
    if (!known.has(ref)) known.set(ref, this.#afdbEntryOf(ref));
    return known.get(ref);
  }

  #afdbEntryOf(ref: string): AfdbEntry | undefined {
    const model = this.#modelOf(ref);
    if (!model) return undefined;
    const texts: (string | undefined)[] = [model.entryId, model.label];
    if (MmcifFormat.is(model.sourceData)) {
      const frame = model.sourceData.data.frame;
      texts.push(frame.header, frame.categories.entry?.getField('id')?.str(0));
    }
    // What the model was loaded from: its ancestors' addresses, file names and labels, nearest first.
    const state = this.ctx.plugin.state.data;
    for (
      let cell = state.cells.get(ref);
      cell && cell.transform.ref !== StateTransform.RootRef;
      cell = state.cells.get(cell.transform.parent)
    ) {
      const p = cell.transform.params as { file?: { name?: string }; label?: string } | undefined;
      texts.push(cell.obj?.label, urlOf(p), p?.file?.name, p?.label);
    }
    return afdbEntryOf(texts);
  }

  /**
   * Whether the PAE of a structure can be fetched: it has none yet, and its model names an associated file with it
   * (Model Archive) or it is an AlphaFold DB entry.
   */
  canFetchPae(ref: string): boolean {
    return !this.hasPae(ref) && (!!this.associatedPae(ref) || !!this.afdbEntry(ref));
  }

  /** Where `model-archive.fetch-pae` would get a structure's PAE from (by state ref): the service's name, if anywhere. */
  paeService(ref: string): 'Model Archive' | 'AlphaFold DB' | undefined {
    return this.associatedPae(ref) ? 'Model Archive' : this.afdbEntry(ref) ? 'AlphaFold DB' : undefined;
  }

  /**
   * A structure's PAE size (by state ref), for its row in the Structure panel; undefined when it has none. The matrix's
   * once it has been read; before, the model's polymer residues when its metric covers them all (counted from the
   * hierarchy and the metric's row count, so a row does not read the N² pairwise rows of a metric nobody has opened),
   * else (a partial metric) the matrix is read for its size.
   */
  paeSize(ref: string): PaeSize | undefined {
    const model = this.#modelOf(ref);
    if (this.#embedded(model)) {
      const m = this.#matrices.get(model);
      if (m) return { residues: m.size, chains: new Set(m.residues.map((r) => r.chain)).size };
      const polymer = polymerSizeOf(model);
      if (coversPolymer(model, polymer)) return polymer;
      const read = this.#matrix(model);
      return read && { residues: read.size, chains: new Set(read.residues.map((r) => r.chain)).size };
    }
    const m = this.#fetched(ref);
    return m && { residues: m.size, chains: new Set(m.residues.map((r) => r.chain)).size };
  }

  /**
   * A structure's selection (by state ref) as a token: the same object while the selection is unchanged (Mol*'s loci of
   * it), another once anything changes it. The PAE card takes its pair off again only while the selection is the one it
   * made.
   */
  selectionToken(ref: string): object | undefined {
    const s = this.ctx.plugin.state.data.cells.get(ref)?.obj;
    if (!PluginStateObject.Molecule.Structure.is(s)) return undefined;
    return this.ctx.plugin.managers.structure.selection.getLoci(s.data);
  }

  /** A structure's PAE matrix (by state ref), for the card; undefined when it has none. */
  paeMatrix(ref: string): PaeMatrix | undefined {
    const model = this.#modelOf(ref);
    return (model && this.#matrix(model)) || this.#fetched(ref);
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
      const m = this.paeMatrix(s.cell.transform.ref);
      if (!m) {
        const ref = s.cell.transform.ref;
        const where = this.associatedPae(ref)
          ? 'the file Model Archive names for it'
          : this.afdbEntry(ref)
            ? `AlphaFold DB (${this.afdbEntry(ref)!.id})`
            : undefined;
        return ApiResult.ok(
          `${label} has no PAE (a ModelCIF's ma_qa_metric_local_pairwise; Model Archive and AlphaFold DB often ship it apart from the model).`,
          {},
          {
            value: { structure: label, hasPae: false },
            ...(where ? { next: `model-archive.fetch-pae gets it from ${where}.` } : {}),
          },
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

  /** Whether the structure still shows this model: the structure and its model are in the tree, unreplaced. */
  #alive(ref: string, modelRef: string, model: Model): boolean {
    const state = this.ctx.plugin.state.data;
    const structure = state.cells.get(ref)?.obj;
    return (
      PluginStateObject.Molecule.Structure.is(structure) &&
      structure.data.model === model &&
      state.cells.get(modelRef)?.obj?.data === model
    );
  }

  /** Where an AlphaFold DB entry's PAE is (its API asked), read as the JSON's matrix over the model. */
  async #afdbSource(entry: AfdbEntry, signal: AbortSignal): Promise<PaeSource> {
    const { url, note } = await resolveAfdbPaeUrl(entry, undefined, signal);
    return {
      id: entry.id,
      service: 'AlphaFold DB',
      url,
      note,
      tag: AfdbPaeTag,
      apply: (b, modelRef) =>
        b
          .to(modelRef)
          .apply(FetchPaeJson, { url, label: `${entry.id} PAE` }, { tags: [AfdbPaeTag, tagOf(AfdbPaeTag, modelRef)] })
          .ref,
      read: (model, obj) => {
        const json = (obj as PluginStateObject.Format.Json).data as object;
        const m = attachedMatrix(model, json);
        return typeof m === 'string' ? m : { key: json, matrix: m };
      },
    };
  }

  /** A model's associated PAE file (Model Archive's), read as its pairwise metric over the model. */
  #associatedSource(id: string, file: AssociatedPae): PaeSource {
    return {
      id,
      service: 'Model Archive',
      url: file.url,
      tag: AssociatedPaeTag,
      apply: (b, modelRef) =>
        b
          .to(modelRef)
          .apply(
            FetchAssociatedPae,
            { url: file.url, path: file.path, label: `${id} PAE` },
            { tags: [AssociatedPaeTag, tagOf(AssociatedPaeTag, modelRef)] },
          ).ref,
      read: (model, obj) => {
        const frame = (obj as PluginStateObject.Format.Cif).data.blocks[0];
        const m = frame && matrixOf(model, frame);
        return m ? { key: frame, matrix: m } : "it has no PAE over the model's residues";
      },
    };
  }

  async #fetch(
    ctx: ApiContext,
    it: {
      ref: string;
      label: string;
      id: string;
      modelRef: string;
      model: Model;
      source: (signal: AbortSignal) => Promise<PaeSource>;
    },
  ): Promise<ApiResult<PaeInfo>> {
    const { ref, label, id, modelRef, model } = it;
    const gone = () => ApiResult.fail('not-found', `${label} was removed or replaced while its PAE was fetched.`);
    const stopped = () => ApiResult.fail('aborted', 'Stopped.');
    let source: PaeSource;
    try {
      source = await it.source(ctx.signal);
    } catch (e) {
      if (ctx.signal.aborted) return stopped();
      return ApiResult.fail(
        'failed',
        e instanceof AfdbError ? e.message : `Could not get the PAE of ${id}: ${(e as Error)?.message ?? e}.`,
        'Check the network, and that AlphaFold DB has a PAE for the entry.',
      );
    }
    if (!this.#alive(ref, modelRef, model)) return gone();
    const { url, note } = source;
    const state = ctx.plugin.state.data;
    // The cell is a child of the model: it goes with the structure, and a saved session restores it with the model.
    const b = state.build();
    const cell = source.apply(b, modelRef);
    const remove = async () => {
      if (!state.cells.has(cell)) return;
      await ctx.plugin.runTask(state.updateTree(state.build().delete(cell)));
    };
    // The stop button's signal reaches the download: the task is asked to abort.
    const task = state.updateTree(b);
    const abort = () => ctx.plugin.managers.task.requestAbort(task.id, 'Stopped');
    ctx.signal.addEventListener('abort', abort, { once: true });
    try {
      await ctx.plugin.runTask(task);
    } catch (e) {
      if (!ctx.signal.aborted) throw e;
    } finally {
      ctx.signal.removeEventListener('abort', abort);
    }
    if (ctx.signal.aborted) {
      await remove();
      return stopped();
    }
    if (!this.#alive(ref, modelRef, model)) {
      await remove();
      return gone();
    }
    const parsed = state.cells.get(cell);
    if (parsed?.status === 'error' || !parsed?.obj) {
      const why = parsed?.errorText ?? 'nothing came back';
      await remove();
      return ApiResult.fail('failed', `Could not get the PAE of ${id} from ${url}: ${why}.`, 'Check the network.');
    }
    const read = source.read(model, parsed.obj);
    if (typeof read === 'string') {
      await remove();
      return ApiResult.fail('failed', `The PAE from ${url} does not fit ${label}: ${read}.`);
    }
    const m = read.matrix;
    this.#attachedMatrices.set(read.key, m);
    const chains = [...new Set(paeChains(m).map((c) => c.label))];
    return ApiResult.ok(
      `Fetched the PAE of ${id} from ${source.service}: ${m.size} × ${m.size} residues, up to ${m.max.toFixed(1)} Å.${note ? ` ${note}` : ''}`,
      { nodes: [ref] },
      {
        value: {
          structure: label,
          hasPae: true,
          ...(note ? { note } : {}),
          metric: m.name,
          size: m.size,
          max: m.max,
          chains,
        },
      },
    );
  }

  /**
   * The addresses a fetch reaches before it runs: the associated files named (AlphaFold DB is a known service). A
   * structure named by its state ref gives its own; any other name (a label, which `run` resolves) every loaded
   * structure's, so the agent's policy never sees fewer addresses than the run reaches.
   */
  #addresses(structure: string): string[] {
    const refs = this.#tree().structures.has(structure)
      ? [structure]
      : this.ctx.plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.transform.ref);
    return refs.flatMap((ref) => {
      const file = this.associatedPae(ref);
      return file ? [file.url] : [];
    });
  }

  // Model Archive names a model's PAE file in the model (an associated file, often a zip); AlphaFold DB ships it as a
  // JSON beside the model. Either is a download and a parse in the state tree (the address is the transform's
  // parameter, so a saved session fetches it again), under the model's tag.
  readonly fetchPae = this.operation({
    id: 'model-archive.fetch-pae',
    title: 'Fetch the PAE of a model',
    params: { structure: StructureParam() },
    modules: ['network'],
    approval: 'network',
    // The structure row's "…" menu, for a model whose PAE is elsewhere.
    structureMenu: {
      label: 'Fetch PAE',
      icon: 'download',
      applies: (_plugin, structure) => this.canFetchPae(structure),
    },
    addresses: ({ structure }) => this.#addresses(structure),
    run: async ({ structure }, ctx): Promise<ApiResult<PaeInfo>> => {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const label = labelOf(ctx, s);
      const ref = s.cell.transform.ref;
      if (this.hasPae(ref)) return ApiResult.already(`${label} has a PAE already.`);
      const file = this.associatedPae(ref);
      const entry = file ? undefined : this.afdbEntry(ref);
      if (!file && !entry) {
        return ApiResult.fail(
          'blocked',
          `${label} names no PAE file and is not an AlphaFold DB entry (an AF-<UniProt>-F<n> id), so there is no PAE to fetch.`,
        );
      }
      const modelRef = this.#modelRefOf(ref);
      const model = this.#modelOf(ref);
      if (!modelRef || !model) return ApiResult.fail('not-found', `${label} has no model.`);
      // A second ask while the first is in flight shares it: the cells are added once.
      const running = this.#fetching.get(modelRef);
      if (running) return running;
      const id = entry?.id ?? model.entryId ?? label;
      const source = file
        ? async () => this.#associatedSource(id, file)
        : (signal: AbortSignal) => this.#afdbSource(entry!, signal);
      const fetching = this.#fetch(ctx, { ref, label, id, modelRef, model, source }).finally(() =>
        this.#fetching.delete(modelRef),
      );
      this.#fetching.set(modelRef, fetching);
      return fetching;
    },
  });
}
