/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import { MmcifFormat } from '@molstar/model/formats/structure/mmcif';
import type { ChainIndex, Model, ResidueIndex } from '@molstar/model/model/structure';
import { AtomicHierarchy } from '@molstar/model/model/structure/model/properties/atomic';
import { type ApiAreaContext, type ApiContext, ApiArea, ApiResult } from '@molstar/prism-viewer/api/area';
import { StructureParam } from '@molstar/prism-viewer/api/params';
import { isFailure, labelOf, oneStructure } from '@molstar/prism-viewer/api/scope';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { type Observable, type Subscription, Subject } from 'rxjs';
import { QualityAssessment } from '../quality-assessment/prop.js';
import {
  type AfdbEntry,
  type AfdbPaeSource,
  AfdbError,
  afdbEntryOf,
  paeValuesFromAfdb,
  resolveAfdbPaeUrl,
} from './afdb.js';
import { type PaeMatrix, type PaeResidue, paeChains } from './data.js';
import { FetchPaeJson } from './fetch.js';

// The `model-archive` area (api.md §2, an extension's own area): the predicted aligned error of a model, read
// from the ModelCIF's `ma_qa_metric_local_pairwise` (Model Archive's models carry it; AlphaFold DB ships its PAE apart
// from the model, and `model-archive.fetch-pae` gets it from there: a download and a parse in one cell under the
// model's cell (`fetch.ts`), tagged with the model: it goes with the structure, is undoable and saved in a session). `model-archive.pae` says whether a structure has it and how big it is; the PAE card reads the matrix
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
/** The parsed JSON's own tag: which model's it is (`<tag>:<model cell ref>`). */
const tagOf = (modelRef: string) => `${AfdbPaeTag}:${modelRef}`;

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

/** What one pass over the state tree finds: the structures' models, the fetched PAE cells, the AlphaFold DB entries. */
interface TreeIndex {
  readonly structures: Map<string, { readonly model: Model; readonly modelRef: string }>;
  /** A model's fetched PAE JSON (parsed), by the model's ref. */
  readonly fetched: Map<string, object>;
  /** By structure ref, read when asked. */
  readonly entries: Map<string, AfdbEntry | undefined>;
}

export class ModelArchiveApi extends ApiArea {
  readonly id = 'model-archive';

  /** By model: a matrix is read once (null: the model has none). */
  readonly #matrices = new WeakMap<Model, PaeMatrix | null>();
  /** By model: whether its file carries a PAE metric (the metric list is read once). */
  readonly #carries = new WeakMap<Model, boolean>();
  /** By parsed JSON (an AlphaFold DB PAE in the state tree): its matrix, read once (null: it does not fit). */
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
    const fetched = new Map<string, object>();
    for (const cell of state.cells.values()) {
      const obj = cell.obj;
      if (PluginStateObject.Molecule.Structure.is(obj)) {
        const modelRef = this.#modelRefOf(cell.transform.ref);
        if (modelRef) structures.set(cell.transform.ref, { model: obj.data.model, modelRef });
      } else if (PluginStateObject.Format.Json.is(obj) && cell.transform.tags?.includes(AfdbPaeTag)) {
        const data: unknown = obj.data;
        if (data && typeof data === 'object') fetched.set(cell.transform.parent, data);
      }
    }
    return (this.#index = { structures, fetched, entries: new Map() });
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

  /** The PAE JSON fetched for a structure's model (its parsed cell), if any. */
  #attachedJson(ref: string): { json: object; model: Model } | undefined {
    const t = this.#tree();
    const s = t.structures.get(ref);
    const json = s && t.fetched.get(s.modelRef);
    return s && json ? { json, model: s.model } : undefined;
  }

  #matrix(model: Model): PaeMatrix | undefined {
    let m = this.#matrices.get(model);
    if (m === undefined) {
      m = matrixOf(model) ?? null;
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
    const a = this.#attachedJson(ref);
    if (!a) return undefined;
    let m = this.#attachedMatrices.get(a.json);
    if (m === undefined) {
      const r = attachedMatrix(a.model, a.json);
      if (typeof r === 'string') this.ctx.plugin.log.warn(`AlphaFold DB's PAE does not fit the model: ${r}.`);
      m = typeof r === 'string' ? null : r;
      this.#attachedMatrices.set(a.json, m);
    }
    return m ?? undefined;
  }

  /** Whether a structure (by state ref) has PAE, without reading the matrix: for the card's toolbar button. */
  hasPae(ref: string): boolean {
    return this.#embedded(this.#modelOf(ref)) || !!this.#attachedJson(ref);
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
      const p = cell.transform.params as { url?: unknown; file?: { name?: string }; label?: string } | undefined;
      const url = typeof p?.url === 'string' ? p.url : (p?.url as { url?: string } | undefined)?.url;
      const file = p?.file as { name?: string } | undefined;
      texts.push(cell.obj?.label, url, file?.name, p?.label);
    }
    return afdbEntryOf(texts);
  }

  /** Whether the PAE of a structure can be fetched from AlphaFold DB: it is an entry of it and has no PAE yet. */
  canFetchPae(ref: string): boolean {
    return !this.hasPae(ref) && !!this.afdbEntry(ref);
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
        const afdb = this.afdbEntry(s.cell.transform.ref);
        return ApiResult.ok(
          `${label} has no PAE (a ModelCIF's ma_qa_metric_local_pairwise; AlphaFold DB ships it apart from the model).`,
          {},
          {
            value: { structure: label, hasPae: false },
            ...(afdb ? { next: `model-archive.fetch-pae gets the PAE of ${afdb.id} from AlphaFold DB.` } : {}),
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

  async #fetch(
    ctx: ApiContext,
    it: { ref: string; label: string; entry: AfdbEntry; modelRef: string; model: Model },
  ): Promise<ApiResult<PaeInfo>> {
    const { ref, label, entry, modelRef, model } = it;
    const gone = () => ApiResult.fail('not-found', `${label} was removed or replaced while its PAE was fetched.`);
    const stopped = () => ApiResult.fail('aborted', 'Stopped.');
    let source: AfdbPaeSource;
    try {
      source = await resolveAfdbPaeUrl(entry, undefined, ctx.signal);
    } catch (e) {
      if (ctx.signal.aborted) return stopped();
      return ApiResult.fail(
        'failed',
        e instanceof AfdbError ? e.message : `Could not get the PAE of ${entry.id}: ${(e as Error)?.message ?? e}.`,
        'Check the network, and that AlphaFold DB has a PAE for the entry.',
      );
    }
    if (!this.#alive(ref, modelRef, model)) return gone();
    const { url, note } = source;
    const state = ctx.plugin.state.data;
    // The cell is a child of the model: it goes with the structure, and a saved session restores it with the model.
    const b = state.build();
    const cell = b
      .to(modelRef)
      .apply(FetchPaeJson, { url, label: `${entry.id} PAE` }, { tags: [AfdbPaeTag, tagOf(modelRef)] });
    const remove = async () => {
      if (!state.cells.has(cell.ref)) return;
      await ctx.plugin.runTask(state.updateTree(state.build().delete(cell.ref)));
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
    const parsed = state.cells.get(cell.ref);
    if (parsed?.status === 'error' || !parsed?.obj) {
      const why = parsed?.errorText ?? 'nothing came back';
      await remove();
      return ApiResult.fail(
        'failed',
        `Could not get the PAE of ${entry.id} from ${url}: ${why}.`,
        'Check the network.',
      );
    }
    const m = attachedMatrix(model, parsed.obj.data);
    if (typeof m === 'string') {
      await remove();
      return ApiResult.fail('failed', `The PAE from ${url} does not fit ${label}: ${m}.`);
    }
    this.#attachedMatrices.set(parsed.obj.data as object, m);
    const chains = [...new Set(paeChains(m).map((c) => c.label))];
    return ApiResult.ok(
      `Fetched the PAE of ${entry.id} from AlphaFold DB: ${m.size} × ${m.size} residues, up to ${m.max.toFixed(1)} Å.${note ? ` ${note}` : ''}`,
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

  // AlphaFold DB ships a model's PAE as a JSON beside it: a download and a parse in the state tree (the address is the
  // download's parameter, so a saved session fetches it again), under the model's tag.
  readonly fetchPae = this.operation({
    id: 'model-archive.fetch-pae',
    title: 'Fetch the PAE of a model from AlphaFold DB',
    params: { structure: StructureParam() },
    modules: ['network'],
    approval: 'network',
    // The structure row's "…" menu, for an AlphaFold DB entry without a PAE.
    structureMenu: {
      label: 'Fetch PAE from AlphaFold DB',
      icon: 'download',
      applies: (_plugin, structure) => this.canFetchPae(structure),
    },
    // AlphaFold DB, a known data service.
    addresses: () => [],
    run: async ({ structure }, ctx): Promise<ApiResult<PaeInfo>> => {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const label = labelOf(ctx, s);
      const ref = s.cell.transform.ref;
      const entry = this.afdbEntry(ref);
      if (!entry) {
        return ApiResult.fail(
          'blocked',
          `${label} is not an AlphaFold DB entry (an AF-<UniProt>-F<n> id), so there is no PAE to fetch.`,
        );
      }
      if (this.hasPae(ref)) return ApiResult.already(`${label} has a PAE already.`);
      const modelRef = this.#modelRefOf(ref);
      const model = this.#modelOf(ref);
      if (!modelRef || !model) return ApiResult.fail('not-found', `${label} has no model.`);
      // A second ask while the first is in flight shares it: the cells are added once.
      const running = this.#fetching.get(modelRef);
      if (running) return running;
      const fetching = this.#fetch(ctx, { ref, label, entry, modelRef, model }).finally(() =>
        this.#fetching.delete(modelRef),
      );
      this.#fetching.set(modelRef, fetching);
      return fetching;
    },
  });
}
