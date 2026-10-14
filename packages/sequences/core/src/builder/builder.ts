/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder (seqspec:format#builders): one builder per place in the grammar, typed from the parameter tables, as
// MVS's is.
// Data methods return handles; a `ref<kind>` param takes a handle or a ref. Sugar writes ordinary nodes: sources,
// main sections, single-layer tracks, a layer's own annotation, row bindings, offsets and gapped pairs as mappings,
// fallbacks of extension types and themes; `getDocument()` generates the refs nobody gave.

import type {
  DocNode,
  KernelLayerParams,
  KernelLayerTypeName,
  LayerCommonParams,
  Metadata,
  SeqViewSpecDoc,
  Snapshot,
  SvsRegistries,
  SvsTreeSchema,
  ThemeRef,
  ValuesFor,
} from "@molstar/seq-core/format";
import { Colors } from "./colors.js";
import {
  cigarFromRows,
  f32,
  isNamespaced,
  notNucleotide,
  SvsFormat,
  SvsVersion,
  translateDna,
  wussPairs,
} from "./helpers.js";

export type BuilderObject = Record<string, unknown>;
export interface BuilderNode {
  kind: string;
  ref?: string;
  params: BuilderObject;
  custom?: BuilderObject;
  children: BuilderNode[];
  parent: BuilderNode | null;
  /** The view a node sits under: its generated ref is `<view ref>.<slug>`. */
  view: BuilderNode | null;
  /** A ref made from another node's: `<ref of node><suffix>` (`.src`, `.main`, `.row-3`). */
  derived?: { node: BuilderNode; suffix: string };
}

/** A node, or a derived ref under it (`R/<suffix>`); stands for its ref in a param that takes its kind (`K`). */
export class NodeHandle<K extends string = string> {
  /** Phantom: the node kind, so that a `ref<kind>` param takes only a handle of one of its kinds. */
  declare readonly __kind?: K;
  readonly node: BuilderNode;
  readonly suffix: string | null;
  protected readonly st: BuilderState;
  constructor(st: BuilderState, node: BuilderNode, suffix: string | null = null) {
    this.st = st;
    this.node = node;
    this.suffix = suffix;
  }
  /** The ref; generated refs are known after `getDocument()`. */
  get ref(): string {
    return this.suffix === null ? (this.node.ref ?? "") : `${this.node.ref}/${this.suffix}`;
  }
}

export type BuilderSchemas = (typeof SvsTreeSchema)["nodes"];
/** A param as the builder takes it: only a `ref<kind>` (a `RefValue`) takes a handle, of one of its kinds.
 *  The mapped type is homomorphic, so lists stay lists and tuples (`domain`, `cds`) keep their length. */
export type BuilderParam<T> = T extends string
  ? T extends { readonly __refKinds?: infer K extends string }
    ? string | NodeHandle<K>
    : T
  : T extends object
    ? { [K in keyof T]: BuilderParam<T[K]> }
    : T;
/** The params a schema accepts, with handles where it has refs. */
export type BuilderValuesFor<P> = BuilderParam<ValuesFor<P>>;
export type BuilderCommon = { ref?: string; custom?: BuilderObject };
export type BuilderParamsOf<K extends keyof BuilderSchemas> = BuilderValuesFor<BuilderSchemas[K]["params"]> &
  BuilderCommon;
export type BuilderSourceParams = BuilderParamsOf<"source">;
export type WithSourceParams<T> = T extends unknown
  ? Omit<T, "source_ref"> & { source_ref?: string | NodeHandle<"source"> | null; source?: BuilderSourceParams }
  : never;
export type WithoutValue<T> = T extends unknown ? Omit<T, "value"> : never;
/** `Omit` over each case of a union. */
export type DistributiveOmit<T, K extends PropertyKey> = T extends unknown ? Omit<T, K> : never;
/** What a mapping's `from` and `to` name. */
export type MappingSide = string | NodeHandle<"sequence" | "ligand" | "unit" | "row">;
/** A layer's data: an annotation or computed node, by handle or ref. */
type Data = string | NodeHandle<"annotation" | "computed">;
export type AnnotationParams = WithSourceParams<WithoutValue<BuilderParamsOf<"annotation">>> & {
  value?: unknown;
  columnar?: boolean;
};
export type MappingParams =
  | BuilderParamsOf<"mapping">
  | (BuilderCommon & { from: MappingSide; to: MappingSide; method: "offset"; offset: number })
  | (BuilderCommon & { from: MappingSide; to: MappingSide; method: "pairwise"; gapped: [string, string] });
/**
 * `split` (seqspec:format#layer-types): a field's name is `{ by, order: null }`; `"alt"` on `variants`, UniProt's grid
 * (shown whole).
 */
export type LayerSplit = string | { by: string; order?: string[] | null; empty?: "show" | "hide" } | null;
/** A layer's own annotation, and the common params of any layer type. */
export type LayerBaseParams = BuilderCommon & { annotation?: AnnotationParams } & BuilderParam<LayerCommonParams>;
/** A kernel type's params: the common ones and the type's own (`letters`' `of`: one polymer, a list, or `"mapped"`). */
export type LayerParams<T extends KernelLayerTypeName = KernelLayerTypeName> = LayerBaseParams &
  ("split" extends keyof KernelLayerParams[T]
    ? Omit<BuilderParam<KernelLayerParams[T]>, "split"> & { split?: LayerSplit }
    : BuilderParam<KernelLayerParams[T]>);
/** An extension layer type (`<namespace>:<name>`), its `options` checked by its registered table. */
export type ExtensionLayerParams = LayerBaseParams & { type: `${string}:${string}` };
/** UniProt's variant rows (Nightingale's residue order), then a deletion and a stop. */
export const VARIANT_ROWS = [..."GAVLISTCMDNEQRKHFYWP", "-", "*"];
/** A kernel provider's params, or an extension's `<namespace>:<name>` with its `options`. */
export type ComputedParams =
  | BuilderParamsOf<"computed">
  | (BuilderCommon & { provider: `${string}:${string}`; options?: BuilderParam<BuilderObject> });
/** A mapping's params but its ends: `from` is the row, `to` the third argument. */
export type BindParams = DistributiveOmit<MappingParams, "from" | "to">;
/** Twin: `KERNEL_LAYER_TYPE_NAMES` in the format's `tables/kernel.ts`; a test keeps them equal. */
const LAYER_TYPES =
  "letters rows ranges points emphasis variants secondary values color ruler pairs logo matrix dot_plot".split(
    " ",
  ) as KernelLayerTypeName[];

export interface BuilderState {
  all: BuilderNode[];
  root: BuilderNode;
  data: BuilderNode;
  sources: Map<string, BuilderNode>;
  /** Each distinct inline `value` of a source, numbered (`sourceKey`). */
  values: Map<unknown, number>;
  snapshots: BuilderObject[];
  registries: SvsRegistries | null;
  /** `reverseComplement`'s sections and maps, written on their own letters layers when the document is made. */
  complements: Map<BuilderNode, string>;
}

function mk(
  st: BuilderState,
  kind: string,
  params: BuilderObject,
  parent: BuilderNode,
  view: BuilderNode | null,
): BuilderNode {
  const { ref, custom, ...rest } = params;
  const n: BuilderNode = { kind, params: rest, children: [], parent, view };
  if (typeof ref === "string") n.ref = ref;
  if (custom) n.custom = custom as BuilderObject;
  parent.children.push(n);
  st.all.push(n);
  return n;
}

/**
 * A value as JSON with its keys sorted at every level (the format's `canonical`; the builder's entry loads only its own
 * files).
 */
function sortedJson(x: unknown): string {
  if (x === null || typeof x !== "object") return JSON.stringify(x) ?? "null";
  if (ArrayBuffer.isView(x)) return sortedJson(Array.from(x as unknown as ArrayLike<number>));
  if (Array.isArray(x)) return `[${x.map(sortedJson).join(",")}]`;
  const keys = Object.keys(x).filter((k) => (x as BuilderObject)[k] !== undefined);
  return `{${keys
    .sort()
    .map((k) => `${JSON.stringify(k)}:${sortedJson((x as BuilderObject)[k])}`)
    .join(",")}}`;
}

/**
 * A source's params as a key (CB-08): every field in name order, nested values too (`sortedJson`, LA-41), so
 * `{uri, format}` and `{format, uri}` are one source, a field left undefined is no field, and an inline `value` is by
 * itself (a number per distinct value), so the key holds no second copy of a file.
 */
function sourceKey(st: BuilderState, p: BuilderObject): string {
  return Object.keys(p)
    .filter((k) => p[k] !== undefined)
    .sort()
    .map((k) => {
      if (k !== "value") return `${k}=${sortedJson(p[k])}`;
      let id = st.values.get(p[k]);
      if (id === undefined) {
        id = st.values.size;
        st.values.set(p[k], id);
      }
      return `value#${id}`;
    })
    .join("\u0000");
}

/** `source: {...}` as a source node `<consumer ref>.src`, placed before the consumer; equal params share one. */
function withSource(st: BuilderState, n: BuilderNode) {
  const p = n.params.source as BuilderObject | undefined;
  if (!p) return n;
  delete n.params.source;
  const key = sourceKey(st, p);
  let src = st.sources.get(key);
  if (!src) {
    src = mk(st, "source", p, st.data, null);
    src.derived = { node: n, suffix: ".src" };
    st.data.children.pop();
    let top: BuilderNode = n;
    while (top.parent && top.parent !== st.data) top = top.parent;
    const at = st.data.children.indexOf(top);
    st.data.children.splice(at < 0 ? st.data.children.length : at, 0, src);
    st.sources.set(key, src);
  }
  n.params.source_ref = new NodeHandle(st, src);
  return n;
}

const KINDS = [Int8Array, Uint8Array, Int16Array, Uint16Array, Int32Array, Uint32Array, Float32Array, Float64Array];
const isObj = (v: unknown): v is BuilderObject =>
  !!v && typeof v === "object" && !Array.isArray(v) && !ArrayBuffer.isView(v);

/**
 * A table column given as a typed array (seqspec:format#values): one of its kinds, kept as it is (`resolveColumns`).
 */
function checkColumn(c: unknown, ref: unknown, k: string) {
  if (!ArrayBuffer.isView(c) || KINDS.some((Kind) => c instanceof Kind)) return;
  const at = `annotation${ref === undefined ? "" : ` "${ref}"`}: column "${k}"`;
  throw new TypeError(`${at} is a ${c.constructor.name}; give an array, or Int8-32, Uint8-32, Float32 or Float64.`);
}

/** `float32`: the columns of numbers with a fraction as `Float32Array`. */
function float32Columns(value: BuilderObject) {
  for (const [k, c] of Object.entries(value)) {
    const xs: unknown[] = c instanceof Float64Array ? Array.from(c) : Array.isArray(c) ? c : [];
    if (xs.some((x) => !Number.isInteger(x ?? 0)) && xs.every((x) => x === null || typeof x === "number"))
      value[k] = Float32Array.from(xs, (x) => (x as number | null) ?? Number.NaN);
  }
}

/** An annotation's `value`: typed columns (and `codes`) checked, kept (`resolveColumns`); `columnar` writes rows as
 *  columns, repeated strings coded; `float32` makes its columns of fractions Float32. */
function annotationValue(p: BuilderObject, ref: unknown) {
  const { columnar, float32, value } = p;
  delete p.columnar;
  delete p.float32;
  if (isObj(value)) for (const [k, c] of Object.entries(value)) checkColumn(isObj(c) ? c.codes : c, ref, k);
  if (columnar && Array.isArray(value)) rowsToColumns(p, value);
  if (float32 && isObj(p.value)) float32Columns(Object.assign(p, { value: { ...p.value } }).value);
}

function rowsToColumns(p: BuilderObject, value: unknown[]) {
  const keys = new Set<string>();
  for (const r of value) for (const k of Object.keys(r as BuilderObject)) keys.add(k);
  p.value = Object.fromEntries(
    [...keys].map((k) => [k, dictionaryColumn(value.map((r) => (r as BuilderObject)[k] ?? null))]),
  );
}

/** One pass: a column of repeated strings (and nulls) as `{ dict, codes }`; all-null, mixed and unique columns plain. */
function dictionaryColumn(col: unknown[]) {
  const codes: (number | null)[] = new Array(col.length);
  const index = new Map<string, number>();
  let strings = 0;
  for (let i = 0; i < col.length; i++) {
    const v = col[i];
    if (v === null) {
      codes[i] = null;
      continue;
    }
    if (typeof v !== "string") return col;
    let code = index.get(v);
    if (code === undefined) {
      code = index.size;
      index.set(v, code);
    }
    codes[i] = code;
    strings++;
  }
  return index.size === strings ? col : { dict: [...index.keys()], codes };
}

export class MoleculeHandle<
  K extends "sequence" | "ligand" | "unit" = "sequence" | "ligand" | "unit",
> extends NodeHandle<K> {
  /** Adds a numbering child. */
  numbering(params: BuilderParamsOf<"numbering">): this {
    mk(this.st, "numbering", params as BuilderObject, this.node, null);
    return this;
  }
}

export class ComplexHandle extends NodeHandle<"complex"> {
  /** A provider's unit `R/<suffix>`; or, given params, a `unit` child (a complex by value). */
  unit(suffix: string): NodeHandle<"unit">;
  unit(params: BuilderParamsOf<"unit">): MoleculeHandle<"unit">;
  unit(x: string | BuilderParamsOf<"unit">) {
    if (typeof x === "string") return new NodeHandle<"unit">(this.st, this.node, x);
    return new MoleculeHandle<"unit">(this.st, mk(this.st, "unit", x as BuilderObject, this.node, null));
  }
  sequence(suffix: string): NodeHandle<"sequence"> {
    return new NodeHandle<"sequence">(this.st, this.node, suffix);
  }
  ligand(suffix: string): NodeHandle<"ligand"> {
    return new NodeHandle<"ligand">(this.st, this.node, suffix);
  }
}

export class AlignmentHandle extends NodeHandle<"alignment"> {
  tree(params: WithSourceParams<BuilderParamsOf<"tree">>): this {
    withSource(this.st, mk(this.st, "tree", params as BuilderObject, this.node, null));
    return this;
  }
  row(k: number): NodeHandle<"row"> {
    return new NodeHandle<"row">(this.st, this.node, String(k));
  }
  /** A mapping from row `k` to `to`, ref `<alignment ref>.row-<k>` unless given. */
  bind(k: number, to: MappingSide, params: BindParams): NodeHandle<"mapping"> {
    const n = mapping(this.st, { ...params, from: this.row(k), to } as MappingParams);
    if (params.ref === undefined) n.derived = { node: this.node, suffix: `.row-${k}` };
    return new NodeHandle<"mapping">(this.st, n);
  }
}

function mapping(st: BuilderState, params: MappingParams): BuilderNode {
  const p = { ...params } as BuilderObject;
  if (p.method === "offset") {
    const o = p.offset as number;
    p.method = "segments";
    p.segments = [o >= 0 ? [0, o, null] : [-o, 0, null]];
    delete p.offset;
  } else if (p.method === "pairwise" && p.gapped) {
    const [a, b] = (p.gapped as string[]).map((s) => s.replaceAll(".", "-"));
    p.cigar = cigarFromRows(a, b);
    delete p.gapped;
  }
  return mk(st, "mapping", p, st.data, null);
}

function dataBuilder(st: BuilderState) {
  const node = (kind: string, params: object) => withSource(st, mk(st, kind, { ...params }, st.data, null));
  return {
    source: (params: BuilderSourceParams) =>
      new NodeHandle<"source">(st, mk(st, "source", { ...params }, st.data, null)),
    sequence: (params: WithSourceParams<BuilderParamsOf<"sequence">>) =>
      new MoleculeHandle<"sequence">(st, node("sequence", params)),
    ligand: (params: WithSourceParams<BuilderParamsOf<"ligand">>) =>
      new MoleculeHandle<"ligand">(st, node("ligand", params)),
    unit: (params: BuilderParamsOf<"unit">) => new MoleculeHandle<"unit">(st, node("unit", params)),
    complex: (params: WithSourceParams<BuilderParamsOf<"complex">>) => new ComplexHandle(st, node("complex", params)),
    alignment: (params: WithSourceParams<BuilderParamsOf<"alignment">>) =>
      new AlignmentHandle(st, node("alignment", params)),
    mapping: (params: MappingParams) => new NodeHandle<"mapping">(st, mapping(st, params)),
    annotation(params: AnnotationParams) {
      const n = node("annotation", params);
      annotationValue(n.params, n.ref);
      return new NodeHandle<"annotation">(st, n);
    },
    computed: (params: ComputedParams) => new NodeHandle<"computed">(st, node("computed", params)),
  };
}
export type DataBuilder = ReturnType<typeof dataBuilder>;

/** An extension theme's fallback from the registered table, when not given. */
function themeFallback(st: BuilderState, t: unknown) {
  const r = t as ThemeRef | null | undefined;
  if (!r || !isNamespaced(r.theme) || r.fallback !== undefined) return r;
  return { ...r, fallback: st.registries?.get("themes", r.theme)?.fallback ?? null };
}

/** A method per kernel layer type, its params that type's; `layer` takes any type, kernel or `<namespace>:<name>`. */
export type TrackBuilder = { [T in KernelLayerTypeName]: (params?: LayerParams<T>) => TrackBuilder } & {
  layer(
    params: { [T in KernelLayerTypeName]: LayerParams<T> & { type: T } }[KernelLayerTypeName] | ExtensionLayerParams,
  ): TrackBuilder;
};
function trackBuilder(st: BuilderState, track: BuilderNode, view: BuilderNode): TrackBuilder {
  const t = {
    layer(params: BuilderObject & { type: string; annotation?: AnnotationParams }) {
      const { annotation, ...rest } = params;
      const n = mk(st, "layer", rest, track, view);
      const p = n.params;
      if (isNamespaced(p.type) && p.fallback == null) {
        p.fallback = st.registries?.get("layerTypes", p.type)?.fallback ?? undefined;
      }
      if (typeof p.split === "string")
        p.split = { by: p.split, order: p.type === "variants" && p.split === "alt" ? [...VARIANT_ROWS] : null };
      if (p.color) p.color = themeFallback(st, p.color);
      if (p.glyph) p.glyph = themeFallback(st, p.glyph);
      if (annotation) {
        const a = withSource(st, mk(st, "annotation", { ...annotation }, n, view));
        annotationValue(a.params, a.ref);
      }
      return t as TrackBuilder;
    },
  } as TrackBuilder;
  for (const type of LAYER_TYPES) {
    (t as Record<string, unknown>)[type] = (params: BuilderObject = {}) => t.layer({ ...params, type } as never);
  }
  return t;
}

function groupBuilder(st: BuilderState, parent: BuilderNode, view: BuilderNode) {
  const g = {
    group: (params: BuilderParamsOf<"group">) =>
      groupBuilder(st, mk(st, "group", params as BuilderObject, parent, view), view),
    track: (params: BuilderParamsOf<"track"> = {}) =>
      trackBuilder(st, mk(st, "track", params as BuilderObject, parent, view), view),
  };
  return g;
}

function sectionBuilder(st: BuilderState, section: BuilderNode, view: BuilderNode) {
  const s = {
    ...groupBuilder(st, section, view),
    /** The section, for a param that names one (`follow`). */
    handle: new NodeHandle<"section">(st, section),
    block(params: BuilderParamsOf<"block">) {
      mk(st, "block", params as BuilderObject, section, view);
      return s;
    },
    focus(params: BuilderParamsOf<"focus">) {
      mk(st, "focus", params as BuilderObject, section, view);
      return s;
    },
  };
  return s;
}
export type SectionBuilder = ReturnType<typeof sectionBuilder>;

/** A view; the section methods go to its main section `<view ref>.main`, made on first use or for section params. */
function viewBuilder(st: BuilderState, params: BuilderParamsOf<"view"> & Partial<BuilderParamsOf<"section">>) {
  const { ref, label, layout, custom, ...sec } = params as BuilderObject;
  const v = mk(st, "view", { ref, label, layout, custom }, st.root, null);
  if (layout === undefined) delete v.params.layout;
  let main: SectionBuilder | null = null;
  const section = (p: BuilderParamsOf<"section">) => sectionBuilder(st, mk(st, "section", p as BuilderObject, v, v), v);
  const first = () => {
    if (main) return main;
    main = section(sec as BuilderParamsOf<"section">);
    v.children[v.children.length - 1].derived = { node: v, suffix: ".main" };
    return main;
  };
  if (Object.keys(sec).length > 0) first();
  return {
    section,
    /** The main section, for a param that names one (`follow`) or a preset (made on first use). */
    get handle() {
      return first().handle;
    },
    block: (p: BuilderParamsOf<"block">) => first().block(p),
    group: (p: BuilderParamsOf<"group">) => first().group(p),
    track: (p?: BuilderParamsOf<"track">) => first().track(p),
    focus: (p: BuilderParamsOf<"focus">) => first().focus(p),
  };
}
export type ViewBuilder = ReturnType<typeof viewBuilder>;

const slug = (s: unknown) =>
  typeof s === "string"
    ? s
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, "-")
        .replace(/^-+|-+$/g, "")
    : "";

/** Generated refs, in building order; explicit refs win and are avoided. */
function generateRefs(st: BuilderState) {
  const used = new Set(st.all.map((n) => n.ref).filter((r) => r !== undefined));
  const free = (base: string) => {
    let ref = base;
    for (let k = 2; used.has(ref); k++) ref = `${base}-${k}`;
    used.add(ref);
    return ref;
  };
  for (const n of st.all) {
    if (n.ref !== undefined) continue;
    if (n.derived) {
      n.ref = free(`${n.derived.node.ref}${n.derived.suffix}`);
      continue;
    }
    const layerType = (l?: BuilderNode) => slug(String(l?.params.type ?? ""));
    const s =
      slug(n.params.label) ||
      (n.kind === "track" ? layerType(n.children[0]) : n.kind === "layer" ? layerType(n) : "") ||
      n.kind;
    n.ref = free(n.view ? `${n.view.ref}.${s}` : s);
  }
}

/** Handles as refs, typed arrays as arrays (float32 as its shortest decimals; an annotation's `value` table aside:
 *  `resolveColumns`), NaN as null. */
function resolve(v: unknown, annotation = false): unknown {
  if (v instanceof NodeHandle) return v.ref;
  if (typeof v === "number") return Number.isNaN(v) ? null : v;
  if (ArrayBuffer.isView(v))
    return Array.from(v as Float64Array, (x) => resolve(v instanceof Float32Array ? f32(x) : x));
  if (Array.isArray(v)) return v.map((x) => resolve(x));
  if (v !== null && typeof v === "object") {
    return Object.fromEntries(
      Object.entries(v)
        .filter(([, x]) => x !== undefined)
        .map(([k, x]) => [k, annotation && k === "value" ? resolveColumns(x) : resolve(x)]),
    );
  }
  return v;
}

/**
 * An annotation's table: typed columns, typed `codes` and a list column's typed cells (a grid's rows) kept as they are
 * (seqspec:format#values, seqspec:format#annotation-tables; a writer makes them JSON arrays: `toJSON`), the rest resolved.
 */
function resolveColumns(value: unknown) {
  if (!isObj(value)) return resolve(value);
  const codes = (c: BuilderObject) =>
    ArrayBuffer.isView(c.codes) ? { ...(resolve(c) as BuilderObject), codes: c.codes } : resolve(c);
  const cells = (c: unknown[]) => c.map((x) => (ArrayBuffer.isView(x) ? x : resolve(x)));
  const col = (c: unknown) =>
    ArrayBuffer.isView(c) ? c : Array.isArray(c) ? cells(c) : isObj(c) ? codes(c) : resolve(c);
  return Object.fromEntries(Object.entries(value).flatMap(([k, c]) => (c === undefined ? [] : [[k, col(c)]])));
}

function toDoc(n: BuilderNode): DocNode {
  const out: DocNode = { kind: n.kind };
  if (n.kind !== "root" && n.kind !== "data") out.ref = n.ref;
  const params = resolve(n.params, n.kind === "annotation") as BuilderObject;
  if (Object.keys(params).length > 0) out.params = params;
  if (n.custom) out.custom = resolve(n.custom) as BuilderObject;
  if (n.children.length > 0) out.children = n.children.map(toDoc);
  return out;
}

/** The layers under a node, in order. */
const layersUnder = (n: BuilderNode): BuilderNode[] =>
  n.children.flatMap((c) => (c.kind === "layer" ? [c] : layersUnder(c)));

/** A document builder. `registries` give extension layer types and themes their fallbacks. */
export function createSeqViewSpec(
  metadata: Metadata & { revision?: string } = {},
  opts: { registries?: SvsRegistries } = {},
) {
  const root: BuilderNode = { kind: "root", params: {}, children: [], parent: null, view: null };
  const data: BuilderNode = { kind: "data", params: {}, children: [], parent: root, view: null };
  root.children.push(data);
  const st: BuilderState = {
    all: [],
    root,
    data,
    sources: new Map(),
    values: new Map(),
    snapshots: [],
    registries: opts.registries ?? null,
    complements: new Map(),
  };
  const { revision, ...meta } = metadata;
  const b = {
    data: dataBuilder(st),
    view: (params: BuilderParamsOf<"view"> & Partial<BuilderParamsOf<"section">>) => viewBuilder(st, params),
    /** Application-level compositions of generic params. */
    preset: {
      /**
       * Several regions on one axis (seqspec:model#pieces-and-layout-maps): a block per `[unit or sequence, start,
       * stop]`, in order, each the piece `range: [start, stop]` of its unit; `direction` on each when given.
       */
      multiRegion<S extends Pick<SectionBuilder, "block">>(
        section: S,
        regions: readonly (readonly [string | NodeHandle<"unit" | "sequence">, number, number])[],
        o: { direction?: 1 | -1 } = {},
      ): S {
        for (const [of, start, stop] of regions) {
          const key = of instanceof NodeHandle && of.node.kind === "sequence" ? "sequence" : "unit";
          section.block({ [key]: of, range: [start, stop], ...(o.direction ? { direction: o.direction } : {}) });
        }
        return section;
      },
      /**
       * Only a section's exons: `positions` over the exon annotation (or a list of them), so each block shows its
       * covered positions, the gaps between them as breaks. Call it on a section from `view.section`.
       */
      exonsOnly<S extends Pick<SectionBuilder, "handle">>(
        section: S,
        exons:
          | string
          | NodeHandle<"annotation" | "computed">
          | readonly (string | NodeHandle<"annotation" | "computed">)[],
      ): S {
        section.handle.node.params.positions = exons;
        return section;
      },
      /**
       * The section's other strand (seqspec:format#letter-maps): `direction: -1` and a complement `letter_map` on the
       * layers that show the blocks' own letters (`letters` and `rows`, without `of`), unless a layer has its own;
       * those layers are the section's when the document is made, so the call may come before its tracks. The map
       * lists the alphabets it is for, so a protein's letters stay as stored; `map` names another registered table
       * (`complement_rna`).
       */
      reverseComplement<S extends Pick<SectionBuilder, "handle">>(section: S, o: { map?: string } = {}): S {
        section.handle.node.params.direction = -1;
        st.complements.set(section.handle.node, o.map ?? "complement");
        return section;
      },
      /**
       * A plasmid map (seqspec:format#keyed-sets, CG2): a circular section of `molecule` (its ruler outermost;
       * `features` as ranges labeled outside the ring with an arrowhead by `direction`'s field; `sites` as labeled
       * ticks on the features' ring, a second layer of their track; `pairs` as chords) and, with `detail` (default), a
       * linear section of the same tracks with its letters, following the ring (whose window shows its view). Only
       * generic params: `style: "circular"`, `labels`, `direction_field`, `display`, size and color themes, `follow`.
       */
      plasmidMap(
        view: Pick<ViewBuilder, "section">,
        molecule: BuilderParamsOf<"block">["sequence"],
        o: {
          ref?: string;
          label?: string;
          features?: Data;
          sites?: Data;
          pairs?: Data;
          direction?: string | null;
          detail?: boolean;
        } = {},
      ): { ring: SectionBuilder; detail: SectionBuilder | null } {
        const direction = o.direction === undefined ? "strand" : o.direction;
        // On the ring the pairs ring is thin and innermost: the chords leave from just inside the features.
        const tracks = (sec: SectionBuilder, linear: boolean) => {
          sec.track({ label: "" }).ruler();
          if (linear) sec.track({ label: "Sequence" }).letters();
          const features = o.features
            ? sec.track({ label: o.sites && !linear ? "Features and sites" : "Features", height: linear ? 2 : 1.4 })
            : null;
          features?.ranges({
            data: o.features,
            labels: true,
            ...(direction ? { direction_field: direction } : {}),
          } as LayerParams<"ranges">);
          // On the ring the sites are a second layer of the features' track (CG3): ticks across its ring, over the
          // features, in one lane, in no data color (a rules theme with no rule: the renderer's faint skin ink, so
          // the skin decides it); linearly a track of their own.
          if (o.sites && features && !linear)
            features.ranges({
              data: o.sites,
              labels: true,
              display: "union",
              radius: { theme: "uniform", params: { size: 12 } },
              color: { theme: "rules", params: { rules: [] } },
            } as LayerParams<"ranges">);
          else if (o.sites)
            sec
              .track({ label: "Sites", height: linear ? 1 : 0.6 })
              .ranges({ data: o.sites, labels: true } as LayerParams<"ranges">);
          if (o.pairs)
            sec.track({ label: "Pairs", height: linear ? 3 : 0.4 }).pairs({ data: o.pairs, display: "both" });
          return sec;
        };
        const label = o.label ?? "Plasmid map";
        const ring = tracks(
          view.section({ ref: o.ref, label, style: "circular" }).block({ sequence: molecule }),
          false,
        );
        const detail =
          o.detail === false
            ? null
            : tracks(
                view.section({ label: `${label}, linear`, follow: ring.handle }).block({ sequence: molecule }),
                true,
              );
        return { ring, detail };
      },
      /**
       * A circular genome (CG2): a circular section of `molecule`, its ruler, its `features` on two rings by
       * `direction`'s field (`+` outside, `-` inside, arrowheads; labeled only with `labels`), each of `values` as bars
       * on a ring of its own; and with `detail` (default) a linear section following it.
       */
      circularGenome(
        view: Pick<ViewBuilder, "section">,
        molecule: BuilderParamsOf<"block">["sequence"],
        o: {
          ref?: string;
          label?: string;
          features?: Data;
          values?: readonly Data[];
          direction?: string;
          labels?: boolean;
          detail?: boolean;
        } = {},
      ): { ring: SectionBuilder; detail: SectionBuilder | null } {
        const direction = o.direction ?? "strand";
        const tracks = (sec: SectionBuilder) => {
          sec.track({ label: "" }).ruler();
          for (const [strand, name] of [
            ["+", "Forward"],
            ["-", "Reverse"],
          ] as const)
            if (o.features)
              sec.track({ label: `${name} strand` }).ranges({
                data: o.features,
                labels: o.labels ?? false,
                direction_field: direction,
                filter: { field: direction, values: [strand] },
              } as LayerParams<"ranges">);
          for (const values of o.values ?? []) sec.track({ height: 1.5 }).values({ data: values });
          return sec;
        };
        const label = o.label ?? "Genome";
        const ring = tracks(view.section({ ref: o.ref, label, style: "circular" }).block({ sequence: molecule }));
        const detail =
          o.detail === false
            ? null
            : tracks(view.section({ label: `${label}, linear`, follow: ring.handle }).block({ sequence: molecule }));
        return { ring, detail };
      },
      /**
       * A sashimi plot (seqspec:format#themes): splice junctions, a `pairs` annotation with a read count per row, as
       * arcs stroked by the count (`width` a `field_scale` into `range` px, 1 to 8) and labeled with it, a track per
       * strand placed one over the other: rows whose `strand` field holds `+` rise from the line between them (`side:
       * "above"`), those with `-` hang from it (`below`); `strand: null` one track above. Only generic params: `side`,
       * `filter`, `width`, `labels`. Tracks of `height` rows each (3), labeled `label` and the strand.
       */
      sashimi<S extends Pick<SectionBuilder, "track">>(
        section: S,
        junctions: string | NodeHandle<"annotation" | "computed">,
        o: {
          label?: string;
          count?: string;
          strand?: string | null;
          range?: readonly [number, number];
          height?: number;
          color?: ThemeRef;
        } = {},
      ): S {
        const width = { theme: "field_scale", params: { field: o.count ?? "count", range: [...(o.range ?? [1, 8])] } };
        const label = o.label ?? "Junctions";
        const strand = o.strand === undefined ? "strand" : o.strand;
        const sides =
          strand === null
            ? ([["above", null]] as const)
            : ([
                ["above", "+"],
                ["below", "-"],
              ] as const);
        for (const [side, value] of sides)
          section.track({ label: value === null ? label : `${label} (${value})`, height: o.height ?? 3 }).pairs({
            data: junctions,
            display: "arcs",
            side,
            labels: true,
            width,
            ...(o.color ? { color: o.color } : {}),
            ...(strand !== null && value !== null ? { filter: { field: strand, values: [value] } } : {}),
          } as LayerParams<"pairs">);
        return section;
      },
      /**
       * A protein alignment over its coding DNA's (seqspec:format#sections): a `codon` mapping from each protein row to
       * its cDNA row (row k to row k, or the `rows` pairs), and a section of each in `view`, the cDNA's following the
       * protein's at a third of its cell. Rows are counted from a `value`; for an alignment read from a source give
       * `rows`.
       */
      codingPair(
        view: Pick<ViewBuilder, "section">,
        protein: AlignmentHandle,
        cdna: AlignmentHandle,
        o: { rows?: readonly (readonly [number, number])[]; labels?: readonly [string, string] } = {},
      ) {
        const value = protein.node.params.value;
        const rows = o.rows ?? (Array.isArray(value) ? value.map((_, k) => [k, k] as const) : null);
        if (!rows)
          throw new Error("codingPair: give rows, [protein row, cDNA row] pairs, for an alignment from a source");
        for (const [p, c] of rows) protein.bind(p, cdna.row(c), { method: "codon" });
        const [pl, cl] = o.labels ?? ["Protein", "Coding DNA"];
        const top = view.section({ label: pl }).block({ alignment: protein });
        const bottom = view.section({ label: cl, follow: top.handle }).block({ alignment: cdna });
        return { protein: top, cdna: bottom };
      },
    },
    /** Color schemes as `color` theme refs: `b.colors.clustalX`, `b.colors.viridis({ domain })`. */
    colors: Colors,
    snapshot(params: Snapshot) {
      st.snapshots.push(params as unknown as BuilderObject);
      return b;
    },
    /** The document; refs not given are generated now. */
    getDocument(more: Metadata = {}): SeqViewSpecDoc {
      generateRefs(st);
      for (const [node, map] of st.complements)
        for (const l of layersUnder(node)) {
          const own = (l.params.type === "letters" || l.params.type === "rows") && l.params.of == null;
          if (own && l.params.letter_map === undefined) l.params.letter_map = map;
        }
      const doc: SeqViewSpecDoc = { format: SvsFormat, version: SvsVersion, root: toDoc(root) };
      if (revision !== undefined) doc.revision = revision;
      if (Object.keys({ ...meta, ...more }).length) doc.metadata = { ...meta, ...more };
      if (st.snapshots.length) doc.snapshots = resolve(st.snapshots) as Snapshot[];
      return doc;
    },
  };
  return b;
}
export type SeqViewSpecBuilder = ReturnType<typeof createSeqViewSpec>;

/** Base pairs from a WUSS (or dot-bracket) string: a `pairs` annotation, category `base_pair`, over `target`. */
export function basePairs(
  data: DataBuilder,
  params: BuilderCommon & {
    target: NonNullable<AnnotationParams["target"]>;
    structure: string;
    notation?: "wuss";
    label?: string;
    provenance?: BuilderObject;
  },
): NodeHandle<"annotation"> {
  const { structure, notation: _, ...rest } = params;
  const partner = wussPairs(structure);
  const start: number[] = [];
  const bStart: number[] = [];
  partner.forEach((j, i) => {
    if (j <= i) return;
    start.push(i);
    bStart.push(j);
  });
  const value = { start, b_start: bStart, category: start.map(() => "base_pair") };
  return data.annotation({ ...rest, geometry: "pairs", value } as AnnotationParams);
}

/** The source a translated protein's provenance names (not exported: the import rule's plugin reads exports). */
const BUILDER = "@molstar/seq-core/builder";

/** `translate`'s params: `from` a DNA or RNA sequence handle the builder holds by `value`. */
export type TranslateParams = BuilderCommon & {
  from: NodeHandle<"sequence">;
  label?: string;
  cds?: [number, number][] | null;
  start?: number;
  strand?: "+" | "-";
  genetic_code?: number;
};

/**
 * A protein `sequence` (`ref`) whose `value` is `from`'s coding range translated now, provenance method `translation`,
 * and its `codon` mapping `<ref>.codon` to `from` (seqspec:format#mappings). Throws when `from` is not a sequence held
 * by value: the builder reads no sources.
 */
export function translate(data: DataBuilder, params: TranslateParams): MoleculeHandle<"sequence"> {
  const { from, label, cds, start, strand, genetic_code, ref, custom } = params;
  const value = from.suffix === null && from.node.kind === "sequence" ? from.node.params.value : undefined;
  if (typeof value !== "string") {
    throw new TypeError(
      "translate: `from` must be a sequence handle given by value; a sequence read from a source is translated by whoever reads it.",
    );
  }
  const why = notNucleotide(value, from.node.params.alphabet);
  if (why) throw new TypeError(`translate: \`from\` must be DNA or RNA; ${why}.`);
  const [s, gc] = [strand ?? "+", genetic_code ?? 1];
  const residues = translateDna(value, cds ?? null, start ?? 0, s, gc);
  const parameters = { from, ...(cds ? { cds } : { start: start ?? 0 }), strand: s, genetic_code: gc };
  const provenance = { source: BUILDER, method: "translation", parameters };
  const protein = data.sequence({ ref, custom, value: residues, alphabet: "protein", label, provenance });
  const codon = { from: protein, to: from, method: "codon", cds, start, strand, genetic_code } as const;
  const m = data.mapping(ref === undefined ? codon : { ...codon, ref: `${ref}.codon` });
  if (ref === undefined) m.node.derived = { node: protein.node, suffix: ".codon" };
  return protein;
}
