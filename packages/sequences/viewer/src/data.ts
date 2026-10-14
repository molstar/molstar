/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Data without a document (seqspec:viewer#data-mode): build inputs as one graph node per entry, or a whole set as one
// node,
// and views written with the format's builder against the set's keys, by default one `tracks` view. Validation in
// data mode resolves refs against the incoming data: the set, or a resolver over the build inputs' keys (`keysOf`).

import { createSeqViewSpec, type SeqViewSpecBuilder } from "@molstar/seq-core/builder";
import { inputSources, isRefSyntax, type SeqViewSpecDoc, streamedSources, withDynamic } from "@molstar/seq-core/format";
import {
  type AnnotationInput,
  ColumnRole,
  ColumnRoles,
  type Geometry,
  KeyKind,
  type KeyResolver,
  type KeyTarget,
  MoleculeKind,
  type SeqSet,
  type SetInput,
  Systems,
} from "@molstar/seq-core/model";
import type { DataEntry, DataMode } from "./graph.ts";

export interface SeqDataInput {
  readonly set?: SeqSet;
  readonly inputs?: SetInput;
}
export interface SetDataOptions {
  readonly views?: (b: SeqViewSpecBuilder) => void;
}

const LISTS = ["molecules", "units", "complexes", "alignments", "mappings", "annotations"] as const;

export function dataMode(data: SeqDataInput): DataMode {
  if (data.set) return { set: data.set };
  const entries = new Map<string, DataEntry>();
  for (const list of LISTS) for (const entry of data.inputs?.[list] ?? []) entries.set(entry.key, { list, entry });
  return { entries };
}

/**
 * What refs name in data mode: the set itself, or the build inputs' keys (row refs `R/<k>` included); and, under the
 * streamed sources the data reads, dynamic nodes not held yet (ST1b-a; seqspec:sources#dynamic-nodes).
 */
export function dataKeys(data: SeqDataInput): KeyResolver {
  if (data.set) return withDynamic(data.set, streamedSources(data.set));
  return withDynamic(inputKeys(data), inputSources(data.inputs?.annotations));
}

function inputKeys(data: SeqDataInput): KeyResolver {
  const i = data.inputs;
  const keys = new Map<string, number>();
  const rows: number[] = [];
  const molecules = (i?.molecules ?? []).map((m, id) => {
    keys.set(m.key, id * 8 + KeyKind.molecule);
    return { kind: m.kind === "ligand" ? MoleculeKind.ligand : MoleculeKind.polymer };
  });
  let unit = 0;
  for (const u of i?.units ?? []) keys.set(u.key, unit++ * 8 + KeyKind.unit);
  for (const m of i?.molecules ?? []) if (m.ownUnit) keys.set(`${m.key}/unit`, unit++ * 8 + KeyKind.unit);
  for (const a of i?.alignments ?? []) rows.push(a.lazy?.rowCount ?? a.gapped.length);
  const lists = [
    ["alignments", KeyKind.alignment],
    ["complexes", KeyKind.complex],
    ["mappings", KeyKind.mapping],
    ["annotations", KeyKind.annotation],
  ] as const;
  for (const [list, kind] of lists)
    for (let id = 0, l = i?.[list] ?? []; id < l.length; id++) keys.set(l[id].key, id * 8 + kind);
  const byName = (key: string, out: KeyTarget) => {
    const k = keys.get(key);
    if (k !== undefined) {
      out.kind = k % 8;
      out.id = Math.floor(k / 8);
      out.row = -1;
      return true;
    }
    // A row ref `<alignment>/<k>`, its index read as the set reads it (digits only).
    const row = Systems.rowOf(key);
    const a = row >= 0 ? keys.get(key.slice(0, key.lastIndexOf("/"))) : undefined;
    if (a === undefined || a % 8 !== KeyKind.alignment || row >= rows[Math.floor(a / 8)]) return false;
    out.kind = KeyKind.row;
    out.id = Math.floor(a / 8);
    out.row = row;
    return true;
  };
  return { byName, molecules };
}

/** A column a `rows` annotation keeps as a field (as the model's build does: not placing, `state` and `node` too). */
const isField = (c: string) =>
  c === "state" ||
  c === "node" ||
  Object.values(ColumnRoles).every((r) => r[c] === undefined || r[c] === ColumnRole.descriptive);

interface Keys {
  units: { key: string; ligand: boolean }[];
  alignments: string[];
  annotations: { key: string; geometry: Geometry }[];
  /** Per alignment key: its `rows` annotations' fields as data columns. */
  rows: Map<string, { data: string; field: string }[]>;
  all: string[];
}

/** Units (polymers, then ligands), alignments, annotations and rows' fields by key. */
function keysOf(data: SeqDataInput): Keys {
  const k: Keys = { units: [], alignments: [], annotations: [], rows: new Map(), all: [] };
  const rowsOf = (alignment: string, data: string, fields: string[]) => {
    if (!isRefSyntax(data, true)) return;
    const list = k.rows.get(alignment) ?? [];
    k.rows.set(alignment, list);
    for (const field of fields) list.push({ data, field });
  };
  const s = data.set;
  if (s) {
    for (let u = 0; u < s.units.count; u++) {
      const ligand = s.molecules[s.units.molecule[u]]?.kind === MoleculeKind.ligand;
      if (s.units.alive[u]) k.units.push({ key: s.strings.get(s.units.key[u]), ligand });
    }
    for (const a of s.alignments) if (a) k.alignments.push(s.strings.get(a.key));
    for (const a of s.annotations) {
      if (!a) continue;
      const key = s.strings.get(a.key);
      const target = s.alignments[a.slots.alignment[0]];
      if (a.geometry !== "rows") k.annotations.push({ key, geometry: a.geometry });
      else if (target)
        rowsOf(
          s.strings.get(target.key),
          key,
          a.fields.map((f) => s.strings.get(f.name)),
        );
    }
    for (const list of [s.molecules, s.complexes, s.alignments, s.mappings, s.annotations])
      for (const x of list) if (x) k.all.push(s.strings.get(x.key));
  } else {
    const i = data.inputs;
    const molecules = i?.molecules ?? [];
    const ligands = new Set(molecules.filter((m) => m.kind === "ligand").map((m) => m.key));
    for (const u of i?.units ?? []) k.units.push({ key: u.key, ligand: ligands.has(u.molecule) });
    for (const m of molecules) if (m.ownUnit) k.units.push({ key: `${m.key}/unit`, ligand: ligands.has(m.key) });
    for (const a of i?.alignments ?? []) k.alignments.push(a.key);
    for (const a of i?.annotations ?? ([] as AnnotationInput[])) {
      const target = a.target ?? "";
      if (a.geometry !== "rows") k.annotations.push({ key: a.key, geometry: a.geometry });
      else rowsOf(target, a.key, Object.keys(a.columns ?? {}).filter(isField));
    }
    for (const list of LISTS) for (const x of i?.[list] ?? []) k.all.push(x.key);
  }
  for (const u of k.units) k.all.push(u.key);
  k.units = [...k.units.filter((u) => !u.ligand), ...k.units.filter((u) => u.ligand)];
  return k;
}

/**
 * `view.<key>` in the ref syntax (seqspec:format#positions), with a numeric suffix past any ref in `used`; added to it.
 */
function slugRef(key: string, used: Set<string>): string {
  const base = `view.${key.replace(/[^A-Za-z0-9_.-]/g, "-")}`;
  let ref = base;
  for (let k = 2; used.has(ref); k++) ref = `${base}-${k}`;
  used.add(ref);
  return ref;
}

/**
 * One `tracks` view, ref `view`, section `view.main`: the blocks, a ruler, letters, a track per annotation, and one
 * `rows` track for the block's alignment when it has `rows` annotations (a `rows` annotation is no track of its own).
 * An annotation whose key is not in the ref syntax gets none: no layer could name it.
 */
function defaultView(data: SeqDataInput) {
  return (b: SeqViewSpecBuilder) => {
    const { units, alignments, annotations, rows, all } = keysOf(data);
    const used = new Set([...all, "view", "view.main", "view.ruler", "view.letters"]);
    const v = b.view({ ref: "view", label: "View" });
    for (const u of units) v.block({ unit: u.key });
    const block = units.length === 0 ? alignments[0] : undefined;
    if (block !== undefined) v.block({ alignment: block });
    v.track({ ref: "view.ruler" }).ruler();
    v.track({ ref: "view.letters" }).letters();
    const columns = block === undefined ? undefined : rows.get(block);
    if (columns?.length) v.track({ ref: slugRef("rows", used) }).rows({ data_columns: columns } as never);
    for (const a of annotations)
      if (isRefSyntax(a.key, true)) {
        // A grid has no type of its name: a dot plot shows it.
        const type = a.geometry === "grid" ? "dot_plot" : a.geometry;
        v.track({ ref: slugRef(a.key, used) }).layer({ type, data: a.key });
      }
  };
}

/** The views document `setData` shows. */
export function dataViews(data: SeqDataInput, opts: SetDataOptions = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  (opts.views ?? defaultView(data))(b);
  return b.getDocument();
}
