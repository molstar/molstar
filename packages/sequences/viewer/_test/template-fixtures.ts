/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b's fixtures: a protein in memory and streamed sources whose manifests carry dynamic nodes (annotations on the
// protein, and nested sources with manifests of their own), a section of the protein's letters with a track template
// (`for_each`) of ranges, and a provider per stream id that records the manifests asked for.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type {
  DataPatch,
  Loaders,
  Manifest,
  SeqViewSpecDoc,
  StreamNode,
  StreamProvider,
} from "@molstar/seq-core/format";
import type { Streams } from "../src/streams.ts";
import { opened } from "./layout-fixtures.ts";

export const P = "MKVLAAGGHHKKLLMMNNPPQQRRSSTTVVWWYY".repeat(6);
export const SECTION = "v.main";
export interface Template {
  kind: "annotation";
  source: string;
  depth?: number | null;
  reaching?: boolean;
}

/** A range on the protein (`at` its start), keyed `key`. */
export const range = (key: string, at = 1, label?: string): StreamNode => ({
  kind: "annotation",
  key,
  params: { geometry: "ranges", target: "p", label, value: { start: [at], stop: [at + 4] } },
});
/** A nested source of `stream`'s items. */
export const nested = (key: string, stream: string): StreamNode => ({
  kind: "source",
  key,
  params: { stream, data_kind: "items" },
});
const manifest = (nodes: StreamNode[]): Manifest => ({
  systems: [{ key: "p", extent: [0, P.length], numbering: null, label: null }],
  levels: [],
  nodes,
});

/** The protein, the `genes` and `notes` sources, and a section of its letters, the template and a ruler. */
export function templateDoc(o: { for_each?: Template; max_depth?: number } = {}) {
  const b = createSeqViewSpec({ title: "templates" });
  const p = b.data.sequence({ ref: "p", value: P });
  b.data.source({ ref: "genes", stream: "genes", data_kind: "items", max_depth: o.max_depth });
  b.data.source({ ref: "notes", stream: "notes", data_kind: "items" });
  const v = b.view({ ref: "v", label: "v" });
  v.block({ sequence: p });
  v.track({ ref: "seq", label: "Seq" }).letters();
  v.track({ ref: "tpl", for_each: o.for_each ?? { kind: "annotation", source: "genes" } }).ranges({ data: "$node" });
  v.track({ ref: "ruler", label: "Ruler" }).ruler();
  return b.getDocument() as SeqViewSpecDoc;
}

/** A provider per stream id over its manifest's nodes; `asked`: the manifests asked for, by stream id. */
export function templateLoaders(manifests: Record<string, StreamNode[]>) {
  const asked: string[] = [];
  const provider = (id: string): StreamProvider => ({
    fetch(r) {
      if (r.system !== null) return Promise.reject(new Error("no data regions"));
      asked.push(id);
      const at = { source: r.source, system: null, start: 0, stop: 0, row_start: 0, row_stop: 0, level: -1 };
      return Promise.resolve({ ...at, kind: "manifest", manifest: manifest(manifests[id] ?? []) } as DataPatch);
    },
  });
  const loaders: Loaders = {
    load: (r) => Promise.reject(new Error(`no file ${r.uri}`)),
    stream: (s) => (s.stream && s.stream in manifests ? provider(s.stream) : null),
  };
  return { loaders, asked };
}

/** The document opened over `manifests` (`genes` and `notes` default to none). */
export async function openTemplates(
  manifests: Record<string, StreamNode[]>,
  o: Parameters<typeof templateDoc>[0] = {},
) {
  const l = templateLoaders({ genes: [], notes: [], ...manifests });
  const t = await opened(templateDoc(o), { loaders: l.loaders, size: { [SECTION]: [1000, 600] } });
  const manager = t.ctx.streams as Streams;
  /** The shown tracks' refs, in order. */
  const tracks = () => {
    const v = t.ctx.graph.assembled.value?.view;
    const rv = v?.resolved;
    const out: string[] = [];
    for (let k = 0; v && rv && k < rv.tracks.count; k++)
      if (rv.tracks.hidden[k] === 0) out.push(v.set.strings.get(rv.tracks.ref[k]));
    return out;
  };
  /** The instances' tracks (`tpl/<node>`), as their nodes. */
  const instances = () => tracks().flatMap((r) => (r.startsWith("tpl/") ? [r.slice(4)] : []));
  /** The prepared state of a layer, by its ref. */
  const state = (ref: string) => {
    const v = t.ctx.graph.assembled.value?.view;
    for (let k = 0; v && k < v.resolved.layers.count; k++)
      if (v.set.strings.get(v.resolved.layers.ref[k]) === ref) return v.prepared[k].state;
    return undefined;
  };
  const frame = async () => {
    t.scheduler.advance(0);
    await t.ctx.graph.update();
  };
  return { ...t, ...l, manager, tracks, instances, state, frame };
}
