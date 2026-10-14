/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A genome-scale table for the validation cost tests, heap guards and bench (VW2): one `values` annotation over 24
// chromosomes whose `target` is dictionary-coded (seqspec:format#positions), as the synthetic genome's GC table with
// `gcTable: "one"`, built as plain data with no builder.

import type { SeqViewSpecDoc } from "../../src/format/index.js";

// biome-ignore lint/suspicious/noExplicitAny: documents are edited freely as JSON.
type J = any;

export const CHROMOSOMES = Array.from({ length: 24 }, (_, i) => `chr${i + 1}`);

export interface CodedOptions {
  /** Trivial snapshots, each hiding or showing the track. */
  snapshots?: number;
  /** The codes as given (a Proxy counting reads); default a `Uint8Array`, chromosome after chromosome. */
  codes?: ArrayLike<number | null>;
  /** The dictionary; default the 24 chromosomes. */
  dict?: string[];
}

/** A values table of `n` rows, its target coded; `gc` is the annotation (data children index 24). */
export function codedDoc(n: number, o: CodedOptions = {}): SeqViewSpecDoc {
  const dict = o.dict ?? CHROMOSOMES;
  let codes = o.codes;
  if (!codes) {
    const c = new Uint8Array(n);
    for (let i = 0; i < n; i++) c[i] = Math.floor((i * dict.length) / n);
    codes = c;
  }
  const value = new Float64Array(n);
  for (let i = 0; i < n; i++) value[i] = (i % 100) / 100;
  const gc = {
    kind: "annotation",
    ref: "gc",
    params: { geometry: "values", value: { target: { dict, codes }, value } },
  };
  const doc: J = {
    format: "seq_view_spec",
    version: "1.0.0",
    root: {
      kind: "root",
      children: [
        {
          kind: "data",
          children: [
            ...CHROMOSOMES.map((ref) => ({ kind: "sequence", ref, params: { alphabet: "dna", value: "ACGT" } })),
            gc,
          ],
        },
        {
          kind: "view",
          ref: "v",
          params: { label: "Genome" },
          children: [
            {
              kind: "section",
              ref: "s",
              children: [
                ...CHROMOSOMES.map((c) => ({ kind: "block", ref: `b.${c}`, params: { sequence: c } })),
                {
                  kind: "track",
                  ref: "t",
                  params: { label: "GC" },
                  children: [{ kind: "layer", ref: "l", params: { type: "values", data: "gc" } }],
                },
              ],
            },
          ],
        },
      ],
    },
    snapshots: Array.from({ length: o.snapshots ?? 0 }, (_, k) => ({
      ref: `s${k}`,
      ops: [{ op: "set", ref: "t", params: { hidden: k % 2 === 0 } }],
    })),
  };
  return doc;
}

/** The annotation node of `codedDoc`. */
export const gcNode = (d: SeqViewSpecDoc): J => (d.root.children as J)[0].children[24];

/** The same document with the annotation's `params` a new object (its table shared): its own checks run again. */
export function freshParams(d: SeqViewSpecDoc): SeqViewSpecDoc {
  const data: J = (d.root.children as J)[0];
  const gc = gcNode(d);
  const children = [...data.children.slice(0, 24), { ...gc, params: { ...gc.params } }];
  return { ...d, root: { ...d.root, children: [{ ...data, children }, ...(d.root.children as J).slice(1)] } };
}

/** An array of codes behind a Proxy that counts the reads of its cells. */
export function countingCodes(codes: number[]): { codes: number[]; reads: () => number } {
  let reads = 0;
  const proxy = new Proxy(codes, {
    get(target, key, receiver) {
      if (typeof key === "string" && /^[0-9]+$/.test(key)) reads++;
      return Reflect.get(target, key, receiver);
    },
  });
  return { codes: proxy, reads: () => reads };
}
