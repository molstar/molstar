/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's node providers (seqspec:model#alignments-and-trees, seqspec:format#mappings; seqspec:format#extensions).
// NP1a: `nj` and `upgma` over the distances of their parent alignment's rows (`@molstar/seq-core/algorithms`). NP1b: `realign`,
// the inputs' residues through a registered MSA aligner, and `codon_rows`, a `codon` mapping per row pair of a protein
// and a cDNA alignment. A node provider runs at load, takes the stage's abort signal (a build that a newer request
// replaced stops at its next slice), and its result is cached by the content of what it read and its normalized options
// (`cached`), so a re-load with the same inputs computes nothing.

import {
  AlignerOptionsError,
  alignmentDistances,
  type GuideTree,
  midpointRoot,
  neighborJoining,
  upgma,
} from "@molstar/seq-core/algorithms";
import type { AlignmentInput, MappingInput, TreeInput } from "@molstar/seq-core/model";
import type { AlignerInput } from "../contract.ts";
import { canonical } from "../patch/tree-patch.ts";
import type { NodeProviderEntry } from "../registry.ts";
import {
  type FieldType,
  int,
  list,
  nullable,
  report as reportField,
  str,
  tuple,
  typeError,
} from "../schema/field-schema.ts";
import { OptionalField, object, SimpleParamsSchema } from "../schema/params-schema.ts";
import { report } from "./issues.ts";

/** Results kept: a small LRU of the latest (a tree of at most `2 * MAX_TREE_ROWS` nodes, a realigned family). */
const KEPT = 16;
const results = new Map<string, unknown>();

/** Two independent 32-bit hashes (FNV-1a and a multiplicative one) of text, folded into `h`. */
function hashInto(h: Uint32Array, text: string) {
  let a = h[0];
  let b = h[1];
  for (let i = 0; i < text.length; i++) {
    const c = text.charCodeAt(i);
    a = Math.imul(a ^ c, 16777619);
    b = Math.imul(b + c + 1, 0x9e3779b1) ^ (b >>> 15);
  }
  h[0] = a;
  h[1] = b;
}

/**
 * A key for what a provider reads: its name and normalized options (canonical: keys sorted at every depth), and the
 * content of `texts` and `more`, each string's length hashed before it, so no two lists of strings share a key by
 * where their boundaries fall (GV-01: labels `["a\nb", "c"]` and `["a", "b\nc"]` gave one tree's labels to the other).
 */
export function contentKey(name: string, options: unknown, texts: readonly string[], more: readonly string[]): string {
  const h = Uint32Array.of(2166136261, 1);
  for (const t of texts) {
    hashInto(h, `${t.length}:`);
    hashInto(h, t);
  }
  hashInto(h, "\0");
  for (const t of more) {
    hashInto(h, `${t.length}:`);
    hashInto(h, t);
  }
  return `${name}|${canonical(options ?? {})}|${texts.length}|${more.length}|${h[0].toString(36)}${h[1].toString(36)}`;
}

/** Each aligner entry's own number: results are keyed by the entry, not its name, which another registry may reuse. */
const entryIds = new WeakMap<object, number>();
let entryCount = 0;
function entryId(e: object): number {
  let id = entryIds.get(e);
  if (id === undefined) {
    id = ++entryCount;
    entryIds.set(e, id);
  }
  return id;
}

/** A key for what a tree provider reads: the provider, its normalized options, and the alignment's rows and labels. */
export function inputKey(name: string, options: unknown, a: AlignmentInput): string {
  return contentKey(name, options, a.gapped, [a.flavor, ...(a.labels ?? [])]);
}

/** The settled result for `key`, or `make()`'s, kept. An aborted or failed run is not kept; null is not kept. */
export async function cached<T>(key: string, make: () => Promise<T | null>): Promise<T | null> {
  const hit = results.get(key) as T | undefined;
  if (hit) {
    results.delete(key);
    results.set(key, hit);
    return hit;
  }
  const r = await make();
  if (r) {
    results.set(key, r);
    if (results.size > KEPT) results.delete(results.keys().next().value as string);
  }
  return r;
}

/** Forgets the kept results (tests). */
export function clearNodeProviderCache() {
  results.clear();
}

/** The alignment's columns as codes, row-major, 0 where a row has no letter: upper case, `-` and `.` are gaps, an A3M's lower case is out. */
export function alignmentCodes(a: AlignmentInput): { codes: Uint8Array; rows: number; columns: number } {
  const a3m = a.flavor === "a3m";
  const rows = a.gapped.length;
  const strip = (s: string) => (a3m ? s.replace(/[a-z.]/g, "") : s);
  const text = a.gapped.map(strip);
  const columns = text.reduce((m, s) => Math.max(m, s.length), 0);
  const codes = new Uint8Array(rows * columns);
  for (let r = 0; r < rows; r++) {
    const s = text[r];
    const o = r * columns;
    for (let c = 0; c < s.length; c++) {
      const k = s.charCodeAt(c);
      codes[o + c] = k === 45 || k === 46 || k > 255 ? 0 : k >= 97 && k <= 122 ? k - 32 : k;
    }
  }
  return { codes, rows, columns };
}

/** A built tree as the model's input: leaf k is row k, named by its label (leaves matched by order, not by label). */
export function treeInput(t: GuideTree, labels: readonly string[] | undefined): TreeInput {
  const n = t.parent.length;
  const leafRow = new Int32Array(n).fill(-1);
  const label: (string | null)[] = new Array(n).fill(null);
  for (let r = 0; r < t.leaves; r++) {
    leafRow[r] = r;
    label[r] = labels?.[r] ?? null;
  }
  return { parent: t.parent, branchLength: t.branchLength, label, leafRow };
}

/** The one distance: `identity` was the same number, and says so. */
const distance: FieldType<"p"> = {
  name: '"p"',
  check(v, path, ctx) {
    if (v === "p") return;
    if (v === "identity")
      reportField(
        ctx,
        "svs.params.type",
        path,
        '"identity" is the same distance as "p": write "p" (mismatches over the columns where both rows have a letter)',
      );
    else typeError(ctx, path, '"p"', v);
  },
};
const treeOptions = SimpleParamsSchema({
  distance: OptionalField(
    distance,
    "p",
    "p: mismatches over the columns where both rows have a letter, as a share of them",
  ),
});

function treeProvider(
  name: string,
  build: (d: Float64Array, n: number, signal: AbortSignal) => Promise<GuideTree>,
): NodeProviderEntry {
  return {
    name,
    kind: "tree",
    options: treeOptions,
    async compute(_ctx, inputs, opts, signal) {
      const a = inputs[0] as AlignmentInput;
      return cached(inputKey(name, opts, a), async () => {
        const { codes, rows, columns } = alignmentCodes(a);
        const tree = await build(await alignmentDistances(codes, rows, columns, signal), rows, signal);
        return treeInput(tree, a.labels);
      });
    },
  };
}

/** One input of `realign`: its residues (what an aligner reads) and the label its row keeps (null: the input's key). */
export type RealignInput = AlignerInput & { readonly label: string | null };

/** An alignment as `codon_rows` reads it: its key and its rows (held ones only: a streamed alignment is not paired). */
export interface AlignmentSummary {
  readonly key: string;
  readonly rowCount: number;
}

const decoder = new TextDecoder();
const realignOptions = SimpleParamsSchema({
  aligner: OptionalField(str, "biomsa", "a registered aligner that aligns many sequences (biomsa)"),
  aligner_options: OptionalField(object(), {}, "the aligner's params, checked by its table; {}: its defaults"),
});

const realign: NodeProviderEntry = {
  name: "realign",
  kind: "alignment",
  options: realignOptions,
  async compute(ctx, inputs, opts, signal) {
    const o = opts as { aligner: string; aligner_options: Record<string, unknown> };
    const ins = inputs as readonly RealignInput[];
    const entry = ctx.registries.get("aligners", o.aligner);
    if (!entry) {
      report(ctx, "svs.provider.unavailable", "", `aligner "${o.aligner}" is not registered`);
      return null;
    }
    if (!entry.alignMany) {
      report(
        ctx,
        "svs.provider.unavailable",
        "",
        `aligner "${o.aligner}" aligns two sequences only: realign needs one that aligns many (biomsa)`,
      );
      return null;
    }
    const texts = ins.map((m) => decoder.decode(m.code));
    const labels = ins.map((m) => m.label ?? String(m.key));
    const key = contentKey(`realign#${entryId(entry)}`, o, texts, [...labels, ...ins.map((m) => m.alphabet ?? "")]);
    return cached<AlignmentInput>(key, async () => {
      let rows: string[] | null;
      try {
        rows = (await entry.alignMany?.(ctx, ins, o.aligner_options, signal)) ?? null;
      } catch (e) {
        signal.throwIfAborted();
        if (!(e instanceof AlignerOptionsError)) throw e;
        report(ctx, "svs.params.type", "/params/options/aligner_options", e.message);
        return null;
      }
      if (rows == null) return null;
      const alphabets = new Set(ins.map((m) => m.alphabet));
      const [alphabet] = alphabets;
      return {
        key: "",
        gapped: rows,
        flavor: "aligned",
        labels,
        ...(alphabets.size === 1 && alphabet ? { alphabet } : {}),
      };
    });
  },
};

const rowPairs = nullable(list(tuple(int, int)));
const codonRows: NodeProviderEntry = {
  name: "codon_rows",
  kind: "mapping",
  options: SimpleParamsSchema({
    rows: OptionalField(rowPairs, null, "[protein row, cDNA row] pairs; null: row k with row k"),
    genetic_code: OptionalField(int, 1, "an NCBI genetic code id (gc.prt; 1 the standard code)"),
  }),
  async compute(ctx, inputs, opts) {
    const [protein, cdna] = inputs as readonly [AlignmentSummary, AlignmentSummary];
    const o = opts as { rows: [number, number][] | null; genetic_code: number };
    let pairs = o.rows;
    if (pairs === null) {
      if (protein.rowCount !== cdna.rowCount) {
        report(
          ctx,
          "svs.mapping.codon",
          "/params/options/rows",
          `"${protein.key}" has ${protein.rowCount} rows and "${cdna.key}" ${cdna.rowCount}: give rows, the [protein row, cDNA row] pairs`,
        );
        return null;
      }
      pairs = Array.from({ length: protein.rowCount }, (_, k) => [k, k]);
    }
    const bad = pairs.find(([p, c]) => p < 0 || p >= protein.rowCount || c < 0 || c >= cdna.rowCount);
    if (bad) {
      report(
        ctx,
        "svs.mapping.codon",
        "/params/options/rows",
        `the pair [${bad}] names a row that does not exist: "${protein.key}" has rows 0 to ${protein.rowCount - 1}, "${cdna.key}" 0 to ${cdna.rowCount - 1}`,
      );
      return null;
    }
    // Keys are suffixes: the compiler prefixes the node's ref.
    return pairs.map(
      ([p, c]): MappingInput => ({
        key: `row-${p}`,
        from: `${protein.key}/${p}`,
        to: `${cdna.key}/${c}`,
        method: "codon",
        geneticCode: o.genetic_code,
      }),
    );
  },
};

export const KernelNodeProviders: readonly NodeProviderEntry[] = [
  treeProvider("nj", async (d, n, signal) => midpointRoot(await neighborJoining(d, n, signal))),
  treeProvider("upgma", upgma),
  realign,
  codonRows,
];
