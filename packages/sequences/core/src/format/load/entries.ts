/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's aligners (`@molstar/seq-core/algorithms`, their plain-data params as a params table) and the four kernel
// providers over the model's (seqspec:model#providers), as registry entries (seqspec:format#extensions).

import { type Aligner, AlignerOptionsError, Aligners } from "@molstar/seq-core/algorithms";
import {
  COVERAGE_PALETTE,
  type Geometry,
  MoleculeKind,
  type ProviderResult,
  Providers,
  type SeqSet,
} from "@molstar/seq-core/model";
import type { SeqQLQuery } from "@molstar/seq-core/query";
import type { AlignerInput, StageContext } from "../contract.ts";
import type { ThemeRef } from "../doc.ts";
import { computeQuery } from "../query.ts";
import type { AlignerEntry, ProviderEntry } from "../registry.ts";
import { type FieldType, float, int, literal, nullable, refined, report, str } from "../schema/field-schema.ts";
import { OptionalField, SimpleParamsSchema, type UnionParamsSchema } from "../schema/params-schema.ts";
import { SvsTreeSchema } from "../tables/svs-tree.ts";
import { report as reportAt } from "./issues.ts";
import { RowProviders } from "./row-providers.ts";

const atLeast = (t: FieldType<number>, min: number): FieldType<number> => ({
  name: `${t.name} (>= ${min})`,
  check(v, path, ctx) {
    t.check(v, path, ctx);
    if (typeof v === "number" && v < min) report(ctx, "svs.params.type", path, `${v} is below the minimum, ${min}`);
  },
});

const TYPE: Readonly<Record<string, string>> = { protein: "amino", dna: "nucleic", rna: "nucleic" };
const decoder = new TextDecoder();

/** An aligner as an entry: its params a table; `type` from the molecules' alphabets; ligands refused. */
export function alignerEntry(al: Aligner): AlignerEntry {
  const fields = Object.fromEntries(
    Object.entries(al.params).map(([k, p]) => {
      let t: FieldType<unknown> =
        p.type === "choice"
          ? literal(...(p.choices ?? []))
          : p.type === "str"
            ? p.pattern
              ? refined("str", (s) => new RegExp(`^(?:${p.pattern})$`).test(s), "svs.params.type", `/${p.pattern}/`)
              : str
            : p.type === "int"
              ? int
              : float;
      if (p.min !== undefined) t = atLeast(t as FieldType<number>, p.min) as FieldType<unknown>;
      return [k, OptionalField(p.nullable ? nullable(t) : t, p.default, p.description)];
    }),
  );
  return {
    name: al.name,
    options: SimpleParamsSchema(fields),
    async align(ctx: StageContext, a: AlignerInput, b: AlignerInput, options: unknown) {
      const ligand = [a, b].find((m) => m.kind !== MoleculeKind.polymer);
      if (ligand) {
        reportAt(
          ctx,
          "svs.ligand.unsupported",
          "",
          `"${ctx.set?.strings.get(ligand.key)}" is a ligand: it cannot be aligned`,
        );
        return null;
      }
      const o = { ...(options as Record<string, string | number | null>) };
      const type = a.alphabet && TYPE[a.alphabet];
      if ("type" in al.params && (o.type == null || o.type === "auto") && type && type === TYPE[b.alphabet ?? ""]) {
        o.type = type;
      }
      try {
        const cigar = await al.align(decoder.decode(a.code), decoder.decode(b.code), o);
        if (cigar === null) reportAt(ctx, "svs.mapping.align", "", `${al.name} cannot align these sequences`);
        return cigar;
      } catch (e) {
        ctx.signal.throwIfAborted();
        // Only invalid options are the options' fault; a size limit or any other failure is the aligner's.
        const message = e instanceof Error ? e.message : String(e);
        if (e instanceof AlignerOptionsError) reportAt(ctx, "svs.params.type", "/params/options", message);
        else reportAt(ctx, "svs.mapping.align", "", message);
        return null;
      }
    },
    ...(al.alignMany && {
      async alignMany(ctx: StageContext, inputs: readonly AlignerInput[], options: unknown, signal: AbortSignal) {
        const ligand = inputs.find((m) => m.kind !== MoleculeKind.polymer);
        if (ligand) {
          reportAt(
            ctx,
            "svs.ligand.unsupported",
            "",
            `"${ctx.set?.strings.get(ligand.key)}" is a ligand: it cannot be aligned`,
          );
          return null;
        }
        const o = { ...(options as Record<string, string | number | null>) };
        const types = new Set(inputs.map((m) => m.alphabet && TYPE[m.alphabet]));
        const [type] = types;
        if ("type" in al.params && (o.type == null || o.type === "auto") && type && types.size === 1) o.type = type;
        const rows = await al.alignMany?.(
          inputs.map((m) => decoder.decode(m.code)),
          o,
          signal,
        );
        if (rows == null) reportAt(ctx, "svs.mapping.align", "", `${al.name} cannot align these sequences`);
        return rows ?? null;
      },
    }),
  };
}

type Opts = Record<string, unknown>;
const cases = (SvsTreeSchema.nodes.computed.params as UnionParamsSchema).cases;

/** A kernel provider: the model's, its options the `computed` node's normalized params. */
function provider(
  name: string,
  geometry: Geometry,
  run: (set: SeqSet, o: Opts) => ProviderResult,
  defaultColor?: ThemeRef,
): ProviderEntry {
  return {
    name,
    geometry,
    defaultColor,
    options: cases[name],
    async compute(ctx, options) {
      if (!ctx.set) return { key: "", geometry, columns: {} };
      const r = run(ctx.set, options as Opts);
      for (const i of r.issues) ctx.report(i);
      // A provider that could not run (a lazy alignment, a generic alphabet) gives none: the node is unavailable.
      return r.annotation;
    },
  };
}

const strs = (v: unknown) => (v == null ? null : (v as string[]));
const stats = (o: Opts) => ({
  key: "",
  alignment: o.alignment as string,
  rowSubset: o.row_subset as number[] | null,
  gaps: o.gaps as "ignore" | "count",
});

/** The query provider (seqspec:seqql#in-documents): its own runner, not a model `Providers` function. */
const QueryProvider: ProviderEntry = {
  name: "query",
  geometry: "ranges",
  options: cases.query,
  compute: (ctx, options) => computeQuery(ctx, options as { query: SeqQLQuery; label: string | null }),
};

export const KernelProviders: readonly ProviderEntry[] = [
  provider(
    "coverage",
    "ranges",
    (set, o) =>
      Providers.coverage(set, {
        key: "",
        of: o.of == null ? null : [o.of as string],
        reference: o.reference as string | null,
      }),
    { theme: "field_palette", params: { palette: COVERAGE_PALETTE } },
  ),
  provider("column_stats", "values", (set, o) =>
    Providers.columnStats(set, {
      ...stats(o),
      measure: o.measure as "conservation",
      referenceRow: o.reference_row as number | null,
    }),
  ),
  provider("profile", "matrix", (set, o) => Providers.profile(set, stats(o))),
  ...RowProviders,
  QueryProvider,
  provider("paired", "ranges", (set, o) =>
    Providers.paired(set, {
      key: "",
      pairs: o.pairs as string[],
      partner: strs(o.partner),
      categories: strs(o.categories),
    }),
  ),
];

export const KernelAligners: readonly AlignerEntry[] = Aligners.map(alignerEntry);
