/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL and Mol* (seqspec:libraries#the-mol-adapter): Mol*'s selection as the `structure` input a run
// reads with `(input structure)`, the extension `molstar.select` (alias `molql`) that runs a PyMOL, VMD, Jmol or
// MolScript selection on the linked structure and gives its residues on their Seq* units, and a result back as a
// Mol* loci. Both directions go through the link (link.ts), so the query and the model know nothing of Mol*.

import { CoordinateSystemKind, Extents } from "@molstar/seq-core/model";
import { QueryAbort, type QueryExtension, SeqSelection, type SymbolRuntime } from "@molstar/seq-core/model/query";
import { defineSymbol, queryIssue } from "@molstar/seq-core/query";
import {
  QueryContext as MolQueryContext,
  type QueryFn as MolQueryFn,
  StructureElement,
  StructureSelection,
} from "molstar/lib/mol-model/structure";
import { compile } from "molstar/lib/mol-script/runtime/query/compiler";
import { Script } from "molstar/lib/mol-script/script";
import { _transpiler } from "molstar/lib/mol-script/transpilers/all";
import type { StructureLink } from "./link.ts";

const LANGS = ["pymol", "vmd", "jmol", "mol-script"] as const;
type Lang = (typeof LANGS)[number];

export const MOLSTAR_SELECT = defineSymbol({
  name: "molstar.select",
  alias: "molql",
  kind: "generator",
  args: [
    { name: "S", type: "str", positional: true, doc: 'a selection, as "polymer and (resn ZN around 3)"' },
    { name: "lang", type: { enum: LANGS }, optional: true, default: "pymol", doc: "its language; default pymol" },
  ],
  result: "query",
  signature: '(molql "S" :lang L)',
  doc: "the residues Mol* selects in the structure shown in 3D, on their sequences, one group",
  examples: [{ text: '(onto (molql "polymer and (resn ZN around 3)") p53)', says: "the zinc site's residues on p53" }],
  hints: [{ code: "seqql.symbol.unavailable", text: "it runs where a structure is shown in 3D and linked" }],
});

/** The definition alone: a host without Mol* checks and documents `molql`, and a run says it is unavailable. */
export const MOLSTAR_QUERY_SYMBOLS: QueryExtension = { namespace: "molstar", symbols: [MOLSTAR_SELECT] };

/** A selection's text as a Mol* query; throws Mol*'s parse or compile error. */
function molQuery(text: string, lang: Lang): MolQueryFn<StructureSelection> {
  const expr = lang === "mol-script" ? Script.toExpression(Script(text, lang)) : _transpiler[lang](text);
  return compile<StructureSelection>(expr);
}

/**
 * `molstar.select` with its runtime, over the link `link()` gives at run time (the structure shown now; null: none).
 * The text is read and compiled once, with the query; a run evaluates it on the link's structure, synchronously.
 */
export function molstarQueryExtension(link: () => StructureLink | null): QueryExtension {
  const runtime: SymbolRuntime = {
    symbol: MOLSTAR_SELECT,
    compile: (c) => {
      const text = c.raw<string>("0") ?? "";
      const lang = c.raw<Lang>("lang") ?? "pymol";
      let query: MolQueryFn<StructureSelection>;
      try {
        query = molQuery(text, lang);
      } catch (e) {
        const why = (e instanceof Error ? e.message : String(e)).split("\n")[0];
        const message = `Mol* cannot read "${text}" as ${lang}: ${why}; write a ${lang} selection, as "polymer and (resn ZN around 3)"`;
        const path = `${c.path}/args/0`;
        throw new QueryAbort(
          queryIssue("molstar.select.syntax", message, { path, signature: MOLSTAR_SELECT.signature }),
        );
      }
      const runs = Extents.create(64);
      return (ctx) => {
        const l = link();
        if (!l || l.set !== ctx.set) {
          return ctx.fail(
            "seqql.symbol.unavailable",
            c.path,
            "molql needs a structure shown in 3D over these sequences",
          );
        }
        const loci = StructureSelection.toLociWithSourceUnits(query(new MolQueryContext(l.structure)));
        runs.count = 0;
        l.seqRunsOf(loci, runs);
        if (runs.count === 0) {
          const n = StructureElement.Loci.size(loci);
          const name = l.structure.model.entryId || l.key;
          ctx.because(c.path, n ? `Mol* selected ${n} atoms of ${name}, none on a sequence here` : `none in ${name}`);
        }
        return SeqSelection.fromExtents(ctx.set, runs);
      };
    },
  };
  return { ...MOLSTAR_QUERY_SYMBOLS, runtimes: [runtime] };
}

/** Mol*'s selection on the link's structure as the `structure` input: its residues on their Seq* units, one group. */
export function structureInput(link: StructureLink, loci: StructureElement.Loci | null): SeqSelection {
  const runs = Extents.create(16);
  if (loci) link.seqRunsOf(loci, runs);
  return SeqSelection.fromExtents(link.set, runs);
}

/** A result back in 3D: its runs as one Mol* loci on the link's structure, through the link's routes. */
export function selectionLoci(link: StructureLink, sel: SeqSelection): StructureElement.Loci {
  const runs = Extents.create(16);
  SeqSelection.toExtents(sel, runs);
  return link.toLoci(runs);
}

/**
 * Which of the link's chains a result reaches, by auth id, and whether it names them (`exact`: every run on a linked
 * chain's unit) or reaches them from a sequence they share: a run on P04637 is any of 1TUP's three p53 chains.
 */
export function structureReach(link: StructureLink, sel: SeqSelection): { chains: string[]; exact: boolean } {
  const runs = Extents.create(16);
  SeqSelection.toExtents(sel, runs);
  const units = new Set(link.chains.flatMap((c) => (c.unit ? [c.unit.id] : [])));
  let exact = runs.count > 0;
  for (let i = 0; i < runs.count; i++) {
    if (runs.kind[i] !== CoordinateSystemKind.unit || !units.has(runs.id[i])) exact = false;
  }
  const members = link.toMembers(runs, Extents.create(16));
  const chains = new Set<string>();
  for (let i = 0; i < members.count; i++) chains.add(link.memberOf(members, i)?.chain.authAsymId ?? "");
  chains.delete("");
  return { chains: [...chains].sort(), exact };
}
