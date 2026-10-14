/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The format's acceptance cases for the Seq* stories: the format's
// test cases and fixture loaders imported as they are (their `node:fs` is main.ts's bundled-fixture module),
// compiled once per case.

import { compile, type DocNode, register, type SeqViewSpecDoc, validate } from "@molstar/seq-core/format";
import type { Issue, SeqSet } from "@molstar/seq-core/model";
import { InteractionsRegistry } from "@molstar/seq-viewer/interactions";
import { CASES } from "../../../../core/_test/format/cases.js";
import CASES_SOURCE from "../../../../core/_test/format/cases.ts?raw";
import { fixtureLoaders, registries } from "../../../../core/_test/format/load-fixtures.js";

/** The fixtures' registries and the interactions extension's format half (case 6's fingerprint). */
export const reg = registries();
register(reg, InteractionsRegistry);

export const CASE_TITLES: Readonly<Record<string, string>> = {
  case1: "Case 1: p53: text and tracks",
  case2: "Case 2: Antibody variable domains, Kabat and IMGT",
  case3: "Case 3: A protein family, its subgroups and its members' features",
  case4: "Case 4: An RNA family with base pairs",
  case5: "Case 5: A DNA construct",
  case6: "Case 6: Imatinib in ABL",
  case7: "Case 7: Barnase–barstar (1BRS)",
  case8: "Case 8: An AlphaFold DB complex, each chain from UniProt",
  case9: "Case 9: A family alignment loaded by view (PF00042)",
  case10: "Case 10: The Stories drafts",
  case12: "Case 12: A designed dimer by value",
  case13: "Case 13: A protein's features, UniProt-style",
  case14: "Case 14: A protein over its coding DNA",
  case15: "Case 15: Indels: protein designs in the wild type's numbering",
  case16: "Case 16: A plasmid map: a ring, its labels and chords, its linear detail and window",
  case17: "Case 17: A genotype matrix: from a pixel row of many samples to a cell each",
  case18: "Case 18: Finding things in p53: queries, search, bookmarks",
  case19: "Case 19: A chromosome by region: streamed, under its ideogram",
  case20: "Case 20: A protein on its genome: streamed through codons, an item coordinate system",
  case21: "Case 21: A gene on its reverse strand: its own direction, splice junctions, the mRNA",
  case22: "Case 22: A family triaged: a guide tree, keyed sets on rows, the vertical overview",
};

const FOCUS_T315 = "!seq-focus=1iep/A&beg=315&end=315&numbering=auth";
/** Descriptions with interactive-markdown links for the cases that tell something. */
const DESCRIPTIONS: Record<string, { text: string; snapshots?: Record<string, string> }> = {
  case1: {
    text: "p53's [DNA-binding domain](!seq-focus=p53&start=101&stop=292) carries [the cancer hotspots](!seq-highlight=hotspots); [four residues](!seq-select=zinc) hold its zinc. The same data [as tracks](!seq-view=tracks).",
  },
  case6: {
    text: `[Imatinib](!seq-highlight=1iep/C) holds ABL's kinase domain in its inactive form; [the gatekeeper, T315](${FOCUS_T315}) guards the back pocket. [Its contacts](!seq-select=contacts) are within 4 Å.`,
  },
  case7: {
    text: "[Barnase](!seq-highlight=barnase) and [barstar](!seq-highlight=barstar) bind through [the interface](!seq-select=iface). Compare [each chain](!seq-view=complex-each) of the asymmetric unit.",
  },
  case10: {
    text: `The drafts' story: [the kinase](!seq-highlight=1iep/A), [imatinib](!seq-highlight=1iep/C) and [the gatekeeper](${FOCUS_T315}). Start with [the first snapshot](!seq-snapshot=s1).`,
    snapshots: {
      s1: "ABL1's kinase domain ([chain A](!seq-highlight=1iep/A)) holds [imatinib](!seq-highlight=1iep/C). [Next](!seq-snapshot=next).",
      s2: "Only [the pocket](!seq-select=pocket): the residues within 4 Å of the drug. [Next](!seq-snapshot=next).",
      s3: `[Threonine 315](${FOCUS_T315}&section=strip.main) guards the back of the pocket; mutated to isoleucine, the drug no longer fits. [Back to the start](!seq-snapshot=s1).`,
    },
  },
};
export function described(name: string, doc: SeqViewSpecDoc): SeqViewSpecDoc {
  const d = DESCRIPTIONS[name];
  if (!d) return doc;
  return {
    ...doc,
    metadata: { ...doc.metadata, description: d.text },
    snapshots: doc.snapshots?.map((x) => ({ ...x, description: d.snapshots?.[x.ref] ?? x.description })),
  };
}

/**
 * The document with the sections of its text views (`layout: "text"`) drawn in the built-in preset `style` (`wrapped`,
 * `compact`).
 */
export function textStyled(doc: SeqViewSpecDoc, style: string): SeqViewSpecDoc {
  const walk = (n: DocNode, text: boolean): DocNode => ({
    ...n,
    params: text && n.kind === "section" ? { ...n.params, style } : n.params,
    children: n.children?.map((c) => walk(c, text || (n.kind === "view" && n.params?.layout === "text"))),
  });
  return { ...doc, root: walk(doc.root as DocNode, false) };
}

/** Case 13 with `split.empty: "hide"` on its split layers: the variants grid keeps only the residues a variant has. */
export function hideEmpty(n: DocNode): DocNode {
  const sp = n.params?.split as object | null | undefined;
  const params = sp ? { ...n.params, split: { ...sp, empty: "hide" } } : n.params;
  return { ...n, params, children: n.children?.map(hideEmpty) };
}

export interface CompiledCase {
  readonly doc: SeqViewSpecDoc;
  readonly set: SeqSet;
  readonly validation: readonly Issue[];
  readonly compile: readonly Issue[];
  readonly views: readonly { ref: string; label: string }[];
}

const cache = new Map<string, Promise<CompiledCase>>();

/** A case compiled on the fixtures; the same promise every time (for React's `use`). */
export function compiledCase(name: string): Promise<CompiledCase> {
  let p = cache.get(name);
  if (!p) {
    const doc = CASES[name]();
    p = compile(doc, reg, fixtureLoaders()).then(({ set, issues }) => ({
      doc,
      set,
      validation: validate(doc, reg),
      compile: issues,
      views: (doc.root.children ?? [])
        .filter((c) => c.kind === "view")
        .map((c) => ({ ref: c.ref as string, label: String(c.params?.label ?? c.ref) })),
    }));
    cache.set(name, p);
  }
  return p;
}

/** JSON for reading: long arrays and strings cut, typed arrays as lists. */
export function preview(value: unknown): string {
  return JSON.stringify(
    value,
    (_, v) => {
      const list = ArrayBuffer.isView(v) ? Array.from(v as unknown as ArrayLike<number>) : v;
      if (Array.isArray(list) && list.length > 12) return [...list.slice(0, 6), `… ${list.length - 6} more`];
      if (typeof list === "string" && list.length > 160) return `${list.slice(0, 120)}… (${list.length} chars)`;
      return list;
    },
    2,
  );
}

/**
 * A case's builder function as cases.ts has it (between `// #region <name>` and `// #endregion`), so the code shown is
 * the code the story runs; and what it uses from outside the function: cases.ts's own helpers and constants, and the
 * fixtures (load-fixtures.ts's imports).
 */
export function caseSource(name: string): { code: string; helpers: string[]; fixtures: string[] } {
  const marker = `// #region ${name}\n`;
  const start = CASES_SOURCE.indexOf(marker);
  const end = CASES_SOURCE.indexOf("// #endregion", start);
  if (start < 0 || end < 0) return { code: "", helpers: [], fixtures: [] };
  const code = CASES_SOURCE.slice(start + marker.length, end).trimEnd();
  const outside = CASES_SOURCE.slice(0, start) + CASES_SOURCE.slice(end);
  const declared = [...outside.matchAll(/^(?:export )?(?:const|function) (?:(\w+)|\[([^\]]+)\])/gm)].flatMap((m) =>
    (m[1] ?? m[2]).split(",").map((x) => x.trim()),
  );
  const imported = /import \{([^}]*)\} from "\.\/load-fixtures\.ts"/.exec(CASES_SOURCE)?.[1] ?? "";
  const used = (names: string[]) => names.filter((n) => n && new RegExp(`\\b${n}\\b`).test(code));
  return {
    code,
    helpers: used(declared.filter((n) => !/^case\d+$/.test(n) && n !== "CASES")),
    fixtures: used(imported.split(",").map((x) => x.trim())),
  };
}
