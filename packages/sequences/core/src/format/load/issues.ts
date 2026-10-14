/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model reports under `model.*` codes that mirror the document's; the compiler reports them as `svs.*`
// (seqspec:model#alignments-and-trees). Codes with the same name map by prefix; the others by this table.

import type { Issue } from "@molstar/seq-core/model";
import type { StageContext } from "../contract.js";
import type { SvsRegistries } from "../registry.js";

const CODES: Readonly<Record<string, string>> = {
  "model.residues.invalid": "svs.sequence.alphabet",
  "model.ligand.element": "svs.params.type",
  "model.molecule.content": "svs.params.missing",
  "model.molecule.length": "svs.numbering.length",
  "model.numbering.name": "svs.ref.duplicate",
  "model.tree.shape": "svs.tree.leaves",
  "model.annotation.target": "svs.ref.derived_unknown",
  "model.annotation.geometry": "svs.params.type",
  "model.annotation.row": "svs.columns.type",
  "model.annotation.matrix": "svs.columns.type",
  "model.annotation.grid": "svs.columns.type",
  "model.annotation.states": "svs.params.missing",
  "model.annotation.link": "svs.params.type",
  "model.position.missing": "svs.params.missing",
  "model.mapping.endpoint": "svs.ref.derived_unknown",
  "model.mapping.numbering": "svs.numbering.unknown",
  "model.mapping.method": "svs.params.type",
  "model.key.duplicate": "svs.ref.duplicate",
  // A removal of a key the set does not hold: only the stream phase removes by key (ST1b-a; seqspec:format#validation).
  "model.key.unknown": "svs.stream.ref",
  "model.complex.unit": "svs.ref.derived_unknown",
  "model.unit.molecule": "svs.ref.derived_unknown",
  "model.unit.range": "svs.position.range",
  "model.unit.observed": "svs.position.range",
  "model.unit.origin": "svs.params.type",
  "model.provider.input": "svs.ref.derived_unknown",
};
/**
 * Model failures a document reads as a missing piece, warnings (seqspec:format#layer-types): the provider's node
 * unavailable.
 */
const MISSING_PIECE: ReadonlySet<string> = new Set([
  "model.provider.lazy",
  "model.provider.alphabet",
  "model.provider.size",
]);

/** The issue under its `svs.*` code; a `model.*` code the table lacks keeps its name under `svs.`. */
export function svsIssue(i: Issue): Issue {
  if (!i.code.startsWith("model.")) return i;
  const code = CODES[i.code] ?? `svs.${i.code.slice(6)}`;
  return MISSING_PIECE.has(i.code) ? { ...i, code, severity: "warning" } : { ...i, code };
}

/**
 * An issue, its severity the registered code's (an error when unregistered): every severity comes from the registry
 * (seqspec:format#validation), so a caller never writes one.
 */
export function issueOf(reg: SvsRegistries, code: string, path: string, message: string): Issue {
  return { code, severity: reg.get("issueCodes", code)?.severity ?? "error", path, message };
}

/** Reports an issue at `path`, its severity the registered code's. */
export function report(ctx: StageContext, code: string, path: string, message: string): void {
  ctx.report(issueOf(ctx.registries, code, path, message));
}

/** A context whose reports are mapped to `svs.*` and placed under `path`; everything else is `ctx`'s. */
export function at(ctx: StageContext, path: string): StageContext {
  const out: StageContext = Object.create(ctx);
  return Object.assign(out, { report: (i: Issue) => ctx.report({ ...svsIssue(i), path: path + i.path }) });
}
