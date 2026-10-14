/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:sifts_coverage` (seqspec:libraries#the-mol-adapter): which positions of a UniProt sequence each structure
// resolves, as one `ranges` annotation on that sequence in its `uniprot:<acc>` numbering. A row per resolved stretch of
// each polymer chain (its observed residues through the entity's SIFTS numbering, or position by position when the
// entity is the sequence itself, as an AlphaFold model's); a chain's rows are one item. Each row carries the JSON field
// `structure` (`StructureRef`: what "Show in 3D" loads), `entry` and `chain` (for filters and `split`) and `label`; the
// annotation carries the `custom.molstar` defaults.

import {
  bool,
  list,
  nullable,
  OptionalField,
  object,
  type ProviderEntry,
  RequiredField,
  ref,
  SimpleParamsSchema,
  str,
} from "@molstar/seq-core/format";
import { KeyKind, MoleculeKind, type NumberingInput, type SeqSet } from "@molstar/seq-core/model";
import type { StructureRef } from "./compose.js";
import { complexUnits } from "./seq-units.js";

export interface CoverageStructure {
  complex: string;
  source: { service: string; id: string };
  label?: string | null;
  assembly?: string | null;
}
export interface CoverageRowsOptions {
  target: string;
  accession: string;
  structures: CoverageStructure[];
  defaults?: Record<string, unknown> | null;
  /** Rows on each chain's own unit, in its indices (instance space), not on the target. */
  instance?: boolean;
}

/** [entity index start, UniProt number start, length] segments of the molecule onto `acc`; null: none. */
function segmentsOf(set: SeqSet, molecule: number, acc: string, target: number): [number, number, number][] | null {
  const m = set.molecules[molecule];
  if (!m) return null;
  const n = m.numberings.byName(set.strings.find(`uniprot:${acc}`));
  if (n?.method === "segments") {
    return Array.from(n.indexStart, (s, k): [number, number, number] => [s, n.labelStart[k], n.segmentLength[k]]);
  }
  if (n?.method === "offset") return [[0, n.first, m.length]];
  // The sequence itself (an AlphaFold model of the accession): equal letters, position by position.
  const t = set.molecules[target];
  if (t && t.length === m.length && Array.from({ length: m.length }).every((_, i) => t.code[i] === m.code[i])) {
    return [[0, 1, m.length]];
  }
  return null;
}

/** The rows `molstar:sifts_coverage` gives (exported for hosts that write the annotation by value). */
export function coverageRows(set: SeqSet, o: CoverageRowsOptions) {
  const t = { kind: 0, id: 0, row: -1 };
  const target = set.byName(o.target, t) && t.kind === KeyKind.molecule ? t.id : -1;
  const rows = {
    beg: [] as number[],
    end: [] as number[],
    id: [] as string[],
    label: [] as string[],
    chain: [] as string[],
    /** The chain's label_asym_id (`chain` is its author ID, what a person reads). */
    asym: [] as string[],
    /** The chain's unit and its row's [start, stop) there (`instance`). */
    unit: [] as string[],
    start: [] as number[],
    stop: [] as number[],
  };
  const entry: string[] = [];
  const structure: StructureRef[] = [];
  for (const s of o.structures) {
    const name = s.label ?? s.source.id.toUpperCase();
    for (const u of complexUnits(set, s.complex)) {
      if (!set.byName(u.ref, t) || t.kind !== KeyKind.unit || !u.mvs.label_asym_id) continue;
      const mol = set.units.molecule[t.id];
      if (set.molecules[mol]?.kind !== MoleculeKind.polymer) continue;
      const segs = segmentsOf(set, mol, o.accession, target);
      if (!segs) continue;
      const obs = set.units.observed[t.id] ?? Int32Array.of(set.units.rangeStart[t.id], set.units.rangeStop[t.id]);
      const asym = u.mvs.label_asym_id;
      const chain = u.mvs.auth_asym_id ?? asym;
      const ref: StructureRef = {
        source: s.source,
        chains: [asym],
        ...(chain !== asym ? { authChains: [chain] } : {}),
        ...(s.assembly != null ? { assembly: s.assembly } : {}),
        ...(s.complex !== s.source.id.toLowerCase() ? { complex: s.complex } : {}),
      };
      for (let k = 0; k + 1 < obs.length; k += 2) {
        for (const [i0, l0, len] of segs) {
          const a = Math.max(obs[k], i0);
          const b = Math.min(obs[k + 1], i0 + len);
          if (a >= b) continue;
          rows.beg.push(l0 + a - i0);
          rows.end.push(l0 + b - 1 - i0);
          rows.id.push(`${s.complex}/${asym}`);
          // The author's chain ID, its label_asym_id beside it where they differ (the hover card's title).
          rows.label.push(`${name} chain ${chain}${chain !== asym ? ` (label_asym_id ${asym})` : ""}`);
          rows.chain.push(chain);
          rows.asym.push(asym);
          rows.unit.push(u.ref);
          rows.start.push(a);
          rows.stop.push(b);
          entry.push(name);
          structure.push(ref);
        }
      }
    }
  }
  return { ...rows, entry, structure };
}

const coverageOptions = SimpleParamsSchema({
  target: RequiredField(ref("sequence"), "the UniProt sequence the rows sit on"),
  accession: RequiredField(str, "its accession: rows are in its uniprot:<acc> numbering"),
  structures: RequiredField(
    list(
      object({
        complex: RequiredField(ref("complex"), "the structure's complex"),
        source: RequiredField(
          object({ service: RequiredField(str, "pdbe or afdb"), id: RequiredField(str, "the entry") }),
          "what Show in 3D loads",
        ),
        label: OptionalField(nullable(str), null, "its name in labels; null: the id in upper case"),
        assembly: OptionalField(nullable(str), null, "the assembly Show in 3D shows; null: the deposited model"),
      }),
    ),
    "the structures, in row order",
  ),
  defaults: OptionalField(nullable(object()), null, "custom.molstar: the 3D defaults (representation, base color)"),
  instance: OptionalField(bool, false, "rows on each chain's own unit (instance space), not the target's"),
});

export const SiftsCoverageProvider: ProviderEntry = {
  name: "molstar:sifts_coverage",
  options: coverageOptions,
  geometry: "ranges",
  async compute(ctx, options) {
    const o = options as CoverageRowsOptions;
    if (!ctx.set) return null;
    const r = coverageRows(ctx.set, o);
    const numbering: NumberingInput["name"] = `uniprot:${o.accession}`;
    const fields = {
      id: r.id,
      label: r.label,
      entry: r.entry,
      chain: r.chain,
      asym: r.asym,
      structure: { json: r.structure },
    };
    // Instance space: each row on its chain's unit, placed on the target through the mappings.
    if (o.instance) {
      const columns = { target: r.unit, start: r.start, stop: r.stop, ...fields };
      return { key: "", geometry: "ranges", label: "Chains", columns, custom: { molstar: o.defaults ?? {} } };
    }
    return {
      key: "",
      geometry: "ranges",
      label: "PDB coverage",
      target: o.target,
      numbering,
      columns: { beg: r.beg, end: r.end, ...fields },
      custom: { molstar: o.defaults ?? {} },
    };
  },
};

/**
 * Coverage depth: runs of positions of the target covered by the same structures, each with how many (`value`), a label
 * and a tooltip naming them ("covered by 3 structures: 1TUP, 2OCJ, 2XWR"), and the JSON field `structures` listing them
 * (so a track over it offers each in 3D and colors none); positions none covers have no row. From the coverage rows;
 * positions in the target's `uniprot:<acc>` numbering, `beg` and `end` inclusive.
 */
export function coverageDepth(rows: ReturnType<typeof coverageRows>) {
  const at = new Map<number, string[]>();
  const refs = new Map<string, StructureRef>();
  for (let k = 0; k < rows.beg.length; k++) {
    refs.set(rows.entry[k], { ...rows.structure[k], chains: undefined });
    for (let p = rows.beg[k]; p <= rows.end[k]; p++) {
      const names = at.get(p) ?? [];
      if (!names.includes(rows.entry[k])) names.push(rows.entry[k]);
      at.set(p, names);
    }
  }
  const out = {
    beg: [] as number[],
    end: [] as number[],
    value: [] as number[],
    label: [] as string[],
    tooltip: [] as string[],
    structures: [] as StructureRef[][],
  };
  const said = (n: readonly string[]) =>
    `covered by ${n.length} structure${n.length === 1 ? "" : "s"}: ${n.join(", ")}`;
  const positions = [...at.keys()].sort((x, y) => x - y);
  for (let k = 0; k < positions.length; ) {
    const names = (at.get(positions[k]) ?? []).slice().sort();
    let j = k + 1;
    while (
      j < positions.length &&
      positions[j] === positions[j - 1] + 1 &&
      (at.get(positions[j]) ?? []).slice().sort().join() === names.join()
    )
      j++;
    out.beg.push(positions[k]);
    out.end.push(positions[j - 1]);
    out.value.push(names.length);
    out.label.push(`${names.length} structure${names.length === 1 ? "" : "s"}`);
    out.tooltip.push(said(names));
    out.structures.push(names.map((n) => refs.get(n) as StructureRef));
    k = j;
  }
  return out;
}

/** `molstar:sifts_depth`: the same options; a `ranges` annotation, the coverage's depth in runs (above). */
export const SiftsDepthProvider: ProviderEntry = {
  name: "molstar:sifts_depth",
  options: coverageOptions,
  geometry: "ranges",
  async compute(ctx, options) {
    const o = options as CoverageRowsOptions;
    if (!ctx.set) return null;
    const d = coverageDepth(coverageRows(ctx.set, o));
    return {
      key: "",
      geometry: "ranges",
      label: "PDB coverage depth",
      target: o.target,
      numbering: `uniprot:${o.accession}`,
      columns: {
        beg: d.beg,
        end: d.end,
        value: d.value,
        label: d.label,
        tooltip: d.tooltip,
        structures: { json: d.structures },
      },
    };
  },
};
