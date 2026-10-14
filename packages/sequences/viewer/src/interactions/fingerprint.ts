/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fingerprints (seqspec:libraries#interactions): ligands or poses (rows) × contacted positions in a shared frame
// (columns), a KindMask per cell, as KLIFS shows them. Built at the edges (once per provider run); `diff` allocates
// nothing.

import {
  type Annotation,
  type ComposedPath,
  type CoordinateSystem,
  CoordinateSystemKind,
  Issue,
  Mappings,
  MoleculeKind,
  type PathStep,
  Paths,
  type SeqSet,
  SortedArray,
  type UnitId,
} from "@molstar/seq-core/model";
import { Interactions } from "./queries.js";
import { KindMask } from "./vocabulary.js";

export interface Fingerprint {
  /** UnitId per ligand or pose. */
  readonly rows: Int32Array;
  /** The columns' coordinate system: a unit, a sequence (its molecule), or an alignment's columns. */
  readonly system: CoordinateSystem;
  /** The contacted positions of the system. */
  readonly columns: SortedArray;
  /** rows × columns, row-major: the kinds per cell. */
  readonly cells: Uint16Array;
}

export interface BuildOptions {
  /** Only these kinds (all by default). */
  kinds?: KindMask;
  /** Receives `interactions.fingerprint.unmapped`, one warning per unit with residues outside the system. */
  issues?: Issue[];
}

/**
 * Each ligand's contacts in every pairs annotation, each contacted polymer residue mapped into the system through
 * `Mappings.path` (an ambiguous path takes its first); whole-unit endpoints and partners that are not polymers (a
 * metal ion, a cofactor, a water) are left out quietly.
 */
function build(
  set: SeqSet,
  pairs: readonly Annotation[],
  ligands: Int32Array,
  system: CoordinateSystem,
  options: BuildOptions = {},
): Fingerprint {
  const to: CoordinateSystem = { kind: system.kind, id: system.id, row: system.row };
  const steps: PathStep[] = [];
  const forward = new Map<number, ComposedPath | null>();
  const unmapped = new Map<number, number>();
  const hits: number[] = [];
  const toSystem = (unit: number) => {
    if (!forward.has(unit)) {
      const from = { kind: CoordinateSystemKind.unit, id: unit, row: -1 };
      const n = Mappings.path(set, from, to, steps);
      forward.set(unit, n === -1 ? null : Mappings.compose(set, from, steps));
    }
    return forward.get(unit);
  };
  ligands.forEach((ligand, r) => {
    for (const a of pairs) {
      Interactions.contactsOf(
        set,
        a,
        ligand as UnitId,
        options.kinds ?? KindMask.all,
        (unit, index, mask) => {
          if (index < 0 || set.molecules[set.units.molecule[unit]]?.kind !== MoleculeKind.polymer) return;
          const path = toSystem(unit);
          const at = path ? Paths.at(path, index) : -1;
          if (at >= 0) hits.push(r, at, mask);
          else unmapped.set(unit, (unmapped.get(unit) ?? 0) + 1);
        },
        null,
      );
    }
  });
  const columns = SortedArray.ofUnsortedArray(Int32Array.from(hits.filter((_, k) => k % 3 === 1)));
  const C = columns.length;
  const cells = new Uint16Array(ligands.length * C);
  for (let k = 0; k < hits.length; k += 3)
    cells[hits[k] * C + SortedArray.indexOf(columns, hits[k + 1])] |= hits[k + 2];
  // One issue for the build, naming the units (GV-07: one per unit).
  if (unmapped.size > 0) {
    const named = [...unmapped].map(([unit, n]) => `${n} of "${set.strings.get(set.units.key[unit])}"`);
    const more = named.length > 5 ? ` and ${named.length - 5} more units` : "";
    const message = `contacted positions with no path into the system, left out: ${named.slice(0, 5).join(", ")}${more}.`;
    options.issues?.push(Issue.warning("interactions.fingerprint.unmapped", "", message));
  }
  return { rows: Int32Array.from(ligands), system: to, columns, cells };
}

/** The kinds row `b` lost and gained against row `a`, per column that differs; `column` is the frame position. */
function diff<C>(
  f: Fingerprint,
  a: number,
  b: number,
  visit: (column: number, lost: KindMask, gained: KindMask, ctx: C) => void,
  ctx: C,
): void {
  const n = f.columns.length;
  for (let c = 0; c < n; c++) {
    const x = f.cells[a * n + c];
    const y = f.cells[b * n + c];
    if (x !== y) visit(f.columns[c], x & ~y, y & ~x, ctx);
  }
}

export const Fingerprints = { build, diff };
