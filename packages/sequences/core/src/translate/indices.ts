/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Extra indices over a system (seqspec:libraries#translate), built by `Translator.index` on first use and cached per
// system. A host passes the data an index reads (its atoms, its labels, its SIFTS rows) as a function of the system;
// any other `IndexKind` is the host's own.

import type { CoordinateSystem } from "@molstar/seq-core/model";
import { firstOf } from "./maps.ts";
import type { IndexKind } from "./translator.ts";

/** Residue → atom ranges: atoms offsets[r] to offsets[r + 1] (each atom's residue index, grouped ascending; Mol*'s). */
const atomRanges = (name: string, residueOfAtom: (s: CoordinateSystem) => Int32Array): IndexKind<Int32Array> => ({
  name,
  build(s, t) {
    const atoms = residueOfAtom(s);
    const count = t.descriptorOf(s)?.stop ?? (atoms.length ? atoms[atoms.length - 1] + 1 : 0);
    const offsets = new Int32Array(count + 1);
    for (let a = 0; a < atoms.length; a++) offsets[atoms[a] + 1]++;
    for (let r = 0; r < count; r++) offsets[r + 1] += offsets[r];
    return offsets;
  },
});

/** Chain → units: a complex's members by their `structure.labelAsymId` (else their key), in member order. */
const chainUnits: IndexKind<ReadonlyMap<string, readonly CoordinateSystem[]>> = {
  name: "chain-units",
  build(s, t) {
    const chains = new Map<string, CoordinateSystem[]>();
    for (const key of t.descriptorOf(s)?.members ?? []) {
      const unit = t.system(key);
      const chain = t.descriptor(key)?.structure?.labelAsymId ?? key;
      const list = chains.get(chain);
      if (!unit) continue;
      if (list) list.push(unit);
      else chains.set(chain, [unit]);
    }
    return chains;
  },
};

/** Numbering label → index: the first position with each label (NO_LABEL: none). */
const labelIndex = (name: string, labelsOf: (s: CoordinateSystem) => Int32Array) =>
  ({ name, build: (s) => firstOf(labelsOf(s)) }) as IndexKind<ReadonlyMap<number, number>>;

/** A SIFTS segment table, (fromStart, toStart, length) triples sorted by fromStart: `TranslateMaps.bySegments`' input. */
const siftsTable = (name: string, rowsOf: (s: CoordinateSystem) => Int32Array): IndexKind<Int32Array> => ({
  name,
  build(s) {
    const rows = rowsOf(s);
    const order = Array.from({ length: Math.floor(rows.length / 3) }, (_, k) => 3 * k).sort(
      (a, b) => rows[a] - rows[b],
    );
    return Int32Array.from(order.flatMap((k) => [rows[k], rows[k + 1], rows[k + 2]]));
  },
});

export const IndexKinds = { atomRanges, chainUnits, labelIndex, siftsTable };
