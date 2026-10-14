/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// 1IEP from the Seq* fixture (lib/seqstar/core/_test/fixtures/1iep.json, from AutoDock Vina's
// basic-docking tutorial, Apache-2.0): ABL kinase chain A (auth 225–498) with imatinib (STI, 37 heavy atoms),
// and per residue the contacts of the crystal pose and of Vina's four docked poses. Here the crystal STI is
// unit "1iep/C" of its own molecule, the poses four units of another (one copy group), atoms named by
// element and ordinal ("C1", "N5"): a `names` numbering "atom" on each ligand molecule.

import { readFileSync } from "node:fs";
import type { AnnotationInput, SetInput } from "../../src/model/index.ts";

interface Pose {
  label: string;
  score: number | null;
  /** [auth residue, kind, atom index, distance]. */
  contacts: [number, string, number, number][];
}
export const IEP = JSON.parse(readFileSync(new URL("../fixtures/1iep.json", import.meta.url), "utf8")) as {
  first: number;
  sequence: string;
  ligand: { code: string; atoms: string[] };
  poses: Pose[];
};

export const ATOM_NAMES = IEP.ligand.atoms.map((e, i) => `${e}${i + 1}`);
export const POSE_UNITS = [0, 1, 2, 3].map((k) => `poses/${k}/ligand`);
/** The ligand unit of each fixture pose: the crystal's, then Vina's. */
export const LIGAND_UNITS = ["1iep/C", ...POSE_UNITS];

export function iep1(): SetInput {
  const ligand = (key: string) => ({
    key,
    kind: "ligand" as const,
    elements: IEP.ligand.atoms,
    chemCompIds: [IEP.ligand.code],
    numberings: [{ name: "atom", method: "names" as const, names: ATOM_NAMES, isDefault: true }],
  });
  return {
    molecules: [
      {
        key: "1iep/entity:1",
        residues: IEP.sequence,
        numberings: [{ name: "auth", method: "offset", first: IEP.first, isDefault: true }],
      },
      ligand("1iep/STI"),
      ligand("poses/STI"),
    ],
    units: [
      { key: "1iep/A", molecule: "1iep/entity:1", origin: "experimental" },
      { key: "1iep/C", molecule: "1iep/STI", origin: "experimental" },
      ...IEP.poses
        .slice(1)
        .map((p, k) => ({ key: POSE_UNITS[k], molecule: "poses/STI", origin: "docked", score: p.score ?? undefined })),
    ],
  };
}

/** Every pose's contacts as one pairs table, labels on both sides: residues by `auth`, atoms by name. */
export function contacts(): AnnotationInput {
  const rows = IEP.poses.flatMap((p, k) => p.contacts.map((c) => ({ c, unit: LIGAND_UNITS[k] })));
  return {
    key: "contacts",
    geometry: "pairs",
    semanticType: "interactions:interaction",
    target: "1iep/A",
    numbering: "auth",
    bNumbering: "atom",
    columns: {
      beg: rows.map((r) => r.c[0]),
      b_target: rows.map((r) => r.unit),
      b_beg: rows.map((r) => ATOM_NAMES[r.c[2]]),
      category: rows.map((r) => r.c[1]),
      distance: rows.map((r) => r.c[3]),
    },
  };
}
