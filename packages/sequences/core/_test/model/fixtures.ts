/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SetInputs in the shape `molstar:structure` yields (seqspec:libraries#the-mol-adapter), written by hand as small
// inline fixtures. They are not recorded provider outputs: those come with V2 (lib/seqstar/core/_test/fixtures/).
//
// 1BRS, from the Seq* fixture (lib/seqstar/core/_test/fixtures/1brs.json, UniProtKB CC BY 4.0,
// wwPDB CC0): entity 1 is barnase's mature chain (P00648 48–157), entity 2 barstar without its initiator
// Met (P11540 2–90) with 1BRS's C40A/C82A. Observed residues are the fixture's for chains A and D, used for
// all three copies; auth numbering equals label numbering in 1BRS. Waters are left out (the default).
//
// 1A3N: hemoglobin α (HBA_HUMAN mature, 141) in chains A and C, β (HBB_HUMAN mature, 146) in B and D, and a
// heme (HEM, 43 heavy atoms, atom names as in the CCD) per chain, label_asym_ids E–H.

import type { SetInput, UnitInput } from "../../src/model/index.ts";

const BARNASE =
  "AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR";
const BARSTAR = "KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS";

function polymerUnit(key: string, molecule: string, length: number, observed?: [number, number][]): UnitInput {
  const auth = { name: "auth", method: "segments" as const, segments: [[0, 1, length]] as [number, number, number][] };
  return { key, molecule, observed, origin: "experimental", numberings: [{ ...auth, isDefault: true }] };
}

export function brs1(): SetInput {
  const label = { name: "label", method: "offset" as const, first: 1 };
  const barnase = (c: string) => polymerUnit(`1brs/${c}`, "1brs/entity:1", 110, [[2, 110]]);
  const barstar = (c: string) =>
    polymerUnit(`1brs/${c}`, "1brs/entity:2", 89, [
      [0, 63],
      [65, 89],
    ]);
  return {
    molecules: [
      { key: "1brs/entity:1", residues: BARNASE, label: "Barnase", numberings: [label] },
      { key: "1brs/entity:2", residues: BARSTAR, label: "Barstar", numberings: [label] },
    ],
    units: [barnase("A"), barnase("B"), barnase("C"), barstar("D"), barstar("E"), barstar("F")],
    complexes: [{ key: "1brs", units: ["A", "B", "C", "D", "E", "F"].map((c) => `1brs/${c}`) }],
  };
}

const HBA =
  "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
const HBB =
  "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
export const HEM_ATOMS = (
  "FE NA NB NC ND CHA CHB CHC CHD C1A C2A C3A C4A CMA CAA CBA CGA O1A O2A C1B C2B C3B C4B CMB CAB CBB " +
  "C1C C2C C3C C4C CMC CAC CBC C1D C2D C3D C4D CMD CAD CBD CGD O1D O2D"
).split(" ");

export function a3n1(): SetInput {
  const heme = (c: string): UnitInput => ({
    key: `1a3n/${c}`,
    molecule: "1a3n/entity:3",
    origin: "experimental",
    numberings: [{ name: "atom", method: "names", names: HEM_ATOMS, isDefault: true }],
  });
  return {
    molecules: [
      { key: "1a3n/entity:1", residues: HBA, label: "Hemoglobin subunit alpha" },
      { key: "1a3n/entity:2", residues: HBB, label: "Hemoglobin subunit beta" },
      {
        key: "1a3n/entity:3",
        kind: "ligand",
        elements: HEM_ATOMS.map((a) => (a === "FE" ? "Fe" : a[0])),
        chemCompIds: ["HEM"],
        label: "Protoporphyrin IX containing Fe",
      },
    ],
    // Entity order, then chain order.
    units: [
      polymerUnit("1a3n/A", "1a3n/entity:1", 141),
      polymerUnit("1a3n/C", "1a3n/entity:1", 141),
      polymerUnit("1a3n/B", "1a3n/entity:2", 146),
      polymerUnit("1a3n/D", "1a3n/entity:2", 146),
      ...["E", "F", "G", "H"].map(heme),
    ],
    complexes: [{ key: "1a3n", units: ["A", "C", "B", "D", "E", "F", "G", "H"].map((c) => `1a3n/${c}`) }],
  };
}
