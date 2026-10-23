/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Writes split-ligand.gro next to this file: a ligand of two residues (LG1: C1-C4 and O1, LG2: C1-C4; the bond
// LG1 C4 - LG2 C1 joins them), a water 3 Å from it and a sodium 2.3 Å from it, in a 2.5 nm cubic box, drifting 3 Å a
// frame along x over 3 frames. Every atom is wrapped into the box, so the bond between the residues crosses the x = 0
// face in the first frame. No bonds in a GRO file: only distances across the box join the residues. Deterministic:
// `node generate-ligand.mjs`.

import { writeFileSync } from 'node:fs';

const Frames = 3;
const L = 2.5;

/** [residue number, residue name, atom name, x, y, z], nm, before wrapping */
function atoms(t) {
  const out = [];
  const x0 = -0.49 + 0.3 * t;
  // a zigzag chain of 8 carbons, 1.49 Å bonds: LG1 has the first four (and O1 on C1), LG2 the last four
  for (let k = 0; k < 8; k++) {
    out.push([
      k < 4 ? 1 : 2,
      k < 4 ? 'LG1' : 'LG2',
      `C${(k % 4) + 1}`,
      x0 + 0.14 * k,
      1.25 + (k % 2 ? 0.025 : -0.025),
      1.25,
    ]);
    if (k === 0) out.push([1, 'LG1', 'O1', x0, 1.25 - 0.025 - 0.143, 1.25]);
  }
  // a water 3 Å above LG2 C2 and a sodium 2.3 Å below LG2 C4
  out.push([3, 'SOL', 'OW', x0 + 0.14 * 5, 1.25 + 0.3, 1.25]);
  out.push([3, 'SOL', 'HW1', x0 + 0.14 * 5 + 0.076, 1.25 + 0.3 + 0.059, 1.25]);
  out.push([3, 'SOL', 'HW2', x0 + 0.14 * 5 - 0.076, 1.25 + 0.3 + 0.059, 1.25]);
  out.push([4, 'NA', 'NA', x0 + 0.14 * 7, 1.25 - 0.23, 1.25]);
  return out;
}

const wrap = (v) => v - Math.floor(v / L) * L;
const f = (v, w, p) => v.toFixed(p).padStart(w);
let gro = '';
for (let t = 0; t < Frames; t++) {
  const list = atoms(t);
  gro += `Split ligand in a box, t= ${t.toFixed(5)}\n${String(list.length).padStart(5)}\n`;
  list.forEach(([rn, res, name, x, y, z], i) => {
    gro += `${String(rn).padStart(5)}${res.padEnd(5)}${name.padStart(5)}${String(i + 1).padStart(5)}${f(wrap(x), 8, 3)}${f(wrap(y), 8, 3)}${f(wrap(z), 8, 3)}\n`;
  });
  gro += `${f(L, 10, 5)}${f(L, 10, 5)}${f(L, 10, 5)}\n`;
}
writeFileSync(new URL('split-ligand.gro', import.meta.url), gro);
