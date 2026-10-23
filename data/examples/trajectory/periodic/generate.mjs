/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Writes split-peptide.gro next to this file: a four-residue peptide and three waters drifting through a rhombic
// dodecahedron (GROMACS's xy-square form, 3 nm), every atom wrapped into the box as an engine writes it, so the peptide
// is split across a face in most frames, the first among them. Deterministic: `node generate.mjs`.

import { writeFileSync } from 'node:fs';

const Frames = 8;
// the box vectors, nm, as columns
const A = [3, 0, 0];
const B = [0, 3, 0];
const C = [1.5, 1.5, 2.12132];

// one residue's backbone and CB, Å, along x; the next residue 3.75 Å on (the C–N bond 1.33 Å)
const Residue = [
  ['N', 0, 0, 0],
  ['CA', 1.21, 0.7, 0],
  ['C', 2.42, 0, 0],
  ['O', 2.42, -1.23, 0],
  ['CB', 1.21, 1.4, 1.2],
];
const Water = [
  ['OW', 0, 0, 0],
  ['HW1', 0.76, 0.59, 0],
  ['HW2', -0.76, 0.59, 0],
];

/** Atoms, nm: [residue number, residue name, atom name, x, y, z] */
function atoms(t) {
  const out = [];
  // the peptide starts across a face (between residues 2 and 3) and drifts along +x, turning slowly about z
  const angle = 0.15 * t;
  const [cos, sin] = [Math.cos(angle), Math.sin(angle)];
  const origin = [-0.05 + 0.45 * t, 1.2 + 0.1 * t, 0.9];
  for (let r = 0; r < 4; r++) {
    for (const [name, x, y, z] of Residue) {
      const px = (x + 3.75 * r) / 10;
      const py = y / 10;
      out.push([
        r + 1,
        'ALA',
        name,
        origin[0] + cos * px - sin * py,
        origin[1] + sin * px + cos * py,
        origin[2] + z / 10,
      ]);
    }
  }
  // waters, two of them near faces
  const waters = [
    [2.95, 0.3, 0.4],
    [1.6, 2.9, 1.2],
    [0.4, 1.9, 2.05],
  ];
  waters.forEach((w, i) => {
    for (const [name, x, y, z] of Water) {
      out.push([5 + i, 'SOL', name, w[0] + 0.05 * t + x / 10, w[1] - 0.03 * t + y / 10, w[2] + z / 10]);
    }
  });
  return out;
}

/** The point moved by a lattice vector into the box (fractional coordinates in [0, 1)), as an engine writes it. */
function wrap([x, y, z]) {
  // C has x and y; A and B are axis-aligned: solve the triangular system for the fractions
  let fc = z / C[2];
  let fb = (y - fc * C[1]) / B[1];
  let fa = (x - fc * C[0] - fb * B[0]) / A[0];
  fa -= Math.floor(fa);
  fb -= Math.floor(fb);
  fc -= Math.floor(fc);
  return [fa * A[0] + fb * B[0] + fc * C[0], fb * B[1] + fc * C[1], fc * C[2]];
}

const f = (v, w, p) => v.toFixed(p).padStart(w);
let gro = '';
for (let t = 0; t < Frames; t++) {
  const list = atoms(t);
  gro += `Split peptide in a dodecahedron, t= ${t.toFixed(5)}\n${String(list.length).padStart(5)}\n`;
  list.forEach(([rn, res, name, ...p], i) => {
    const [x, y, z] = wrap(p);
    gro += `${String(rn).padStart(5)}${res.padEnd(5)}${name.padStart(5)}${String(i + 1).padStart(5)}${f(x, 8, 3)}${f(y, 8, 3)}${f(z, 8, 3)}\n`;
  });
  // v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
  gro += [A[0], B[1], C[2], A[1], A[2], B[0], B[2], C[0], C[1]].map((v) => f(v, 10, 5)).join('') + '\n';
}
writeFileSync(new URL('split-peptide.gro', import.meta.url), gro);
