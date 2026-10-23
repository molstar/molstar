/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Writes long-chain.lammpstrj next to this file: a LAMMPS dump (`id mol type x y z ix iy iz`) of a 12-bead chain,
// 17.6 Å long, in a 20 Å box from -10 to 10, drifting 6 Å a frame along x across its faces over 4 frames. Every atom
// is wrapped into the box with its image flags, as LAMMPS writes it. The atom ids interleave the chain (beads 0, 7,
// 1, 8, …), so atoms with consecutive ids are up to 11.2 Å apart, more than half the box: only the image flags make the
// chain whole. Deterministic: `node generate-lammps.mjs`.

import { writeFileSync } from 'node:fs';

const Frames = 4;
const Lo = -10;
const L = 20;
const Spacing = 1.6;
// the chain's beads in id order
const IdOrder = [0, 7, 1, 8, 2, 9, 3, 10, 4, 11, 5, 6];

let dump = '';
for (let t = 0; t < Frames; t++) {
  dump += `ITEM: TIMESTEP\n${t * 100}\nITEM: NUMBER OF ATOMS\n${IdOrder.length}\n`;
  dump += `ITEM: BOX BOUNDS pp pp pp\n${Lo} ${Lo + L}\n${Lo} ${Lo + L}\n${Lo} ${Lo + L}\n`;
  dump += 'ITEM: ATOMS id mol type x y z ix iy iz\n';
  IdOrder.forEach((bead, i) => {
    // along x, a zigzag in y; the chain starts at x = -12 (across the lower face) and drifts +6 Å a frame
    const p = [-12 + 6 * t + bead * Spacing, bead % 2 ? 0.5 : -0.5, 2];
    const image = p.map((v) => Math.floor((v - Lo) / L));
    const wrapped = p.map((v, k) => v - image[k] * L);
    dump += `${i + 1} 1 1 ${wrapped.map((v) => v.toFixed(4)).join(' ')} ${image.join(' ')}\n`;
  });
}
writeFileSync(new URL('long-chain.lammpstrj', import.meta.url), dump);
