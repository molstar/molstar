/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small ensemble for tests of the trajectory and the animation (data/examples has no multi-model file).

/** An ensemble: one residue, three models, the atoms moving by a tenth of an angstrom in each. */
export function ensemble(models = 3): File {
  const lines: string[] = [];
  for (let m = 1; m <= models; m++) {
    lines.push(`MODEL     ${String(m).padStart(4)}`);
    const x = (m * 0.1).toFixed(3).padStart(8);
    lines.push(`ATOM      1  N   GLY A   1    ${x}   1.000   1.000  1.00  0.00           N  `);
    lines.push(`ATOM      2  CA  GLY A   1    ${x}   2.000   1.000  1.00  0.00           C  `);
    lines.push(`ATOM      3  C   GLY A   1    ${x}   2.000   2.000  1.00  0.00           C  `);
    lines.push('ENDMDL');
  }
  lines.push('END');
  return new File([lines.join('\n')], 'ensemble.pdb');
}
