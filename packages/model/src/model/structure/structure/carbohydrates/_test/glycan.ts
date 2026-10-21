/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A synthetic polysaccharide as mmCIF `atom_site`, for carbohydrate tests and benchmarks: strands of hexopyranoses
// (C1-C5, O5 ring, O4, C6, O6) in one chain, as MD files without chain ids and polysaccharide entries give. Rings are
// planar hexagons with bonds of 1.43-1.52 Å, so bonds come from distances. Residue i's C1 bonds to residue i+1's O4
// (1-4 links: rings 2 bonds apart); every `branchEvery`-th residue carries a residue on its O6 (1-6 links: 3 bonds).

export interface GlycanParams {
  /** Strands in chain A, 6 Å apart in z. */
  strands: number;
  /** Residues per strand. */
  length: number;
  /** A 1-6 branch on every `branchEvery`-th residue of a strand (from the second); 0 for none. */
  branchEvery: number;
  /** The first branch of strand 0 in its own chain C (an inter-unit link), and an ASN in chain B on strand 0's last C1. */
  otherChains?: boolean;
  /** A carbon (C7) between a branch point's O6 and the branch's C1: the rings are then 4 bonds apart, not linked. */
  spacer?: boolean;
  /** Every `altEvery`-th residue of a strand (from the second) in two conformations, alt ids A and B (0.2 Å apart). */
  altEvery?: number;
}

const Comps = ['GLC', 'BGC', 'MAN', 'GAL'];
const RingR = 1.45;

function dir(deg: number): [number, number] {
  const a = (deg * Math.PI) / 180;
  return [Math.cos(a), Math.sin(a)];
}

export function glycanCif(p: GlycanParams): string {
  const rows: string[] = [];
  let id = 0;
  let alt = '.';
  const atom = (
    symbol: string,
    name: string,
    comp: string,
    asym: string,
    entity: number,
    seq: number,
    x: number,
    y: number,
    z: number,
  ) => {
    rows.push(
      `HETATM ${++id} ${symbol} ${name} ${alt} ${comp} ${asym} ${entity} . ${seq} ${asym} ${x.toFixed(3)} ${y.toFixed(3)} ${z.toFixed(3)} 1 10 1`,
    );
  };
  /** A residue centred at (cx, cy, z), rotated by `rot` degrees: C1 at `rot`, C4 at `rot` + 180, C6 and O6 at + 240. */
  const sugar = (comp: string, asym: string, seq: number, cx: number, cy: number, z: number, rot: number) => {
    const ring = ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'];
    for (let k = 0; k < 6; k++) {
      const [dx, dy] = dir(rot + 60 * k);
      atom(ring[k][0], ring[k], comp, asym, 1, seq, cx + RingR * dx, cy + RingR * dy, z);
    }
    const [ax, ay] = dir(rot + 180);
    atom('O', 'O4', comp, asym, 1, seq, cx + (RingR + 1.43) * ax, cy + (RingR + 1.43) * ay, z);
    const [bx, by] = dir(rot + 240);
    atom('C', 'C6', comp, asym, 1, seq, cx + (RingR + 1.52) * bx, cy + (RingR + 1.52) * by, z);
    atom('O', 'O6', comp, asym, 1, seq, cx + (RingR + 2.95) * bx, cy + (RingR + 2.95) * by, z);
  };

  // C1 of residue i to O4 of residue i + 1: 1.43 Å along x
  const step = 2 * RingR + 2 * 1.43;
  // a branch's ring centre from its parent's: O6 at 4.40 Å, then 1.43 Å to the branch's C1, rotated to face it
  const [bx, by] = dir(240);
  const branch = RingR + 2.95 + 1.43 + RingR + (p.spacer ? 1.43 : 0);

  const deferred: (() => void)[] = [];
  for (let s = 0; s < p.strands; s++) {
    const z = 6 * s;
    let seq = 1000 * s;
    let firstBranch = true;
    for (let r = 0; r < p.length; r++) {
      const cx = r * step;
      if (p.altEvery && r % p.altEvery === 1) {
        ++seq;
        alt = 'A';
        sugar(Comps[r % Comps.length], 'A', seq, cx, 0, z, 0);
        alt = 'B';
        sugar(Comps[r % Comps.length], 'A', seq, cx, 0, z + 0.2, 0);
        alt = '.';
      } else {
        sugar(Comps[r % Comps.length], 'A', ++seq, cx, 0, z, 0);
      }
      if (p.branchEvery > 0 && r > 0 && r % p.branchEvery === 0) {
        if (p.spacer) {
          const d = RingR + 2.95 + 1.43;
          atom('C', 'C7', Comps[r % Comps.length], 'A', 1, seq, cx + d * bx, d * by, z);
        }
        const comp = Comps[(r + 1) % Comps.length];
        if (p.otherChains && s === 0 && firstBranch) {
          deferred.push(() => sugar(comp, 'C', 1, cx + branch * bx, branch * by, z, 60));
        } else {
          sugar(comp, 'A', ++seq, cx + branch * bx, branch * by, z, 60);
        }
        firstBranch = false;
      }
    }
    if (p.otherChains && s === 0) {
      const x = (p.length - 1) * step + RingR;
      deferred.unshift(() => {
        atom('N', 'ND2', 'ASN', 'B', 2, 1, x + 1.45, 0, z);
        atom('C', 'CG', 'ASN', 'B', 2, 1, x + 2.95, 0, z);
      });
    }
  }
  for (const f of deferred) f();

  return `data_glycan
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
${rows.join('\n')}
`;
}
