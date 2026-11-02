/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/util';
import { chemCompHash } from '../perception/chemcomp-cache.js';
import { getInterUnitNeighbours } from '../perception/inter-unit-neighbours.js';
import { perceiveBondOrders } from '../perception/perceive.js';
import { bondRows, ccdBondRows, compareBonds, modelFromPdb, perceived, readFixture, unitOf } from './util.js';

// Hard cases with their fixtures (_spec §4): perceived orders against the component's CCD definition.

describe('bond-order perception: hard cases', () => {
  it('1nco-chr: neocarzinostatin chromophore (double, triple, aromatic, non-aromatic 5-ring) matches its CCD', async () => {
    const model = await modelFromPdb(readFixture('1nco-chr.pdb'), '1NCO');
    const unit = unitOf(perceived(model), 'CHR');
    const expected = await ccdBondRows(readFixture('ccd-chr.cif'));
    const { missed, spurious } = compareBonds(bondRows(unit), expected);
    expect({ missed, spurious }).toEqual({ missed: [], spurious: [] });
  });

  it('5kre-ser124-6wg: the link C12-SER.OG is seen from the ligand unit, and in its pattern key', async () => {
    const model = await modelFromPdb(readFixture('5kre-ser124-6wg.pdb'), '5KRE');
    const s = perceived(model);
    const ligand = unitOf(s, '6WG');
    const serine = unitOf(s, 'SER');
    expect(ligand).not.toBe(serine);

    const { label_atom_id } = model.atomicHierarchy.atoms;
    const c12 = ligand.elements.find((e) => label_atom_id.value(e) === 'C12')!;
    const interUnit = getInterUnitNeighbours(model, () => Structure.ofModel(model));
    expect((interUnit.get(c12) ?? []).map((e) => label_atom_id.value(e))).toEqual(['OG']);

    const hash = chemCompHash(
      ligand,
      ligand.bonds,
      interUnit,
      0 as UnitIndex,
      ligand.elements.length as UnitIndex,
      'auto|',
    );
    expect(hash).toContain('C12-SER.OG');

    const rows = bondRows(ligand);
    const c12o1 = rows.find((r) => r.a === 'C12' && r.b === 'O1')!;
    expect(c12o1.order).toBe(2);
  });
});

// Basic cases (_spec §4).

const Benzene = [
  'HETATM    1  C1  BNZ A   1       1.390   0.000   0.000  1.00  0.00           C  ',
  'HETATM    2  C2  BNZ A   1       0.695   1.204   0.000  1.00  0.00           C  ',
  'HETATM    3  C3  BNZ A   1      -0.695   1.204   0.000  1.00  0.00           C  ',
  'HETATM    4  C4  BNZ A   1      -1.390   0.000   0.000  1.00  0.00           C  ',
  'HETATM    5  C5  BNZ A   1      -0.695  -1.204   0.000  1.00  0.00           C  ',
  'HETATM    6  C6  BNZ A   1       0.695  -1.204   0.000  1.00  0.00           C  ',
];

async function intraBondOrders(pdb: string) {
  const s = perceived(await modelFromPdb(pdb));
  const rows = bondRows(s.units[0] as Unit.Atomic);
  return {
    rows,
    orders: rows.map((r) => r.order),
    edgeCount: rows.length,
    aromaticCount: rows.filter((r) => r.aromatic).length,
  };
}

describe('bond-order perception: basic cases', () => {
  it('perceives an aromatic ring (benzene) as a Kekule structure', async () => {
    // planar regular hexagon, C-C = 1.39 A, residue BNZ (not a known residue,
    // no chem_comp_bond) -> orders must be perceived from coordinates
    const { orders, edgeCount, aromaticCount } = await intraBondOrders([...Benzene, 'END'].join('\n'));
    expect(edgeCount).toBe(6); // six ring bonds detected by distance
    // Kekule benzene: three double + three single
    expect(orders.filter((o) => o === 2).length).toBe(3);
    expect(orders.filter((o) => o === 1).length).toBe(3);
    // all ring bonds flagged aromatic
    expect(aromaticCount).toBe(6);
  });

  it('perceives orders for CONECT-derived connectivity (no explicit orders)', async () => {
    // benzene whose bonds come from CONECT records (basic connectivity only): struct_conn covalent bonds without
    // pdbx_value_order, flagged OrderUnknown, assigned a Kekule structure.
    const { orders, edgeCount, aromaticCount } = await intraBondOrders(
      [
        ...Benzene,
        'CONECT    1    2    6                                                            ',
        'CONECT    2    1    3                                                            ',
        'CONECT    3    2    4                                                            ',
        'CONECT    4    3    5                                                            ',
        'CONECT    5    4    6                                                            ',
        'CONECT    6    5    1                                                            ',
        'END',
      ].join('\n'),
    );
    expect(edgeCount).toBe(6);
    expect(orders.filter((o) => o === 2).length).toBe(3);
    expect(aromaticCount).toBe(6);
  });

  it('writes orders and the aromatic flag only: the input flags are kept (_spec §1.2)', async () => {
    const conect = [1, 2, 3, 4, 5, 6].map((i) => `CONECT${String(i).padStart(5)}${String((i % 6) + 1).padStart(5)}`);
    const model = await modelFromPdb([...Benzene, ...conect, 'END'].join('\n'));
    const builtIn = (Structure.ofModel(model).units[0] as Unit.Atomic).bonds.edgeProps.flags;
    const out = (perceived(model).units[0] as Unit.Atomic).bonds.edgeProps.flags;
    // CONECT bonds: OrderUnknown, never Computed (a user excluding `computed` bonds keeps them)
    expect(Array.from(out, (f) => f & ~BondType.Flag.Aromatic)).toEqual(Array.from(builtIn));
    expect(Array.from(out).some((f) => BondType.is(f, BondType.Flag.Computed))).toBe(false);
  });

  it('does not push a double onto an exocyclic amino N of an aromatic ring carbon', async () => {
    // 2-aminopyrimidine-like: ring carbon C2 is bonded to two ring nitrogens (N1, N3)
    // and one exocyclic amino nitrogen (N7). C2's pi bond belongs to the ring, so the
    // exocyclic C2-N7 bond must stay single (it must not be read as guanidinium).
    const { rows, orders } = await intraBondOrders(
      [
        'HETATM    1  N1  APM A   1       1.390   0.000   0.000  1.00  0.00           N  ',
        'HETATM    2  C2  APM A   1       0.695   1.204   0.000  1.00  0.00           C  ',
        'HETATM    3  N3  APM A   1      -0.695   1.204   0.000  1.00  0.00           N  ',
        'HETATM    4  C4  APM A   1      -1.390   0.000   0.000  1.00  0.00           C  ',
        'HETATM    5  C5  APM A   1      -0.695  -1.204   0.000  1.00  0.00           C  ',
        'HETATM    6  C6  APM A   1       0.695  -1.204   0.000  1.00  0.00           C  ',
        'HETATM    7  N7  APM A   1       1.370   2.373   0.000  1.00  0.00           N  ',
        'END',
      ].join('\n'),
    );
    expect(rows.find((r) => r.a === 'C2' && r.b === 'N7')!.order).toBe(1); // exocyclic amino bond stays single
    expect(orders.filter((o) => o === 2).length).toBe(3); // three ring doubles only
  });

  it('perceives a carboxylate (one C=O, one C-O)', async () => {
    // acetate-like: C(methyl)-C(=O)(-O), planar; residue ACX
    const { orders } = await intraBondOrders(
      [
        'HETATM    1  C   ACX A   1       0.000   0.000   0.000  1.00  0.00           C  ',
        'HETATM    2  CT  ACX A   1      -1.520   0.000   0.000  1.00  0.00           C  ',
        'HETATM    3  O1  ACX A   1       0.640   1.060   0.000  1.00  0.00           O  ',
        'HETATM    4  O2  ACX A   1       0.620  -1.080   0.000  1.00  0.00           O  ',
        'END',
      ].join('\n'),
    );
    // exactly one double bond (the C=O), the rest single
    expect(orders.filter((o) => o === 2).length).toBe(1);
  });

  it('leaves a canonical residue untouched (no spurious orders)', async () => {
    // glycine: only the backbone C=O is double (the core table's amino-acid C-O rule); perception adds nothing.
    const { orders } = await intraBondOrders(
      [
        'ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N  ',
        'ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00  0.00           C  ',
        'ATOM      3  C   GLY A   1       2.000   1.420   0.000  1.00  0.00           C  ',
        'ATOM      4  O   GLY A   1       1.250   2.390   0.000  1.00  0.00           O  ',
        'END',
      ].join('\n'),
    );
    // N-CA, CA-C single; C=O double (from the order table, not perception)
    expect(orders.filter((o) => o === 2).length).toBe(1);
  });

  it('handles altlocs', async () => {
    // PBD 6FUX has 1 altloc for a oxygen carbonyl (OG2). If altlocs are ignored,
    // the Carbon is perceived as sp3, connected to 3 single bonds.
    const { orders, edgeCount, aromaticCount } =
      await intraBondOrders(`HETATM 2139  C11 SRY A 303      33.664  71.619  31.003  1.00 51.05           C
HETATM 2140  N11 SRY A 303      33.655  72.167  32.286  1.00 56.10           N
HETATM 2141  CA1 SRY A 303      34.566  72.694  33.110  1.00 58.34           C
HETATM 2142  NB1 SRY A 303      34.052  73.088  34.310  1.00 63.09           N
HETATM 2143  NC1 SRY A 303      35.812  72.823  32.807  1.00 58.10           N
HETATM 2144  C21 SRY A 303      34.860  71.388  30.029  1.00 48.03           C
HETATM 2145  O21 SRY A 303      36.174  71.595  30.489  1.00 50.40           O
HETATM 2146  C31 SRY A 303      34.765  70.234  29.094  1.00 44.45           C
HETATM 2147  N31 SRY A 303      35.634  70.116  27.955  1.00 44.31           N
HETATM 2148  CD1 SRY A 303      36.819  69.542  27.774  1.00 46.65           C
HETATM 2149  NE1 SRY A 303      37.382  69.600  26.595  1.00 49.92           N
HETATM 2150  NF1 SRY A 303      37.471  68.930  28.635  1.00 49.37           N
HETATM 2151  C41 SRY A 303      33.402  69.700  28.905  1.00 42.20           C
HETATM 2152  O41 SRY A 303      33.523  68.452  28.242  1.00 36.53           O
HETATM 2153  C51 SRY A 303      32.319  69.698  30.034  1.00 45.87           C
HETATM 2154  O51 SRY A 303      31.699  68.366  30.386  1.00 45.09           O
HETATM 2155  C61 SRY A 303      32.490  70.642  31.207  1.00 52.73           C
HETATM 2156  O61 SRY A 303      32.360  70.076  32.562  1.00 58.33           O
HETATM 2157  C12 SRY A 303      33.161  68.470  26.811  1.00 33.80           C
HETATM 2158  C22 SRY A 303      33.048  67.123  26.293  1.00 32.15           C
HETATM 2159  C32 SRY A 303      31.547  66.774  26.426  1.00 31.57           C
HETATM 2160  O32 SRY A 303      30.941  66.628  25.099  1.00 30.00           O
HETATM 2161  CG2 SRY A 303      31.377  65.493  27.278  1.00 32.12           C
HETATM 2162  OG2ASRY A 303      31.727  64.432  26.692  0.50 30.76           O
HETATM 2163  OG2BSRY A 303      30.170  65.084  27.156  0.50 32.03           O
HETATM 2164  C42 SRY A 303      31.014  68.073  27.034  1.00 33.38           C
HETATM 2165  CH2 SRY A 303      29.651  68.451  26.425  1.00 34.28           C
HETATM 2166  O42 SRY A 303      31.910  69.142  26.745  1.00 34.07           O
HETATM 2167  O13 SRY A 303      33.480  66.820  24.969  1.00 32.07           O
HETATM 2168  C13 SRY A 303      34.864  66.982  24.694  1.00 33.00           C
HETATM 2169  C23 SRY A 303      35.167  66.104  23.483  1.00 33.44           C
HETATM 2170  N23 SRY A 303      34.699  64.715  23.476  1.00 34.30           N
HETATM 2171  CI3 SRY A 303      35.053  63.964  24.721  1.00 36.28           C
HETATM 2172  C33 SRY A 303      34.871  66.901  22.206  1.00 32.82           C
HETATM 2173  O33 SRY A 303      34.998  66.199  21.104  1.00 32.38           O
HETATM 2174  C43 SRY A 303      35.183  68.418  22.241  1.00 34.43           C
HETATM 2175  O43 SRY A 303      34.580  69.116  21.102  1.00 34.09           O
HETATM 2176  C53 SRY A 303      35.037  69.027  23.492  1.00 35.62           C
HETATM 2177  O53 SRY A 303      35.172  68.323  24.700  1.00 35.46           O
HETATM 2178  C63 SRY A 303      35.283  70.514  23.598  1.00 37.81           C
HETATM 2179  O63 SRY A 303      35.420  70.895  24.913  1.00 42.50           O
END`);
    expect(edgeCount).toBe(43);
    expect(orders.filter((o) => o === 2).length).toBe(4);
    expect(aromaticCount).toBe(0);
  });
});

describe('bond-order perception: residue gate (_spec §1.9)', () => {
  it('auto skips a residue whose orders are all stated, and never builds the inter-unit neighbours', async () => {
    const model = await modelFromPdb([...Benzene, 'END'].join('\n'));
    const unit = unitOf(Structure.ofModel(model), 'BNZ');
    // the distance bonds as if a file stated them single (no Computed flag)
    const { offset, a, b, edgeCount, edgeProps } = unit.bonds;
    const flags = edgeProps.flags.map((f) => f & ~BondType.Flag.Computed);
    const stated = IntAdjacencyGraph.create(
      offset,
      a,
      b,
      edgeCount,
      { ...edgeProps, flags },
      unit.bonds.props,
    ) as IntraUnitBonds;
    const getInterUnit = vi.fn(() => getInterUnitNeighbours(model, () => Structure.ofModel(model)));

    const out = perceiveBondOrders(unit, stated, getInterUnit, 'auto');
    expect(Array.from(out.edgeProps.order)).toEqual(Array.from(stated.edgeProps.order));
    expect(Array.from(out.edgeProps.flags)).toEqual(Array.from(flags));
    expect(out.props.canRemap).toBe(stated.props.canRemap);
    expect(getInterUnit).not.toHaveBeenCalled();

    // the same ring with unknown orders is perceived
    const perceivedBonds = perceiveBondOrders(unit, unit.bonds, getInterUnit, 'auto');
    expect(Array.from(perceivedBonds.edgeProps.order).filter((o) => o === 2).length).toBe(6); // 3 bonds, both directions
    expect(getInterUnit).toHaveBeenCalled();
  });
});
