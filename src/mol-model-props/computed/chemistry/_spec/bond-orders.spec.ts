/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Tokenizer } from '../../../../mol-io/reader/common/text/tokenizer';
import { PdbFile } from '../../../../mol-io/reader/pdb/schema';
import { pdbToMmCif } from '../../../../mol-model-formats/structure/pdb/to-cif';
import { trajectoryFromMmCIF } from '../../../../mol-model-formats/structure/mmcif';
import { Task } from '../../../../mol-task';
import { Structure } from '../../../../mol-model/structure';
import { Unit } from '../../../../mol-model/structure';
import { BondType } from '../../../../mol-model/structure/model/types';
import { computeUnitBondOrders } from '../bond-orders';

function makePdb(pdbText: string): PdbFile {
    const lines = Tokenizer.readAllLines(pdbText);
    return { lines, variant: 'pdb' };
}

async function structureFromPdb(pdbText: string) {
    const cif = await pdbToMmCif(makePdb(pdbText));
    const trajectory = await trajectoryFromMmCIF(cif).run();
    const model = await Task.resolveInContext(trajectory.getFrameAtIndex(0));
    return Structure.ofModel(model);
}

/**
 * Resolve the perceived per-edge order/flags for the structure's first unit. Perception is now an
 * opt-in computed property, so the bond graph no longer carries perceived orders — we run
 * `computeUnitBondOrders` ('auto') and read its per-edge override arrays (parallel to the unit's bond edges).
 */
function perceivedEdges(structure: Structure) {
    const unit = structure.units[0] as Unit.Atomic;
    const { edgeCount, offset, b } = unit.bonds;
    const ov = computeUnitBondOrders(structure, unit, 'auto');
    const order = ov ? ov.order : unit.bonds.edgeProps.order;
    const flags = ov ? ov.flags : unit.bonds.edgeProps.flags;
    return { unit, edgeCount, offset, b, order, flags };
}

function intraBondOrders(structure: Structure) {
    const { unit, edgeCount, offset, b, order, flags } = perceivedEdges(structure);
    const orders: number[] = [];
    const aromatic: number[] = [];
    for (let a = 0; a < unit.elements.length; a++) {
        for (let t = offset[a]; t < offset[a + 1]; t++) {
            if (a < b[t]) {
                orders.push(order[t]);
                if (flags[t] & (BondType.Flag.Aromatic | BondType.Flag.AromaticHuckel)) aromatic.push(1);
            }
        }
    }
    return { orders, edgeCount, aromaticCount: aromatic.length };
}

describe('bond-order perception (Sayle)', () => {
    it('perceives an aromatic ring (benzene) as a Kekule structure', async () => {
        // planar regular hexagon, C-C = 1.39 A, residue BNZ (not a known residue,
        // no chem_comp_bond) -> orders must be perceived from coordinates
        const structure = await structureFromPdb([
            'HETATM    1  C1  BNZ A   1       1.390   0.000   0.000  1.00  0.00           C  ',
            'HETATM    2  C2  BNZ A   1       0.695   1.204   0.000  1.00  0.00           C  ',
            'HETATM    3  C3  BNZ A   1      -0.695   1.204   0.000  1.00  0.00           C  ',
            'HETATM    4  C4  BNZ A   1      -1.390   0.000   0.000  1.00  0.00           C  ',
            'HETATM    5  C5  BNZ A   1      -0.695  -1.204   0.000  1.00  0.00           C  ',
            'HETATM    6  C6  BNZ A   1       0.695  -1.204   0.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const { orders, edgeCount, aromaticCount } = intraBondOrders(structure);
        expect(edgeCount).toBe(6); // six ring bonds detected by distance
        // Kekule benzene: three double + three single
        expect(orders.filter(o => o === 2).length).toBe(3);
        expect(orders.filter(o => o === 1).length).toBe(3);
        // all ring bonds flagged aromatic
        expect(aromaticCount).toBe(6);
    });

    it('perceives orders for CONECT-derived connectivity (no explicit orders)', async () => {
        // benzene whose bonds come from CONECT records (basic connectivity only).
        // These become struct_conn covalent bonds without pdbx_value_order, which must
        // be marked perceivable and assigned a Kekule structure.
        const structure = await structureFromPdb([
            'HETATM    1  C1  BNZ A   1       1.390   0.000   0.000  1.00  0.00           C  ',
            'HETATM    2  C2  BNZ A   1       0.695   1.204   0.000  1.00  0.00           C  ',
            'HETATM    3  C3  BNZ A   1      -0.695   1.204   0.000  1.00  0.00           C  ',
            'HETATM    4  C4  BNZ A   1      -1.390   0.000   0.000  1.00  0.00           C  ',
            'HETATM    5  C5  BNZ A   1      -0.695  -1.204   0.000  1.00  0.00           C  ',
            'HETATM    6  C6  BNZ A   1       0.695  -1.204   0.000  1.00  0.00           C  ',
            'CONECT    1    2    6                                                            ',
            'CONECT    2    1    3                                                            ',
            'CONECT    3    2    4                                                            ',
            'CONECT    4    3    5                                                            ',
            'CONECT    5    4    6                                                            ',
            'CONECT    6    5    1                                                            ',
            'END                                                                             ',
        ].join('\n'));

        const { orders, edgeCount, aromaticCount } = intraBondOrders(structure);
        expect(edgeCount).toBe(6);
        expect(orders.filter(o => o === 2).length).toBe(3);
        expect(aromaticCount).toBe(6);
    });

    it('does not push a double onto an exocyclic amino N of an aromatic ring carbon', async () => {
        // 2-aminopyrimidine-like: ring carbon C2 is bonded to two ring nitrogens (N1, N3)
        // and one exocyclic amino nitrogen (N7). C2's pi bond belongs to the ring, so the
        // exocyclic C2-N7 bond must stay single (it must not be read as guanidinium).
        const structure = await structureFromPdb([
            'HETATM    1  N1  APM A   1       1.390   0.000   0.000  1.00  0.00           N  ',
            'HETATM    2  C2  APM A   1       0.695   1.204   0.000  1.00  0.00           C  ',
            'HETATM    3  N3  APM A   1      -0.695   1.204   0.000  1.00  0.00           N  ',
            'HETATM    4  C4  APM A   1      -1.390   0.000   0.000  1.00  0.00           C  ',
            'HETATM    5  C5  APM A   1      -0.695  -1.204   0.000  1.00  0.00           C  ',
            'HETATM    6  C6  APM A   1       0.695  -1.204   0.000  1.00  0.00           C  ',
            'HETATM    7  N7  APM A   1       1.370   2.373   0.000  1.00  0.00           N  ',
            'END                                                                             ',
        ].join('\n'));

        const unit = structure.units[0] as Unit.Atomic;
        const { label_atom_id } = unit.model.atomicHierarchy.atoms;
        const local = new Map<string, number>();
        for (let i = 0; i < unit.elements.length; i++) local.set(label_atom_id.value(unit.elements[i]), i);
        const c2 = local.get('C2')!, n7 = local.get('N7')!;
        const { offset, b, edgeProps } = unit.bonds;
        let exocyclicOrder = -1;
        for (let t = offset[c2]; t < offset[c2 + 1]; t++) if (b[t] === n7) exocyclicOrder = edgeProps.order[t];
        expect(exocyclicOrder).toBe(1); // exocyclic amino bond stays single

        const { orders } = intraBondOrders(structure);
        expect(orders.filter(o => o === 2).length).toBe(3); // three ring doubles only
    });

    it('perceives a carboxylate (one C=O, one C-O)', async () => {
        // acetate-like: C(methyl)-C(=O)(-O), planar; residue ACX
        const structure = await structureFromPdb([
            'HETATM    1  C   ACX A   1       0.000   0.000   0.000  1.00  0.00           C  ',
            'HETATM    2  CT  ACX A   1      -1.520   0.000   0.000  1.00  0.00           C  ',
            'HETATM    3  O1  ACX A   1       0.640   1.060   0.000  1.00  0.00           O  ',
            'HETATM    4  O2  ACX A   1       0.620  -1.080   0.000  1.00  0.00           O  ',
            'END                                                                             ',
        ].join('\n'));

        const { orders } = intraBondOrders(structure);
        // exactly one double bond (the C=O), the rest single
        expect(orders.filter(o => o === 2).length).toBe(1);
    });

    it('leaves a residue with a table template untouched (no spurious orders)', async () => {
        // glycine: only the backbone C=O is double (from the table's AminoAcidNames C-O
        // special case); perception must add nothing.
        const structure = await structureFromPdb([
            'ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N  ',
            'ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00  0.00           C  ',
            'ATOM      3  C   GLY A   1       2.000   1.420   0.000  1.00  0.00           C  ',
            'ATOM      4  O   GLY A   1       1.250   2.390   0.000  1.00  0.00           O  ',
            'END                                                                             ',
        ].join('\n'));

        const { orders } = intraBondOrders(structure);
        // N-CA, CA-C single; C=O double (from the order table, not perception)
        expect(orders.filter(o => o === 2).length).toBe(1);
    });

    it('handles altlocs', async () => {
        // PBD 6FUX has 1 altloc for a oxygen carbonyl (OG2). If altlocs are ignored,
        // the Carbon is perceived as sp3, connected to 3 single bonds.
        const structure = await structureFromPdb(`HETATM 2139  C11 SRY A 303      33.664  71.619  31.003  1.00 51.05           C  
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

        const { orders, edgeCount, aromaticCount } = intraBondOrders(structure);
        expect(edgeCount).toBe(43);
        expect(orders.filter(o => o === 2).length).toBe(4);
        expect(aromaticCount).toBe(0);
    });
});
