/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import * as fs from 'fs';
import * as path from 'path';
import { Tokenizer } from '../../../../../../mol-io/reader/common/text/tokenizer';
import { PdbFile } from '../../../../../../mol-io/reader/pdb/schema';
import { pdbToMmCif } from '../../../../../../mol-model-formats/structure/pdb/to-cif';
import { trajectoryFromMmCIF } from '../../../../../../mol-model-formats/structure/mmcif';
import { Task } from '../../../../../../mol-task';
import { Structure } from '../../../structure';
import { Unit } from '../../../unit';
import { BondType } from '../../../../model/types';

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

function intraBondOrders(structure: Structure) {
    const unit = structure.units[0] as Unit.Atomic;
    const { edgeCount, offset, b, edgeProps } = unit.bonds;
    const orders: number[] = [];
    const aromatic: number[] = [];
    for (let a = 0; a < unit.elements.length; a++) {
        for (let t = offset[a]; t < offset[a + 1]; t++) {
            if (a < b[t]) {
                orders.push(edgeProps.order[t]);
                if (edgeProps.flags[t] & BondType.Flag.Aromatic) aromatic.push(1);
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

    it('does not add a second C=O when the table already placed one (C-terminal carboxylate)', async () => {
        // ALA is an amino acid (so the table gives backbone C=O order 2) but is NOT in
        // IntraBondOrderTable, so it is no longer skipped wholesale - perception runs on
        // it. Its C-terminal carboxylate (C bonded to O + OXT) matches the carboxyl group,
        // but since C already has its table C=O, no second double may be assigned.
        const structure = await structureFromPdb([
            'ATOM      1  N   ALA A   1      -0.500   1.400   0.000  1.00  0.00           N  ',
            'ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  ',
            'ATOM      3  CB  ALA A   1      -0.500  -1.400   0.000  1.00  0.00           C  ',
            'ATOM      4  C   ALA A   1       1.520   0.000   0.000  1.00  0.00           C  ',
            'ATOM      5  O   ALA A   1       2.140   1.050   0.000  1.00  0.00           O  ',
            'ATOM      6  OXT ALA A   1       2.140  -1.050   0.000  1.00  0.00           O  ',
            'END                                                                             ',
        ].join('\n'));

        const { orders } = intraBondOrders(structure);
        // exactly one double bond (the backbone/table C=O); C-OXT must stay single
        expect(orders.filter(o => o === 2).length).toBe(1);
    });
});

// --- data-driven validation suite -------------------------------------------

interface ExpectedBond {
    a: string;
    /** Single atom name, or an array when any of the listed partners is equally valid
     *  (e.g. interchangeable nitro O, carboxylate O detected from CCD by build-data.mjs). */
    b: string | string[];
    order?: 2 | 3;
    aromatic?: true;
}
/** A bond (atom-name pair, unordered) whose perceived order the suite knowingly gets wrong; `reason`
 *  documents the open issue so it stays tracked. Exempts the pair from BOTH the forward check (an
 *  expected multiple that is missed) and the reverse check (a spurious multiple, e.g. when a missed
 *  double leaves an atom's valence to be matched elsewhere). The structure still runs and every
 *  other bond is checked. */
interface SkipBond { a: string; b: string; reason: string; }
interface ManifestEntry { pdbId: string; compId: string; skipBonds?: SkipBond[]; }

const dataDir = path.join(__dirname, 'data');

async function structureFromDataFile(filename: string) {
    const text = fs.readFileSync(path.join(dataDir, filename), 'utf8');
    return structureFromPdb(text);
}

function checkPerceived(structure: Structure, expected: ExpectedBond[], skipBonds: SkipBond[] = []) {
    const unit = structure.units[0] as Unit.Atomic;
    const { label_atom_id, label_comp_id } = unit.model.atomicHierarchy.atoms;
    const nameToLocal = new Map<string, number>();
    const resname = label_comp_id.value(unit.elements[0]);
    for (let i = 0; i < unit.elements.length; i++) {
        nameToLocal.set(label_atom_id.value(unit.elements[i]), i);
    }
    const { offset, b, edgeProps } = unit.bonds;

    // Set of every atom-pair the CCD considers non-single (double / triple / aromatic),
    // as a sorted "min|max" local-index key. Interchangeable partners contribute all
    // their variants. Used by the reverse check to reject spurious perceived multiples.
    const expectedKey = (i: number, j: number) => `${Math.min(i, j)}|${Math.max(i, j)}`;
    const expectedPairs = new Set<string>();

    // Known-unperceivable bonds whose forward check is exempted (see SkipBond / manifest.json).
    const localOf = (name: string) => {
        const i = nameToLocal.get(name);
        if (i === undefined) throw new Error(`skipBonds atom not found: ${name}`);
        return i;
    };
    const skipKeys = new Set(skipBonds.map(s => expectedKey(localOf(s.a), localOf(s.b))));

    for (const exp of expected) {
        const bNames = Array.isArray(exp.b) ? exp.b : [exp.b];
        const bondLabel = `${resname} ${exp.a}-(${bNames.join('|')})`;
        const u = nameToLocal.get(exp.a);
        if (u === undefined) throw new Error(`atom not found: ${exp.a}`);
        const vs = bNames.map(name => {
            const v = nameToLocal.get(name);
            if (v === undefined) throw new Error(`atom not found: ${name}`);
            return v;
        });
        for (const v of vs) expectedPairs.add(expectedKey(u, v));
        // For interchangeable partners, succeed as soon as any partner has the expected bond.
        let found = false;
        for (const v of vs) {
            for (let t = offset[u]; t < offset[u + 1]; t++) {
                if (b[t] !== v) continue;
                if (exp.order !== undefined && edgeProps.order[t] !== exp.order) continue;
                if (exp.aromatic && !(edgeProps.flags[t] & BondType.Flag.Aromatic)) continue;
                found = true;
            }
            if (found) break;
        }
        // A known-unperceivable bond is exempt from the forward check; flag it if it is now
        // perceived correctly so the obsolete skip can be removed from the manifest.
        if (vs.some(v => skipKeys.has(expectedKey(u, v)))) {
            if (found) {
                console.warn(`obsolete skipBonds entry: ${bondLabel} is now perceived correctly`);
            }
            continue;
        }
        expect({ bond: bondLabel, found }).toEqual({ bond: bondLabel, found: true });
    }

    // Reverse check: perception must not invent multiple bonds the CCD doesn't have.
    // (A perceived order>1 inside an aromatic ring is allowed against the expected aromatic
    // pair — our Kekule may differ from the CCD's, but it stays within the same ring bonds.)
    const nameOf = (i: number) => label_atom_id.value(unit.elements[i]);
    for (let u = 0; u < unit.elements.length; u++) {
        for (let t = offset[u]; t < offset[u + 1]; t++) {
            const v = b[t];
            if (u >= v) continue;
            if (edgeProps.order[t] <= 1) continue;
            if (skipKeys.has(expectedKey(u, v))) continue; // known-wrong pair, exempt both ways
            const spurious = !expectedPairs.has(expectedKey(u, v));
            const bondLabel = `${resname} ${nameOf(u)}=${nameOf(v)} (order ${edgeProps.order[t]})`;
            expect({ bond: bondLabel, spurious }).toEqual({ bond: bondLabel, spurious: false });
        }
    }
}

const manifestPath = path.join(dataDir, 'manifest.json');
const manifest: ManifestEntry[] = fs.existsSync(manifestPath)
    ? JSON.parse(fs.readFileSync(manifestPath, 'utf8'))
    : [];

describe('bond-order perception (data)', () => {
    for (const { pdbId, compId, skipBonds } of manifest) {
        const stem = `${pdbId}_${compId}`;
        const pdbPath = path.join(dataDir, `${stem}.pdb`);
        if (!fs.existsSync(pdbPath)) {
            it.skip(`${pdbId} / ${compId} — run build-data.mjs to generate data`, () => {});
            continue;
        }
        it(`${pdbId} / ${compId}`, async () => {
            const structure = await structureFromDataFile(`${stem}.pdb`);
            const expected: ExpectedBond[] = JSON.parse(
                fs.readFileSync(path.join(dataDir, `${stem}_expected.json`), 'utf8')
            );
            checkPerceived(structure, expected, skipBonds ?? []);
        });
    }
});
