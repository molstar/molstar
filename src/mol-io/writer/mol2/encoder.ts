/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Category } from '../cif/encoder';
import { LigandEncoder } from '../ligand-encoder';
import { StringBuilder } from '../../../mol-util';
import { getCategoryInstanceData } from '../cif/encoder/util';
import { BondType } from '../../../mol-model/structure/model/types';
// import { ComponentBond } from '../../../mol-model-formats/structure/property/bonds/comp';

// // type MOL_TYPE = 'SMALL' | 'BIOPOLYMER' | 'PROTEIN' | 'NUCLEIC_ACID' | 'SACCHARIDE';
// // type CHARGE_TYPE = 'NO_CHARGES' | 'DEL_RE' | 'GASTEIGER' | 'GAST_HUCK' | 'HUCKEL' | 'PULLMAN' | 'GAUSS80_CHARGES' | 'AMPAC_CHARGES' | 'MULLIKEN_CHARGES' | 'DICT_ CHARGES' | 'MMFF94_CHARGES' | 'USER_CHARGES';
// const NON_METAL_ATOMS = 'H D B C N O F Si P S Cl As Se Br Te I At He Ne Ar Kr Xe Rn'.split(' ');
// type BondData = { label_atom_id: string, order: number, aromatic: boolean };

// specification: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
export class Mol2Encoder extends LigandEncoder {
    private out: StringBuilder;

    _writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void {
        const a = StringBuilder.create();
        const b = StringBuilder.create();
        const { instance, source } = getCategoryInstanceData(category, context);

        // write header
        const name = this.getName(instance, source);
        StringBuilder.writeSafe(this.builder, `# Name: ${name}\n# Created by ${this.encoder}\n\n`);

        const bondMap = this.componentBondData.entries.get(name)!;
        let bondCount = 0;

        const atoms = this.getAtoms(instance, source);
        StringBuilder.writeSafe(a, '@<TRIPOS>ATOM\n');
        StringBuilder.writeSafe(b, '@<TRIPOS>BOND\n');
        atoms.forEach((atom1, label_atom_id1) => {
            const { index: i1 } = atom1;
            bondMap.map.get(label_atom_id1)!.forEach((bond, label_atom_id2) => {
                const atom2 = atoms.get(label_atom_id2);
                if (!atom2) return;

                const { index: i2, type_symbol: type_symbol2 } = atom2;
                if (i1 < i2 && !this.skipHydrogen(type_symbol2)) {
                    const { order, flags } = bond;
                    const ar = BondType.is(BondType.Flag.Aromatic, flags);
                    StringBuilder.writeSafe(b, `${++bondCount} ${i1 + 1} ${i2 + 1} ${ar ? 'ar' : order}`);
                    StringBuilder.newline(b);
                }
            });

            // const sybyl = this.mapToSybyl(label_atom_id1, atom1.type_symbol, bondMap);
            const sybyl = '?';
            StringBuilder.writeSafe(a, `${i1 + 1} ${label_atom_id1} ${atom1.Cartn_x.toFixed(3)} ${atom1.Cartn_y.toFixed(3)} ${atom1.Cartn_z.toFixed(3)} ${sybyl} 1 ${name} 0.000\n`);
        });

        // could write something like 'SMALL\nNO_CHARGES', for now let's write **** indicating non-optional, yet missing, string values
        StringBuilder.writeSafe(this.out, `@<TRIPOS>MOLECULE\n${name}\n${atoms.size} ${bondCount} 1\n****\n****\n\n`);
        StringBuilder.writeSafe(this.out, StringBuilder.getString(a));
        StringBuilder.writeSafe(this.out, StringBuilder.getString(b));
        StringBuilder.writeSafe(this.out, `@<TRIPOS>SUBSTRUCTURE\n1 ${name} 1\n`);
    }

    // private toArray(map: Map<string, { order: number, flags: number }>): BondData[] {
    //     const array: BondData[] = [];
    //     map.forEach((v, k) => array.push({ label_atom_id: k, order: v.order, aromatic: BondType.is(BondType.Flag.Aromatic, v.flags) }));
    //     return array;
    // }

    // see https://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
    // cannot account for covalently bound amino acids etc
    // private mapToSybyl(label_atom_id: string, type_symbol: string, bondMap: ComponentBond.Entry) {
    //     const partialBondMap = bondMap.map.get(label_atom_id)!;
    //     const bond = this.toArray(partialBondMap);

    //     const num_bond = bond.length;
    //     const nonmet = bond.filter(b => this.isNonMetalBond(b));
    //     const num_nonmet = nonmet.length;
    //     const arom = bond.filter(e => e.aromatic);
    //     const num_arom = arom.length;

    //     // TODO if altLoc: 'Du' // 1.1
    //     // TODO if end of polymeric bond: 'Du' // 1.2
    //     if (type_symbol === 'D') return 'H'; // 1.3
    //     if (type_symbol === 'P') return 'P.3'; // 1.4, 4mpo/ligand?encoding=mol2&auth_seq_id=203 (PO4)
    //     if (type_symbol === 'Co' || type_symbol === 'Ru') return type_symbol + '.oh'; // 1.5
    //     if (type_symbol === 'C') { // 1.6
    //         if (num_bond >= 4 && bond.every(b => b.order === 1)) return 'C.3'; // 1.6.1, 3rga/ligand?encoding=mol2&auth_seq_id=307 (MOH)
    //         if (num_bond === 3 && this.isCat(label_atom_id, bond, bondMap)) return 'C.cat'; // 1.6.2, 1acj/ligand?encoding=mol2&auth_seq_id=44 (ARG), 5vjb/ligand?encoding=mol2&auth_seq_id=101 (GAI)
    //         if (num_bond >= 2 && num_arom >= 2) return 'C.ar'; // 1.6.3, 1acj/ligand?encoding=mol2&auth_seq_id=30 (PHE), 1acj/ligand?encoding=mol2&auth_seq_id=63 (TYR), 1acj/ligand?encoding=mol2&auth_seq_id=84 (TRP), 1acj/ligand?encoding=mol2&auth_seq_id=999 (THA)
    //         if ((num_bond === 1 || num_bond === 2) && bond.filter(b => b.order === 3).length === 1) return 'C.1'; // 1.6.4, 3i04/ligand?encoding=mol2&auth_asym_id=C&auth_seq_id=900 (CYN)
    //         return 'C.2'; // 1.6.5
    //     }
    //     if (type_symbol === 'O') { // 1.7
    //         if (num_nonmet === 1) { // 1.7.1
    //             if (this.isOC(nonmet[0], bondMap)) return 'O.co2'; // 1.7.1.1, 4h2v/ligand?encoding=mol2&auth_seq_id=403 (ACT)
    //             if (this.isOP(nonmet[0], bondMap)) return 'O.co2'; // 1.7.1.2, 4mpo/ligand?encoding=mol2&auth_seq_id=203 (PO4)
    //         }
    //         if (num_nonmet >= 2 && bond.every(b => b.order === 1)) return 'O.3'; // 1.7.2, 1acj/ligand?encoding=mol2&auth_seq_id=601 (HOH), 3rga/ligand?encoding=mol2&auth_seq_id=307 (MOH)
    //         return 'O.2'; // 1.7.3, 1acj/ligand?encoding=mol2&auth_seq_id=4 (SER)
    //     }
    //     if (type_symbol === 'N') { // 1.8
    //         if (num_nonmet === 4 && bond.every(b => b.order === 1)) return 'N.4'; // 1.8.1, 4ikf/ligand?encoding=mol2&auth_seq_id=403 (NH4)
    //         if (num_bond >= 2 && num_arom === 2) return 'N.ar'; // 1.8.2, 1acj/ligand?encoding=mol2&auth_seq_id=84 (TRP), 1acj/ligand?encoding=mol2&auth_seq_id=999 (THA)
    //         if (num_nonmet === 1 && nonmet.some(b => b.order === 3)) return 'N.1'; // 1.8.3, 3i04/ligand?encoding=mol2&auth_asym_id=C&auth_seq_id=900 (CYN)
    //         if (num_nonmet === 2 && (nonmet[0].order + nonmet[1].order === 4)) return 'N.1'; // 1.8.4, 3sbr/ligand?encoding=mol2&auth_seq_id=640&auth_asym_id=D (N2O)
    //         if (num_nonmet === 3 && this.hasCOCS(nonmet, bondMap)) return 'N.am'; // 1.8.5, 3zfz/ligand?encoding=mol2&auth_seq_id=1669 (1W8)
    //         if (num_nonmet === 3) { // 1.8.6
    //             if (nonmet.filter(b => b.order > 1).length === 1) return 'N.pl3'; // 1.8.6.1, 4hon/ligand?encoding=mol2&auth_seq_id=407 (NO3)
    //             if (nonmet.every(b => b.order === 1)) {
    //                 if (this.isNpl3(nonmet, bondMap)) return 'N.pl3'; // 1.8.6.1.1 & 1.8.6.1.2, 1acj/ligand?encoding=mol2&auth_seq_id=44 (ARG), 5vjb/ligand?encoding=mol2&auth_seq_id=101 (GAI)
    //             }
    //             return 'N.3';
    //         }
    //         return 'N.2'; // 1.8.7, 1acj/ligand?encoding=mol2&auth_seq_id=4 (SER)
    //     }
    //     if (type_symbol === 'S') { // 1.9
    //         if (num_nonmet === 3 && this.countOfOxygenWithSingleNonmet(nonmet, bondMap) === 1) return 'S.o'; // 1.9.1, 4i03/ligand?encoding=mol2&auth_seq_id=312 (DMS)
    //         if (num_nonmet === 4 && this.countOfOxygenWithSingleNonmet(nonmet, bondMap) === 2) return 'S.o2'; // 1.9.2, 1udt/ligand?encoding=mol2&auth_seq_id=1000 (VIA)
    //         if (num_nonmet >= 2 && bond.every(b => b.order === 1)) return 'S.3'; // 1.9.3, 3zfz/ligand?encoding=mol2&auth_seq_id=1669 (1W8)
    //         return 'S.2'; // 1.9.4, 4gpc/ligand?encoding=mol2&auth_seq_id=902 (SO4)
    //     }
    //     if (type_symbol === 'Ti' || type_symbol === 'Cr') { // 1.10
    //         return type_symbol + (num_bond <= 4 ? '.th' : '.oh'); // 1.10.1 & 1.10.2
    //     }
    //     return type_symbol; // 1.11
    // }

    // private isNonMetalBond(b: BondData): boolean {
    //     // would be nice to have type_symbol here...
    //     return NON_METAL_ATOMS.some(a => this.getLabel(b.label_atom_id).startsWith(a));
    // }

    // // 1.8.6.2.1: If one single bond is to an atom that forms a bond of type double, triple, aromatic or
    // // delocalised .AND. one other single bond is to H then atom_type is N.pl3
    // // 1.8.6.2.2: If one single bond is to an atom that forms a bond of type double, triple, aromatic or
    // // delocalised .AND. neither of the other single bonds are to H .AND. sum_of_angles around N .ge. 350 deg then atom_type is N.pl3
    // // TODO cannot check accurately for delocalized bonds
    // // TODO cannot check accurately for 2nd criterion without coordinates
    // private isNpl3(nonmet: BondData[], bondMap: ComponentBond.Entry): boolean {
    //     for (let i = 0, il = nonmet.length; i < il; i++) {
    //         const consumed = nonmet[i];
    //         // determine index that fulfills 1st criterion
    //         if (this.toArray(bondMap.map.get(consumed.label_atom_id)!).some(b => b.order > 1 || b.aromatic)) {
    //             if (nonmet.filter(b => b !== consumed).filter(b => this.getLabel(b.label_atom_id).startsWith('H')).length === 1) return true; // 1.8.6.2.1
    //             if (nonmet.filter(b => b !== consumed).every(b => !this.getLabel(b.label_atom_id).startsWith('H'))) return true; // 1.8.6.2.2
    //         }
    //     }
    //     return false;
    // }

    // // If bond is to carbon .AND. carbon forms a total of 3 bonds, 2 of which are to an oxygen
    // // forming only 1 non-metal bond then atom_type is O.co2
    // private isOC(nonmet: BondData, bondMap: ComponentBond.Entry): boolean {
    //     if (!this.getLabel(nonmet.label_atom_id).startsWith('C')) return false;
    //     const carbonBonds = this.toArray(bondMap.map.get(nonmet.label_atom_id)!);
    //     if (carbonBonds.length !== 3) return false;
    //     return carbonBonds.filter(b => this.getLabel(b.label_atom_id).startsWith('O') &&
    //         this.toArray(bondMap.map.get(b.label_atom_id)!).filter(ob => this.isNonMetalBond(ob)).length === 1).length === 2;
    // }

    // // If bond is to phosphorus .AND. phosphorus forms at least 2 bonds to an oxygen forming
    // // only 1 non-metal bond then atom_type is O.co2
    // private isOP(nonmet: BondData, bondMap: ComponentBond.Entry): boolean {
    //     if (!this.getLabel(nonmet.label_atom_id).startsWith('P')) return false;
    //     const phosphorusBonds = this.toArray(bondMap.map.get(nonmet.label_atom_id)!);
    //     if (phosphorusBonds.length < 2) return false;
    //     return phosphorusBonds.filter(b => this.getLabel(b.label_atom_id).startsWith('O') &&
    //         this.toArray(bondMap.map.get(b.label_atom_id)!).filter(ob => this.isNonMetalBond(ob)).length === 1).length >= 2;
    // }

    // // If num_bond .eq. 3 .AND. all bonds are acyclic .AND. all bonds are to nitrogen .AND. each
    // // nitrogen forms bonds to 2 other atoms both of which are not oxygen then atom_type is C.cat.
    // private isCat(root: string, bond: BondData[], bondMap: ComponentBond.Entry): boolean {
    //     if (bond.some(b => !this.getLabel(b.label_atom_id).startsWith('N'))) return false;
    //     const nitrogenBonds = bond.map(b => b.label_atom_id).map(label_atom_id => this.toArray(bondMap.map.get(label_atom_id)!));

    //     // ensure no cycles
    //     const all = [];
    //     const unique = new Set();
    //     nitrogenBonds.forEach(a => a.map(b => b.label_atom_id).filter(lai => lai !== root).forEach(lai => { all.push(lai); unique.add(lai); }));
    //     if (all.length !== unique.size) return false;

    //     return nitrogenBonds.every(a => a.length >= 2 && a.every(b => b.label_atom_id !== 'O'));
    // }

    // private countOfOxygenWithSingleNonmet(nonmet: BondData[], bondMap: ComponentBond.Entry): number {
    //     return nonmet.map(b => b.label_atom_id)
    //         .filter(label_atom_id => this.getLabel(label_atom_id).startsWith('O'))
    //         .map(label_atom_id => this.toArray(bondMap.map.get(label_atom_id)!)
    //             .filter(b => this.isNonMetalBond(b)).length === 1)
    //         .length;
    // }

    // // If num_nonmet .eq. 3 .AND. one bond is to C=O or C=S then atom_type is N.am
    // private hasCOCS(nonmet: BondData[], bondMap: ComponentBond.Entry): boolean {
    //     return nonmet.map(b => b.label_atom_id)
    //         .filter(label_atom_id => this.getLabel(label_atom_id).startsWith('C'))
    //         .filter(label_atom_id => this.toArray(bondMap.map.get(label_atom_id)!)
    //             .filter(b => b.order === 2)
    //             .filter(b => this.getLabel(b.label_atom_id).startsWith('O') || this.getLabel(b.label_atom_id).startsWith('S')))
    //         .length === 1;
    // }

    protected writeFullCategory<Ctx>(sb: StringBuilder, category: Category<Ctx>, context?: Ctx) {
        const { instance, source } = getCategoryInstanceData(category, context);
        const fields = instance.fields;
        const src = source[0];
        const data = src.data;

        const it = src.keys();
        const key = it.move();
        for (let _f = 0; _f < fields.length; _f++) {
            const f = fields[_f]!;

            StringBuilder.writeSafe(sb, `# ${category.name}.${f.name}: `);
            const val = f.value(key, data, 0);
            StringBuilder.writeSafe(sb, val as string);
            StringBuilder.newline(sb);
        }
        StringBuilder.newline(sb);
    }

    encode(): void {
        // write meta-information, do so after ctab
        if (this.error || this.metaInformation) {
            StringBuilder.writeSafe(this.builder, StringBuilder.getString(this.meta));
        }
        StringBuilder.writeSafe(this.builder, StringBuilder.getString(this.out));

        this.encoded = true;
    }

    constructor(encoder: string, metaInformation: boolean, hydrogens: boolean) {
        super(encoder, metaInformation, hydrogens);
        this.out = StringBuilder.create();
    }
}