// TODO: make this work when the time comes.
// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { Column, Table } from 'mol-data/db';
// import { Interval, Segmentation } from 'mol-data/int';
// import { mmCIF_Schema as mmCIF } from 'mol-io/reader/cif/schema/mmcif';
// import { Atoms } from 'mol-io/reader/gro/schema';
// import UUID from 'mol-util/uuid';
// import Format from '../format';
// import Model from '../model';
// import { AtomicConformation, AtomicData, AtomicSegments, AtomsSchema, ChainsSchema, ResiduesSchema } from '../properties/atomic';
// import { CoarseHierarchy } from '../properties/coarse';
// import { Entities } from '../properties/common';
// import Sequence from '../properties/sequence';
// import { ModelSymmetry } from '../properties/symmetry';
// import { guessElement } from '../properties/utils/guess-element';
// import { getAtomicKeys } from '../properties/utils/keys';
// import { ElementSymbol } from '../types';

// import gro_Format = Format.gro

// type HierarchyOffsets = { residues: ArrayLike<number>, chains: ArrayLike<number> }

// function findHierarchyOffsets(atomsData: Atoms, bounds: Interval) {
//     const start = Interval.start(bounds), end = Interval.end(bounds);
//     const residues = [start], chains = [start];

//     const { residueName, residueNumber } = atomsData;

//     for (let i = start + 1; i < end; i++) {
//         const newResidue = !residueNumber.areValuesEqual(i - 1, i)
//             || !residueName.areValuesEqual(i - 1, i);
//         console.log(residueName.value(i - 1), residueName.value(i), residueNumber.value(i - 1), residueNumber.value(i), newResidue)
//         if (newResidue) residues[residues.length] = i;
//     }
//     console.log(residues, residues.length)
//     return { residues, chains };
// }

// function guessElementSymbol (value: string) {
//     return ElementSymbol(guessElement(value));
// }

// function createHierarchyData(atomsData: Atoms, offsets: HierarchyOffsets): AtomicData {
//     console.log(atomsData.atomName)
//     const atoms = Table.ofColumns(AtomsSchema, {
//         type_symbol: Column.ofArray({ array: Column.mapToArray(atomsData.atomName, guessElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
//         label_atom_id: atomsData.atomName,
//         auth_atom_id: atomsData.atomName,
//         label_alt_id: Column.Undefined(atomsData.count, Column.Schema.str),
//         pdbx_formal_charge: Column.Undefined(atomsData.count, Column.Schema.int)
//     });

//     const residues = Table.view(Table.ofColumns(ResiduesSchema, {
//         group_PDB: Column.Undefined(atomsData.count, Column.Schema.Aliased<'ATOM' | 'HETATM'>(Column.Schema.str)),
//         label_comp_id: atomsData.residueName,
//         auth_comp_id: atomsData.residueName,
//         label_seq_id: atomsData.residueNumber,
//         auth_seq_id: atomsData.residueNumber,
//         pdbx_PDB_ins_code: Column.Undefined(atomsData.count, Column.Schema.str),
//     }), ResiduesSchema, offsets.residues);
//     // Optimize the numeric columns
//     Table.columnToArray(residues, 'label_seq_id', Int32Array);
//     Table.columnToArray(residues, 'auth_seq_id', Int32Array);

//     // const chains = Table.ofColumns(Hierarchy.ChainsSchema, {
//     //     label_asym_id: Column.ofConst('A', atomsData.count, Column.Schema.str),
//     //     auth_asym_id: Column.ofConst('A', atomsData.count, Column.Schema.str),
//     //     label_entity_id: Column.Undefined(atomsData.count, Column.Schema.str)
//     // });

//     const chains = Table.ofUndefinedColumns(ChainsSchema, 0);

//     return { atoms, residues, chains };
// }

// function getConformation(atoms: Atoms): AtomicConformation {
//     return {
//         id: UUID.create(),
//         atomId: atoms.atomNumber,
//         occupancy: Column.Undefined(atoms.count, Column.Schema.int),
//         B_iso_or_equiv: Column.Undefined(atoms.count, Column.Schema.float),
//         x: Column.mapToArray(atoms.x, x => x * 10, Float32Array),
//         y: Column.mapToArray(atoms.y, y => y * 10, Float32Array),
//         z: Column.mapToArray(atoms.z, z => z * 10, Float32Array)
//     }
// }

// function isHierarchyDataEqual(a: AtomicData, b: AtomicData) {
//     // need to cast because of how TS handles type resolution for interfaces https://github.com/Microsoft/TypeScript/issues/15300
//     return Table.areEqual(a.residues as Table<ResiduesSchema>, b.residues as Table<ResiduesSchema>)
//         && Table.areEqual(a.atoms as Table<AtomsSchema>, b.atoms as Table<AtomsSchema>)
// }

// function createModel(format: gro_Format, modelNum: number, previous?: Model): Model {
//     const structure = format.data.structures[modelNum];
//     const bounds = Interval.ofBounds(0, structure.atoms.count);

//     const hierarchyOffsets = findHierarchyOffsets(structure.atoms, bounds);
//     const hierarchyData = createHierarchyData(structure.atoms, hierarchyOffsets);

//     if (previous && isHierarchyDataEqual(previous.atomicHierarchy, hierarchyData)) {
//         return {
//             ...previous,
//             atomicConformation: getConformation(structure.atoms)
//         };
//     }

//     const hierarchySegments: AtomicSegments = {
//         residueSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, bounds),
//         chainSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, bounds),
//     }

//     // TODO: create a better mock entity
//     const entityTable = Table.ofRows<mmCIF['entity']>(mmCIF.entity, [{
//         id: '0',
//         src_method: 'syn',
//         type: 'polymer',
//         pdbx_number_of_molecules: 1
//     }]);

//     const entities: Entities = { data: entityTable, getEntityIndex: Column.createIndexer(entityTable.id) };

//     const hierarchyKeys = getAtomicKeys(hierarchyData, entities, hierarchySegments);
//     const atomicHierarchy = { ...hierarchyData, ...hierarchyKeys, ...hierarchySegments };
//     return {
//         id: UUID.create(),
//         sourceData: format,
//         modelNum,
//         atomicHierarchy,
//         entities,
//         sequence: Sequence.fromAtomicHierarchy(atomicHierarchy),
//         atomicConformation: getConformation(structure.atoms),
//         coarseHierarchy: CoarseHierarchy.Empty,
//         coarseConformation: void 0 as any,
//         symmetry: ModelSymmetry.Default
//     };
// }

// function buildModels(format: gro_Format): ReadonlyArray<Model> {
//     const models: Model[] = [];

//     format.data.structures.forEach((_, i) => {
//         const model = createModel(format, i, models.length > 0 ? models[models.length - 1] : void 0);
//         models.push(model);
//     });
//     return models;
// }

// export default buildModels;
