// /**
//  * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import { Model, Unit, PropertyGetterProvider } from '../model'
// import * as Data from '../data'
// //import { Vec3 } from '../../utils/linear-algebra'

// // function provider<T>(p: (model: Model) => Property<T>): PropertyProvider<T> {
// //     return p;
// // }

// function atomPropertyProvider<T>(property: (ps: Model['structure']['properties']) => ArrayLike<T>): PropertyGetterProvider<T> {
//     return m => {
//         const a = m.structure.properties.atoms;
//         const p = property(m.structure.properties);
//         return ({ unit, atom }) => p[a[unit][atom]];
//     };
// }

// function unitStructureProvider<T>(p: (structure: Unit.Structure, data: Data.Structure, atom: number) => T): PropertyGetterProvider<T> {
//     return m => {
//         const units = m.structure.units;
//         const data = m.structure.properties.data;
//         return ({ unit, atom }) => p(units[unit], data[unit], atom);
//     };
// }

// export const atom_site = {
//     label_atom_id: atomPropertyProvider(ps => ps.),
//     label_comp_id: unitStructureProvider((u, d, a) => d.residues.name[u.atomResidue[a]])
// }