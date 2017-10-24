// /**
//  * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import { Selector } from '../model'
// //import * as Data from '../data'
// //import { Vec3 } from '../../utils/linear-algebra'

// export const atom = {
//     name: Selector((m, u) => {
//         const unit = m.units[u].structure;
//         const atoms = unit.atoms;
//         const name = unit.data.atoms.name;
//         return a => name[atoms[a]];
//     }),
//     // position: Selector<Vec3>((m, u) => {
//     //     const unit = m.units[u].structure;
//     //     const atoms = unit.atoms;
//     //     const { positions: { x, y, z } } = m.conformation[u];
//     //     const { operator: { transform } } = unit;
//     //     return (atom, position) => {
//     //         const a = atoms[atom];
//     //         const p = position || Vec3.zero();
//     //         Vec3.set(p, x[a], y[a], z[a]);
//     //         return Vec3.transformMat4(p, p, transform);
//     //     };
//     // }),
//     // inversePosition: Selector<Vec3>((m, u) => {
//     //     const unit = m.structure[u];
//     //     const atoms = unit.atoms;
//     //     const { positions: { x, y, z } } = m.conformation[u];
//     //     const { operator: { inverse } } = unit;

//     //     return (atom, position) => {
//     //         const a = atoms[atom];
//     //         const p = position || Vec3.zero();
//     //         Vec3.set(p, x[a], y[a], z[a]);
//     //         return Vec3.transformMat4(p, p, inverse);
//     //     };
//     // })
// }

// export const residue = {

// }