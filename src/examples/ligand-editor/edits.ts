/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { JSONCifLigandGraph, JSONCifLigandGraphAtom, JSONCifLigandGraphBondProps } from '../../extensions/json-cif/ligand-graph';
import { Vec3 } from '../../mol-math/linear-algebra';
import { VdwRadius } from '../../mol-model/structure/model/properties/atomic';
import { ElementSymbol } from '../../mol-model/structure/model/types';
import { attachRGroup, RGroupName } from './r-groups';

export const TopologyEdits = {
    setElement: async (graph: JSONCifLigandGraph, atomIds: number[], type_symbol: string) => {
        for (const id of atomIds) {
            graph.modifyAtom(id, { type_symbol });
        }
    },
    addElement: async (graph: JSONCifLigandGraph, parentId: number, type_symbol: string) => {
        const p = graph.getAtom(parentId);
        if (!p) return;

        const c = graph.getAtomCoords(p);
        const dir = approximateAddAtomDirection(graph, p);
        const r = 2 / 5 * (VdwRadius(ElementSymbol(p.row.type_symbol ?? 'C')) + VdwRadius(ElementSymbol(type_symbol)));
        const newAtom = graph.addAtom({
            ...p.row,
            // NOTE: this is not correct for editing protein atoms
            // as they should have atom names from CCD, or at least the should be
            // unique. This should be fine for small ligand editing.
            auth_atom_id: type_symbol,
            label_atom_id: type_symbol,
            type_symbol,
            Cartn_x: c[0] + dir[0] * r,
            Cartn_y: c[1] + dir[1] * r,
            Cartn_z: c[2] + dir[2] * r
        });
        graph.addOrUpdateBond(p, newAtom, { value_order: 'sing', type_id: 'covale' });
        return newAtom;
    },
    removeAtoms: async (graph: JSONCifLigandGraph, atomIds: number[]) => {
        for (const id of atomIds) {
            graph.removeAtom(id);
        }
    },
    removeBonds: async (graph: JSONCifLigandGraph, atomIds: number[]) => {
        for (let i = 0; i < atomIds.length; ++i) {
            for (let j = i + 1; j < atomIds.length; ++j) {
                graph.removeBond(atomIds[i], atomIds[j]);
            }
        }
    },
    updateBonds: async (graph: JSONCifLigandGraph, atomIds: number[], props: JSONCifLigandGraphBondProps) => {
        // TODO: iterate on the all-pairs behavior
        // e.g. only add bonds if there is no path connecting them,
        // or by a distance threshold, ...
        for (let i = 0; i < atomIds.length; ++i) {
            for (let j = i + 1; j < atomIds.length; ++j) {
                graph.addOrUpdateBond(atomIds[i], atomIds[j], props);
            }
        }
    },
    attachRgroup: async (graph: JSONCifLigandGraph, atomId: number, name: RGroupName) => {
        await attachRGroup(graph, name, atomId);
    }
};

function approximateAddAtomDirection(graph: JSONCifLigandGraph, parent: JSONCifLigandGraphAtom) {
    let deltas: Vec3[] = [];
    const bonds = graph.bondByKey.get(parent.key);
    if (!bonds?.length) return Vec3.create(1, 0, 0);

    const c = graph.getAtomCoords(parent);
    for (const b of bonds) {
        const delta = Vec3.sub(Vec3(), graph.getAtomCoords(b.atom_2), c);
        deltas.push(delta);
    }

    if (deltas.length === 1) {
        const ret = Vec3.negate(Vec3(), deltas[0]);
        Vec3.normalize(ret, ret);
        return ret;
    }

    if (deltas.length === 2) {
        const ret = Vec3.add(Vec3(), deltas[0], deltas[1]);
        Vec3.normalize(ret, ret);
        Vec3.negate(ret, ret);
        return ret;
    }

    // Take the first three deltas and cross-product them
    deltas = deltas.slice(0, 3);
    const crossProducts: Vec3[] = [];
    for (let i = 0; i < deltas.length; ++i) {
        for (let j = i + 1; j < deltas.length; ++j) {
            const cross = Vec3.cross(Vec3(), deltas[i], deltas[j]);
            Vec3.normalize(cross, cross);
            crossProducts.push(cross);
        }
    }
    for (let i = 1; i < crossProducts.length; ++i) {
        Vec3.matchDirection(crossProducts[i], crossProducts[i], crossProducts[0]);
    }

    const avg = Vec3.create(0, 0, 0);
    for (const cp of crossProducts) {
        Vec3.add(avg, avg, cp);
    }
    Vec3.normalize(avg, avg);
    return avg;
}