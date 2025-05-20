/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { JSONCifLigandGraph, JSONCifLigandGraphAtom, JSONCifLigandGraphBondProps } from '../../extensions/json-cif/ligand-graph';
import { Quat, Vec3 } from '../../mol-math/linear-algebra';
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

export type GeometryEditFn = (param: number) => JSONCifLigandGraph;

export const GeometryEdits = {
    twist: (graph: JSONCifLigandGraph, atomIds: number[]): GeometryEditFn => {
        if (atomIds.length !== 2) {
            throw new Error('Twist requires exactly two atoms.');
        }

        const { left, right } = splitGraph(graph, atomIds[0], atomIds[1]);

        const active = left.length <= right.length ? left : right;
        const a = left.length <= right.length ? atomIds[0] : atomIds[1];
        const b = left.length <= right.length ? atomIds[1] : atomIds[0];

        const pivot = graph.getAtomCoords(a);
        const axis = Vec3.sub(Vec3(), pivot, graph.getAtomCoords(b));
        Vec3.normalize(axis, axis);

        const basePositions = active.map(a => graph.getAtomCoords(a));
        const xform = Quat();
        const p = Vec3();

        return (angle: number) => {
            Quat.setAxisAngle(xform, axis, angle);
            for (let i = 0; i < active.length; ++i) {
                Vec3.copy(p, basePositions[i]);
                Vec3.sub(p, p, pivot);
                Vec3.transformQuat(p, p, xform);
                Vec3.add(p, p, pivot);
                graph.modifyAtom(active[i], {
                    Cartn_x: p[0],
                    Cartn_y: p[1],
                    Cartn_z: p[2]
                });
            }
            return graph;
        };
    },
    stretch: (graph: JSONCifLigandGraph, atomIds: number[]): GeometryEditFn => {
        if (atomIds.length !== 2) {
            throw new Error('Stretch requires exactly two atoms.');
        }

        const { left, right } = splitGraph(graph, atomIds[0], atomIds[1]);

        const a = graph.getAtomCoords(atomIds[0]);
        const b = graph.getAtomCoords(atomIds[1]);
        const center = Vec3.add(Vec3(), b, a);
        Vec3.scale(center, center, 0.5);
        const baseDelta = Vec3.sub(Vec3(), a, center);
        const baseLeft = left.map(a => graph.getAtomCoords(a));
        const baseRight = right.map(a => graph.getAtomCoords(a));

        const p = Vec3();
        const delta = Vec3();

        return (factor: number) => {
            Vec3.scale(delta, baseDelta, factor);
            for (let i = 0; i < left.length; ++i) {
                Vec3.copy(p, baseLeft[i]);
                Vec3.add(p, p, delta);
                graph.modifyAtom(left[i], {
                    Cartn_x: p[0],
                    Cartn_y: p[1],
                    Cartn_z: p[2]
                });
            }
            for (let i = 0; i < right.length; ++i) {
                Vec3.copy(p, baseRight[i]);
                Vec3.sub(p, p, delta);
                graph.modifyAtom(right[i], {
                    Cartn_x: p[0],
                    Cartn_y: p[1],
                    Cartn_z: p[2]
                });
            }
            return graph;
        };
    },
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

function getAtomDepths(graph: JSONCifLigandGraph, atomId: number) {
    return graph.traverse(atomId, 'bfs', new Map<string, number>(), (a, depths, pred) => {
        depths.set(a.key, pred ? depths.get(pred.atom_1.key)! + 1 : 0);
    });
}

function splitGraph(graph: JSONCifLigandGraph, leftId: number, rightId: number) {
    const xs = getAtomDepths(graph, leftId);
    const ys = getAtomDepths(graph, rightId);

    const l: JSONCifLigandGraphAtom[] = [];
    const r: JSONCifLigandGraphAtom[] = [];
    for (const a of graph.atoms) {
        if (xs.has(a.key) && ys.has(a.key)) {
            if (xs.get(a.key)! < ys.get(a.key)!) l.push(a);
            else r.push(a);
        } else if (xs.has(a.key)) {
            l.push(a);
        } else if (ys.has(a.key)) {
            r.push(a);
        }
    }

    return { left: l, right: r };
}