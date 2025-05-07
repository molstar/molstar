/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { JSONCifLigandGraph, LigandGraphBondProps } from '../../extensions/json-cif/ligand-graph';
import { attachRGroup, RGroupName } from './r-groups';

export const TopologyEdits = {
    setElement: async (graph: JSONCifLigandGraph, atomIds: number[], type_symbol: string) => {
        for (const id of atomIds) {
            graph.modifyAtom(id, { type_symbol });
        }
    },
    addElement: async (graph: JSONCifLigandGraph, parentId: number, type_symbol: string) => {
        const newAtom = graph.attachAtom(parentId, { type_symbol });
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
    updateBonds: async (graph: JSONCifLigandGraph, atomIds: number[], props: LigandGraphBondProps) => {
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