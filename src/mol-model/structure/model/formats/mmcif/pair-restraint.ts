/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// ihm_predicted_contact_restraint: {
//     id: int,
//     entity_id_1: str,
//     entity_id_2: str,
//     asym_id_1: str,
//     asym_id_2: str,
//     comp_id_1: str,
//     comp_id_2: str,
//     seq_id_1: int,
//     seq_id_2: int,
//     atom_id_1: str,
//     atom_id_2: str,
//     distance_upper_limit: float,
//     probability: float,
//     restraint_type: Aliased<'lower bound' | 'upper bound' | 'lower and upper bound'>(str),
//     model_granularity: Aliased<'by-residue' | 'by-feature' | 'by-atom'>(str),
//     dataset_list_id: int,
//     software_id: int,
// },
// ihm_cross_link_restraint: {
//     id: int,
//     group_id: int,
//     entity_id_1: str,
//     entity_id_2: str,
//     asym_id_1: str,
//     asym_id_2: str,
//     comp_id_1: str,
//     comp_id_2: str,
//     seq_id_1: int,
//     seq_id_2: int,
//     atom_id_1: str,
//     atom_id_2: str,
//     restraint_type: Aliased<'harmonic' | 'upper bound' | 'lower bound'>(str),
//     conditional_crosslink_flag: Aliased<'ALL' | 'ANY'>(str),
//     model_granularity: Aliased<'by-residue' | 'by-feature' | 'by-atom'>(str),
//     distance_threshold: float,
//     psi: float,
//     sigma_1: float,
//     sigma_2: float,
// },

import Model from '../../model'
import { Table } from 'mol-data/db'
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';
import { findAtomIndexByLabelName } from './util';
import { Element, Unit } from '../../../structure';

function findAtomIndex(model: Model, entityId: string, asymId: string, compId: string, seqId: number, atomId: string) {
    if (!model.atomicHierarchy.atoms.auth_atom_id.isDefined) return -1
    const residueIndex = model.atomicHierarchy.findResidueKey(entityId, compId, asymId, seqId, '')
    if (residueIndex < 0) return -1
    return findAtomIndexByLabelName(model, residueIndex, atomId, '') as Element
}

// function findElementIndex(model: Model, entityId: string, asymId: string, compId: string, seqId: number, atomId: string) {

// }

// function key(entityId: string, asymId: string, compId: string, seqId: number, atomId: string) {
//     return `${entityId}|${asymId}|${compId}|${seqId}|${atomId}`
// }

export interface IHMCrossLinkRestraint {
    getIndicesByElement: (element: Element, kind: Unit.Kind) => number[]
    data: Table<mmCIF_Schema['ihm_cross_link_restraint']>
}

export namespace IHMCrossLinkRestraint {
    export const PropName = '__CrossLinkRestraint__';
    export function fromModel(model: Model): IHMCrossLinkRestraint | undefined {
        if (model.properties[PropName]) return model.properties[PropName]

        if (model.sourceData.kind !== 'mmCIF') return
        const { ihm_cross_link_restraint } = model.sourceData.data;
        if (!ihm_cross_link_restraint._rowCount) return

        const p1 = {
            entity_id: ihm_cross_link_restraint.entity_id_1,
            asym_id: ihm_cross_link_restraint.asym_id_1,
            comp_id: ihm_cross_link_restraint.comp_id_1,
            seq_id: ihm_cross_link_restraint.seq_id_1,
            atom_id: ihm_cross_link_restraint.atom_id_1,
        }

        const p2: typeof p1 = {
            entity_id: ihm_cross_link_restraint.entity_id_2,
            asym_id: ihm_cross_link_restraint.asym_id_2,
            comp_id: ihm_cross_link_restraint.comp_id_2,
            seq_id: ihm_cross_link_restraint.seq_id_2,
            atom_id: ihm_cross_link_restraint.atom_id_2,
        }

        function _add(map: Map<Element, number[]>, element: Element, row: number) {
            const indices = map.get(element)
            if (indices) indices.push(row)
            else map.set(element, [ row ])
        }

        function add(row: number, ps: typeof p1) {
            const entityId = ps.entity_id.value(row)
            const asymId = ps.asym_id.value(row)
            const seqId = ps.seq_id.value(row)

            if (ihm_cross_link_restraint.model_granularity.value(row) === 'by-atom') {
                const atomicElement = findAtomIndex(model, entityId, asymId, ps.comp_id.value(row), seqId, ps.atom_id.value(row))
                if (atomicElement >= 0) _add(atomicElementMap, atomicElement as Element, row)
            } else if (model.coarseHierarchy.isDefined) {
                const sphereElement = model.coarseHierarchy.spheres.findSequenceKey(entityId, asymId, seqId)
                if (sphereElement >= 0) {
                    _add(sphereElementMap, sphereElement as Element, row)
                } else {
                    const gaussianElement = model.coarseHierarchy.gaussians.findSequenceKey(entityId, asymId, seqId)
                    if (gaussianElement >= 0) _add(gaussianElementMap, gaussianElement as Element, row)
                }
            }
        }

        function getMapByKind(kind: Unit.Kind) {
            switch (kind) {
                case Unit.Kind.Atomic: return atomicElementMap;
                case Unit.Kind.Spheres: return sphereElementMap;
                case Unit.Kind.Gaussians: return gaussianElementMap;
            }
        }

        /** map from atomic element to cross link indices */
        const atomicElementMap: Map<Element, number[]> = new Map()
        /** map from sphere element to cross link indices */
        const sphereElementMap: Map<Element, number[]> = new Map()
        /** map from gaussian element to cross link indices */
        const gaussianElementMap: Map<Element, number[]> = new Map()

        const emptyIndexArray: number[] = [];

        for (let i = 0; i < ihm_cross_link_restraint._rowCount; ++i) {
            add(i, p1)
            add(i, p2)
        }

        const crossLinkRestraint = {
            getIndicesByElement: (element: Element, kind: Unit.Kind) => {
                const map = getMapByKind(kind)
                const idx = map.get(element)
                return idx !== undefined ? idx : emptyIndexArray
            },
            data: ihm_cross_link_restraint
        }
        model.properties[PropName] = crossLinkRestraint
        return crossLinkRestraint
    }
}