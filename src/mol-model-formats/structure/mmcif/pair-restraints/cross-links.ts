/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../../mol-model/structure/model/model'
import { Table } from '../../../../mol-data/db'
import { mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';
import { findAtomIndexByLabelName } from '../util';
import { Unit } from '../../../../mol-model/structure';
import { ElementIndex } from '../../../../mol-model/structure/model/indexing';

function findAtomIndex(model: Model, entityId: string, asymId: string, seqId: number, atomId: string) {
    if (!model.atomicHierarchy.atoms.auth_atom_id.isDefined) return -1
    const residueIndex = model.atomicHierarchy.index.findResidue(entityId, asymId, seqId)
    if (residueIndex < 0) return -1
    return findAtomIndexByLabelName(model, residueIndex, atomId, '') as ElementIndex
}

export interface IHMCrossLinkRestraint {
    getIndicesByElement: (element: ElementIndex, kind: Unit.Kind) => number[]
    data: Table<mmCIF_Schema['ihm_cross_link_restraint']>
}

export namespace IHMCrossLinkRestraint {
    export const PropName = '__CrossLinkRestraint__';
    export function fromModel(model: Model): IHMCrossLinkRestraint | undefined {
        if (model._staticPropertyData[PropName]) return model._staticPropertyData[PropName]

        if (model.sourceData.kind !== 'mmCIF') return
        const { ihm_cross_link_restraint } = model.sourceData.data;
        if (!ihm_cross_link_restraint._rowCount) return

        const p1 = {
            entity_id: ihm_cross_link_restraint.entity_id_1,
            asym_id: ihm_cross_link_restraint.asym_id_1,
            seq_id: ihm_cross_link_restraint.seq_id_1,
            atom_id: ihm_cross_link_restraint.atom_id_1,
        }

        const p2: typeof p1 = {
            entity_id: ihm_cross_link_restraint.entity_id_2,
            asym_id: ihm_cross_link_restraint.asym_id_2,
            seq_id: ihm_cross_link_restraint.seq_id_2,
            atom_id: ihm_cross_link_restraint.atom_id_2,
        }

        function _add(map: Map<ElementIndex, number[]>, element: ElementIndex, row: number) {
            const indices = map.get(element)
            if (indices) indices.push(row)
            else map.set(element, [ row ])
        }

        function add(row: number, ps: typeof p1) {
            const entityId = ps.entity_id.value(row)
            const asymId = ps.asym_id.value(row)
            const seqId = ps.seq_id.value(row)

            if (ihm_cross_link_restraint.model_granularity.value(row) === 'by-atom') {
                const atomicElement = findAtomIndex(model, entityId, asymId, seqId, ps.atom_id.value(row))
                if (atomicElement >= 0) _add(atomicElementMap, atomicElement as ElementIndex, row)
            } else if (model.coarseHierarchy.isDefined) {
                const sphereElement = model.coarseHierarchy.spheres.findSequenceKey(entityId, asymId, seqId)
                if (sphereElement >= 0) {
                    _add(sphereElementMap, sphereElement, row)
                } else {
                    const gaussianElement = model.coarseHierarchy.gaussians.findSequenceKey(entityId, asymId, seqId)
                    if (gaussianElement >= 0) _add(gaussianElementMap, gaussianElement, row)
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
        const atomicElementMap: Map<ElementIndex, number[]> = new Map()
        /** map from sphere element to cross link indices */
        const sphereElementMap: Map<ElementIndex, number[]> = new Map()
        /** map from gaussian element to cross link indices */
        const gaussianElementMap: Map<ElementIndex, number[]> = new Map()

        const emptyIndexArray: number[] = [];

        for (let i = 0; i < ihm_cross_link_restraint._rowCount; ++i) {
            add(i, p1)
            add(i, p2)
        }

        const crossLinkRestraint = {
            getIndicesByElement: (element: ElementIndex, kind: Unit.Kind) => {
                const map = getMapByKind(kind)
                const idx = map.get(element)
                return idx !== undefined ? idx : emptyIndexArray
            },
            data: ihm_cross_link_restraint
        }
        model._staticPropertyData[PropName] = crossLinkRestraint
        return crossLinkRestraint
    }
}