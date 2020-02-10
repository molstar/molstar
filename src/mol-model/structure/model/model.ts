/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import UUID from '../../../mol-util/uuid';
import StructureSequence from './properties/sequence';
import { AtomicHierarchy, AtomicConformation, AtomicRanges } from './properties/atomic';
import { ModelSymmetry } from './properties/symmetry';
import { CoarseHierarchy, CoarseConformation } from './properties/coarse';
import { Entities, ChemicalComponentMap, MissingResidues } from './properties/common';
import { CustomProperties } from '../common/custom-property';
import { SecondaryStructure } from './properties/seconday-structure';
import { SaccharideComponentMap } from '../structure/carbohydrates/constants';
import { ModelFormat } from '../../../mol-model-formats/structure/format';
import { calcModelCenter } from './util';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Mutable } from '../../../mol-util/type-helpers';
import { Coordinates } from '../coordinates';
import { Topology } from '../topology';
import { _parse_mmCif } from '../../../mol-model-formats/structure/mmcif/parser';
import { Task } from '../../../mol-task';
import { IndexPairBonds } from '../../../mol-model-formats/structure/mmcif/bonds/index-pair';

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
export interface Model extends Readonly<{
    id: UUID,
    entryId: string,
    label: string,

    /** the name of the entry/file/collection the model is part of */
    entry: string,

    /**
     * corresponds to
     * - for IHM: `ihm_model_list.model_id`
     * - for standard mmCIF: `atom_site.pdbx_PDB_model_num`
     * - for models from coordinates: frame index
     */
    modelNum: number,

    sourceData: ModelFormat,

    symmetry: ModelSymmetry,
    entities: Entities,
    sequence: StructureSequence,

    atomicHierarchy: AtomicHierarchy,
    atomicConformation: AtomicConformation,
    atomicRanges: AtomicRanges,

    properties: {
        /** secondary structure provided by the input file */
        readonly secondaryStructure: SecondaryStructure,
        /** maps modified residue name to its parent */
        readonly modifiedResidues: Readonly<{
            parentId: ReadonlyMap<string, string>,
            details: ReadonlyMap<string, string>
        }>,
        /** map that holds details about unobserved or zero occurrence residues */
        readonly missingResidues: MissingResidues,
        /** maps residue name to `ChemicalComponent` data */
        readonly chemicalComponentMap: ChemicalComponentMap
        /** maps residue name to `SaccharideComponent` data */
        readonly saccharideComponentMap: SaccharideComponentMap
    },

    customProperties: CustomProperties,

    /**
     * Not to be accessed directly, each custom property descriptor
     * defines property accessors that use this field to store the data.
     */
    _staticPropertyData: { [name: string]: any },
    _dynamicPropertyData: { [name: string]: any },

    coarseHierarchy: CoarseHierarchy,
    coarseConformation: CoarseConformation
}> {

} { }

export namespace Model {
    // TODO: is this enough?
    export type Trajectory = ReadonlyArray<Model>

    export function trajectoryFromModelAndCoordinates(model: Model, coordinates: Coordinates): Trajectory {
        const trajectory: Mutable<Model.Trajectory> = []
        const { frames } = coordinates
        for (let i = 0, il = frames.length; i < il; ++i) {
            const f = frames[i]
            const m = {
                ...model,
                id: UUID.create22(),
                modelNum: i,
                atomicConformation: Coordinates.getAtomicConformation(f, model.atomicConformation.atomId),
                // TODO: add support for supplying sphere and gaussian coordinates in addition to atomic coordinates?
                // coarseConformation: coarse.conformation,
                customProperties: new CustomProperties(),
                _staticPropertyData: Object.create(null),
                _dynamicPropertyData: Object.create(null)
            }
            trajectory.push(m)
        }
        return trajectory
    }

    export function trajectoryFromTopologyAndCoordinates(topology: Topology, coordinates: Coordinates): Task<Trajectory> {
        return Task.create('Create Trajectory', async ctx => {
            const model = (await _parse_mmCif(topology.format, ctx))[0];
            if (!model) throw new Error('found no model')
            const trajectory = trajectoryFromModelAndCoordinates(model, coordinates)
            const bondData = { pairs: topology.bonds, count: model.atomicHierarchy.atoms._rowCount }
            for (const m of trajectory) {
                IndexPairBonds.attachFromData(m, bondData)
            }
            return trajectory
        })
    }

    const CenterProp = '__Center__'
    export function getCenter(model: Model): Vec3 {
        if (model._dynamicPropertyData[CenterProp]) return model._dynamicPropertyData[CenterProp]
        const center = calcModelCenter(model.atomicConformation, model.coarseConformation)
        model._dynamicPropertyData[CenterProp] = center
        return center
    }
}