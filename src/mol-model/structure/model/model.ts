/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import UUID from '../../../mol-util/uuid';
import StructureSequence from './properties/sequence';
import { AtomicHierarchy, AtomicConformation, AtomicRanges } from './properties/atomic';
import { CoarseHierarchy, CoarseConformation } from './properties/coarse';
import { Entities, ChemicalComponentMap, MissingResidues, StructAsymMap } from './properties/common';
import { CustomProperties } from '../common/custom-property';
import { SaccharideComponentMap } from '../structure/carbohydrates/constants';
import { ModelFormat } from '../../../mol-model-formats/structure/format';
import { calcModelCenter } from './util';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Mutable } from '../../../mol-util/type-helpers';
import { Coordinates } from '../coordinates';
import { Topology } from '../topology';
import { Task } from '../../../mol-task';
import { IndexPairBonds } from '../../../mol-model-formats/structure/property/bonds/index-pair';
import { createModels } from '../../../mol-model-formats/structure/basic/parser';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { ChainIndex } from './indexing';
import { SymmetryOperator } from '../../../mol-math/geometry';

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

    /**
     * This is a hack to allow "model-index coloring"
     */
    trajectoryInfo: {
        index: number,
        size: number
    },

    sourceData: ModelFormat,

    entities: Entities,
    sequence: StructureSequence,

    atomicHierarchy: AtomicHierarchy,
    atomicConformation: AtomicConformation,
    atomicRanges: AtomicRanges,
    atomicChainOperatorMappinng: Map<ChainIndex, SymmetryOperator>,

    properties: {
        /** map that holds details about unobserved or zero occurrence residues */
        readonly missingResidues: MissingResidues,
        /** maps residue name to `ChemicalComponent` data */
        readonly chemicalComponentMap: ChemicalComponentMap
        /** maps residue name to `SaccharideComponent` data */
        readonly saccharideComponentMap: SaccharideComponentMap
        /** maps label_asym_id name to `StructAsym` data */
        readonly structAsymMap: StructAsymMap
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
            const model = (await createModels(topology.basic, topology.sourceData, ctx))[0];
            if (!model) throw new Error('found no model')
            const trajectory = trajectoryFromModelAndCoordinates(model, coordinates)
            const bondData = { pairs: topology.bonds, count: model.atomicHierarchy.atoms._rowCount }
            const indexPairBonds = IndexPairBonds.fromData(bondData)

            let index = 0;
            for (const m of trajectory) {
                IndexPairBonds.Provider.set(m, indexPairBonds)
                m.trajectoryInfo.index = index++;
                m.trajectoryInfo.size = trajectory.length;
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

    //

    export function isFromPdbArchive(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        return (
            db.database_2.database_id.isDefined
        )
    }

    export function hasSecondaryStructure(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        return (
            db.struct_conf.id.isDefined ||
            db.struct_sheet_range.id.isDefined
        )
    }

    export function hasCrystalSymmetry(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        return (
            db.symmetry._rowCount === 1 && db.cell._rowCount === 1 && !(
                db.symmetry.Int_Tables_number.value(0) === 1 &&
                db.cell.angle_alpha.value(0) === 90 &&
                db.cell.angle_beta.value(0) === 90 &&
                db.cell.angle_gamma.value(0) === 90 &&
                db.cell.length_a.value(0) === 1 &&
                db.cell.length_b.value(0) === 1 &&
                db.cell.length_c.value(0) === 1
            )
        )
    }

    export function isFromXray(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        for (let i = 0; i < db.exptl.method.rowCount; i++) {
            const v = db.exptl.method.value(i).toUpperCase()
            if (v.indexOf('DIFFRACTION') >= 0) return true
        }
        return false
    }

    export function isFromEm(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        for (let i = 0; i < db.exptl.method.rowCount; i++) {
            const v = db.exptl.method.value(i).toUpperCase()
            if (v.indexOf('MICROSCOPY') >= 0) return true
        }
        return false
    }

    export function isFromNmr(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        for (let i = 0; i < db.exptl.method.rowCount; i++) {
            const v = db.exptl.method.value(i).toUpperCase()
            if (v.indexOf('NMR') >= 0) return true
        }
        return false
    }

    export function hasXrayMap(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        // Check exprimental method to exclude models solved with
        // 'ELECTRON CRYSTALLOGRAPHY' which also have structure factors
        if (!isFromXray(model)) return false
        const { db } = model.sourceData.data
        const { status_code_sf } = db.pdbx_database_status
        return status_code_sf.isDefined && status_code_sf.value(0) === 'REL'
    }

    /**
     * Also checks for `content_type` of 'associated EM volume' to exclude cases
     * like 6TEK which are solved with 'X-RAY DIFFRACTION' but have an related
     * EMDB entry of type 'other EM volume'.
     */
    export function hasEmMap(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false
        const { db } = model.sourceData.data
        const { db_name, content_type } = db.pdbx_database_related
        for (let i = 0, il = db.pdbx_database_related._rowCount; i < il; ++i) {
            if (db_name.value(i).toUpperCase() === 'EMDB' && content_type.value(i) === 'associated EM volume') {
                return true
            }
        }
        return false
    }

    export function hasDensityMap(model: Model) {
        return hasXrayMap(model) || hasEmMap(model)
    }
}