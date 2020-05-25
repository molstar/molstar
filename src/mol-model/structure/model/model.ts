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
import { CustomProperties } from '../../custom-property';
import { SaccharideComponentMap } from '../structure/carbohydrates/constants';
import { ModelFormat } from '../../../mol-model-formats/format';
import { calcModelCenter, getAsymIdCount } from './util';
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
import { ModelSymmetry } from '../../../mol-model-formats/structure/property/symmetry';
import { Column } from '../../../mol-data/db';
import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';

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

    function _trajectoryFromModelAndCoordinates(model: Model, coordinates: Coordinates) {
        const trajectory: Mutable<Model.Trajectory> = [];
        const { frames } = coordinates;

        const srcIndex = model.atomicHierarchy.atoms.sourceIndex;
        const isIdentity = Column.isIdentity(srcIndex);
        const srcIndexArray = isIdentity ? void 0 : srcIndex.toArray({ array: Int32Array });

        for (let i = 0, il = frames.length; i < il; ++i) {
            const f = frames[i];
            const m = {
                ...model,
                id: UUID.create22(),
                modelNum: i,
                atomicConformation: isIdentity
                    ? Coordinates.getAtomicConformation(f, model.atomicConformation.atomId)
                    : Coordinates.getAtomicConformationReordered(f, model.atomicConformation.atomId, srcIndexArray!),
                // TODO: add support for supplying sphere and gaussian coordinates in addition to atomic coordinates?
                // coarseConformation: coarse.conformation,
                customProperties: new CustomProperties(),
                _staticPropertyData: Object.create(null),
                _dynamicPropertyData: Object.create(null)
            };
            trajectory.push(m);
        }
        return { trajectory, srcIndexArray };
    }

    export function trajectoryFromModelAndCoordinates(model: Model, coordinates: Coordinates): Trajectory {
        return _trajectoryFromModelAndCoordinates(model, coordinates).trajectory;
    }

    function invertIndex(xs: ArrayLike<number>) {
        const ret = new Int32Array(xs.length);
        for (let i = 0, _i = xs.length; i < _i; i++) {
            ret[xs[i]] = i;
        }
        return ret;
    }

    export function trajectoryFromTopologyAndCoordinates(topology: Topology, coordinates: Coordinates): Task<Trajectory> {
        return Task.create('Create Trajectory', async ctx => {
            const model = (await createModels(topology.basic, topology.sourceData, ctx))[0];
            if (!model) throw new Error('found no model');
            const { trajectory, srcIndexArray } = _trajectoryFromModelAndCoordinates(model, coordinates);

            // TODO: cache the inverted index somewhere?
            const invertedIndex = srcIndexArray ? invertIndex(srcIndexArray) : void 0;
            const pairs = srcIndexArray
                ? {
                    indexA: Column.ofIntArray(Column.mapToArray(topology.bonds.indexA, i => invertedIndex![i], Int32Array)),
                    indexB: Column.ofIntArray(Column.mapToArray(topology.bonds.indexB, i => invertedIndex![i], Int32Array)),
                    order: topology.bonds.order
                }
                : topology.bonds;

            const bondData = { pairs, count: model.atomicHierarchy.atoms._rowCount };
            const indexPairBonds = IndexPairBonds.fromData(bondData);

            let index = 0;
            for (const m of trajectory) {
                IndexPairBonds.Provider.set(m, indexPairBonds);
                TrajectoryInfo.set(m, { index: index++, size: trajectory.length });
            }
            return trajectory;
        });
    }

    const CenterProp = '__Center__';
    export function getCenter(model: Model): Vec3 {
        if (model._dynamicPropertyData[CenterProp]) return model._dynamicPropertyData[CenterProp];
        const center = calcModelCenter(model.atomicConformation, model.coarseConformation);
        model._dynamicPropertyData[CenterProp] = center;
        return center;
    }

    const TrajectoryInfoProp = '__TrajectoryInfo__';
    export type TrajectoryInfo = { readonly index: number, readonly size: number }
    export const TrajectoryInfo = {
        get(model: Model): TrajectoryInfo {
            return model._staticPropertyData[TrajectoryInfoProp] || { index: 0, size: 1 };
        },
        set(model: Model, trajectoryInfo: TrajectoryInfo) {
            return model._staticPropertyData[TrajectoryInfoProp] = trajectoryInfo;
        }
    };

    const AsymIdCountProp = '__AsymIdCount__';
    export type AsymIdCount = { readonly auth: number, readonly label: number }
    export const AsymIdCount = {
        get(model: Model): AsymIdCount {
            if (model._staticPropertyData[AsymIdCountProp]) return model._staticPropertyData[AsymIdCountProp];
            const asymIdCount = getAsymIdCount(model);
            model._staticPropertyData[AsymIdCountProp] = asymIdCount;
            return asymIdCount;
        },
    };

    export type AsymIdOffset = { auth: number, label: number };
    export const AsymIdOffset = CustomModelProperty.createSimple<AsymIdOffset>('asym_id_offset', 'static');

    export type Index = number;
    export const Index = CustomModelProperty.createSimple<Index>('index', 'static');

    //

    export function hasCarbohydrate(model: Model): boolean {
        return model.properties.saccharideComponentMap.size > 0;
    }

    export function hasProtein(model: Model): boolean {
        const { subtype } = model.entities;
        for (let i = 0, il = subtype.rowCount; i < il; ++i) {
            if (subtype.value(i).startsWith('polypeptide')) return true;
        }
        return false;
    }

    export function hasNucleic(model: Model): boolean {
        const { subtype } = model.entities;
        for (let i = 0, il = subtype.rowCount; i < il; ++i) {
            const s = subtype.value(i);
            if (s.endsWith('ribonucleotide hybrid') || s.endsWith('ribonucleotide')) return true;
        }
        return false;
    }

    export function isFromPdbArchive(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        return (
            db.database_2.database_id.isDefined ||
            // 4 character PDB id
            model.entryId.match(/^[1-9][a-z0-9]{3,3}$/i) !== null ||
            // long PDB id
            model.entryId.match(/^pdb_[0-9]{4,4}[1-9][a-z0-9]{3,3}$/i) !== null
        );
    }

    export function hasSecondaryStructure(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        return (
            db.struct_conf.id.isDefined ||
            db.struct_sheet_range.id.isDefined
        );
    }

    const tmpAngles90 = Vec3.create(1.5707963, 1.5707963, 1.5707963); // in radians
    const tmpLengths1 = Vec3.create(1, 1, 1);
    export function hasCrystalSymmetry(model: Model): boolean {
        const spacegroup = ModelSymmetry.Provider.get(model)?.spacegroup;
        return !!spacegroup && !(
            spacegroup.num === 1 &&
            Vec3.equals(spacegroup.cell.anglesInRadians, tmpAngles90) &&
            Vec3.equals(spacegroup.cell.size, tmpLengths1)
        );
    }

    export function isFromXray(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        for (let i = 0; i < db.exptl.method.rowCount; i++) {
            const v = db.exptl.method.value(i).toUpperCase();
            if (v.indexOf('DIFFRACTION') >= 0) return true;
        }
        return false;
    }

    export function isFromEm(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        for (let i = 0; i < db.exptl.method.rowCount; i++) {
            const v = db.exptl.method.value(i).toUpperCase();
            if (v.indexOf('MICROSCOPY') >= 0) return true;
        }
        return false;
    }

    export function isFromNmr(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        for (let i = 0; i < db.exptl.method.rowCount; i++) {
            const v = db.exptl.method.value(i).toUpperCase();
            if (v.indexOf('NMR') >= 0) return true;
        }
        return false;
    }

    export function hasXrayMap(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        // Check exprimental method to exclude models solved with
        // 'ELECTRON CRYSTALLOGRAPHY' which also have structure factors
        if (!isFromXray(model)) return false;
        const { db } = model.sourceData.data;
        const { status_code_sf } = db.pdbx_database_status;
        return status_code_sf.isDefined && status_code_sf.value(0) === 'REL';
    }

    /**
     * Also checks for `content_type` of 'associated EM volume' to exclude cases
     * like 6TEK which are solved with 'X-RAY DIFFRACTION' but have an related
     * EMDB entry of type 'other EM volume'.
     */
    export function hasEmMap(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        const { db_name, content_type } = db.pdbx_database_related;
        for (let i = 0, il = db.pdbx_database_related._rowCount; i < il; ++i) {
            if (db_name.value(i).toUpperCase() === 'EMDB' && content_type.value(i) === 'associated EM volume') {
                return true;
            }
        }
        return false;
    }

    export function hasDensityMap(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        return hasXrayMap(model) || hasEmMap(model);
    }

    export function probablyHasDensityMap(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        return hasDensityMap(model) || (
            // check if from pdb archive but missing relevant meta data
            isFromPdbArchive(model) && (
                !db.exptl.method.isDefined ||
                (isFromXray(model) && (
                    !db.pdbx_database_status.status_code_sf.isDefined ||
                    db.pdbx_database_status.status_code_sf.valueKind(0) === Column.ValueKind.Unknown
                )) ||
                (isFromEm(model) && (
                    !db.pdbx_database_related.db_name.isDefined
                ))
            )
        );
    }
}