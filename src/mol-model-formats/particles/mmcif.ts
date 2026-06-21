/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { ParticleList } from '../../mol-model/particles/particle-list';
import { CustomProperties } from '../../mol-model/custom-property';
import { CifBlock, CifFile } from '../../mol-io/reader/cif/data-model';
import { toDatabase } from '../../mol-io/reader/cif/schema';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { getMatrices, parseOperatorList, expandOperators } from '../structure/property/assembly';
import { ModelFormat } from '../format';

export type MmcifVariant = 'auto' | 'cellpack' | 'petworld' | 'standard';

export interface MmcifAssemblyParticleListOptions {
    /** Assembly identifier from `_pdbx_struct_assembly.id`. */
    readonly assemblyId: string
    /**
     * Individual chain IDs (`label_asym_id`) to include.
     * An empty or omitted array means "include all chains for this assembly".
     */
    readonly asymIds?: ReadonlyArray<string>
    readonly label?: string
    /**
     * Which mmCIF variant to use for building the particle list.
     * - `'auto'` (default): detect automatically from the CIF block header and categories.
     * - `'cellpack'`: interpret `_entity.pdbx_description` as a dot-separated `compartment.entityName` path.
     * - `'standard'`: like cellpack but treat `_entity.pdbx_description` as a plain entity name only.
     * - `'petworld'`: read the non-standard `_pdbx_struct_assembly_gen.PDB_model_num` field;
     *   one particle per operator combo per gen row.
     */
    readonly variant?: MmcifVariant
}

/**
 * Return all unique assembly IDs found in `_pdbx_struct_assembly_gen` of the first CIF block.
 * Results are sorted lexicographically.
 */
export function getAssemblyIdsFromMmcif(cifFile: CifFile): string[] {
    const block = cifFile.blocks[0];
    if (!block) return [];
    const db = toDatabase(mmCIF_Schema, block);
    const { assembly_id } = db.pdbx_struct_assembly_gen;
    const ids = new Set<string>();
    for (let i = 0, il = db.pdbx_struct_assembly_gen._rowCount; i < il; i++) {
        ids.add(assembly_id.value(i));
    }
    return Array.from(ids).sort();
}

/**
 * Return the sorted unique individual chain IDs for the given assembly in the first CIF block.
 */
export function getAsymIdsFromMmcif(cifFile: CifFile, assemblyId: string): string[] {
    const block = cifFile.blocks[0];
    if (!block) return [];
    const db = toDatabase(mmCIF_Schema, block);
    const { assembly_id, asym_id_list } = db.pdbx_struct_assembly_gen;
    const ids = new Set<string>();
    for (let i = 0, il = db.pdbx_struct_assembly_gen._rowCount; i < il; i++) {
        if (assembly_id.value(i) !== assemblyId) continue;
        for (const chain of asym_id_list.value(i)) ids.add(chain);
    }
    return Array.from(ids).sort();
}

export function createParticleListFromMmcifAssembly(cifFile: CifFile, options: MmcifAssemblyParticleListOptions): ParticleList {
    const block = cifFile.blocks[0];
    if (!block) throw new Error('CIF file contains no data blocks.');
    const variant = resolveVariant(block, options.variant);
    if (variant === 'petworld') return buildPetworldParticleList(cifFile, block, options);
    return buildCellpackStandardParticleList(cifFile, block, options, variant);
}

function detectMmcifVariant(block: CifBlock): 'cellpack' | 'petworld' | 'standard' {
    if (block.header.toLowerCase().includes('cellpack')) return 'cellpack';
    if (block.categories['pdbx_struct_assembly_gen']?.getField('PDB_model_num')?.isDefined) return 'petworld';
    return 'standard';
}

function resolveVariant(block: CifBlock, requested?: MmcifVariant): 'cellpack' | 'petworld' | 'standard' {
    if (!requested || requested === 'auto') return detectMmcifVariant(block);
    return requested;
}

function buildCellpackStandardParticleList(
    cifFile: CifFile,
    block: CifBlock,
    options: MmcifAssemblyParticleListOptions,
    variant: 'cellpack' | 'standard'
): ParticleList {
    const db = toDatabase(mmCIF_Schema, block);
    const { pdbx_struct_assembly_gen, pdbx_struct_oper_list } = db;

    if (pdbx_struct_oper_list._rowCount === 0) {
        throw new Error('CIF file contains no _pdbx_struct_oper_list entries.');
    }

    const matrices = getMatrices(pdbx_struct_oper_list);

    const { assembly_id, oper_expression, asym_id_list } = pdbx_struct_assembly_gen;
    const asymFilter = options.asymIds?.length
        ? new Set<string>(options.asymIds)
        : undefined;

    // Collect matching gen rows and their expanded operator combinations.
    interface GenEntry {
        genIndex: number
        operCombinations: string[][] // each combination is an array of operator IDs
    }

    const entries: GenEntry[] = [];
    let totalCount = 0;

    for (let i = 0, il = pdbx_struct_assembly_gen._rowCount; i < il; i++) {
        if (assembly_id.value(i) !== options.assemblyId) continue;

        const operList = parseOperatorList(oper_expression.value(i));
        const combinations = expandOperators(operList);
        entries.push({ genIndex: i, operCombinations: combinations });
        const chainCount = asymFilter
            ? asym_id_list.value(i).filter(c => asymFilter.has(c)).length
            : asym_id_list.value(i).length;
        totalCount += chainCount * combinations.length;
    }

    if (totalCount === 0) {
        if (entries.length === 0) {
            throw new Error(`No _pdbx_struct_assembly_gen rows found for assembly '${options.assemblyId}'.`);
        }
        throw new Error(
            asymFilter !== undefined
                ? `No chains matched assembly '${options.assemblyId}' and asym IDs '${options.asymIds!.join(', ')}'.`
                : `Assembly '${options.assemblyId}' has no expanded operators.`
        );
    }

    // Build entity_id → compartment name and entity name from _entity.pdbx_description.
    // In the 'cellpack' variant, descriptions are dot-separated paths:
    //   compartment = all segments except the last  (e.g. "root.mge.surface.proteins")
    //   entity name = last segment                  (e.g. "MG_191_192_NAP")
    // In the 'standard' variant, the full description is used as the entity name with no compartment.
    const entityToCompartment = new Map<string, string>();
    const entityToName = new Map<string, string>();
    for (let i = 0, il = db.entity._rowCount; i < il; i++) {
        const entityId = db.entity.id.value(i);
        const desc = db.entity.pdbx_description.value(i).join(',');
        if (variant === 'cellpack') {
            const dotIdx = desc.lastIndexOf('.');
            if (dotIdx > 0) {
                entityToCompartment.set(entityId, desc.substring(0, dotIdx));
                entityToName.set(entityId, desc.substring(dotIdx + 1));
            } else if (desc) {
                entityToName.set(entityId, desc);
            }
        } else {
            if (desc) entityToName.set(entityId, desc);
        }
    }

    // Build unique compartment name → index map (populated lazily below).
    const compartmentNameToIdx = new Map<string, number>();
    // Build unique entity name → index map (populated lazily below).
    const entityNameToIdx = new Map<string, number>();

    const keys = new Int32Array(totalCount);
    const targets = new Int32Array(totalCount);
    const compartments = new Int32Array(totalCount).fill(-1);
    const entities = new Int32Array(totalCount).fill(-1);
    const coordinates = new Float32Array(totalCount * 3);
    const rotations = new Float32Array(totalCount * 4);
    const radii = new Float32Array(totalCount);

    // Per-particle metadata for labels (stored as flat arrays to avoid object allocation).
    const labelChainId = new Array<string>(totalCount);
    const labelOpCombo = new Array<string>(totalCount);

    // Build per-chain centroid from _atom_site Cartesian coordinates and collect
    // chain → entity_id mapping at the same time (avoids relying on _struct_asym).
    // Applying the full operator to the centroid of the reference chains gives a meaningful
    // particle position even when the operator has zero translation (pure rotation).
    const { label_asym_id: siteAsymId, label_entity_id: siteEntityId, Cartn_x, Cartn_y, Cartn_z } = db.atom_site;
    interface ChainAccum { x: number; y: number; z: number; n: number; radiusSq: number }
    const chainAccum = new Map<string, ChainAccum>();
    const chainToEntityId = new Map<string, string>();
    for (let i = 0, il = db.atom_site._rowCount; i < il; i++) {
        const chain = siteAsymId.value(i);
        let acc = chainAccum.get(chain);
        if (!acc) {
            acc = { x: 0, y: 0, z: 0, n: 0, radiusSq: 0 };
            chainAccum.set(chain, acc);
            chainToEntityId.set(chain, siteEntityId.value(i));
        }
        acc.x += Cartn_x.value(i);
        acc.y += Cartn_y.value(i);
        acc.z += Cartn_z.value(i);
        acc.n++;
    }
    // Second pass: compute per-chain max squared distance from centroid (= bounding sphere radius²).
    for (let i = 0, il = db.atom_site._rowCount; i < il; i++) {
        const chain = siteAsymId.value(i);
        const acc = chainAccum.get(chain)!;
        const cx = acc.x / acc.n, cy = acc.y / acc.n, cz = acc.z / acc.n;
        const dx = Cartn_x.value(i) - cx;
        const dy = Cartn_y.value(i) - cy;
        const dz = Cartn_z.value(i) - cz;
        const d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > acc.radiusSq) acc.radiusSq = d2;
    }

    const combined = Mat4();
    const quaternion = Quat();
    const centroid = Vec3();
    const position = Vec3();

    // Assign a unique target index per distinct chain ID.
    const chainToTargetIdx = new Map<string, number>();
    let nextTargetIdx = 0;
    for (let ei = 0, eil = entries.length; ei < eil; ei++) {
        for (const chain of asym_id_list.value(entries[ei].genIndex)) {
            if (asymFilter && !asymFilter.has(chain)) continue;
            if (!chainToTargetIdx.has(chain)) chainToTargetIdx.set(chain, nextTargetIdx++);
        }
    }

    let count = 0;
    for (let ei = 0, eil = entries.length; ei < eil; ei++) {
        const entry = entries[ei];
        const chainGroup = asym_id_list.value(entry.genIndex);

        for (const chain of chainGroup) {
            if (asymFilter && !asymFilter.has(chain)) continue;
            // Per-chain centroid.
            const acc = chainAccum.get(chain);
            Vec3.set(centroid,
                acc && acc.n > 0 ? acc.x / acc.n : 0,
                acc && acc.n > 0 ? acc.y / acc.n : 0,
                acc && acc.n > 0 ? acc.z / acc.n : 0
            );

            for (const combo of entry.operCombinations) {
                // Multiply matrices left-to-right: first op in combo is outermost.
                Mat4.setIdentity(combined);
                for (let k = 0; k < combo.length; k++) {
                    const m = matrices.get(combo[k]);
                    if (!m) throw new Error(`Operator '${combo[k]}' not found in _pdbx_struct_oper_list.`);
                    Mat4.mul(combined, combined, m);
                }

                Vec3.transformMat4(position, centroid, combined);
                Quat.normalize(quaternion, Quat.fromMat4(quaternion, combined));

                const cOffset = count * 3;
                coordinates[cOffset + 0] = position[0];
                coordinates[cOffset + 1] = position[1];
                coordinates[cOffset + 2] = position[2];

                const qOffset = count * 4;
                rotations[qOffset + 0] = quaternion[0];
                rotations[qOffset + 1] = quaternion[1];
                rotations[qOffset + 2] = quaternion[2];
                rotations[qOffset + 3] = quaternion[3];

                radii[count] = Math.sqrt(acc ? acc.radiusSq : 0) * Mat4.getMaxScaleOnAxis(combined);

                keys[count] = count;
                targets[count] = chainToTargetIdx.get(chain)!;

                const entityId = chainToEntityId.get(chain);
                const compartmentName = entityId !== undefined ? entityToCompartment.get(entityId) : undefined;
                if (compartmentName !== undefined) {
                    if (!compartmentNameToIdx.has(compartmentName)) {
                        compartmentNameToIdx.set(compartmentName, compartmentNameToIdx.size);
                    }
                    compartments[count] = compartmentNameToIdx.get(compartmentName)!;
                }

                const entityName = entityId !== undefined ? entityToName.get(entityId) : undefined;
                if (entityName !== undefined) {
                    if (!entityNameToIdx.has(entityName)) {
                        entityNameToIdx.set(entityName, entityNameToIdx.size);
                    }
                    entities[count] = entityNameToIdx.get(entityName)!;
                }

                labelChainId[count] = chain;
                labelOpCombo[count] = combo.join('×');
                count++;
            }
        }
    }

    // Build targetMapping: target index → single chain ID.
    const targetMapping = new Map<number, ReadonlyArray<string>>();
    for (const [chain, idx] of chainToTargetIdx) {
        targetMapping.set(idx, [chain]);
    }

    // Build compartmentInfo: compartment index → compartment name.
    const compartmentInfo = new Map<number, string>();
    for (const [name, idx] of compartmentNameToIdx) {
        compartmentInfo.set(idx, name);
    }

    // Build entityInfo: entity index → entity name.
    const entityInfo = new Map<number, string>();
    for (const [name, idx] of entityNameToIdx) {
        entityInfo.set(idx, name);
    }

    const assemblyId = options.assemblyId;

    return {
        label: buildMmcifLabel(options.label, assemblyId),
        count,
        keys,
        targets,
        targetMapping,
        compartments: compartmentInfo.size > 0 ? compartments : undefined,
        compartmentInfo: compartmentInfo.size > 0 ? compartmentInfo : undefined,
        entities: entityInfo.size > 0 ? entities : undefined,
        entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
        coordinates,
        rotations,
        radii,
        getParticleLabel: (index: number) => {
            const chain = labelChainId[index];
            const opCombo = labelOpCombo[index];
            return `#${index + 1} | assembly ${assemblyId} | chain ${chain} | ops ${opCombo}`;
        },
        sourceData: MmcifParticleFormat.create(cifFile),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

function buildPetworldParticleList(
    cifFile: CifFile,
    block: CifBlock,
    options: MmcifAssemblyParticleListOptions
): ParticleList {
    const db = toDatabase(mmCIF_Schema, block);
    const { pdbx_struct_assembly_gen, pdbx_struct_oper_list } = db;

    if (pdbx_struct_oper_list._rowCount === 0) {
        throw new Error('CIF file contains no _pdbx_struct_oper_list entries.');
    }

    const matrices = getMatrices(pdbx_struct_oper_list);
    const { assembly_id, oper_expression } = pdbx_struct_assembly_gen;

    // PDB_model_num is a non-standard field absent from the typed mmCIF_Schema; read it via raw categories.
    const pdbModelNumField = block.categories['pdbx_struct_assembly_gen'].getField('PDB_model_num')!;

    // Read _pdbx_model.name for entity names (non-standard category; rows are in model order, 1-indexed).
    const pdbxModelName = block.categories['pdbx_model']?.getField('name');
    const modelNumToName = new Map<number, string>();
    if (pdbxModelName) {
        for (let i = 0, il = pdbxModelName.rowCount; i < il; i++) {
            modelNumToName.set(i + 1, pdbxModelName.str(i));
        }
    }

    // Collect matching gen rows.
    // One particle per operator combo per gen row = one particle per complex instance.
    interface GenEntry {
        genIndex: number
        modelNum: number
        operCombinations: string[][]
    }

    const entries: GenEntry[] = [];
    let totalCount = 0;

    for (let i = 0, il = pdbx_struct_assembly_gen._rowCount; i < il; i++) {
        if (assembly_id.value(i) !== options.assemblyId) continue;
        const operList = parseOperatorList(oper_expression.value(i));
        const combinations = expandOperators(operList);
        const modelNum = pdbModelNumField.int(i);
        entries.push({ genIndex: i, modelNum, operCombinations: combinations });
        totalCount += combinations.length;
    }

    if (totalCount === 0) {
        if (entries.length === 0) {
            throw new Error(`No _pdbx_struct_assembly_gen rows found for assembly '${options.assemblyId}'.`);
        }
        throw new Error(`Assembly '${options.assemblyId}' has no expanded operators.`);
    }

    // Map distinct PDB_model_num values (ascending) → 0-based trajectory model index.
    // mol* `TrajectoryFromMmCif` creates one model per distinct pdbx_PDB_model_num in ascending order,
    // so the rank of a model number equals its trajectory frame index.
    const distinctModelNums = Array.from(new Set(entries.map(e => e.modelNum))).sort((a, b) => a - b);
    const modelNumToIndex = new Map<number, number>();
    distinctModelNums.forEach((mn, idx) => modelNumToIndex.set(mn, idx));

    // targetModels: target ID (= trajectory model index) → trajectory model index (identity).
    // Each target's reference structure is the full structure of that trajectory model.
    const targetModels = new Map<number, number>();
    for (const idx of modelNumToIndex.values()) {
        targetModels.set(idx, idx);
    }

    // Build per-model centroid from _atom_site (aggregate all atoms with matching pdbx_PDB_model_num).
    const { Cartn_x, Cartn_y, Cartn_z, pdbx_PDB_model_num: siteModelNum } = db.atom_site;
    interface ModelAccum { x: number; y: number; z: number; n: number; radiusSq: number }
    const modelAccum = new Map<number, ModelAccum>();
    for (let i = 0, il = db.atom_site._rowCount; i < il; i++) {
        const modelNum = siteModelNum.value(i);
        let acc = modelAccum.get(modelNum);
        if (!acc) {
            acc = { x: 0, y: 0, z: 0, n: 0, radiusSq: 0 };
            modelAccum.set(modelNum, acc);
        }
        acc.x += Cartn_x.value(i);
        acc.y += Cartn_y.value(i);
        acc.z += Cartn_z.value(i);
        acc.n++;
    }
    // Second pass: compute per-model max squared distance from centroid (= bounding sphere radius²).
    for (let i = 0, il = db.atom_site._rowCount; i < il; i++) {
        const modelNum = siteModelNum.value(i);
        const acc = modelAccum.get(modelNum)!;
        const cx = acc.x / acc.n, cy = acc.y / acc.n, cz = acc.z / acc.n;
        const dx = Cartn_x.value(i) - cx;
        const dy = Cartn_y.value(i) - cy;
        const dz = Cartn_z.value(i) - cz;
        const d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > acc.radiusSq) acc.radiusSq = d2;
    }

    const keys = new Int32Array(totalCount);
    const targets = new Int32Array(totalCount);
    const entities = new Int32Array(totalCount).fill(-1);
    const coordinates = new Float32Array(totalCount * 3);
    const rotations = new Float32Array(totalCount * 4);
    const radii = new Float32Array(totalCount);

    const labelModelNum = new Int32Array(totalCount);
    const labelOpCombo = new Array<string>(totalCount);

    const entityNameToIdx = new Map<string, number>();

    const combined = Mat4();
    const quaternion = Quat();
    const centroid = Vec3();
    const position = Vec3();

    let count = 0;
    for (let ei = 0, eil = entries.length; ei < eil; ei++) {
        const entry = entries[ei];
        const targetIdx = modelNumToIndex.get(entry.modelNum)!; // 0-based trajectory model index
        const modelName = modelNumToName.get(entry.modelNum);

        const acc = modelAccum.get(entry.modelNum);
        Vec3.set(centroid,
            acc && acc.n > 0 ? acc.x / acc.n : 0,
            acc && acc.n > 0 ? acc.y / acc.n : 0,
            acc && acc.n > 0 ? acc.z / acc.n : 0
        );

        for (const combo of entry.operCombinations) {
            // Multiply matrices left-to-right: first op in combo is outermost.
            Mat4.setIdentity(combined);
            for (let k = 0; k < combo.length; k++) {
                const m = matrices.get(combo[k]);
                if (!m) throw new Error(`Operator '${combo[k]}' not found in _pdbx_struct_oper_list.`);
                Mat4.mul(combined, combined, m);
            }

            Vec3.transformMat4(position, centroid, combined);
            Quat.normalize(quaternion, Quat.fromMat4(quaternion, combined));

            const cOffset = count * 3;
            coordinates[cOffset + 0] = position[0];
            coordinates[cOffset + 1] = position[1];
            coordinates[cOffset + 2] = position[2];

            const qOffset = count * 4;
            rotations[qOffset + 0] = quaternion[0];
            rotations[qOffset + 1] = quaternion[1];
            rotations[qOffset + 2] = quaternion[2];
            rotations[qOffset + 3] = quaternion[3];

            radii[count] = Math.sqrt(acc ? acc.radiusSq : 0) * Mat4.getMaxScaleOnAxis(combined);

            keys[count] = count;
            targets[count] = targetIdx;

            if (modelName !== undefined) {
                if (!entityNameToIdx.has(modelName)) entityNameToIdx.set(modelName, entityNameToIdx.size);
                entities[count] = entityNameToIdx.get(modelName)!;
            }

            labelModelNum[count] = entry.modelNum;
            labelOpCombo[count] = combo.join('×');
            count++;
        }
    }

    // Build entityInfo: entity index → model name.
    const entityInfo = new Map<number, string>();
    for (const [name, idx] of entityNameToIdx) {
        entityInfo.set(idx, name);
    }

    const assemblyId = options.assemblyId;

    return {
        label: buildMmcifLabel(options.label, assemblyId),
        count,
        keys,
        targets,
        targetModels,
        entities: entityInfo.size > 0 ? entities : undefined,
        entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
        coordinates,
        rotations,
        radii,
        getParticleLabel: (index: number) => {
            const modelNum = labelModelNum[index];
            const opCombo = labelOpCombo[index];
            return `#${index + 1} | assembly ${assemblyId} | model ${modelNum} | ops ${opCombo}`;
        },
        sourceData: MmcifParticleFormat.create(cifFile),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

function buildMmcifLabel(label?: string, assemblyId?: string): string {
    const base = label || 'Particles';
    return assemblyId ? `${base} (assembly ${assemblyId})` : base;
}

//

export { MmcifParticleFormat };

type MmcifParticleFormat = ModelFormat<CifFile>

namespace MmcifParticleFormat {
    export function is(x?: ModelFormat): x is MmcifParticleFormat {
        return x?.kind === 'mmcif-particle';
    }

    export function create(cifFile: CifFile): MmcifParticleFormat {
        return { kind: 'mmcif-particle', name: 'mmcif-particle', data: cifFile };
    }
}
