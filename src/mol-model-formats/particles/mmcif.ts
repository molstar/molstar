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
import { mmCIF_Database, mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { parseOperatorList } from '../structure/property/assembly';
import { OperatorCombinations, composeOperatorCombination, expandOperatorCombinations, formatOperatorCombination, readParticleOperators } from './operators';
import { ModelFormat } from '../format';
import { binaryCifHasColumn, getBinaryCifHeader } from '../../mol-io/common/binary-cif';
import { StringLike } from '../../mol-io/common/string-like';
import { FileNameInfo } from '../../mol-util/file-info';
import { RuntimeContext } from '../../mol-task';

export type MmcifVariant = 'auto' | 'cellpack' | 'petworld' | 'standard';

const PETWORLD_CATEGORY = '_pdbx_struct_assembly_gen';
const PETWORLD_COLUMN = 'PDB_model_num';
const PETWORLD_FIELD = `${PETWORLD_CATEGORY}.${PETWORLD_COLUMN}`;

/**
 * Content-sniff a `.cif`/`.bcif` file (before full parsing) to detect whether it looks like
 * a CellPack or PetWorld mmCIF particle assembly, for use as a `DataFormatProvider.isApplicable`.
 * Returns `false` for ordinary/standard mmCIF structure files so those keep loading as
 * regular trajectories by default.
 */
export function looksLikeMmcifParticles(info: FileNameInfo, data: StringLike | Uint8Array): boolean {
    if (info.ext === 'bcif') {
        try {
            const bcifHeader = getBinaryCifHeader(data as Uint8Array);
            const header = bcifHeader.dataBlocks[0]?.header;
            if (header && header.toLowerCase().includes('cellpack')) return true;
            if (binaryCifHasColumn(bcifHeader, PETWORLD_CATEGORY, PETWORLD_COLUMN)) return true;
        } catch (e) {
            console.error(e);
        }
        return false;
    } else if (info.ext === 'cif') {
        const str = data as StringLike;
        const header = str.substring(0, 400);
        if (header.toLowerCase().includes('cellpack')) return true;
        if (str.includes(PETWORLD_FIELD)) return true;
    }
    return false;
}

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
    * - `'cellpack'`: interpret `_entity.pdbx_description` as a dot-separated `compartment.entityName` path;
    *   entities below a `.fibers` compartment use assembly-gen rows as ordered fibers.
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

export async function createParticleListFromMmcifAssembly(cifFile: CifFile, options: MmcifAssemblyParticleListOptions, ctx: RuntimeContext): Promise<ParticleList> {
    const block = cifFile.blocks[0];
    if (!block) throw new Error('CIF file contains no data blocks.');
    const variant = resolveVariant(block, options.variant);
    if (variant === 'petworld') return await buildPetworldParticleList(cifFile, block, options, ctx);
    return await buildCellpackStandardParticleList(cifFile, block, options, variant, ctx);
}

/** The centroid of a group of `_atom_site` rows and their maximum distance from it. */
interface AtomSiteBounds {
    readonly center: Vec3
    readonly radius: number
}

/**
 * Applying the operator of a particle to the centroid of its reference atoms gives a meaningful
 * position even when the operator has zero translation (pure rotation).
 */
async function getAtomSiteBounds<K>(db: mmCIF_Database, key: (row: number) => K, ctx: RuntimeContext): Promise<Map<K, AtomSiteBounds>> {
    const { Cartn_x, Cartn_y, Cartn_z, _rowCount } = db.atom_site;
    const updateChunk = 10000;

    interface Accum { x: number, y: number, z: number, n: number, radiusSq: number }
    const accums = new Map<K, Accum>();

    for (let i = 0; i < _rowCount; i++) {
        const k = key(i);
        let acc = accums.get(k);
        if (!acc) {
            acc = { x: 0, y: 0, z: 0, n: 0, radiusSq: 0 };
            accums.set(k, acc);
        }
        acc.x += Cartn_x.value(i);
        acc.y += Cartn_y.value(i);
        acc.z += Cartn_z.value(i);
        acc.n++;
        if (i % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Calculating particle bounds', current: i, max: _rowCount * 2 });
        }
    }

    for (const acc of accums.values()) {
        acc.x /= acc.n;
        acc.y /= acc.n;
        acc.z /= acc.n;
    }

    for (let i = 0; i < _rowCount; i++) {
        const acc = accums.get(key(i))!;
        const dx = Cartn_x.value(i) - acc.x;
        const dy = Cartn_y.value(i) - acc.y;
        const dz = Cartn_z.value(i) - acc.z;
        const d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > acc.radiusSq) acc.radiusSq = d2;
        if (i % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Calculating particle bounds', current: _rowCount + i, max: _rowCount * 2 });
        }
    }

    const bounds = new Map<K, AtomSiteBounds>();
    for (const [k, acc] of accums) {
        bounds.set(k, { center: Vec3.create(acc.x, acc.y, acc.z), radius: Math.sqrt(acc.radiusSq) });
    }
    return bounds;
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

async function buildCellpackStandardParticleList(
    cifFile: CifFile,
    block: CifBlock,
    options: MmcifAssemblyParticleListOptions,
    variant: 'cellpack' | 'standard',
    ctx: RuntimeContext
): Promise<ParticleList> {
    const db = toDatabase(mmCIF_Schema, block);
    const { pdbx_struct_assembly_gen, pdbx_struct_oper_list } = db;

    if (pdbx_struct_oper_list._rowCount === 0) {
        throw new Error('CIF file contains no _pdbx_struct_oper_list entries.');
    }

    const operators = readParticleOperators(block);

    const { assembly_id, oper_expression, asym_id_list } = pdbx_struct_assembly_gen;
    const asymFilter = options.asymIds?.length
        ? new Set<string>(options.asymIds)
        : undefined;

    // Collect matching gen rows and their expanded operator combinations.
    interface GenEntry {
        genIndex: number
        combinations: OperatorCombinations
    }

    const entries: GenEntry[] = [];
    let totalCount = 0;

    for (let i = 0, il = pdbx_struct_assembly_gen._rowCount; i < il; i++) {
        if (assembly_id.value(i) !== options.assemblyId) continue;

        const operList = parseOperatorList(oper_expression.value(i));
        const combinations = expandOperatorCombinations(operList, operators);
        entries.push({ genIndex: i, combinations });
        const chainCount = asymFilter
            ? asym_id_list.value(i).filter(c => asymFilter.has(c)).length
            : asym_id_list.value(i).length;
        totalCount += chainCount * combinations.count;
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
    // CellPack entities whose penultimate segment is "fibers" use each assembly-gen row as
    // a separate ordered fiber. In the 'standard' variant, the full description is used as
    // the entity name with no compartment.
    const entityToCompartment = new Map<string, string>();
    const entityToName = new Map<string, string>();
    const fiberEntityIds = new Set<string>();
    for (let i = 0, il = db.entity._rowCount; i < il; i++) {
        const entityId = db.entity.id.value(i);
        const desc = db.entity.pdbx_description.value(i).join(',');
        if (variant === 'cellpack') {
            const segments = desc.split('.');
            if (segments.length > 1) {
                entityToCompartment.set(entityId, segments.slice(0, -1).join('.'));
                entityToName.set(entityId, segments[segments.length - 1]);
                if (segments[segments.length - 2] === 'fibers') fiberEntityIds.add(entityId);
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

    const { label_asym_id: siteAsymId, label_entity_id: siteEntityId } = db.atom_site;
    const chainBounds = await getAtomSiteBounds(db, i => siteAsymId.value(i), ctx);
    const chainToEntityId = new Map<string, string>();
    const updateChunk = 10000;
    for (let i = 0, il = db.atom_site._rowCount; i < il; i++) {
        const chain = siteAsymId.value(i);
        if (!chainToEntityId.has(chain)) chainToEntityId.set(chain, siteEntityId.value(i));
        if (i % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Mapping particle chains', current: i, max: il });
        }
    }

    let fiberCount = 0;
    let fiberPointCount = 0;
    for (const entry of entries) {
        if (entry.combinations.count < 2) continue;
        for (const chain of asym_id_list.value(entry.genIndex)) {
            if (asymFilter && !asymFilter.has(chain)) continue;
            const entityId = chainToEntityId.get(chain);
            if (entityId !== undefined && fiberEntityIds.has(entityId)) {
                fiberCount++;
                fiberPointCount += entry.combinations.count;
            }
        }
    }

    const keys = new Int32Array(totalCount);
    const targets = new Int32Array(totalCount);
    const compartments = new Int32Array(totalCount).fill(-1);
    const entities = new Int32Array(totalCount).fill(-1);
    const coordinates = new Float32Array(totalCount * 3);
    const rotations = new Float32Array(totalCount * 4);
    const radii = new Float32Array(totalCount);
    const fiberOffsets = new Int32Array(fiberCount + 1);
    const fiberIndices = new Int32Array(fiberPointCount);

    // Per-particle metadata for labels, kept as indices so no string is allocated per particle.
    const labelEntry = new Int32Array(totalCount);
    const labelCombo = new Int32Array(totalCount);

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
    let fiberIndex = 0;
    let fiberPosition = 0;
    for (let ei = 0, eil = entries.length; ei < eil; ei++) {
        const entry = entries[ei];
        const chainGroup = asym_id_list.value(entry.genIndex);

        for (const chain of chainGroup) {
            if (asymFilter && !asymFilter.has(chain)) continue;
            const bounds = chainBounds.get(chain);
            Vec3.copy(centroid, bounds ? bounds.center : Vec3.origin);
            const entityId = chainToEntityId.get(chain);
            const isFiber = entry.combinations.count >= 2 && entityId !== undefined && fiberEntityIds.has(entityId);

            for (let ci = 0, cil = entry.combinations.count; ci < cil; ci++) {
                composeOperatorCombination(operators, entry.combinations, ci, combined);

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

                radii[count] = (bounds ? bounds.radius : 0) * Mat4.getMaxScaleOnAxis(combined);

                keys[count] = count;
                targets[count] = chainToTargetIdx.get(chain)!;

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

                labelEntry[count] = ei;
                labelCombo[count] = ci;
                if (isFiber) fiberIndices[fiberPosition++] = count;
                count++;
                if (count % updateChunk === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Building particles', current: count, max: totalCount });
                }
            }
            if (isFiber) fiberOffsets[++fiberIndex] = fiberPosition;
        }
    }
    if (ctx.shouldUpdate) {
        await ctx.update({ message: 'Building particles', current: count, max: totalCount });
    }

    const chainMapping = new Map<number, string>();
    for (const [chain, idx] of chainToTargetIdx) chainMapping.set(idx, chain);

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
    // Only the ids are captured below, so the label closure does not retain the matrix array.
    const operatorIds = operators.ids;
    const combinations = entries.map(e => e.combinations);

    return {
        label: buildMmcifLabel(options.label, assemblyId),
        count,
        keys,
        targets,
        compartments: compartmentInfo.size > 0 ? compartments : undefined,
        compartmentInfo: compartmentInfo.size > 0 ? compartmentInfo : undefined,
        entities: entityInfo.size > 0 ? entities : undefined,
        entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
        coordinates,
        rotations,
        radii,
        fibers: fiberCount > 0 ? { count: fiberCount, offsets: fiberOffsets, indices: fiberIndices } : undefined,
        getParticleLabel: (index: number) => {
            const entity = entityInfo.get(entities[index]);
            const chain = chainMapping.get(targets[index]);
            const opCombo = formatOperatorCombination(operatorIds, combinations[labelEntry[index]], labelCombo[index]);
            return `#${index + 1} | ${entity} | chain ${chain} | ops ${opCombo}`;
        },
        sourceData: MmcifParticleFormat.create(cifFile),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

async function buildPetworldParticleList(
    cifFile: CifFile,
    block: CifBlock,
    options: MmcifAssemblyParticleListOptions,
    ctx: RuntimeContext
): Promise<ParticleList> {
    const db = toDatabase(mmCIF_Schema, block);
    const { pdbx_struct_assembly_gen, pdbx_struct_oper_list } = db;

    if (pdbx_struct_oper_list._rowCount === 0) {
        throw new Error('CIF file contains no _pdbx_struct_oper_list entries.');
    }

    const operators = readParticleOperators(block);
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
        combinations: OperatorCombinations
    }

    const entries: GenEntry[] = [];
    let totalCount = 0;

    for (let i = 0, il = pdbx_struct_assembly_gen._rowCount; i < il; i++) {
        if (assembly_id.value(i) !== options.assemblyId) continue;
        const operList = parseOperatorList(oper_expression.value(i));
        const combinations = expandOperatorCombinations(operList, operators);
        const modelNum = pdbModelNumField.int(i);
        entries.push({ genIndex: i, modelNum, combinations });
        totalCount += combinations.count;
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

    const keys = new Int32Array(totalCount);
    const targets = new Int32Array(totalCount);
    const entities = new Int32Array(totalCount).fill(-1);
    const coordinates = new Float32Array(totalCount * 3);
    const rotations = new Float32Array(totalCount * 4);
    const radii = new Float32Array(totalCount);

    const labelEntry = new Int32Array(totalCount);
    const labelCombo = new Int32Array(totalCount);

    const entityNameToIdx = new Map<string, number>();

    const { pdbx_PDB_model_num: siteModelNum } = db.atom_site;
    const modelBounds = await getAtomSiteBounds(db, i => siteModelNum.value(i), ctx);

    const combined = Mat4();
    const quaternion = Quat();
    const centroid = Vec3();
    const position = Vec3();
    const updateChunk = 10000;

    let count = 0;
    for (let ei = 0, eil = entries.length; ei < eil; ei++) {
        const entry = entries[ei];
        const targetIdx = modelNumToIndex.get(entry.modelNum)!; // 0-based trajectory model index
        const modelName = modelNumToName.get(entry.modelNum) || `Model ${entry.modelNum}`;

        const bounds = modelBounds.get(entry.modelNum);
        Vec3.copy(centroid, bounds ? bounds.center : Vec3.origin);

        for (let ci = 0, cil = entry.combinations.count; ci < cil; ci++) {
            composeOperatorCombination(operators, entry.combinations, ci, combined);

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

            radii[count] = (bounds ? bounds.radius : 0) * Mat4.getMaxScaleOnAxis(combined);

            keys[count] = count;
            targets[count] = targetIdx;

            if (!entityNameToIdx.has(modelName)) entityNameToIdx.set(modelName, entityNameToIdx.size);
            entities[count] = entityNameToIdx.get(modelName)!;

            labelEntry[count] = ei;
            labelCombo[count] = ci;
            count++;
            if (count % updateChunk === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'Building particles', current: count, max: totalCount });
            }
        }
    }
    if (ctx.shouldUpdate) {
        await ctx.update({ message: 'Building particles', current: count, max: totalCount });
    }

    // Build entityInfo: entity index → model name.
    const entityInfo = new Map<number, string>();
    for (const [name, idx] of entityNameToIdx) {
        entityInfo.set(idx, name);
    }

    const assemblyId = options.assemblyId;
    // Only the ids are captured below, so the label closure does not retain the matrix array.
    const operatorIds = operators.ids;
    const combinations = entries.map(e => e.combinations);

    return {
        label: buildMmcifLabel(options.label, assemblyId),
        count,
        keys,
        targets,
        entities: entityInfo.size > 0 ? entities : undefined,
        entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
        coordinates,
        rotations,
        radii,
        getParticleLabel: (index: number) => {
            const entity = entityInfo.get(entities[index])!;
            const opCombo = formatOperatorCombination(operatorIds, combinations[labelEntry[index]], labelCombo[index]);
            return `#${index + 1} | ${entity} | ops ${opCombo}`;
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
