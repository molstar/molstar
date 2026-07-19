/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Quat, Vec3 } from '../../mol-math/linear-algebra';
import { Euler } from '../../mol-math/linear-algebra/3d/euler';
import { Particle, ParticleList } from '../../mol-model/particles/particle-list';
import { ParticleTrajectory } from '../../mol-model/particles/particle-trajectory';
import { CustomProperties } from '../../mol-model/custom-property';
import { SimulariumFile, SimulariumAgentBuffer as AB, SimulariumMinValuesPerAgent as MIN_VALUES_PER_AGENT, SimulariumVisType } from '../../mol-io/reader/simularium/schema';
import { ModelFormat } from '../format';
import { objectForEach } from '../../mol-util/object';
import { Box3D } from '../../mol-math/geometry/primitives/box3d';
import { Sphere3D } from '../../mol-math/geometry/primitives/sphere3d';
import { Boundary } from '../../mol-math/geometry/boundary';
import { Cell } from '../../mol-math/geometry/spacegroup/cell';

export interface SimulariumParticleListOptions {
    /** Index of the trajectory frame to convert into a particle list. */
    readonly frameIndex: number
    /**
     * Spatial scale factor applied to all coordinates (to angstrom). A value `<= 0`
     * (or omitted) auto-derives the factor from `trajectoryInfo.spatialUnits`.
     */
    readonly scale?: number
    /** Agent type ids to include. An empty/omitted list includes all types. */
    readonly typeFilter?: ReadonlyArray<number>
}

/** Conversion factor from a Simularium spatial unit name to angstrom. */
function spatialUnitToAngstrom(name?: string): number {
    if (!name) return 1;
    switch (name.toLowerCase().trim()) {
        case 'm': case 'meter': case 'meters': return 1e10;
        case 'cm': case 'centimeter': case 'centimeters': return 1e8;
        case 'mm': case 'millimeter': case 'millimeters': return 1e7;
        case 'um': case 'µm': case 'micron': case 'microns': case 'micrometer': case 'micrometers': return 1e4;
        case 'nm': case 'nanometer': case 'nanometers': return 10;
        case 'a': case 'å': case 'ang': case 'angstrom': case 'angstroms': return 1;
        case 'pm': case 'picometer': case 'picometers': return 0.01;
        default: return 1;
    }
}

function resolveScale(file: SimulariumFile, scale?: number): number {
    if (scale && scale > 0) return scale;
    const units = file.trajectoryInfo.spatialUnits;
    return (units?.magnitude ?? 1) * spatialUnitToAngstrom(units?.name);
}

/** Number of particles a single agent contributes (and whether it forms a fiber). */
function agentParticleCount(visType: number, nSubpoints: number): { particles: number, fiberPoints: number } {
    if (visType === SimulariumVisType.FIBER) {
        const fiberPoints = Math.floor(nSubpoints / 3);
        if (fiberPoints >= 1) return { particles: fiberPoints, fiberPoints };
    }
    return { particles: 1, fiberPoints: 0 };
}

export function createParticleListFromSimularium(file: SimulariumFile, options: SimulariumParticleListOptions): ParticleList {
    const { trajectoryInfo, frames } = file;
    if (frames.length === 0) throw new Error('Simularium file contains no frames.');

    const frameIndex = Math.max(0, Math.min(options.frameIndex | 0, frames.length - 1));
    const frame = frames[frameIndex];

    const scale = resolveScale(file, options.scale);
    const typeFilter = options.typeFilter && options.typeFilter.length > 0
        ? new Set<number>(options.typeFilter)
        : undefined;

    const buf = frame.data;
    const bufLen = buf.length;

    // Pass 1: determine particle, fiber, and fiber-point counts.
    let particleCount = 0;
    let fiberCount = 0;
    let fiberPointCount = 0;
    for (let j = 0; j + MIN_VALUES_PER_AGENT <= bufLen;) {
        const visType = buf[j + AB.VIS_TYPE];
        const typeId = buf[j + AB.TYPE_ID] | 0;
        const nSub = buf[j + AB.N_SUBPOINTS] | 0;
        if (!typeFilter || typeFilter.has(typeId)) {
            const { particles, fiberPoints } = agentParticleCount(visType, nSub);
            particleCount += particles;
            if (fiberPoints > 0) { fiberCount += 1; fiberPointCount += fiberPoints; }
        }
        j += MIN_VALUES_PER_AGENT + nSub;
    }

    if (particleCount === 0) {
        throw new Error(typeFilter
            ? 'No Simularium agents matched the selected agent types.'
            : 'No Simularium agents were found in the selected frame.');
    }

    const keys = new Int32Array(particleCount);
    const targets = new Int32Array(particleCount);
    const coordinates = new Float32Array(particleCount * 3);
    const rotations = new Float32Array(particleCount * 4);
    const radii = new Float32Array(particleCount);
    const entities = new Int32Array(particleCount);

    const fiberOffsets = new Int32Array(fiberCount + 1);
    const fiberIndices = new Int32Array(fiberPointCount);

    // Per-particle source info for labels.
    const particleTypeId = new Int32Array(particleCount);
    const particleInstanceId = new Float64Array(particleCount);

    const typeIdToEntityIdx = new Map<number, number>();
    const entityInfo = new Map<number, string>();

    if (trajectoryInfo.typeMapping) {
        objectForEach(trajectoryInfo.typeMapping, (value, key) => {
            const typeId = Number(key);
            const entityIdx = typeIdToEntityIdx.size;
            typeIdToEntityIdx.set(typeId, entityIdx);
            entityInfo.set(entityIdx, value.name);
        });
    }

    const typeName = (t: number) => trajectoryInfo.typeMapping?.[String(t)]?.name ?? `type ${t}`;
    const entityIndexOf = (t: number) => {
        let idx = typeIdToEntityIdx.get(t);
        if (idx === undefined) {
            idx = typeIdToEntityIdx.size;
            typeIdToEntityIdx.set(t, idx);
            entityInfo.set(idx, typeName(t));
        }
        return idx;
    };

    const euler = Euler();
    const quat = Quat();
    const local = Vec3();

    let count = 0;
    let fiberIdx = 0;
    let fiberPos = 0;

    for (let j = 0; j + MIN_VALUES_PER_AGENT <= bufLen;) {
        const visType = buf[j + AB.VIS_TYPE];
        const instanceId = buf[j + AB.INSTANCE_ID];
        const typeId = buf[j + AB.TYPE_ID] | 0;
        const px = buf[j + AB.POS_X], py = buf[j + AB.POS_Y], pz = buf[j + AB.POS_Z];
        const rx = buf[j + AB.ROT_X], ry = buf[j + AB.ROT_Y], rz = buf[j + AB.ROT_Z];
        const radius = buf[j + AB.RADIUS];
        const nSub = buf[j + AB.N_SUBPOINTS] | 0;
        const subStart = j + MIN_VALUES_PER_AGENT;

        if (!typeFilter || typeFilter.has(typeId)) {
            const entityIdx = entityIndexOf(typeId);
            const { fiberPoints } = agentParticleCount(visType, nSub);

            if (fiberPoints > 0) {
                // FIBER agent: explode each subpoint into its own particle (rendered as a sphere).
                euler[0] = rx; euler[1] = ry; euler[2] = rz;
                Quat.fromEuler(quat, euler, 'XYZ');

                fiberOffsets[fiberIdx] = fiberPos;
                for (let s = 0; s < fiberPoints; ++s) {
                    const so = subStart + s * 3;
                    Vec3.set(local, buf[so + 0], buf[so + 1], buf[so + 2]);
                    Vec3.transformQuat(local, local, quat);

                    const cOffset = count * 3;
                    coordinates[cOffset + 0] = (px + local[0]) * scale;
                    coordinates[cOffset + 1] = (py + local[1]) * scale;
                    coordinates[cOffset + 2] = (pz + local[2]) * scale;

                    const qOffset = count * 4;
                    rotations[qOffset + 3] = 1; // identity for subpoint spheres

                    radii[count] = radius * scale;
                    entities[count] = entityIdx;
                    targets[count] = typeId; // target id = agent type id (per-type reference geometry)
                    keys[count] = count;
                    particleTypeId[count] = typeId;
                    particleInstanceId[count] = instanceId;

                    fiberIndices[fiberPos++] = count;
                    ++count;
                }
                fiberOffsets[++fiberIdx] = fiberPos;
            } else {
                // DEFAULT agent (or fiber without subpoints): a single particle at the agent position.
                euler[0] = rx; euler[1] = ry; euler[2] = rz;
                Quat.fromEuler(quat, euler, 'XYZ');

                const cOffset = count * 3;
                coordinates[cOffset + 0] = px * scale;
                coordinates[cOffset + 1] = py * scale;
                coordinates[cOffset + 2] = pz * scale;

                const qOffset = count * 4;
                rotations[qOffset + 0] = quat[0];
                rotations[qOffset + 1] = quat[1];
                rotations[qOffset + 2] = quat[2];
                rotations[qOffset + 3] = quat[3];

                radii[count] = radius * scale;
                entities[count] = entityIdx;
                targets[count] = typeId; // target id = agent type id (per-type reference geometry)
                keys[count] = count;
                particleTypeId[count] = typeId;
                particleInstanceId[count] = instanceId;
                ++count;
            }
        }

        j += MIN_VALUES_PER_AGENT + nSub;
    }

    const title = trajectoryInfo.trajectoryTitle;
    const label = title ? `${title} (frame ${frame.frameNumber})` : `Simularium particles (frame ${frame.frameNumber})`;

    // the simulation box is an orthogonal (rectangular) cell centered on the origin
    const cell = trajectoryInfo.size
        ? Cell.create(
            Vec3.create(trajectoryInfo.size.x * scale, trajectoryInfo.size.y * scale, trajectoryInfo.size.z * scale),
            Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2)
        )
        : undefined;

    const particleList: ParticleList = {
        label,
        count,
        keys,
        targets,
        entities,
        entityInfo,
        coordinates,
        rotations,
        radii,
        cell,
        fibers: fiberCount > 0 ? { count: fiberCount, offsets: fiberOffsets, indices: fiberIndices } : undefined,
        getParticleLabel: (index: number) => {
            const t = particleTypeId[index];
            const id = particleInstanceId[index];
            return `#${index + 1} | ${typeName(t)} | id ${id}`;
        },
        sourceData: SimulariumFormat.create(file),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };

    // set the particle boundary to the size of the simulation box if available
    if (trajectoryInfo.size) {
        const { x, y, z } = trajectoryInfo.size;
        let maxRadius = 0;
        for (let i = 0; i < count; ++i) {
            if (radii[i] > maxRadius) maxRadius = radii[i];
        }
        // box around origin with the given size, scaled to angstrom
        const box = Box3D.create(
            Vec3.create(-x * scale / 2 - maxRadius, -y * scale / 2 - maxRadius, -z * scale / 2 - maxRadius),
            Vec3.create(x * scale / 2 + maxRadius, y * scale / 2 + maxRadius, z * scale / 2 + maxRadius)
        );
        const sphere = Sphere3D.fromBox3D(Sphere3D(), box);
        const boundary: Boundary = { box, sphere };
        Particle.setBoundary(particleList, boundary);
    }

    return particleList;
}

export interface SimulariumParticleTrajectoryOptions {
    /**
     * Spatial scale factor applied to all coordinates (to angstrom). A value `<= 0`
     * (or omitted) auto-derives the factor from `trajectoryInfo.spatialUnits`.
     */
    readonly scale?: number
    /** Agent type ids to include. An empty/omitted list includes all types. */
    readonly typeFilter?: ReadonlyArray<number>
}

/** Lazy `ParticleTrajectory` backed by a `SimulariumFile`. Converts frames on demand. */
class SimulariumParticleTrajectory implements ParticleTrajectory {
    readonly frameCount: number;
    readonly representative: ParticleList;

    getFrameAtIndex(i: number): ParticleList {
        return createParticleListFromSimularium(this.file, {
            frameIndex: i,
            scale: this.options.scale && this.options.scale > 0 ? this.options.scale : void 0,
            typeFilter: this.options.typeFilter,
        });
    }

    constructor(private file: SimulariumFile, private options: SimulariumParticleTrajectoryOptions) {
        this.frameCount = file.frames.length;
        this.representative = this.getFrameAtIndex(0);
    }
}

/** Create a lazy `ParticleTrajectory` from a `SimulariumFile`. Frames are converted on demand. */
export function createSimulariumParticleTrajectory(file: SimulariumFile, options: SimulariumParticleTrajectoryOptions = {}): ParticleTrajectory {
    return new SimulariumParticleTrajectory(file, options);
}

/** Return the number of frames available in the file. */
export function getSimulariumFrameCount(file: SimulariumFile): number {
    return file.frames.length;
}

/** Return the agent types declared in `typeMapping`, sorted by id. */
export function getSimulariumAgentTypeNames(file: SimulariumFile): { id: number, name: string }[] {
    const out: { id: number, name: string }[] = [];
    const typeMapping = file.trajectoryInfo.typeMapping ?? {};
    for (const key of Object.keys(typeMapping)) {
        const id = Number(key);
        if (!Number.isFinite(id)) continue;
        out.push({ id, name: typeMapping[key]?.name ?? `type ${id}` });
    }
    out.sort((a, b) => a.id - b.id);
    return out;
}

/** An external geometry (PDB structure or OBJ mesh) referenced by an agent type. */
export interface SimulariumAgentGeometry {
    readonly id: number
    readonly name: string
    /** Normalized to upper case, e.g. `PDB` or `OBJ`. */
    readonly displayType: string
    readonly url: string
    readonly color?: string
}

/**
 * Return the agent types whose `geometry` references an external PDB structure or
 * OBJ mesh (i.e. `displayType` is `PDB` or `OBJ` and a `url` is given), sorted by id.
 */
export function getSimulariumAgentGeometries(file: SimulariumFile): SimulariumAgentGeometry[] {
    const out: SimulariumAgentGeometry[] = [];
    const typeMapping = file.trajectoryInfo.typeMapping ?? {};
    for (const key of Object.keys(typeMapping)) {
        const id = Number(key);
        if (!Number.isFinite(id)) continue;
        const entry = typeMapping[key];
        const geometry = entry?.geometry;
        if (!geometry?.url) continue;
        const displayType = (geometry.displayType ?? '').toUpperCase();
        if (displayType !== 'PDB' && displayType !== 'OBJ') continue;
        out.push({ id, name: entry?.name ?? `type ${id}`, displayType, url: geometry.url, color: geometry.color });
    }
    out.sort((a, b) => a.id - b.id);
    return out;
}

/**
 * Return the agent types whose `geometry.displayType` is `SPHERE`, sorted by id.
 * These are agents that are rendered as spheres and can be associated with an
 * external structure file matched by file base name against the agent type name.
 */
export function getSimulariumSphereAgentTypes(file: SimulariumFile): { id: number, name: string }[] {
    const out: { id: number, name: string }[] = [];
    const typeMapping = file.trajectoryInfo.typeMapping ?? {};
    for (const key of Object.keys(typeMapping)) {
        const id = Number(key);
        if (!Number.isFinite(id)) continue;
        const entry = typeMapping[key];
        const displayType = (entry?.geometry?.displayType ?? '').toUpperCase();
        if (displayType !== 'SPHERE') continue;
        out.push({ id, name: entry?.name ?? `type ${id}` });
    }
    out.sort((a, b) => a.id - b.id);
    return out;
}

/**
 * Resolve the spatial scale factor (to angstrom) for a Simularium file, mirroring the
 * conversion applied to particle coordinates and radii. A `scale <= 0` (or omitted)
 * auto-derives the factor from `trajectoryInfo.spatialUnits`.
 */
export function getSimulariumSpatialScale(file: SimulariumFile, scale?: number): number {
    return resolveScale(file, scale);
}

//

export { SimulariumFormat };

type SimulariumFormat = ModelFormat<SimulariumFile>

namespace SimulariumFormat {
    export function is(x?: ModelFormat): x is SimulariumFormat {
        return x?.kind === 'simularium';
    }

    export function create(simularium: SimulariumFile): SimulariumFormat {
        return { kind: 'simularium', name: 'simularium', data: simularium };
    }
}
