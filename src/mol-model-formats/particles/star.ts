/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../mol-data/db';
import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { ParticleList } from '../../mol-model/particles/particle-list';
import { CustomProperties } from '../../mol-model/custom-property';
import { RelionStarFile } from '../../mol-io/reader/relion/star';
import { RelionStar_Database } from '../../mol-io/reader/relion/schema';
import { degToRad } from '../../mol-math/misc';
import { ModelFormat } from '../format';

type Particles = RelionStar_Database['particles'];
type Optics = RelionStar_Database['optics'];

type CoordinateKind = 'pixel' | 'angstrom';

type TripletColumns = {
    x: Column<number>
    y: Column<number>
    z?: Column<number>
    rowCount: number
    kind: CoordinateKind
};

type TripletSpec = {
    x: keyof Particles
    y: keyof Particles
    z: keyof Particles
    kind: CoordinateKind
};

type Angles = { rot: number, tilt: number, psi: number };

type AngleSpec = {
    rot: 'rlnAngleRot' | 'rlnTomoSubtomogramRot'
    tilt: 'rlnAngleTilt' | 'rlnTomoSubtomogramTilt'
    psi: 'rlnAnglePsi' | 'rlnTomoSubtomogramPsi'
};

const RelionCoordinateSpecs: TripletSpec[] = [
    { x: 'rlnCenteredCoordinateXAngst', y: 'rlnCenteredCoordinateYAngst', z: 'rlnCenteredCoordinateZAngst', kind: 'angstrom' },
    { x: 'rlnCoordinateX', y: 'rlnCoordinateY', z: 'rlnCoordinateZ', kind: 'pixel' },
];

const RelionOriginSpecs: TripletSpec[] = [
    { x: 'rlnOriginXAngst', y: 'rlnOriginYAngst', z: 'rlnOriginZAngst', kind: 'angstrom' },
    { x: 'rlnOriginX', y: 'rlnOriginY', z: 'rlnOriginZ', kind: 'pixel' },
];

const ParticleAngleSpec: AngleSpec = {
    rot: 'rlnAngleRot', tilt: 'rlnAngleTilt', psi: 'rlnAnglePsi',
};

const SubtomogramAngleSpec: AngleSpec = {
    rot: 'rlnTomoSubtomogramRot', tilt: 'rlnTomoSubtomogramTilt', psi: 'rlnTomoSubtomogramPsi',
};

const RelionPixelSizeFields: ReadonlyArray<keyof Optics> = [
    'rlnTomoTiltSeriesPixelSize',
    'rlnImagePixelSize',
    'rlnMicrographPixelSize',
    'rlnDetectorPixelSize',
];

function getTripletColumns(particles: Particles, spec: TripletSpec): TripletColumns | undefined {
    const x = particles[spec.x] as Column<number>;
    const y = particles[spec.y] as Column<number>;
    const z = particles[spec.z] as Column<number>;
    if (!x.isDefined || !y.isDefined) return;
    return {
        x, y,
        z: z.isDefined ? z : void 0,
        rowCount: x.rowCount,
        kind: spec.kind,
    };
}

function getTripletColumnsFromSpecs(particles: Particles, specs: readonly TripletSpec[]) {
    for (const spec of specs) {
        const cols = getTripletColumns(particles, spec);
        if (cols) return cols;
    }
}

function readTriplet(out: Vec3, cols: TripletColumns, row: number) {
    return Vec3.set(out, cols.x.value(row), cols.y.value(row), cols.z ? cols.z.value(row) : 0);
}

function readAngles(out: Angles, particles: Particles, row: number, spec: AngleSpec) {
    out.rot = particles[spec.rot].value(row);
    out.tilt = particles[spec.tilt].value(row);
    out.psi = particles[spec.psi].value(row);
    return out;
}

function hasAngleColumns(particles: Particles, spec: AngleSpec) {
    return particles[spec.rot].isDefined
        && particles[spec.tilt].isDefined
        && particles[spec.psi].isDefined;
}

function relionEulerToRotation(out: Mat4, rot: number, tilt: number, psi: number) {
    // RELION's `Euler_angles2matrix(rot, tilt, psi)` (a.k.a. pyEM `euler2matrix`)
    // is the rotation that places the reference at the particle's orientation
    // in the tomogram. As a composition: A = (Rz(rot) Ry(tilt) Rz(psi))^T,
    // i.e. the transpose of the standard intrinsic ZYZ active matrix.
    const rotZ = Mat4.fromRotation(Mat4(), degToRad(rot), Vec3.unitZ);
    const tiltY = Mat4.fromRotation(Mat4(), degToRad(tilt), Vec3.unitY);
    const psiZ = Mat4.fromRotation(Mat4(), degToRad(psi), Vec3.unitZ);

    Mat4.mul(out, tiltY, psiZ);
    Mat4.mul(out, rotZ, out);
    Mat4.transpose(out, out);
    return out;
}

export interface RelionParticleListOptions {
    readonly tomograms?: ReadonlyArray<string>
    readonly micrographs?: ReadonlyArray<string>
    /** Override pixel size (Å/pixel) used to convert pixel-space coordinates to angstrom. */
    readonly pixelSize?: number
    /** Uniform particle radius in angstrom assigned to every particle. Leave 0 or undefined to omit radii. */
    readonly particleRadius?: number
}

function buildRelionLabel(particleBlockHeader: string, tomograms?: ReadonlyArray<string>, micrographs?: ReadonlyArray<string>) {
    const base = particleBlockHeader || 'RELION';
    const filters: string[] = [];
    if (tomograms !== void 0 && tomograms.length > 0) filters.push(tomograms.join(', '));
    if (micrographs !== void 0 && micrographs.length > 0) filters.push(micrographs.join(', '));
    if (filters.length > 0) {
        return `${base} particles (${filters.join('; ')})`;
    }
    return `${base} particles`;
}

/** Find the first defined pixel-size column in the optics block. */
function getOpticsPixelSizeColumn(optics: Optics | undefined): Column<number> | undefined {
    if (!optics) return;
    for (const name of RelionPixelSizeFields) {
        const col = optics[name] as Column<number>;
        if (col.isDefined && col.rowCount > 0) return col;
    }
}

type ResolvedPixelSize = {
    fixed: number | undefined
    perRow: Column<number> | undefined
    perGroup: Map<number, number> | undefined
};

function resolvePixelSize(particles: Particles, optics: Optics | undefined, options: RelionParticleListOptions): ResolvedPixelSize {
    const overrideValid = options.pixelSize !== void 0 && Number.isFinite(options.pixelSize) && options.pixelSize > 0;
    if (overrideValid) return { fixed: options.pixelSize, perRow: undefined, perGroup: undefined };

    // Per-row pixel size in the particles block takes precedence over the optics block.
    const particlePixelSizeCol = particles.rlnPixelSize;
    if (particlePixelSizeCol.isDefined && particlePixelSizeCol.rowCount > 0) {
        return { fixed: undefined, perRow: particlePixelSizeCol, perGroup: undefined };
    }

    const opticsCol = getOpticsPixelSizeColumn(optics);
    if (!opticsCol || !optics) return { fixed: undefined, perRow: undefined, perGroup: undefined };

    const opticsGroupCol = optics.rlnOpticsGroup;
    const particleGroupCol = particles.rlnOpticsGroup;
    if (opticsGroupCol.isDefined && particleGroupCol.isDefined) {
        const perGroup = new Map<number, number>();
        for (let i = 0, n = opticsGroupCol.rowCount; i < n; ++i) {
            perGroup.set(opticsGroupCol.value(i), opticsCol.value(i));
        }
        return { fixed: undefined, perRow: undefined, perGroup };
    }

    return { fixed: opticsCol.value(0), perRow: undefined, perGroup: undefined };
}

function buildStarAttrMaps(count: number, rowCount: number, defs: Array<[string, string, Float32Array | undefined]>) {
    const attributes = new Map<string, Float32Array>();
    const attributeInfo = new Map<string, { label: string, min: number, max: number }>();
    for (const [key, label, buf] of defs) {
        if (!buf) continue;
        const arr = count === rowCount ? buf : buf.slice(0, count);
        let min = Infinity, max = -Infinity;
        for (let i = 0; i < count; i++) {
            const v = arr[i];
            if (isFinite(v)) {
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        if (!isFinite(min)) continue;
        attributes.set(key, arr);
        attributeInfo.set(key, { label, min, max });
    }
    return attributes.size > 0
        ? { attributes, attributeInfo }
        : { attributes: undefined, attributeInfo: undefined };
}

export function createParticleListFromRelionStar(data: RelionStarFile, options: RelionParticleListOptions = {}): ParticleList {
    const { particleBlock, particles, optics } = data;
    const coordinateColumns = getTripletColumnsFromSpecs(particles, RelionCoordinateSpecs);
    if (!coordinateColumns) throw new Error(`Block '${particleBlock.header}' does not define supported particle coordinates.`);
    const originColumns = getTripletColumnsFromSpecs(particles, RelionOriginSpecs);

    const rowCount = coordinateColumns.rowCount;
    const tomoNameCol = particles.rlnTomoName;
    const micrographNameCol = particles.rlnMicrographName;
    const tomoFilter = options.tomograms !== void 0 && options.tomograms.length > 0
        ? new Set<string>(options.tomograms)
        : void 0;
    const micrographFilter = options.micrographs !== void 0 && options.micrographs.length > 0
        ? new Set<string>(options.micrographs)
        : void 0;

    const hasParticleAngles = hasAngleColumns(particles, ParticleAngleSpec);
    const hasSubtomogramAngles = hasAngleColumns(particles, SubtomogramAngleSpec);
    const hasRotations = hasParticleAngles || hasSubtomogramAngles;

    const { fixed: fixedPixelSize, perRow: perRowPixelSize, perGroup: perGroupPixelSize } = resolvePixelSize(particles, optics, options);
    const particleGroupCol = particles.rlnOpticsGroup;

    const _keys = new Int32Array(rowCount);
    const _coordinates = new Float32Array(rowCount * 3);
    const _rotations = hasRotations ? new Float32Array(rowCount * 4) : undefined;

    const _classCol = particles.rlnClassNumber;
    const _scoreCol = particles.rlnMaxValueProbDistribution;
    const _logCol = particles.rlnLogLikeliContribution;
    const _normCol = particles.rlnNormCorrection;
    const _classAttr = _classCol.isDefined ? new Float32Array(rowCount) : undefined;
    const _scoreAttr = _scoreCol.isDefined ? new Float32Array(rowCount) : undefined;
    const _logAttr = _logCol.isDefined ? new Float32Array(rowCount) : undefined;
    const _normAttr = _normCol.isDefined ? new Float32Array(rowCount) : undefined;

    const position = Vec3();
    const originShift = Vec3();
    const originRotation = Mat4();
    const rotation = Mat4();
    const subtomogram = Mat4();
    const quaternion = Quat();
    const particleAngles: Angles = { rot: 0, tilt: 0, psi: 0 };
    const subtomogramAngles: Angles = { rot: 0, tilt: 0, psi: 0 };

    let count = 0;
    for (let row = 0; row < rowCount; ++row) {
        if (tomoFilter !== void 0) {
            if (!tomoNameCol.isDefined) break;
            if (!tomoFilter.has(tomoNameCol.value(row))) continue;
        }
        if (micrographFilter !== void 0) {
            if (!micrographNameCol.isDefined) break;
            if (!micrographFilter.has(micrographNameCol.value(row))) continue;
        }

        const pixelScale = fixedPixelSize
            ?? perRowPixelSize?.value(row)
            ?? (perGroupPixelSize && particleGroupCol.isDefined
                ? perGroupPixelSize.get(particleGroupCol.value(row))
                : undefined)
            ?? 1;

        readTriplet(position, coordinateColumns, row);
        if (coordinateColumns.kind === 'pixel') Vec3.scale(position, position, pixelScale);

        if (hasSubtomogramAngles) {
            readAngles(subtomogramAngles, particles, row, SubtomogramAngleSpec);
        }

        if (originColumns) {
            readTriplet(originShift, originColumns, row);
            if (originColumns.kind === 'pixel') Vec3.scale(originShift, originShift, pixelScale);
            if (hasSubtomogramAngles) {
                relionEulerToRotation(originRotation, subtomogramAngles.rot, subtomogramAngles.tilt, subtomogramAngles.psi);
                Vec3.transformMat4(originShift, originShift, originRotation);
            }
            Vec3.sub(position, position, originShift);
        }

        _keys[count] = row;

        if (_classAttr) _classAttr[count] = _classCol.value(row);
        if (_scoreAttr) _scoreAttr[count] = _scoreCol.value(row);
        if (_logAttr) _logAttr[count] = _logCol.value(row);
        if (_normAttr) _normAttr[count] = _normCol.value(row);

        const cOffset = count * 3;
        _coordinates[cOffset + 0] = position[0];
        _coordinates[cOffset + 1] = position[1];
        _coordinates[cOffset + 2] = position[2];

        if (_rotations) {
            if (hasParticleAngles) {
                readAngles(particleAngles, particles, row, ParticleAngleSpec);
                relionEulerToRotation(rotation, particleAngles.rot, particleAngles.tilt, particleAngles.psi);
            } else {
                Mat4.setIdentity(rotation);
            }
            if (hasSubtomogramAngles) {
                relionEulerToRotation(subtomogram, subtomogramAngles.rot, subtomogramAngles.tilt, subtomogramAngles.psi);
                Mat4.mul(rotation, subtomogram, rotation);
            }
            Quat.normalize(quaternion, Quat.fromMat4(quaternion, rotation));
            const qOffset = count * 4;
            _rotations[qOffset + 0] = quaternion[0];
            _rotations[qOffset + 1] = quaternion[1];
            _rotations[qOffset + 2] = quaternion[2];
            _rotations[qOffset + 3] = quaternion[3];
        }

        ++count;
    }

    if (count === 0) {
        const filterDesc: string[] = [];
        if (tomoFilter !== void 0) filterDesc.push(`tomograms '${options.tomograms!.join(', ')}'`);
        if (micrographFilter !== void 0) filterDesc.push(`micrographs '${options.micrographs!.join(', ')}'`);
        throw new Error(filterDesc.length > 0
            ? `No RELION particle rows matched ${filterDesc.join(' and ')}.`
            : `Block '${particleBlock.header}' does not contain any readable particle rows.`);
    }

    const keys = count === rowCount ? _keys : _keys.slice(0, count);
    const coordinates = count === rowCount ? _coordinates : _coordinates.slice(0, count * 3);
    const rotations = _rotations && (count === rowCount ? _rotations : _rotations.slice(0, count * 4));

    const { attributes, attributeInfo } = buildStarAttrMaps(count, rowCount, [
        ['class', 'Class', _classAttr],
        ['score', 'Score', _scoreAttr],
        ['logLikelihood', 'Log-Likelihood', _logAttr],
        ['normCorrection', 'Norm Correction', _normAttr],
    ]);

    const radii = options.particleRadius && options.particleRadius > 0
        ? new Float32Array(count).fill(options.particleRadius)
        : undefined;

    return {
        label: buildRelionLabel(particleBlock.header, options.tomograms, options.micrographs),
        count,
        keys,
        targets: new Int32Array(count),
        coordinates,
        rotations,
        radii,
        attributes,
        attributeInfo,
        getParticleLabel: (index: number) => {
            const row = keys[index];
            const tomoName = tomoNameCol.isDefined ? tomoNameCol.value(row) : undefined;
            const micrographName = micrographNameCol.isDefined ? micrographNameCol.value(row) : undefined;
            const label: string[] = [`#${row + 1}`];
            if (tomoName) label.push(tomoName);
            if (micrographName) label.push(micrographName);
            return label.join(' | ');
        },
        sourceData: RelionStarFormat.create(data),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

//

export { RelionStarFormat };

type RelionStarFormat = ModelFormat<RelionStarFile>

namespace RelionStarFormat {
    export function is(x?: ModelFormat): x is RelionStarFormat {
        return x?.kind === 'relion-star';
    }

    export function create(relionStar: RelionStarFile): RelionStarFormat {
        return { kind: 'relion-star', name: 'relion-star', data: relionStar };
    }
}
