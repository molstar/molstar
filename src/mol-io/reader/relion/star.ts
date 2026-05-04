/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author OpenAI
 */

import { CifBlock, CifField, CifFile } from '../cif';
import { Column } from '../../../mol-data/db';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { ParticleDistanceUnit, ParticleList, ParticleListParticle } from '../particle-list';

export type RelionDistanceUnit = ParticleDistanceUnit;

export interface RelionStarEulerAngles {
    rot: number
    tilt: number
    psi: number
}

export type RelionStarParticle = ParticleListParticle;
export type RelionStarParticleList = ParticleList;

type TripletFieldSpec = {
    x: string[]
    y: string[]
    z: string[]
    unit: RelionDistanceUnit
}

const CoordinateSpecs: TripletFieldSpec[] = [
    { x: ['rlnCenteredCoordinateXAngst', 'rlnCenteredCoordinateXAngstrom'], y: ['rlnCenteredCoordinateYAngst', 'rlnCenteredCoordinateYAngstrom'], z: ['rlnCenteredCoordinateZAngst', 'rlnCenteredCoordinateZAngstrom'], unit: 'angstrom' },
    { x: ['rlnCoordinateX'], y: ['rlnCoordinateY'], z: ['rlnCoordinateZ'], unit: 'pixel' },
];

const OriginSpecs: TripletFieldSpec[] = [
    { x: ['rlnOriginXAngst', 'rlnOriginXAngstrom'], y: ['rlnOriginYAngst', 'rlnOriginYAngstrom'], z: ['rlnOriginZAngst', 'rlnOriginZAngstrom'], unit: 'angstrom' },
    { x: ['rlnOriginX'], y: ['rlnOriginY'], z: ['rlnOriginZ'], unit: 'pixel' },
];

const ParticleAngleSpecs = [
    { rot: ['rlnAngleRot'], tilt: ['rlnAngleTilt'], psi: ['rlnAnglePsi'] },
];

const SubtomogramAngleSpecs = [
    { rot: ['rlnTomoSubtomogramRot'], tilt: ['rlnTomoSubtomogramTilt'], psi: ['rlnTomoSubtomogramPsi'] },
];

const PixelSizeFields = [
    'rlnTomoTiltSeriesPixelSize',
    'rlnImagePixelSize',
    'rlnMicrographPixelSize',
    'rlnDetectorPixelSize',
];

function getFlatField(block: CifBlock, name: string) {
    return block.categories[name]?.getField('');
}

function getFirstFlatField(block: CifBlock, names: readonly string[]) {
    for (const name of names) {
        const field = getFlatField(block, name);
        if (field) return field;
    }
}

function hasPresentValue(field: CifField | undefined, row: number) {
    return !!field && field.valueKind(row) === Column.ValueKinds.Present;
}

function getOptionalNumber(field: CifField | undefined, row: number) {
    return hasPresentValue(field, row) ? field!.float(row) : void 0;
}

function getOptionalString(field: CifField | undefined, row: number) {
    return hasPresentValue(field, row) ? field!.str(row) : void 0;
}

function getFieldRowCount(field: CifField | undefined, block: CifBlock, name: string) {
    return field ? block.categories[name].rowCount : 0;
}

function getTripletFields(block: CifBlock, spec: TripletFieldSpec) {
    const xName = spec.x.find(name => !!getFlatField(block, name));
    const yName = spec.y.find(name => !!getFlatField(block, name));
    const zName = spec.z.find(name => !!getFlatField(block, name));

    if (!xName || !yName) return;

    return {
        xName,
        yName,
        zName,
        x: getFlatField(block, xName)!,
        y: getFlatField(block, yName)!,
        z: zName ? getFlatField(block, zName) : void 0,
        rowCount: block.categories[xName].rowCount,
        unit: spec.unit,
    };
}

function getTripletFieldsFromSpecs(block: CifBlock, specs: readonly TripletFieldSpec[]) {
    for (const spec of specs) {
        const fields = getTripletFields(block, spec);
        if (fields) return fields;
    }
}

function getAngles(block: CifBlock, row: number, specs: ReadonlyArray<{ rot: readonly string[], tilt: readonly string[], psi: readonly string[] }>) {
    for (const spec of specs) {
        const rot = getOptionalNumber(getFirstFlatField(block, spec.rot), row);
        const tilt = getOptionalNumber(getFirstFlatField(block, spec.tilt), row);
        const psi = getOptionalNumber(getFirstFlatField(block, spec.psi), row);
        if (rot === void 0 && tilt === void 0 && psi === void 0) continue;
        return {
            rot: rot ?? 0,
            tilt: tilt ?? 0,
            psi: psi ?? 0,
        };
    }
}

function getTripletValue(block: CifBlock, row: number, specs: readonly TripletFieldSpec[]) {
    for (const spec of specs) {
        const fields = getTripletFields(block, spec);
        if (!fields) continue;

        const x = getOptionalNumber(fields.x, row);
        const y = getOptionalNumber(fields.y, row);
        const z = getOptionalNumber(fields.z, row);
        if (x === void 0 || y === void 0) continue;

        return {
            value: Vec3.create(x, y, z ?? 0),
            unit: fields.unit,
        };
    }
}

function findParticlesBlock(file: CifFile) {
    let fallback: CifBlock | undefined;
    for (const block of file.blocks) {
        if (getTripletFields(block, CoordinateSpecs[0]) || getTripletFields(block, CoordinateSpecs[1])) {
            if (!fallback) fallback = block;
            if (block.header.toLowerCase().includes('particle')) return block;
        }
    }
    return fallback;
}

function findOpticsBlock(file: CifFile) {
    for (const block of file.blocks) {
        if (block.header.toLowerCase().includes('optics')) return block;
    }
    for (const block of file.blocks) {
        if (PixelSizeFields.some(name => !!getFlatField(block, name))) return block;
    }
}

function getNumericFieldUniqueValue(block: CifBlock | undefined, names: readonly string[]) {
    if (!block) return;

    for (const name of names) {
        const field = getFlatField(block, name);
        if (!field) continue;

        const rowCount = getFieldRowCount(field, block, name);
        let value: number | undefined = void 0;
        let isUnique = true;

        for (let i = 0; i < rowCount; i++) {
            const current = getOptionalNumber(field, i);
            if (current === void 0) continue;
            if (value === void 0) {
                value = current;
            } else if (Math.abs(value - current) > 1e-6) {
                isUnique = false;
                break;
            }
        }

        if (isUnique && value !== void 0 && Number.isFinite(value) && value > 0) return value;
    }
}

function degToRad(value: number) {
    return value * Math.PI / 180;
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

export function parseRelionStarParticleList(file: CifFile): RelionStarParticleList {
    const particleBlock = findParticlesBlock(file);
    if (!particleBlock) throw new Error('No RELION particle data block with coordinates was found.');

    const coordinateFields = getTripletFieldsFromSpecs(particleBlock, CoordinateSpecs);
    if (!coordinateFields) throw new Error(`Block '${particleBlock.header}' does not define supported particle coordinates.`);

    const rowCount = coordinateFields.rowCount;

    const opticsBlock = findOpticsBlock(file);
    const opticsGroupField = getFlatField(particleBlock, 'rlnOpticsGroup');
    const tomogramField = getFlatField(particleBlock, 'rlnTomoName');
    const micrographField = getFlatField(particleBlock, 'rlnMicrographName');
    const tomogramParticleField = getFlatField(particleBlock, 'rlnTomoParticleName');
    const imageField = getFlatField(particleBlock, 'rlnImageName');
    const groupNumberField = getFlatField(particleBlock, 'rlnGroupNumber');
    const classNumberField = getFlatField(particleBlock, 'rlnClassNumber');
    const warnings: string[] = [];
    const particles: RelionStarParticle[] = [];

    for (let row = 0; row < rowCount; row++) {
        const coordinate = getTripletValue(particleBlock, row, CoordinateSpecs);
        if (!coordinate) continue;

        const origin = getTripletValue(particleBlock, row, OriginSpecs) ?? { value: Vec3.create(0, 0, 0), unit: 'pixel' as const };
        const particleAngles = getAngles(particleBlock, row, ParticleAngleSpecs) ?? { rot: 0, tilt: 0, psi: 0 };
        const subtomogramAngles = getAngles(particleBlock, row, SubtomogramAngleSpecs);
        const subtomogram = subtomogramAngles
            ? relionEulerToRotation(Mat4(), subtomogramAngles.rot, subtomogramAngles.tilt, subtomogramAngles.psi)
            : void 0;
        const rotation = relionEulerToRotation(Mat4(), particleAngles.rot, particleAngles.tilt, particleAngles.psi);
        if (subtomogram) {
            Mat4.mul(rotation, subtomogram, rotation);
        }
        const opticsGroup = getOptionalString(opticsGroupField, row);
        const tomogram = getOptionalString(tomogramField, row);
        const micrograph = getOptionalString(micrographField, row);
        const tomogramParticle = getOptionalString(tomogramParticleField, row);
        const image = getOptionalString(imageField, row);
        const groupNumber = getOptionalNumber(groupNumberField, row);
        const classNumber = getOptionalNumber(classNumberField, row);

        particles.push({
            index: row,
            coordinate: coordinate.value,
            coordinateUnit: coordinate.unit,
            origin: origin.value,
            originUnit: origin.unit,
            originRotation: subtomogram ? Mat4.copy(Mat4(), subtomogram) : void 0,
            rotation,
            metadata: {
                opticsGroup,
                groupNumber,
                classNumber,
                tomogram,
                tomoName: tomogram,
                micrograph,
                micrographName: micrograph,
                tomoParticleName: tomogramParticle,
                imageName: image,
                particleRot: particleAngles.rot,
                particleTilt: particleAngles.tilt,
                particlePsi: particleAngles.psi,
                subtomogramRot: subtomogramAngles?.rot,
                subtomogramTilt: subtomogramAngles?.tilt,
                subtomogramPsi: subtomogramAngles?.psi,
            }
        });
    }

    if (particles.length === 0) {
        throw new Error(`Block '${particleBlock.header}' does not contain any readable particle rows.`);
    }

    if (!getAngles(particleBlock, 0, ParticleAngleSpecs) && !getAngles(particleBlock, 0, SubtomogramAngleSpecs)) {
        warnings.push('No RELION rotation columns were found; particle rotations default to identity.');
    }

    return {
        format: 'relion-star',
        particleBlockHeader: particleBlock.header,
        opticsBlockHeader: opticsBlock?.header,
        particles,
        suggestedScale: coordinateFields.unit === 'pixel'
            ? getNumericFieldUniqueValue(opticsBlock, PixelSizeFields) ?? getNumericFieldUniqueValue(particleBlock, PixelSizeFields) ?? 1
            : 1,
        warnings,
    };
}
