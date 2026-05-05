/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Column } from '../../../mol-data/db';
import { CifAliases } from '../cif/data-model';

import Schema = Column.Schema;

const str = Schema.str;
const float = Schema.float;
const int = Schema.int;

// from https://github.com/3dem/relion/blob/master/src/metadata_label.h

export const RelionStar_Schema = {
    particles: {
        /** X-coordinate (in Angstrom) for the origin of rotation */
        rlnOriginXAngst: float,
        /** Y-coordinate (in Angstrom) for the origin of rotation */
        rlnOriginYAngst: float,
        /** Z-coordinate (in Angstrom) for the origin of rotation */
        rlnOriginZAngst: float,
        /** X-coordinate (in pixels) for the origin of rotation */
        rlnOriginX: float,
        /** Y-coordinate (in pixels) for the origin of rotation */
        rlnOriginY: float,
        /** Z-coordinate (in pixels) for the origin of rotation */
        rlnOriginZ: float,

        /** First Euler angle (rot, in degrees) */
        rlnAngleRot: float,
        /** Second Euler angle (tilt, in degrees) */
        rlnAngleTilt: float,
        /** Third Euler, or in-plane angle (psi, in degrees) */
        rlnAnglePsi: float,

        /** ID (i.e. a unique number) for a micrograph */
        rlnMicrographId: int,
        /** Name of the micrograph */
        rlnMicrographName: str,

        /** First Euler angle of a subtomogram (rot, in degrees) */
        rlnTomoSubtomogramRot: float,
        /** Second Euler angle of a subtomogram (tilt, in degrees) */
        rlnTomoSubtomogramTilt: float,
        /** Third Euler angle of a subtomogram (psi, in degrees) */
        rlnTomoSubtomogramPsi: float,

        /** X-Position of an image in a micrograph (in pixels) */
        rlnCoordinateX: float,
        /** Y-Position of an image in a micrograph (in pixels) */
        rlnCoordinateY: float,
        /** Z-Position of an image in a 3D micrograph, i.e. tomogram (in pixels) */
        rlnCoordinateZ: float,
        /** X-Position of an image in a micrograph (in Angstroms, with the center being 0,0) */
        rlnCenteredCoordinateXAngst: float,
        /** Y-Position of an image in a micrograph (in Angstroms, with the center being 0,0) */
        rlnCenteredCoordinateYAngst: float,
        /** Z-Position of an image in a 3D micrograph, i.e. tomogram (in Angstroms, with the center being 0,0,0) */
        rlnCenteredCoordinateZAngst: float,

        /** Class number for which a particle has its highest probability */
        rlnClassNumber: int,
        /** ID (i.e. a unique number) for a particle */
        rlnParticleId: int,
        /** Name for a particle */
        rlnParticleName: str,
        /** Original name for a particle */
        rlnOriginalParticleName: str,
        /** Group of particles with identical optical properties */
        rlnOpticsGroup: int,
        /** Tomogram name */
        rlnTomoName: str,
        /** Size of the pixels in the references and images (in Angstroms) */
        rlnPixelSize: float,

        /** The number of a group of images */
        rlnGroupNumber: int,

        // pixel-size labels can also appear inline in particle blocks
        /** Pixel size of the original tilt series */
        rlnTomoTiltSeriesPixelSize: float,
        /** Pixel size (in Angstrom) */
        rlnImagePixelSize: float,
        /** Pixel size of (averaged) micrographs after binning in Angstrom/pixel */
        rlnMicrographPixelSize: float,
        /** Pixel size of the detector (in micrometers) */
        rlnDetectorPixelSize: float,
    },
    optics: {
        /** Group of particles with identical optical properties */
        rlnOpticsGroup: int,
        /** Pixel size of the original tilt series */
        rlnTomoTiltSeriesPixelSize: float,
        /** Pixel size (in Angstrom) */
        rlnImagePixelSize: float,
        /** Pixel size of (averaged) micrographs after binning in Angstrom/pixel */
        rlnMicrographPixelSize: float,
        /** Pixel size of the detector (in micrometers) */
        rlnDetectorPixelSize: float,
    },
};

export const RelionStar_Aliases: CifAliases = {
    'particles.rlnOriginXAngst': ['rlnOriginXAngst', 'rlnOriginXAngstrom'],
    'particles.rlnOriginYAngst': ['rlnOriginYAngst', 'rlnOriginYAngstrom'],
    'particles.rlnOriginZAngst': ['rlnOriginZAngst', 'rlnOriginZAngstrom'],
    'particles.rlnOriginX': ['rlnOriginX'],
    'particles.rlnOriginY': ['rlnOriginY'],
    'particles.rlnOriginZ': ['rlnOriginZ'],

    'particles.rlnAngleRot': ['rlnAngleRot', 'wrpAngleRot1'],
    'particles.rlnAngleTilt': ['rlnAngleTilt', 'wrpAngleTilt1'],
    'particles.rlnAnglePsi': ['rlnAnglePsi', 'wrpAnglePsi1'],

    'particles.rlnMicrographId': ['rlnMicrographId'],
    'particles.rlnMicrographName': ['rlnMicrographName'],

    'particles.rlnTomoSubtomogramRot': ['rlnTomoSubtomogramRot'],
    'particles.rlnTomoSubtomogramTilt': ['rlnTomoSubtomogramTilt'],
    'particles.rlnTomoSubtomogramPsi': ['rlnTomoSubtomogramPsi'],

    'particles.rlnCoordinateX': ['rlnCoordinateX', 'wrpCoordinateX1'],
    'particles.rlnCoordinateY': ['rlnCoordinateY', 'wrpCoordinateY1'],
    'particles.rlnCoordinateZ': ['rlnCoordinateZ', 'wrpCoordinateZ1'],
    'particles.rlnCenteredCoordinateXAngst': ['rlnCenteredCoordinateXAngst', 'rlnCenteredCoordinateXAngstrom'],
    'particles.rlnCenteredCoordinateYAngst': ['rlnCenteredCoordinateYAngst', 'rlnCenteredCoordinateYAngstrom'],
    'particles.rlnCenteredCoordinateZAngst': ['rlnCenteredCoordinateZAngst', 'rlnCenteredCoordinateZAngstrom'],

    'particles.rlnClassNumber': ['rlnClassNumber'],
    'particles.rlnParticleId': ['rlnParticleId'],
    'particles.rlnParticleName': ['rlnParticleName'],
    'particles.rlnOriginalParticleName': ['rlnOriginalParticleName'],
    'particles.rlnOpticsGroup': ['rlnOpticsGroup'],
    'particles.rlnTomoName': ['rlnTomoName'],
    'particles.rlnPixelSize': ['rlnPixelSize'],

    'particles.rlnGroupNumber': ['rlnGroupNumber'],

    'particles.rlnTomoTiltSeriesPixelSize': ['rlnTomoTiltSeriesPixelSize'],
    'particles.rlnImagePixelSize': ['rlnImagePixelSize'],
    'particles.rlnMicrographPixelSize': ['rlnMicrographPixelSize'],
    'particles.rlnDetectorPixelSize': ['rlnDetectorPixelSize'],

    'optics.rlnOpticsGroup': ['rlnOpticsGroup'],
    'optics.rlnTomoTiltSeriesPixelSize': ['rlnTomoTiltSeriesPixelSize'],
    'optics.rlnImagePixelSize': ['rlnImagePixelSize'],
    'optics.rlnMicrographPixelSize': ['rlnMicrographPixelSize'],
    'optics.rlnDetectorPixelSize': ['rlnDetectorPixelSize'],
};

export type RelionStar_Schema = typeof RelionStar_Schema;
export interface RelionStar_Database extends Database<RelionStar_Schema> {};