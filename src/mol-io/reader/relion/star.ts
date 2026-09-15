/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifBlock, CifFile } from '../cif';
import { toDatabase } from '../cif/schema';
import { ReaderResult as Result } from '../result';
import { Column } from '../../../mol-data/db';
import { RelionStar_Aliases, RelionStar_Database, RelionStar_Schema } from './schema';

const CoordinateFieldAliases = [
    'rlnCenteredCoordinateXAngst',
    'rlnCenteredCoordinateXAngstrom',
    'rlnCoordinateX',
    'wrpCoordinateX1',
];

export interface RelionStarFile {
    readonly source: CifFile
    readonly particleBlock: CifBlock
    readonly opticsBlock?: CifBlock
    readonly particles: RelionStar_Database['particles']
    readonly optics?: RelionStar_Database['optics']
}

function hasCoordinateFields(block: CifBlock) {
    return CoordinateFieldAliases.some(name => !!block.categories[name]);
}

function findParticlesBlock(file: CifFile) {
    let fallback: CifBlock | undefined;
    for (const block of file.blocks) {
        if (!hasCoordinateFields(block)) continue;
        if (!fallback) fallback = block;
        if (block.header.toLowerCase().includes('particle')) return block;
    }
    return fallback;
}

function findOpticsBlock(file: CifFile) {
    for (const block of file.blocks) {
        if (block.header.toLowerCase().includes('optics')) return block;
    }
}

export function parseRelionStar(file: CifFile) {
    const particleBlock = findParticlesBlock(file);
    if (!particleBlock) {
        return Result.error<RelionStarFile>('No RELION particle data block with coordinates was found.');
    }

    const opticsBlock = findOpticsBlock(file);
    const particles = toDatabase(RelionStar_Schema, particleBlock, RelionStar_Aliases).particles;
    const optics = opticsBlock
        ? toDatabase(RelionStar_Schema, opticsBlock, RelionStar_Aliases).optics
        : undefined;

    return Result.success<RelionStarFile>({
        source: file,
        particleBlock,
        opticsBlock,
        particles,
        optics,
    });
}

export function getRelionStarTomogramNames(file: CifFile): string[] {
    const result = parseRelionStar(file);
    if (result.isError) return [];
    const tomoName = result.result.particles.rlnTomoName;
    if (!tomoName.isDefined) return [];
    const tomograms = new Set<string>();
    for (let row = 0, n = tomoName.rowCount; row < n; ++row) {
        if (tomoName.valueKind(row) !== Column.ValueKinds.Present) continue;
        const v = tomoName.value(row);
        if (v) tomograms.add(v);
    }
    return Array.from(tomograms).sort();
}

export function getRelionStarMicrographNames(file: CifFile): string[] {
    const result = parseRelionStar(file);
    if (result.isError) return [];
    const micrographName = result.result.particles.rlnMicrographName;
    if (!micrographName.isDefined) return [];
    const micrographs = new Set<string>();
    for (let row = 0, n = micrographName.rowCount; row < n; ++row) {
        if (micrographName.valueKind(row) !== Column.ValueKinds.Present) continue;
        const v = micrographName.value(row);
        if (v) micrographs.add(v);
    }
    return Array.from(micrographs).sort();
}
