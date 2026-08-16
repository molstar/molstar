/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParticleList } from '../../../mol-model/particles/particle-list';
import { Asset } from '../../../mol-util/assets';
import { matchParticleTargetFiles } from '../particle-targets';

function particles(targets: number[], entities?: number[], names?: ReadonlyArray<[number, string]>): ParticleList {
    return {
        count: targets.length,
        targets: Int32Array.from(targets),
        entities: entities ? Int32Array.from(entities) : undefined,
        entityInfo: names ? new Map(names.map(([id, name]) => [id, { name }])) : undefined,
    } as ParticleList;
}

function file(name: string): Asset.File {
    return { kind: 'file', id: name as any, name };
}

describe('matchParticleTargetFiles', () => {
    it('matches exact entity names and strips regular and compressed extensions', () => {
        const data = particles([3, 4], [0, 1], [[0, 'ProteinA'], [1, 'Density']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb'), file('Density.map.gz')]);

        expect(result.matches.map(m => [m.file.name, m.targetIds])).toEqual([
            ['ProteinA.pdb', [3]],
            ['Density.map.gz', [4]],
        ]);
        expect(result.warnings).toEqual([]);
    });

    it('is case-sensitive and warns for unmatched files', () => {
        const data = particles([3], [0], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('proteina.pdb')]);

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual(["Could not match particle target file 'proteina.pdb' to an entity name."]);
    });

    it('loads one entity file for multiple target ids', () => {
        const data = particles([3, 4], [0, 0], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb')]);

        expect(result.matches[0].targetIds).toEqual([3, 4]);
    });

    it('skips target ids associated with multiple or missing entity names', () => {
        const data = particles([3, 3, 4], [0, 1, -1], [[0, 'A'], [1, 'B']]);
        const result = matchParticleTargetFiles(data, [file('A.pdb')]);

        expect(result.matches).toEqual([]);
        expect(result.warnings).toContain('Cannot match particle target 3 because it does not map to exactly one entity name.');
        expect(result.warnings).toContain('Cannot match particle target 4 because it has no entity name.');
    });

    it('skips duplicate basenames', () => {
        const data = particles([3], [0], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb'), file('ProteinA.cif')]);

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual(["Cannot match particle target files named 'ProteinA' because the filename is not unique."]);
    });

    it('skips explicitly excluded target ids', () => {
        const data = particles([3, 4], [0, 0], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb')], new Set([3]));

        expect(result.matches[0].targetIds).toEqual([4]);
    });

    it('silently ignores files shadowed entirely by explicit targets', () => {
        const data = particles([3], [0], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb')], new Set([3]));

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual([]);
    });

    it('warns when entity metadata is unavailable', () => {
        const result = matchParticleTargetFiles(particles([3]), [file('ProteinA.pdb')]);

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual(['Cannot match particle target files because the particle list has no entity metadata.']);
    });
});