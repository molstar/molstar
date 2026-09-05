/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomProperties } from '../../../mol-model/custom-property';
import { ParticleList } from '../../../mol-model/particles/particle-list';
import { Asset } from '../../../mol-util/assets';
import { matchParticleTargetFiles } from '../particle-targets';

/** One particle per given target id, with the given (optional) entity index and entity names. */
function particles(targetEntities: ReadonlyArray<[targetId: number, entity: number | undefined]>, names?: ReadonlyArray<[number, string]>): ParticleList {
    const count = targetEntities.length;
    return {
        count,
        keys: Int32Array.from({ length: count }, (_, i) => i),
        targets: Int32Array.from(targetEntities.map(([id]) => id)),
        targetInfo: new Map(targetEntities.map(([id, entity]) => [id, { entity }])),
        entityInfo: names ? new Map(names.map(([id, name]) => [id, { name }])) : undefined,
        coordinates: new Float32Array(count * 3),
        getParticleLabel: (index: number) => `#${index}`,
        sourceData: { kind: 'test', name: 'test', data: {} },
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}


function file(name: string): Asset.File {
    return { kind: 'file', id: name as any, name };
}

describe('matchParticleTargetFiles', () => {
    it('matches exact entity names and strips regular and compressed extensions', () => {
        const data = particles([[3, 0], [4, 1]], [[0, 'ProteinA'], [1, 'Density']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb'), file('Density.map.gz')]);

        expect(result.matches.map(m => [m.file.name, m.targetIds])).toEqual([
            ['ProteinA.pdb', [3]],
            ['Density.map.gz', [4]],
        ]);
        expect(result.warnings).toEqual([]);
    });

    it('is case-sensitive and warns for unmatched files', () => {
        const data = particles([[3, 0]], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('proteina.pdb')]);

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual(["Could not match particle target file 'proteina.pdb' to an entity name."]);
    });

    it('loads one entity file for multiple target ids', () => {
        const data = particles([[3, 0], [4, 0]], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb')]);

        expect(result.matches[0].targetIds).toEqual([3, 4]);
    });

    it('skips target ids with no entity name', () => {
        const data = particles([[3, 0], [4, undefined]], [[0, 'A']]);
        const result = matchParticleTargetFiles(data, [file('A.pdb')]);

        expect(result.matches[0].targetIds).toEqual([3]);
        expect(result.warnings).toContain('Cannot match particle target 4 because it has no entity name.');
    });

    it('skips duplicate basenames', () => {
        const data = particles([[3, 0]], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb'), file('ProteinA.cif')]);

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual(["Cannot match particle target files named 'ProteinA' because the filename is not unique."]);
    });

    it('skips explicitly excluded target ids', () => {
        const data = particles([[3, 0], [4, 0]], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb')], new Set([3]));

        expect(result.matches[0].targetIds).toEqual([4]);
    });

    it('silently ignores files shadowed entirely by explicit targets', () => {
        const data = particles([[3, 0]], [[0, 'ProteinA']]);
        const result = matchParticleTargetFiles(data, [file('ProteinA.pdb')], new Set([3]));

        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual([]);
    });

    it('warns when entity metadata is unavailable', () => {
        const result = matchParticleTargetFiles(particles([[3, undefined]]), [file('ProteinA.pdb')]);


        expect(result.matches).toEqual([]);
        expect(result.warnings).toEqual(['Cannot match particle target files because the particle list has no entity metadata.']);
    });
});