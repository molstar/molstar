/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData } from '../../../mvs-data';
import { builderDemo, createMVSBuilder } from '../mvs-builder';


describe('mvs-builder', () => {
    it('mvs-builder demo works', () => {
        const mvsData = builderDemo();
        expect(typeof mvsData.metadata.version).toEqual('string');
        expect(typeof mvsData.metadata.timestamp).toEqual('string');
        expect(MVSData.validationIssues(mvsData)).toEqual(undefined);
    });

    it('volume builder works', () => {
        const builder = createMVSBuilder();
        builder
            .download({ url: 'https://www.ebi.ac.uk/pdbe/densities/x-ray/1tqn/box/-22.367,-33.367,-21.634/-7.106,-10.042,-0.937?detail=3' })
            .parse({ format: 'bcif' })
            .volume({ channel_id: '2FO-FC' })
            .representation({ type: 'isosurface', absolute_isovalue: 1, show_wireframe: true, show_faces: false });
        const state = builder.getState();
        expect(state).toBeDefined();
    });
});
