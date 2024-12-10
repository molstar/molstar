/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { MVSData } from '../../../mvs-data';
import { builderDemo } from '../mvs-builder';


describe('mvs-builder', () => {
    it('mvs-builder demo works', async () => {
        const mvsData = builderDemo();
        expect(typeof mvsData.metadata.version).toEqual('string');
        expect(typeof mvsData.metadata.timestamp).toEqual('string');
        expect(MVSData.validationIssues(mvsData)).toEqual(undefined);
    });
});
