/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { treeValidationIssues } from '../../generic/tree-schema';
import { builderDemo } from '../mvs-builder';
import { MVSTreeSchema } from '../mvs-tree';


describe('mvs-builder', () => {
    it('mvs-builder demo works', async () => {
        const mvsData = builderDemo();
        expect(typeof mvsData.metadata.version).toEqual('string');
        expect(typeof mvsData.metadata.timestamp).toEqual('string');
        expect(treeValidationIssues(MVSTreeSchema, mvsData.root)).toEqual(undefined);
    });
});
