/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import fs from 'fs';
import { MVSData } from '../mvs-data';


describe('MVSData', () => {
    it('MVSData functions work', async () => {
        const data = fs.readFileSync('examples/mvs/1cbs.mvsj', { encoding: 'utf8' });
        const mvsData = MVSData.fromMVSJ(data);
        expect(mvsData).toBeTruthy();

        expect(MVSData.validationIssues(mvsData)).toEqual(undefined);

        expect(MVSData.isValid(mvsData)).toEqual(true);

        const reencoded = MVSData.toMVSJ(mvsData);
        expect(reencoded.replace(/\s/g, '')).toEqual(data.replace(/\s/g, ''));

        const prettyString = MVSData.toPrettyString(mvsData);
        expect(typeof prettyString).toEqual('string');
        expect(prettyString.length).toBeGreaterThan(0);
    });

    it('MVSData builder works', async () => {
        const builder = MVSData.createBuilder();
        expect(builder).toBeTruthy();

        const mvsData = builder.getState();
        expect(MVSData.validationIssues(mvsData)).toEqual(undefined);

        builder
            .download({ url: 'http://example.com' })
            .parse({ format: 'mmcif' })
            .assemblyStructure({ assembly_id: '1' })
            .component({ selector: 'polymer' })
            .representation()
            .color({ color: 'green', selector: { label_asym_id: 'A' } });
        const mvsData2 = builder.getState();
        expect(MVSData.validationIssues(mvsData2)).toEqual(undefined);
    });
});
