/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parse3DG } from '../3dg/parser';

const basic3dgString = `1(mat)	1420000	0.791377837067	10.9947291355	-13.1882897693
1(mat)	1440000	-0.268241283699	10.5200875887	-13.0896257278
1(mat)	1460000	-1.3853075236	10.5513787498	-13.1440142173
1(mat)	1480000	-1.55984101733	11.4340829129	-13.6026301209
1(mat)	1500000	-0.770991778399	11.4758488546	-14.5881137222
1(mat)	1520000	-0.0848245107875	12.2624690808	-14.354289628
1(mat)	1540000	-0.458643807046	12.5985791771	-13.4701149287
1(mat)	1560000	-0.810322906201	12.2461643989	-12.3172933413
1(mat)	1580000	-2.08211172035	12.8886838656	-12.8742007778
1(mat)	1600000	-3.52093948201	13.1850935438	-12.4118684428`;

describe('3dg reader', () => {
    it('basic', async () => {
        const parsed = await parse3DG(basic3dgString).run();
        expect(parsed.isError).toBe(false);

        if (parsed.isError) return;
        const { chromosome, position, x, y, z } = parsed.result.table;
        expect(chromosome.value(0)).toBe('1(mat)');
        expect(position.value(1)).toBe(1440000);
        expect(x.value(5)).toBe(-0.0848245107875);
        expect(y.value(5)).toBe(12.2624690808);
        expect(z.value(5)).toBe(-14.354289628);
    });
});