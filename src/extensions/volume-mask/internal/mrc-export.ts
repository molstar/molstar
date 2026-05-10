/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Grid } from '../../../mol-model/volume';
import { CCP4Writer } from '../../../mol-io/writer/ccp4/ccp4';

export function downloadMrc(grid: Grid, maskData: Uint8Array | Float32Array, filename = 'mask.mrc') {
    const buf = CCP4Writer.writeMrc(grid, maskData);
    const blob = new Blob([buf], { type: 'application/octet-stream' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}
