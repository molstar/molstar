/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../mol-task';
import { DcdFile } from '../../mol-io/reader/dcd/parser';
import { Coordinates, Frame, Time } from '../../mol-model/structure/coordinates';

const charmmTimeUnitFactor = 20.45482949774598;

export function coordinatesFromDcd(dcdFile: DcdFile): Task<Coordinates> {
    return Task.create('Parse DCD', async ctx => {
        await ctx.update('Converting to coordinates');

        const { header } = dcdFile;

        let deltaTime = Time(1, 'step');
        if (header.DELTA) {
            deltaTime = Time(header.DELTA * charmmTimeUnitFactor, 'ps');
        }

        let offsetTime = Time(0, deltaTime.unit);
        if (header.ISTART >= 1) {
            offsetTime = Time((header.ISTART - 1) * deltaTime.value, deltaTime.unit);
        }

        const frames: Frame[] = [];
        for (let i = 0, il = dcdFile.frames.length; i < il; ++i) {
            frames.push({
                ...dcdFile.frames[i],
                time: Time(offsetTime.value + deltaTime.value * i, deltaTime.unit)
            });
        }

        return Coordinates.create(frames, deltaTime, offsetTime);
    });
}
