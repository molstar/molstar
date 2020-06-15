/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../../mol-task';
import { XtcFile } from '../../mol-io/reader/xtc/parser';
import { Coordinates, Frame, Time } from '../../mol-model/structure/coordinates';
import { Cell } from '../../mol-math/geometry/spacegroup/cell';

export function coordinatesFromXtc(file: XtcFile): Task<Coordinates> {
    return Task.create('Parse XTC', async ctx => {
        await ctx.update('Converting to coordinates');

        const deltaTime = Time(file.deltaTime, 'step');
        const offsetTime = Time(file.timeOffset, deltaTime.unit);

        const frames: Frame[] = [];
        for (let i = 0, il = file.frames.length; i < il; ++i) {
            frames.push({
                elementCount: file.frames[i].count,
                // TODO:
                cell: Cell.empty(),
                x: file.frames[i].x,
                y: file.frames[i].y,
                z: file.frames[i].z,
                time: Time(offsetTime.value + deltaTime.value * i, deltaTime.unit)
            });
        }

        return Coordinates.create(frames, deltaTime, offsetTime);
    });
}
