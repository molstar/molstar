/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../mol-task';
import { NctrajFile } from '../../mol-io/reader/nctraj/parser';
import { Coordinates, Frame, Time } from '../../mol-model/structure/coordinates';
import { Cell } from '../../mol-math/geometry/spacegroup/cell';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Mutable } from '../../mol-util/type-helpers';

export function coordinatesFromNctraj(file: NctrajFile): Task<Coordinates> {
    return Task.create('Parse NCTRAJ', async ctx => {
        await ctx.update('Converting to coordinates');

        const deltaTime = Time(file.deltaTime, 'step');
        const offsetTime = Time(file.timeOffset, deltaTime.unit);

        const frames: Frame[] = [];
        for (let i = 0, il = file.coordinates.length; i < il; ++i) {
            const c = file.coordinates[i];
            const elementCount = c.length / 3;
            const x = new Float32Array(elementCount);
            const y = new Float32Array(elementCount);
            const z = new Float32Array(elementCount);
            for (let j = 0, jl = c.length; j < jl; j += 3) {
                x[j / 3] = c[j];
                y[j / 3] = c[j + 1];
                z[j / 3] = c[j + 2];
            }
            const frame: Mutable<Frame> = {
                elementCount,
                x, y, z,
                xyzOrdering: { isIdentity: true },
                time: Time(offsetTime.value + deltaTime.value * i, deltaTime.unit)
            };
            // TODO: handle case where cell_lengths and cell_angles are set, i.e., angles not 90deg
            if (file.cell_lengths) {
                const lengths = file.cell_lengths[i];
                const x = Vec3.scale(Vec3(), Vec3.unitX, lengths[0]);
                const y = Vec3.scale(Vec3(), Vec3.unitY, lengths[1]);
                const z = Vec3.scale(Vec3(), Vec3.unitZ, lengths[2]);
                frame.cell = Cell.fromBasis(x, y, z);
            }
            frames.push(frame);
        }

        return Coordinates.create(frames, deltaTime, offsetTime);
    });
}
