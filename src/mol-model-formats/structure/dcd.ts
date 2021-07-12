/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../mol-task';
import { DcdFile } from '../../mol-io/reader/dcd/parser';
import { Coordinates, Frame, Time } from '../../mol-model/structure/coordinates';
import { Vec3 } from '../../mol-math/linear-algebra';
import { degToRad, halfPI } from '../../mol-math/misc';
import { Cell } from '../../mol-math/geometry/spacegroup/cell';
import { Mutable } from '../../mol-util/type-helpers';
import { EPSILON, equalEps } from '../../mol-math/linear-algebra/3d/common';

const charmmTimeUnitFactor = 20.45482949774598;

export function coordinatesFromDcd(dcdFile: DcdFile): Task<Coordinates> {
    return Task.create('Parse DCD', async ctx => {
        await ctx.update('Converting to coordinates');

        const { header } = dcdFile;

        const deltaTime = header.DELTA
            ? Time(header.DELTA * charmmTimeUnitFactor, 'ps')
            : Time(1, 'step');

        const offsetTime = header.ISTART >= 1
            ? Time((header.ISTART - 1) * deltaTime.value, deltaTime.unit)
            : Time(0, deltaTime.unit);

        const frames: Frame[] = [];
        for (let i = 0, il = dcdFile.frames.length; i < il; ++i) {
            const dcdFrame = dcdFile.frames[i];
            const frame: Mutable<Frame> = {
                elementCount: dcdFrame.elementCount,
                time: Time(offsetTime.value + deltaTime.value * i, deltaTime.unit),

                x: dcdFrame.x,
                y: dcdFrame.y,
                z: dcdFrame.z,

                xyzOrdering: { isIdentity: true }
            };

            if (dcdFrame.cell) {
                // this is not standardized, using heuristics to handle variants
                const c = dcdFrame.cell;
                if (c[1] >= -1 && c[1] <= 1 && c[3] >= -1 && c[3] <= 1 && c[4] >= -1 && c[4] <= 1) {
                    frame.cell = Cell.create(
                        Vec3.create(c[0], c[2], c[5]),
                        Vec3.create(
                            degToRad(90 - Math.asin(c[1]) * 90 / halfPI),
                            degToRad(90 - Math.asin(c[3]) * 90 / halfPI),
                            degToRad(90 - Math.asin(c[4]) * 90 / halfPI)
                        )
                    );
                } else if (
                    c[0] < 0 || c[1] < 0 || c[2] < 0 || c[3] < 0 || c[4] < 0 || c[5] < 0 ||
                    c[3] > 180 || c[4] > 180 || c[5] > 180
                ) {
                    frame.cell = Cell.fromBasis(
                        Vec3.create(c[0], c[1], c[3]),
                        Vec3.create(c[1], c[2], c[4]),
                        Vec3.create(c[3], c[4], c[5])
                    );
                } else {
                    frame.cell = Cell.create(
                        Vec3.create(c[0], c[2], c[5]),
                        // interpret angles very close to 0 as 90 deg
                        Vec3.create(
                            degToRad(equalEps(c[1], 0, EPSILON) ? 90 : c[1]),
                            degToRad(equalEps(c[3], 0, EPSILON) ? 90 : c[3]),
                            degToRad(equalEps(c[4], 0, EPSILON) ? 90 : c[4])
                        )
                    );
                }
            }
            frames.push(frame);
        }

        return Coordinates.create(frames, deltaTime, offsetTime);
    });
}
