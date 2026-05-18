/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author OpenAI
 */

import { Column } from '../../../mol-data/db';
import { parseCifText } from '../cif/text/parser';
import { parseRelionStar } from '../relion/star';

test('parses RELION STAR blocks and keeps particle and optics blocks', async () => {
    const data = `data_optics
loop_
_rlnOpticsGroup
_rlnTomoTiltSeriesPixelSize
1 4.5

data_particles
loop_
_rlnCoordinateX
_rlnCoordinateY
_rlnCoordinateZ
_rlnTomoName
10 20 30 tomo-a
40 50 60 tomo-b
`;

    const parsed = await parseCifText(data).run();
    if (parsed.isError) throw new Error(parsed.message);

    const relion = parseRelionStar(parsed.result);
    if (relion.isError) throw new Error(relion.message);

    expect(relion.result.particleBlock.header).toBe('particles');
    expect(relion.result.opticsBlock?.header).toBe('optics');
    expect(relion.result.source.blocks).toHaveLength(2);
});

test('exposes typed RELION particle and optics tables', async () => {
    const data = `data_optics
loop_
_rlnOpticsGroup
_rlnImagePixelSize
1 1.25

data_particles
loop_
_rlnCoordinateX
_rlnCoordinateY
_rlnCoordinateZ
_rlnAngleRot
_rlnAngleTilt
_rlnAnglePsi
_rlnOpticsGroup
_rlnTomoName
10 20 30 11 22 33 1 tomo-a
40 50 60 44 55 66 1 tomo-b
`;

    const parsed = await parseCifText(data).run();
    if (parsed.isError) throw new Error(parsed.message);

    const relion = parseRelionStar(parsed.result);
    if (relion.isError) throw new Error(relion.message);

    const { particles, optics } = relion.result;

    expect(particles.rlnCoordinateX.isDefined).toBe(true);
    expect(particles.rlnCoordinateX.rowCount).toBe(2);
    expect(particles.rlnCoordinateX.value(0)).toBeCloseTo(10);
    expect(particles.rlnCoordinateY.value(1)).toBeCloseTo(50);
    expect(particles.rlnCoordinateZ.value(1)).toBeCloseTo(60);

    expect(particles.rlnAngleRot.isDefined).toBe(true);
    expect(particles.rlnAngleRot.value(1)).toBeCloseTo(44);

    expect(particles.rlnOpticsGroup.value(0)).toBe(1);
    expect(particles.rlnTomoName.value(0)).toBe('tomo-a');

    expect(particles.rlnCenteredCoordinateXAngst.isDefined).toBe(false);

    expect(optics).toBeDefined();
    expect(optics!.rlnOpticsGroup.isDefined).toBe(true);
    expect(optics!.rlnImagePixelSize.value(0)).toBeCloseTo(1.25);
});

test('aliases multi-variant centered coordinate field names', async () => {
    const data = `data_particles
loop_
_rlnCenteredCoordinateXAngstrom
_rlnCenteredCoordinateYAngstrom
_rlnCenteredCoordinateZAngstrom
1.5 2.5 3.5
`;

    const parsed = await parseCifText(data).run();
    if (parsed.isError) throw new Error(parsed.message);

    const relion = parseRelionStar(parsed.result);
    if (relion.isError) throw new Error(relion.message);

    const { particles } = relion.result;
    expect(particles.rlnCenteredCoordinateXAngst.isDefined).toBe(true);
    expect(particles.rlnCenteredCoordinateXAngst.value(0)).toBeCloseTo(1.5);
    expect(particles.rlnCenteredCoordinateYAngst.value(0)).toBeCloseTo(2.5);
    expect(particles.rlnCenteredCoordinateZAngst.value(0)).toBeCloseTo(3.5);
    // pixel-coordinate alias should not be present
    expect(particles.rlnCoordinateX.valueKind(0)).toBe(Column.ValueKinds.NotPresent);
});

