/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { PositionLocation } from '../../../../mol-geo/util/location-iterator';
import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { Color } from '../../../../mol-util/color';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { BodyLabels } from '../labels';
import { BodyLabelColorTheme, BodyLabelColorThemeParams, BodyLabelColorThemeProvider, makeLabelAtPosition } from '../theme';
import { CanonicalOrder, SwappedOrder, createTestVolume, offsetOf } from './test-volume';

const props = PD.getDefaultValues(BodyLabelColorThemeParams);
const Red = Color(0xff0000), Blue = Color(0x0000ff);

describe('BodyLabelColorTheme', () => {
    it('falls back to a uniform theme without labels', () => {
        const volume = createTestVolume([2, 2, 2], () => 1);
        const theme = BodyLabelColorTheme({ volume }, props);
        expect(theme.granularity).toBe('uniform');
        expect(theme.color(PositionLocation(Vec3.create(0.5, 0.5, 0.5)), false)).toBe(props.unassignedColor);
        expect(BodyLabelColorThemeProvider.isApplicable({ volume })).toBe(true);
        expect(BodyLabelColorThemeProvider.isApplicable({})).toBe(false);
    });

    for (const order of [CanonicalOrder, SwappedOrder]) {
        it(`colors a vertex by the densest corner of its cell (axis order ${order})`, () => {
            // density grows with i; a 2-voxel step at i = 2 belongs to body 2
            const matrix = Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2));
            Mat4.setTranslation(matrix, Vec3.create(100, 0, 0));
            const volume = createTestVolume([4, 4, 4], i => i >= 2 ? 1 : 0, order, matrix);
            const store = BodyLabels.ensure(volume);
            store.bodies.push({ id: 2, name: 'Two', color: Red, views: [], remainder: false, voxelCount: 32 }, { id: 3, name: 'Three', color: Blue, views: [], remainder: false, voxelCount: 0 });
            for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++) {
                store.labels[offsetOf(volume, 2, j, k)] = 2;
                store.labels[offsetOf(volume, 3, j, k)] = 2;
            }
            store.version = 7;

            const theme = BodyLabelColorTheme({ volume }, props);
            expect(theme.granularity).toBe('vertex');
            expect(theme.contextHash).toBe(7);

            // vertex on the edge between i = 1 (outside) and i = 2 (inside), world x = 100 + 2 * 1.5
            const onEdge = PositionLocation(Vec3.create(103, 2 * 1.2, 2 * 2.7));
            expect(theme.color(onEdge, false)).toBe(Red);
            // deep in the unassigned low-density region
            expect(theme.color(PositionLocation(Vec3.create(100.4, 1, 1)), false)).toBe(props.unassignedColor);
            // outside the grid: clamped to the boundary cell
            expect(theme.color(PositionLocation(Vec3.create(150, 1, 1)), false)).toBe(Red);
            expect(theme.color(PositionLocation(Vec3.create(-50, 1, 1)), false)).toBe(props.unassignedColor);
        });
    }

    it('prefers the densest corner even when the vertex is nearer the other one', () => {
        const volume = createTestVolume([3, 1, 1], i => i === 1 ? 5 : 0);
        const store = BodyLabels.ensure(volume);
        store.labels[offsetOf(volume, 1, 0, 0)] = 1;
        const labelAt = makeLabelAtPosition(volume, store.labels);
        expect(labelAt(Vec3.create(0.1, 0, 0))).toBe(1);
        expect(labelAt(Vec3.create(1.9, 0, 0))).toBe(1);
    });
});
