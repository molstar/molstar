/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { Camera } from '../../../mol-canvas3d/camera';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { RendererParams } from '../../renderer';
import { createLightingData, MaxWebGPULights } from '../lighting';

describe('native WebGPU lighting capacity', () => {
    it('accepts the portable 1024-light uniform budget and rejects only larger inputs', () => {
        const props = PD.getDefaultValues(RendererParams);
        const light = { inclination: 90, azimuth: 0, color: Color(0xffffff), intensity: 1 };
        props.light = Array.from({ length: MaxWebGPULights }, () => ({ ...light }));
        const camera = new Camera({}, { x: 0, y: 0, width: 1, height: 1 });
        const data = createLightingData(camera, props, 1, 1);
        expect(data.length).toBe(12 + MaxWebGPULights * 8);
        props.light = [...props.light, { ...light }];
        expect(() => createLightingData(camera, props, 1, 1)).toThrow(`at most ${MaxWebGPULights}`);
    });
});
