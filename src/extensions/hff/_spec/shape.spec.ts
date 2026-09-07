/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Asserts the SFF shape params: the mesh defaults that make thin segmentation
 * surfaces readable, and the procedural-animation controls inherited from
 * `Mesh.Params` so HFF surfaces can be animated from the standard panel.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { SffShapeParams } from '../shape';

describe('hff shape params', () => {
    it('renders double-sided and keeps the back face in the segment colour', () => {
        const defaults = PD.getDefaultValues(SffShapeParams);
        // SFF segments are thin oriented surfaces; the mol* default is false
        expect(defaults.doubleSided).toBe(true);
        // 0 = do not tint the interior grey, which is the mol* default of 1
        expect(defaults.interior.colorStrength).toBe(0);
    });

    it('exposes the animation PD.Group with wiggle and tumble controls', () => {
        const animParam = (SffShapeParams as any).animation;
        expect(animParam).toBeTruthy();
        expect(animParam.type).toBe('group');

        const params = animParam.params;
        for (const k of [
            'wiggleMode',
            'wiggleSpeed',
            'wiggleAmplitude',
            'wiggleFrequency',
            'tumbleSpeed',
            'tumbleAmplitude',
            'tumbleFrequency',
        ]) {
            expect(params[k]).toBeTruthy();
        }
    });

    it('default props produce a usable animation block (amplitudes 0 = animation off)', () => {
        const defaults = PD.getDefaultValues(SffShapeParams);
        expect(defaults.animation).toBeTruthy();
        expect(defaults.animation.wiggleAmplitude).toBe(0);
        expect(defaults.animation.tumbleAmplitude).toBe(0);
    });
});
