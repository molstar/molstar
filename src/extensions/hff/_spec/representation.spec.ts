/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Asserts that the SFF representation params expose the procedural-animation
 * controls (wiggle / tumble) so users can animate HFF surfaces from the
 * standard parameter panel.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { SffRepresentationParams } from '../representation';

describe('hff representation params', () => {
    it('exposes the animation PD.Group with wiggle and tumble controls', () => {
        const animParam = (SffRepresentationParams as any).animation;
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
        const defaults = PD.getDefaultValues(SffRepresentationParams);
        expect(defaults.animation).toBeTruthy();
        expect(defaults.animation.wiggleAmplitude).toBe(0);
        expect(defaults.animation.tumbleAmplitude).toBe(0);
    });

    it('exposes the visuals MultiSelect with mesh and wire options', () => {
        const v = (SffRepresentationParams as any).visuals;
        expect(v.type).toBe('multi-select');
        const optionNames = v.options.map((o: any[]) => o[0]);
        expect(optionNames).toEqual(expect.arrayContaining(['mesh', 'wire']));
    });
});
