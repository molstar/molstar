/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

export type AnimationData = {
    uWiggleSpeed: ValueCell<number>,
    uWiggleAmplitude: ValueCell<number>,
    uWiggleFrequency: ValueCell<number>,
    uWiggleMode: ValueCell<number>,
    uTumbleSpeed: ValueCell<number>,
    uTumbleAmplitude: ValueCell<number>,
    uTumbleFrequency: ValueCell<number>,
}

export function getAnimationParam() {
    return PD.Group({
        wiggleMode: PD.Select('position', [['position', 'Position'], ['group', 'Group']] as const, { description: 'Noise seeding mode. Position: spatially correlated (nearby atoms move together). Group: per-group independent noise.' }),
        wiggleSpeed: PD.Numeric(7, { min: 0, max: 10, step: 0.1 }, { description: 'Speed of vertex wiggle animation.' }),
        wiggleAmplitude: PD.Numeric(0, { min: 0, max: 5, step: 0.01 }, { description: 'Amplitude of vertex wiggle animation.' }),
        wiggleFrequency: PD.Numeric(0.2, { min: 0.01, max: 2, step: 0.01 }, { description: 'Spatial frequency of vertex wiggle noise (position mode). Lower values correlate nearby atoms more.' }),
        tumbleSpeed: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, { description: 'Speed of instance tumble animation.' }),
        tumbleAmplitude: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, { description: 'Amplitude of instance tumble animation. In Ångströms of implied surface displacement.' }),
        tumbleFrequency: PD.Numeric(0.2, { min: 0, max: 2, step: 0.01 }, { description: 'Spatial frequency multiplier for tumble noise.' }),
    });
}
export type AnimationParam = ReturnType<typeof getAnimationParam>
export type AnimationProps = AnimationParam['defaultValue'];

export function areAnimationPropsEqual(a: AnimationProps, b: AnimationProps): boolean {
    return a.wiggleMode === b.wiggleMode
        && a.wiggleSpeed === b.wiggleSpeed
        && a.wiggleAmplitude === b.wiggleAmplitude
        && a.wiggleFrequency === b.wiggleFrequency
        && a.tumbleSpeed === b.tumbleSpeed
        && a.tumbleAmplitude === b.tumbleAmplitude
        && a.tumbleFrequency === b.tumbleFrequency;
}

export function createAnimationValues(props: AnimationProps) {
    return {
        uWiggleSpeed: ValueCell.create(props.wiggleSpeed),
        uWiggleAmplitude: ValueCell.create(props.wiggleAmplitude),
        uWiggleFrequency: ValueCell.create(props.wiggleFrequency),
        uWiggleMode: ValueCell.create(props.wiggleMode === 'position' ? 0 : 1),
        uTumbleSpeed: ValueCell.create(props.tumbleSpeed),
        uTumbleAmplitude: ValueCell.create(props.tumbleAmplitude),
        uTumbleFrequency: ValueCell.create(props.tumbleFrequency),
    };
}

export function updateAnimationValues(values: AnimationData, props: AnimationProps) {
    ValueCell.updateIfChanged(values.uWiggleSpeed, props.wiggleSpeed);
    ValueCell.updateIfChanged(values.uWiggleAmplitude, props.wiggleAmplitude);
    ValueCell.updateIfChanged(values.uWiggleFrequency, props.wiggleFrequency);
    ValueCell.updateIfChanged(values.uWiggleMode, props.wiggleMode === 'position' ? 0 : 1);
    ValueCell.updateIfChanged(values.uTumbleSpeed, props.tumbleSpeed);
    ValueCell.updateIfChanged(values.uTumbleAmplitude, props.tumbleAmplitude);
    ValueCell.updateIfChanged(values.uTumbleFrequency, props.tumbleFrequency);
}
