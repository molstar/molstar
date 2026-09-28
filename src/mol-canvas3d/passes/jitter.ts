/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Camera } from '../camera';
import { StereoCamera } from '../camera/stereo';

export const JitterVectors = [
    [
        [0, 0]
    ],
    [
        [0, 0], [-4, -4]
    ],
    [
        [0, 0], [6, -2], [-6, 2], [2, 6]
    ],
    [
        [0, 0], [-1, 3], [5, 1], [-3, -5],
        [-5, 5], [-7, -1], [3, 7], [7, -7]
    ],
    [
        [0, 0], [-1, -3], [-3, 2], [4, -1],
        [-5, -2], [2, 5], [5, 3], [3, -5],
        [-2, 6], [0, -7], [-4, -6], [-6, 4],
        [-8, 0], [7, -4], [6, 7], [-7, -8]
    ],
    [
        [0, 0], [-7, -5], [-3, -5], [-5, -4],
        [-1, -4], [-2, -2], [-6, -1], [-4, 0],
        [-7, 1], [-1, 2], [-6, 3], [-3, 3],
        [-7, 6], [-3, 6], [-5, 7], [-1, 7],
        [5, -7], [1, -6], [6, -5], [4, -4],
        [2, -3], [7, -2], [1, -1], [4, -1],
        [2, 1], [6, 2], [0, 4], [4, 4],
        [2, 5], [7, 5], [5, 6], [3, 7]
    ]
];

JitterVectors.forEach(offsetList => {
    offsetList.forEach(offset => {
        // 0.0625 = 1 / 16
        offset[0] *= 0.0625;
        offset[1] *= 0.0625;
    });
});

export function getJitterOffsets(sampleLevel: number) {
    return JitterVectors[Math.max(0, Math.min(sampleLevel, 5))];
}

export function getTemporalSamplesPerFrame(sampleLevel: number) {
    return Math.pow(2, Math.max(0, sampleLevel - 2));
}

export function setJitter(camera: Camera | StereoCamera, offset: number[]) {
    const { width, height } = camera.viewport;
    camera.viewOffset.enabled = true;
    Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height);
    camera.update();
}

export function clearJitter(camera: Camera | StereoCamera) {
    camera.viewOffset.enabled = false;
    camera.update();
}
