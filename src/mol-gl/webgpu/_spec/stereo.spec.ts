/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../../mol-canvas3d/camera';
import { DefaultStereoCameraProps, StereoCamera } from '../../../mol-canvas3d/camera/stereo';
import { Mat4, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { WebGPUStereoCamera } from '../stereo';

describe('native stereo cameras', () => {
    it('retains established eye transforms, asymmetric projections and odd viewport splits', () => {
        for (const width of [256, 255]) {
            const parent = new Camera({ position: Vec3.create(4, 5, 20), radius: 10, radiusMax: 10 }, { x: 13, y: 17, width, height: 197 });
            parent.scale = 2; parent.update();
            const legacy = new StereoCamera(parent), native = new WebGPUStereoCamera();
            legacy.update(); native.update(parent, DefaultStereoCameraProps);
            for (const side of ['left', 'right'] as const) {
                expect(native[side].viewport).toEqual(legacy[side].viewport);
                expect(Mat4.areEqual(native[side].view, legacy[side].view, 1e-6)).toBe(true);
                expect(Mat4.areEqual(native[side].projection, legacy[side].projection, 1e-6)).toBe(true);
                expect(native[side].scale).toBe(2);
                expect(native[side].isAsymmetricProjection).toBe(true);
                expect(native[side].near).toBe(parent.near);
            }
        }
    });
    it('applies independent pixel jitter without replacing an eye projection', () => {
        const parent = new Camera({ position: Vec3.create(0, 0, 20), radius: 10, radiusMax: 10 }, { x: 0, y: 0, width: 256, height: 256 }); parent.update();
        const native = new WebGPUStereoCamera(); native.update(parent, DefaultStereoCameraProps);
        for (const eye of [native.left, native.right]) {
            const before = eye.project(Vec4(), Vec3.create(1, 2, 0)), projection = Mat4.clone(eye.projection);
            eye.viewOffset.enabled = true;
            Camera.setViewOffset(eye.viewOffset, eye.viewport.width, eye.viewport.height, 0.25, -0.375, eye.viewport.width, eye.viewport.height); eye.update();
            const jittered = eye.project(Vec4(), Vec3.create(1, 2, 0));
            expect(jittered[0]).toBeCloseTo(before[0] - 0.25, 6);
            expect(jittered[1]).toBeCloseTo(before[1] - 0.375, 6);
            eye.viewOffset.enabled = false; eye.update();
            expect(Mat4.areEqual(eye.projection, projection, 1e-6)).toBe(true);
        }
        expect(parent.viewOffset.enabled).toBe(false);
    });
    it('preserves crop offsets and round-trips points with the correct eye', () => {
        const parent = new Camera({ position: Vec3.create(0, 0, 20), radius: 10, radiusMax: 10 }, { x: 0, y: 0, width: 320, height: 240 }); parent.update();
        parent.viewOffset.enabled = true;
        Camera.setViewOffset(parent.viewOffset, 640, 480, 80, 40, 320, 240);
        const offset = { ...parent.viewOffset }, native = new WebGPUStereoCamera(); native.update(parent, DefaultStereoCameraProps);
        const point = Vec3.create(1, 2, 0);
        for (const eye of [native.left, native.right]) {
            expect(eye.viewOffset).toEqual(offset);
            expect(Vec3.distance(eye.unproject(Vec3(), eye.project(Vec4(), point)), point)).toBeLessThan(1e-6);
        }
        expect(parent.viewOffset).toEqual(offset);
    });
});
