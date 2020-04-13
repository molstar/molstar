import { Vec3 } from '../3d';

describe('vec3', () => {
    const vec1 = [ 1, 2, 3 ] as Vec3;
    const vec2 = [ 2, 3, 1 ] as Vec3;
    const orthVec1 = [ 0, 1, 0 ] as Vec3;
    const orthVec2 = [ 1, 0, 0 ] as Vec3;

    it('angle calculation', () => {
        expect(Vec3.angle(vec1, vec1) * 360 / (2 * Math.PI)).toBe(0.0);
        expect(Vec3.angle(orthVec1, orthVec2) * 360 / (2 * Math.PI)).toBe(90.0);
        expect(Vec3.angle(vec1, vec2)).toBeCloseTo(0.666946);
    });
});