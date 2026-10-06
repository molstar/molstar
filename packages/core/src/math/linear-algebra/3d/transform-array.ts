import { Mat4 } from './mat4.js';

const identityTransform = new Float32Array(16);
Mat4.toArray(Mat4.identity(), identityTransform, 0);

export function fillIdentityTransform(transform: Float32Array, count: number) {
  for (let i = 0; i < count; i++) {
    transform.set(identityTransform, i * 16);
  }
  return transform;
}
