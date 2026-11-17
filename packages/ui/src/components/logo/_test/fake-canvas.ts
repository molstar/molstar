/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { vi } from 'vitest';

/** A Path2D that records nothing (jsdom and Node have none): the glass only builds and draws paths. */
export class FakePath2D {
  moveTo(): void {}
  lineTo(): void {}
  closePath(): void {}
}

/** A 2D context that counts what the glass asks of it, and draws nothing. */
export function fakeContext() {
  const counts = { stroke: 0, drawImage: 0, gradients: 0 };
  // Its numbers start fractional (a field that held an integer and is then written a fraction changes shape: an
  // allocation that is the fake's, not the glass's).
  const ctx = {
    globalAlpha: 0.5,
    lineWidth: 0.5,
    lineJoin: 'miter',
    strokeStyle: '' as unknown,
    save() {},
    restore() {},
    translate() {},
    scale() {},
    setTransform() {},
    clearRect() {},
    clip() {},
    stroke() {
      counts.stroke++;
    },
    drawImage() {
      counts.drawImage++;
    },
    createLinearGradient() {
      counts.gradients++;
      return { addColorStop() {} };
    },
  };
  return { ctx: ctx as unknown as CanvasRenderingContext2D, counts };
}

/** Path2D where there is none; restored after the file by the setup (`vi.stubGlobal`). */
export function stubPath2D(): void {
  if (typeof (globalThis as any).Path2D !== 'function') vi.stubGlobal('Path2D', FakePath2D);
}
