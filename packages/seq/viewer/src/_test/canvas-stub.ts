/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';

/**
 * jsdom has no canvas: its `getContext` logs "Not implemented" through the virtual console with an error and a stack,
 * every call. A mounted viewer asks for a context once a canvas, so a file that mounts one without reading its drawing
 * stubs it with this: a recording context kept per canvas. Undone by `vi.restoreAllMocks`.
 */
export function stubCanvasContext(): void {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as { __ctx?: RecordingContext }).__ctx ??= new RecordingContext()) as never;
  });
}
