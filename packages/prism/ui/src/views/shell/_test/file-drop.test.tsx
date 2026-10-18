/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { flush } from '../../../dom/scope.js';
import { FileDropOverlay } from '../file-drop.js';

/** A drag event carrying `types` and `files` (jsdom has no DataTransfer). */
function drag(type: string, files: File[]) {
  const e = new Event(type, { bubbles: true, cancelable: true });
  Object.defineProperty(e, 'dataTransfer', { value: { types: files.length ? ['Files'] : ['text/plain'], files } });
  return e;
}

describe('FileDropOverlay', () => {
  it('shows only for drags carrying files, hands dropped files over, hides on leave and drop', () => {
    const root = document.createElement('div');
    document.body.append(root);
    const dropped: File[][] = [];
    render(() => <FileDropOverlay target={root} onDrop={(f) => dropped.push([...f])} />, root);
    const overlay = root.lastElementChild as HTMLElement;
    expect(overlay.hidden).toBe(true);

    root.dispatchEvent(drag('dragenter', []));
    flush();
    expect(overlay.hidden).toBe(true);

    const file = new File(['data_1TQN'], '1tqn.cif');
    root.dispatchEvent(drag('dragenter', [file]));
    flush();
    expect(overlay.hidden).toBe(false);
    overlay.dispatchEvent(drag('dragleave', [file]));
    flush();
    expect(overlay.hidden).toBe(true);

    root.dispatchEvent(drag('dragenter', [file]));
    flush();
    const drop = drag('drop', [file]);
    overlay.dispatchEvent(drop);
    flush();
    expect(drop.defaultPrevented).toBe(true);
    expect(overlay.hidden).toBe(true);
    expect(dropped).toEqual([[file]]);
  });
});
