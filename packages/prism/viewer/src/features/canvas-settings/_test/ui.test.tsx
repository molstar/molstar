/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DParams, type Canvas3DProps, type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it } from 'vitest';
import { expectOk, until } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { canvasSettingsUI } from '../ui.js';

// The Canvas card's sections over a stand-in canvas: an edit is debounced into one step, and the form keeps what the
// person is editing while the canvas reports older values (a drag is not pulled back).

let open: PrismViewer | undefined;
afterEach(async () => {
  document.body.innerHTML = '';
  await open?.dispose();
  open = undefined;
});

async function withCanvas() {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
  open = v;
  const props = structuredClone(PD.getDefaultValues(Canvas3DParams)) as Canvas3DProps;
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: {
      get props() {
        return structuredClone(props);
      },
      setProps: () => {},
      dispose() {},
    },
    configurable: true,
  });
  v.plugin.setCanvas3DProps = (change: PartialCanvas3DProps) => {
    if (change.transparentBackground !== undefined) props.transparentBackground = change.transparentBackground;
    v.plugin.events.canvas3d.settingsUpdated.next(void 0);
  };
  // The canvas source reads once the canvas is there.
  (v.plugin.behaviors.canvas3d.initialized as unknown as { next(v: boolean): void }).next(true);
  return { v, props };
}

describe('the Canvas card', () => {
  it('sends an edit after its debounce as one step, and keeps the form while the edit is on its way', async () => {
    const { v, props } = await withCanvas();
    const ui = canvasSettingsUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    const host = document.createElement('div');
    document.body.append(host);
    ui.sections![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    const toggle = () => host.querySelector('[role=switch]') as HTMLButtonElement;
    expect(toggle().getAttribute('aria-checked')).toBe('false');

    toggle().click();
    expect(toggle().getAttribute('aria-checked')).toBe('true');
    // The canvas reports an older value while the edit waits: the form keeps the person's.
    v.plugin.events.canvas3d.settingsUpdated.next(void 0);
    expect(props.transparentBackground).toBe(false);
    expect(toggle().getAttribute('aria-checked')).toBe('true');

    await until(() => props.transparentBackground === true);
    await until(() => v.managers.history.running.value === 0);
    expect(v.managers.history.steps.value).toHaveLength(1);
    expect(v.managers.history.steps.value[0].label).toMatch(/Transparent background on/);

    // Landed: the form follows the canvas again (Undo, a command).
    props.transparentBackground = false;
    v.plugin.events.canvas3d.settingsUpdated.next(void 0);
    await until(() => toggle().getAttribute('aria-checked') === 'false');
    expectOk(await v.api.canvas.illumination!({ on: true }));
  });
});
