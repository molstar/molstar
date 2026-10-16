/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { ModelExportExtension } from '../prism-viewer.js';
import { modelExportUI } from '../prism-viewer-ui.js';

// The Models section of the Export card over a headless viewer, in jsdom: Download says nothing of its own, the run
// path's message (the viewer's) is the one. No network.

const cbs = () =>
  new File([readFileSync(resolve(import.meta.dirname, '../../../../data/examples/1cbs_updated.cif'))], '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  document.body.innerHTML = '';
});

async function until(check: () => boolean, ms = 5000) {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error('timed out');
    await new Promise((r) => setTimeout(r, 5));
  }
}

describe('the Models section of the Export card', () => {
  it('leaves the status to the run path: one message for a failure, a warning, and for a success', async () => {
    const base = DefaultPrismViewerSpec();
    const v = await createPrismViewerHeadless({
      spec: { ...base, extensions: [...(base.extensions ?? []), ModelExportExtension] },
    });
    viewers.push(v);
    const loaded = await v.api.load.files({ files: [cbs()] });
    if (!loaded.ok) throw new Error(loaded.message);

    const ui = modelExportUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    const host = document.createElement('div');
    document.body.append(host);
    ui.sections![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    const download = () => [...host.querySelectorAll('button')].find((b) => b.textContent?.includes('Download'))!;
    const messages = () => v.managers.notifications.history.value;
    const settled = async () => {
      await until(() => !download().disabled && !download().getAttribute('aria-busy'));
      await new Promise((r) => setTimeout(r, 20));
    };

    // Nothing is selected: the run path warns, once.
    const only = [...host.querySelectorAll('[role=switch]')][0] as HTMLButtonElement;
    only.click();
    await until(() => only.getAttribute('aria-checked') === 'true');
    download().click();
    await until(() => messages().length > 0);
    await settled();
    expect(messages().map((m) => m.kind)).toEqual(['warning']);

    // It worked: one message again (the run path's), not a second 'success' from the section.
    only.click();
    await until(() => only.getAttribute('aria-checked') === 'false');
    (globalThis.URL as { createObjectURL?: unknown }).createObjectURL ??= () => 'blob:test';
    download().click();
    await until(() => messages().length > 1);
    await settled();
    expect(messages().map((m) => m.kind)).toEqual(['info', 'warning']);
  });
});
