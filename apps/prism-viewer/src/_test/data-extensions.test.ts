/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { DataOnlyExtensions } from '../data-extensions.js';

// Mol*'s data-only extensions (plan 5.2) on a headless viewer: they register what the v5 behaviors did, and what a
// structure carries of them is listed by `annotations.list`. No network: nothing here applies a validation theme (it
// would fetch the report).

const cbs = () =>
  new File([readFileSync(new URL('../../../../data/examples/1cbs_updated.cif', import.meta.url))], '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(load = false) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ...DataOnlyExtensions] },
  });
  viewers.push(v);
  if (load) {
    const r = await v.api.load.files({ files: [cbs()] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

describe('the data-only extensions in the Prism viewer', () => {
  it('lists one extension each, with unique ids', () => {
    expect(DataOnlyExtensions.map((e) => e.id)).toEqual([
      'dnatco',
      'rcsb',
      'pdbe',
      'wwpdb',
      'sb-ncbr',
      'g3d',
      'backgrounds',
    ]);
  });

  it('registers their color themes, formats and background styles in the plugin', async () => {
    const { PluginConfig } = await import('@molstar/plugin/config');
    const v = await viewer();
    const themes = v.plugin.representation.structure.themes.colorThemeRegistry;
    for (const name of [
      'rcsb-geometry-quality',
      'rcsb-density-fit',
      'rcsb-random-coil-index',
      'pdbe-structure-quality-report',
      'sb-ncbr-partial-charges',
    ]) {
      expect(themes.has(name), name).toBe(true);
    }
    expect(v.plugin.dataFormats.list.map((f) => f.name)).toContain('g3d');
    expect(v.plugin.config.get(PluginConfig.Background.Styles)?.length).toBe(4);
  });

  it('lists the validation reports for an entry with a PDB id, and no partial charges for it', async () => {
    const v = await viewer(true);
    const r = await v.api.annotations.list({});
    const text = JSON.stringify(r.ok && r.value);
    expect(text).toContain('rcsb-geometry-quality');
    expect(text).toContain('pdbe-structure-quality');
    // 1CBS is an X-ray entry: density fit may apply, the random coil index (NMR) and partial charges do not.
    expect(text).not.toContain('rcsb-random-coil-index');
    expect(text).not.toContain('partial-charges');
    const none = await v.api.annotations.apply({ annotation: 'partial-charges' });
    expect(!none.ok && none.code).toBe('not-found');
  });
});
