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
import { ModelArchiveExtension, PlddtTrack } from '../prism-viewer.js';

// Model Archive's quality in the Prism viewer (plan 5.1) on a headless viewer: AlphaFold DB's P04637 model (recorded,
// Seq*'s fixture; pLDDT in `ma_qa_metric_local`) and 1CBS (an experimental structure). No network.

const fixture = (path: string, name: string) => new File([readFileSync(new URL(path, import.meta.url))], name);
const alphafold = () =>
  fixture(
    '../../../../packages/sequences/molstar-poc/_test/fixtures/structures/AF-P04637-F1-model_v6.bcif',
    'AF-P04637-F1.bcif',
  );
const cbs = () => fixture('../../../../data/examples/1cbs_updated.cif', '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(files: File[], withExtension = true) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ...(withExtension ? [ModelArchiveExtension] : [])] },
  });
  viewers.push(v);
  for (const f of files) {
    const r = await v.api.load.files({ files: [f] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

async function until(check: () => boolean, ms = 5000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

const polymerThemes = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) =>
    s.components.flatMap((c) =>
      c.representations.map((r) => (r.cell.transform.params as { colorTheme: { name: string } }).colorTheme.name),
    ),
  );

describe('the model-archive extension in the Prism viewer', () => {
  it('lists pLDDT for a predicted model and colors the polymer by it, as one step', async () => {
    const v = await viewer([alphafold()]);
    const list = await v.api.annotations.list({});
    expect(list.ok && JSON.stringify(list.value)).toMatch(/plddt/);
    const steps = v.managers.history.steps.value.length;
    const r = await v.api.annotations.apply({ annotation: 'plddt' });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(polymerThemes(v)).toContain('plddt-confidence');
    expect(v.managers.history.steps.value.length).toBe(steps + 1);
    // QMEAN is not what an AlphaFold model carries.
    expect((await v.api.annotations.apply({ annotation: 'qmean' })).ok).toBe(false);
  });

  it('adds the pLDDT track to the sequence, with its menu and its bands, only for a model that carries it', async () => {
    const v = await viewer([alphafold()]);
    const seq = (await v.views.sequence!.load()).viewers.get();
    await until(() => JSON.stringify(seq.seq.document.shown.value).includes(`"ref":"${PlddtTrack}"`));
    expect(seq.legend.value.sections?.map((s) => s.title)).toEqual(['Secondary structure', 'Model confidence']);
    expect(seq.trackActions(PlddtTrack).map((a) => a.label)).toEqual(['Color by pLDDT']);
    seq.trackActions(PlddtTrack)[0].run();
    await until(() => polymerThemes(v).includes('plddt-confidence'));
    // An experimental structure: no track, no bands.
    const other = await viewer([cbs()]);
    const plain = (await other.views.sequence!.load()).viewers.get();
    await until(() => plain.shown.value !== undefined);
    expect(JSON.stringify(plain.seq.document.shown.value)).not.toContain(PlddtTrack);
    expect(plain.legend.value.sections?.map((s) => s.title)).toEqual(['Secondary structure']);
  });

  it('brings its commands, and without it the viewer has no pLDDT at all', async () => {
    const v = await viewer([]);
    expect(v.managers.commands.get('model-archive.color-plddt')).toBeDefined();
    const bare = await viewer([alphafold()], false);
    expect((await bare.api.annotations.apply({ annotation: 'plddt' })).ok).toBe(false);
  });
});
