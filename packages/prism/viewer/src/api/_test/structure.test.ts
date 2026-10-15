/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';

// `structure` (operations.md §2, plan 2.2) on a headless viewer.

const dataDir = path.resolve(import.meta.dirname, '../../../../../../data/examples');
const file = (name: string, as: string) => new File([fs.readFileSync(path.join(dataDir, name))], as);

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

async function setup(...files: [string, string][]) {
  viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
  for (const [name, as] of files) await viewer.api.load.files({ files: [file(name, as)] as never });
  return viewer;
}

const reprTypes = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) =>
    s.components.flatMap((c) => c.representations.map((r) => (r.cell.transform.params as any).type.name as string)),
  );

describe('structure', () => {
  it('applies a preset by id, names the valid ones when wrong, and undoes', async () => {
    const v = await setup(['1cbs_updated.cif', '1cbs.cif']);
    const before = reprTypes(v).sort();
    const wrong = await v.api.structure.preset({ preset: 'nonsense' });
    expect(wrong).toMatchObject({ ok: false, code: 'invalid' });
    if (!wrong.ok) expect(wrong.hint).toMatch(/polymer-cartoon/);

    const r = await v.api.structure.preset({ preset: 'atomic-detail' });
    expect(r.ok).toBe(true);
    expect(reprTypes(v)).not.toEqual(before);
    await v.api.history.undo();
    expect(reprTypes(v).sort()).toEqual(before);
  });

  it('switches to an assembly and back, and says when it already shows one', async () => {
    const v = await setup(['1cbs_updated.cif', '1cbs.cif']);
    const missing = await v.api.structure.assembly({ id: '99' });
    expect(missing).toMatchObject({ ok: false, code: 'invalid' });
    const s = v.plugin.managers.structure.hierarchy.current.structures[0];
    const type = (s.cell.transform.params as any).type.name;
    const r = await v.api.structure.assembly({ id: '1' });
    expect(r.ok).toBe(true);
    const again = await v.api.structure.assembly({ id: '1' });
    expect(again.ok && again.did).toMatch(/already/);
    expect(type).toBeDefined();
  });

  it('refuses a model index past the end, with the range', async () => {
    const v = await setup(['1crn.cif', '1crn.cif']);
    const r = await v.api.structure.model({ index: 3 });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    const same = await v.api.structure.model({ index: 0 });
    expect(same.ok && same.did).toMatch(/already/);
  });

  it('asks which structure when several are loaded and none is named', async () => {
    const v = await setup(['1crn.cif', '1crn.cif'], ['1cbs_updated.cif', '1cbs.cif']);
    const r = await v.api.structure.model({ index: 0 });
    expect(r).toMatchObject({ ok: false, code: 'ambiguous' });
    if (!r.ok) expect(r.hint).toMatch(/1CRN/i);
    const named = await v.api.structure.model({ index: 0, structure: '1cbs' });
    expect(named.ok).toBe(true);
  });

  it('hides, shows and removes, each a step', async () => {
    const v = await setup(['1crn.cif', '1crn.cif'], ['1cbs_updated.cif', '1cbs.cif']);
    const hidden = await v.api.structure.visibility({ visible: false, structure: '1crn' });
    expect(hidden.ok).toBe(true);
    const ref = v.api.structure.list.value.find((x) => x.label === '1CRN')!.ref;
    expect(v.plugin.state.data.cells.get(ref)!.state.isHidden).toBe(true);
    const again = await v.api.structure.visibility({ visible: false, structure: '1crn' });
    expect(again.ok && again.did).toMatch(/already/);

    const removed = await v.api.structure.remove({ structure: '1crn' });
    expect(removed).toMatchObject({ ok: true });
    expect(v.api.structure.list.value.map((x) => x.label)).toEqual(['1CBS']);
    await v.api.history.undo();
    expect(v.api.structure.list.value).toHaveLength(2);
  });
});
