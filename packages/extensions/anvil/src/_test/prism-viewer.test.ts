/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { Color } from '@molstar/core/util/color';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { MembraneInfo } from '../operations.js';
import { AnvilExtension } from '../prism-viewer.js';
import { MembraneOrientation } from '../prop.js';
import type { MembraneOrientationProps } from '../representation.js';

// ANVIL in the Prism viewer on a headless viewer (parity.md §6): 1GRM (gramicidin A, a channel in the membrane, two
// 16-residue chains), 1CBS (a soluble protein: ANVIL still places a membrane, as it does for any protein) and a .gro
// (no protein chains it can use). Computed locally; no network.

const fixture = (name: string) =>
  new File([readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url))], name);

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(files: string[]) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), AnvilExtension] },
  });
  viewers.push(v);
  for (const f of files) {
    const r = await v.api.load.files({ files: [fixture(f)] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

const membranes = (v: PrismViewer) =>
  [...v.plugin.state.data.cells.values()].filter((c) =>
    c.transform.tags?.includes(MembraneOrientation.Tag.Representation),
  );

const look = (v: PrismViewer) => membranes(v)[0].transform.params as MembraneOrientationProps;

describe('the anvil extension in the Prism viewer', () => {
  it('places the membrane of gramicidin across the channel and shows it, as one step that Undo removes', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const steps = v.managers.history.steps.value.length;
    const r = await v.api.membrane!.show({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    const info = (r.ok ? r.value : undefined) as MembraneInfo;
    // The channel runs along z: the membrane's normal does too, and most of its 32 residues are in it.
    expect(Math.abs(info.normal[2])).toBeGreaterThan(0.99);
    expect(info.thickness).toBeGreaterThanOrEqual(20);
    expect(info.thickness).toBeLessThanOrEqual(40);
    expect(info.embeddedResidues).toBeGreaterThan(16);
    expect(info.shown).toBe(true);
    expect(r.ok && r.did).toMatch(/^Showed the membrane of 1GRM: \d+(\.\d)? Å thick, normal \(/);
    expect(membranes(v).length).toBe(1);
    expect(membranes(v)[0].obj?.label).toBe('Membrane Orientation');
    expect(v.managers.history.steps.value.length).toBe(steps + 1);

    const again = await v.api.membrane!.show({});
    expect(again.ok && again.unchanged).toBe(true);

    await v.api.history.undo();
    expect(membranes(v).length).toBe(0);
  });

  it('says where the membrane is without showing it', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const r = await v.api.membrane!.get({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    const info = (r.ok ? r.value : undefined) as MembraneInfo;
    expect(info.shown).toBe(false);
    expect(info.settings.adjust).toBe(14);
    expect(r.ok && r.did).toMatch(/; not shown\.$/);
    expect(membranes(v).length).toBe(0);
  });

  it("fits it again with other settings, kept in the structure's custom properties, and Undo restores them", async () => {
    const v = await viewer(['1cbs_updated.cif']);
    const m = v.api.membrane!;
    const first = await m.show({});
    expect(first.ok, first.ok ? '' : first.message).toBe(true);
    const before = (first.ok ? first.value : undefined) as MembraneInfo;
    const update = vi.spyOn(membranes(v)[0].transform.transformer.definition, 'update');

    const r = await m.fit({ minThickness: 10, maxThickness: 30, adjust: 0 });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    const after = (r.ok ? r.value : undefined) as MembraneInfo;
    expect(after.settings).toMatchObject({ minThickness: 10, maxThickness: 30, adjust: 0 });
    expect(after.thickness).not.toBe(before.thickness);
    // The planes follow: the membrane's cell is updated, not left with the old placement.
    expect(update).toHaveBeenCalled();
    expect(membranes(v).length).toBe(1);
    const props = v.plugin.managers.structure.hierarchy.current.structures[0].properties!.cell.transform.params as {
      autoAttach: string[];
      properties: Record<string, { adjust: number }>;
    };
    expect(props.autoAttach).toContain('anvil_computed_membrane_orientation');
    expect(props.properties.anvil_computed_membrane_orientation.adjust).toBe(0);

    const same = await m.fit({ minThickness: 10, maxThickness: 30, adjust: 0 });
    expect(same.ok && same.unchanged).toBe(true);

    await v.api.history.undo();
    const back = await m.get({});
    expect(back.ok && (back.value as MembraneInfo).settings.adjust).toBe(14);
    expect(back.ok && (back.value as MembraneInfo).thickness).toBe(before.thickness);
    update.mockRestore();
  });

  it('fits and shows a membrane that is not shown yet', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const r = await v.api.membrane!.fit({ tmdetDefinition: true });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(membranes(v).length).toBe(1);
  });

  it('sets how it looks, keeping what is not given', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const m = v.api.membrane!;
    const blocked = await m.options({ opacity: 0.2 });
    expect(!blocked.ok && blocked.code).toBe('blocked');
    await m.show({});
    const r = await m.options({ rims: false, color: 'steelblue', opacity: 0.2 });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(look(v)).toMatchObject({ visuals: ['bilayer-planes'], sectorOpacity: 0.2, radiusFactor: 1.2 });
    expect(Color.toHexStyle(look(v).color)).toBe('#4682b4');
    const same = await m.options({ opacity: 0.2 });
    expect(same.ok && same.unchanged).toBe(true);
    const none = await m.options({ planes: false });
    expect(!none.ok && none.code).toBe('invalid');
    const bad = await m.options({ color: 'not-a-color' });
    expect(!bad.ok && bad.code).toBe('invalid');
  });

  it('hides it, and says so when there is none', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const m = v.api.membrane!;
    const none = await m.hide({});
    expect(none.ok && none.unchanged).toBe(true);
    await m.show({});
    const r = await m.hide({});
    expect(r.ok && r.did).toBe('Hid the membrane of 1GRM.');
    expect(membranes(v).length).toBe(0);
  });

  it('refuses a structure without protein chains, and brings its commands and skill', async () => {
    const v = await viewer(['test.gro']);
    const r = await v.api.membrane!.show({});
    expect(!r.ok && r.code).toBe('blocked');
    expect(membranes(v).length).toBe(0);
    for (const id of ['show', 'hide', 'fit', 'options', 'get']) {
      expect(v.managers.commands.get(`membrane.${id}`)).toBeDefined();
      expect(v.api.describe(`membrane.${id}`)!.definition.approval).toBeUndefined();
    }
    expect(v.api.describe('membrane.get')!.definition.read).toBe(true);
    expect(AnvilExtension.skills?.[0].tools.map((t) => t.operation)).toEqual([
      'membrane.get',
      'membrane.show',
      'membrane.hide',
      'membrane.fit',
      'membrane.options',
    ]);
  });
});
