/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { afterEach, describe, expect, it } from 'vitest';
import { KinemageExtension } from '../prism-viewer.js';
import { Kin } from './kin-fixture.js';

// The kinemage extension in the Prism viewer on a headless viewer (parity.md §6), with a small .kin (kin-fixture.ts; no
// network, no canvas: the view operation says so).

const kin = (name = 'test.kin') => new File([Kin], name);

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(open = true) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), KinemageExtension] },
  });
  viewers.push(v);
  if (open) {
    const r = await v.api.load.files({ files: [kin()] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

type Info = {
  groups: { name: string; visible: boolean; subgroups: { name: string; visible: boolean }[] }[];
  masters: { name: string; visible: boolean }[];
  views: { name: string }[];
  animate: string[];
};
async function listing(v: PrismViewer): Promise<Info> {
  const r = await v.api.kinemage!.get({});
  if (!r.ok) throw new Error(r.message);
  return (r.value as unknown as { kinemages: Info[] }).kinemages[0];
}

const groupShown = async (v: PrismViewer, name: string) =>
  (await listing(v)).groups.find((g) => g.name === name)?.visible;

describe('the kinemage extension in the Prism viewer', () => {
  it('opens a .kin through load.files and lists its groups, subgroups, masters and views', async () => {
    const v = await viewer();
    const k = await listing(v);
    expect(k.groups.map((g) => g.name)).toEqual(['Contacts', 'Structure', 'Alternate']);
    expect(k.groups[1].subgroups.map((s) => s.name)).toEqual(['Backbone', 'Sidechains']);
    expect(k.masters.map((m) => m.name)).toEqual(['bad overlap', 'clashes']);
    expect(k.views.length).toBeGreaterThan(0);
    expect(k.animate).toEqual(['Structure', 'Alternate']);
    // Only the first animate frame shows at the start; the groups after it are off (the reader's rule).
    expect(k.groups.map((g) => [g.name, g.visible])).toEqual([
      ['Contacts', true],
      ['Structure', true],
      ['Alternate', false],
    ]);
    // Its shapes are in the scene (a representation each for lines, spheres and points).
    expect(v.plugin.state.data.selectQ((q) => q.ofType(PluginStateObject.Shape.Representation3D)).length).toBe(3);
  });

  it('lists nothing before a kinemage is loaded, and its operations say why they cannot run', async () => {
    const v = await viewer(false);
    const r = await v.api.kinemage!.get({});
    expect(r.ok && r.did).toBe('No kinemage is loaded.');
    expect(v.api.describe('kinemage.visibility')?.available).toBe('No kinemage is loaded.');
    const hide = await v.api.kinemage!.visibility({ name: 'Structure', visible: false });
    expect(hide.ok).toBe(false);
  });

  it('hides and shows a group as one step each, already-so changes nothing, and Undo restores it', async () => {
    const v = await viewer();
    const steps = v.managers.history.steps.value.length;
    const hide = await v.api.kinemage!.visibility({ kind: 'group', name: 'contacts', visible: false });
    expect(hide.ok, hide.ok ? '' : hide.message).toBe(true);
    expect(await groupShown(v, 'Contacts')).toBe(false);
    expect(v.managers.history.steps.value.length).toBe(steps + 1);

    const again = await v.api.kinemage!.visibility({ kind: 'group', name: 'Contacts', visible: false });
    expect(again.ok && again.unchanged).toBe(true);
    expect(v.managers.history.steps.value.length).toBe(steps + 1);

    await v.api.history.undo();
    expect(await groupShown(v, 'Contacts')).toBe(true);
  });

  it('hides a subgroup by its name and a master, and refuses unknown names', async () => {
    const v = await viewer();
    expect((await v.api.kinemage!.visibility({ kind: 'subgroup', name: 'Sidechains', visible: false })).ok).toBe(true);
    expect((await listing(v)).groups[1].subgroups.find((s) => s.name === 'Sidechains')?.visible).toBe(false);
    expect((await v.api.kinemage!.visibility({ kind: 'master', name: 'clashes', visible: false })).ok).toBe(true);
    expect((await listing(v)).masters.find((m) => m.name === 'clashes')?.visible).toBe(false);
    const bad = await v.api.kinemage!.visibility({ kind: 'master', name: 'nope', visible: false });
    expect(!bad.ok && bad.code).toBe('not-found');
    const none = await v.api.kinemage!.visibility({ kind: 'group', name: '' });
    expect(!none.ok && none.code).toBe('invalid');
  });

  it('steps the animation frames, showing one animate group at a time', async () => {
    const v = await viewer();
    const first = await v.api.kinemage!.animate({});
    expect(first.ok, first.ok ? '' : first.message).toBe(true);
    const a = await listing(v);
    expect(a.groups.filter((g) => g.visible).map((g) => g.name)).toEqual(['Contacts', 'Alternate']);
    await v.api.kinemage!.animate({});
    const b = await listing(v);
    expect(b.groups.filter((g) => g.visible).map((g) => g.name)).toEqual(['Contacts', 'Structure']);
    const two = await v.api.kinemage!.animate({ sequence: '2animate' });
    expect(!two.ok && two.code).toBe('invalid');
  });

  it('needs a canvas for a view, and says which kinemage when several are loaded', async () => {
    const v = await viewer();
    const r = await v.api.kinemage!.view({});
    expect(!r.ok && r.code).toBe('unavailable');
    expect((await v.api.load.files({ files: [kin('second.kin')] })).ok).toBe(true);
    const which = await v.api.kinemage!.visibility({ name: 'Structure', visible: false });
    expect(!which.ok && which.code).toBe('ambiguous');
    const second = await v.api.kinemage!.visibility({ name: 'Structure', visible: false, kinemage: '2' });
    expect(second.ok).toBe(true);
  });

  it('registers its commands and its skill', async () => {
    const v = await viewer(false);
    const ids = v.managers.commands.commands.value.map((c) => c.id);
    expect(ids).toEqual(expect.arrayContaining(['kinemage.visibility', 'kinemage.view', 'kinemage.animate']));
  });
});
