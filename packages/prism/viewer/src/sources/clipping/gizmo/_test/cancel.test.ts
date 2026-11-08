/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Representation } from '@molstar/graphics/repr/representation';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import { Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { cbs, headlessViewer, until } from '../../../../_test/harness.js';
import type { PrismViewer } from '../../../../viewer.js';
import { type Cut, type CutObject, findCuts } from '../../cuts.js';
import { ClipBinding, type ClipItem } from '../binding.js';

// A drag's ends on the clipping binding, over a headless viewer: a preview is on the live representation only, a cancel
// puts back what the state holds, a change under a drag invalidates it.

async function setup() {
  const v = await headlessViewer({ files: [cbs()] });
  const added = await v.api.clipping.add({ shape: 'cube', target: 'polymer' });
  if (!added.ok) throw new Error(added.message);
  const cut = added.value as Cut;
  const ref = cut.targets[0].ref;
  const binding = new ClipBinding(v.plugin, new Subject<readonly ClipItem[]>());
  const moved: CutObject = { ...cut.object, position: Vec3.create(11, 22, 33) };
  return { v, cut, ref, binding, moved };
}

const repr = (v: PrismViewer, ref: string) =>
  (v.plugin.state.data.cells.get(ref)!.obj!.data as { repr: Representation.Any }).repr;
const shown = (v: PrismViewer, ref: string) =>
  (repr(v, ref).props as unknown as { clip: { objects: CutObject[] } }).clip.objects[0].position;
const stored = (v: PrismViewer, ref: string) =>
  RepresentationClip.get(v.plugin.state.data.cells.get(ref)!.transform.params)!.objects[0].position;
const atMoved = (p: Vec3) => Vec3.exactEquals(p, Vec3.create(11, 22, 33));

describe('a drag on a cut', () => {
  it('previews on the representation only, and a cancel puts back what the state holds', async () => {
    const { v, cut, ref, binding, moved } = await setup();
    binding.begin('a', cut.object);
    await binding.preview('a', cut.object, moved);
    expect(atMoved(shown(v, ref))).toBe(true);
    expect(atMoved(stored(v, ref))).toBe(false);

    binding.cancel();
    await until(() => !atMoved(shown(v, ref)));
    expect(Vec3.exactEquals(shown(v, ref), cut.object.position)).toBe(true);
    binding.dispose();
  });

  it('keeps a commit not taken cancelable, and a commit without a commit function is not taken', async () => {
    const { v, cut, ref, binding, moved } = await setup();
    binding.begin('a', cut.object);
    await binding.preview('a', cut.object, moved);
    expect(await binding.commit('a', cut.object, moved)).toBe(false);
    binding.setCommit(async () => false);
    expect(await binding.commit('a', cut.object, moved)).toBe(false);
    expect(atMoved(stored(v, ref))).toBe(false);
    binding.cancel();
    await until(() => !atMoved(shown(v, ref)));
    binding.dispose();
  });

  it('commits as one step, which a later cancel leaves alone', async () => {
    const { v, cut, ref, binding, moved } = await setup();
    binding.setCommit(
      async (_, next) =>
        (await v.api.clipping.update({ cut: 1, position: [...next.position] as [number, number, number] })).ok,
    );
    binding.begin('a', cut.object);
    await binding.preview('a', cut.object, moved);
    const steps = v.managers.history.steps.value.length;
    expect(await binding.commit('a', cut.object, moved)).toBe(true);
    expect(v.managers.history.steps.value.length).toBe(steps + 1);
    expect(atMoved(stored(v, ref))).toBe(true);
    // the commit ended the drag: a cancel touches no representation (a revert would show what the state now holds
    // anyway, so the screen alone cannot tell)
    const update = vi.spyOn(repr(v, ref), 'createOrUpdate');
    binding.cancel();
    expect(update).not.toHaveBeenCalled();
    expect(atMoved(shown(v, ref))).toBe(true);
    binding.dispose();
  });

  it('is invalidated by a change of what it is on, and its cancel keeps that change', async () => {
    const { v, cut, ref, binding, moved } = await setup();
    const invalidated: string[] = [];
    binding.invalidated.subscribe((k) => invalidated.push(k));
    binding.begin('a', cut.object);
    await binding.preview('a', cut.object, moved);
    expect((await v.api.clipping.update({ cut: 1, position: [1, 2, 3] })).ok).toBe(true);
    expect(invalidated).toContain('a');
    binding.cancel();
    await until(() => Vec3.exactEquals(shown(v, ref), Vec3.create(1, 2, 3)));
    binding.dispose();
  });

  it('is not invalidated by a change elsewhere, or once it ended', async () => {
    const { v, cut, binding } = await setup();
    const invalidated: string[] = [];
    binding.invalidated.subscribe((k) => invalidated.push(k));
    binding.begin('a', { ...cut.object, tag: 'another' });
    expect((await v.api.clipping.update({ cut: 1, position: [1, 2, 3] })).ok).toBe(true);
    expect(invalidated).toEqual([]);
    binding.begin('a', findCuts(v.plugin)[0].object);
    binding.cancel();
    expect((await v.api.clipping.update({ cut: 1, position: [4, 5, 6] })).ok).toBe(true);
    expect(invalidated).toEqual([]);
    binding.dispose();
  });
});
