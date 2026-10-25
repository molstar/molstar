/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import { ValueCell } from '@molstar/core/util/value-cell';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { SegmentedArray } from '@molstar/graphics/gl/segmented-array';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { WebGPU } from './webgpu-module.js';

/** a compact grid of atoms, filling the view */
const Pdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, (Math.floor(i / 4) % 4) * 2.5, Math.floor(i / 16) * 2.5];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

/** two chains of the grid's atoms, side by side: an object each */
const TwoChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const chain = i < 24 ? 'A' : 'B';
  const residue = String((i % 24) + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5 + (i < 24 ? 0 : 12), (Math.floor(i / 4) % 6) * 2.5, 0];
  return `ATOM  ${serial}  CA  ALA ${chain}${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 64;
/** pick positions over the scene and beside it */
const Points = [
  [32, 32],
  [20, 40],
  [44, 24],
  [28, 12],
  [2, 2],
];

// the pick passes replayed from the renderer's draw lists find what they find drawn directly
describe.skipIf(!WebGPU.module)('draw lists on WebGPU', () => {
  it('pick as the direct passes do, replayed', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'cartoon' });
      const canvas = plugin.canvas3d!;
      const { stats } = canvas.graphics;
      const pixelRatio = canvas.graphics.pixelRatio;

      /** the picks at `Points`, each after a drawn frame (the pick passes render again after one) */
      let time = now();
      /** the draw lists replayed by the picks (not by the frames between them) */
      let pickReplays = 0;
      const pickAll = async () => {
        const ids: string[] = [];
        for (const [x, y] of Points) {
          canvas.commit(true);
          canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
          const replayed = stats.drawLists.replayed;
          const pick = await canvas.asyncIdentify(Vec2.create(x / pixelRatio, y / pixelRatio))?.wait();
          pickReplays += stats.drawLists.replayed - replayed;
          const id = pick?.id;
          ids.push(id ? `${id.objectId}/${id.instanceId}/${id.groupId}` : 'none');
        }
        return ids;
      };

      // compiled, then drawn directly
      await pickAll();
      const direct = await pickAll();
      expect(direct.filter((id) => id !== 'none').length).toBeGreaterThan(1);

      canvas.setProps({ renderer: { ...canvas.props.renderer, drawLists: true } });
      // recorded, then replayed twice
      expect(await pickAll()).toEqual(direct);
      pickReplays = 0;
      expect(await pickAll()).toEqual(direct);
      expect(await pickAll()).toEqual(direct);
      expect(pickReplays).toBeGreaterThan(0);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('record again where the shared quad indices were replaced (grown)', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        renderer: { ...canvas.props.renderer, drawLists: true },
      });
      const graphics = canvas.graphics as WebGPUContext;
      const { stats } = graphics;
      const passes = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
        return Array.from(await graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
      };
      // compiled, recorded, then replayed with the first quad indices
      let before = await frame();
      for (let i = 0; i < 8 && stats.drawLists.replayed === 0; ++i) before = await frame();
      expect(stats.drawLists.replayed).toBeGreaterThan(0);

      // replaced (as by more spheres than they hold), the scene unchanged: no recording binds the released buffer
      const { recorded, replayed } = stats.drawLists;
      const id = graphics.quadIndices.id;
      graphics.quadIndices.get(graphics.quadIndices.get(1)!.itemCount + 1);
      expect(graphics.quadIndices.id).not.toBe(id);
      const after = await frame();
      expect(stats.drawLists.replayed).toBe(replayed);
      expect(stats.drawLists.recorded).toBeGreaterThan(recorded);
      expect(after).toEqual(before);
      // replayed again with the new buffer
      expect(await frame()).toEqual(before);
      expect(stats.drawLists.replayed).toBeGreaterThan(replayed);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('replay with changed markers written, not recorded again (hover)', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      // one frame per state: no accumulating samples
      canvas.setProps({ multiSample: { ...canvas.props.multiSample, mode: 'off' } });
      const { stats } = canvas.graphics;
      const passes = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        return Array.from(await canvas.graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
      };
      const EveryLoci = { kind: 'every-loci' as const };
      /** the frames of the scene highlighted and not, alternating */
      const toggle = async () => {
        const frames: number[][] = [];
        for (let i = 0; i < 4; ++i) {
          canvas.mark({ loci: EveryLoci }, i % 2 === 0 ? MarkerAction.Highlight : MarkerAction.RemoveHighlight);
          frames.push(await frame());
        }
        return frames;
      };

      // compiled, then drawn directly
      await toggle();
      const direct = await toggle();
      expect(direct[0]).not.toEqual(direct[1]);

      canvas.setProps({ renderer: { ...canvas.props.renderer, drawLists: true } });
      // recorded
      await toggle();
      const { recorded, flushed } = stats.drawLists;
      // the marker changes are data: replayed with it written
      expect(await toggle()).toEqual(direct);
      expect(stats.drawLists.recorded).toBe(recorded);
      expect(stats.drawLists.flushed).toBeGreaterThan(flushed);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('replay the marking passes of marking-only redraws (the camera resting, no new frame)', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        renderer: { ...canvas.props.renderer, drawLists: true },
      });
      const { stats } = canvas.graphics;
      let time = now();
      canvas.commit(true);
      canvas.requestDraw();
      canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
      await canvas.graphics.resources.pending.whenReady();
      const EveryLoci = { kind: 'every-loci' as const };
      /** highlighted and not, alternating, each redrawn as marking only (no draw requested) */
      const toggle = async () => {
        for (let i = 0; i < 4; ++i) {
          canvas.mark({ loci: EveryLoci }, i % 2 === 0 ? MarkerAction.Highlight : MarkerAction.RemoveHighlight);
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
        }
      };

      // compiled, then recorded
      await toggle();
      await toggle();
      const { recorded, replayed } = stats.drawLists;
      await toggle();
      await toggle();
      expect(stats.drawLists.replayed).toBeGreaterThan(replayed);
      expect(stats.drawLists.recorded).toBe(recorded);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('hover single groups and whole objects: changing chunks of few objects drawn directly, nothing recorded again', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        renderer: { ...canvas.props.renderer, drawLists: true },
      });
      const { stats } = canvas.graphics;
      let time = now();
      canvas.commit(true);
      canvas.requestDraw();
      canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
      await canvas.graphics.resources.pending.whenReady();
      const objects = canvas.getRenderObjects().filter((o) => o.state.visible);
      expect(objects.length).toBe(2);
      /** a hover moving over the groups (atoms): the previous highlight removed, the next one added, marking only */
      let previous: ReturnType<typeof canvas.getLoci> | undefined;
      const hover = async (objectIndex: number, group: number) => {
        const loci = canvas.getLoci({ objectId: objects[objectIndex].id, instanceId: 0, groupId: group });
        if (previous) canvas.mark(previous, MarkerAction.RemoveHighlight);
        canvas.mark(loci, MarkerAction.Highlight);
        previous = loci;
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
      };
      const counts = () => ({ ...stats.drawLists });
      const delta = (a: ReturnType<typeof counts>, b: ReturnType<typeof counts>) => ({
        replayed: b.replayed - a.replayed,
        flushed: b.flushed - a.flushed,
        recorded: b.recorded - a.recorded,
        direct: b.chunksDirect - a.chunksDirect,
      });

      // compiled, then recorded
      for (let g = 0; g < 4; ++g) await hover(0, g);
      for (let g = 0; g < 4; ++g) await hover(0, g);
      const within0 = counts();
      for (let g = 4; g < 12; ++g) await hover(0, g);
      const within = delta(within0, counts());
      const across0 = counts();
      for (let g = 0; g < 8; ++g) await hover(g % 2, g);
      const across = delta(across0, counts());

      /** a hover marking a whole object (every group of it), marking only */
      let previousObject: ReturnType<typeof canvas.getLoci>[] = [];
      const hoverObject = async (objectIndex: number) => {
        for (const loci of previousObject) canvas.mark(loci, MarkerAction.RemoveHighlight);
        previousObject = Array.from({ length: 24 }, (_, g) =>
          canvas.getLoci({ objectId: objects[objectIndex].id, instanceId: 0, groupId: g }),
        );
        for (const loci of previousObject) canvas.mark(loci, MarkerAction.Highlight);
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
      };
      if (previous) canvas.mark(previous, MarkerAction.RemoveHighlight);
      for (let i = 0; i < 4; ++i) await hoverObject(i % 2);
      const whole0 = counts();
      for (let i = 0; i < 8; ++i) await hoverObject(i % 2);
      const whole = delta(whole0, counts());

      // the hovered object's chunk changes every frame and holds few objects: drawn directly in both marking passes (a
      // flush and a replay cost more there), the other passes replayed, nothing recorded again
      expect(within.recorded).toBe(0);
      expect(within.replayed).toBe(16);
      expect(within.direct).toBe(16);
      // across objects and whole objects: the depth pass replays; the mask pass draws the marked objects only, recorded
      // again when they change
      expect(across.replayed).toBeGreaterThanOrEqual(8);
      expect(across.recorded).toBeLessThanOrEqual(8);
      expect(whole.replayed).toBeGreaterThanOrEqual(8);
      expect(whole.recorded).toBeLessThanOrEqual(8);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('draw the marker ranges uploaded as a whole upload of the layer draws them', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      // highlights tint the color of the draw pass (not drawn by the marking pass)
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        marking: { ...canvas.props.marking, enabled: false },
      });
      const passes = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        return Array.from(await canvas.graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
      };
      await frame();
      const plain = await frame();
      const objects = canvas.getRenderObjects().filter((o) => o.state.visible);
      // single atoms highlighted and unhighlighted: the marker layer uploaded by its changed ranges
      let previous: ReturnType<typeof canvas.getLoci> | undefined;
      for (const group of [3, 17, 8, 20, 11]) {
        const loci = canvas.getLoci({ objectId: objects[0].id, instanceId: 0, groupId: group });
        if (previous) canvas.mark(previous, MarkerAction.RemoveHighlight);
        canvas.mark(loci, MarkerAction.Highlight);
        previous = loci;
        await frame();
      }
      const ranged = await frame();
      // the same markers uploaded whole (an update without logged ranges)
      for (const o of objects) ValueCell.update(o.values.tMarker, o.values.tMarker.ref.value);
      const whole = await frame();
      expect(ranged).not.toEqual(plain);
      expect(ranged).toEqual(whole);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it.each([
    ['cartoon', true],
    ['cartoon', false],
    ['spacefill', true],
  ] as const)(
    'hover residue after residue (%s, marking pass %s): no earlier highlight stays',
    async (type, marking) => {
      const plugin = await HeadlessPluginContext.create(
        { webgpu: WebGPU.module! },
        DefaultPluginSpec(),
        { width: Size, height: Size },
        { backend: 'webgpu' },
      );
      try {
        await plugin.init();
        const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        await plugin.builders.structure.representation.addRepresentation(structure, { type });
        const canvas = plugin.canvas3d!;
        canvas.setProps({
          multiSample: { ...canvas.props.multiSample, mode: 'off' },
          marking: { ...canvas.props.marking, enabled: marking },
        });
        const passes = plugin.renderer.canvas3dContext!.passes;
        let time = now();
        const frame = async () => {
          canvas.commit(true);
          canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
          return Array.from(await canvas.graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
        };
        await frame();
        const objects = canvas.getRenderObjects().filter((o) => o.state.visible);
        const loci = (group: number) => canvas.getLoci({ objectId: objects[0].id, instanceId: 0, groupId: group });
        // hovers as the viewer issues them: the previous one cleared (all clear for a moment), the next one marked
        let previous: ReturnType<typeof canvas.getLoci> | undefined;
        for (const group of [3, 17, 8, 20]) {
          if (previous) canvas.mark(previous, MarkerAction.RemoveHighlight);
          canvas.mark(loci(group), MarkerAction.Highlight);
          previous = loci(group);
          await frame();
        }
        const hovered = await frame();
        // the last residue alone, its markers uploaded whole
        canvas.mark(previous!, MarkerAction.RemoveHighlight);
        await frame();
        canvas.mark(previous!, MarkerAction.Highlight);
        for (const o of objects) ValueCell.update(o.values.tMarker, o.values.tMarker.ref.value);
        const alone = await frame();
        expect(hovered).toEqual(alone);
      } finally {
        plugin.dispose();
      }
    },
    60_000,
  );

  it('draw the marker parts of a merged object uploaded as a whole upload of its layer draws them', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { renderMerged: true },
      });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        marking: { ...canvas.props.marking, enabled: false },
      });
      const passes = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        return Array.from(await canvas.graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
      };
      await frame();
      const plain = await frame();
      const merged = canvas.getRenderObjects().find((o) => o.state.visible && isMergedRenderObject(o));
      if (!merged) throw new Error('expected a merged render object');
      // single atoms of both members (an instance each) highlighted and unhighlighted: the changed members' parts written
      let previous: ReturnType<typeof canvas.getLoci> | undefined;
      for (const [instanceId, groupId] of [
        [0, 3],
        [1, 5],
        [0, 8],
        [1, 2],
      ]) {
        const loci = canvas.getLoci({ objectId: merged.id, instanceId, groupId });
        if (previous) canvas.mark(previous, MarkerAction.RemoveHighlight);
        canvas.mark(loci, MarkerAction.Highlight);
        previous = loci;
        await frame();
      }
      const parts = await frame();
      // the same markers uploaded whole: a new segmented array of the same parts
      const image = merged.values.tMarker.ref.value;
      const array = image.array as unknown as SegmentedArray;
      const copy = SegmentedArray.create(
        array.ctor,
        array.parts.map((p) => ({ ...p })),
        array.length,
      );
      ValueCell.update(merged.values.tMarker, { ...image, array: copy as unknown as Uint8Array });
      const whole = await frame();
      expect(parts).not.toEqual(plain);
      expect(parts).toEqual(whole);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('replay without planning only while nothing the plan read changed (state, camera)', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const repr = await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'cartoon' });
      const canvas = plugin.canvas3d!;
      canvas.setProps({ multiSample: { ...canvas.props.multiSample, mode: 'off' } });
      const { stats } = canvas.graphics;
      const passes = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        return Array.from(await canvas.graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
      };
      const reprObject = plugin.state.data.select(repr!.ref)[0].obj!.data.repr;
      const snapshot = canvas.camera.getSnapshot();
      /** still, the spacefill hidden, the camera moved, the spacefill shown again: each frame twice */
      const sequence = async () => {
        const frames: number[][] = [];
        frames.push(await frame(), await frame());
        reprObject.setState({ visible: false });
        frames.push(await frame(), await frame());
        canvas.camera.setState({ position: Vec3.add(Vec3(), snapshot.position, Vec3.create(4, 0, 0)) }, 0);
        frames.push(await frame(), await frame());
        reprObject.setState({ visible: true });
        frames.push(await frame(), await frame());
        canvas.camera.setState(snapshot, 0);
        return frames;
      };

      // compiled, then drawn directly
      await sequence();
      const direct = await sequence();
      // hidden, and shown with the camera moved
      expect(direct[0]).not.toEqual(direct[2]);
      expect(direct[0]).not.toEqual(direct[6]);

      canvas.setProps({ renderer: { ...canvas.props.renderer, drawLists: true } });
      await sequence();
      const unplanned = stats.drawLists.unplanned;
      expect(await sequence()).toEqual(direct);
      // the second frame of each state replays without planning
      expect(stats.drawLists.unplanned).toBeGreaterThan(unplanned);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
