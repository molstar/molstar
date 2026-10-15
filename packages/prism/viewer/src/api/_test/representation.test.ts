/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection, type StateTransformer } from '@molstar/core/state';
import { ClippingStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/clipping';
import { EmissiveStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/emissive';
import { OverpaintStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { TransparencyStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/transparency';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// `representation` (operations.md §2, plan 2.2) on a headless viewer with 1CBS.

const cells = (v: PrismViewer, t: StateTransformer) =>
  v.plugin.state.data.select(StateSelection.Generators.ofTransformer(t));

const components = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0].components;

const reprs = (v: PrismViewer) => components(v).flatMap((c) => c.representations);

const themeNames = (v: PrismViewer) => reprs(v).map((r) => (r.cell.transform.params as any).colorTheme.name as string);

const nameOf = (c: { cell: { obj?: { label?: string } } }) => c.cell.obj?.label ?? '';

/** The label of a component that has a representation of a type, or the first component. */
const polymerLabel = (v: PrismViewer) =>
  nameOf(components(v).find((c) => /polymer/i.test(nameOf(c))) ?? components(v)[0]);

describe('representation', () => {
  it('colors a target with an overpaint layer, and undo removes it', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(cells(v, OverpaintStructureRepresentation3DFromBundle)).toHaveLength(0);

    const r = await v.api.representation.color({ target: 'ligand', color: 'red' });
    expectOk(r);
    expect(r.changed.nodes?.length).toBeGreaterThan(0);
    expect(cells(v, OverpaintStructureRepresentation3DFromBundle).length).toBeGreaterThan(0);

    expectOk(await v.api.history.undo());
    expect(cells(v, OverpaintStructureRepresentation3DFromBundle)).toHaveLength(0);
  });

  it('says when a color is already painted, and refuses a bad color or neither/both of color and theme', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.representation.color({ target: 'ligand', color: '#ff0000' }));
    const again = await v.api.representation.color({ target: 'ligand', color: 'red' });
    expectOk(again);
    expect(again.unchanged).toBe(true);

    expect(await v.api.representation.color({ target: 'ligand', color: 'nonsense' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.representation.color({ target: 'ligand' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.representation.color({ target: 'ligand', color: 'red', theme: 'element-symbol' })).toMatchObject(
      {
        ok: false,
        code: 'invalid',
      },
    );
  });

  it('applies a color theme to all representations, and names the valid themes when wrong', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(themeNames(v).every((n) => n === 'element-symbol')).toBe(false);

    const r = await v.api.representation.color({ target: 'all', theme: 'element-symbol' });
    expectOk(r);
    expect(themeNames(v).length).toBeGreaterThan(0);
    expect(themeNames(v).every((n) => n === 'element-symbol')).toBe(true);

    const again = await v.api.representation.color({ target: 'all', theme: 'element-symbol' });
    expectOk(again);
    expect(again.unchanged).toBe(true);

    const wrong = await v.api.representation.color({ target: 'all', theme: 'no-such-theme' });
    expect(wrong).toMatchObject({ ok: false, code: 'invalid' });
    if (!wrong.ok) expect(wrong.hint).toMatch(/element-symbol/);

    expectOk(await v.api.history.undo());
    expect(themeNames(v).every((n) => n === 'element-symbol')).toBe(false);
  });

  it('sets transparency on a target, resets it, and says when there is nothing to reset', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(await v.api.representation.reset({ layer: 'transparency' })).toMatchObject({ ok: true, unchanged: true });
    // Opaque where nothing is transparent: nothing to do.
    expect(await v.api.representation.transparency({ target: 'polymer', alpha: 0 })).toMatchObject({
      ok: true,
      unchanged: true,
    });

    const r = await v.api.representation.transparency({ target: 'polymer', alpha: 0.5 });
    expectOk(r);
    expect(r.unchanged).toBeUndefined();
    expect(cells(v, TransparencyStructureRepresentation3DFromBundle).length).toBeGreaterThan(0);
    const same = await v.api.representation.transparency({ target: 'polymer', alpha: 0.5 });
    expect(same).toMatchObject({ ok: true, unchanged: true });
    expect(await v.api.representation.transparency({ target: 'polymer', alpha: 2 })).toMatchObject({
      ok: false,
      code: 'invalid',
    });

    expectOk(await v.api.representation.reset({ layer: 'transparency' }));
    expect(cells(v, TransparencyStructureRepresentation3DFromBundle)).toHaveLength(0);
  });

  it('creates the emissive and clipping layers, and resets each', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.representation.emissive({ target: 'ligand', strength: 0.4 }));
    expect(cells(v, EmissiveStructureRepresentation3DFromBundle).length).toBeGreaterThan(0);
    expectOk(await v.api.representation.clipping({ target: 'ligand' }));
    expect(cells(v, ClippingStructureRepresentation3DFromBundle).length).toBeGreaterThan(0);

    const wrong = await v.api.representation.clipping({ target: 'ligand', groups: ['seven'] as never });
    expect(wrong).toMatchObject({ ok: false, code: 'invalid' });

    expectOk(await v.api.representation.reset({ layer: 'emissive' }));
    expect(cells(v, EmissiveStructureRepresentation3DFromBundle)).toHaveLength(0);
    expect(cells(v, ClippingStructureRepresentation3DFromBundle).length).toBeGreaterThan(0);
    expectOk(await v.api.representation.reset({ layer: 'clipping' }));
    expect(cells(v, ClippingStructureRepresentation3DFromBundle)).toHaveLength(0);
  });

  it('says a target that selects nothing selects nothing', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.representation.color({ target: 'A:9999', color: 'red' });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.message).toMatch(/selects nothing/);
    expect(cells(v, OverpaintStructureRepresentation3DFromBundle)).toHaveLength(0);
  });

  it('sets a representation type on a component, once', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const component = polymerLabel(v);
    const shown = () => reprs(v).filter((x) => (x.cell.transform.params as any).type.name === 'ball-and-stick').length;
    const before = shown();

    const r = await v.api.representation.set({ component, type: 'ball-and-stick' });
    expectOk(r);
    expect(shown()).toBe(before + 1);
    expect(r.changed.nodes?.length).toBe(1);

    // The same again: nothing changes, so no history step.
    const again = await v.api.representation.set({ component, type: 'ball-and-stick' });
    expectOk(again);
    expect(again.unchanged).toBe(true);

    // Another size: the params of the existing one are updated, not another added.
    const wider = await v.api.representation.set({ component, type: 'ball-and-stick', params: { sizeFactor: 0.3 } });
    expectOk(wider);
    expect(shown()).toBe(before + 1);
    const bs = components(v)
      .find((c) => nameOf(c) === component)!
      .representations.find((x) => (x.cell.transform.params as any).type.name === 'ball-and-stick')!;
    expect((bs.cell.transform.params as any).type.params.sizeFactor).toBe(0.3);
    const same = await v.api.representation.set({ component, type: 'ball-and-stick', params: { sizeFactor: 0.3 } });
    expect(same).toMatchObject({ ok: true, unchanged: true });

    // Two steps (the add, the update): undo twice.
    expectOk(await v.api.history.undo());
    expect(shown()).toBe(before + 1);
    expectOk(await v.api.history.undo());
    expect(shown()).toBe(before);
  });

  it('refuses a wrong type, a wrong param, and a wrong component, with the valid ones', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const component = polymerLabel(v);
    const type = await v.api.representation.set({ component, type: 'nonsense' });
    expect(type).toMatchObject({ ok: false, code: 'invalid' });
    if (!type.ok) expect(type.hint ?? '').toMatch(/ball-and-stick/);

    const param = await v.api.representation.set({ component, type: 'ball-and-stick', params: { nope: 1 } });
    expect(param).toMatchObject({ ok: false, code: 'invalid' });
    if (!param.ok) expect(param.hint ?? '').toMatch(/sizeFactor/);

    const missing = await v.api.representation.set({ component: 'no such component', type: 'cartoon' });
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
    if (!missing.ok) expect(missing.hint ?? '').toContain(component);
  });

  it('removes a type, or all, from a component', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const polymer = components(v).find((c) => c.representations.length > 0)!;
    const component = nameOf(polymer);
    const type = (polymer.representations[0].cell.transform.params as any).type.name as string;
    const count = reprs(v).length;

    expect(await v.api.representation.remove({ component, type: 'spacefill' })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expectOk(await v.api.representation.remove({ component, type }));
    expect(reprs(v).length).toBeLessThan(count);
    expect(
      components(v).find((c) => c.cell.transform.ref === polymer.cell.transform.ref)?.representations ?? [],
    ).toHaveLength(0);

    expectOk(await v.api.history.undo());
    expect(reprs(v)).toHaveLength(count);
  });

  it('sets the size factor and the material on representations', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const component = polymerLabel(v);
    expectOk(await v.api.representation.set({ component, type: 'ball-and-stick' }));
    const bs = () => reprs(v).find((x) => (x.cell.transform.params as any).type.name === 'ball-and-stick')!;
    const base = (bs().cell.transform.params as any).type.params.sizeFactor as number;

    expectOk(await v.api.representation.size({ factor: 2 }));
    expect((bs().cell.transform.params as any).type.params.sizeFactor).toBeCloseTo(base * 2);
    expect(await v.api.representation.size({ factor: 2 })).toMatchObject({ ok: true, unchanged: true });
    expect(await v.api.representation.size({})).toMatchObject({ ok: false, code: 'invalid' });

    expectOk(await v.api.representation.material({ material: { metalness: 0.5, roughness: 0.3, bumpiness: 0 } }));
    expect((bs().cell.transform.params as any).type.params.material).toMatchObject({ metalness: 0.5, roughness: 0.3 });
    expect(
      await v.api.representation.material({ material: { metalness: 0.5, roughness: 0.3, bumpiness: 0 } }),
    ).toMatchObject({ ok: true, unchanged: true });
  });
});
