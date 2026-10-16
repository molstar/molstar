/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection, StateTransform } from '@molstar/core/state';
import { setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { beforeAll, describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import { measurementLocis, readMeasurements } from '../../../sources/measurements.js';
import { createTargetResolver } from '../../../sources/targets.js';
import type { PrismViewer } from '../../../viewer.js';
import { MeasurementGroupTag } from '@molstar/plugin/state/manager/structure/measurement';
import { MeasurementsFeature } from '../index.js';
import { shortLabel } from '../source.js';

// The `measure` area (operations.md §2) on a headless viewer with 1CBS: the measurement manager's cells, as results.

const steps = (v: PrismViewer) => v.api.history.steps.value.length;
const manager = (v: PrismViewer) => v.plugin.managers.structure.measurement;

// Measurement texts are drawn into a font atlas, which Node has no canvas for: a stand-in context that draws nothing.
beforeAll(() => {
  const context = {
    measureText: () => ({ width: 8 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});

async function setup() {
  return headlessViewer({ files: [cbs()], extensions: [MeasurementsFeature] });
}

const Å = 'Å';
const A = (n: number) => `A:${n}@CA`;

describe('measure', () => {
  it('measures a distance between two targets, as a step Undo removes', async () => {
    const v = await setup();
    const before = steps(v);
    const r = await v.api.measure!.distance({ a: A(10), b: A(20) });
    expectOk(r);
    expect(r.value?.kind).toBe('distance');
    expect(r.value?.value).toMatch(new RegExp(`^\\d+\\.\\d\\d ${Å}$`));
    expect(r.value?.atoms).toHaveLength(2);
    expect(r.value?.label).toContain(' – ');
    expect(r.did).toContain(r.value!.value);
    expect(steps(v)).toBe(before + 1);
    expect(manager(v).state.distances).toHaveLength(1);
    expect(manager(v).state.distances[0].transform.ref).toBe(r.value?.ref);

    // The cells are under the measurements group (tagged), as the classic UI makes them.
    const state = v.plugin.state.data;
    const group = StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, MeasurementGroupTag);
    expect(group).toBeDefined();
    const selection = manager(v).state.distances[0].transform.parent;
    expect(state.tree.transforms.get(selection)?.parent).toBe(group);

    expectOk(await v.api.history.undo());
    expect(manager(v).state.distances).toHaveLength(0);
    expect(StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, MeasurementGroupTag)).toBeUndefined();
  });

  it('says it is already measured, unless a copy is asked for', async () => {
    const v = await setup();
    expectOk(await v.api.measure!.distance({ a: A(10), b: A(20) }));
    const before = steps(v);
    const again = await v.api.measure!.distance({ a: A(10), b: A(20) });
    expectOk(again);
    expect(again.unchanged).toBe(true);
    expect(again.value?.kind).toBe('distance');
    expect(steps(v)).toBe(before);
    expect(manager(v).state.distances).toHaveLength(1);

    // The order is the identity: B to A is another.
    expectOk(await v.api.measure!.distance({ a: A(20), b: A(10) }));
    expect(manager(v).state.distances).toHaveLength(2);

    const copy = await v.api.measure!.distance({ a: A(10), b: A(20) }, { copy: true });
    expectOk(copy);
    expect(copy.unchanged).toBeUndefined();
    expect(manager(v).state.distances).toHaveLength(3);
  });

  it('measures a residue as its center, and angle, dihedral and label points', async () => {
    const v = await setup();
    const residues = await v.api.measure!.distance({ a: 'A:10', b: 'ligand' });
    expectOk(residues);
    expect(residues.value?.value).toMatch(/\d/);

    const angle = await v.api.measure!.angle({ a: A(10), b: A(11), c: A(12) });
    expectOk(angle);
    expect(angle.value?.value).toMatch(/^\d+\.\d\d°$/);
    expect(angle.value?.atoms).toHaveLength(3);

    const dihedral = await v.api.measure!.dihedral({ a: A(10), b: A(11), c: A(12), d: A(13) });
    expectOk(dihedral);
    expect(dihedral.value?.value).toMatch(/^-?\d+\.\d\d°$/);
    expect(dihedral.value?.atoms).toHaveLength(4);

    const label = await v.api.measure!.label({ a: A(10), text: 'Hello' });
    expectOk(label);
    expect(label.value?.value).toBe('Hello');
    // The same atom with another text is another label.
    const other = await v.api.measure!.label({ a: A(10), text: 'World' });
    expectOk(other);
    expect(other.unchanged).toBeUndefined();
    const same = await v.api.measure!.label({ a: A(10), text: 'World' });
    expect(same).toMatchObject({ ok: true, unchanged: true });
    expect(manager(v).state.labels).toHaveLength(2);
  });

  it('makes an orientation and a plane of a target, or of the selection', async () => {
    const v = await setup();
    const orientation = await v.api.measure!.orientation({ target: 'ligand' });
    expectOk(orientation);
    expect(orientation.value?.kind).toBe('orientation');
    expect(manager(v).state.orientations).toHaveLength(1);
    expect(orientation.unchanged).toBeUndefined();
    expect(await v.api.measure!.orientation({ target: 'ligand' })).toMatchObject({ ok: true, unchanged: true });

    const plane = await v.api.measure!.plane({ target: 'A:10-20' });
    expectOk(plane);
    expect(manager(v).state.planes).toHaveLength(1);

    // Nothing selected, no target: it asks for one.
    const none = await v.api.measure!.plane({});
    expect(none).toMatchObject({ ok: false, code: 'invalid' });
    expectOk(await v.api.select.set({ target: 'A:30-35' }));
    const selected = await v.api.measure!.plane({});
    expectOk(selected);
    expect(manager(v).state.planes).toHaveLength(2);
    expect(selected.value?.label).toMatch(/Residues/);
  });

  it('takes the latest picks when no targets are given', async () => {
    const v = await setup();
    const none = await v.api.measure!.distance({});
    expect(none).toMatchObject({ ok: false, code: 'invalid' });
    if (!none.ok) expect(none.hint).toMatch(/pick/);

    const picks = v.plugin.managers.structure.selection;
    const loci = async (t: string) => {
      const r = await createTargetResolver(v.plugin).resolve(t);
      if (!r.ok) throw new Error(r.message);
      return r.parts[0].loci;
    };
    picks.fromLoci('add', await loci(A(10)), false);
    expect(await v.api.measure!.distance({})).toMatchObject({ ok: false, code: 'invalid' });
    picks.fromLoci('add', await loci(A(20)), false);
    // The newest is the first point.
    const r = await v.api.measure!.distance({});
    expectOk(r);
    const explicit = await v.api.measure!.distance({ a: A(20), b: A(10) });
    expect(explicit).toMatchObject({ ok: true, unchanged: true });
  });

  it('fails clearly when a target selects nothing, or the points are not all given', async () => {
    const v = await setup();
    const before = steps(v);
    const nothing = await v.api.measure!.distance({ a: A(10), b: 'A:9999@CA' });
    expect(nothing).toMatchObject({ ok: false, code: 'invalid' });
    if (!nothing.ok) expect(nothing.message).toMatch(/second point selects nothing/);

    const some = await v.api.measure!.angle({ a: A(10), b: A(11) });
    expect(some).toMatchObject({ ok: false, code: 'invalid' });
    if (!some.ok) expect(some.message).toMatch(/needs 3 points, 2 given/);
    expect(manager(v).state.distances).toHaveLength(0);
    expect(steps(v)).toBe(before);

    const empty = await headlessViewer({ extensions: [MeasurementsFeature] });
    expect(await empty.api.measure!.distance({ a: A(10), b: A(20) })).toMatchObject({ ok: false });
  });

  it("reads a measurement's loci by its ref, as the list does, without describing the others", async () => {
    const v = await setup();
    const d = await v.api.measure!.distance({ a: A(10), b: A(20) });
    expectOk(d);
    expectOk(await v.api.measure!.label({ a: A(10), text: 'x' }));
    const listed = readMeasurements(v.plugin).find((m) => m.ref === d.value!.ref)!;
    expect(listed.locis).toHaveLength(2);
    expect(measurementLocis(v.plugin, d.value!.ref)).toEqual(listed.locis);
    expect(measurementLocis(v.plugin, 'nope')).toEqual([]);
  });

  it('describes the measurements, shows or hides and removes them', async () => {
    const v = await setup();
    const none = await v.api.measure!.list();
    expectOk(none);
    expect(none.value).toEqual([]);
    expect(await v.api.measure!.remove({})).toMatchObject({ ok: false, code: 'not-found' });

    const d = await v.api.measure!.distance({ a: A(10), b: A(20) });
    expectOk(d);
    const l = await v.api.measure!.label({ a: A(10), text: 'x' });
    expectOk(l);
    const got = await v.api.measure!.list();
    expectOk(got);
    expect(got.value?.map((m) => m.kind)).toEqual(['label', 'distance']);
    expect(got.value?.[1]).toMatchObject({ ref: d.value?.ref, hidden: false });

    const hidden = await v.api.measure!.visibility({ measurement: d.value!.ref, visible: false });
    expectOk(hidden);
    expect((await v.api.measure!.list()).ok && (await v.api.measure!.list())).toMatchObject({
      value: [{ hidden: false }, { hidden: true }],
    });
    expect(await v.api.measure!.visibility({ measurement: d.value!.ref, visible: false })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expectOk(await v.api.measure!.visibility({ measurement: d.value!.ref, visible: true }));

    // Hidden through its representation cell alone (the State tree's eye): showing it shows it.
    expectOk(await v.api.state!.visibility({ node: d.value!.ref, visible: false }));
    const viaRepresentation = await v.api.measure!.list();
    expectOk(viaRepresentation);
    expect(viaRepresentation.value?.[1]).toMatchObject({ hidden: true });
    expectOk(await v.api.measure!.visibility({ measurement: d.value!.ref, visible: true }));
    const shownAgain = await v.api.measure!.list();
    expectOk(shownAgain);
    expect(shownAgain.value?.[1]).toMatchObject({ hidden: false });

    // Which one? Several exist.
    expect(await v.api.measure!.remove({})).toMatchObject({ ok: false, code: 'ambiguous' });
    expect(await v.api.measure!.remove({ measurement: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });

    const before = steps(v);
    const removed = await v.api.measure!.remove({ measurement: d.value!.ref });
    expectOk(removed);
    expect(manager(v).state.distances).toHaveLength(0);
    expect(manager(v).state.labels).toHaveLength(1);
    expect(steps(v)).toBe(before + 1);
    expectOk(await v.api.history.undo());
    expect(manager(v).state.distances).toHaveLength(1);

    // By its label, too.
    const byLabel = await v.api.measure!.remove({ measurement: d.value!.label });
    expectOk(byLabel);
    expect(manager(v).state.distances).toHaveLength(0);
  });

  it('asks before removing, from the agent', async () => {
    const v = await setup();
    expect(v.api.measure!.remove.definition.approval).toBe('destructive');
  });

  it('sets the options, updating the measurements made', async () => {
    const v = await setup();
    const d = await v.api.measure!.distance({ a: A(10), b: A(20) });
    expectOk(d);
    const r = await v.api.measure!.options({ distanceUnitLabel: 'nm', textColor: '#ff8800' });
    expectOk(r);
    expect(manager(v).state.options.distanceUnitLabel).toBe('nm');
    expect(manager(v).state.options.textColor).toBe(0xff8800);
    const got = await v.api.measure!.list();
    expectOk(got);
    expect(got.value?.[0].value).toMatch(/ nm$/);

    expect(await v.api.measure!.options({ distanceUnitLabel: 'nm', textColor: '#ff8800' })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    // Unset values keep what is.
    expect(await v.api.measure!.options({})).toMatchObject({ ok: true, unchanged: true });
    expect(await v.api.measure!.options({ textColor: 'red' })).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('the Measure section model', () => {
  it('measures from the picks and follows the measurements made', async () => {
    const v = await setup();
    const model = v.api.measure!.model(v.api);
    // A view binds to what it shows; the nodes follow the measurements while it does.
    const watching = [model.points.subscribe(), model.added.subscribe()];
    expect(model.points.value).toHaveLength(0);
    expect(model.kinds.map((k) => k.id)).toEqual(['distance', 'angle', 'dihedral', 'label', 'orientation', 'plane']);

    const picks = v.plugin.managers.structure.selection;
    for (const t of [A(10), A(20), A(30)]) {
      const r = await createTargetResolver(v.plugin).resolve(t);
      if (!r.ok) throw new Error(r.message);
      picks.fromLoci('add', r.parts[0].loci, false);
    }
    expect(model.points.value).toHaveLength(3);
    const ids = model.points.value.map((p) => p.id);
    model.movePoint(ids[0], 1);
    expect(model.points.value.map((p) => p.id)).toEqual([ids[1], ids[0], ids[2]]);
    model.movePointTo(ids[2], 0);
    expect(model.points.value[0].id).toBe(ids[2]);
    model.removePoint(ids[2]);
    expect(model.points.value).toHaveLength(2);

    model.add('distance');
    await until(() => model.added.value.length === 1);
    const node = model.added.value[0];
    expect(node.kind).toBe('measurement');
    expect(node.description.value).toMatch(/\d/);
    node.toggleVisibility();
    await until(() => node.isHidden.value);
    node.remove();
    await until(() => model.added.value.length === 0);
    for (const w of watching) w.unsubscribe();
  });

  it('listens to the plugin only while something watches, and reads once per frame', async () => {
    const v = await setup();
    const model = v.api.measure!.model(v.api);
    const events = v.plugin.state.data.events as unknown as Record<string, { observers: unknown[] }>;
    const history = (v.plugin.managers.structure.selection.events as any).additionsHistoryUpdated;
    const observers = () => [events.changed.observers.length, history.observers.length];
    const idle = observers();
    // Reading `value` is a read, not a subscription.
    expect(model.added.value).toHaveLength(0);
    expect(model.points.value).toHaveLength(0);
    expect(observers()).toEqual(idle);

    const watching = [model.added.subscribe(), model.points.subscribe()];
    expect(observers()).not.toEqual(idle);
    expectOk(await v.api.measure!.distance({ a: A(10), b: A(20) }));
    // `value` is current at once, without waiting for the frame.
    expect(model.added.value).toHaveLength(1);
    for (const w of watching) w.unsubscribe();
    expect(observers()).toEqual(idle);
  });

  it('shows the short form of a point and the full label for the tooltip', async () => {
    const v = await setup();
    const model = v.api.measure!.model(v.api);
    const r = await createTargetResolver(v.plugin).resolve('A:10@CA');
    if (!r.ok) throw new Error(r.message);
    v.plugin.managers.structure.selection.fromLoci('add', r.parts[0].loci, false);
    const [point] = model.points.value;
    expect(point.label).toContain('|');
    expect(point.short).toMatch(/^.+ · .+$/);
    expect(point.short!.split(' · ')[0]).toBe(point.label.split('|')[0].trim());
    expect(point.short!.endsWith(point.label.split('|').pop()!.trim())).toBe(true);
    expect(shortLabel('CYS 442 SG')).toBe('CYS 442 SG');
    expect(shortLabel('REA 200 | B [auth A] | Instance ASM-1 | Model 1 | 1CBS')).toBe('REA 200 · 1CBS');
  });

  it('keeps the options form while an edit is on its way, and sends one step for a drag', async () => {
    const v = await setup();
    const model = v.api.measure!.model(v.api);
    const steps0 = steps(v);
    // Several changes of one drag: one step.
    model.options.set('textColor', 0x112233 as never);
    model.options.set('textColor', 0x223344 as never);
    model.options.set('textColor', 0x334455 as never);
    // Mid-drag, an older value arriving does not pull the form back.
    model.loadOptions!({ distanceUnitLabel: 'Å', textColor: 0xaaaaaa });
    expect((model.options.values.value as any).textColor).toBe(0x334455);
    await until(() => manager(v).state.options.textColor === 0x334455);
    expect(steps(v)).toBe(steps0 + 1);
    // Once it landed, the form follows again.
    await until(() => {
      model.loadOptions!({ distanceUnitLabel: 'nm', textColor: 0x334455 });
      return (model.options.values.value as any).distanceUnitLabel === 'nm';
    });
  });

  it('restores the options with Undo (they are custom state of the checkpoint)', async () => {
    const v = await setup();
    expectOk(await v.api.measure!.options({ distanceUnitLabel: 'nm', textColor: '#ff8800' }));
    expect(manager(v).state.options.distanceUnitLabel).toBe('nm');
    expectOk(await v.api.history.undo());
    expect(manager(v).state.options.distanceUnitLabel).toBe('\u212B');
    expect(manager(v).state.options.textColor).not.toBe(0xff8800);
  });
});

describe('the feature', () => {
  it('brings commands for its operations and one that opens the Measure section', async () => {
    const v = await setup();
    const ids = ['distance', 'angle', 'dihedral', 'label', 'orientation', 'plane', 'remove', 'options', 'open'].map(
      (n) => `measure.${n}`,
    );
    for (const id of ids) expect(v.managers.commands.get(id), id).toBeDefined();
    // Every operation has a command (operations.md §4).
    const defined = v.api
      .area('measure')!
      .operations()
      .map((o) => o.definition.id);
    for (const id of defined) expect(v.managers.commands.get(id), id).toBeDefined();
  });

  it('offers Measure… only where there is a Selection tab to open', async () => {
    const context = { empty: false };
    // Headless has no tabs.
    const v = await setup();
    expect(v.managers.commands.get('measure.open')!.when?.(context)).toMatch(/no Selection panel/);
    v.managers.viewState.setDock('right', { tabs: ['selection'] });
    expect(v.managers.commands.get('measure.open')!.when?.(context)).toBe(true);
  });

  it('puts the Measure section in the Selection tab', async () => {
    const v = await setup();
    const ui = await MeasurementsFeature.ui!({ api: v.api, sources: undefined as never, managers: v.managers });
    expect(ui.sections?.map((s) => [s.id, s.slot, s.title, s.icon, s.order])).toEqual([
      ['measurements.section.measure', 'selection', 'Measure', 'ruler', 2],
    ]);
  });
});
