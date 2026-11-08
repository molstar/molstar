/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { Mat4, Vec3, type Vec4 } from '@molstar/core/math/linear-algebra';
import { Binding } from '@molstar/core/util/binding';
import { Color } from '@molstar/core/util/color';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import { BehaviorSubject, Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import type { Canvas3D } from '../../../canvas3d/canvas3d.js';
import { Mesh } from '../../../geo/geometry/mesh/mesh.js';
import { LineStride, type Lines } from '../../../geo/geometry/lines/lines.js';
import { LinesBuilder } from '../../../geo/geometry/lines/lines-builder.js';
import { PointsBuilder } from '../../../geo/geometry/points/points-builder.js';
import type { GraphicsContext } from '../../../gl/context.js';
import type { Representation } from '../../../repr/representation.js';
import { ShapeRepresentation } from '../../../repr/shape/representation.js';
import {
  DefaultOccludedAlpha,
  type GizmoBinding,
  GizmoController,
  type GizmoItem,
  type GizmoPart,
  type GuideFrame,
  guideLines,
} from '../controller.js';
import { transformHandles } from '../transform.js';

// The gizmo controller's lifecycle on a fake canvas: hover takes the trackball's drags and gives them back, a drag
// previews and commits once, and every other end (Esc, a change under it, a commit not taken, nothing changed)
// cancels and restores.

type Value = { readonly x: number };

const graphics = { capabilities: { flatInterpolation: true } } as unknown as GraphicsContext;
const drags = {
  dragRotate: Binding([]),
  dragPan: Binding([]),
} as unknown as Canvas3D['attribs']['trackball']['bindings'];

function fakeCanvas() {
  const added: Representation.Any[] = [];
  const attribs = { trackball: { bindings: drags } };
  const input = {
    pixelRatio: 1,
    // the element's size: the viewport's
    width: 100,
    height: 100,
    drag: new Subject<{ isStart: boolean; x: number; y: number; dx: number; dy: number }>(),
    interactionEnd: new Subject<void>(),
    keyDown: new Subject<{ key: string }>(),
  };
  const canvas = {
    graphics,
    add: vi.fn((r: Representation.Any) => void added.push(r)),
    remove: vi.fn((r: Representation.Any) => void added.splice(added.indexOf(r), 1)),
    update: vi.fn(),
    requestDraw: vi.fn(),
    attribs,
    props: { trackball: { rotateSpeed: 5 } },
    setAttribs: vi.fn((a: { trackball: { bindings: typeof drags } }) => {
      attribs.trackball.bindings = a.trackball.bindings;
    }),
    interaction: { hover: new Subject<{ current: Representation.Loci }>() },
    input,
    didDraw: new Subject<number>(),
    camera: {
      state: {
        position: Vec3.create(0, 0, 10),
        target: Vec3.create(0, 0, 0),
        up: Vec3.create(0, 1, 0),
        fov: Math.PI / 4,
        mode: 'perspective',
      },
      viewport: { x: 0, y: 0, width: 100, height: 100 },
      project: (out: Vec4, p: Vec3) => {
        out[0] = 50 + p[0] * 10;
        out[1] = 50 - p[1] * 10;
        return out;
      },
    },
  };
  return { canvas: canvas as unknown as Canvas3D, fake: canvas, added };
}

function fakeBinding(commitResult = true) {
  const items = new BehaviorSubject<readonly GizmoItem<Value>[]>([{ key: 'a', value: { x: 0 }, label: 'A' }]);
  const invalidated = new Subject<string>();
  const binding = {
    items,
    invalidated,
    // the handles alone: one representation a gizmo (the hit area's own test adds it)
    layout: () => ({
      ...headsOnly(transformHandles(1, { planes: false, rotate: false })),
      frame: Mat4.identity(),
    }),
    sameGeometry: () => true,
    equal: (a: Value, b: Value) => a.x === b.x,
    apply: (from: Value, _h: unknown, _f: unknown, dx: number) => ({ x: from.x + dx }),
    begin: vi.fn(),
    preview: vi.fn(async () => {}),
    commit: vi.fn(async () => commitResult),
    cancel: vi.fn(),
  } satisfies GizmoBinding<Value>;
  return binding;
}

const settle = () => new Promise((r) => setTimeout(r, 0));

/** Transform handles as one part: their heads and squares (`handles`), all their groups */
function headsOnly<T extends { parts: readonly GizmoPart[]; handles: Mesh }>(h: T) {
  return { ...h, parts: [{ role: 'handle' as const, geometry: h.handles }] };
}

/** The fake canvas' handle representation (the gizmo's only one: no body) and a loci on its group `g` */
function hoverOn(fake: ReturnType<typeof fakeCanvas>, g: number, at = 0) {
  const repr = fake.added[at];
  const shape = (repr.getAllLoci()[0] as Shape.Loci).shape;
  const loci = ShapeGroup.Loci(shape, [{ ids: OrderedSet.ofSingleton(g), instance: 0 }]);
  fake.fake.interaction.hover.next({ current: { repr, loci } });
}

async function drag(fake: ReturnType<typeof fakeCanvas>, steps: number[]) {
  fake.fake.input.drag.next({ isStart: true, x: 50, y: 50, dx: 0, dy: 0 });
  for (const x of steps) fake.fake.input.drag.next({ isStart: false, x: 50 + x, y: 50, dx: 0, dy: 0 });
  await settle();
}

describe('GizmoController', () => {
  it('takes the trackball drags over a handle and gives them back whole', async () => {
    const fake = fakeCanvas();
    const c = new GizmoController(fake.canvas, fakeBinding());
    await settle();
    expect(fake.added).toHaveLength(1);
    hoverOn(fake, 1);
    expect(fake.fake.attribs.trackball.bindings.dragRotate).toBe(Binding.Empty);
    fake.fake.interaction.hover.next({ current: { repr: undefined, loci: { kind: 'empty-loci' } } as never });
    expect(fake.fake.attribs.trackball.bindings).toEqual(drags);
    c.dispose();
    expect(fake.added).toHaveLength(0);
  });

  it('previews a drag and commits it once, as the last move left it', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    hoverOn(fake, 1);
    await drag(fake, [5, 10, 20]);
    expect(binding.begin).toHaveBeenCalledWith('a', { x: 0 });
    expect(binding.preview).toHaveBeenLastCalledWith('a', { x: 0 }, { x: 20 });
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(binding.commit).toHaveBeenCalledTimes(1);
    expect(binding.commit).toHaveBeenCalledWith('a', { x: 0 }, { x: 20 });
    expect(binding.cancel).not.toHaveBeenCalled();
    expect(c.dragging.value).toBeUndefined();
    c.dispose();
  });

  it('measures a move from where the camera zoomed in a frame made there; commits and cancels from the press', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const calls: { from: number; dx: number; perPixel: number }[] = [];
    const apply = binding.apply;
    binding.apply = (from, h, frame, dx, dy) => {
      calls.push({ from: from.x, dx, perPixel: frame.perPixel });
      return apply(from, h, frame, dx, dy);
    };
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    hoverOn(fake, 1);
    await drag(fake, [10]);
    const before = calls[calls.length - 1].perPixel;
    // zoomed in to half the distance
    Vec3.set(fake.fake.camera.state.position, 0, 0, 5);
    fake.fake.input.drag.next({ isStart: false, x: 65, y: 50, dx: 0, dy: 0 });
    await settle();
    const after = calls[calls.length - 1];
    expect(after.from).toBe(10);
    expect(after.dx).toBe(5);
    expect(after.perPixel).toBeCloseTo(before / 2, 6);
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(binding.commit).toHaveBeenCalledWith('a', { x: 0 }, { x: 15 });
    c.dispose();
  });

  it('cancels on Esc, on a change under the drag, on a commit not taken, and when nothing changed', async () => {
    for (const end of ['esc', 'invalidated', 'refused', 'unchanged'] as const) {
      const fake = fakeCanvas();
      const binding = fakeBinding(end !== 'refused');
      const c = new GizmoController(fake.canvas, binding);
      await settle();
      hoverOn(fake, 1);
      await drag(fake, end === 'unchanged' ? [0] : [10]);
      if (end === 'esc') fake.fake.input.keyDown.next({ key: 'Escape' });
      if (end === 'invalidated') binding.invalidated.next('a');
      fake.fake.input.interactionEnd.next();
      await settle();
      expect(binding.cancel, end).toHaveBeenCalledWith('a', { x: 0 });
      if (end !== 'refused') expect(binding.commit, end).not.toHaveBeenCalled();
      expect(c.dragging.value, end).toBeUndefined();
      // the hover left: the bindings back whole
      fake.fake.interaction.hover.next({ current: { repr: undefined, loci: { kind: 'empty-loci' } } as never });
      expect(fake.fake.attribs.trackball.bindings, end).toEqual(drags);
      c.dispose();
    }
  });

  it('follows its items: a new one drawn, a gone one removed', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    binding.items.next([...binding.items.value, { key: 'b', value: { x: 1 }, label: 'B' }]);
    await settle();
    expect(fake.added).toHaveLength(2);
    binding.items.next([]);
    await settle();
    expect(fake.added).toHaveLength(0);
    expect(c.owns(undefined)).toBe(false);
    c.dispose();
  });

  it('removes the gizmo of an item its binding no longer lays out', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const layout = binding.layout;
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    expect(fake.added).toHaveLength(1);
    // not shown at its new value
    binding.layout = (() => undefined) as unknown as typeof layout;
    binding.items.next([{ key: 'a', value: { x: 1 }, label: 'A' }]);
    await settle();
    expect(fake.added).toHaveLength(0);
    binding.layout = layout;
    binding.items.next([{ key: 'a', value: { x: 2 }, label: 'A' }]);
    await settle();
    expect(fake.added).toHaveLength(1);
    c.dispose();
  });

  it('shares the trackball drags with another controller on the canvas: back whole when the last lets go', async () => {
    const fake = fakeCanvas();
    const a = new GizmoController(fake.canvas, fakeBinding());
    const b = new GizmoController(fake.canvas, fakeBinding());
    await settle();
    expect(fake.added).toHaveLength(2);
    // straight from b's handle to a's (one hover each), then off both
    hoverOn(fake, 1, 1);
    hoverOn(fake, 1, 0);
    expect(fake.fake.attribs.trackball.bindings.dragRotate).toBe(Binding.Empty);
    fake.fake.interaction.hover.next({ current: { repr: undefined, loci: { kind: 'empty-loci' } } as never });
    expect(fake.fake.attribs.trackball.bindings).toEqual(drags);
    a.dispose();
    b.dispose();
  });

  it('cancels a commit that fails, and starts no drag while a commit is in flight', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    let fail: (e: Error) => void = () => {};
    binding.commit.mockImplementationOnce(() => new Promise<boolean>((_, reject) => (fail = reject)));
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    hoverOn(fake, 1);
    await drag(fake, [10]);
    fake.fake.input.interactionEnd.next();
    await settle();
    // pressed again while the commit runs: no drag
    await drag(fake, [5]);
    expect(binding.begin).toHaveBeenCalledTimes(1);
    fail(new Error('no'));
    await settle();
    expect(binding.cancel).toHaveBeenCalledWith('a', { x: 0 });
    expect(warn).toHaveBeenCalled();
    c.dispose();
  });

  it("starts no drag while a cancel's revert waits on a preview in flight", async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    let finish: () => void = () => {};
    binding.preview.mockImplementationOnce(() => new Promise<void>((resolve) => (finish = resolve)));
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    hoverOn(fake, 1);
    await drag(fake, [10]);
    // cancelled while the preview runs: the revert waits on it
    fake.fake.input.keyDown.next({ key: 'Escape' });
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(binding.cancel).not.toHaveBeenCalled();
    // pressed again before the revert: no drag (it would start from the preview it reverts)
    await drag(fake, [5]);
    expect(binding.begin).toHaveBeenCalledTimes(1);
    finish();
    await settle();
    expect(binding.cancel).toHaveBeenCalledWith('a', { x: 0 });
    // after it: a drag from the reverted value
    hoverOn(fake, 1);
    await drag(fake, [5]);
    expect(binding.begin).toHaveBeenCalledTimes(2);
    c.dispose();
  });

  it('cancels the drag of an item that goes away', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const c = new GizmoController(fake.canvas, binding);
    await settle();
    hoverOn(fake, 1);
    await drag(fake, [10]);
    binding.items.next([]);
    await settle();
    expect(binding.cancel).toHaveBeenCalledWith('a', { x: 0 });
    expect(c.dragging.value).toBeUndefined();
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(binding.commit).not.toHaveBeenCalled();
    c.dispose();
  });

  it('adds its gizmos to the overlay: the body x-ray, the handles on top, at the size the layout says', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const c = new GizmoController(fake.canvas, {
      ...binding,
      layout: () => {
        const h = transformHandles(1);
        return { ...h, parts: [...h.parts, { role: 'body', geometry: h.handles }], frame: Mat4.identity() };
      },
    });
    await settle();
    const [handles] = fake.added;
    const body = fake.added.find((r) => r.state.pickable === false)!;
    expect(body.state.pickable).toBe(false);
    expect(c.owns(body)).toBe(true);
    expect(handles.state.layer).toBe('overlay');
    expect(body.state.layer).toBe('overlay');
    expect(handles.props.occludedAlpha).toBe(DefaultOccludedAlpha.handles);
    expect(body.props.occludedAlpha).toBe(DefaultOccludedAlpha.body);
    expect(handles.props.screenSize).toBe(0);
    c.dispose();

    const sized = fakeCanvas();
    const d = new GizmoController(sized.canvas, {
      ...binding,
      layout: () => ({
        ...transformHandles(1),
        frame: Mat4.identity(),
        screenSize: 80,
        occludedAlpha: { handles: 0.5 },
      }),
    });
    await settle();
    expect(sized.added[0].props.screenSize).toBe(80);
    expect(sized.added[0].props.occludedAlpha).toBe(0.5);
    d.dispose();
  });

  // the camera looks down -z: the z arrow points at the viewer, the planes across x and y are edge-on
  it('fades the handles turned edge-on to the view, hidden past the pick threshold', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const c = new GizmoController(fake.canvas, {
      ...binding,
      layout: () => ({ ...transformHandles(1, { rotate: false }), frame: Mat4.identity() }),
    });
    await settle();
    const t = transparencyOf(fake.added[0]);
    // move, slide x, y, z; planes across x, y, z
    expect([0, 1, 2, 3, 4, 5, 6].map((g) => t.get(g) ?? 0)).toEqual([0, 0, 0, 1, 1, 1, 0]);

    // turned: the camera looks along the diagonal of x and z; the z arrow is back, the plane across y still edge-on
    Vec3.set(fake.fake.camera.state.position, 10, 0, 10);
    fake.fake.didDraw.next(0);
    const u = transparencyOf(fake.added[0]);
    expect([0, 1, 2, 3, 4, 5, 6].map((g) => u.get(g) ?? 0)).toEqual([0, 0, 0, 0, 0, 1, 0]);
    c.dispose();
  });

  it('keeps a slide drawn end-on with fadeSlides off; the planes still fade', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const c = new GizmoController(
      fake.canvas,
      { ...binding, layout: () => ({ ...transformHandles(1, { rotate: false }), frame: Mat4.identity() }) },
      { fadeSlides: false },
    );
    await settle();
    const t = transparencyOf(fake.added[0]);
    // move, slide x, y, z; planes across x, y, z: looking down -z, no slide faded, the planes across x and y are
    expect([0, 1, 2, 3, 4, 5, 6].map((g) => t.get(g) ?? 0)).toEqual([0, 0, 0, 0, 1, 1, 0]);

    // on and off while shown
    c.setFadeSlides(true);
    expect(transparencyOf(fake.added[0]).get(3)).toBe(1);
    c.setFadeSlides(false);
    expect(transparencyOf(fake.added[0]).get(3) ?? 0).toBe(0);
    c.dispose();
  });

  it('applies a turn by steps, each from the last value, and cancels to the value at the press', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const groups = [{ handle: { kind: 'turn' as const }, label: 'turn', color: Color(0) }];
    const seen: [number, number][] = [];
    const c = new GizmoController(fake.canvas, {
      ...binding,
      layout: () => ({ ...transformHandles(1, { rotate: false }), groups, frame: Mat4.identity() }),
      apply: (from, _h, _f, dx) => {
        seen.push([from.x, dx]);
        return { x: from.x + 100 };
      },
    });
    await settle();
    hoverOn(fake, 0);
    await drag(fake, [5, 8]);
    // the press, then each step from the last value
    expect(seen).toEqual([
      [0, 0],
      [100, 5],
      [200, 3],
    ]);
    c.cancel();
    await settle();
    expect(binding.cancel).toHaveBeenCalledWith('a', { x: 0 });
    c.dispose();
  });

  // the origin at (1, 0, 0): on the fake screen at (60, 50); the press at (50, 50), 10 px left of it
  it('applies a roll by steps, its grab the last cursor from the origin', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const groups = [{ handle: { kind: 'roll' as const }, label: 'roll', color: Color(0) }];
    const seen: { from: number; grab: readonly number[] }[] = [];
    const frame = Mat4.fromTranslation(Mat4(), Vec3.create(1, 0, 0));
    const c = new GizmoController(fake.canvas, {
      ...binding,
      layout: () => ({ ...headsOnly(transformHandles(1, { rotate: false })), groups, frame }),
      apply: (from, _h, f) => {
        seen.push({ from: from.x, grab: [...f.grab!] });
        return { x: from.x + 1 };
      },
    });
    await settle();
    hoverOn(fake, 0);
    await drag(fake, [5, 8]);
    expect(seen).toEqual([
      { from: 0, grab: [-10, 0] },
      { from: 1, grab: [-10, 0] },
      { from: 2, grab: [-5, 0] },
    ]);
    c.dispose();
  });

  it('moves on a press on a slide kept end-on, and slides on one seen from the side', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const apply = vi.fn(binding.apply);
    const c = new GizmoController(
      fake.canvas,
      { ...binding, apply, layout: () => ({ ...transformHandles(1, { rotate: false }), frame: Mat4.identity() }) },
      { fadeSlides: false },
    );
    await settle();
    // looking down -z: the z arrow's head covers the knob
    hoverOn(fake, 3);
    await drag(fake, [5]);
    expect(apply.mock.lastCall![1].kind).toBe('move');
    fake.fake.input.interactionEnd.next();
    await settle();

    // the x arrow lies across the view
    hoverOn(fake, 1);
    await drag(fake, [5]);
    expect(apply.mock.lastCall![1].kind).toBe('slide');
    c.dispose();
  });

  // the camera looks down -z: the z arrow faded out, in its hit area too
  it('adds a hit area: picked only, its picks the handles, faded with them', async () => {
    const fake = fakeCanvas();
    const binding = fakeBinding();
    const h = transformHandles(1, { planes: false, rotate: false });
    const c = new GizmoController(fake.canvas, {
      ...binding,
      layout: () => ({
        ...h,
        parts: h.parts.filter((p) => p.geometry === h.handles || p.geometry === h.hitArea),
        frame: Mat4.identity(),
      }),
    });
    await settle();
    expect(fake.added).toHaveLength(2);
    const hit = fake.added.find((r) => r.state.pickOnly)!;
    const handles = fake.added.find((r) => !r.state.pickOnly)!;
    expect(hit.state.layer).toBe('overlay');
    expect(hit.renderObjects[0].state.pickOnly).toBe(true);
    expect(handles.renderObjects[0].state.pickOnly).toBe(false);
    expect(c.owns(hit)).toBe(true);

    // a pick on it is a pick on the handles: their shape's loci (hover, highlight, clicks)
    const id = hit.renderObjects[0].id;
    const loci = hit.getLoci({ objectId: id, groupId: 1, instanceId: 0 }) as ShapeGroup.Loci;
    expect(loci.shape).toBe((handles.getAllLoci()[0] as Shape.Loci).shape);
    expect(OrderedSet.start(loci.groups[0].ids)).toBe(1);

    // hovered through it: the handle's drag
    fake.fake.interaction.hover.next({ current: { repr: hit, loci } });
    expect(fake.fake.attribs.trackball.bindings.dragRotate).toBe(Binding.Empty);

    expect(transparencyOf(hit)).toEqual(transparencyOf(handles));
    expect(transparencyOf(hit).get(3)).toBe(1);
    c.dispose();
  });

  // the camera looks down -z: the z arrow faded out, its shaft (a line) too
  it('draws thin handles as lines: their picks and marks the handles, faded with them', async () => {
    const fake = fakeCanvas();
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => {
        const h = transformHandles(1, { planes: false, rotate: false });
        // the shafts as lines: from the center to each arrow's head (the slides' groups 1 to 3)
        const shafts = LinesBuilder.create(4, 4);
        for (const a of [0, 1, 2]) {
          const end = [0, 0, 0];
          end[a] = 0.8;
          shafts.add(0, 0, 0, end[0], end[1], end[2], a + 1);
        }
        return {
          ...h,
          parts: [
            { role: 'handle', geometry: h.handles },
            { role: 'handle', geometry: shafts.getLines() },
          ],
          frame: Mat4.identity(),
        };
      },
    });
    await settle();
    expect(fake.added).toHaveLength(2);
    const lines = fake.added.find((r) => r.renderObjects[0].type === 'lines')!;
    const handles = fake.added.find((r) => r.renderObjects[0].type === 'mesh')!;
    expect(lines.state.layer).toBe('overlay');
    expect(lines.props.lineSizeAttenuation).toBe(false);

    // a pick on a line is the handle's
    const loci = lines.getLoci({ objectId: lines.renderObjects[0].id, groupId: 1, instanceId: 0 }) as ShapeGroup.Loci;
    const handleShape = (handles.getAllLoci()[0] as Shape.Loci).shape;
    expect(loci.shape).toBe(handleShape);
    // a highlight of the handle marks its line
    expect(lines.mark(ShapeGroup.Loci(handleShape, [{ ids: OrderedSet.ofSingleton(1), instance: 0 }]), 1)).toBe(true);
    // hovered through it: the handle's drag
    fake.fake.interaction.hover.next({ current: { repr: lines, loci } });
    expect(fake.fake.attribs.trackball.bindings.dragRotate).toBe(Binding.Empty);

    expect(transparencyOf(lines).get(3)).toBe(1);
    expect(transparencyOf(lines)).toEqual(transparencyOf(handles));
    c.dispose();
  });

  // the handles in two parts: in world units, and at a size of their own on screen (a clip's move handles beside its
  // outline), both the handles' loci
  it('draws fixed-size handles beside the world ones, their picks the handles', async () => {
    const fake = fakeCanvas();
    const { handles, groups, hitArea } = transformHandles(1, { planes: false, rotate: false });
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({
        parts: [
          { role: 'handle', geometry: transformHandles(4, { move: false, rotate: false }).handles },
          { role: 'handle', geometry: handles, screenSize: 90 },
          { role: 'hit', geometry: hitArea, screenSize: 90 },
        ],
        groups,
        frame: Mat4.identity(),
      }),
    });
    await settle();
    // the world handles, the fixed ones and their hit area
    expect(fake.added).toHaveLength(3);
    const world = fake.added.find((r) => r.props.screenSize === 0)!;
    const fixed = fake.added.find((r) => r.props.screenSize === 90 && !r.state.pickOnly)!;
    const fixedHit = fake.added.find((r) => r.props.screenSize === 90 && r.state.pickOnly)!;
    expect([world, fixed, fixedHit].every((r) => r.state.layer === 'overlay')).toBe(true);

    const handleShape = (world.getAllLoci()[0] as Shape.Loci).shape;
    for (const r of [fixed, fixedHit]) {
      const loci = r.getLoci({ objectId: r.renderObjects[0].id, groupId: 1, instanceId: 0 }) as ShapeGroup.Loci;
      expect(loci.shape).toBe(handleShape);
    }
    // hovered through the fixed part: the handle's drag
    const loci = fixed.getLoci({ objectId: fixed.renderObjects[0].id, groupId: 1, instanceId: 0 });
    fake.fake.interaction.hover.next({ current: { repr: fixed, loci } });
    expect(fake.fake.attribs.trackball.bindings.dragRotate).toBe(Binding.Empty);
    expect(transparencyOf(fixed)).toEqual(transparencyOf(world));
    c.dispose();
  });

  // knobs as points: round, a set size on screen, their picks the handles'
  it('draws knobs as points of a set size, their picks the handles', async () => {
    const fake = fakeCanvas();
    const builder = PointsBuilder.create(4, 4);
    builder.add(1, 0, 0, 1);
    builder.add(0, 1, 0, 2);
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => {
        const h = headsOnly(transformHandles(1, { planes: false, rotate: false }));
        return { ...h, parts: [...h.parts, { role: 'handle', geometry: builder.getPoints() }], frame: Mat4.identity() };
      },
    });
    await settle();
    const points = fake.added.find((r) => r.renderObjects[0].type === 'points')!;
    expect(points.props.pointSizeAttenuation).toBe(false);
    expect(points.props.pointStyle).toBe('circle');
    const loci = points.getLoci({ objectId: points.renderObjects[0].id, groupId: 2, instanceId: 0 }) as ShapeGroup.Loci;
    const handles = fake.added.find((r) => r.renderObjects[0].type === 'mesh')!;
    expect(loci.shape).toBe((handles.getAllLoci()[0] as Shape.Loci).shape);
    expect(OrderedSet.start(loci.groups[0].ids)).toBe(2);
    c.dispose();
  });

  // a part facing the camera (a circle about the view axis): a point of its own size and look, its radius and depth in
  // units of the screen size
  it('draws a point part at its size, style and depth on screen', async () => {
    const fake = fakeCanvas();
    const h = headsOnly(transformHandles(1, { move: false, free: false }));
    const ring = PointsBuilder.create(1, 1);
    ring.add(0, 0, 0, 0);
    const point = { radius: 0.5, style: 'ring' as const, ringWidth: 0.2, depth: 0.25 };
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({
        ...h,
        parts: [...h.parts, { role: 'handle' as const, geometry: ring.getPoints(), point }],
        frame: Mat4.identity(),
        screenSize: 80,
      }),
    });
    await settle();
    const drawn = fake.added.find((r) => r.renderObjects[0].type === 'points')!;
    expect(drawn.props.pointStyle).toBe('ring');
    expect(drawn.props.pointRingWidth).toBe(0.2);
    expect(drawn.props.pointDepthOffset).toBe(20);
    expect(drawn.props.pointFront).toBe(false);
    expect(drawn.props.pointSizeAttenuation).toBe(false);
    c.dispose();
  });

  // an extra of the handles' loci is one of them (picked, faded); one of its own loci is drawn only
  it("adds a layout's extras in its frame, removed with it; a handle where its loci are the handles'", async () => {
    const fake = fakeCanvas();
    const h = headsOnly(transformHandles(1, { planes: false, rotate: false }));
    const frame = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    const extra = (standIn: boolean) => ({
      build: (graphics: GraphicsContext, handles: Shape) =>
        ShapeRepresentation((_, s: Shape<Mesh>) => s, Mesh.Utils, {
          graphics,
          ...(standIn ? { lociShape: () => handles } : {}),
        }),
      data: Shape.create(
        'x',
        {},
        h.handles,
        () => Color(0),
        () => 1,
        () => '',
        4,
      ),
    });
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({ ...h, extras: [extra(true), extra(false)], frame }),
    });
    await settle();
    expect(fake.added).toHaveLength(3);
    const [, handleExtra, otherExtra] = fake.added;
    for (const r of fake.added) {
      expect(r.state.layer).toBe('overlay');
      expect(Mat4.areEqual(r.state.transform, frame, 1e-9)).toBe(true);
    }
    const pick = (r: Representation.Any) => r.getLoci({ objectId: r.renderObjects[0].id, groupId: 1, instanceId: 0 });
    fake.fake.interaction.hover.next({ current: { repr: otherExtra, loci: pick(otherExtra) } });
    expect(fake.fake.attribs.trackball.bindings).toEqual(drags);
    fake.fake.interaction.hover.next({ current: { repr: handleExtra, loci: pick(handleExtra) } });
    expect(fake.fake.attribs.trackball.bindings.dragRotate).toBe(Binding.Empty);
    // the z arrow toward the viewer, faded in the extra as in the parts
    expect(transparencyOf(handleExtra).get(3)).toBeGreaterThan(0);
    expect(transparencyOf(handleExtra)).toEqual(transparencyOf(fake.added[0]));
    c.dispose();
    expect(fake.added).toHaveLength(0);
  });

  it('dims the other handles during a drag, and shows them again after', async () => {
    const fake = fakeCanvas();
    Vec3.set(fake.fake.camera.state.position, 6, 5, 8);
    const c = new GizmoController(fake.canvas, fakeBinding());
    await settle();
    expect(transparencyOf(fake.added[0]).size).toBe(0);
    hoverOn(fake, 1);
    await drag(fake, [10]);
    const t = transparencyOf(fake.added[0]);
    // the slide along x dragged; move, slide y and z dimmed
    expect([0, 1, 2, 3].map((g) => t.get(g) ?? 0)).toEqual([0.75, 0, 0.75, 0.75]);
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(transparencyOf(fake.added[0]).size).toBe(0);
    c.dispose();
  });

  // move 0, slides 1-3, rings 4-6, roll 7, trackball 8
  const rings = () => headsOnly(transformHandles(1, { planes: false }));

  it('hides every other handle during a turn, and the turning ones during a move (the other moving ones dimmed)', async () => {
    const fake = fakeCanvas();
    Vec3.set(fake.fake.camera.state.position, 6, 5, 8);
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({ ...rings(), frame: Mat4.identity() }),
    });
    await settle();
    const groups = [0, 1, 2, 3, 4, 5, 6, 7, 8];
    hoverOn(fake, 5);
    await drag(fake, [10]);
    const t = transparencyOf(fake.added[0]);
    expect(groups.map((g) => t.get(g) ?? 0)).toEqual([1, 1, 1, 1, 1, 0, 1, 1, 1]);
    fake.fake.input.interactionEnd.next();
    await settle();
    hoverOn(fake, 1);
    await drag(fake, [10]);
    const u = transparencyOf(fake.added[0]);
    expect(groups.map((g) => u.get(g) ?? 0)).toEqual([0.75, 0, 0.75, 0.75, 1, 1, 1, 1, 1]);
    c.dispose();
  });

  // a fill shown while its group is hovered, a shaft while a slide is dragged
  it('shows the parts of a state only in it: on hover, during a drag of their kinds', async () => {
    const fake = fakeCanvas();
    Vec3.set(fake.fake.camera.state.position, 6, 5, 8);
    const h = rings();
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({
        ...h,
        parts: [
          ...h.parts,
          { role: 'handle', geometry: h.handles, when: { hover: 8 }, alpha: 0.2 },
          { role: 'handle', geometry: h.handles, when: { drag: ['slide', 'move'] } },
        ],
        frame: Mat4.identity(),
      }),
    });
    await settle();
    const [, fill, shafts] = fake.added;
    expect(fill.props.alpha).toBe(0.2);
    expect([fill.state.visible, shafts.state.visible]).toEqual([false, false]);
    hoverOn(fake, 8);
    expect([fill.state.visible, shafts.state.visible]).toEqual([true, false]);
    // during a drag of the free rotation: no fill
    await drag(fake, [10]);
    expect(fill.state.visible).toBe(false);
    fake.fake.input.interactionEnd.next();
    await settle();
    // hovered again after: the fill back
    expect(fill.state.visible).toBe(true);
    hoverOn(fake, 1);
    expect(fill.state.visible).toBe(false);
    await drag(fake, [10]);
    expect(shafts.state.visible).toBe(true);
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(shafts.state.visible).toBe(false);
    c.dispose();
  });

  // the plane through the origin across the view, its normal away from the viewer (the far half cut); whole when turned
  it('cuts a front part to its half toward the viewer, as the camera moves, and shows it whole while turned', async () => {
    const fake = fakeCanvas();
    const h = rings();
    const frame = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({ ...h, parts: [{ ...h.parts[0], front: true }], frame }),
    });
    await settle();
    const [ringPart] = fake.added;
    const plane = () => ringPart.props.clip.objects[0];
    const normal = () => {
      const { axis, angle } = plane().rotation;
      return Vec3.transformDirection(
        Vec3(),
        Vec3.create(0, 1, 0),
        Mat4.fromRotation(Mat4(), (angle * Math.PI) / 180, axis),
      );
    };
    const away = (eye: Vec3) => Vec3.normalize(Vec3(), Vec3.sub(Vec3(), Vec3.create(1, 2, 3), eye));
    expect(plane().type).toBe('plane');
    expect([...plane().position]).toEqual([1, 2, 3]);
    normal().forEach((x, i) => expect(x).toBeCloseTo(away(Vec3.create(0, 0, 10))[i], 6));
    // the camera moved
    Vec3.set(fake.fake.camera.state.position, 10, 2, 3);
    fake.fake.didDraw.next(0);
    await settle();
    normal().forEach((x, i) => expect(x).toBeCloseTo(away(Vec3.create(10, 2, 3))[i], 6));
    // a ring turned: the plane far behind, nothing cut
    hoverOn(fake, 4);
    await drag(fake, [10]);
    await settle();
    expect(Vec3.distance(plane().position, Vec3.create(1, 2, 3))).toBeGreaterThan(1e5);
    fake.fake.input.interactionEnd.next();
    await settle();
    expect([...plane().position]).toEqual([1, 2, 3]);
    c.dispose();
  });

  // a shape representation keeps its geometry's buffers: the guides drawn follow the cursor only built into them
  it('updates the guides drawn as the cursor moves', async () => {
    const fake = fakeCanvas();
    Vec3.set(fake.fake.camera.state.position, 6, 5, 8);
    const c = new GizmoController(fake.canvas, {
      ...fakeBinding(),
      layout: () => ({ ...rings(), frame: Mat4.identity() }),
    });
    await settle();
    hoverOn(fake, 5);
    fake.fake.input.drag.next({ isStart: true, x: 70, y: 50, dx: 0, dy: 0 });
    await settle();
    const guides = fake.added.find((r) => r.renderObjects[0].type === 'lines')!;
    const drawn = () => Array.from(guides.renderObjects[0].values.tLine.ref.value.array as Float32Array);
    const before = drawn();
    // far to the left
    fake.fake.input.drag.next({ isStart: false, x: 5, y: 50, dx: 0, dy: 0 });
    await settle();
    expect(drawn()).not.toEqual(before);
    c.dispose();
  });

  it('draws a moved handle grey where it started and the guides along its axis, removed when it ends', async () => {
    const fake = fakeCanvas();
    Vec3.set(fake.fake.camera.state.position, 6, 5, 8);
    const frame = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    const c = new GizmoController(fake.canvas, { ...fakeBinding(), layout: () => ({ ...rings(), frame }) });
    await settle();
    const parts = [...fake.added];
    expect(parts).toHaveLength(1);
    hoverOn(fake, 1);
    await drag(fake, [10]);
    await settle();
    // a ghost of the part, the guides
    expect(fake.added).toHaveLength(3);
    const ghost = fake.added.find((r) => !parts.includes(r) && r.renderObjects[0].type === 'mesh')!;
    expect(ghost.state.pickable).toBe(false);
    expect(Mat4.areEqual(ghost.state.transform, frame, 1e-9)).toBe(true);
    // the other groups hidden
    const t = transparencyOf(ghost);
    expect([0, 1, 2, 7].map((g) => t.get(g) ?? 0)).toEqual([1, 0, 1, 1]);
    const guides = fake.added.find((r) => r.renderObjects[0].type === 'lines')!;
    expect(guides.state.pickable).toBe(false);
    expect(guides.state.layer).toBe('overlay');
    fake.fake.input.interactionEnd.next();
    await settle();
    expect(fake.added).toEqual(parts);
    c.dispose();
  });
});

/** Each line's start and end of a lines geometry */
function linesOf(lines: Lines) {
  const a = lines.lineData.ref.value.array;
  const out: { from: number[]; to: number[]; group: number }[] = [];
  for (let i = 0; i < lines.lineCount; i++) {
    const o = i * LineStride;
    out.push({ from: [a[o], a[o + 1], a[o + 2]], to: [a[o + 4], a[o + 5], a[o + 6]], group: a[o + 3] });
  }
  return out;
}

describe('guideLines', () => {
  const origin = Vec3.create(1, 0, 0);
  const frame = Mat4.mul(
    Mat4(),
    Mat4.fromTranslation(Mat4(), origin),
    Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.create(0, 0, 1)),
  );
  // looking down -Z from 10 Å, a tenth of an Å a pixel; the cursor 3 Å above the origin
  const at = (o: Partial<GuideFrame> = {}): GuideFrame => ({
    frame,
    origin,
    reach: 5,
    depth: 10,
    cursor: Vec3.create(1, 3, 0),
    right: Vec3.create(1, 0, 0),
    up: Vec3.create(0, 1, 0),
    perPixel: 0.1,
    ...o,
  });

  it("draws a slide's axis as the frame turns it, through the origin either way", () => {
    const [axis, ...rest] = linesOf(guideLines({ kind: 'slide', axis: Vec3.create(1, 0, 0) }, at()));
    expect(rest).toHaveLength(0);
    // local x turned onto world y
    [1, -5, 0].forEach((x, i) => expect(axis.from[i]).toBeCloseTo(x, 5));
    [1, 5, 0].forEach((x, i) => expect(axis.to[i]).toBeCloseTo(x, 5));
    expect(axis.group).toBe(0);
  });

  it("draws a plane's two axes", () => {
    expect(linesOf(guideLines({ kind: 'plane', normal: Vec3.create(0, 0, 1) }, at()))).toHaveLength(2);
  });

  it("draws a roll's radius at the start and dashes to the cursor, no axis", () => {
    const lines = linesOf(guideLines({ kind: 'roll' }, at({ start: Vec3.create(2, 0, 0) })));
    expect(lines.filter((l) => l.group === 0)).toHaveLength(0);
    expect(lines[0].to).toEqual([2, 0, 0]);
    expect(lines.length).toBeGreaterThan(2);
  });

  it("draws a ring's axis, its radius at the grab and now, dashes and a double arrow at the cursor", () => {
    const lines = linesOf(
      guideLines(
        { kind: 'rotate', axis: Vec3.create(0, 1, 0) },
        at({ start: Vec3.create(2, 0, 0), current: Vec3.create(1, 1, 0) }),
      ),
    );
    // the axis (local y turned onto world -x, in the screen: `reach` either way), the radius now
    const [axis, now] = lines.filter((l) => l.group === 0);
    expect(axis.from[0]).toBeCloseTo(6, 5);
    expect(axis.to[0]).toBeCloseTo(-4, 5);
    expect(now.to).toEqual([1, 1, 0]);
    const [radius, ...rest] = lines.filter((l) => l.group === 1);
    expect(radius.from).toEqual([1, 0, 0]);
    expect(radius.to).toEqual([2, 0, 0]);
    // dashes up from the origin to the cursor (6 px a dash: 0.6 Å)
    const dashes = rest.filter((l) => Math.abs(l.from[0] - 1) < 1e-5 && Math.abs(l.to[0] - 1) < 1e-5);
    expect(dashes.length).toBeGreaterThan(1);
    // the double arrow at the cursor, across the dashes: 9 px (0.9 Å) either way along x
    const across = rest.find((l) => Math.abs(l.from[1] - 3) < 1e-5 && Math.abs(l.to[1] - 3) < 1e-5)!;
    expect(Math.abs(across.from[0] - across.to[0])).toBeCloseTo(1.8, 5);
  });

  it("ends a ring's axis toward the viewer at half the origin's depth, not past the eye", () => {
    // 20° off the view (outside the end-on cone, `Fade.slide`): local (0, -sin, cos) turned onto world (sin, 0, cos)
    const t = (20 * Math.PI) / 180;
    const lines = linesOf(
      guideLines({ kind: 'rotate', axis: Vec3.create(0, -Math.sin(t), Math.cos(t)) }, at({ reach: 50 })),
    );
    const [axis] = lines.filter((l) => l.group === 0);
    // toward the viewer: 5 Å in front of the eye (at z 10); away: the whole reach
    expect(axis.to[2]).toBeCloseTo(5, 5);
    expect(axis.to[0] - origin[0]).toBeCloseTo(5 * Math.tan(t), 5);
    expect(axis.from[2]).toBeCloseTo(-50 * Math.cos(t), 5);
  });

  it("leaves out a ring's axis end-on to the view", () => {
    // 5° off the view
    const t = (5 * Math.PI) / 180;
    const axis = Vec3.create(0, -Math.sin(t), Math.cos(t));
    const lines = linesOf(guideLines({ kind: 'rotate', axis }, at({ current: Vec3.create(1, 1, 0) })));
    // only the radius now
    expect(lines.filter((l) => l.group === 0).map((l) => l.to)).toEqual([[1, 1, 0]]);
  });

  it("draws a cross where the trackball's turn is about", () => {
    const lines = linesOf(guideLines({ kind: 'turn' }, at()));
    expect(lines).toHaveLength(2);
    for (const l of lines) {
      [0, 1, 2].forEach((i) => expect((l.from[i] + l.to[i]) / 2).toBeCloseTo(origin[i], 5));
      expect(l.group).toBe(1);
    }
  });
});

/** The transparency of each group of a gizmo's handles (its state's layers) */
function transparencyOf(repr: Representation.Any) {
  const out = new Map<number, number>();
  for (const { loci, value } of repr.state.transparency.layers) {
    OrderedSet.forEach((loci as ShapeGroup.Loci).groups[0].ids, (g) => void out.set(g, value));
  }
  return out;
}
