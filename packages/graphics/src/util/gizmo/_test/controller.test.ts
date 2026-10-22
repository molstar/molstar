/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { Mat4, Vec3, type Vec4 } from '@molstar/core/math/linear-algebra';
import { Binding } from '@molstar/core/util/binding';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import { BehaviorSubject, Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import type { Canvas3D } from '../../../canvas3d/canvas3d.js';
import type { GraphicsContext } from '../../../gl/context.js';
import type { Representation } from '../../../repr/representation.js';
import { type GizmoBinding, GizmoController, type GizmoItem } from '../controller.js';
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
    layout: () => ({ ...transformHandles(1, { planes: false, rotate: false }), frame: Mat4.identity() }),
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
});
