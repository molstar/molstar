/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { now } from '@molstar/core/util/now';
import { LinesBuilder } from '@molstar/graphics/geo/geometry/lines/lines-builder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { PointsBuilder } from '@molstar/graphics/geo/geometry/points/points-builder';
import {
  type GizmoBinding,
  GizmoController,
  type GizmoControllerOptions,
} from '@molstar/graphics/util/gizmo/controller';
import { transformHandles } from '@molstar/graphics/util/gizmo/transform';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { BehaviorSubject } from 'rxjs';
import { WebGPU } from './webgpu-module.js';

const Pdb = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n';
const Size = 256;
/** CSS pixels a unit of the handles is: the arrows 80 long */
const ScreenSize = 80;

/** A binding of one gizmo at the origin (inside the atom): move and the slides, at a fixed size on screen */
function binding(): GizmoBinding<number> {
  return {
    items: new BehaviorSubject([{ key: 'a', value: 0, label: 'A' }]),
    layout: () => ({
      ...transformHandles(1, { planes: false, rotate: false }),
      frame: Mat4.identity(),
      screenSize: ScreenSize,
    }),
    sameGeometry: () => true,
    equal: (a, b) => a === b,
    apply: (from) => from,
    preview: async () => {},
    commit: async () => true,
    cancel: () => {},
  };
}

/** One point at the origin of `group` */
function pointAt(group: number) {
  const points = PointsBuilder.create(1, 1);
  points.add(0, 0, 0, group);
  return points.getPoints();
}

/** A circle facing the camera (a `ring` point: a band from 0.92 to 1.08), 1 unit a `ScreenSize` pixels */
function ringBinding(): GizmoBinding<number> {
  return {
    ...binding(),
    layout: () => ({
      parts: [{ role: 'handle', geometry: pointAt(0), point: { radius: 1.08, style: 'ring', ringWidth: 0.16 / 1.08 } }],
      groups: [{ handle: { kind: 'scale' }, label: 'Drag to size it', color: Color(0x9db8e0) }],
      frame: Mat4.identity(),
      screenSize: ScreenSize,
    }),
  };
}

/** A knob at the origin (group 0) and, as points facing the camera, a circle with a disc far behind it (group 1) */
function facingPartBinding(): GizmoBinding<number> {
  const knob = MeshBuilder.createState(1024, 512);
  knob.currentGroup = 0;
  addSphere(knob, Vec3(), 0.25, 2);
  return {
    ...binding(),
    layout: () => ({
      parts: [
        { role: 'handle', geometry: MeshBuilder.getMesh(knob), screenSize: ScreenSize },
        {
          role: 'handle',
          geometry: pointAt(1),
          screenSize: ScreenSize,
          point: { radius: 1.08, style: 'ring', ringWidth: 0.16 / 1.08 },
        },
        { role: 'hit', geometry: pointAt(1), screenSize: ScreenSize, point: { radius: 1, depth: -1.25 } },
      ],
      groups: [
        { handle: { kind: 'move' }, label: 'Drag to move', color: Color(0x9db8e0) },
        { handle: { kind: 'turn' }, label: 'Drag to turn', color: Color(0x9db8e0) },
      ],
      frame: Mat4.identity(),
    }),
  };
}

/** A knob at the origin, `ScreenSize` pixels a unit, and a square outline `2 * Reach` wide in world units */
const Reach = 15;
function twoPartBinding(): GizmoBinding<number> {
  const knob = MeshBuilder.createState(1024, 512);
  knob.currentGroup = 0;
  addSphere(knob, Vec3(), 0.25, 2);
  const outline = LinesBuilder.create(4, 1);
  const corners = [
    [-Reach, -Reach],
    [Reach, -Reach],
    [Reach, Reach],
    [-Reach, Reach],
  ];
  corners.forEach(([x, y], i) => {
    const [nx, ny] = corners[(i + 1) % 4];
    outline.add(x, y, 0, nx, ny, 0, 1);
  });
  return {
    ...binding(),
    layout: () => ({
      parts: [
        { role: 'handle', geometry: outline.getLines() },
        { role: 'handle', geometry: MeshBuilder.getMesh(knob), screenSize: ScreenSize },
      ],
      groups: [
        { handle: { kind: 'move' }, label: 'Drag to move', color: Color(0x9db8e0) },
        { handle: { kind: 'scale' }, label: 'Drag to size it', color: Color(0x9db8e0) },
      ],
      frame: Mat4.identity(),
    }),
  };
}

/** A knob (a point, group 0) at `KnobAt` and a line (group 1) from it toward the camera at (0, 0, 80) and the screen's center */
const KnobAt = Vec3.create(15, 0, 0);
function knobAndLineBinding(): GizmoBinding<number> {
  const knob = PointsBuilder.create(1, 1);
  knob.add(KnobAt[0], KnobAt[1], KnobAt[2], 0);
  const line = LinesBuilder.create(1, 1);
  const toward = Vec3.create(0, 0, 40);
  line.add(KnobAt[0], KnobAt[1], KnobAt[2], toward[0], toward[1], toward[2], 1);
  return {
    ...binding(),
    layout: () => ({
      parts: [
        { role: 'handle', geometry: line.getLines() },
        { role: 'handle', geometry: knob.getPoints() },
      ],
      groups: [
        { handle: { kind: 'scale' }, label: 'Drag to size it', color: Color(0x9db8e0) },
        { handle: { kind: 'faces', faces: [{ axis: 0, sign: 1 }] }, label: 'Drag to move it', color: Color(0x9db8e0) },
      ],
      frame: Mat4.identity(),
    }),
  };
}

// the gizmo library's handles are in the overlay layer: drawn and picked on top of the structure they are inside, at
// their size on screen
/**
 * A plugin with an atom's spacefill, the camera at `position` looking at it (y up), and a controller of `gizmos`
 * (its gizmo made); picks in CSS pixels
 */
async function setup(gizmos: GizmoBinding<number>, position: Vec3, options?: GizmoControllerOptions) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
  const canvas = plugin.canvas3d!;
  canvas.setProps({ camera: { ...canvas.props.camera, helper: { axes: { name: 'off', params: {} } } } });

  let time = now();
  const pickAt = async (x: number, y: number) => {
    canvas.commit(true);
    canvas.requestDraw();
    canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
    await canvas.graphics.resources.pending.whenReady();
    const pick = await canvas.asyncIdentify(Vec2.create(x, y))?.wait();
    return pick ? canvas.getLoci(pick.id) : undefined;
  };
  // the camera reset to the atom first
  await pickAt(0, 0);
  canvas.camera.setState({ position, target: Vec3(), up: Vec3.unitY }, 0);
  const controller = new GizmoController(canvas, gizmos, options);
  const groupAt = async (x: number, y: number) => {
    const picked = await pickAt(x, y);
    if (!picked || !controller.owns(picked.repr) || picked.loci.kind !== 'group-loci') return undefined;
    // a shape's group loci (`ShapeGroup.Loci`)
    return OrderedSet.start((picked.loci as { groups: { ids: OrderedSet }[] }).groups[0].ids);
  };
  // made asynchronously
  for (let i = 0; i < 20 && !(await isMade()); i++);
  async function isMade() {
    for (const [x, y] of [
      [Size / 2, Size / 2],
      [Size / 2 + ScreenSize, Size / 2],
    ])
      if ((await groupAt(x, y)) !== undefined) return true;
    return false;
  }
  return { plugin, canvas, controller, pickAt, groupAt };
}

describe.skipIf(!WebGPU.module)('gizmos in the overlay layer on WebGPU', () => {
  it('picks a handle inside the structure, at its size on screen, and not one faded edge-on', async () => {
    // looking at the atom along the diagonal of x and z
    const { plugin, canvas, controller, pickAt, groupAt } = await setup(binding(), Vec3.create(30, 0, 30));
    try {
      // the y arrow (group 2) straight up from the center: its head, no shaft at rest, nothing past its tip
      const c = Size / 2;
      expect(await groupAt(c, c - 0.92 * ScreenSize)).toBe(2);
      expect(await groupAt(c, c - 0.5 * ScreenSize)).toBe(undefined);
      expect(await groupAt(c, c - 1.3 * ScreenSize)).toBe(undefined);
      // the center, on top of the atom it is inside
      expect(await groupAt(c, c)).toBe(0);
      // zoomed in: the same size on screen
      canvas.camera.setState({ position: Vec3.create(10, 0, 10) }, 0);
      expect(await groupAt(c, c - 0.92 * ScreenSize)).toBe(2);

      // looking down y: the y arrow points at the viewer, in front of the knob, but faded out and not picked
      canvas.camera.setState({ position: Vec3.create(0, 30, 0.001), target: Vec3(), up: Vec3.unitZ }, 0);
      await pickAt(0, 0);
      expect(await groupAt(c, c)).toBe(0);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // a circle facing the camera as a point: round from any side, not edge-on (a line across the center), in the first
  // frame of every view, the pick's included
  it('draws a circle facing the camera round from any side', async () => {
    const { plugin, canvas, controller, groupAt } = await setup(ringBinding(), Vec3.create(30, 0, 0));
    try {
      const c = Size / 2;
      for (const position of [Vec3.create(30, 0, 0), Vec3.create(0, 5, 30), Vec3.create(-20, -20, 5)]) {
        canvas.camera.setState({ position }, 0);
        // above and right of the center at its radius, nothing at the center
        expect(await groupAt(c, c - ScreenSize)).toBe(0);
        expect(await groupAt(c + ScreenSize, c)).toBe(0);
        expect(await groupAt(c, c)).toBe(undefined);
      }
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // the transform handles' trackball: the area inside the rings picked as it, behind every other handle; the center in
  // front of every other; a ring's front half only
  it('picks the trackball inside the rings, the center over an edge-on ring, a ring by its front half', async () => {
    const freeBinding: GizmoBinding<number> = {
      ...binding(),
      layout: () => ({
        ...transformHandles(1, { planes: false }),
        frame: Mat4.identity(),
        screenSize: ScreenSize,
      }),
    };
    // groups: center 0, slides 1-3, rings 4-6, roll 7, trackball 8; looking along the diagonal of x and z (the axes
    // on the screen's vertical and horizontal, none on its diagonals)
    const { plugin, controller, groupAt } = await setup(freeBinding, Vec3.create(30, 0, 30));
    try {
      const c = Size / 2;
      const S = ScreenSize;
      for (let i = 0; i < 20 && (await groupAt(c + 0.25 * S, c - 0.25 * S)) === undefined; i++);
      expect(await groupAt(c + 0.25 * S, c - 0.25 * S)).toBe(8);
      expect(await groupAt(c - 0.25 * S, c + 0.25 * S)).toBe(8);
      // the center over the y ring, edge-on across it
      expect(await groupAt(c, c)).toBe(0);
      // the y arrow's head past the rings
      expect(await groupAt(c, c - 0.92 * S)).toBe(2);
      // the roll circle on the diagonal
      const r = (0.8 * S) / Math.SQRT2;
      expect(await groupAt(c + r, c - r)).toBe(7);
      // past the roll circle
      expect(await groupAt(c + 0.75 * S, c - 0.75 * S)).toBe(undefined);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // the x ring alone (about x, in yz; seen along the diagonal the z ring would lie on its far half): 45° from its top,
  // its near half left of the center, its far half right of it cut (the trackball there)
  it('picks a ring by its front half only', async () => {
    const ringOnly: GizmoBinding<number> = {
      ...binding(),
      layout: () => ({
        ...transformHandles(1, { planes: false, rings: [0] }),
        frame: Mat4.identity(),
        screenSize: ScreenSize,
      }),
    };
    // groups: center 0, slides 1-3, the x ring 4, roll 5, trackball 6
    const { plugin, controller, groupAt } = await setup(ringOnly, Vec3.create(30, 0, 30));
    try {
      const c = Size / 2;
      const [x, y] = [(0.65 * ScreenSize * Math.SQRT1_2) / Math.SQRT2, 0.65 * ScreenSize * Math.SQRT1_2];
      for (let i = 0; i < 20 && (await groupAt(c - x, c - y)) === undefined; i++);
      expect(await groupAt(c - x, c - y)).toBe(4);
      expect(await groupAt(c + x, c - y)).toBe(6);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // a knob drawn in front of the lines ending in it (`Points.Params.pointFront`): the line running from it toward the
  // camera is nearer than the knob's center over the knob's pixels, yet the knob is picked there
  it('picks a knob over a line ending in it', async () => {
    // a knob 24 px across: the pick is drawn at a fraction of the pixels
    const { plugin, canvas, controller, groupAt } = await setup(knobAndLineBinding(), Vec3.create(0, 0, 80), {
      pointSize: 24,
    });
    try {
      const p = canvas.camera.project(Vec4(), KnobAt);
      const [x, y] = [p[0], Size - p[1]];
      // the lines' pipelines may compile after the points'
      for (let i = 0; i < 40 && ((await groupAt(x, y)) === undefined || (await groupAt(x - 18, y)) === undefined); i++);
      // past the knob, along the line
      expect(await groupAt(x - 18, y)).toBe(1);
      // toward the center of the screen (where the line goes), within the knob's radius (12 px)
      expect(await groupAt(x - 6, y)).toBe(0);
      expect(await groupAt(x, y)).toBe(0);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // parts of a fixed part facing the camera, the frame not: the disk behind everything takes what the knob leaves
  it('faces a part of the fixed handles to the camera, its disk behind the others', async () => {
    const { plugin, canvas, controller, groupAt } = await setup(facingPartBinding(), Vec3.create(30, 0, 0));
    try {
      const c = Size / 2;
      for (const position of [Vec3.create(30, 0, 0), Vec3.create(0, 5, 30), Vec3.create(-20, -20, 5)]) {
        canvas.camera.setState({ position }, 0);
        for (let i = 0; i < 20 && (await groupAt(c + 0.6 * ScreenSize, c)) === undefined; i++);
        expect(await groupAt(c, c)).toBe(0);
        expect(await groupAt(c + 0.6 * ScreenSize, c)).toBe(1);
        expect(await groupAt(c, c - 0.6 * ScreenSize)).toBe(1);
        expect(await groupAt(c + 1.3 * ScreenSize, c)).toBe(undefined);
      }
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // a slide kept end-on: drawn over the center, the center's hit disk still in front of it
  it('picks the center over a slide kept drawn when it points at the viewer', async () => {
    const { plugin, controller, groupAt } = await setup(binding(), Vec3.create(0, 30, 0.001), { fadeSlides: false });
    try {
      const c = Size / 2;
      expect(await groupAt(c, c)).toBe(0);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // a gizmo in two parts: the outline at the scene's scale, the knob at a fixed size on screen, one gizmo's handles
  it('keeps fixed handles at their size on screen while the outline scales, picked along it', async () => {
    const { plugin, canvas, controller, groupAt } = await setup(twoPartBinding(), Vec3.create(0, 0, 80));
    try {
      const c = Size / 2;
      const widths: number[] = [];
      for (const distance of [80, 50]) {
        canvas.camera.setState({ position: Vec3.create(0, 0, distance), target: Vec3(), up: Vec3.unitY }, 0);
        await groupAt(0, 0);
        /** where a point of the outline's right edge is on the canvas, `y` of the reach above the center */
        const edgeAt = (y: number) => {
          const p = canvas.camera.project(Vec4(), Vec3.create(Reach, y * Reach, 0));
          return [p[0], Size - p[1]] as const;
        };
        // the lines' pipelines may compile after the knob's (what `setup` waits for)
        for (let i = 0; i < 20 && (await groupAt(...edgeAt(0))) === undefined; i++);
        // the knob: its radius a quarter of a unit at any zoom
        expect(await groupAt(c + 0.2 * ScreenSize, c)).toBe(0);
        expect(await groupAt(c + 0.2 * ScreenSize, c + 0.4 * ScreenSize)).toBe(undefined);
        // the right edge where the world puts it, anywhere along it; nothing past it
        for (const y of [-0.7, 0, 0.7]) expect(await groupAt(...edgeAt(y))).toBe(1);
        const [x] = edgeAt(0);
        expect(await groupAt(x + 12, c)).toBe(undefined);
        widths.push(x - c);
      }
      // nearer: the outline larger on screen, the knob as it was
      expect(widths[1]).toBeGreaterThan(widths[0] * 1.4);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });

  // lines are at least a pixel of the target they are drawn in (`max(1.0, linewidth)`), the pick's too: a thin line
  // handle is picked anywhere along it without a hit area (a mesh shaft as thin is not: it falls between the pixels)
  it('picks a line handle at its middle without a hit area', async () => {
    const lineBinding: GizmoBinding<number> = {
      ...binding(),
      layout: () => {
        const { groups } = transformHandles(1, { planes: false, rotate: false });
        // the shafts as lines in the world (a line's width set on screen, its length the world's): the slides' groups
        const shafts = LinesBuilder.create(4, 4);
        for (const a of [0, 1, 2]) {
          const end = [0, 0, 0];
          end[a] = 40;
          shafts.add(0, 0, 0, end[0], end[1], end[2], a + 1);
        }
        return { parts: [{ role: 'handle', geometry: shafts.getLines() }], groups, frame: Mat4.identity() };
      },
    };
    const { plugin, controller, groupAt } = await setup(lineBinding, Vec3.create(30, 0, 30));
    try {
      const c = Size / 2;
      // the lines' pipelines may compile after the knob's (what `setup` waits for)
      for (let i = 0; i < 20 && (await groupAt(c, c - 50)) === undefined; i++);
      // the y shaft (group 2), a line, along it (past the atom at the origin) and a pixel aside
      for (const dy of [40, 50, 60]) expect(await groupAt(c, c - dy)).toBe(2);
      expect(await groupAt(c - 3, c - 50)).toBe(2);
      controller.dispose();
    } finally {
      plugin.dispose();
    }
  });
});
