/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The hover the allocation cases measure (hover.quick.ts): a Canvas3D without canvas on the GL shim or the stub device, a
// row of spheres, and the step moving a highlight from one sphere to the next with a tick.

import { BehaviorSubject, Subject } from 'rxjs';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { InputObserver } from '@molstar/core/util/input/input-observer';
import { MarkerAction } from '@molstar/core/util/marker-action';
import type { now } from '@molstar/core/util/now';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import { Canvas3D, Canvas3DContext, Canvas3DParams } from '../../canvas3d/canvas3d.js';
import { Passes } from '../../canvas3d/passes/passes.js';
import { addSphere } from '../../geo/geometry/mesh/builder/sphere.js';
import { Mesh } from '../../geo/geometry/mesh/mesh.js';
import { MeshBuilder } from '../../geo/geometry/mesh/mesh-builder.js';
import { ShapeRepresentation } from '../../repr/shape/representation.js';
import type { GraphicsContext } from '../context.js';
import { createContext } from '../webgl/context.js';
import { createWebGPUContext } from '../webgpu/context.js';
import { createGl } from './gl.shim.js';
import { createStubDevice } from './stub-gpu.js';

const Size = 64;
const SphereCount = 8;

/** A Canvas3D context without canvas or input, on `graphics` (as the headless plugin's) */
function createHeadlessContext(graphics: GraphicsContext): Canvas3DContext {
  const assetManager = new DefaultAssetManager();
  const props = { ...Canvas3DContext.DefaultProps };
  const passes = new Passes(graphics, assetManager, props);
  const input = InputObserver.create();
  return {
    graphics,
    input,
    passes,
    attribs: { ...Canvas3DContext.DefaultAttribs },
    props,
    contextLost: new Subject<now.Timestamp>(),
    contextRestored: new Subject<now.Timestamp>(),
    contextChanged: new Subject<Canvas3DContext.ContextChange>(),
    assetManager,
    changed: new BehaviorSubject<undefined>(undefined),
    pixelScale: 1,
    syncPixelScale: () => {},
    setProps: () => {},
    release: () => false,
    restore: async () => false,
    isReleased: false,
    dispose: () => {
      input.dispose();
      graphics.destroy();
    },
  };
}

/** A row of spheres, a group each, as one shape representation */
async function createSpheres() {
  const state = MeshBuilder.createState(4096, 4096);
  for (let g = 0; g < SphereCount; ++g) {
    state.currentGroup = g;
    addSphere(state, Vec3.create(g * 3 - 10, 0, 0), 1, 2);
  }
  const shape = Shape.create(
    'spheres',
    {},
    MeshBuilder.getMesh(state),
    () => Color(0x3f7fbf),
    () => 1,
    () => '',
    SphereCount,
  );
  const repr = ShapeRepresentation((_: unknown, s: Shape<Mesh>) => s, Mesh.Utils);
  await repr.createOrUpdate({}, shape).run();
  return { repr, shape };
}

export async function createHoverCase(backend: 'webgl' | 'webgpu') {
  const graphics: GraphicsContext =
    backend === 'webgl'
      ? createContext(createGl(Size, Size, { preserveDrawingBuffer: true }))
      : createWebGPUContext(createStubDevice(), {
          width: Size,
          height: Size,
          // the device probes these by drawing and reading back, which a stub can't: as a GPU has them
          fragDepthInBundles: true,
          viewportOutsideTarget: true,
        });
  const context = createHeadlessContext(graphics);
  // no scale line: its label's font atlas needs a 2d canvas
  const { camera } = PD.getDefaultValues(Canvas3DParams);
  const canvas3d = Canvas3D.create(context, {
    camera: { ...camera, helper: { ...camera.helper, scale: { name: 'off', params: {} } } },
  });
  const { repr, shape } = await createSpheres();
  canvas3d.add(repr);
  canvas3d.commit(true);
  canvas3d.requestCameraReset({ durationMs: 0 });
  // the loci a pick gives, one per sphere (made once: picking is not the tick's)
  const loci = Array.from({ length: SphereCount }, (_, g) =>
    ShapeGroup.Loci(shape, [{ ids: OrderedSet.ofSingleton(g), instance: 0 }]),
  );
  let time = 0;
  const tick = () => canvas3d.tick((time += 16) as now.Timestamp, { isSynchronous: true });
  // the first frames drawn, programs or pipelines ready
  for (let k = 0; k < 20; ++k) {
    const skipped = graphics.stats.skippedDraws;
    canvas3d.requestDraw();
    tick();
    await graphics.resources.pending.whenReady();
    if (k > 2 && graphics.stats.skippedDraws === skipped) break;
  }
  let current = 0;
  const hover = (i: number) => {
    const next = i % SphereCount;
    canvas3d.mark({ loci: loci[current], repr }, MarkerAction.RemoveHighlight);
    canvas3d.mark({ loci: loci[next], repr }, MarkerAction.Highlight);
    current = next;
    tick();
  };
  // every sphere hovered once: its marking programs or pipelines ready too
  for (let i = 0; i < 2 * SphereCount; ++i) {
    hover(i);
    await graphics.resources.pending.whenReady();
  }
  return { canvas3d, context, hover };
}
