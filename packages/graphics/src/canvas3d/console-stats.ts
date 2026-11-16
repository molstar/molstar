/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getVertexArenaStats } from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import { isTimingMode } from '@molstar/core/util/debug';
import { toFixed } from '@molstar/core/util/number';
import type { Canvas3DParts } from './parts.js';

/** Logs a canvas's render items, memory and scene to the console (debug mode, the console stats providers) */
export function logCanvas3DStats(parts: Canvas3DParts) {
  const { graphics, scene, passes, hiZ } = parts;
  const items = scene.renderables.map((r) => ({
    drawCount: r.values.drawCount.ref.value,
    instanceCount: r.values.instanceCount.ref.value,
    materialId: r.materialId,
    renderItemId: r.id,
    geometryType: r.values.dGeometryType.ref.value,
    'byteCount [MiB]': toFixed(r.getByteCount() / 1024 / 1024, 3),
    radius: toFixed(r.values.boundingSphere.ref.value.radius, 3),
  }));

  console.groupCollapsed(`${items.length} RenderItems`);

  if (items.length <= 64) {
    console.table(items);
  } else {
    console.log(items);
  }
  console.log(JSON.stringify(graphics.stats, undefined, 4));

  const { texture, cubeTexture, attribute, elements, pixelPack, renderbuffer, storageBuffer } =
    graphics.resources.getByteCounts();
  const arena = graphics.compute ? getVertexArenaStats(graphics.compute) : undefined;
  console.log(
    JSON.stringify(
      {
        texture: `${(texture / 1024 / 1024).toFixed(3)} MiB`,
        cubeTexture: `${(cubeTexture / 1024 / 1024).toFixed(3)} MiB`,
        attribute: `${(attribute / 1024 / 1024).toFixed(3)} MiB`,
        elements: `${(elements / 1024 / 1024).toFixed(3)} MiB`,
        pixelPack: `${(pixelPack / 1024 / 1024).toFixed(3)} MiB`,
        renderbuffer: `${(renderbuffer / 1024 / 1024).toFixed(3)} MiB`,
        storageBuffer: `${(storageBuffer / 1024 / 1024).toFixed(3)} MiB`,
        ...(arena
          ? {
              vertexArena: `${(arena.usedByteCount / 1024 / 1024).toFixed(3)} of ${(arena.byteCount / 1024 / 1024).toFixed(3)} MiB used (${arena.usedPages} of ${arena.pages} pages, ${arena.arenas} arenas${arena.exhausted ? ', out of memory' : ''})`,
            }
          : {}),
      },
      undefined,
      4,
    ),
  );

  console.log(
    JSON.stringify(
      {
        renderables: `${(scene.renderables.reduce((sum, r) => sum + r.getByteCount(), 0) / 1024 / 1024).toFixed(3)} MiB`,
        passes: {
          draw: `${(passes.draw.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
          illumination: `${(passes.illumination.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
          pick: `${(passes.pick.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
          hiZ: `${(hiZ.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
        },
      },
      undefined,
      4,
    ),
  );

  console.log(
    JSON.stringify(
      {
        scene: {
          radius: toFixed(scene.boundingSphere.radius, 3),
          radiusVisible: toFixed(scene.boundingSphereVisible.radius, 3),
          transparency: scene.transparency,
          markerAverage: toFixed(scene.markerAverage, 2),
          emissiveAverage: toFixed(scene.emissiveAverage, 2),
          wiggleAverage: toFixed(scene.wiggleAverage, 2),
          opacityAverage: toFixed(scene.opacityAverage, 2),
          transparencyMin: toFixed(scene.transparencyMin, 2),
          hasOpaque: scene.hasOpaque,
          hasAnimation: scene.hasAnimation,
          counts: {
            all: scene.renderables.length,
            primitives: scene.primitives.renderables.length,
            volumes: scene.volumes.renderables.length,
          },
        },
      },
      undefined,
      4,
    ),
  );

  if (isTimingMode) {
    console.log(JSON.stringify(graphics.timer.formatedStats(), undefined, 4));
  }

  console.groupEnd();
}
