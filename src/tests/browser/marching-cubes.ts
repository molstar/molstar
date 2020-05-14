/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Canvas3DParams, Canvas3D } from '../../mol-canvas3d/canvas3d';
import { ColorNames } from '../../mol-util/color/names';
import { PositionData, Box3D, Sphere3D } from '../../mol-math/geometry';
import { OrderedSet } from '../../mol-data/int';
import { Vec3 } from '../../mol-math/linear-algebra';
import { computeGaussianDensityTexture2d, computeGaussianDensity } from '../../mol-math/geometry/gaussian-density';
import { calcActiveVoxels } from '../../mol-gl/compute/marching-cubes/active-voxels';
import { createHistogramPyramid } from '../../mol-gl/compute/histogram-pyramid/reduction';
import { createIsosurfaceBuffers } from '../../mol-gl/compute/marching-cubes/isosurface';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { Color } from '../../mol-util/color';
import { createRenderObject } from '../../mol-gl/render-object';
import { Representation } from '../../mol-repr/representation';
import { computeMarchingCubesMesh } from '../../mol-geo/util/marching-cubes/algorithm';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

const parent = document.getElementById('app')!;
parent.style.width = '100%';
parent.style.height = '100%';

const canvas = document.createElement('canvas');
parent.appendChild(canvas);
resizeCanvas(canvas, parent);

const canvas3d = Canvas3D.fromCanvas(canvas, PD.merge(Canvas3DParams, PD.getDefaultValues(Canvas3DParams), {
    renderer: { backgroundColor: ColorNames.white },
    camera: { mode: 'orthographic' }
}));
canvas3d.animate();

async function init() {
    const { webgl } = canvas3d;

    const position: PositionData = {
        x: [0, 2],
        y: [0, 2],
        z: [0, 2],
        indices: OrderedSet.ofSortedArray([0, 1]),
    };
    const box = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(2, 2, 2));
    const radius = () => 1.8;
    const props = {
        resolution: 0.1,
        radiusOffset: 0,
        smoothness: 1.5
    };
    const isoValue = Math.exp(-props.smoothness);

    if (true) {
        console.time('gpu gaussian2');
        const densityTextureData2 = await computeGaussianDensityTexture2d(position, box, radius, props, webgl).run();
        webgl.waitForGpuCommandsCompleteSync();
        console.timeEnd('gpu gaussian2');

        console.time('gpu mc2');
        console.time('gpu mc active2');
        const activeVoxelsTex2 = calcActiveVoxels(webgl, densityTextureData2.texture, densityTextureData2.gridDim, densityTextureData2.gridTexDim, isoValue, densityTextureData2.gridTexScale);
        webgl.waitForGpuCommandsCompleteSync();
        console.timeEnd('gpu mc active2');

        console.time('gpu mc pyramid2');
        const compacted2 = createHistogramPyramid(webgl, activeVoxelsTex2, densityTextureData2.gridTexScale);
        webgl.waitForGpuCommandsCompleteSync();
        console.timeEnd('gpu mc pyramid2');

        console.time('gpu mc vert2');
        createIsosurfaceBuffers(webgl, activeVoxelsTex2, densityTextureData2.texture, compacted2, densityTextureData2.gridDim, densityTextureData2.gridTexDim, densityTextureData2.transform, isoValue);
        webgl.waitForGpuCommandsCompleteSync();
        console.timeEnd('gpu mc vert2');
        console.timeEnd('gpu mc2');
    }

    console.time('gpu gaussian');
    const densityTextureData = await computeGaussianDensityTexture2d(position, box, radius, props, webgl).run();
    webgl.waitForGpuCommandsCompleteSync();
    console.timeEnd('gpu gaussian');

    console.time('gpu mc');
    console.time('gpu mc active');
    const activeVoxelsTex = calcActiveVoxels(webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, isoValue, densityTextureData.gridTexScale);
    webgl.waitForGpuCommandsCompleteSync();
    console.timeEnd('gpu mc active');

    console.time('gpu mc pyramid');
    const compacted = createHistogramPyramid(webgl, activeVoxelsTex, densityTextureData.gridTexScale);
    webgl.waitForGpuCommandsCompleteSync();
    console.timeEnd('gpu mc pyramid');

    console.time('gpu mc vert');
    const gv = createIsosurfaceBuffers(webgl, activeVoxelsTex, densityTextureData.texture, compacted, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.transform, isoValue);
    webgl.waitForGpuCommandsCompleteSync();
    console.timeEnd('gpu mc vert');
    console.timeEnd('gpu mc');

    console.log({ ...webgl.stats, programCount: webgl.stats.resourceCounts.program, shaderCount: webgl.stats.resourceCounts.shader });

    const mcBoundingSphere = Sphere3D.fromBox3D(Sphere3D(), densityTextureData.bbox);
    const mcIsosurface = TextureMesh.create(gv.vertexCount, 1, gv.vertexGroupTexture, gv.normalTexture, mcBoundingSphere);
    const mcIsoSurfaceProps = { doubleSided: true, flatShaded: true, alpha: 1.0 };
    const mcIsoSurfaceValues = TextureMesh.Utils.createValuesSimple(mcIsosurface, mcIsoSurfaceProps, Color(0x112299), 1);
    // console.log('mcIsoSurfaceValues', mcIsoSurfaceValues)
    const mcIsoSurfaceState = TextureMesh.Utils.createRenderableState(mcIsoSurfaceProps);
    const mcIsoSurfaceRenderObject = createRenderObject('texture-mesh', mcIsoSurfaceValues, mcIsoSurfaceState, -1);
    const mcIsoSurfaceRepr = Representation.fromRenderObject('texture-mesh', mcIsoSurfaceRenderObject);

    canvas3d.add(mcIsoSurfaceRepr);
    canvas3d.requestCameraReset();

    //

    console.time('cpu gaussian');
    const densityData = await computeGaussianDensity(position, box, radius, { ...props, useGpu: false }, webgl).run();
    console.timeEnd('cpu gaussian');
    console.log({ densityData });

    const params = {
        isoLevel: isoValue,
        scalarField: densityData.field,
        idField: densityData.idField
    };

    console.time('cpu mc');
    const surface = await computeMarchingCubesMesh(params).run();
    console.timeEnd('cpu mc');
    console.log('surface', surface);
    Mesh.transform(surface, densityData.transform);
    const meshProps = { doubleSided: true, flatShaded: false, alpha: 1.0 };
    const meshValues = Mesh.Utils.createValuesSimple(surface, meshProps, Color(0x995511), 1);
    const meshState = Mesh.Utils.createRenderableState(meshProps);
    const meshRenderObject = createRenderObject('mesh', meshValues, meshState, -1);
    const meshRepr = Representation.fromRenderObject('mesh', meshRenderObject);

    canvas3d.add(meshRepr);
    canvas3d.requestCameraReset();
}

init();