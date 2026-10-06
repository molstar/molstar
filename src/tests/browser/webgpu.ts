/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { WebGPUImagePass } from '../../mol-canvas3d/passes/webgpu-image';
import { fromHalfFloat, toHalfFloat } from '../../mol-util/number-conversion';
import { JitterVectors } from '../../mol-canvas3d/passes/jitter';
import { MultiSampleParams } from '../../mol-canvas3d/passes/multi-sample';
import { Camera } from '../../mol-canvas3d/camera';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { GPUTextureUsage } from '../../mol-gl/webgpu/compat';
import { WebGPUTextureData } from '../../mol-gl/webgpu/texture-data';
import { DirectVolumeValues } from '../../mol-gl/renderable/direct-volume';
import { MeshValues } from '../../mol-gl/renderable/mesh';
import { packIntToRGBArray } from '../../mol-util/number-packing';
import { Sphere3D } from '../../mol-math/geometry';
import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Points } from '../../mol-geo/geometry/points/points';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';
import { CylindersBuilder } from '../../mol-geo/geometry/cylinders/cylinders-builder';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Text } from '../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../mol-geo/geometry/text/text-builder';
import { createTransform, updateTransformData } from '../../mol-geo/geometry/transform-data';
import { createRenderObject, GraphicsRenderObject } from '../../mol-gl/render-object';
import { RendererParams, RendererProps } from '../../mol-gl/renderer';
import { WebGPUContext } from '../../mol-gl/webgpu/context';
import { WebGPURenderer } from '../../mol-gl/webgpu/renderer';
import { WebGPURayPick } from '../../mol-gl/webgpu/ray-pick';
import { WebGPUIlluminationCompose } from '../../mol-gl/webgpu/illumination-compose';
import { IlluminationParams } from '../../mol-canvas3d/passes/illumination';
import { WebGPUTracing } from '../../mol-gl/webgpu/tracing';
import { TracingParams } from '../../mol-canvas3d/passes/tracing';
import { WebGPUDepthPyramid } from '../../mol-gl/webgpu/depth-pyramid';
import { createWebGPUGeometry, createWebGPUThemes, value } from '../../mol-gl/webgpu/geometry';
import { computeMarchingCubesTextureMeshWebGPU, WebGPUTextureMeshGeometry } from '../../mol-gl/webgpu/texture-mesh';
import { computeMarchingCubesWebGPU } from '../../mol-gl/webgpu/marching-cubes';
import { computeMarchingCubesMesh } from '../../mol-geo/util/marching-cubes/algorithm';
import { createWrappedTensor } from '../../mol-repr/volume/util';
import { CubeVertices } from '../../mol-geo/util/marching-cubes/tables';
import { StaticBasisAndOrbitals, CreateOrbitalVolume, CreateOrbitalDensityVolume, CreateOrbitalRepresentation3D } from '../../extensions/alpha-orbitals/transforms';
import { computeOrbitalGridWebGPU } from '../../extensions/alpha-orbitals/gpu/webgpu';
import { initCubeGrid, CubeGridComputationParams, AlphaOrbital } from '../../extensions/alpha-orbitals/data-model';
import { sphericalCollocation } from '../../extensions/alpha-orbitals/collocation';
import { createSphericalCollocationGrid } from '../../extensions/alpha-orbitals/orbitals';
import { createSphericalCollocationDensityGrid } from '../../extensions/alpha-orbitals/density';
import { calcMeshColorSmoothingWebGPU, calcMeshColorSmoothingTextureWebGPU } from '../../mol-gl/webgpu/color-smoothing';
import { calcMeshColorSmoothing, ColorSmoothingInput } from '../../mol-geo/geometry/mesh/color-smoothing';
import { GaussianDensityWebGPU } from '../../mol-gl/webgpu/gaussian-density';
import { GaussianDensityCPU } from '../../mol-math/geometry/gaussian-density/cpu';
import { OrderedSet } from '../../mol-data/int/ordered-set';
import { Box3D } from '../../mol-math/geometry/primitives/box3d';
import { CommonSurfaceParams, getStructureConformationAndRadius } from '../../mol-repr/structure/visual/util/common';
import { Mat3, Mat4, Vec2, Vec3, Vec4, Tensor } from '../../mol-math/linear-algebra';
import { createWebGPUHandleHelper, HandleGroup, HandleHelperParams, isHandleLoci } from '../../mol-canvas3d/helper/handle-helper';
import { createWebGPUPointerHelper } from '../../mol-canvas3d/helper/pointer-helper';
import { SsaoParams } from '../../mol-canvas3d/passes/ssao';
import { MarkingParams } from '../../mol-canvas3d/passes/marking';
import { BackgroundParams } from '../../mol-canvas3d/passes/background';
import { Asset } from '../../mol-util/assets';
import { FxaaParams } from '../../mol-canvas3d/passes/fxaa';
import { OutlineParams } from '../../mol-canvas3d/passes/outline';
import { CameraHelperAxis, CameraHelperParams, createWebGPUCameraHelper, isCameraAxesLoci } from '../../mol-canvas3d/helper/camera-helper';
import { DefaultStereoCameraProps } from '../../mol-canvas3d/camera/stereo';
import { WebGPUStereoCamera } from '../../mol-gl/webgpu/stereo';
import { BloomParams } from '../../mol-canvas3d/passes/bloom';
import { DofParams } from '../../mol-canvas3d/passes/dof';
import { CasParams } from '../../mol-canvas3d/passes/cas';
import { SmaaParams } from '../../mol-canvas3d/passes/smaa';
import { ShadowParams } from '../../mol-canvas3d/passes/shadow';
import { PostprocessingParams } from '../../mol-canvas3d/passes/postprocessing';
import { Clip } from '../../mol-util/clip';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueCell } from '../../mol-util/value-cell';
import { PluginContext } from '../../mol-plugin/context';
import { PluginConfig } from '../../mol-plugin/config';
import { DefaultPluginSpec, PluginSpec } from '../../mol-plugin/spec';
import { DebugHelpers } from '../../extensions/debug-helpers';
import { WebGPUDebugRegistry } from '../../mol-canvas3d/helper/debug-registry';
import { RuntimeContext } from '../../mol-task';
import { Substance } from '../../mol-theme/substance';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Emissive } from '../../mol-theme/emissive';
import { isEmptyLoci, EveryLoci } from '../../mol-model/loci';
import { now } from '../../mol-util/now';
import { Volume } from '../../mol-model/volume';
import { CustomProperties } from '../../mol-model/custom-property';
import { Theme } from '../../mol-theme/theme';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { MarkerAction } from '../../mol-util/marker-action';

function assert(condition: unknown, message: string): asserts condition {
    if (!condition) throw new Error(message);
}

function showImage(label: string, image: ImageData) {
    const figure = document.createElement('figure');
    figure.style.cssText = 'display:inline-block;vertical-align:top;margin:12px;color:white;font:14px sans-serif';
    const caption = document.createElement('figcaption');
    caption.textContent = label;
    const canvas = document.createElement('canvas');
    canvas.width = image.width; canvas.height = image.height;
    canvas.style.cssText = 'display:block;background:#333;margin-top:8px';
    canvas.getContext('2d')!.putImageData(image, 0, 0);
    figure.append(caption, canvas); document.body.appendChild(figure);
}

function showPixels(label: string, pixels: { width: number, height: number, array: Uint8Array }) {
    const array = new Uint8ClampedArray(pixels.array);
    for (let i = 0; i < array.length; i += 4) {
        if (array[i + 3] > 0 && array[i + 3] < 255) for (let c = 0; c < 3; c++) array[i + c] = array[i + c] * 255 / array[i + 3];
    }
    showImage(label, new ImageData(array, pixels.width, pixels.height));
}

async function verify() {
    const canvas = document.querySelector('canvas')!;
    canvas.style.cssText = 'position:absolute;visibility:hidden';
    const context = await WebGPUContext.create(canvas);
    const errors: string[] = [];
    context.errors.subscribe(error => errors.push(error.message));
    const renderer = await WebGPURenderer.create(context);
    const camera = new Camera({ position: Vec3.create(0, 0, 20), target: Vec3(), radius: 10, radiusMax: 10, fog: 0 }, { x: 0, y: 0, width: canvas.width, height: canvas.height });
    const props = PD.getDefaultValues(RendererParams);
    props.backgroundColor = Color(0x000000);
    const results: string[] = [];
    try {
        await verifyDepthPyramid(context);
        results.push('native SSAO depth pyramids, odd and single-pixel dimensions, independent transparent alpha, resizing and analytical mip readbacks');
        await verifyMarchingCubesCompute(context);
        results.push('native GPU-generated texture meshes, direct storage-buffer rendering, CPU theme mirrors, affine normals and native marching cubes, all 256 cube cases, CPU triangle/normal parity, negative iso-levels, axis orders, GPU prefix-scan compaction, periodic boundaries, floodfill views, cropped regions and batching');
        await verifyColorSmoothingCompute(context);
        results.push('native weighted color/overlay grid accumulation, CPU references, logical instances, transforms, strides, empty cells and batching');
        await verifyOrbitalCompute(context);
        results.push('native orbital and electron-density compute, L=0–4, spherical orders, contractions, cutoffs, occupancies, task integration and batching');
        await verifyGaussianCompute(context);
        results.push('native Gaussian compute grids, analytical density/atom IDs, sparse subsets, radius offsets, smoothness, spatial bins/full-scan equivalence, batching and CPU-reference comparison');
        camera.update();
        const mesh = Mesh.create(new Float32Array([-4, -4, 0, 4, -4, 0, 0, 4, 0]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array([7, 7, 7]), 3, 1);
        const meshProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
        const values = Mesh.Utils.createValuesSimple(mesh, meshProps, Color(0xff0000), 1);
        const object = createRenderObject('mesh', values, Mesh.Utils.createRenderableState(meshProps), -1);
        renderer.render([object], camera, props);
        let pixels = await renderer.readPixels();
        const offset = (128 * pixels.width + 128) * 4;
        assert(pixels.array[offset] > 240 && pixels.array[offset + 1] < 10, `Native triangle must be red at the center. ${errors.join('; ')}`);
        const picked = await renderer.pick(128, 128);
        assert(picked?.id.objectId === object.id && picked.id.groupId === 7 && picked.id.instanceId === 0, 'Integer picking must retain object, group and instance IDs.');
        assert(picked.depth > 0 && picked.depth < 1, 'Picking must return WebGPU depth.');
        renderer.render([object], camera, props, true);
        const lodReference = await renderer.readPixels(), lodThemeBuilds = renderer.themeResourceStats.builds;
        const cameraDistance = -camera.view[14] / camera.scale;
        ValueCell.update(values.uLod, Vec4.create(cameraDistance + 1, cameraDistance + 10, 0, 0));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => i % 4 !== 3 || v === 0) && !(await renderer.pick(128, 128)), 'Distance bounds must remove geometry colors and canonical picking.');
        ValueCell.update(values.uLod, Vec4.create(cameraDistance - 10, cameraDistance + 10, 0, 0));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === lodReference.array[i]), 'Zero-overlap distance bounds must retain the complete in-range frame.');
        ValueCell.update(values.uLod, Vec4.create(cameraDistance - 2, cameraDistance + 10, 4, 0));
        renderer.render([object], camera, props, true);
        const fadedLod = await renderer.readPixels();
        const bayerRows = [[1, 13, 4, 16], [9, 5, 12, 8], [3, 15, 2, 14], [11, 7, 10, 6]];
        let lodCovered = 0, lodRemoved = 0;
        for (let y = 0; y < 256; y++) for (let x = 0; x < 256; x++) {
            const o = (y * 256 + x) * 4;
            if (!lodReference.array[o + 3]) continue;
            const keep = bayerRows[(255 - y) % 4][x % 4] <= 8;
            assert(fadedLod.array[o + 3] === (keep ? lodReference.array[o + 3] : 0), 'Half-distance fade must match the independent four-by-four coverage mask.');
            if (keep) lodCovered++; else lodRemoved++;
        }
        assert(lodCovered > 100 && lodRemoved > 100, 'Distance fading must preserve and remove populated molecular pixels.');
        ValueCell.update(values.uLod, Vec4.create(cameraDistance - 10, cameraDistance + 2, 4, 0));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === fadedLod.array[i]), 'Far-distance overlap must produce the same analytical half fade as the near overlap.');
        const savedPickThreshold = props.pickingAlphaThreshold;
        props.pickingAlphaThreshold = 0.75;
        renderer.render([object], camera, props, true);
        assert(!(await renderer.pick(128, 128)), 'Distance fading below the picking threshold must suppress selection.');
        props.pickingAlphaThreshold = 0.25;
        renderer.render([object], camera, props, true);
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Distance fading above the picking threshold must retain the molecular group.');
        props.pickingAlphaThreshold = savedPickThreshold;
        ValueCell.update(values.uLod, Vec4.create(0, 0, 0, 0));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === lodReference.array[i]) && renderer.themeResourceStats.builds === lodThemeBuilds, 'Disabling distance LOD must restore the frame without rebuilding geometry themes.');
        camera.scale = 2; camera.update();
        const scaledDistance = -camera.view[14] / camera.scale;
        ValueCell.update(values.uLod, Vec4.create(0, scaledDistance - 1, 0, 0));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => i % 4 !== 3 || v === 0) && !(await renderer.pick(128, 128)), 'Far-distance exclusion must account for runtime camera scaling.');
        camera.scale = 1; camera.update();
        Object.assign(camera.viewport, { x: 7, y: 9, width: 242, height: 238 });
        for (const jitter of [[0, 0], [0.375, -0.375]]) {
            camera.viewOffset.enabled = true;
            Camera.setViewOffset(camera.viewOffset, 242, 238, jitter[0], jitter[1], 242, 238); camera.update();
            ValueCell.update(values.uLod, Vec4.create(0, 0, 0, 0));
            renderer.render([object], camera, props, true);
            const viewportReference = await renderer.readPixels();
            ValueCell.update(values.uLod, Vec4.create(cameraDistance - 2, cameraDistance + 10, 4, 0));
            renderer.render([object], camera, props, true);
            const viewportFade = await renderer.readPixels();
            for (let y = 0; y < 256; y++) for (let x = 0; x < 256; x++) {
                const o = (y * 256 + x) * 4;
                if (!viewportReference.array[o + 3]) continue;
                const column = ((Math.floor(x + 0.5 + jitter[0] * 4) % 4) + 4) % 4;
                const row = ((Math.floor(255 - y + 0.5 + jitter[1] * 4) % 4) + 4) % 4;
                assert(viewportFade.array[o + 3] === (bayerRows[row][column] <= 8 ? viewportReference.array[o + 3] : 0), 'Offset viewports and signed sampling jitter must retain the independently computed distance coverage mask.');
            }
        }
        camera.viewOffset.enabled = false;
        Object.assign(camera.viewport, { x: 0, y: 0, width: 256, height: 256 }); camera.update();
        ValueCell.update(values.uLod, Vec4.create(0, 0, 0, 0));
        renderer.render([object], camera, props);
        results.push('native distance bounds, zero-overlap visibility, analytical ordered fading, picking thresholds, camera scale, offset viewports, signed sampling jitter and exact restoration');
        await verifyGeometryResources(renderer, camera, object, props);
        results.push('geometry and transform buffer retention across theme/material updates, changed transforms/positions, ValueCell replacement and failed resource rebuilds');
        await verifyTextureResources(renderer, camera, object, props);
        results.push('unchanged image/palette/spatial texture retention, source-revision invalidation, independent uploads, atomic failed rebuilds and canonical picking');
        await verifyIncrementalBufferSizes(renderer, camera, props);
        renderer.render([object], camera, props);
        results.push('incremental GPU buffer updates, retained capacity, packed byte-channel shrinking/clearing, growth and stable picking');
        const rayPicker = new WebGPURayPick(context);
        try {
            for (const origin of [Vec3.create(0, 0, 20), Vec3.create(10, 0, 20)]) {
                const ray = Ray3D.targetTo(Ray3D(), Ray3D.create(origin, Vec3()), Vec3());
                const hit = await rayPicker.pick(ray, camera, [object], props, 0);
                assert(hit?.id.objectId === object.id && hit.id.groupId === 7 && hit.id.instanceId === 0, 'Native straight and oblique rays must retain mesh object/group/instance IDs.');
                assert(Vec3.magnitude(hit.position) < 0.0001, 'Native mesh rays must return the analytical triangle intersection.');
            }
            object.state.pickable = false;
            const ray = Ray3D.create(Vec3.create(0, 0, 20), Vec3.create(0, 0, -1));
            assert(await rayPicker.pick(ray, camera, [object], props, 3) === undefined, 'Ray queries must respect mesh pickability.');
            object.state.pickable = true;
            assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Ray queries must retain the displayed canonical selection target.');
            const preserved = await renderer.readPixels();
            assert(preserved.array.every((v, i) => v === pixels.array[i]), 'Mesh ray queries must retain the displayed color frame.');
        } finally { object.state.pickable = true; rayPicker.dispose(); }
        results.push('native ray-based mesh picking, analytical straight/oblique intersections, pickability and independent selection/color targets');
        results.push('mesh rendering, color, integer picking and depth');
        await verifyTracingInput(context, renderer, camera, props);
        await verifyTracingBounces(context, renderer, props);
        results.push('native illumination G-buffers, inward view normals, raw material color, density/emission, front/back thickness depth, fog independence, transparency exclusion, resizing, native ray tracing, primary emission, analytical accumulation and restart');
        const positions = new WebGPUTextureData(), normals = new WebGPUTextureData(), groups = new WebGPUTextureData();
        const positionArray = new Float32Array([-4, -4, 0, 1, 4, -4, 0, 1, 0, 4, 0, 1, 100, 100, 100, 1]);
        positions.load({ width: 2, height: 2, array: positionArray });
        normals.load({ width: 2, height: 2, array: new Float32Array([0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]) });
        const packedGroups = new Uint8Array(16);
        for (let i = 0; i < 3; i++) packIntToRGBArray(7, packedGroups, i * 4);
        groups.load({ width: 2, height: 2, array: packedGroups });
        const textureGeometry = TextureMesh.create(3, 8, positions, groups, normals, Sphere3D.create(Vec3(), 6));
        const textureProps = { ...PD.getDefaultValues(TextureMesh.Params), ignoreLight: true, doubleSided: true };
        const textureValues = TextureMesh.Utils.createValuesSimple(textureGeometry, textureProps, Color(0xff0000), 1);
        const textureObject = createRenderObject('texture-mesh', textureValues, TextureMesh.Utils.createRenderableState(textureProps), -1);
        renderer.render([textureObject], camera, props);
        const texturePixels = await renderer.readPixels();
        assert(texturePixels.array.every((v, i) => v === pixels.array[i]), 'Native texture-mesh geometry must match ordinary triangle mesh rendering.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Packed texture-mesh group IDs must remain pickable.');
        const colors = new Uint8Array(8 * 3); colors.set([0, 255, 0], 7 * 3);
        ValueCell.update(textureValues.dColorType, 'group'); ValueCell.update(textureValues.tColor, { width: 8, height: 1, array: colors });
        renderer.render([textureObject], camera, props);
        const themedTexture = await renderer.readPixels();
        assert(themedTexture.array[offset + 1] > 240 && themedTexture.array[offset] < 10, 'Texture-mesh group themes must decode packed IDs.');
        const textureTransforms = new Float32Array(32); textureTransforms.set(Mat4.fromTranslation(Mat4(), Vec3.create(100, 0, 0)), 0); textureTransforms.set(Mat4.identity(), 16);
        ValueCell.update(textureValues.aTransform, textureTransforms); ValueCell.update(textureValues.aInstance, new Float32Array([0, 1])); ValueCell.update(textureValues.instanceCount, 2);
        renderer.render([textureObject], camera, props);
        assert((await renderer.pick(128, 128))?.id.instanceId === 1, 'Texture-mesh picking must retain native instance identities.');
        const moved = new Float32Array(positionArray); for (let i = 0; i < 3; i++) moved[i * 4] += 100;
        positions.load({ width: 2, height: 2, array: moved }); ValueCell.update(textureValues.tPosition, positions);
        renderer.render([textureObject], camera, props);
        assert(!(await renderer.pick(128, 128)), 'Texture-mesh position updates must rebuild native geometry.');
        positions.load({ width: 2, height: 2, array: positionArray }); ValueCell.update(textureValues.tPosition, positions);
        groups.load({ width: 2, height: 2, array: Float32Array.from(packedGroups, v => v / 255) }); ValueCell.update(textureValues.tGroup, groups);
        renderer.render([textureObject], camera, props);
        const restoredTexture = await renderer.readPixels();
        assert(restoredTexture.array.every((v, i) => v === themedTexture.array[i]) && (await renderer.pick(128, 128))?.id.groupId === 7, 'Normalized float groups and geometry restoration must preserve native texture meshes.');
        positions.destroy(); normals.destroy(); groups.destroy();
        renderer.render([object], camera, props);
        results.push('native texture meshes, packed byte/float groups, themes, instance picking and geometry updates');
        const originalSnapshot = camera.getSnapshot();
        for (const mode of ['perspective', 'orthographic'] as const) {
            camera.setState({ mode }, 0); camera.scale = 1; camera.update();
            renderer.render([object], camera, props);
            const baseline = await renderer.readPixels();
            for (const scale of [0.25, 2]) {
                camera.scale = scale; camera.update(); renderer.render([object], camera, props);
                const scaled = await renderer.readPixels();
                assert(scaled.array.every((v, i) => v === baseline.array[i]), 'Camera scale must preserve molecular geometry coverage in both projection modes.');
                assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Camera scale must preserve native picking IDs.');
                const pass = new WebGPUImagePass(context, camera, () => [object], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: false, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' }, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false } });
                const exported = await pass.getImageData(RuntimeContext.Synchronous, canvas.width, canvas.height);
                assert(exported.data.every((v, i) => v === scaled.array[i]), 'Screenshot cameras must preserve runtime camera scale.');
                await pass.dispose();
            }
        }
        camera.scale = 1; camera.setState(originalSnapshot, 0); camera.update();
        ValueCell.update(values.dIgnoreLight, false);
        Mat4.fromRotation(camera.headRotation, Math.PI / 3, Vec3.create(0, 1, 0));
        renderer.render([object], camera, props);
        const rotatedLighting = await renderer.readPixels();
        const rotatedPass = new WebGPUImagePass(context, camera, () => [object], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: false, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' }, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false } });
        const rotatedExport = await rotatedPass.getImageData(RuntimeContext.Synchronous, canvas.width, canvas.height);
        assert(rotatedExport.data.every((v, i) => v === rotatedLighting.array[i]), 'Exports must preserve head-rotated lighting.');
        await rotatedPass.dispose();
        Mat4.copy(camera.headRotation, Mat4.zero());
        renderer.render([object], camera, props);
        const unrotatedLighting = await renderer.readPixels();
        assert(rotatedLighting.array.some((v, i) => v !== unrotatedLighting.array[i]), 'Head rotation must change native light directions.');
        ValueCell.update(values.dIgnoreLight, true);
        renderer.render([object], camera, props);
        results.push('camera scale in perspective/orthographic rendering, picking and screenshots; head-rotated lighting exports');
        await verifyMultiSample(renderer, context, camera, object, props);
        results.push('native weighted supersampling, temporal convergence, all sample levels, canonical picking, premultiplied alpha, camera restoration and screenshot sampling');
        await verifyIlluminationPipeline(context, renderer, camera, object, props);
        results.push('native progressive illumination in the renderer, canonical picking, idle retention, material updates, full screenshot iteration completion and disable restoration');
        await verifyIlluminationSamples(context, renderer, camera, object, props);
        results.push('progressive illumination supersampling, all jitter levels, weighted GPU frame references, canonical picking, camera restoration, cropping and full screenshot convergence');
        await verifyIlluminationEffects(context, renderer, camera, object, props);
        results.push('illumination effect combinations, outline suppression, transparent SSAO, absence of duplicate opaque SSAO/shadows, picking and screenshot exports');
        await verifyCameraAxes(context, renderer, camera, props);
        results.push('native camera axes, all four corners, orientation, display scaling, labels, highlighting, picking, illumination exclusion and screenshot exports');
        await verifyStereo(context, renderer, camera, object, props);
        results.push('native stereo eye composition, odd/offset viewports, independent full/temporal sampling, per-eye picking and world positions, stereo updates and mono restoration');
        await verifyWorldHelpers(context, renderer, camera, props);
        results.push('native handle/pointer helpers, transformed picking, highlighting, occlusion, theme updates, display sizing, screenshots and illumination exclusion');
        showPixels('Native mesh', pixels);
        const post = PD.getDefaultValues(PostprocessingParams);
        post.occlusion = { name: 'off', params: {} }; post.bloom = { name: 'off', params: {} };
        post.antialiasing = { name: 'fxaa', params: PD.getDefaultValues(FxaaParams) };
        renderer.render([object], camera, props, false, 1, undefined, post);
        const smoothed = await renderer.readPixels();
        assert(smoothed.array.some((v, i) => i % 4 === 0 && v > 0 && v < 240), 'FXAA must produce partial coverage along geometry edges.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Postprocessing must preserve original integer picking.');
        showPixels('Native FXAA edges', smoothed);
        post.antialiasing = { name: 'smaa', params: { edgeThreshold: 0.1, maxSearchSteps: 16 } };
        renderer.render([object], camera, props, false, 1, undefined, post);
        const smaa = await renderer.readPixels();
        assert(smaa.array.some((v, i) => i % 4 === 0 && v > 0 && v < 240), 'SMAA must smooth geometric silhouette edges.');
        assert(smaa.array.some((v, i) => v !== smoothed.array[i]), 'SMAA must use its own morphology and lookup tables.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'SMAA must preserve the original picking attachment.');
        showPixels('Native SMAA edges', smaa);
        renderer.render([object], camera, props, true, 1, undefined, post); const smaaTransparent = await renderer.readPixels();
        assert(smaaTransparent.array.some((v, i) => i % 4 === 3 && v > 0 && v < 240), 'SMAA must filter transparent silhouette coverage.');
        for (let i = 0; i < smaaTransparent.array.length; i += 4) assert(smaaTransparent.array[i] <= smaaTransparent.array[i + 3], 'SMAA colors must remain premultiplied by the filtered alpha.');
        ValueCell.update(values.uColor, Vec3.create(0.1, 0.1, 0.1));
        post.antialiasing.params.edgeThreshold = 0.15;
        renderer.render([object], camera, props, true, 1, undefined, post); const belowThreshold = await renderer.readPixels();
        renderer.render([object], camera, props, true); const dimBaseline = await renderer.readPixels();
        assert(belowThreshold.array.every((v, i) => v === dimBaseline.array[i]), 'SMAA edge threshold must leave low contrast silhouettes unchanged.');
        post.antialiasing.params.edgeThreshold = 0.05;
        renderer.render([object], camera, props, true, 1, undefined, post); const lowThreshold = await renderer.readPixels();
        assert(lowThreshold.array.some((v, i) => v !== dimBaseline.array[i]), 'Lowering SMAA threshold must detect low contrast edges.');
        post.antialiasing.params.maxSearchSteps = 0;
        renderer.render([object], camera, props, true, 1, undefined, post); const zeroSearch = await renderer.readPixels();
        assert(zeroSearch.array[offset + 3] === 255 && zeroSearch.array[3] === 0, 'Zero SMAA search steps must produce valid alpha.');
        ValueCell.update(values.uColor, Vec3.create(1, 0, 0));
        post.antialiasing.params.maxSearchSteps = 16;
        Object.assign(camera.viewport, { x: 32, y: 24, width: 192, height: 160 }); camera.update();
        renderer.render([object], camera, props, true, 1, undefined, post); const offsetSmaa = await renderer.readPixels();
        assert(offsetSmaa.array[3] === 0 && (await renderer.pick(128, 152))?.id.groupId === 7, 'SMAA must respect offset viewport coordinates.');
        Object.assign(camera.viewport, { x: 0, y: 0, width: 256, height: 256 }); camera.update();
        post.antialiasing = { name: 'fxaa', params: PD.getDefaultValues(FxaaParams) };
        results.push('native three-pass SMAA, color threshold, search steps, premultiplied alpha, offset viewports and picking');
        ValueCell.update(values.uColor, Vec3.create(0.5, 0.5, 0.5));
        renderer.render([object], camera, props, false, 1, undefined, post);
        const beforeSharpen = await renderer.readPixels();
        post.sharpening = { name: 'on', params: { sharpness: 1, denoise: false } };
        renderer.render([object], camera, props, false, 1, undefined, post);
        const sharpened = await renderer.readPixels();
        assert(sharpened.array.some((v, i) => v !== beforeSharpen.array[i]), 'Native sharpening must adjust smoothed edges.');
        post.sharpening.params.denoise = true;
        renderer.render([object], camera, props, false, 1, undefined, post);
        const denoised = await renderer.readPixels();
        assert(denoised.array.some((v, i) => v !== sharpened.array[i]), 'Sharpening denoise controls must affect the filtered image.');
        showPixels('Native sharpening', denoised);
        post.sharpening = { name: 'off', params: {} };
        ValueCell.update(values.uColor, Vec3.create(1, 0, 0));
        post.outline = { name: 'on', params: { ...PD.getDefaultValues(OutlineParams), color: Color(0x0088ff) } };
        post.antialiasing = { name: 'off', params: {} };
        renderer.render([object], camera, props, false, 1, undefined, post);
        const outlined = await renderer.readPixels();
        let outsideEdge = -1;
        for (let i = 0; i < outlined.array.length; i += 4) if (outlined.array[i + 2] > 200 && pixels.array[i] === 0) { outsideEdge = i / 4; break; }
        assert(outsideEdge >= 0, 'Depth outlines must add the configured outline color around geometry.');
        assert(!(await renderer.pick(outsideEdge % 256, Math.floor(outsideEdge / 256))), 'Outline pixels outside geometry must remain unpickable.');
        showPixels('Native depth outline', outlined);
        post.enabled = false;
        renderer.render([object], camera, props, false, 1, undefined, post);
        const disabled = await renderer.readPixels();
        assert(disabled.array.every((v, i) => v === pixels.array[i]), 'Disabling postprocessing must restore the unprocessed frame.');
        post.enabled = true; post.outline.params.includeTransparent = false;
        ValueCell.update(values.alpha, 0.5);
        renderer.render([object], camera, props, true, 1, undefined, post);
        const transparentExcluded = await renderer.readPixels();
        assert(!transparentExcluded.array.some((v, i) => i % 4 === 2 && v > 100), 'Transparent objects must be excluded from outlines when requested.');
        post.outline.params.includeTransparent = true;
        renderer.render([object], camera, props, true, 1, undefined, post);
        const transparentIncluded = await renderer.readPixels();
        assert(transparentIncluded.array.some((v, i) => i % 4 === 2 && v > 100 && transparentIncluded.array[i + 1] === 255), 'Transparent outline alpha must double and clamp the source opacity, as in the existing compositor.');
        const outlineBackground = props.backgroundColor; props.backgroundColor = Color(0x204060);
        const background = Color.toRgb(props.backgroundColor);
        for (const alpha of [0.1, 0.25, 0.5]) {
            ValueCell.update(values.alpha, alpha);
            renderer.render([object], camera, props, true, 1, undefined, post); const clearOutline = await renderer.readPixels();
            renderer.render([object], camera, props, false, 1, undefined, post); const solidOutline = await renderer.readPixels();
            for (let i = 0; i < solidOutline.array.length; i += 4) {
                for (let c = 0; c < 3; c++) assert(Math.abs(solidOutline.array[i + c] - (clearOutline.array[i + c] + background[c] * (1 - clearOutline.array[i + 3] / 255))) <= 2, 'Transparent outlines over opaque backgrounds must match premultiplied compositing of the transparent frame.');
                assert(solidOutline.array[i + 3] === 255, 'Opaque background outlines must preserve canvas opacity.');
            }
            assert(!(await renderer.pick(outsideEdge % 256, Math.floor(outsideEdge / 256))), 'Transparent outline coverage must remain unpickable outside the source surface.');
            if (alpha === 0.25) showPixels('Transparent outline over opaque background', solidOutline);
            post.bloom = { name: 'on', params: { mode: 'emissive', strength: 1, radius: 0, threshold: 0, transparency: false } };
            renderer.render([object], camera, props, false, 1, undefined, post); const noTransparentEmission = await renderer.readPixels();
            assert(noTransparentEmission.array.every((v, i) => v === solidOutline.array[i]), 'Excluding transparent bloom must retain background-independent outline coverage.');
            post.bloom = { name: 'off', params: {} };
        }
        props.backgroundColor = outlineBackground;
        for (const opacity of [0.25, 1]) {
            ValueCell.update(values.alpha, opacity);
            renderer.render([object], camera, props, true, 1, undefined, post); const noFogOutline = await renderer.readPixels();
            const originalOutlineFog = camera.state.fog;
            for (const fog of [25, 75]) {
                camera.setState({ fog }, 0); camera.update();
                renderer.render([object], camera, props, true, 1, undefined, post); const foggedOutline = await renderer.readPixels();
                const distance = Math.abs(Vec3.transformMat4(Vec3(), Vec3(), camera.view)[2]);
                const t = Math.max(0, Math.min(1, (distance - camera.fogNear) / (camera.fogFar - camera.fogNear)));
                const expected = Math.round(noFogOutline.array[outsideEdge * 4 + 3] * (1 - t * t * (3 - 2 * t)));
                assert(Math.abs(foggedOutline.array[outsideEdge * 4 + 3] - expected) <= 1, 'Transparent outline alpha must follow the camera fog smoothstep at its source view depth.');
                assert(!(await renderer.pick(outsideEdge % 256, Math.floor(outsideEdge / 256))), 'Fogged outline pixels must remain outside the selectable geometry.');
            }
            camera.setState({ fog: originalOutlineFog }, 0); camera.update();
            renderer.render([object], camera, props, true, 1, undefined, post);
            assert((await renderer.readPixels()).array.every((v, i) => v === noFogOutline.array[i]), 'Disabling outline fog must restore the exact original frame.');
        }
        ValueCell.update(values.alpha, 1);
        await verifyOutlineLayers(renderer, camera, props);
        await verifyMarking(renderer, camera, props);
        await verifyBackground(renderer, camera, props);
        results.push('native gradient/image/cubemap backgrounds, analytical colors/alpha, coverage, aspect cropping, blur mips, rotation, fog, picking, asset replacement and screenshot exports');
        results.push('native selection/hover edges, depth occlusion, ghost opacity, edge thickness, highlight priority, fog, offset viewports, picking and screenshot exports');
        results.push('opaque/transparent outline depth separation, overlap compositing, nearest source alpha and independent picking');
        results.push('native FXAA, outlines, transparent outlines, effect disable and picking preservation');
        results.push('native contrast-adaptive sharpening and denoise controls');
        post.outline = { name: 'off', params: {} };
        post.dof = { name: 'on', params: { blurSize: 9, blurSpread: 1, inFocus: 0, PPM: 1, center: 'camera-target', mode: 'plane' } };
        renderer.render([object], camera, props, false, 1, undefined, post);
        const focused = await renderer.readPixels();
        assert(focused.array[offset] > 240, 'DOF must preserve the focused center plane.');
        post.dof.params.inFocus = 10;
        renderer.render([object], camera, props, false, 1, undefined, post);
        const defocused = await renderer.readPixels();
        assert(defocused.array.some((v, i) => v !== focused.array[i]), 'Changing focus distance must blur out-of-focus edges.');
        showPixels('Native depth of field', defocused);
        post.dof.params.mode = 'sphere'; post.dof.params.center = 'scene-center'; post.dof.params.inFocus = 0;
        renderer.render([object], camera, props, true, 1, undefined, post);
        const spherical = await renderer.readPixels();
        assert(spherical.array[3] === 0 && spherical.array[offset + 3] === 255, 'Spherical DOF must preserve transparent background and geometry alpha.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'DOF must retain geometry picking identifiers.');
        post.dof.params.PPM = 0;
        renderer.render([object], camera, props, true, 1, undefined, post);
        const zeroRange = await renderer.readPixels();
        assert(zeroRange.array[offset + 3] === 255, 'Zero focus range must produce valid pixels.');
        results.push('native planar/spherical DOF, scene-center focus, alpha and picking preservation');
        post.dof = { name: 'off', params: {} };
        renderer.render([object], camera, props, true); const beforeBloom = await renderer.readPixels();
        post.bloom = { name: 'on', params: { mode: 'luminosity', strength: 1, radius: 0, threshold: 0.1, transparency: true } };
        renderer.render([object], camera, props, true, 1, undefined, post);
        const glowing = await renderer.readPixels();
        let halo = -1;
        for (let i = 0; i < glowing.array.length; i += 4) if (beforeBloom.array[i + 3] === 0 && glowing.array[i] > 5 && glowing.array[i + 3] > 0) { halo = i / 4; break; }
        assert(halo >= 0, 'Luminosity bloom must spread color into a transparent halo.');
        assert(!(await renderer.pick(halo % 256, Math.floor(halo / 256))), 'Bloom halos must not become pickable geometry.');
        showPixels('Native luminosity bloom', glowing);
        post.bloom.params.radius = 1;
        renderer.render([object], camera, props, true, 1, undefined, post); const wider = await renderer.readPixels();
        assert(wider.array.some((v, i) => v !== glowing.array[i]), 'Bloom radius must change the mixture of five blur levels.');
        post.bloom.params.threshold = 1;
        renderer.render([object], camera, props, true, 1, undefined, post); const thresholded = await renderer.readPixels();
        assert(thresholded.array.every((v, i) => v === beforeBloom.array[i]), 'Luminosity threshold must reject dimmer geometry.');
        post.bloom.params.mode = 'emissive';
        renderer.render([object], camera, props, true, 1, undefined, post); const noEmission = await renderer.readPixels();
        assert(noEmission.array.every((v, i) => v === beforeBloom.array[i]), 'Non-emissive geometry must not seed emissive bloom.');
        ValueCell.update(values.uEmissive, 1);
        renderer.render([object], camera, props, true, 1, undefined, post); const emitted = await renderer.readPixels();
        assert(emitted.array.some((v, i) => v !== beforeBloom.array[i]), 'Uniform emissive material must seed native bloom independently of the luminosity threshold.');
        showPixels('Native emissive bloom', emitted);
        ValueCell.update(values.uColor, Vec3.create(0.2, 0.2, 0.2)); ValueCell.update(values.dIgnoreLight, false);
        const originalLights = props.light, originalAmbient = props.ambientIntensity;
        props.light = props.light.map(l => ({ ...l, intensity: 0 })); props.ambientIntensity = 0;
        renderer.render([object], camera, props, true, 1, undefined, post); const darkLighting = await renderer.readPixels();
        props.light = props.light.map(l => ({ ...l, intensity: 10 }));
        renderer.render([object], camera, props, true, 1, undefined, post); const brightLighting = await renderer.readPixels();
        assert(darkLighting.array[offset] !== brightLighting.array[offset], 'Lighting must change the visible material.');
        for (let i = 0; i < brightLighting.array.length; i += 4) if (beforeBloom.array[i + 3] === 0) {
            assert(darkLighting.array[i] === brightLighting.array[i], 'Emissive halo color must remain independent of material lighting.');
        }
        props.light = originalLights; props.ambientIntensity = originalAmbient;
        ValueCell.update(values.uColor, Vec3.create(1, 0, 0)); ValueCell.update(values.dIgnoreLight, true);
        camera.setState({ fog: 100 }, 0); camera.update();
        renderer.render([object], camera, props, true, 1, undefined, post); const foggedEmission = await renderer.readPixels();
        assert(foggedEmission.array.some((v, i) => i % 4 === 0 && beforeBloom.array[i + 3] === 0 && emitted.array[i] > 5 && v < emitted.array[i]), 'Fog must attenuate emitted light before bloom spreads it into the halo.');
        camera.setState({ fog: 0 }, 0); camera.update();
        ValueCell.update(values.alpha, 0.5); post.bloom.params.transparency = false;
        renderer.render([object], camera, props, true); const translucentBaseline = await renderer.readPixels();
        renderer.render([object], camera, props, true, 1, undefined, post); const translucentExcluded = await renderer.readPixels();
        assert(translucentExcluded.array.every((v, i) => v === translucentBaseline.array[i]), 'Disabling transparent emissive bloom must exclude transparent geometry.');
        post.bloom.params.transparency = true;
        renderer.render([object], camera, props, true, 1, undefined, post); const translucentGlow = await renderer.readPixels();
        assert(translucentGlow.array.some((v, i) => v !== translucentBaseline.array[i]), 'Transparent emission must participate when enabled.');
        post.bloom.params.strength = 0;
        renderer.render([object], camera, props, true, 1, undefined, post); const noGlow = await renderer.readPixels();
        assert(noGlow.array.every((v, i) => v === translucentBaseline.array[i]), 'Zero bloom strength must restore the source frame.');
        ValueCell.update(values.alpha, 1); post.bloom.params.strength = 1;
        Object.assign(camera.viewport, { x: 32, y: 24, width: 192, height: 160 }); camera.update();
        renderer.render([object], camera, props, true, 1, undefined, post); const framedGlow = await renderer.readPixels();
        assert(framedGlow.array[3] === 0 && (await renderer.pick(128, 152))?.id.groupId === 7, 'Bloom must respect offset viewports without changing geometry picking.');
        Object.assign(camera.viewport, { x: 0, y: 0, width: 256, height: 256 }); camera.update();
        ValueCell.update(values.uEmissive, 0); ValueCell.update(values.alpha, 1);
        results.push('five-level luminosity/emissive bloom, lighting independence, fog, radius, threshold, transparent emission, alpha and picking');

        const savedLighting = { light: props.light, ambientColor: props.ambientColor, ambientIntensity: props.ambientIntensity, exposure: props.exposure, celSteps: props.celSteps };
        ValueCell.update(values.uColor, Vec3.create(0.5, 0.5, 0.5)); ValueCell.update(values.dIgnoreLight, false);
        ValueCell.update(values.uMetalness, 0); ValueCell.update(values.uRoughness, 1);
        props.ambientIntensity = 0; props.exposure = 1;
        props.light = [{ inclination: 180, azimuth: 0, color: Color(0xff0000), intensity: 1 }];
        renderer.render([object], camera, props, true); const redLight = await renderer.readPixels();
        assert(redLight.array[offset] >= 128 && redLight.array[offset] <= 132 && redLight.array[offset + 1] < 5, 'Head-on red light must reproduce the matte Lambert/GGX reference intensity and color.');
        props.light = [{ ...props.light[0], inclination: 0 }];
        renderer.render([object], camera, props, true); const reverseLight = await renderer.readPixels();
        assert(reverseLight.array[offset] < 5, 'Moving a light behind the surface must remove its direct illumination.');
        props.light = [{ ...props.light[0], inclination: 180 }, { ...props.light[0], inclination: 180, color: Color(0x00ff00) }];
        renderer.render([object], camera, props, true); const twoLights = await renderer.readPixels();
        assert(twoLights.array[offset] > 120 && twoLights.array[offset + 1] > 120 && twoLights.array[offset + 2] < 5, 'Multiple colored directional lights must all illuminate the material.');
        props.light = []; props.ambientColor = Color(0x0000ff); props.ambientIntensity = 0.5;
        renderer.render([object], camera, props, true); const ambient = await renderer.readPixels();
        assert(ambient.array[offset] < 5 && ambient.array[offset + 2] >= 63 && ambient.array[offset + 2] <= 65, 'Ambient light must use its configured color and intensity.');
        props.exposure = 2;
        renderer.render([object], camera, props, true); const exposed = await renderer.readPixels();
        assert(Math.abs(exposed.array[offset + 2] - ambient.array[offset + 2] * 2) <= 1 && exposed.array[offset + 3] === ambient.array[offset + 3], 'Exposure must scale material brightness while preserving alpha.');
        props.light = [{ inclination: 180, azimuth: 0, color: Color(0xffffff), intensity: 1 }]; props.ambientIntensity = 0; props.exposure = 1;
        ValueCell.update(values.uMetalness, 1);
        renderer.render([object], camera, props, true); const roughMetal = await renderer.readPixels();
        ValueCell.update(values.uRoughness, 0.1);
        renderer.render([object], camera, props, true); const glossyMetal = await renderer.readPixels();
        assert(roughMetal.array[offset] < 40 && glossyMetal.array[offset] > 240, 'Roughness and metalness must control the GGX specular highlight.');
        showPixels('Native glossy metal', glossyMetal);
        ValueCell.update(values.uMetalness, 0); ValueCell.update(values.uRoughness, 1);
        const substance = new Uint8Array(32); substance.set([255, 0, 0, 255], 28);
        ValueCell.update(values.tSubstance, { array: substance, width: 8, height: 1 });
        ValueCell.update(values.dSubstance, true);
        renderer.render([object], camera, props, true); const substanceMetal = await renderer.readPixels();
        ValueCell.update(values.dSubstance, false);
        ValueCell.update(values.uMetalness, 0.99); ValueCell.update(values.uRoughness, 0.01);
        renderer.render([object], camera, props, true); const substanceReference = await renderer.readPixels();
        assert(substanceMetal.array.every((v, i) => Math.abs(v - substanceReference.array[i]) <= 1), 'Substance overlays must reproduce the clamped material mixture in native GGX lighting.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Material overlays must preserve atom/group picking.');
        ValueCell.update(values.uMetalness, 0); ValueCell.update(values.uRoughness, 1);
        ValueCell.update(values.dSubstance, true); ValueCell.update(values.uSubstanceStrength, 0.5);
        renderer.render([object], camera, props, true); const partialSubstance = await renderer.readPixels();
        ValueCell.update(values.dSubstance, false);
        ValueCell.update(values.uMetalness, 0.25); ValueCell.update(values.uRoughness, 0.5);
        renderer.render([object], camera, props, true); const partialReference = await renderer.readPixels();
        assert(partialSubstance.array.every((v, i) => Math.abs(v - partialReference.array[i]) <= 1), 'Partial substance strength must preserve the existing pre-mix and interpolation behavior.');
        ValueCell.update(values.uMetalness, 0); ValueCell.update(values.uRoughness, 1);
        ValueCell.update(values.dSubstance, true); ValueCell.update(values.uSubstanceStrength, 0);
        renderer.render([object], camera, props, true); const disabledSubstance = await renderer.readPixels();
        ValueCell.update(values.dSubstance, false); ValueCell.update(values.uSubstanceStrength, 1);
        renderer.render([object], camera, props, true); const matteReference = await renderer.readPixels();
        assert(disabledSubstance.array.every((v, i) => v === matteReference.array[i]), 'Zero substance strength must restore the original material.');
        results.push('native substance material overlays, strength updates, clamped mixing and picking preservation');
        const spatialMaterial = new WebGPUTextureData(), spatialBytes = new Uint8Array(24 * 24 * 4);
        for (let z = 0; z < 3; z++) for (let y = 0; y < 12; y++) for (let x = 0; x < 12; x++) {
            const index = ((Math.floor(z / 2) * 12 + y) * 24 + z % 2 * 12 + x) * 4;
            spatialBytes.set([40 + x * 8 + y * 2 + z * 16, 180 - x * 5 + z * 10, 0, 128 + x * 4], index);
        }
        spatialMaterial.load({ width: 24, height: 24, array: spatialBytes });
        ValueCell.update(values.tSubstanceGrid, spatialMaterial);
        ValueCell.update(values.uSubstanceGridDim, Vec3.create(12, 12, 3));
        ValueCell.update(values.uSubstanceGridTransform, Vec4.create(-6.5, -6.5, -1.5, 1));
        ValueCell.update(values.dSubstance, true);
        for (const translation of [0, 1]) for (const scale of [0.5, 1, 2]) {
            const matrix = Mat4.fromTranslation(Mat4(), Vec3.create(translation, 0, 0));
            ValueCell.update(values.aTransform, new Float32Array(matrix)); camera.scale = scale; camera.update();
            ValueCell.update(values.dSubstanceType, 'volumeInstance');
            renderer.render([object], camera, props, true); const spatialPixels = await renderer.readPixels();
            const reference = new Uint8Array(12);
            for (let i = 0; i < 3; i++) {
                const x = mesh.vertexBuffer.ref.value[i * 3] + 6 + translation, y = mesh.vertexBuffer.ref.value[i * 3 + 1] + 6;
                reference.set([40 + x * 8 + y * 2 + 24, 195 - x * 5, 0, 128 + x * 4], i * 4);
            }
            ValueCell.update(values.dSubstanceType, 'vertexInstance');
            ValueCell.update(values.tSubstance, { array: reference, width: 3, height: 1 });
            renderer.render([object], camera, props, true); const referencePixels = await renderer.readPixels();
            assert(spatialPixels.array.every((v, i) => Math.abs(v - referencePixels.array[i]) <= 1), 'Spatial material sampling must match bilinear slice/z interpolation in molecular coordinates, including instance transforms and camera scale.');
            assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Spatial material grids must preserve group picking.');
            if (translation === 0 && scale === 1) showPixels('Native spatial material', spatialPixels);
        }
        camera.scale = 1; camera.update(); ValueCell.update(values.aTransform, new Float32Array(Mat4.identity()));
        ValueCell.update(values.dSubstance, false); ValueCell.update(values.dSubstanceType, 'groupInstance');
        spatialMaterial.destroy();
        renderer.render([object], camera, props, true); const afterSpatial = await renderer.readPixels();
        assert(afterSpatial.array.every((v, i) => v === matteReference.array[i]), 'Disabling spatial materials must restore the original frame after grid disposal.');
        results.push('native spatial material grids, slice interpolation, instance transforms, camera scale, picking and restoration');
        const spatialOverlay = new WebGPUTextureData(), overpaintBytes = new Uint8Array(spatialBytes);
        for (let i = 3; i < overpaintBytes.length; i += 4) overpaintBytes[i] = 255;
        spatialOverlay.load({ width: 24, height: 24, array: overpaintBytes });
        const scalarOverlay = new WebGPUTextureData(); scalarOverlay.load({ width: 24, height: 24, array: spatialBytes });
        ValueCell.update(values.tOverpaintGrid, spatialOverlay); ValueCell.update(values.tTransparencyGrid, scalarOverlay); ValueCell.update(values.tEmissiveGrid, scalarOverlay);
        ValueCell.update(values.uOverpaintGridDim, Vec3.create(12, 12, 3)); ValueCell.update(values.uTransparencyGridDim, Vec3.create(12, 12, 3)); ValueCell.update(values.uEmissiveGridDim, Vec3.create(12, 12, 3));
        ValueCell.update(values.uOverpaintGridTransform, Vec4.create(-6.5, -6.5, -1.5, 1)); ValueCell.update(values.uTransparencyGridTransform, Vec4.create(-6.5, -6.5, -1.5, 1)); ValueCell.update(values.uEmissiveGridTransform, Vec4.create(-6.5, -6.5, -1.5, 1));
        ValueCell.update(values.dOverpaint, true); ValueCell.update(values.dTransparency, true); ValueCell.update(values.dEmissive, true);
        const overlayPickThreshold = props.pickingAlphaThreshold; props.pickingAlphaThreshold = 0.1;
        for (const scale of [0.5, 1, 2]) {
            camera.scale = scale; camera.update();
            ValueCell.update(values.dOverpaintType, 'volumeInstance'); ValueCell.update(values.dTransparencyType, 'volumeInstance'); ValueCell.update(values.dEmissiveType, 'volumeInstance');
            renderer.render([object], camera, props, true); const spatialPixels = await renderer.readPixels();
            const rgba = new Uint8Array(12), scalar = new Uint8Array(3);
            for (let i = 0; i < 3; i++) {
                const x = mesh.vertexBuffer.ref.value[i * 3] + 6, y = mesh.vertexBuffer.ref.value[i * 3 + 1] + 6;
                rgba.set([40 + x * 8 + y * 2 + 24, 195 - x * 5, 0, 255], i * 4); scalar[i] = 128 + x * 4;
            }
            ValueCell.update(values.dOverpaintType, 'vertexInstance'); ValueCell.update(values.dTransparencyType, 'vertexInstance'); ValueCell.update(values.dEmissiveType, 'vertexInstance');
            ValueCell.update(values.tOverpaint, { array: rgba, width: 3, height: 1 }); ValueCell.update(values.tTransparency, { array: scalar, width: 3, height: 1 }); ValueCell.update(values.tEmissive, { array: scalar, width: 3, height: 1 });
            renderer.render([object], camera, props, true); const referencePixels = await renderer.readPixels();
            assert(spatialPixels.array.every((v, i) => Math.abs(v - referencePixels.array[i]) <= 1), 'Spatial overpaint/transparency/emission must match the sampled per-vertex reference and retain transparent draw classification.');
            assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Combined spatial overlays must preserve selection IDs.');
            if (scale === 1) showPixels('Native spatial overlays', spatialPixels);
        }
        ValueCell.update(values.dTransparency, false); ValueCell.update(values.dEmissive, false);
        spatialOverlay.load({ width: 1, height: 1, array: new Uint8Array([0, 0, 255, 128]) });
        ValueCell.update(values.tOverpaintGrid, spatialOverlay); ValueCell.update(values.uOverpaintGridDim, Vec3.create(1, 1, 1));
        ValueCell.update(values.dOverpaintType, 'volumeInstance'); ValueCell.update(values.uOverpaintStrength, 0.5);
        camera.scale = 1; camera.update(); renderer.render([object], camera, props, true); const partialOverpaint = await renderer.readPixels();
        const overpaintAlpha = 128 / 255, overpaintWeight = overpaintAlpha * 0.5;
        ValueCell.update(values.dOverpaint, false);
        ValueCell.update(values.uColor, Vec3.create(0.5 * (1 - overpaintWeight) + 0.5 * (1 - overpaintAlpha) * 0.5 * overpaintWeight,
            0.5 * (1 - overpaintWeight) + 0.5 * (1 - overpaintAlpha) * 0.5 * overpaintWeight,
            0.5 * (1 - overpaintWeight) + (0.5 * (1 - overpaintAlpha) + overpaintAlpha) * 0.5 * overpaintWeight));
        renderer.render([object], camera, props, true); const overpaintReference = await renderer.readPixels();
        assert(partialOverpaint.array.every((v, i) => Math.abs(v - overpaintReference.array[i]) <= 1), 'Partial spatial overpaint strength must preserve pre-mixing before fragment material blending.');
        ValueCell.update(values.uColor, Vec3.create(0.5, 0.5, 0.5)); ValueCell.update(values.uOverpaintStrength, 1);
        ValueCell.update(values.dOverpaintType, 'groupInstance'); ValueCell.update(values.dTransparencyType, 'groupInstance'); ValueCell.update(values.dEmissiveType, 'groupInstance');
        camera.scale = 1; camera.update(); spatialOverlay.destroy(); scalarOverlay.destroy();
        props.pickingAlphaThreshold = overlayPickThreshold;
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === matteReference.array[i]), 'Clearing spatial overlays must restore the original frame.');
        results.push('native spatial overpaint/transparency/emission, combined interpolation, camera scale, transparent draw classification, picking and restoration');
        const spatialColors = new WebGPUTextureData(); spatialColors.load({ width: 24, height: 24, array: spatialBytes });
        ValueCell.update(values.tColorGrid, spatialColors); ValueCell.update(values.uColorGridDim, Vec3.create(12, 12, 3));
        ValueCell.update(values.uColorGridTransform, Vec4.create(-6.5, -6.5, -1.5, 1)); ValueCell.update(values.dIgnoreLight, true);
        for (const type of ['volume', 'volumeInstance']) for (const scale of [0.5, 1, 2]) {
            ValueCell.update(values.aTransform, new Float32Array(Mat4.fromTranslation(Mat4(), Vec3.create(1, 0, 0))));
            camera.scale = scale; camera.update(); ValueCell.update(values.dColorType, type);
            renderer.render([object], camera, props, true); const gridFrame = await renderer.readPixels();
            const reference = new Uint8Array(9);
            for (let i = 0; i < 3; i++) {
                const x = mesh.vertexBuffer.ref.value[i * 3] + 6 + (type === 'volumeInstance' ? 1 : 0), y = mesh.vertexBuffer.ref.value[i * 3 + 1] + 6;
                reference.set([40 + x * 8 + y * 2 + 24, 195 - x * 5, 0], i * 3);
            }
            ValueCell.update(values.dColorType, 'vertex'); ValueCell.update(values.tColor, { array: reference, width: 3, height: 1 });
            renderer.render([object], camera, props, true); const vertexFrame = await renderer.readPixels();
            assert(gridFrame.array.every((v, i) => Math.abs(v - vertexFrame.array[i]) <= 1), 'Invariant and instance color grids must match sampled vertex colors across camera scales.');
            if (type === 'volumeInstance' && scale === 1) showPixels('Native spatial color grid', gridFrame);
        }
        camera.scale = 1; camera.update(); ValueCell.update(values.aTransform, new Float32Array(Mat4.identity()));
        ValueCell.update(values.dColorType, 'volumeInstance'); ValueCell.update(values.dOverpaint, true);
        const colorOverpaint = new Uint8Array(32); colorOverpaint.set([255, 0, 128, 255], 28);
        ValueCell.update(values.tOverpaint, { array: colorOverpaint, width: 8, height: 1 });
        renderer.render([object], camera, props, true); const colorPainted = await renderer.readPixels();
        ValueCell.update(values.dOverpaint, false); ValueCell.update(values.dColorType, 'uniform'); ValueCell.update(values.uColor, Vec3.create(1, 0, 128 / 255));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - colorPainted.array[i]) <= 1), 'Group overpaint must be applied after native spatial color sampling.');
        spatialColors.load({ width: 1, height: 1, array: new Uint8Array([128, 0, 0, 255]) });
        ValueCell.update(values.tColorGrid, spatialColors); ValueCell.update(values.uColorGridDim, Vec3.create(1, 1, 1));
        for (const filter of ['nearest', 'linear'] as const) {
            ValueCell.update(values.dColorType, 'volume'); ValueCell.update(values.dUsePalette, true);
            ValueCell.update(values.tPalette, { array: new Uint8Array([255, 0, 0, 0, 255, 0, 0, 0, 255, 255, 255, 0]), width: 4, height: 1, filter });
            renderer.render([object], camera, props, true); const paletteFrame = await renderer.readPixels();
            ValueCell.update(values.dColorType, 'uniform'); ValueCell.update(values.dUsePalette, false);
            ValueCell.update(values.uColor, filter === 'nearest' ? Vec3.create(0, 0, 1) : Vec3.create(0, 0.5, 0.5));
            renderer.render([object], camera, props, true);
            assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - paletteFrame.array[i]) <= 1), 'Native spatial palette colors must preserve nearest/linear texture filtering.');
        }
        spatialColors.destroy(); ValueCell.update(values.uColor, Vec3.create(0.5, 0.5, 0.5)); ValueCell.update(values.dIgnoreLight, false);
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === matteReference.array[i]), 'Disabling spatial colors must restore the original material frame.');
        results.push('native invariant/instance color grids, camera scale, group overpaint, nearest/linear palettes and restoration');
        const encodedColors = new Uint8Array(9), encodedScalars: number[] = [];
        for (let i = 0; i < 3; i++) {
            const scalar = [0.125, 0.875, 0.5][i], packed = Math.round(scalar * 16777214);
            packIntToRGBArray(packed, encodedColors, i * 3); encodedScalars.push(packed / 16777214);
        }
        ValueCell.update(values.dColorType, 'vertex'); ValueCell.update(values.tColor, { array: encodedColors, width: 3, height: 1 });
        ValueCell.update(values.dUsePalette, true); ValueCell.update(values.dIgnoreLight, true);
        const projected = [0, 1, 2].map(i => {
            const p = camera.project(Vec4(), Vec3.fromArray(Vec3(), mesh.vertexBuffer.ref.value, i * 3));
            return [p[0], 256 - p[1]];
        });
        const [[ax, ay], [bx, by], [cx, cy]] = projected;
        const determinant = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
        for (const filter of ['nearest', 'linear'] as const) {
            ValueCell.update(values.tPalette, { array: new Uint8Array([0, 0, 0, 255, 255, 255]), width: 2, height: 1, filter });
            renderer.render([object], camera, props, true); const interpolatedPalette = await renderer.readPixels();
            for (const [x, y] of [[110, 140], [145, 150], [128, 95]]) {
                const a = ((by - cy) * (x + 0.5 - cx) + (cx - bx) * (y + 0.5 - cy)) / determinant;
                const b = ((cy - ay) * (x + 0.5 - cx) + (ax - cx) * (y + 0.5 - cy)) / determinant;
                const scalar = a * encodedScalars[0] + b * encodedScalars[1] + (1 - a - b) * encodedScalars[2];
                const expected = filter === 'nearest' ? (scalar >= 0.5 ? 255 : 0) : Math.round(Math.max(0, Math.min(1, scalar * 2 - 0.5)) * 255);
                const offset = (y * interpolatedPalette.width + x) * 4;
                assert(Math.abs(interpolatedPalette.array[offset] - expected) <= 1 && interpolatedPalette.array[offset + 3] === 255, 'Palette lookup must happen after interpolation, matching analytical pixel barycentric values and filtering.');
            }
            assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Fragment palette lookup must preserve group picking.');
            if (filter === 'linear') showPixels('Native interpolated palette', interpolatedPalette);
        }
        ValueCell.update(values.dUsePalette, false); ValueCell.update(values.dColorType, 'uniform'); ValueCell.update(values.uColor, Vec3.create(0.5, 0.5, 0.5));
        const partialPaint = new Uint8Array(32); partialPaint.set([0, 0, 255, 128], 28);
        ValueCell.update(values.dOverpaint, true); ValueCell.update(values.tOverpaint, { array: partialPaint, width: 8, height: 1 }); ValueCell.update(values.uOverpaintStrength, 0.5);
        renderer.render([object], camera, props, true); const groupPartial = await renderer.readPixels();
        ValueCell.update(values.dOverpaint, false);
        ValueCell.update(values.uColor, Vec3.create(0.5 * (1 - overpaintWeight) + 0.5 * (1 - overpaintAlpha) * 0.5 * overpaintWeight,
            0.5 * (1 - overpaintWeight) + 0.5 * (1 - overpaintAlpha) * 0.5 * overpaintWeight,
            0.5 * (1 - overpaintWeight) + (0.5 * (1 - overpaintAlpha) + overpaintAlpha) * 0.5 * overpaintWeight));
        renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - groupPartial.array[i]) <= 1), 'Group overpaint must preserve the same pre-mixing/fragment strength formula as spatial grids.');
        ValueCell.update(values.uColor, Vec3.create(0.5, 0.5, 0.5)); ValueCell.update(values.uOverpaintStrength, 1); ValueCell.update(values.dIgnoreLight, false);
        results.push('ordinary fragment palette interpolation, analytical barycentric references, nearest/linear filters, picking and partial group overpaint');
        const cachedEmission = new Uint8Array(8); cachedEmission[7] = 255;
        ValueCell.update(values.tEmissive, { array: cachedEmission, width: 8, height: 1 });
        ValueCell.update(values.dEmissiveType, 'groupInstance'); ValueCell.update(values.dEmissive, true); ValueCell.update(values.uEmissiveStrength, 0.5);
        renderer.render([object], camera, props, true); const enabledEmission = await renderer.readPixels();
        assert(enabledEmission.array.some((v, i) => v !== matteReference.array[i]), 'A cached group-emission overlay must visibly affect surface brightness when enabled.');
        ValueCell.update(values.dEmissive, false); renderer.render([object], camera, props, true); const disabledEmission = await renderer.readPixels();
        assert(disabledEmission.array.every((v, i) => v === matteReference.array[i]), 'Disabling group emission must restore the frame while retaining its cached texture data.');
        assert(values.tEmissive.ref.value.array === cachedEmission, 'Emission enable/disable must not require clearing its cached data.');
        ValueCell.update(values.dEmissive, true); renderer.render([object], camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === enabledEmission.array[i]), 'Re-enabling cached emission must reproduce the same surface brightness.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Emission toggles must preserve surface picking.');
        ValueCell.update(values.dEmissive, false); ValueCell.update(values.uEmissiveStrength, 1);
        results.push('cached group emission enable/disable/re-enable, frame restoration and picking');
        ValueCell.update(values.uMetalness, 0); ValueCell.update(values.uRoughness, 1); ValueCell.update(values.dCelShaded, true);
        props.celSteps = 5;
        renderer.render([object], camera, props, true); const celFive = await renderer.readPixels();
        props.celSteps = 2;
        renderer.render([object], camera, props, true); const celTwo = await renderer.readPixels();
        assert(celTwo.array[offset] > celFive.array[offset] + 20, 'Cel shading must quantize lighting using the configured step count.');
        ValueCell.update(values.dCelShaded, false);
        const originalNormals = values.aNormal.ref.value;
        ValueCell.update(values.aNormal, new Float32Array([0.8, 0, 0.6, 0.8, 0, 0.6, 0.8, 0, 0.6]));
        renderer.render([object], camera, props, true); const smoothNormal = await renderer.readPixels();
        ValueCell.update(values.dFlatShaded, true);
        renderer.render([object], camera, props, true); const flatNormal = await renderer.readPixels();
        assert(flatNormal.array[offset] > smoothNormal.array[offset] + 20, 'Flat shading must reconstruct face normals instead of using vertex normals.');
        ValueCell.update(values.aNormal, originalNormals); ValueCell.update(values.dFlatShaded, false);
        renderer.render([object], camera, props, true); const bumpBaseline = await renderer.readPixels();
        ValueCell.update(values.dFlipSided, true);
        renderer.render([object], camera, props, true); const flipped = await renderer.readPixels();
        assert(flipped.array[offset] < 5, 'Flipping a double-sided surface must invert its lighting normal.');
        ValueCell.update(values.uDoubleSided, false);
        renderer.render([object], camera, props, true); const frontCulled = await renderer.readPixels();
        assert(frontCulled.array[offset + 3] === 0 && !(await renderer.pick(128, 128)), 'Flipped single-sided surfaces must cull their front faces and picking fragments.');
        ValueCell.update(values.dFlipSided, false); ValueCell.update(values.uDoubleSided, true);
        ValueCell.update(values.uBumpiness, 1); ValueCell.update(values.uBumpFrequency, 1); ValueCell.update(values.uBumpAmplitude, 0.5);
        renderer.render([object], camera, props, true); const bumped = await renderer.readPixels();
        assert(bumped.array.some((v, i) => v !== bumpBaseline.array[i]), 'Bump frequency, amplitude and material bumpiness must perturb native shading.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Material effects must preserve molecular geometry picking.');
        showPixels('Native procedural bump shading', bumped);
        ValueCell.update(values.uBumpiness, 0); ValueCell.update(values.uBumpFrequency, 0); ValueCell.update(values.uBumpAmplitude, 0);
        ValueCell.update(values.dIgnoreLight, true); Object.assign(props, savedLighting);
        results.push('native GGX materials, colored directional/ambient lights, exposure, metalness, roughness, cel/flat/flipped shading and procedural bumps');

        ValueCell.update(values.uColor, Vec3.create(1, 0, 0));
        ValueCell.update(values.aNormal, new Float32Array([0.8, 0, 0.6, 0.8, 0, 0.6, 0.8, 0, 0.6]));
        const savedFalloff = props.xrayEdgeFalloff; props.xrayEdgeFalloff = 1;
        ValueCell.update(values.dXrayShaded, 'on');
        renderer.render([object], camera, props, true); const xray = await renderer.readPixels();
        assert(Math.abs(xray.array[offset + 3] - 102) <= 1, 'X-ray opacity must be one minus the normal-facing factor.');
        assert(xray.array[offset] === xray.array[offset + 3], 'X-ray colors must remain premultiplied by their opacity.');
        ValueCell.update(values.dXrayShaded, 'inverted');
        renderer.render([object], camera, props, true); const invertedXray = await renderer.readPixels();
        assert(Math.abs(invertedXray.array[offset + 3] - 153) <= 1, 'Inverted x-ray opacity must use the normal-facing factor.');
        props.xrayEdgeFalloff = 2; ValueCell.update(values.dXrayShaded, 'on');
        renderer.render([object], camera, props, true); const steepXray = await renderer.readPixels();
        assert(steepXray.array[offset + 3] > xray.array[offset + 3] + 50, 'X-ray edge falloff must change the opacity curve.');
        showPixels('Native x-ray surface', steepXray);
        ValueCell.update(values.dXrayShaded, 'off'); ValueCell.update(values.aNormal, originalNormals); props.xrayEdgeFalloff = savedFalloff;
        const originalElements = values.elements.ref.value;
        const originalInteriorColor = Vec4.clone(values.uInteriorColor.ref.value);
        const originalInteriorSubstance = Vec4.clone(values.uInteriorSubstance.ref.value);
        ValueCell.update(values.elements, new Uint32Array([0, 2, 1]));
        ValueCell.update(values.uInteriorColor, Vec4.create(0, 1, 0, 1));
        renderer.render([object], camera, props, true); const interior = await renderer.readPixels();
        assert(interior.array[offset + 1] > 240 && interior.array[offset] < 5, 'Back faces must use their configured interior color.');
        ValueCell.update(values.alpha, 0.5);
        renderer.render([object], camera, props, true); const hiddenBack = await renderer.readPixels();
        assert(hiddenBack.array[offset + 3] === 0, 'Transparent back faces must be discarded when disabled.');
        ValueCell.update(values.dTransparentBackfaces, 'on');
        renderer.render([object], camera, props, true); const translucentBack = await renderer.readPixels();
        assert(Math.abs(translucentBack.array[offset + 3] - 128) <= 1, 'Enabled transparent back faces must retain surface opacity.');
        ValueCell.update(values.dTransparentBackfaces, 'opaque');
        renderer.render([object], camera, props, true); const opaqueBack = await renderer.readPixels();
        assert(opaqueBack.array[offset + 3] === 255 && opaqueBack.array[offset + 1] > 240, 'Opaque back faces must override surface opacity.');
        showPixels('Native interior surface', opaqueBack);
        ValueCell.update(values.alpha, 1); ValueCell.update(values.dIgnoreLight, false);
        ValueCell.update(values.uInteriorSubstance, Vec4.create(0, 1, 0, 1));
        renderer.render([object], camera, props, true); const matteInterior = await renderer.readPixels();
        ValueCell.update(values.uInteriorSubstance, Vec4.create(1, 0.2, 0, 1));
        renderer.render([object], camera, props, true); const metalInterior = await renderer.readPixels();
        assert(metalInterior.array.some((v, i) => v !== matteInterior.array[i]), 'Interior material strength must blend its own roughness and metalness.');
        ValueCell.update(values.elements, originalElements); ValueCell.update(values.uInteriorColor, originalInteriorColor);
        ValueCell.update(values.uInteriorSubstance, originalInteriorSubstance); ValueCell.update(values.dTransparentBackfaces, 'off');
        ValueCell.update(values.dIgnoreLight, true);
        results.push('native x-ray/inverted opacity and falloff, interior colors/materials, and transparent back-face modes');

        const behindMesh = Mesh.create(new Float32Array([-4, -4, -2, 4, -4, -2, 0, 4, -2]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array([9, 9, 9]), 3, 1);
        const behind = createRenderObject('mesh', Mesh.Utils.createValuesSimple(behindMesh, meshProps, Color(0x0000ff), 1), Mesh.Utils.createRenderableState(meshProps), -1);
        const savedThreshold = props.pickingAlphaThreshold; props.pickingAlphaThreshold = 0.5;
        ValueCell.update(values.alpha, 0.2);
        renderer.render([object, behind], camera, props, true); const faint = await renderer.readPixels();
        assert(faint.array[offset] > 40 && faint.array[offset + 2] > 180, 'Faint foreground surfaces must remain visible over geometry behind them.');
        assert((await renderer.pick(128, 128))?.id.objectId === behind.id, 'Below-threshold foreground opacity must allow picking the geometry behind it.');
        props.pickingAlphaThreshold = 0.1;
        renderer.render([object, behind], camera, props, true); const lowerThreshold = await renderer.readPixels();
        assert((await renderer.pick(128, 128))?.id.objectId === object.id, 'Lowering picking opacity threshold must select the foreground surface.');
        assert(lowerThreshold.array.every((v, i) => v === faint.array[i]), 'Picking threshold must not change the color rendering.');
        object.state.pickable = false;
        renderer.render([object, behind], camera, props, true);
        assert((await renderer.pick(128, 128))?.id.objectId === behind.id, 'Unpickable geometry must not erase pickable geometry behind it.');
        object.state.pickable = true; object.state.colorOnly = true;
        renderer.render([object, behind], camera, props, true);
        assert((await renderer.pick(128, 128))?.id.objectId === behind.id, 'Color-only objects must be excluded from selection.');
        object.state.colorOnly = false; ValueCell.update(values.alpha, 1); props.pickingAlphaThreshold = 0.5;
        ValueCell.update(values.aNormal, new Float32Array([0.8, 0, 0.6, 0.8, 0, 0.6, 0.8, 0, 0.6])); ValueCell.update(values.dXrayShaded, 'on');
        renderer.render([object, behind], camera, props, true);
        assert((await renderer.pick(128, 128))?.id.objectId === behind.id, 'X-ray facing opacity must participate in selection thresholding.');
        ValueCell.update(values.dXrayShaded, 'inverted');
        renderer.render([object, behind], camera, props, true);
        assert((await renderer.pick(128, 128))?.id.objectId === object.id, 'Inverted x-ray opacity above threshold must remain selectable.');
        ValueCell.update(values.dXrayShaded, 'off'); ValueCell.update(values.aNormal, originalNormals); props.pickingAlphaThreshold = savedThreshold;
        const layeredMesh = Mesh.create(new Float32Array([-4, -4, 1, 4, -4, 1, 0, 4, 1, -4, -4, -1, 4, -4, -1, 0, 4, -1]), new Uint32Array([0, 1, 2, 3, 4, 5]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array([3, 3, 3, 4, 4, 4]), 6, 2);
        const layeredProps = { ...meshProps, alpha: 0.8 };
        const layered = createRenderObject('mesh', Mesh.Utils.createValuesSimple(layeredMesh, layeredProps, Color(0xff0000), 1), Mesh.Utils.createRenderableState(layeredProps), -1);
        renderer.render([layered], camera, props, true);
        assert((await renderer.pick(128, 128))?.id.groupId === 3, 'Selection must choose the nearest transparent triangle even if a farther triangle is drawn last in the same object.');
        ValueCell.update(layered.values.alpha, 1); ValueCell.update(layered.values.uGroupCount, 5);
        ValueCell.update(layered.values.dColorType, 'group');
        const layerColors = new Uint8Array(15); layerColors[9] = 255; layerColors[14] = 255;
        ValueCell.update(layered.values.tColor, { array: layerColors, width: 5, height: 1 });
        ValueCell.update(layered.values.dTransparency, true);
        ValueCell.update(layered.values.tTransparency, { array: new Uint8Array([0, 0, 0, 0, 128]), width: 5, height: 1 });
        renderer.render([layered], camera, props, true); const mixedOpacity = await renderer.readPixels();
        assert(mixedOpacity.array[offset] > 240 && mixedOpacity.array[offset + 2] < 10, 'Opaque fragments inside a mixed-opacity object must occlude farther transparent fragments.');
        assert((await renderer.pick(128, 128))?.id.groupId === 3, 'Mixed opacity must retain nearest opaque group selection.');
        results.push('per-fragment opaque/transparent surface depth classification');
        results.push('independent opacity-threshold picking, selection through faint/unpickable objects, x-ray opacity and nearest transparent depth');
        await verifyWeightedTransparency(context, renderer, camera, props);
        results.push('weighted intersecting transparency, triangle/object order independence, analytic color/alpha, opaque occlusion, independent picking, mode updates, supersampling and screenshot parity');
        await verifyDepthPeeling(context, renderer, camera, props);
        results.push('native dual-depth peeling, analytic layered color/alpha, iteration limits, triangle order independence, nearest picking, opaque occlusion, supersampling and screenshot parity');

        ValueCell.update(values.uColor, Vec3.create(0, 1, 0));
        renderer.render([object], camera, props);
        pixels = await renderer.readPixels();
        assert(pixels.array[offset + 1] > 240 && pixels.array[offset] < 10, 'ValueCell updates must replace GPU color data.');
        results.push('dynamic color update');
        const clipDefaults = Clip.Params.objects.ctor();
        for (const type of ['plane', 'sphere', 'cube', 'cylinder', 'infiniteCone'] as const) {
            const clipObject = { ...clipDefaults, type, scale: Vec3.create(100, 100, 100), position: Vec3.create(0, 10, 0) };
            // Cone clips negative z; move its apex above the surface.
            if (type === 'infiniteCone') clipObject.position = Vec3.create(0, 0, 100);
            const clippedProps = { ...meshProps, clip: { variant: 'pixel' as const, objects: [clipObject] } };
            const clippedValues = Mesh.Utils.createValuesSimple(mesh, clippedProps, Color(0xff0000), 1);
            const clippedObject = createRenderObject('mesh', clippedValues, Mesh.Utils.createRenderableState(clippedProps), -1);
            renderer.render([clippedObject], camera, props);
            assert(!(await renderer.pick(128, 128)), `Native ${type} clipping must discard center fragments and picking.`);
            Mesh.Utils.updateValues(clippedValues, { ...clippedProps, clip: { ...clippedProps.clip, objects: [{ ...clipObject, invert: true }] } });
            renderer.render([clippedObject], camera, props);
            assert(await renderer.pick(128, 128), `Native ${type} clip inversion must preserve outside fragments.`);
            const outsidePosition = type === 'plane' ? Vec3.create(0, -1000, 0) : type === 'infiniteCone' ? Vec3.create(0, 0, -1000) : Vec3.create(1000, 0, 0);
            Mesh.Utils.updateValues(clippedValues, { ...clippedProps, clip: { ...clippedProps.clip, objects: [{ ...clipObject, position: outsidePosition }] } });
            renderer.render([clippedObject], camera, props);
            assert(await renderer.pick(128, 128), `Native ${type} clipping must retain fragments outside its shape.`);
            const transformed = Mat4.fromTranslation(Mat4(), Vec3.create(1000, 1000, 1000));
            Mesh.Utils.updateValues(clippedValues, { ...clippedProps, clip: { ...clippedProps.clip, objects: [{ ...clipObject, transform: transformed }] } });
            renderer.render([clippedObject], camera, props);
            assert(await renderer.pick(128, 128), 'Clip transforms must map sample points before the signed-distance test.');
            Mesh.Utils.updateValues(clippedValues, { ...clippedProps, clip: { ...clippedProps.clip, variant: 'instance' } });
            renderer.render([clippedObject], camera, props);
            assert(!(await renderer.pick(128, 128)), 'Instance clipping must remove the entire clipped instance.');
            ValueCell.update(clippedValues.dClipping, true); ValueCell.update(clippedValues.dClippingType, 'instance');
            ValueCell.update(clippedValues.tClipping, { array: new Uint8Array([1]), width: 1, height: 1 });
            renderer.render([clippedObject], camera, props);
            assert(await renderer.pick(128, 128), 'Clipping masks must exempt designated clip objects.');
        }
        results.push('five native clip shapes, inversion, instance clipping and masks');
        for (const mode of [0, 1]) {
            ValueCell.update(values.uWiggleAmplitude, 2); ValueCell.update(values.uWiggleMode, mode);
            renderer.setTime(0); renderer.render([object], camera, props); const before = await renderer.readPixels();
            renderer.setTime(1); renderer.render([object], camera, props); const after = await renderer.readPixels();
            assert(before.array.some((v, i) => v !== after.array[i]), 'Native wiggle animation must change pixels as time advances.');
        }
        ValueCell.update(values.uWiggleAmplitude, 0);
        ValueCell.update(values.dWiggle, true); ValueCell.update(values.dWiggleType, 'instance'); ValueCell.update(values.uWiggleStrength, 2);
        ValueCell.update(values.tWiggle, { array: new Uint8Array([255]), width: 1, height: 1 });
        renderer.setTime(0); renderer.render([object], camera, props); const overlayBefore = await renderer.readPixels();
        renderer.setTime(1); renderer.render([object], camera, props); const overlayAfter = await renderer.readPixels();
        assert(overlayBefore.array.some((v, i) => v !== overlayAfter.array[i]), 'Per-instance wiggle weights must animate independently of base amplitude.');
        ValueCell.update(values.dWiggle, false); ValueCell.update(values.uTumbleAmplitude, 8);
        renderer.setTime(0); renderer.render([object], camera, props); const tumbleBefore = await renderer.readPixels();
        renderer.setTime(1); renderer.render([object], camera, props); const tumbleAfter = await renderer.readPixels();
        assert(tumbleBefore.array.some((v, i) => v !== tumbleAfter.array[i]), 'Native instance tumble must change pixels as time advances.');
        props.enableAnimation = false;
        renderer.setTime(0); renderer.render([object], camera, props); const frozenBefore = await renderer.readPixels();
        renderer.setTime(1); renderer.render([object], camera, props); const frozenAfter = await renderer.readPixels();
        assert(frozenBefore.array.every((v, i) => v === frozenAfter.array[i]), 'Disabling renderer animation must freeze geometry.');
        props.enableAnimation = true;
        ValueCell.update(values.uTumbleAmplitude, 0);
        results.push('position/group wiggle, instance tumble and animation disable');

        const transforms = new Float32Array(32);
        transforms.set(Mat4.identity()); transforms.set(Mat4.identity(), 16);
        transforms[12] = -4; transforms[28] = 4;
        const instanceValues = Mesh.Utils.createValuesSimple(mesh, meshProps, Color(0x0088ff), 1, createTransform(transforms, 2, undefined, 0, 0));
        const instanced = createRenderObject('mesh', instanceValues, Mesh.Utils.createRenderableState(meshProps), -1);
        renderer.render([instanced], camera, props);
        const left = await renderer.pick(65, 128), right = await renderer.pick(190, 128);
        assert(left?.id.instanceId === 0 && right?.id.instanceId === 1, 'Instanced geometry must preserve distinct instance IDs.');
        results.push('instanced transforms and identity');
        showPixels('Instanced transforms', await renderer.readPixels());

        object.state.visible = false;
        renderer.render([object], camera, props);
        assert(await renderer.pick(128, 128) === undefined, 'Hidden objects must not be pickable.');
        results.push('visibility');
        object.state.visible = true;

        const points = Points.create(new Float32Array([0, 0, 0]), new Float32Array([3]), 1);
        const pointProps = { ...PD.getDefaultValues(Points.Params), pointSizeAttenuation: false, pointStyle: 'circle' as const };
        const pointValues = Points.Utils.createValuesSimple(points, pointProps, Color(0xffff00), 16);
        const pointObject = createRenderObject('points', pointValues, Points.Utils.createRenderableState(pointProps), -1);
        renderer.render([pointObject], camera, props);
        assert((await renderer.pick(128, 128))?.id.groupId === 3, 'Expanded WebGPU point sprites must render and pick.');
        results.push('point sprites');
        showPixels('Point sprite', await renderer.readPixels());

        const sphere = Spheres.create(new Float32Array([0, 0, 0]), new Float32Array([4]), 1);
        const sphereProps = { ...PD.getDefaultValues(Spheres.Params), ignoreLight: true };
        const sphereValues = Spheres.Utils.createValuesSimple(sphere, sphereProps, Color(0xff8800), 2);
        const sphereObject = createRenderObject('spheres', sphereValues, Spheres.Utils.createRenderableState(sphereProps), -1);
        renderer.render([sphereObject], camera, props);
        assert((await renderer.pick(128, 128))?.id.groupId === 4, 'Native spheres must render and pick.');
        for (const [x, y] of [[128, 128], [138, 128], [128, 138], [118, 118]]) {
            const hit = await renderer.pick(x, y);
            assert(hit?.id.groupId === 4, 'Compact sphere geometry must retain analytical interior picking.');
            const point = camera.unproject(Vec3(), Vec3.create(x + 0.5, 256 - y - 0.5, hit!.depth));
            assert(Math.abs(Vec3.magnitude(point) - 2) < 0.0001, 'Sphere depth must lie on the analytical radius rather than a tessellated face.');
        }
        renderer.render([sphereObject], camera, props, true);
        const sphereLodReference = await renderer.readPixels();
        ValueCell.update(sphereValues.uLod, Vec4.create(cameraDistance - 2, cameraDistance + 10, 4, 0));
        renderer.render([sphereObject], camera, props, true);
        const sphereLodFade = await renderer.readPixels();
        for (let y = 0; y < 256; y++) for (let x = 0; x < 256; x++) {
            const o = (y * 256 + x) * 4;
            if (!sphereLodReference.array[o + 3]) continue;
            assert(sphereLodFade.array[o + 3] === (bayerRows[(255 - y) % 4][x % 4] <= 8 ? sphereLodReference.array[o + 3] : 0), 'Sphere distance fading must use its center rather than vary across the tessellated surface.');
        }
        ValueCell.update(sphereValues.uLod, Vec4.create(0, 0, 0, 0));
        renderer.render([sphereObject], camera, props);
        results.push('native sphere-center distance fading with independently checked coverage');
        const sphereCamera = camera.getSnapshot(), ellipsoidTransforms = new Float32Array(32);
        const ellipsoids = [-4, 4].map(x => {
            const transform = Mat4.fromScaling(Mat4(), Vec3.create(1, 0.5, 1.5)); transform[12] = x; return transform;
        });
        ellipsoids.forEach((matrix, i) => ellipsoidTransforms.set(matrix, i * 16));
        const ellipsoidValues = Spheres.Utils.createValuesSimple(sphere, sphereProps, Color(0xff0000), 2, createTransform(ellipsoidTransforms, 2, undefined, 0, 0));
        ValueCell.update(ellipsoidValues.aInstance, new Float32Array([1, 0]));
        ValueCell.update(ellipsoidValues.dColorType, 'instance');
        ValueCell.update(ellipsoidValues.tColor, { array: new Uint8Array([255, 0, 0, 0, 0, 255]), width: 2, height: 1 });
        const ellipsoidObject = createRenderObject('spheres', ellipsoidValues, Spheres.Utils.createRenderableState(sphereProps), -1);
        for (const mode of ['perspective', 'orthographic'] as const) {
            camera.setState({ mode }, 0); camera.update(); renderer.render([ellipsoidObject], camera, props, true);
            const image = await renderer.readPixels();
            for (let i = 0; i < ellipsoids.length; i++) {
                const projected = camera.project(Vec4(), Vec3.create(ellipsoids[i][12], 0, 0));
                const x = Math.floor(projected[0]), y = Math.floor(256 - projected[1]), hit = await renderer.pick(x, y);
                assert(hit?.id.instanceId === 1 - i && hit.id.groupId === 4, 'Analytical ellipsoids must retain reordered logical instance and molecular group IDs.');
                const world = camera.unproject(Vec3(), Vec3.create(x + 0.5, 256 - y - 0.5, hit!.depth));
                const local = Vec3.transformMat4(Vec3(), world, Mat4.invert(Mat4(), ellipsoids[i]));
                assert(Math.abs(Vec3.magnitude(local) - 2) < 0.0001, 'Perspective and orthographic sphere depth must remain analytical under nonuniform instance transforms.');
                assert(image.array[(y * 256 + x) * 4 + (i === 0 ? 2 : 0)] > 240, 'Compact sphere colors must use reordered instance themes.');
            }
        }
        camera.setState(sphereCamera, 0); camera.update(); renderer.render([sphereObject], camera, props);
        results.push('compact analytical spheres, nonuniform ellipsoids, perspective/orthographic depth, reordered instance colors and molecular picking');
        await verifySphereLods(renderer, camera, props);
        await verifySphereLodInstances(renderer, camera, props);
        renderer.render([sphereObject], camera, props);
        results.push('native sphere distance levels, prepared subsets, reduced draws, analytical radius scaling, overlaps and exact restoration');
        results.push('native sphere instance-level culling, reordered grid/physical instances, reduced draws, colors, molecular picking and exact restoration');
        results.push('spheres');
        showPixels('Sphere', await renderer.readPixels());

        Spheres.Utils.updateValues(sphereValues, { ...sphereProps, alpha: 0.8, alphaThickness: 4 });
        renderer.render([sphereObject], camera, props, true);
        const thinSphere = await renderer.readPixels();
        assert(Math.abs(thinSphere.array[offset + 3] - 102) <= 1, 'Transparent spheres must scale alpha by radius divided by thickness.');
        assert((await renderer.pick(128, 128))?.id.groupId === 4, 'Sphere thickness must preserve theme-alpha selection thresholds.');
        Spheres.Utils.updateValues(sphereValues, { ...sphereProps, alpha: 0.8, alphaThickness: 1 });
        renderer.render([sphereObject], camera, props, true);
        assert(Math.abs((await renderer.readPixels()).array[offset + 3] - 204) <= 1, 'Sphere thickness must clamp its alpha multiplier to one.');
        Spheres.Utils.updateValues(sphereValues, { ...sphereProps, alphaThickness: 4 });
        renderer.render([sphereObject], camera, props, true);
        assert((await renderer.readPixels()).array[offset + 3] === 255, 'Opaque spheres must ignore alpha thickness.');
        results.push('sphere radius-dependent transparency updates with independent selection thresholds');

        const linesBuilder = LinesBuilder.create(4, 4);
        linesBuilder.add(-4, 0, 0, 4, 0, 0, 8);
        const linesProps = { ...PD.getDefaultValues(Lines.Params), lineSizeAttenuation: false };
        const lineObject = createRenderObject('lines', Lines.Utils.createValuesSimple(linesBuilder.getLines(), linesProps, Color(0xffffff), 4), Lines.Utils.createRenderableState(linesProps), -1);
        renderer.render([lineObject], camera, props);
        assert((await renderer.pick(128, 128))?.id.groupId === 8, 'Native thick lines must render and pick.');
        results.push('thick lines');
        showPixels('Thick line', await renderer.readPixels());

        const cylinderBuilder = CylindersBuilder.create(6, 6);
        cylinderBuilder.add(-4, 0, 0, 4, 0, 0, 1, true, true, 2, 5);
        const cylinderProps = { ...PD.getDefaultValues(Cylinders.Params), ignoreLight: true };
        const cylinderObject = createRenderObject('cylinders', Cylinders.Utils.createValuesSimple(cylinderBuilder.getCylinders(), cylinderProps, Color(0xff00ff), 0.5), Cylinders.Utils.createRenderableState(cylinderProps), -1);
        renderer.render([cylinderObject], camera, props);
        assert((await renderer.pick(128, 128))?.id.groupId === 5, 'Native cylinder meshes must render and pick.');
        results.push('cylinders');
        showPixels('Cylinder', await renderer.readPixels());
        const bondValues = cylinderObject.values;
        const bondColors = new Uint8Array(36); bondColors.set([0, 255, 0], 15); bondColors.set([255, 0, 0, 0, 0, 255], 30);
        ValueCell.update(bondValues.dColorType, 'group'); ValueCell.update(bondValues.dDualColor, true);
        ValueCell.update(bondValues.tColor, { array: bondColors, width: 12, height: 1 });
        ValueCell.update(bondValues.aColorMode, new Float32Array(6).fill(3));
        renderer.render([cylinderObject], camera, props, true); const gradientBond = await renderer.readPixels();
        assert(Math.abs(gradientBond.array[offset] - 191) <= 2 && Math.abs(gradientBond.array[offset + 2] - 64) <= 2, 'Interpolated half-bonds must blend toward the midpoint color.');
        assert(gradientBond.array[offset - 12 * 4] > gradientBond.array[offset + 12 * 4] && gradientBond.array[offset - 12 * 4 + 2] < gradientBond.array[offset + 12 * 4 + 2], 'Bond endpoint interpolation must vary along the cylinder axis.');
        assert((await renderer.pick(128, 128))?.id.groupId === 5, 'Bond color interpolation must preserve molecular group IDs.');
        showPixels('Interpolated bond', gradientBond);
        ValueCell.update(bondValues.aColorMode, new Float32Array(6).fill(1));
        renderer.render([cylinderObject], camera, props, true); const blueBond = await renderer.readPixels();
        assert(blueBond.array[offset] < 5 && blueBond.array[offset + 2] > 250, 'Fixed endpoint modes must use the second bond color.');
        ValueCell.update(bondValues.aColorMode, new Float32Array(6).fill(2));
        renderer.render([cylinderObject], camera, props, true); const singleColorBond = await renderer.readPixels();
        assert(singleColorBond.array[offset + 1] > 250 && singleColorBond.array[offset + 2] < 5, 'Single-color cylinder mode must retain undoubled group indexing.');
        results.push('dual bond colors, half-bond gradients, fixed endpoint modes and preserved picking');
        const unclippedCamera = camera.getSnapshot();
        try {
            for (const mode of ['perspective', 'orthographic'] as const) {
                camera.setState({ mode, radius: 0.25, fog: 0, clipFar: false }, 0); camera.update();
                assert(camera.near > 19 && camera.near < 20, 'Interior fixture must cut the sphere and cylinder at the near plane.');
                for (const [shape, values, update, group] of [
                    [sphereObject, sphereValues, (solidInterior: boolean, alpha: number, transparentBackfaces: 'off' | 'on' | 'opaque' = 'on') => Spheres.Utils.updateValues(sphereValues, { ...sphereProps, solidInterior, alpha, alphaThickness: 0, transparentBackfaces }), 4],
                    [cylinderObject, bondValues, (solidInterior: boolean, alpha: number, transparentBackfaces: 'off' | 'on' | 'opaque' = 'on') => Cylinders.Utils.updateValues(bondValues, { ...cylinderProps, solidInterior, alpha, transparentBackfaces }), 5]
                ] as const) {
                    update(true, 0.5);
                    ValueCell.update(values.uInteriorColor, Vec4.create(0, 0, 1, 1));
                    renderer.render([shape], camera, props, true);
                    const interior = await renderer.readPixels(), hit = await renderer.pick(128, 128);
                    assert(interior.array[offset] < 5 && interior.array[offset + 2] > 120 && Math.abs(interior.array[offset + 3] - 128) <= 1, `${mode} ${shape.type} solid interiors must fill near-plane cuts once with interior color and theme alpha: ${Array.from(interior.array.subarray(offset, offset + 4))}.`);
                    assert(hit?.id.groupId === group && hit.depth < 0.00001, `${mode} solid interiors must retain near-plane depth and molecular picking IDs.`);
                    for (const backfaces of ['off', 'opaque'] as const) {
                        update(true, 0.5, backfaces);
                        renderer.render([shape], camera, props, true);
                        assert((await renderer.readPixels()).array[offset + 3] === (backfaces === 'off' ? 0 : 255), 'Solid interiors must honor hidden and opaque transparent-backface settings.');
                        assert((await renderer.pick(128, 128))?.id.groupId === group, 'Backface color settings must preserve independent theme-alpha picking.');
                        if (backfaces === 'opaque') {
                            const input = renderer.renderTracingInput([shape], camera, props);
                            const front = new Float32Array((await context.readTexture(input.depth, 0, 0, input.depth.width, input.depth.height)).buffer)[128 * input.depth.width + 128];
                            const back = new Float32Array((await context.readTexture(input.backDepth, 0, 0, input.backDepth.width, input.backDepth.height)).buffer)[128 * input.backDepth.width + 128];
                            assert(front < 0.00001 && back > front, 'Opaque interiors of translucent impostors must retain illumination thickness.');
                        }
                    }
                    update(true, 1); ValueCell.update(values.uInteriorColor, Vec4.create(0, 0, 1, 1));
                    const input = renderer.renderTracingInput([shape], camera, props);
                    const depth = async (texture: GPUTexture) => new Float32Array((await context.readTexture(texture, 0, 0, texture.width, texture.height)).buffer)[128 * texture.width + 128];
                    const frontDepth = await depth(input.depth), backDepth = await depth(input.backDepth);
                    assert(frontDepth < 0.00001 && backDepth > frontDepth && backDepth < 1, `${mode} ${shape.type} clipped interiors must retain near/far tracing depth for physical thickness: ${frontDepth}, ${backDepth}.`);
                    const normal = Array.from(new Uint16Array((await context.readTexture(input.normal, 128, 128, 1, 1, 8)).buffer), fromHalfFloat);
                    const albedo = Array.from(new Uint16Array((await context.readTexture(input.albedo, 128, 128, 1, 1, 8)).buffer), fromHalfFloat);
                    assert(normal[2] < -0.999 && albedo[0] === 0 && albedo[2] === 1, 'Clipped interiors must supply their plane normal and interior material to illumination.');
                    renderer.render([object, shape], camera, props, true);
                    assert((await renderer.pick(128, 128))?.id.groupId === group, 'Near-plane interiors must occlude molecular geometry behind them.');
                    const pass = new WebGPUImagePass(context, camera, () => [shape], { renderer: props, transparentBackground: true, cameraHelper: { axes: { name: 'off', params: {} } }, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false }, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } });
                    try {
                        const exported = await pass.getImageData(RuntimeContext.Synchronous, canvas.width, canvas.height);
                        assert(exported.data[offset] === 0 && exported.data[offset + 2] === 255 && exported.data[offset + 3] === 255, 'Screenshot exports must retain solid clipped interiors and their material.');
                    } finally { await pass.dispose(); }
                    update(false, 1);
                    renderer.render([shape], camera, props, true);
                    assert((await renderer.readPixels()).array[offset + 3] === 0, 'Disabling solid interiors must leave a near-plane cut open.');
                }
                Spheres.Utils.updateValues(sphereValues, { ...sphereProps, solidInterior: true });
                const stretched = Mat4.fromScaling(Mat4(), Vec3.create(2, 0.5, 1));
                ValueCell.update(sphereValues.aTransform, new Float32Array(stretched));
                renderer.render([sphereObject], camera, props, true);
                for (const [point, inside] of [[Vec3.create(3, 0, 0.25), true], [Vec3.create(0, 1.3, 0.25), false]] as const) {
                    const screen = camera.project(Vec4(), point);
                    const hit = await renderer.pick(Math.floor(screen[0]), Math.floor(canvas.height - screen[1]));
                    assert(Boolean(hit) === inside, 'Near-plane sphere interiors must preserve nonuniform instance scaling.');
                }
                ValueCell.update(sphereValues.aTransform, new Float32Array(Mat4.identity()));
            }
        } finally {
            camera.setState(unclippedCamera, 0); camera.update();
            ValueCell.update(sphereValues.uInteriorColor, Vec4()); ValueCell.update(bondValues.uInteriorColor, Vec4());
            Spheres.Utils.updateValues(sphereValues, sphereProps); Cylinders.Utils.updateValues(bondValues, cylinderProps);
        }
        results.push('sphere/cylinder solid near-plane interiors, perspective/orthographic projection, interior color, single-layer transparency, depth and picking');
        const frontClip = { ...Clip.Params.objects.ctor(), type: 'plane' as const, rotation: { axis: Vec3.create(1, 0, 0), angle: -90 } };
        for (const [shape, values, update, group] of [
            [sphereObject, sphereValues, (solidInterior: boolean) => Spheres.Utils.updateValues(sphereValues, { ...sphereProps, solidInterior, clip: { variant: 'pixel', objects: [frontClip] } }), 4],
            [cylinderObject, bondValues, (solidInterior: boolean) => Cylinders.Utils.updateValues(bondValues, { ...cylinderProps, solidInterior, clip: { variant: 'pixel', objects: [frontClip] } }), 5]
        ] as const) {
            update(true); ValueCell.update(values.uInteriorColor, Vec4.create(0, 0, 1, 1));
            renderer.render([shape], camera, props, true);
            const interior = await renderer.readPixels();
            assert(interior.array[offset + 2] > 250 && interior.array[offset] < 5 && (await renderer.pick(128, 128))?.id.groupId === group, 'Pixel clipping must expose the retained back intersection with interior color and molecular picking.');
            update(false); renderer.render([shape], camera, props, true);
            assert((await renderer.readPixels()).array[offset + 3] === 0, 'Single-sided impostors without solid interiors must discard clipped front surfaces.');
        }
        Spheres.Utils.updateValues(sphereValues, sphereProps); Cylinders.Utils.updateValues(bondValues, cylinderProps);
        const openBondBuilder = CylindersBuilder.create(6, 6);
        openBondBuilder.add(-0.5, 0, 2, 0.5, 0, -2, 1, false, false, 2, 5);
        const openBondValues = Cylinders.Utils.createValuesSimple(openBondBuilder.getCylinders(), { ...cylinderProps, transparentBackfaces: 'on', alpha: 0.5 }, Color(0xff0000), 1);
        const openBond = createRenderObject('cylinders', openBondValues, Cylinders.Utils.createRenderableState(cylinderProps), -1);
        ValueCell.update(openBondValues.uInteriorColor, Vec4.create(0, 0, 1, 1));
        renderer.render([openBond], camera, props, true);
        const endInterior = await renderer.readPixels();
        assert(endInterior.array[offset] < 5 && endInterior.array[offset + 2] > 120 && Math.abs(endInterior.array[offset + 3] - 128) <= 1, 'Solid uncapped cylinder ends must render one interior layer along the bond axis.');
        assert((await renderer.pick(128, 128))?.id.groupId === 5, 'Solid uncapped bond ends must retain picking IDs.');
        Cylinders.Utils.updateValues(openBondValues, { ...cylinderProps, transparentBackfaces: 'on' });
        const openInput = renderer.renderTracingInput([openBond], camera, props);
        const openNormal = Array.from(new Uint16Array((await context.readTexture(openInput.normal, 128, 128, 1, 1, 8)).buffer), fromHalfFloat);
        assert(Math.abs(openNormal[0]) < 0.01 && openNormal[2] < -0.999, 'Artificial solid bond ends must use the inward view-ray normal for illumination.');
        Cylinders.Utils.updateValues(openBondValues, { ...cylinderProps, solidInterior: false });
        renderer.render([openBond], camera, props, true);
        assert(!(await renderer.pick(128, 128)), 'Disabling solid interiors must restore uncapped bond ends.');
        results.push('pixel-clipped sphere/bond back intersections, solid uncapped bond ends and single-layer end-on transparency');

        const textProps = { ...PD.getDefaultValues(Text.Params), attachment: 'middle-center' as const, background: false };
        const textBuilder = TextBuilder.create(textProps, 8, 8);
        textBuilder.add('C', 0, 0, 0, 0, 1, 6);
        const textObject = createRenderObject('text', Text.Utils.createValuesSimple(textBuilder.getText(), textProps, Color(0xffffff), 4), Text.Utils.createRenderableState(textProps), -1);
        renderer.render([textObject], camera, props);
        pixels = await renderer.readPixels();
        assert(pixels.array.some((value, i) => i % 4 === 0 && value > 200), 'Native SDF labels must produce visible glyph pixels.');
        const labelReference = pixels;
        for (const scale of [0.25, 2]) {
            camera.scale = scale; camera.update(); renderer.render([textObject], camera, props);
            const scaled = await renderer.readPixels();
            assert(scaled.array.every((v, i) => v === labelReference.array[i]), 'Text size and offsets must scale with molecular coordinates.');
        }
        camera.scale = 1; camera.update(); renderer.render([textObject], camera, props);
        results.push('SDF text labels');
        showPixels('SDF label', pixels);

        canvas.width = 320; canvas.height = 192;
        Object.assign(camera.viewport, { width: 320, height: 192 });
        camera.update();
        renderer.render([sphereObject], camera, props, true);
        pixels = await renderer.readPixels();
        assert(pixels.width === 320 && pixels.height === 192 && pixels.array[3] === 0, 'Resize and transparent image readback must work.');
        assert((await renderer.pick(160, 96))?.id.groupId === 4, 'Picking attachments must resize with the canvas.');
        results.push('resize, transparent background and RGBA readback');
        results.push(...await verifyPlugin());
        await context.device.queue.onSubmittedWorkDone();
        assert(errors.length === 0, `WebGPU validation errors: ${errors.join('; ')}`);
        const info = context.adapter.info;
        return { results, adapter: { vendor: info.vendor, architecture: info.architecture, device: info.device, description: info.description, isFallbackAdapter: info.isFallbackAdapter }, errors };
    } catch (error) {
        await context.device.queue.onSubmittedWorkDone();
        if (errors.length) throw new Error(`${String(error)}\nWebGPU validation errors: ${errors.join('; ')}`);
        throw error;
    } finally {
        renderer.dispose();
        context.dispose();
    }
}

async function verifySphereLodInstances(renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const snapshot = camera.getSnapshot();
    const params = { ...PD.getDefaultValues(Spheres.Params), ignoreLight: true, lodLevels: [
        { minDistance: 0, maxDistance: 25, overlap: 0, stride: 1, scaleBias: 1 },
        { minDistance: 25, maxDistance: 100, overlap: 0, stride: 2, scaleBias: 1 },
    ] };
    const sphere = Spheres.create(new Float32Array(3), new Float32Array([7]), 1);
    const matrices = new Float32Array(48);
    for (let i = 0; i < 3; i++) matrices.set(Mat4.fromTranslation(Mat4(), Vec3.create([-4, 0, 4][i], 0, [0, -20, -120][i])), i * 16);
    let reference: Uint8Array | undefined;
    try {
        for (const grid of [false, true]) {
            const values = Spheres.Utils.createValuesSimple(sphere, params, Color(0xff0000), 0.5, createTransform(matrices, 3, undefined, 0, 0));
            if (grid) updateTransformData(values, values.invariantBoundingSphere.ref.value, 4, 20);
            else ValueCell.update(values.aInstance, new Float32Array([2, 0, 1]));
            ValueCell.update(values.dColorType, 'instance');
            ValueCell.update(values.tColor, { array: new Uint8Array(grid ? [0, 0, 255, 255, 0, 0, 0, 255, 0] : [255, 0, 0, 0, 255, 0, 0, 0, 255]), width: 3, height: 1 });
            const object = createRenderObject('spheres', values, Spheres.Utils.createRenderableState(params), -1);
            // Keep the complete assembly inside the camera clip range while the
            // LOD distance logic independently culls the out-of-range instance.
            camera.setState({ position: Vec3.create(0, 0, 20), target: Vec3(), radius: 50, radiusMax: 50, fog: 0 }, 0); camera.update();
            renderer.render([object], camera, props, true);
            const pixels = await renderer.readPixels();
            assert(renderer.stats.instanceCount === 2 && renderer.stats.triangleCount === 28, 'Distance levels must submit only one near and one far instance, excluding the third instance.');
            for (const [center, id, channel] of [[Vec3.create(-4, 0, 0), grid ? 0 : 2, 2], [Vec3.create(0, 0, -20), grid ? 1 : 0, 0]] as const) {
                const projected = camera.project(Vec4(), center), x = Math.floor(projected[0]), y = Math.floor(256 - projected[1]);
                const pick = await renderer.pick(x, y);
                assert(pick?.id.instanceId === id && pick.id.groupId === 7, 'Culled sphere draws must preserve physical base instances and logical molecular picking IDs.');
                assert(pixels.array[(y * 256 + x) * 4 + channel] > 240, 'Culled sphere draws must preserve reordered instance colors.');
            }
            if (reference) assert(pixels.array.every((v, i) => v === reference![i]), 'Grid reordering must preserve the exact visible instance frame.');
            else reference = pixels.array;
            camera.setState({ position: Vec3.create(0, 0, 120) }, 0); camera.update();
            renderer.render([object], camera, props, true);
            assert(renderer.stats.drawCount === 0 && (await renderer.readPixels()).array.every((v, i) => i % 4 !== 3 || v === 0), 'Fully out-of-range assemblies must submit no sphere draws.');
            camera.setState({ position: Vec3.create(0, 0, 20) }, 0); camera.update(); renderer.render([object], camera, props, true);
            assert((await renderer.readPixels()).array.every((v, i) => v === pixels.array[i]), 'Restoring an assembly camera must restore the exact sphere instance frame.');
        }
    } finally { camera.setState(snapshot, 0); camera.update(); }
}

async function verifySphereLods(renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const snapshot = camera.getSnapshot(), scale = camera.scale;
    const centers = [-6, -2, 2, 6];
    const geometry = Spheres.create(new Float32Array(centers.flatMap(x => [x, 0, 0])), new Float32Array([0, 1, 2, 3]), 4);
    const levels = [
        { minDistance: 0, maxDistance: 25, overlap: 0, stride: 1, scaleBias: 1 },
        { minDistance: 25, maxDistance: 100, overlap: 0, stride: 2, scaleBias: 1 },
    ];
    const params = { ...PD.getDefaultValues(Spheres.Params), ignoreLight: true, lodLevels: levels };
    const values = Spheres.Utils.createValuesSimple(geometry, params, Color(0xff8800), 0.5);
    const object = createRenderObject('spheres', values, Spheres.Utils.createRenderableState(params), -1);
    const draw = async (distance: number, mode: 'perspective' | 'orthographic' = 'perspective') => {
        camera.setState({ position: Vec3.create(0, 0, distance), target: Vec3(), radius: 10, radiusMax: 10, fog: 0, mode }, 0); camera.update();
        const before = renderer.sphereLodStats.indices;
        renderer.render([object], camera, props, true);
        return { pixels: await renderer.readPixels(), indices: renderer.sphereLodStats.indices - before, triangles: renderer.stats.triangleCount };
    };
    const check = async (groups: number[], radius: number) => {
        for (let group = 0; group < centers.length; group++) {
            const projected = camera.project(Vec4(), Vec3.create(centers[group] * camera.scale, 0, 0));
            const x = Math.floor(projected[0]), y = Math.floor(256 - projected[1]), hit = await renderer.pick(x, y);
            if (!groups.includes(group)) { assert(!hit, 'Coarse sphere distance levels must omit spheres outside their prepared subset.'); continue; }
            assert(hit?.id.groupId === group, `Sphere distance levels must retain the original molecular group (${camera.state.mode}, distance ${camera.state.position[2]}, scale ${camera.scale}, group ${group}, hit ${hit?.id.groupId}).`);
            const point = Vec3.scale(Vec3(), camera.unproject(Vec3(), Vec3.create(x + 0.5, 256 - y - 0.5, hit!.depth)), 1 / camera.scale);
            assert(Math.abs(Vec3.distance(point, Vec3.create(centers[group], 0, 0)) - radius) < 0.0001, 'Sphere distance levels must apply their analytical radius scale.');
        }
    };
    try {
        camera.scale = 1;
        const near = await draw(17); await check([0, 1, 2, 3], 0.5);
        const far = await draw(40); await check([0, 2], 1);
        assert(far.indices > 0 && far.indices * 2 === near.indices, 'Far sphere levels must submit exactly half the near-level indices for a stride-two subset.');
        assert(far.triangles * 2 === near.triangles, 'Renderer geometry statistics must count the submitted sphere level instead of all retained levels.');
        const orthographic = await draw(40, 'orthographic'); await check([0, 2], 1);
        assert(orthographic.indices === far.indices, 'Orthographic sphere levels must retain the same prepared subset and draw work.');
        camera.scale = 2;
        await draw(60); await check([0, 2], 1);
        camera.scale = 1;
        const beyond = await draw(120);
        assert(beyond.indices === 0 && beyond.pixels.array.every((v, i) => i % 4 !== 3 || v === 0), 'Out-of-range sphere levels must submit no indices and produce no visible fragments.');
        assert(renderer.stats.drawCount === 0 && renderer.stats.triangleCount === 0, 'Culled sphere levels must report no rendered geometry.');
        const restored = await draw(17);
        assert(restored.pixels.array.every((v, i) => v === near.pixels.array[i]), 'Returning to a near distance level must restore the exact frame.');
        const overlap = { ...params, lodLevels: [
            { ...levels[0], maxDistance: 26, overlap: 4 },
            { ...levels[1], minDistance: 22, overlap: 4 },
        ] };
        Spheres.Utils.updateValues(values, overlap);
        await draw(24);
        for (const [group, radius] of [[0, 0.5], [1, 0.25], [2, 0.5], [3, 0.25]]) {
            const projected = camera.project(Vec4(), Vec3.create(centers[group], 0, 0));
            const x = Math.floor(projected[0]), y = Math.floor(256 - projected[1]), hit = await renderer.pick(x, y);
            assert(hit?.id.groupId === group, 'Overlapping sphere distance levels must preserve molecular picking.');
            const point = camera.unproject(Vec3(), Vec3.create(x + 0.5, 256 - y - 0.5, hit!.depth));
            assert(Math.abs(Vec3.distance(point, Vec3.create(centers[group], 0, 0)) - radius) < 0.0001, 'Both sphere distance overlaps must shrink radii by the analytical smoothstep factor.');
        }
        Spheres.Utils.updateValues(values, params);
        const restoredLevels = await draw(17);
        assert(restoredLevels.pixels.array.every((v, i) => v === near.pixels.array[i]), 'Changing and restoring sphere level settings must restore their original frame.');
    } finally { camera.scale = scale; camera.setState(snapshot, 0); camera.update(); }
}

async function verifyGeometryResources(renderer: WebGPURenderer, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    assert(object.type === 'mesh', 'Geometry resource verification requires a mesh.');
    const values = object.values as MeshValues;
    const originalMarker = values.uMarker.ref.value, originalLight = values.dIgnoreLight.ref.value, originalMetalness = values.uMetalness.ref.value;
    const before = await renderer.readPixels();
    const initial = { ...renderer.geometryResourceStats };
    const textureUploads = renderer.textureResourceStats.uploads;
    const bufferAllocations = renderer.bufferResourceStats.allocations;
    const bufferUpdates = renderer.bufferResourceStats.updates;
    ValueCell.update(values.uMarker, 1);
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.some((v, i) => v !== before.array[i]), 'Marker updates must visibly highlight retained geometry.');
    ValueCell.update(values.uMarker, originalMarker);
    ValueCell.update(values.dIgnoreLight, false);
    ValueCell.update(values.uMetalness, 1);
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.some((v, i) => v !== before.array[i]), 'Material updates must visibly shade retained geometry.');
    ValueCell.update(values.dIgnoreLight, originalLight);
    ValueCell.update(values.uMetalness, originalMetalness);
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Restoring marker and material values must restore the exact frame.');
    assert(Object.keys(initial).every(key => renderer.geometryResourceStats[key as keyof typeof initial] === initial[key as keyof typeof initial]), 'Theme and material updates must retain geometry construction, vertex/index uploads and transform buffers.');
    assert(renderer.bufferResourceStats.allocations === bufferAllocations && renderer.bufferResourceStats.updates > bufferUpdates, 'Marker and material changes must update existing GPU buffers without allocations.');

    const transforms = values.aTransform.ref.value;
    const themeBeforeTransforms = { ...renderer.themeResourceStats };
    const shifted = transforms.slice(); shifted[12] = 2;
    ValueCell.update(values.aTransform, shifted);
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.some((v, i) => v !== before.array[i]), 'Changing an instance transform must move retained geometry.');
    assert(renderer.geometryResourceStats.transformUploads === initial.transformUploads + 1 && renderer.geometryResourceStats.vertexUploads === initial.vertexUploads && renderer.geometryResourceStats.indexUploads === initial.indexUploads, 'Transform updates must upload only the transform buffer.');
    ValueCell.update(values.aTransform, transforms); renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Restoring the transform must reproduce the original frame.');
    assert(renderer.bufferResourceStats.allocations === bufferAllocations, 'Changing and restoring transforms must update existing buffers without allocations.');
    assert(renderer.themeResourceStats.builds === themeBeforeTransforms.builds && renderer.themeResourceStats.uploads === themeBeforeTransforms.uploads, 'Transform-only updates must skip CPU theme expansion and retain the GPU theme buffer.');

    const originalPositionBox = values.aPosition.ref, positions = originalPositionBox.value;
    const moved = positions.slice(); for (let i = 0; i < moved.length; i += 3) moved[i] += 2;
    const beforePositionUpdate = { ...renderer.geometryResourceStats };
    ValueCell.update(values.aPosition, moved); renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.some((v, i) => v !== before.array[i]), 'Position updates must rebuild visible geometry.');
    assert(renderer.geometryResourceStats.builds === beforePositionUpdate.builds + 1 && renderer.geometryResourceStats.vertexUploads === beforePositionUpdate.vertexUploads + 1 && renderer.geometryResourceStats.indexUploads === beforePositionUpdate.indexUploads + 1 && renderer.geometryResourceStats.transformUploads === beforePositionUpdate.transformUploads, 'Changed geometry must upload new vertices and indices while retaining transforms.');
    // A replacement cell can have the same revision as its predecessor.
    const replacement = ValueCell.create(positions);
    ValueCell.update(replacement, positions);
    assert(replacement.ref.version === values.aPosition.ref.version, 'Replacement attributes must exercise equal revisions with distinct cell identities.');
    ValueCell.set(values.aPosition, replacement.ref);
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Replacing an attribute ValueCell must invalidate geometry by cell identity.');
    const beforeFailure = { ...renderer.geometryResourceStats };
    const writesBeforeFailure = renderer.bufferResourceStats.updates;
    ValueCell.update(values.dColorType, 'volume');
    let rejected = false;
    try {
        renderer.render([object], camera, props);
    } catch (error) {
        rejected = error instanceof Error && error.message.includes('grid');
    } finally {
        ValueCell.update(values.dColorType, 'uniform');
    }
    assert(rejected, 'Invalid spatial color data must reject the resource rebuild.');
    assert(renderer.bufferResourceStats.updates === writesBeforeFailure, 'Rejected resource rebuilds must not commit pending writes into live GPU buffers.');
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Failed resource rebuilds must preserve retained vertex/index/transform buffers.');
    assert(Object.keys(beforeFailure).every(key => renderer.geometryResourceStats[key as keyof typeof beforeFailure] === beforeFailure[key as keyof typeof beforeFailure]), 'Failed theme updates must preserve uploaded geometry and transforms.');
    assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Retained geometry must preserve canonical group picking.');
    assert(renderer.textureResourceStats.uploads === textureUploads, 'Geometry, transform, marker and material updates must retain every unchanged image/palette/grid texture.');
    ValueCell.set(values.aPosition, originalPositionBox);
    renderer.render([object], camera, props);
}

async function verifyIncrementalBufferSizes(renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const mesh = Mesh.create(new Float32Array([-4, -4, 0, 4, -4, 0, 0, 4, 0]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array([7, 7, 7]), 3, 1);
    const meshProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
    const values = Mesh.Utils.createValuesSimple(mesh, meshProps, Color(0xff0000), 1);
    const object = createRenderObject('mesh', values, Mesh.Utils.createRenderableState(meshProps), -1);
    renderer.render([object], camera, props);
    const before = await renderer.readPixels();
    ValueCell.update(values.dClipObjectCount, 1); ValueCell.update(values.dClipping, true);
    ValueCell.update(values.uClipObjectType, [3]); ValueCell.update(values.uClipObjectInvert, [false]);
    ValueCell.update(values.uClipObjectPosition, [0, 0, 0]); ValueCell.update(values.uClipObjectRotation, [0, 0, 0, 1]);
    ValueCell.update(values.uClipObjectScale, [100, 100, 100]); ValueCell.update(values.uClipObjectTransform, [...Mat4.identity()]);
    const masks = new Uint8Array(8); masks[7] = 1;
    ValueCell.update(values.tClipping, { array: masks, width: 8, height: 1 });
    renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'A packed group mask must exempt its triangle from clipping.');
    const allocations = renderer.bufferResourceStats.allocations;
    ValueCell.update(values.tClipping, { array: new Uint8Array(1), width: 1, height: 1 });
    renderer.render([object], camera, props);
    assert(!(await renderer.pick(128, 128)), 'Shrinking a packed channel must clear stale bytes beyond its new length.');
    assert(renderer.bufferResourceStats.allocations === allocations, 'Shrinking packed channels must retain their GPU buffer capacity.');
    ValueCell.update(values.tClipping, { array: masks, width: 8, height: 1 }); renderer.render([object], camera, props);
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Growing a channel within its capacity must restore its original pixels.');
    assert(renderer.bufferResourceStats.allocations === allocations, 'Growth within capacity must use an incremental buffer update.');
    const larger = new Uint8Array(12); larger[7] = 1;
    ValueCell.update(values.tClipping, { array: larger, width: 12, height: 1 }); renderer.render([object], camera, props);
    assert(renderer.bufferResourceStats.allocations === allocations + 1, 'Growth beyond capacity must allocate only the affected buffer.');
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Buffer growth must preserve rendered group masks.');
    ValueCell.update(values.tClipping, { array: new Uint8Array(0), width: 1, height: 1 }); renderer.render([object], camera, props);
    assert(!(await renderer.pick(128, 128)), 'Empty packed channels must clear all previously uploaded masks.');
    assert(renderer.bufferResourceStats.allocations === allocations + 1, 'Clearing a channel must retain its enlarged capacity.');
}

async function verifyTextureResources(renderer: WebGPURenderer, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    const values = object.values as MeshValues;
    const saved = {
        colorType: values.dColorType.ref, colorGrid: values.tColorGrid.ref, dimensions: values.uColorGridDim.ref,
        transform: values.uColorGridTransform.ref, marker: values.uMarker.ref, palette: values.tPalette.ref,
        overpaint: values.dOverpaint.ref, overpaintType: values.dOverpaintType.ref,
    };
    const grid = new WebGPUTextureData();
    const before = await renderer.readPixels();
    const initialUploads = renderer.textureResourceStats.uploads;
    try {
        grid.load({ width: 1, height: 1, array: new Uint8Array([0, 255, 0, 255]) });
        ValueCell.update(values.tColorGrid, grid);
        ValueCell.update(values.dColorType, 'volume');
        ValueCell.update(values.uColorGridDim, Vec3.create(1, 1, 1));
        ValueCell.update(values.uColorGridTransform, Vec4.create(0, 0, 0, 1));
        renderer.render([object], camera, props);
        const green = await renderer.readPixels(), offset = (128 * green.width + 128) * 4;
        assert(green.array[offset + 1] > 240 && green.array[offset] < 10 && green.array[offset + 2] < 10, 'Spatial texture updates must visibly replace uniform colors.');
        assert(renderer.textureResourceStats.uploads === initialUploads + 1, 'Activating a spatial grid must upload only the changed texture.');
        grid.load({ width: 1, height: 1, array: new Uint8Array([0, 0, 255, 255]) });
        renderer.render([object], camera, props);
        const blue = await renderer.readPixels();
        assert(blue.array[offset + 2] > 240 && blue.array[offset + 1] < 10, 'Reloading a grid without updating its ValueCell must refresh GPU colors.');
        assert(renderer.textureResourceStats.uploads === initialUploads + 2, 'A grid source revision must upload exactly one texture.');
        ValueCell.update(values.uMarker, 1); renderer.render([object], camera, props);
        assert((await renderer.readPixels()).array.some((v, i) => v !== blue.array[i]), 'Highlighting must affect retained spatial grid colors.');
        ValueCell.set(values.uMarker, saved.marker); renderer.render([object], camera, props);
        assert((await renderer.readPixels()).array.every((v, i) => v === blue.array[i]), 'Removing a marker must restore the exact spatial frame.');
        assert(renderer.textureResourceStats.uploads === initialUploads + 2, 'Marker updates must retain the uploaded spatial grid.');

        ValueCell.update(values.tPalette, { ...saved.palette.value, array: saved.palette.value.array.slice() });
        ValueCell.update(values.dOverpaint, true); ValueCell.update(values.dOverpaintType, 'unsupported');
        let rejected = false;
        try {
            renderer.render([object], camera, props);
        } catch (error) {
            rejected = error instanceof Error && error.message.includes('granularity');
        }
        assert(rejected, 'Invalid overpaint granularity must reject a rebuild after texture allocation.');
        ValueCell.set(values.dOverpaint, saved.overpaint); ValueCell.set(values.dOverpaintType, saved.overpaintType);
        ValueCell.set(values.tPalette, saved.palette);
        renderer.render([object], camera, props);
        assert((await renderer.readPixels()).array.every((v, i) => v === blue.array[i]), 'Failed rebuilds must preserve borrowed image/palette/grid textures for the corrected frame.');
        assert(renderer.textureResourceStats.uploads === initialUploads + 3, 'Failed rebuilds must release the new palette while retaining the previous palette and grid.');
        assert((await renderer.pick(128, 128))?.id.groupId === 7, 'Texture reuse and failed updates must retain canonical picking.');
    } finally {
        ValueCell.set(values.dColorType, saved.colorType); ValueCell.set(values.tColorGrid, saved.colorGrid);
        ValueCell.set(values.uColorGridDim, saved.dimensions); ValueCell.set(values.uColorGridTransform, saved.transform);
        ValueCell.set(values.uMarker, saved.marker); ValueCell.set(values.tPalette, saved.palette);
        ValueCell.set(values.dOverpaint, saved.overpaint); ValueCell.set(values.dOverpaintType, saved.overpaintType);
        renderer.render([object], camera, props); grid.destroy();
    }
    assert((await renderer.readPixels()).array.every((v, i) => v === before.array[i]), 'Restoring uniform themes must recover the original frame.');
}

async function verifyIlluminationSamples(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    const post = { ...PD.getDefaultValues(PostprocessingParams), enabled: false };
    const illumination = { ...PD.getDefaultValues(IlluminationParams), enabled: true, maxIterations: 2, denoise: false, rendersPerFrame: [1, 1] as [number, number] };
    const sampling = { ...PD.getDefaultValues(MultiSampleParams), mode: 'on' as const };
    const original = { ...camera.viewOffset };
    const draw = (changed: boolean) => renderer.render([object], camera, props, true, 1, undefined, post, undefined, sampling, changed, false, illumination);
    let checkedEdges = 0;
    try {
        for (const exponent of [2, 3]) for (const level of [0, 1, 2, 3, 4, 5]) {
            illumination.maxIterations = exponent; sampling.sampleLevel = level;
            const max = Math.pow(2, exponent), offsets = JitterVectors[level];
            renderer.render([object], camera, props, true, 1, undefined, post);
            const baseline = await renderer.readPixels(), pick = await renderer.pick(128, 128);
            const expected = Float32Array.from(baseline.array, value => value / max);
            let lastReference = baseline.array;
            for (let iteration = 1; iteration < max; iteration++) {
                const offset = offsets[Math.min(offsets.length - 1, Math.floor(iteration * offsets.length / max))];
                camera.viewOffset.enabled = true;
                Camera.setViewOffset(camera.viewOffset, camera.viewport.width, camera.viewport.height, offset[0], offset[1], camera.viewport.width, camera.viewport.height); camera.update();
                renderer.render([object], camera, props, true, 1, undefined, post);
                const pixels = await renderer.readPixels(); lastReference = pixels.array;
                for (let i = 0; i < expected.length; i++) expected[i] += pixels.array[i] / max;
            }
            Object.assign(camera.viewOffset, original); camera.update(); draw(true);
            assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Supersampled illumination must start with a complete canonical frame.');
            for (let iteration = 1; iteration < max; iteration++) {
                draw(false);
                assert(JSON.stringify(camera.viewOffset) === JSON.stringify(original), 'Every illumination jitter iteration must restore the original camera offset.');
            }
            const pixels = await renderer.readPixels();
            assert(!renderer.illuminationNeedsFrame && Number(renderer.illuminationProgress) === max && pixels.array.every((v, i) => Math.abs(v - expected[i]) <= 2), `Illumination sampling level ${level} with ${max} iterations must match an independently averaged set of GPU frames.`);
            assert(JSON.stringify(await renderer.pick(128, 128)) === JSON.stringify(pick), 'Illumination supersampling must retain canonical picking IDs and depth.');
            let edges = 0;
            for (let i = 0; i < baseline.array.length && edges < 6; i += 4) if (baseline.array[i + 3] !== lastReference[i + 3]) {
                const pixel = i / 4, edge = await renderer.pick(pixel % baseline.width, Math.floor(pixel / baseline.width));
                assert(baseline.array[i + 3] ? JSON.stringify(edge?.id) === JSON.stringify(pick?.id) && Math.abs(edge!.depth - pick!.depth) < 0.000001 : !edge, 'Edge picking must follow the canonical frame even where the final jitter changes geometry coverage.');
                edges++;
            }
            checkedEdges += edges;

            if (level > 0) assert(pixels.array.some((v, i) => i % 4 === 3 && v > 0 && v < 255), 'Illumination jitter must produce fractional edge coverage.');
        }
        assert(checkedEdges > 0, 'Jitter verification must include changed geometry coverage at actual raster edges.');
        illumination.maxIterations = 2; sampling.sampleLevel = 2;
        draw(true); for (let i = 1; i < 4; i++) draw(false);
        const pixels = await renderer.readPixels(), straight = new Uint8ClampedArray(pixels.array);
        for (let i = 0; i < straight.length; i += 4) if (straight[i + 3] > 0 && straight[i + 3] < 255) for (let c = 0; c < 3; c++) straight[i + c] = straight[i + c] * 255 / straight[i + 3];
        const pass = new WebGPUImagePass(context, camera, () => [object], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: post, illumination, multiSample: sampling });
        const exported = await pass.getImageData(RuntimeContext.Synchronous, pixels.width, pixels.height);
        assert(exported.data.every((v, i) => v === straight[i]), 'Supersampled illumination exports must finish all iterations and convert premultiplied edge colors to straight alpha.');
        await pass.dispose();
        sampling.sampleLevel = 0;
        camera.viewOffset.enabled = true;
        Camera.setViewOffset(camera.viewOffset, camera.viewport.width * 2, camera.viewport.height * 2, camera.viewport.width / 2, camera.viewport.height / 2, camera.viewport.width, camera.viewport.height); camera.update();
        const crop = { ...camera.viewOffset };
        renderer.render([object], camera, props, true, 1, undefined, post); const cropped = await renderer.readPixels();
        draw(true); for (let i = 1; i < 4; i++) draw(false);
        assert(JSON.stringify(camera.viewOffset) === JSON.stringify(crop) && (await renderer.readPixels()).array.every((v, i) => Math.abs(v - cropped.array[i]) <= 1), 'Illumination sampling must preserve an existing cropped camera view.');
    } finally { Object.assign(camera.viewOffset, original); camera.update(); renderer.render([object], camera, props); }
}

async function verifyIlluminationPipeline(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    const post = { ...PD.getDefaultValues(PostprocessingParams), enabled: false }, samples = { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' as const };
    const illumination = { ...PD.getDefaultValues(IlluminationParams), enabled: true, maxIterations: 2, denoise: false, rendersPerFrame: [1, 1] as [number, number] };
    renderer.render([object], camera, props, true, 1, undefined, post);
    const baseline = await renderer.readPixels(), picked = await renderer.pick(128, 128);
    const draw = (changed: boolean) => renderer.render([object], camera, props, true, 1, undefined, post, undefined, samples, changed, false, illumination);
    const encoderPrototype = Object.getPrototypeOf(context.device.createCommandEncoder()), beginRenderPass = encoderPrototype.beginRenderPass;
    const countInputs = (action: () => void) => {
        let inputs = 0, thickness = 0;
        encoderPrototype.beginRenderPass = function (descriptor: GPURenderPassDescriptor) {
            if (descriptor.label === 'molstar-tracing-gbuffer') inputs++;
            if (descriptor.label === 'molstar-tracing-thickness') thickness++;
            return beginRenderPass.call(this, descriptor);
        };
        try { action(); } finally { encoderPrototype.beginRenderPass = beginRenderPass; }
        return { inputs, thickness };
    };
    const dispatches = context.stats.computeDispatches;
    const stable = countInputs(() => {
        draw(true);
        assert(Number(renderer.illuminationProgress) === 1 && renderer.illuminationNeedsFrame, 'A changed illuminated scene must start progressive tracing.');
        const quality = renderer.illuminationQuality;
        assert(quality.steps === 16 && quality.refineSteps === 1 && quality.rendersPerFrame === 1, 'The first native illumination iteration must use the existing inexpensive workload.');
        for (let i = 1; i < 4; i++) draw(false);
    });
    assert(stable.inputs === 1 && stable.thickness === 1 && context.stats.computeDispatches - dispatches === 4, 'Stable illumination must reuse opaque and thickness inputs while tracing every iteration.');
    const fixed = countInputs(() => renderer.render([object], camera, props, true, 1, undefined, post, undefined, samples, true, false, { ...illumination, thicknessMode: 'fixed', maxIterations: 0 }));
    assert(fixed.inputs === 1 && fixed.thickness === 0, 'Fixed thickness must skip the unused back-depth geometry pass.');
    countInputs(() => { draw(true); for (let i = 1; i < 4; i++) draw(false); });
    assert(Number(renderer.illuminationProgress) === 4 && !renderer.illuminationNeedsFrame, 'Native illumination must stop at its configured iteration count.');
    const lit = await renderer.readPixels();
    assert(lit.array.every((v, i) => Math.abs(v - baseline.array[i]) <= 1), 'Illumination must retain the analytical color and coverage of isolated diffuse geometry.');
    assert(JSON.stringify(await renderer.pick(128, 128)) === JSON.stringify(picked), 'Illumination must preserve canonical picking IDs and depth.');
    draw(false); assert((await renderer.readPixels()).array.every((v, i) => v === lit.array[i]), 'Converged illumination must retain its output when idle.');
    const pass = new WebGPUImagePass(context, camera, () => [object], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: post, illumination, multiSample: samples });
    const exported = await pass.getImageData(RuntimeContext.Synchronous, baseline.width, baseline.height);
    assert(exported.data.every((v, i) => Math.abs(v - lit.array[i]) <= 1), 'Illuminated screenshots must finish all tracing iterations before readback.');
    await pass.dispose();
    const overlayMesh = Mesh.create(new Float32Array([-4, -4, 1, 4, -4, 1, 0, 4, 1]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
    const overlayProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true, alpha: 0.5 };
    const overlay = createRenderObject('mesh', Mesh.Utils.createValuesSimple(overlayMesh, overlayProps, Color(0x0000ff), 1), Mesh.Utils.createRenderableState(overlayProps), -1);
    for (const objects of [[object, overlay], [overlay]]) {
        renderer.render(objects, camera, props, true, 1, undefined, post);
        const transparentReference = await renderer.readPixels();
        renderer.render(objects, camera, props, true, 1, undefined, post, undefined, samples, true, false, { ...illumination, maxIterations: 0 });
        assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - transparentReference.array[i]) <= 1), 'Illumination must preserve premultiplied transparent overlays both over opaque geometry and over an empty background.');
    }
    illumination.denoise = true; draw(false);
    assert(Number(renderer.illuminationProgress) === 1 && renderer.illuminationNeedsFrame, 'Illumination setting changes must restart accumulation.');
    if (object.type === 'mesh') {
        const originalColor = Vec3.clone(object.values.uColor.ref.value);
        ValueCell.update(object.values.uColor, Vec3.create(0, 1, 0)); const refreshed = countInputs(() => draw(true));
        assert(refreshed.inputs === 1 && refreshed.thickness === 1, 'Material changes must refresh both cached illumination input passes.');
        const updated = await renderer.readPixels(), center = (128 * updated.width + 128) * 4;
        assert(Number(renderer.illuminationProgress) === 1 && updated.array[center + 1] > 240 && updated.array[center] < 10, 'Material updates must restart illumination and appear in the newly traced frame.');
        ValueCell.update(object.values.uColor, originalColor); draw(true);
    }

    renderer.render([object], camera, props, true, 1, undefined, post, undefined, samples, true, false, { ...illumination, enabled: false });
    assert(!renderer.illuminationNeedsFrame && (await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Disabling illumination must restore normal rendering and stop progressive frames.');
}

async function verifyWorldHelpers(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    let pixelRatio = 1;
    const on = { name: 'on' as const, params: HandleHelperParams.handle.map('on').defaultValue };
    const handle = createWebGPUHandleHelper(() => pixelRatio, { handle: on });
    const pointer = createWebGPUPointerHelper({ enabled: 'on', color: Color(0x00ff00), hitColor: Color(0xff0000) });
    const post = PD.getDefaultValues(PostprocessingParams);
    post.occlusion = { name: 'off', params: {} }; post.antialiasing = { name: 'off', params: {} };
    const draw = (objects: GraphicsRenderObject[] = []) => renderer.render(objects, camera, props, true, 1, undefined, post);
    const project = (position: Vec3) => { const p = camera.project(Vec4(), position); return [Math.floor(p[0]), Math.floor(context.canvas.height - p[1])] as const; };
    const savedScale = camera.scale;
    try {
        renderer.setWorldHelpers(handle, undefined);
        handle.update(camera, Vec3.create(2, 1, 0), Mat3.identity()); draw();
        const endpoint = project(Vec3.create(5.3, 1, 0));
        const picked = await renderer.pick(...endpoint);
        assert(picked?.id.groupId === HandleGroup.TranslateObjectX && isHandleLoci(handle.getLoci(picked.id)), 'Native handle endpoints must retain helper loci after translation.');
        const baseline = await renderer.readPixels();
        handle.mark(handle.getLoci(picked.id), MarkerAction.Highlight); draw();
        assert((await renderer.readPixels()).array.some((v, i) => v !== baseline.array[i]), 'Handle highlighting must update native colors.');
        handle.mark(handle.getLoci(picked.id), MarkerAction.RemoveHighlight); draw();
        assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Removing handle highlights must restore the original output.');
        const rotation = Mat3.fromMat4(Mat3(), Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.unitZ));
        handle.update(camera, Vec3.create(2, 1, 0), rotation); draw();
        assert((await renderer.pick(...project(Vec3.create(2, 4.3, 0))))?.id.groupId === HandleGroup.TranslateObjectX, 'Handle rotation must preserve world translation.');
        const sphere = handle.getBoundingSphere(Sphere3D(), 0);
        pixelRatio = 2; draw();
        assert(Math.abs(handle.getBoundingSphere(Sphere3D(), 0).radius - sphere.radius * 2) < 0.001, 'Handle geometry must follow display scaling.');
        const transform = handle.getRenderObjects()[0].values.aTransform.ref.value;
        assert(transform[12] === 2 && transform[13] === 1, 'Display scaling must preserve the handle pose.');
        pixelRatio = 1; handle.getRenderObjects(); handle.update(camera, Vec3(), Mat3.identity());
        const planeMesh = Mesh.create(new Float32Array([-8, -8, 5, 8, -8, 5, 8, 8, 5, -8, 8, 5]), new Uint32Array([0, 1, 2, 0, 2, 3]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(4), 4, 2);
        const planeProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
        const plane = createRenderObject('mesh', Mesh.Utils.createValuesSimple(planeMesh, planeProps, Color(0x0000ff), 1), Mesh.Utils.createRenderableState(planeProps), -1);
        draw([plane]);
        assert((await renderer.pick(128, 128))?.id.objectId === plane.id, 'Opaque scene geometry must occlude handle picking.');
        handle.update(camera, Vec3.create(0, 0, 6), Mat3.identity()); draw([plane]);
        assert((await renderer.pick(128, 128))?.id.objectId === handle.getRenderObjects()[0].id, 'Foreground handles must win depth-tested picking.');
        renderer.setWorldHelpers(undefined, undefined); draw([plane]);
        const plain = await renderer.readPixels(), center = (128 * plain.width + 128) * 4;
        renderer.setWorldHelpers(undefined, pointer); pointer.update([], [], Vec3.create(0, 0, 6)); draw([plane]);
        const red = await renderer.readPixels();
        assert(red.array[center] > 240 && red.array[center + 2] < 10, 'Pointer hit markers must render in front of geometry.');
        assert((await renderer.pick(128, 128))?.id.objectId === plane.id, 'Pointer overlays must preserve molecular picking.');
        pointer.setProps({ hitColor: Color(0xff00ff) }); draw([plane]);
        const magenta = await renderer.readPixels();
        assert(magenta.array[center] > 240 && magenta.array[center + 2] > 240, 'Pointer color updates must refresh the group theme.');
        pointer.setProps({ alpha: 0.5 }); draw([plane]);
        const blended = await renderer.readPixels();
        assert(Math.abs(blended.array[center] - 128) < 3 && blended.array[center + 2] > 250, 'Transparent pointer markers must blend over the scene.');
        pointer.update([], [], undefined); draw([plane]);
        assert((await renderer.readPixels()).array.every((v, i) => v === plain.array[i]), 'Clearing pointer input must remove stale geometry.');
        pointer.setProps({ alpha: 1 }); pointer.update([], [], Vec3()); draw();
        const coverage = (pixels: Uint8Array) => pixels.reduce((sum, v, i) => sum + (i % 4 === 3 && v > 0 ? 1 : 0), 0);
        const area = coverage((await renderer.readPixels()).array);
        camera.scale = 3; camera.update(); draw();
        assert(Math.abs(coverage((await renderer.readPixels()).array) - area) <= 4, 'Pointer markers must preserve display size under model scaling.');
        camera.scale = savedScale; camera.update();
        pointer.update([], [Vec3.create(-2, 0, 0)], Vec3.create(2, 0, 0));
        handle.update(camera, Vec3.create(2, 1, 0), rotation); renderer.setWorldHelpers(handle, pointer); draw();
        const annotated = await renderer.readPixels();
        const image = new WebGPUImagePass(context, camera, () => [], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: post, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } }, undefined, { handle, pointer });
        try {
            const exported = await image.getImageData(RuntimeContext.Synchronous, annotated.width, annotated.height);
            assert(exported.data.every((v, i) => v === annotated.array[i]), 'Screenshot exports must retain handle and pointer geometry.');
        } finally { await image.dispose(); }
        assert(handle.getRenderObjects().length > 0 && pointer.getRenderObjects().length > 0, 'Disposing an export must preserve borrowed helper scenes.');
        const input = renderer.renderTracingInput([], camera, props);
        const normal = await context.readTexture(input.normal, 0, 0, annotated.width, annotated.height, 8);
        assert(normal.every(v => v === 0), 'World helper overlays must stay out of illumination inputs.');
    } finally {
        renderer.setWorldHelpers(undefined, undefined); handle.scene.clear(); pointer.scene.clear(); camera.scale = savedScale; camera.update(); draw();
    }
}

async function verifyStereo(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    const viewport = { ...camera.viewport }, snapshot = camera.getSnapshot();
    const stereoProps = { ...DefaultStereoCameraProps, eyeSeparation: 0.1, focus: 1 };
    const eyes = new WebGPUStereoCamera();
    const post = PD.getDefaultValues(PostprocessingParams);
    post.occlusion = { name: 'off', params: {} }; post.outline = { name: 'off', params: {} };
    post.antialiasing = { name: 'fxaa', params: PD.getDefaultValues(FxaaParams) };
    const output = { width: context.canvas.width, height: context.canvas.height, present: false };
    const readIds = async () => {
        const texture = (renderer as unknown as { outputPicking: GPUTexture }).outputPicking;
        return new Uint32Array((await context.readTexture(texture, 0, 0, output.width, output.height, 16)).buffer);
    };
    try {
        for (const v of [viewport, { x: 13, y: 17, width: 221, height: 193 }, { x: 0, y: 0, width: 257, height: 191 }]) {
            output.width = Math.max(context.canvas.width, v.x + v.width);
            output.height = v.width === 257 ? v.height : context.canvas.height;
            Object.assign(camera.viewport, v); camera.update(); eyes.update(camera, stereoProps);
            for (const mode of ['off', 'on', 'temporal'] as const) {
                const sampling = { ...PD.getDefaultValues(MultiSampleParams), mode, sampleLevel: 3, reuseOcclusion: false };
                const references = [];
                for (const eye of [eyes.left, eyes.right]) {
                    renderer.render([object], eye, props, true, 1, output, post, undefined, sampling, true);
                    for (let i = 0; renderer.multiSampleNeedsFrame && i < 16; i++) renderer.render([object], eye, props, true, 1, output, post, undefined, sampling, false);
                    assert(!renderer.multiSampleNeedsFrame, 'Independent reference eyes must finish their configured sampling.');
                    references.push({ pixels: (await renderer.readPixels()).array, ids: await readIds() });
                }
                renderer.renderStereo([object], camera, stereoProps, props, true, 1, output, post, undefined, sampling, true);
                for (let i = 0; renderer.multiSampleNeedsFrame && i < 16; i++) renderer.renderStereo([object], camera, stereoProps, props, true, 1, output, post, undefined, sampling, false);
                assert(!renderer.multiSampleNeedsFrame, 'Stereo temporal sampling must converge independently for both eyes.');
                const combined = await renderer.readPixels(), ids = await readIds();
                for (let y = 0; y < output.height; y++) for (let x = 0; x < output.width; x++) {
                    const right = x >= eyes.right.viewport.x && x < eyes.right.viewport.x + eyes.right.viewport.width && y >= output.height - v.y - v.height && y < output.height - v.y;
                    const reference = references[right ? 1 : 0], i = (y * output.width + x) * 4;
                    for (let c = 0; c < 4; c++) {
                        assert(Math.abs(combined.array[i + c] - reference.pixels[i + c]) <= 1, `Stereo ${mode} colors must match independently rendered eyes at (${x}, ${y}).`);
                        assert(ids[i + c] === reference.ids[i + c], 'Stereo selection must preserve each canonical eye ID and packed depth.');
                    }
                }
                for (const eye of [eyes.left, eyes.right]) {
                    const point = eye.project(Vec4(), Vec3()), x = Math.floor(point[0]), y = Math.floor(output.height - point[1]);
                    const picked = await renderer.pick(x, y);
                    assert(picked?.id.objectId === object.id, 'Both stereo eyes must retain molecular object identity.');
                    const active = renderer.getPickingCamera(x, y, camera);
                    assert(Mat4.areEqual(active.projection, eye.projection, 1e-6) && Mat4.areEqual(active.view, eye.view, 1e-6), 'Stereo picking must select the corresponding asymmetric eye camera.');
                    const world = active.unproject(Vec3(), Vec3.create(x, output.height - y, picked!.depth));
                    assert(Math.abs(world[2]) < 1e-4, 'Stereo picking must reconstruct the known mesh plane instead of using the parent projection.');
                }
                const retained = combined.array;
                if (mode === 'temporal') {
                    renderer.renderStereo([object], camera, stereoProps, props, true, 1, output, post, undefined, sampling, false);
                    assert((await renderer.readPixels()).array.every((value, i) => value === retained[i]), 'Converged stereo frames must remain stable when idle.');
                    renderer.renderStereo([object], camera, { ...stereoProps, focus: 10 }, props, true, 1, output, post, undefined, sampling, false);
                    assert(renderer.multiSampleNeedsFrame, 'Changed stereo settings must restart both temporal eye accumulators.');
                }
            }
        }
        output.width = context.canvas.width; output.height = context.canvas.height;
        Object.assign(camera.viewport, viewport); camera.update();
        renderer.render([object], camera, props, true, 1, output, post);
        const mono = await renderer.readPixels();
        renderer.renderStereo([object], camera, stereoProps, props, true, 1, output, post);
        assert((await renderer.readPixels()).array.some((value, i) => value !== mono.array[i]), 'Enabling stereo must replace the single view with two projected views.');
        renderer.render([object], camera, props, true, 1, output, post);
        assert((await renderer.readPixels()).array.every((value, i) => value === mono.array[i]) && renderer.getPickingCamera(128, 128, camera) === camera, 'Disabling stereo must restore mono pixels and picking matrices.');
    } finally { Object.assign(camera.viewport, viewport); camera.setState(snapshot, 0); camera.update(); renderer.render([object], camera, props); }
}

async function verifyCameraAxes(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const helperProps = PD.getDefaultValues(CameraHelperParams);
    assert(helperProps.axes.name === 'on', 'Orientation axes must be enabled by default.');
    const axes = helperProps.axes.params;
    axes.alpha = 1; axes.showPlanes = false;
    let pixelRatio = 1;
    const helper = createWebGPUCameraHelper(() => pixelRatio, helperProps);
    const snapshot = camera.getSnapshot(), scale = camera.scale;
    const post = { ...PD.getDefaultValues(PostprocessingParams), enabled: false };
    const draw = () => renderer.render([], camera, props, true, pixelRatio, undefined, post);
    const ids = async () => {
        const texture = (renderer as unknown as { selectionTexture: GPUTexture }).selectionTexture;
        const packed = new Uint32Array((await context.readTexture(texture, 0, 0, texture.width, texture.height, 16)).buffer);
        const points = new Map<number, { x: number, y: number }>();
        for (let i = 0; i < packed.length; i += 4) if (packed[i]) points.set(packed[i + 2], { x: i / 4 % texture.width, y: Math.floor(i / 4 / texture.width) });
        return points;
    };
    renderer.setCameraHelper(helper);
    try {
        for (const location of ['bottom-left', 'bottom-right', 'top-left', 'top-right'] as const) {
            helper.setProps({ axes: { name: 'on', params: { ...axes, location } } }); draw();
            const pixels = await renderer.readPixels(), points = await ids();
            const x = points.get(CameraHelperAxis.X), y = points.get(CameraHelperAxis.Y);
            assert(x && y, `${location} axes must provide both X and Y picking groups.`);
            for (const point of [x, y]) {
                assert(location.endsWith('left') ? point.x < pixels.width / 2 : point.x > pixels.width / 2, 'Axis placement must follow its horizontal corner setting.');
                assert(location.startsWith('top') ? point.y < pixels.height / 2 : point.y > pixels.height / 2, 'Axis placement must follow its vertical corner setting.');
            }
            const picked = await renderer.pick(x.x, x.y);
            assert(picked && isCameraAxesLoci(helper.getLoci(picked.id)), 'Native helper picks must resolve existing camera-axis loci.');
            assert(!(await renderer.pick(128, 128)), 'Orientation axes must leave the empty scene center unpickable.');
            assert(pixels.array.some((v, i) => i % 4 === 0 && v > 200) && pixels.array.some((v, i) => i % 4 === 1 && v > 100), 'Configured red and green axis colors must appear in the native output.');
        }
        helper.setProps({ axes: { name: 'on', params: { ...axes, location: 'bottom-left' } } }); draw();
        const baseline = await renderer.readPixels(), points = await ids(), x = points.get(CameraHelperAxis.X)!;
        const picked = (await renderer.pick(x.x, x.y))!;
        helper.mark(helper.getLoci(picked.id), MarkerAction.Highlight); draw();
        assert((await renderer.readPixels()).array.some((v, i) => v !== baseline.array[i]), 'Hover highlighting must change the selected native axis.');
        assert(JSON.stringify(await renderer.pick(x.x, x.y)) === JSON.stringify(picked), 'Axis highlighting must preserve IDs and depth.');
        helper.mark(helper.getLoci(picked.id), MarkerAction.RemoveHighlight); draw();
        assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Removing an axis highlight must restore its original frame.');
        const antialiased = PD.getDefaultValues(PostprocessingParams);
        antialiased.occlusion = { name: 'off', params: {} }; antialiased.outline = { name: 'off', params: {} };
        antialiased.antialiasing = { name: 'fxaa', params: PD.getDefaultValues(FxaaParams) };
        renderer.render([], camera, props, true, 1, undefined, antialiased);
        assert((await renderer.readPixels()).array.some((v, i) => v !== baseline.array[i]), 'Native antialiasing must filter camera-helper edges.');
        assert(JSON.stringify(await renderer.pick(x.x, x.y)) === JSON.stringify(picked), 'Antialiased helper edges must preserve canonical picking.');
        helper.setProps({ axes: { name: 'on', params: { ...axes, location: 'bottom-left', alpha: props.pickingAlphaThreshold / 2 } } }); draw();
        const faint = await renderer.readPixels();
        assert(faint.array.some((v, i) => i % 4 === 3 && v > 0 && v < 255) && !(await renderer.pick(x.x, x.y)), 'Faint axes must render with transparent alpha while respecting the independent picking threshold.');
        helper.setProps({ axes: { name: 'on', params: { ...axes, location: 'bottom-left' } } }); draw();
        camera.scale = 3; camera.update(); draw();
        assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Orientation axes must remain invariant to molecular model scaling.');
        camera.scale = scale; camera.setState({ position: Vec3.create(15, 10, 20) }, 0); camera.update(); draw();
        assert((await renderer.readPixels()).array.some((v, i) => v !== baseline.array[i]) && (await ids()).has(CameraHelperAxis.Z), 'Camera rotation must rotate native axes and expose their Z group.');
        camera.setState(snapshot, 0); camera.update();
        pixelRatio = 2; draw();
        const scaled = await renderer.readPixels();
        const coverage = (array: Uint8Array) => array.filter((v, i) => i % 4 === 3 && v > 0).length;
        assert(coverage(scaled.array) > 2 * coverage(baseline.array), 'Axes must rebuild geometry when display pixel ratio changes.');
        pixelRatio = 1;
        helper.setProps({ axes: { name: 'on', params: { ...axes, location: 'bottom-left', showLabels: true, labelX: 'a', labelY: 'b', labelZ: 'c' } } }); draw();
        const labeled = await renderer.readPixels();
        assert(labeled.array.some((v, i) => v !== baseline.array[i]) && helper.getRenderObjects().length === 2, 'Configured labels must render with the native text path.');
        const pass = new WebGPUImagePass(context, camera, () => [], { renderer: props, cameraHelper: helper.props, transparentBackground: true, postprocessing: post, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } });
        try {
            const exported = await pass.getImageData(RuntimeContext.Synchronous, labeled.width, labeled.height), straight = new Uint8ClampedArray(labeled.array);
            for (let i = 0; i < straight.length; i += 4) if (straight[i + 3] > 0 && straight[i + 3] < 255) for (let c = 0; c < 3; c++) straight[i + c] = straight[i + c] * 255 / straight[i + 3];
            assert(exported.data.every((v, i) => v === straight[i]), 'Transparent exports must preserve camera axes and custom labels.');
            pass.setProps({ cameraHelper: { axes: { name: 'off', params: {} } } });
            assert((await pass.getImageData(RuntimeContext.Synchronous, labeled.width, labeled.height)).data.every(v => v === 0), 'Disabling screenshot axes must remove their geometry and labels.');
        } finally { await pass.dispose(); }
        renderer.renderTracingInput([], camera, props);
        const normal = (renderer as unknown as { tracingInput: { textures: { normal: GPUTexture } } }).tracingInput.textures.normal;
        assert((await context.readTexture(normal, 0, 0, normal.width, normal.height, 8)).every(v => v === 0), 'Orientation axes must never enter the molecular tracing G-buffer.');
        helper.setProps({ axes: { name: 'off', params: {} } }); draw();
        assert((await renderer.readPixels()).array.every(v => v === 0) && !(await renderer.pick(x.x, x.y)), 'Disabling axes must remove their color and picking IDs.');
    } finally { renderer.setCameraHelper(undefined); helper.scene.clear(); camera.scale = scale; camera.setState(snapshot, 0); camera.update(); }
}

async function verifyIlluminationEffects(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    const illumination = { ...PD.getDefaultValues(IlluminationParams), enabled: true, maxIterations: 0, denoise: false, ignoreOutline: false, rendersPerFrame: [1, 1] as [number, number] };
    const samples = { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' as const };
    const base = PD.getDefaultValues(PostprocessingParams);
    base.occlusion = { name: 'off', params: {} }; base.shadow = { name: 'off', params: {} };
    base.outline = { name: 'off', params: {} }; base.antialiasing = { name: 'off', params: {} };
    const outline = { ...base, outline: { name: 'on' as const, params: { ...PD.getDefaultValues(OutlineParams), color: Color(0x0000ff) } } };
    const bloom = { ...base, bloom: { name: 'on' as const, params: { ...PD.getDefaultValues(BloomParams), mode: 'luminosity' as const, threshold: 0, strength: 0.3 } } };
    const dof = { ...base, dof: { name: 'on' as const, params: { ...PD.getDefaultValues(DofParams), inFocus: 4, PPM: 1, blurSize: 3 } } };
    const sharpen = { ...base, sharpening: { name: 'on' as const, params: PD.getDefaultValues(CasParams) } };
    const fxaa = { ...base, antialiasing: { name: 'fxaa' as const, params: PD.getDefaultValues(FxaaParams) } };
    const smaa = { ...base, antialiasing: { name: 'smaa' as const, params: PD.getDefaultValues(SmaaParams) } };
    const combined = { ...outline, bloom: bloom.bloom, dof: dof.dof, sharpening: sharpen.sharpening, antialiasing: smaa.antialiasing };
    renderer.render([object], camera, props, true, 1, undefined, base);
    const canonical = await renderer.pick(128, 128), unprocessed = await renderer.readPixels();
    for (const [name, post] of [['outline', outline], ['bloom', bloom], ['DOF', dof], ['sharpening', sharpen], ['FXAA', fxaa], ['SMAA', smaa], ['combined', combined]] as const) {
        renderer.render([object], camera, props, true, 1, undefined, post);
        const reference = await renderer.readPixels();
        // DOF and CAS preserve a uniformly colored surface over transparent pixels;
        // their combined case below also filters the gradients produced by bloom.
        if (name !== 'DOF' && name !== 'sharpening') assert(reference.array.some((v, i) => Math.abs(v - unprocessed.array[i]) > 2), `${name} must visibly affect the reference geometry.`);
        renderer.render([object], camera, props, true, 1, undefined, post, undefined, samples, true, false, illumination);
        const lit = await renderer.readPixels();
        assert(lit.array.every((v, i) => Math.abs(v - reference.array[i]) <= 2), `${name} must compose with isolated illuminated geometry like the ordinary effect path.`);
        assert(JSON.stringify(await renderer.pick(128, 128)) === JSON.stringify(canonical), `${name} with illumination must preserve canonical picking.`);
        if (name === 'combined') {
            const pass = new WebGPUImagePass(context, camera, () => [object], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: post, illumination, multiSample: samples });
            try {
                const exported = await pass.getImageData(RuntimeContext.Synchronous, lit.width, lit.height);
                const straight = new Uint8ClampedArray(lit.array);
                for (let i = 0; i < straight.length; i += 4) if (straight[i + 3] > 0 && straight[i + 3] < 255) for (let c = 0; c < 3; c++) straight[i + c] = straight[i + c] * 255 / straight[i + 3];
                assert(exported.data.every((v, i) => Math.abs(v - straight[i]) <= 2), 'Combined illumination effects must survive transparent screenshot exports.');
            } finally { await pass.dispose(); }
        }
    }
    renderer.render([object], camera, props, true, 1, undefined, outline, undefined, samples, true, false, { ...illumination, ignoreOutline: true });
    assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - unprocessed.array[i]) <= 1), 'ignoreOutline must suppress configured outlines during illumination.');
    assert(outline.outline.name === 'on', 'Illumination must preserve caller-owned outline settings.');

    const shaded = { ...base, occlusion: { name: 'on' as const, params: PD.getDefaultValues(SsaoParams) }, shadow: { name: 'on' as const, params: PD.getDefaultValues(ShadowParams) } };
    const prototype = Object.getPrototypeOf(context.device.createCommandEncoder()), begin = prototype.beginRenderPass;
    const passes: string[] = [];
    prototype.beginRenderPass = function (descriptor: GPURenderPassDescriptor) { passes.push(descriptor.label ?? ''); return begin.call(this, descriptor); };
    try { renderer.render([object], camera, props, true, 1, undefined, shaded, undefined, samples, true, false, illumination); } finally { prototype.beginRenderPass = begin; }
    assert(!passes.includes('molstar-postprocess-ssao') && !passes.includes('molstar-postprocess-shadowCompose'), 'Opaque illumination must skip duplicate SSAO and screen-space shadow passes.');
    assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - unprocessed.array[i]) <= 1), 'Opaque illumination must retain its traced lighting when SSAO and screen-space shadows are configured.');
    assert(shaded.shadow.name === 'on' && shaded.occlusion.name === 'on', 'Illumination must preserve caller-owned shadow and SSAO settings.');

    const overlayMesh = Mesh.create(new Float32Array([-4, -4, 1, 4, -4, 1, 0, 4, 1]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
    const overlayProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true, alpha: 0.5 };
    const overlay = createRenderObject('mesh', Mesh.Utils.createValuesSimple(overlayMesh, overlayProps, Color(0x0000ff), 1), Mesh.Utils.createRenderableState(overlayProps), -1);
    shaded.occlusion.params.transparentThreshold = 1;
    passes.length = 0;
    prototype.beginRenderPass = function (descriptor: GPURenderPassDescriptor) { passes.push(descriptor.label ?? ''); return begin.call(this, descriptor); };
    try { renderer.render([object, overlay], camera, props, true, 1, undefined, shaded, undefined, samples, true, false, illumination); } finally { prototype.beginRenderPass = begin; }
    assert(passes.includes('molstar-postprocess-ssao') && passes.includes('molstar-postprocess-ssaoTransparentColor') && !passes.includes('molstar-postprocess-shadowCompose'), 'Illumination must retain transparent SSAO without reapplying opaque screen-space shadows.');
    const mixed = await renderer.readPixels();
    assert(mixed.array[(128 * mixed.width + 128) * 4 + 2] > 0 && (await renderer.pick(128, 128))?.id.objectId === overlay.id, 'Transparent SSAO must retain illuminated overlay color and independent picking.');
    assert(Math.abs(mixed.array[(128 * mixed.width + 128) * 4] - 128) <= 2, 'Transparent SSAO must leave the traced red opaque layer intact under the half-opacity blue overlay.');
    renderer.render([object], camera, props);
}

async function verifyTracingBounces(context: WebGPUContext, renderer: WebGPURenderer, rendererProps: RendererProps) {
    const camera = new Camera({ position: Vec3.create(0, 0, 20), target: Vec3(), radius: 10, radiusMax: 10, fog: 0 }, { x: 0, y: 0, width: 64, height: 64 }); camera.update();
    const props = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
    const make = (vertices: number[], indices: number[], normals: number[], color: Color) => {
        const mesh = Mesh.create(new Float32Array(vertices), new Uint32Array(indices), new Float32Array(normals), new Float32Array(vertices.length / 3), vertices.length / 3, indices.length / 3);
        return createRenderObject('mesh', Mesh.Utils.createValuesSimple(mesh, props, color, 1), Mesh.Utils.createRenderableState(props), -1);
    };
    const floor = make([-6, -6, 0, 6, -6, 0, 6, 6, 0, -6, 6, 0], [0, 1, 2, 0, 2, 3], [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1], Color(0x808080));
    const wall = make([1, -6, 0, 1, 6, 0, 1, 6, 8, 1, -6, 8], [0, 2, 1, 0, 3, 2], [-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0], Color(0x00ff00));
    const input = renderer.renderTracingInput([floor, wall], camera, rendererProps, 64, 64);
    const tracer = await WebGPUTracing.create(context), tracing = PD.getDefaultValues(TracingParams);
    tracing.rendersPerFrame = [16, 16]; tracing.steps = 64; tracing.refineSteps = 4; tracing.bounces = 4;
    const lights = { ...rendererProps, light: [], ambientColor: Color(0xffffff), ambientIntensity: 0.5, exposure: 1 };
    try {
        const read = async (iteration: number) => new Float32Array((await context.readTexture(tracer.render(camera, lights, tracing, input, iteration), 0, 0, 64, 64, 16)).buffer);
        const first = await read(0), second = await read(1);
        assert(first.every(Number.isFinite) && second.every(Number.isFinite), 'Multiple diffuse bounces and Russian roulette must produce finite GPU values.');
        let red = 0, green = 0, samples = 0;
        for (let y = 22; y < 42; y++) for (let x = 28; x < 34; x++) {
            const i = (y * 64 + x) * 4; red += second[i]; green += second[i + 1]; samples++;
        }
        assert(green / samples > red / samples + 0.005 && red > 0, `A green wall must transfer green irradiance onto the neighboring neutral diffuse surface: red=${red / samples}, green=${green / samples}, center=${Array.from(second.slice((32 * 64 + 30) * 4, (32 * 64 + 30) * 4 + 4))}.`);
        assert(first.some((v, i) => i % 4 < 3 && Math.abs(v - second[i]) > 0.001), 'Progressive frames must use distinct PCG ray samples.');
        assert(second.every((v, i) => i % 4 !== 3 || v === 0 || v === 0.5), 'A second tracing frame must retain the analytical accumulation weight.');
        assert((await read(0)).every((v, i) => v === first[i]), 'Restarting the same scene must replay deterministic native PCG samples.');
        ValueCell.update(wall.values.uColor, Vec3.create(1, 0, 0));
        renderer.renderTracingInput([floor, wall], camera, rendererProps, 64, 64);
        await read(0); const recolored = await read(1);
        let redTransfer = 0, greenTransfer = 0;
        for (let y = 22; y < 42; y++) for (let x = 28; x < 34; x++) {
            const i = (y * 64 + x) * 4; redTransfer += recolored[i]; greenTransfer += recolored[i + 1];
        }
        assert(redTransfer / samples > greenTransfer / samples + 0.005, 'Changing the wall to red must reverse the neighboring surface color transfer.');

        tracing.thicknessMode = 'fixed'; tracing.thickness = 2;
        const fixed = await read(0);
        assert(fixed.every(Number.isFinite) && fixed.some((v, i) => i % 4 < 3 && v > 0), 'Fixed-thickness tracing must support diffuse bounce transport.');
        tracing.shadowEnable = true; tracing.shadowThickness = 0; tracing.shadowSoftness = 0.1;
        const shadowTexture = tracer.render(camera, rendererProps, tracing, input, 0);
        const shadowed = new Float32Array((await context.readTexture(shadowTexture, 0, 0, 64, 64, 16)).buffer);
        assert(shadowed.every(Number.isFinite), 'Automatic shadow thickness and soft-shadow sampling must retain finite ray values.');
    } finally { tracer.dispose(); }
}

async function verifyWeightedTransparency(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const vertices = new Float32Array([-4, -4, -2, 4, -4, 2, 0, 4, 0, -4, -4, 2, 4, -4, -2, 0, 4, 0]);
    const mesh = Mesh.create(vertices, new Uint32Array([0, 1, 2, 3, 4, 5]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array([0, 0, 0, 1, 1, 1]), 6, 2);
    const shapeProps = { ...PD.getDefaultValues(Mesh.Params), alpha: 0.5, ignoreLight: true, doubleSided: true, transparentBackfaces: 'on' as const };
    const values = Mesh.Utils.createValuesSimple(mesh, shapeProps, Color(0xffffff), 1);
    ValueCell.update(values.dColorType, 'group'); ValueCell.update(values.tColor, { array: new Uint8Array([255, 0, 0, 0, 0, 255]), width: 2, height: 1 });
    const shape = createRenderObject('mesh', values, Mesh.Utils.createRenderableState(shapeProps), -1);
    const offset = (128 * 256 + 128) * 4;
    const draw = async () => { renderer.render([shape], camera, props, true); return renderer.readPixels(); };
    try {
        renderer.setTransparency('wboit');
        const original = await draw();
        assert(Math.abs(original.array[offset] - 96) <= 2 && Math.abs(original.array[offset + 2] - 96) <= 2 && Math.abs(original.array[offset + 3] - 191) <= 1, 'Equal-opacity intersecting surfaces must match weighted color and multiplicative revealage.');
        ValueCell.update(values.elements, new Uint32Array([3, 4, 5, 0, 1, 2]));
        const reversed = await draw();
        assert(reversed.array.every((v, i) => v === original.array[i]), 'Weighted transparency must be independent of triangle order within the same object.');
        const project = (x: number) => { const p = camera.project(Vec4(), Vec3.create(x, 0, 0)); return [Math.floor(p[0]), Math.floor(context.canvas.height - p[1])] as const; };
        for (const x of [-1, 1]) {
            const [sx, sy] = project(x), i = (sy * 256 + sx) * 4, hit = await renderer.pick(sx, sy);
            assert(x > 0 ? original.array[i] > original.array[i + 2] : original.array[i + 2] > original.array[i], 'Depth weighting must favor the locally nearer intersecting surface.');
            assert(hit?.id.groupId === (x > 0 ? 0 : 1), 'Weighted color accumulation must preserve nearest molecular selection.');
        }
        ValueCell.update(values.dTransparency, true); ValueCell.update(values.tTransparency, { array: new Uint8Array([0, 128]), width: 2, height: 1 });
        const uneven = await draw(), blueAlpha = 0.5 * (1 - 128 / 255), expectedAlpha = 1 - 0.5 * (1 - blueAlpha);
        const redFraction = 0.25 / (0.25 + blueAlpha * blueAlpha);
        assert(Math.abs(uneven.array[offset] - expectedAlpha * redFraction * 255) <= 2 && Math.abs(uneven.array[offset + 2] - expectedAlpha * (1 - redFraction) * 255) <= 2 && Math.abs(uneven.array[offset + 3] - expectedAlpha * 255) <= 1, 'Unequal opacity must retain alpha-weighted color and multiplicative coverage.');
        const [faintX, faintY] = project(-1);
        renderer.render([shape], camera, { ...props, pickingAlphaThreshold: 0.4 }, true);
        assert((await renderer.pick(faintX, faintY))?.id.groupId === 0, 'Weighted faint foreground layers must allow selection of a stronger layer behind them.');
        renderer.render([shape], camera, { ...props, pickingAlphaThreshold: 0.1 }, true);
        assert((await renderer.pick(faintX, faintY))?.id.groupId === 1 && (await renderer.readPixels()).array.every((v, i) => v === uneven.array[i]), 'Changing selection thresholds must preserve weighted colors.');
        ValueCell.update(values.dTransparency, false);
        const split = [0, 1].map(group => {
            const part = Mesh.create(vertices.slice(group * 9, group * 9 + 9), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3).fill(group), 3, 1);
            return createRenderObject('mesh', Mesh.Utils.createValuesSimple(part, shapeProps, Color(group ? 0x0000ff : 0xff0000), 1), Mesh.Utils.createRenderableState(shapeProps), -1);
        });
        renderer.render(split, camera, props, true); const splitFrame = await renderer.readPixels();
        renderer.render([...split].reverse(), camera, props, true);
        assert((await renderer.readPixels()).array.every((v, i) => v === splitFrame.array[i]) && splitFrame.array.every((v, i) => v === original.array[i]), 'Weighted colors must agree across render-object grouping and object order.');
        const opaqueProps = { ...shapeProps, alpha: 1 };
        const opaqueMesh = Mesh.create(new Float32Array([-4, -4, 3, 4, -4, 3, 0, 4, 3]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
        const opaque = createRenderObject('mesh', Mesh.Utils.createValuesSimple(opaqueMesh, opaqueProps, Color(0x00ff00), 1), Mesh.Utils.createRenderableState(opaqueProps), -1);
        renderer.render([shape, opaque], camera, props, true); const occluded = await renderer.readPixels();
        assert(occluded.array[offset] === 0 && occluded.array[offset + 1] === 255 && occluded.array[offset + 2] === 0 && occluded.array[offset + 3] === 255, 'Weighted fragments behind opaque geometry must not contribute color or alpha.');
        renderer.setTransparency('blended'); const blended = await draw();
        assert(blended.array.some((v, i) => v !== original.array[i]), 'Explicit blended mode must retain ordinary ordered alpha composition.');
        renderer.setTransparency('wboit'); assert((await draw()).array.every((v, i) => v === original.array[i]), 'Changing transparency modes must invalidate and rebuild the rendered frame.');
        const pass = new WebGPUImagePass(context, camera, () => [shape], { renderer: props, transparentBackground: true, cameraHelper: { axes: { name: 'off', params: {} } }, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false }, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } }, undefined, { transparency: () => renderer.transparencyMode });
        try {
            for (const mode of ['wboit', 'blended', 'wboit'] as const) {
                renderer.setTransparency(mode);
                const exported = await pass.getImageData(RuntimeContext.Synchronous, 256, 256);
                const straight = new Uint8ClampedArray(mode === 'wboit' ? original.array : blended.array);
                for (let i = 0; i < straight.length; i += 4) if (straight[i + 3] > 0 && straight[i + 3] < 255) for (let c = 0; c < 3; c++) straight[i + c] = straight[i + c] * 255 / straight[i + 3];
                assert(exported.data.every((v, i) => v === straight[i]), 'Existing screenshot passes must follow live transparency modes and preserve straight alpha.');
            }
        } finally { await pass.dispose(); }
        const sampling = { ...PD.getDefaultValues(MultiSampleParams), mode: 'on' as const, sampleLevel: 2 };
        renderer.render([shape], camera, props, true, 1, undefined, undefined, undefined, sampling);
        const sampled = await renderer.readPixels();
        assert(Math.abs(sampled.array[offset] - original.array[offset]) <= 1 && Math.abs(sampled.array[offset + 3] - original.array[offset + 3]) <= 1, 'Supersampling must retain weighted color and alpha in constant-coverage regions.');
    } finally { renderer.setTransparency('blended'); }
}

async function verifyDepthPeeling(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const palette = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 1]];
    const positions = new Float32Array(54), normals = new Float32Array(54), groups = new Float32Array(18), indices = new Uint32Array(18);
    for (let layer = 0; layer < 6; layer++) {
        positions.set([-4, -4, 6 - layer * 2, 4, -4, 6 - layer * 2, 0, 4, 6 - layer * 2], layer * 9);
        normals.set([0, 0, 1, 0, 0, 1, 0, 0, 1], layer * 9); groups.fill(layer, layer * 3, layer * 3 + 3); indices.set([layer * 3, layer * 3 + 1, layer * 3 + 2], layer * 3);
    }
    const shapeProps = { ...PD.getDefaultValues(Mesh.Params), alpha: 0.5, ignoreLight: true, doubleSided: true, transparentBackfaces: 'on' as const };
    const geometry = Mesh.create(positions, indices, normals, groups, 18, 6), values = Mesh.Utils.createValuesSimple(geometry, shapeProps, Color(0xffffff), 1);
    ValueCell.update(values.dColorType, 'group'); ValueCell.update(values.tColor, { array: Uint8Array.from(palette.flat(), v => v * 255), width: 6, height: 1 });
    const shape = createRenderObject('mesh', values, Mesh.Utils.createRenderableState(shapeProps), -1), offset = (128 * 256 + 128) * 4;
    const reverse = new Uint32Array(18); for (let i = 0; i < 6; i++) reverse.set(indices.subarray((5 - i) * 3, (6 - i) * 3), i * 3);
    const draw = async (iterations: number) => { renderer.setDpoitIterations(iterations); renderer.render([shape], camera, props, true); return renderer.readPixels(); };
    let complete: Awaited<ReturnType<WebGPURenderer['readPixels']>> | undefined;
    try {
        renderer.setTransparency('dpoit');
        for (const iterations of [1, 2, 3, 10]) {
            const selected = [0, 1, 2, 3, 4, 5].filter(i => i < iterations || i >= 6 - iterations);
            const expected = [0, 0, 0, 0];
            for (const i of selected) {
                const weight = 0.5 * (1 - expected[3]);
                for (let c = 0; c < 3; c++) expected[c] += palette[i][c] * weight;
                expected[3] += weight;
            }
            ValueCell.update(values.elements, indices); const frame = await draw(iterations);
            assert(expected.every((v, c) => Math.abs(frame.array[offset + c] - v * 255) <= 1), `Dual peeling with ${iterations} iterations must match analytical front-to-back alpha composition: ${Array.from(frame.array.subarray(offset, offset + 4))}, expected ${expected.map(v => v * 255)}.`);
            ValueCell.update(values.elements, reverse); const reversed = await draw(iterations);
            assert(reversed.array.every((v, i) => v === frame.array[i]), 'Dual peeling must not depend on triangle drawing order.');
            assert((await renderer.pick(128, 128))?.id.groupId === 0, 'Peeling iteration limits must preserve nearest molecular selection.');
            if (iterations === 3) complete = frame;
            if (iterations === 10) assert(frame.array.every((v, i) => v === complete!.array[i]), 'Empty later peeling iterations must preserve converged color and alpha.');
        }
        for (const separation of [0, 0.0001]) {
            const pairPositions = positions.slice(0, 18); for (let i = 11; i < 18; i += 3) pairPositions[i] = 6 - separation;
            const pairMesh = Mesh.create(pairPositions, new Uint32Array([0, 1, 2, 3, 4, 5]), normals.slice(0, 18), groups.slice(0, 6), 6, 2);
            const pairValues = Mesh.Utils.createValuesSimple(pairMesh, shapeProps, Color(0xffffff), 1);
            ValueCell.update(pairValues.dColorType, 'group'); ValueCell.update(pairValues.tColor, { array: new Uint8Array([255, 0, 0, 0, 255, 0]), width: 2, height: 1 });
            const pair = createRenderObject('mesh', pairValues, Mesh.Utils.createRenderableState(shapeProps), -1);
            renderer.setDpoitIterations(1); renderer.render([pair], camera, props, true); const first = await renderer.readPixels();
            ValueCell.update(pairValues.elements, new Uint32Array([3, 4, 5, 0, 1, 2])); renderer.render([pair], camera, props, true);
            assert((await renderer.readPixels()).array.every((v, i) => v === first.array[i]), 'Coplanar and closely spaced peeled layers must remain independent of triangle order.');
            const expected = separation === 0 ? [128, 128, 0, 128] : [128, 64, 0, 191];
            assert(expected.every((v, c) => Math.abs(first.array[offset + c] - v) <= 1), 'Peeling must preserve MAX-blended coplanar ties and distinguish full-precision adjacent depths.');
        }
        const opaqueProps = { ...shapeProps, alpha: 1 };
        const opaqueMesh = Mesh.create(new Float32Array([-4, -4, 8, 4, -4, 8, 0, 4, 8]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
        const opaque = createRenderObject('mesh', Mesh.Utils.createValuesSimple(opaqueMesh, opaqueProps, Color(0x00ff00), 1), Mesh.Utils.createRenderableState(opaqueProps), -1);
        renderer.render([shape, opaque], camera, props, true); const occluded = await renderer.readPixels();
        assert(occluded.array[offset] === 0 && occluded.array[offset + 1] === 255 && occluded.array[offset + 2] === 0 && occluded.array[offset + 3] === 255, 'Peeling must exclude every transparent fragment behind opaque geometry.');
        const pass = new WebGPUImagePass(context, camera, () => [shape], { dpoitIterations: 3, renderer: props, transparentBackground: true, cameraHelper: { axes: { name: 'off', params: {} } }, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false }, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } }, undefined, { transparency: () => renderer.transparencyMode });
        try {
            const exported = await pass.getImageData(RuntimeContext.Synchronous, 256, 256), straight = new Uint8ClampedArray(complete!.array);
            for (let i = 0; i < straight.length; i += 4) if (straight[i + 3] > 0 && straight[i + 3] < 255) for (let c = 0; c < 3; c++) straight[i + c] = straight[i + c] * 255 / straight[i + 3];
            assert(exported.data.every((v, i) => v === straight[i]), 'Depth-peeled screenshot exports must retain layer count, color and straight alpha.');
        } finally { await pass.dispose(); }
        renderer.setDpoitIterations(3);
        renderer.render([shape], camera, props, true, 1, undefined, undefined, undefined, { ...PD.getDefaultValues(MultiSampleParams), mode: 'on', sampleLevel: 2 });
        const sampled = await renderer.readPixels();
        assert([0, 1, 2, 3].every(c => Math.abs(sampled.array[offset + c] - complete!.array[offset + c]) <= 1), 'Supersampling must retain depth-peeled colors and transparency.');
    } finally { renderer.setTransparency('blended'); renderer.setDpoitIterations(2); }
}

async function verifyTracingInput(context: WebGPUContext, renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const mesh = Mesh.create(new Float32Array([-4, -4, 0, 4, -4, 0, 0, 4, 0, -4, -4, -2, 4, -4, -2, 0, 4, -2]), new Uint32Array([0, 1, 2, 3, 5, 4]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, -1]), new Float32Array(6), 6, 2);
    const meshProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true, density: 0.7, emissive: 0.25 };
    const values = Mesh.Utils.createValuesSimple(mesh, meshProps, Color(0x804020), 1);
    const object = createRenderObject('mesh', values, Mesh.Utils.createRenderableState(meshProps), -1);
    const read = async (texture: GPUTexture, x = 128, y = 128) => {
        const bytes = await context.readTexture(texture, x, y, 1, 1, 8);
        return Array.from(new Uint16Array(bytes.buffer), fromHalfFloat);
    };
    const readDepth = async (texture: GPUTexture, x = 128, y = 128) => {
        const bytes = await context.readTexture(texture, 0, 0, texture.width, texture.height);
        return new Float32Array(bytes.buffer)[y * texture.width + x];
    };
    const original = camera.getSnapshot();
    try {
        let input = renderer.renderTracingInput([object], camera, props);
        const shaded = await read(input.shaded), normal = await read(input.normal), albedo = await read(input.albedo);
        const expected = [128 / 255, 64 / 255, 32 / 255];
        assert(albedo.every((v, i) => Math.abs(v - (i === 3 ? 0.7 : expected[i])) < 0.001), 'Tracing albedo and density must retain their native material values.');
        assert(shaded.every((v, i) => Math.abs(v - (i === 3 ? 1 : expected[i] * 1.25 * props.exposure)) < 0.001), 'Tracing shaded colors must include material emission and exposure.');
        assert(Math.abs(normal[0]) < 0.001 && Math.abs(normal[1]) < 0.001 && Math.abs(normal[2] + 1) < 0.001 && Math.abs(normal[3] - 0.25) < 0.001, 'Tracing normals must retain the inward view-space convention and independent emission.');
        const front = await readDepth(input.depth), back = await readDepth(input.backDepth);
        assert(front > 0 && front < back && back < 1, 'Tracing must capture the nearest and farthest opaque depth for automatic thickness.');
        const tracer = await WebGPUTracing.create(context), tracing = PD.getDefaultValues(TracingParams);
        tracing.rendersPerFrame = [2, 2]; tracing.shadowEnable = false;
        try {
            const tracePixel = async (iteration: number) => {
                const output = tracer.render(camera, props, tracing, input, iteration);
                const bytes = await context.readTexture(output, 128, 128, 1, 1, 16);
                return Array.from(new Float32Array(bytes.buffer));
            };
            for (let iteration = 0; iteration < 3; iteration++) {
                const traced = await tracePixel(iteration);
                assert(traced.every((v, i) => Math.abs(v - (i === 3 ? 1 / (iteration + 1) : shaded[i] + albedo[i] * normal[3])) < 0.00001), 'Isolated diffuse surfaces must retain direct lighting plus primary emission and accumulate analytically.');
            }
            const restarted = await tracePixel(0);
            assert(Math.abs(restarted[3] - 1) < 0.000001, 'Restarting tracing must discard previous accumulation.');
            const background = await context.readTexture(tracer.render(camera, props, tracing, input, 1), 0, 0, 1, 1, 16);
            assert(new Float32Array(background.buffer).every(v => v === 0), 'Tracing background pixels must remain empty.');
            tracing.shadowEnable = true;
            const shadowRenderer = { ...props, light: [], ambientIntensity: 1 };
            const ambientOutput = tracer.render(camera, shadowRenderer, tracing, input, 0);
            const ambient = new Float32Array((await context.readTexture(ambientOutput, 128, 128, 1, 1, 16)).buffer);
            assert(ambient.every((v, i) => Math.abs(v - restarted[i]) < 0.00001), 'Ambient-only illumination must preserve ray color when soft shadows are enabled.');
            const directOutput = tracer.render(camera, props, tracing, input, 0);
            const direct = new Float32Array((await context.readTexture(directOutput, 128, 128, 1, 1, 16)).buffer);
            assert(direct.every((v, i) => Math.abs(v - restarted[i]) < 0.00001), 'An isolated plane must retain its unobstructed directional irradiance with soft shadows enabled.');
            const compose = await WebGPUIlluminationCompose.create(context), illumination = PD.getDefaultValues(IlluminationParams);
            const traced = tracer.render(camera, props, tracing, input, 0), composeProps = { ...props, backgroundColor: Color(0x123456) };
            const composed = async (transparent: boolean, x = 128, y = 128, source = traced, iteration = 0, fullSampling = false) => {
                const encoder = context.device.createCommandEncoder();
                const texture = compose.render(encoder, source, input, camera, illumination, composeProps, transparent, iteration, fullSampling);
                context.device.queue.submit([encoder.finish()]);
                const bytes = await context.readTexture(texture, x, y, 1, 1);
                if (context.format.startsWith('bgra')) { const r = bytes[0]; bytes[0] = bytes[2]; bytes[2] = r; }
                return bytes;
            };
            try {
                illumination.denoise = false;
                const color = await composed(true);
                assert(color.every((v, i) => Math.abs(v - Math.round(Math.min(1, restarted[i]) * 255)) <= 1), 'Illumination composition must convert traced HDR color to opaque display pixels.');
                assert((await composed(true, 0, 0)).every(v => v === 0), 'Transparent illumination composition must retain empty background pixels.');
                illumination.denoise = true;
                assert((await composed(true)).every((v, i) => Math.abs(v - color[i]) <= 1), 'Normal-guided denoising must preserve constant lighting on flat surfaces.');
                const synthetic = context.device.createTexture({ size: [traced.width, traced.height], format: 'rgba32float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
                const normalTexture = context.device.createTexture({ size: [traced.width, traced.height], format: 'rgba16float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
                const pixels = new Float32Array(traced.width * traced.height * 4), normalPixels = new Uint16Array(pixels.length), originalNormal = input.normal;
                for (let i = 0; i < pixels.length; i += 4) { pixels.set([0.4, 0.2, 0.3, 1], i); normalPixels.set([0, 0, toHalfFloat(-1), 0], i); }
                pixels[(128 * traced.width + 128) * 4] = 1;
                context.device.queue.writeTexture({ texture: synthetic }, pixels, { bytesPerRow: traced.width * 16 }, [traced.width, traced.height]);
                const reference = (threshold: number, split = false) => {
                    let spatial = 0;
                    for (let x = -6; x <= 6; x++) for (let y = -6; y <= 6; y++) if ((x !== 0 || y !== 0) && (!split || x >= 0)) spatial += Math.exp(-(x * x + y * y) / 18);
                    const neighbors = spatial * Math.exp(-0.36 / (2 * threshold * threshold));
                    return Math.round((1 + 0.4 * neighbors) / (1 + neighbors) * 255);
                };
                try {
                    const early = await composed(true, 128, 128, synthetic, 0), final = await composed(true, 128, 128, synthetic, 16);
                    assert(Math.abs(early[0] - reference(illumination.denoiseThreshold[1])) <= 1 && Math.abs(final[0] - reference(illumination.denoiseThreshold[0])) <= 1, 'Illumination denoising must match the independent Gaussian/color-weight reference and adapt its threshold with convergence.');
                    assert(final[0] > early[0] && (await composed(true, 128, 128, synthetic, 0, true))[0] === final[0], 'Full supersampling must use the final denoise threshold.');
                    for (let y = 0; y < traced.height; y++) for (let x = 0; x < 128; x++) normalPixels.set([0, toHalfFloat(-1), 0, 0], (y * traced.width + x) * 4);
                    context.device.queue.writeTexture({ texture: normalTexture }, normalPixels, { bytesPerRow: traced.width * 8 }, [traced.width, traced.height]);
                    input.normal = normalTexture;
                    assert(Math.abs((await composed(true, 128, 128, synthetic, 0))[0] - reference(illumination.denoiseThreshold[1], true)) <= 1, 'Normal-guided denoising must exclude samples across a perpendicular surface boundary.');
                    illumination.denoise = false;
                    assert((await composed(true, 128, 128, synthetic, 0))[0] === 255, 'Disabling illumination denoising must preserve the original pixel.');
                    illumination.denoise = true;
                } finally { input.normal = originalNormal; synthetic.destroy(); normalTexture.destroy(); }

                camera.setState({ fog: 100 }, 0); camera.update();
                const t = Math.max(0, Math.min(1, (20 * camera.scale - camera.fogNear) / (camera.fogFar - camera.fogNear)));
                const fog = t * t * (3 - 2 * t), fogged = await composed(true);
                assert(fog > 0.1 && fog < 0.9 && fogged.every((v, i) => Math.abs(v - Math.round(Math.min(1, restarted[i] * (1 - fog)) * 255)) <= 1), 'Transparent traced lighting must apply fog once to both color and alpha.');
                const solid = await composed(false), background = Color.toRgbNormalized(composeProps.backgroundColor);
                assert(solid.every((v, i) => Math.abs(v - (i === 3 ? 255 : Math.round(Math.min(1, restarted[i] * (1 - fog) + background[i] * fog) * 255))) <= 1), 'Solid traced lighting must mix with the configured fog color while retaining opaque alpha.');
                camera.setState(original, 0); camera.update();
                const backgroundPixel = await composed(false, 0, 0);
                assert(backgroundPixel.every((v, i) => v === (i === 3 ? 255 : Math.round(background[i] * 255))), 'Solid illumination composition must retain its configured background.');
            } finally { compose.dispose(); }


        } finally { tracer.dispose(); }

        const projectedFront = camera.project(Vec4(), Vec3()), projectedBack = camera.project(Vec4(), Vec3.create(0, 0, -2));
        assert(Math.abs(front - projectedFront[2]) < 0.000001 && Math.abs(back - projectedBack[2]) < 0.000001, 'Tracing front/back depth must match analytic camera projection.');
        camera.setState({ fog: 100 }, 0); camera.update(); input = renderer.renderTracingInput([object], camera, props);
        assert((await read(input.shaded)).every((v, i) => v === shaded[i]), 'Illumination inputs must be independent of fog so composition applies it once.');
        camera.setState(original, 0); camera.update();
        ValueCell.update(values.alpha, 0.5); input = renderer.renderTracingInput([object], camera, props);
        assert((await read(input.shaded)).every(v => v === 0) && (await readDepth(input.depth)) === 1 && (await readDepth(input.backDepth)) === 0, 'Translucent geometry must leave opaque tracing buffers empty.');
        ValueCell.update(values.alpha, 1); object.state.visible = false; input = renderer.renderTracingInput([object], camera, props);
        assert((await read(input.albedo)).every(v => v === 0), 'Hidden geometry must be excluded from illumination inputs.');
        object.state.visible = true;
        const viewport = { ...camera.viewport };
        Object.assign(camera.viewport, { x: 10, y: 15, width: 96, height: 80 }); camera.update();
        input = renderer.renderTracingInput([object], camera, props, 128, 112);
        assert(input.shaded.width === 128 && input.depth.height === 112 && (await read(input.albedo, 58, 57))[0] > 0.49, 'Tracing buffers must resize and render offset viewports.');
        assert((await read(input.shaded, 0, 0)).every(v => v === 0) && await readDepth(input.depth, 0, 0) === 1, 'Tracing buffers must retain empty pixels outside the viewport.');
        Object.assign(camera.viewport, viewport); camera.update();
    } finally { camera.setState(original, 0); camera.update(); }
}

async function verifyMultiSample(renderer: WebGPURenderer, context: WebGPUContext, camera: Camera, object: GraphicsRenderObject, props: RendererProps) {
    const post = { ...PD.getDefaultValues(PostprocessingParams), enabled: false };
    const options = { ...PD.getDefaultValues(MultiSampleParams), mode: 'on' as const, reuseOcclusion: false };
    const draw = (changed = true, temporal = false, force = false) => renderer.render([object], camera, props, true, 1, undefined, post, undefined, { ...options, mode: temporal ? 'temporal' : 'on' }, changed, force);
    renderer.render([object], camera, props, true, 1, undefined, post);
    const baseline = await renderer.readPixels(), pick = await renderer.pick(128, 128);
    const original = { ...camera.viewOffset };
    for (const level of [0, 1, 2, 3, 4, 5]) {
        options.sampleLevel = level;
        const reference = new Float32Array(baseline.array.length), offsets = JitterVectors[level];
        for (let i = 0; i < offsets.length; i++) {
            camera.viewOffset.enabled = true;
            Camera.setViewOffset(camera.viewOffset, camera.viewport.width, camera.viewport.height, offsets[i][0], offsets[i][1], camera.viewport.width, camera.viewport.height); camera.update();
            renderer.render([object], camera, props, true, 1, undefined, post);
            const pixels = await renderer.readPixels(), weight = 1 / offsets.length + (1 / 32) * (-0.5 + (i + 0.5) / offsets.length);
            for (let j = 0; j < reference.length; j++) reference[j] += pixels.array[j] * weight;
        }
        Object.assign(camera.viewOffset, original); camera.update(); draw();
        const sampled = await renderer.readPixels();
        assert(sampled.array.every((v, i) => Math.abs(v - reference[i]) <= 2), `Native sample level ${level} must match independently accumulated jittered frames.`);
        assert(!renderer.multiSampleNeedsFrame && JSON.stringify(camera.viewOffset) === JSON.stringify(original), 'Full sampling must complete and restore the camera.');
        assert(JSON.stringify((await renderer.pick(128, 128))?.id) === JSON.stringify(pick?.id), 'Supersampling must preserve canonical picking.');
        assert(sampled.array.every((v, i) => i % 4 === 3 || v <= sampled.array[i - i % 4 + 3]), 'Sampled transparent colors must remain premultiplied.');
        if (level === 0) assert(sampled.array.every((v, i) => Math.abs(v - baseline.array[i]) <= 1), 'Level zero must reproduce the unjittered frame.');
        else assert(sampled.array.some((v, i) => i % 4 === 3 && v > 0 && v < 255), 'Subpixel jitter must create fractional edge coverage.');
    }
    options.sampleLevel = 2;
    draw(true, true);
    assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Temporal sampling must begin with the canonical frame.');
    for (let i = 0; i < 4; i++) { assert(renderer.multiSampleNeedsFrame, 'Temporal refinement must remain pending until all samples finish.'); draw(false, true); }
    assert(!renderer.multiSampleNeedsFrame, 'Temporal refinement must converge.');
    const temporal = await renderer.readPixels(); draw(false, true);
    assert((await renderer.readPixels()).array.every((v, i) => v === temporal.array[i]), 'Completed temporal sampling must retain its frame.');
    draw(true, true); assert(renderer.multiSampleNeedsFrame, 'Scene changes must restart temporal sampling.');
    draw(true, true, true); assert(!renderer.multiSampleNeedsFrame, 'Flicker reduction must support completing temporal sampling immediately.');
    const full = await renderer.readPixels();
    const image = new WebGPUImagePass(context, camera, () => [object], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: post, multiSample: { ...options, mode: 'temporal' } });
    const exported = await image.getImageData(RuntimeContext.Synchronous, baseline.width, baseline.height);
    const straight = new Uint8ClampedArray(full.array);
    for (let i = 0; i < straight.length; i += 4) if (straight[i + 3] > 0 && straight[i + 3] < 255) for (let c = 0; c < 3; c++) straight[i + c] = straight[i + c] * 255 / straight[i + 3];
    assert(exported.data.every((v, i) => v === straight[i]), 'Temporal screenshot exports must complete all samples before converting to straight alpha.');
    await image.dispose();
    const ao = PD.getDefaultValues(PostprocessingParams);
    ao.antialiasing = { name: 'off', params: {} }; ao.bloom = { name: 'off', params: {} };
    const encoderPrototype = Object.getPrototypeOf(context.device.createCommandEncoder());
    const beginRenderPass = encoderPrototype.beginRenderPass;
    let aoPasses = 0;
    encoderPrototype.beginRenderPass = function (descriptor: GPURenderPassDescriptor) {
        if (descriptor.label === 'molstar-postprocess-ssao') aoPasses++;
        return beginRenderPass.call(this, descriptor);
    };
    try {
        renderer.render([object], camera, props, true, 1, undefined, ao, undefined, { ...options, reuseOcclusion: false });
        assert(aoPasses === 5, 'Disabling occlusion reuse must compute SSAO for the baseline and every jitter sample.');
        const independentAo = await renderer.readPixels(); aoPasses = 0;
        renderer.render([object], camera, props, true, 1, undefined, ao, undefined, { ...options, reuseOcclusion: true });
        assert(aoPasses === 1, 'Occlusion reuse must compute SSAO once per supersampled frame.');
        const reusedAo = await renderer.readPixels();
        assert(reusedAo.array.every((v, i) => Math.abs(v - independentAo.array[i]) <= 2), 'Reused SSAO must preserve flat-surface color and edge coverage.');
        aoPasses = 0;
        for (let i = 0; i < 5; i++) renderer.render([object], camera, props, true, 1, undefined, ao, undefined, { ...options, mode: 'temporal', reuseOcclusion: true }, i === 0);
        assert(aoPasses === 1 && !renderer.multiSampleNeedsFrame, 'Temporal refinement must retain SSAO across successive frames.');
    } finally { encoderPrototype.beginRenderPass = beginRenderPass; }
    camera.viewOffset.enabled = true;
    Camera.setViewOffset(camera.viewOffset, camera.viewport.width * 2, camera.viewport.height * 2, camera.viewport.width / 2, camera.viewport.height / 2, camera.viewport.width, camera.viewport.height); camera.update();
    const crop = { ...camera.viewOffset }; options.sampleLevel = 0;
    renderer.render([object], camera, props, true, 1, undefined, post); const cropped = await renderer.readPixels(); draw();
    assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - cropped.array[i]) <= 1) && JSON.stringify(camera.viewOffset) === JSON.stringify(crop), 'Sampling must preserve an existing cropped view offset.');
    Object.assign(camera.viewOffset, original); camera.update();
    renderer.render([object], camera, props);
}

async function verifyPlugin() {
    const container = document.createElement('div');
    container.style.cssText = 'position:relative;width:320px;height:240px';
    const canvas = document.createElement('canvas');
    canvas.width = 320; canvas.height = 240;
    canvas.style.cssText = 'width:320px;height:240px';
    container.appendChild(canvas); document.body.appendChild(container);
    const spec = DefaultPluginSpec();
    spec.behaviors.push(PluginSpec.Behavior(DebugHelpers));
    const plugin = new PluginContext(spec);
    assert(plugin.config.get(PluginConfig.General.RenderingBackend) === 'webgpu', 'New plugins must select WebGPU by default.');
    const errors: string[] = [];
    try {
        await plugin.init();
        assert(await plugin.initViewerAsync(canvas, container), 'WebGPU plugin initialization must succeed.');
        assert(plugin.canvas3dContext?.props.transparency === 'wboit' && plugin.canvas3dContext.webgpuRenderer?.transparencyMode === 'wboit', 'Default native canvas contexts must use weighted blended transparency.');
        plugin.canvas3dContext.setProps({ transparency: 'blended' });
        assert(plugin.canvas3dContext.webgpuRenderer!.transparencyMode as string === 'blended', 'Canvas context updates must change the actual native transparency pass.');
        plugin.canvas3dContext.setProps({ transparency: 'wboit' });
        for (let i = 0; i < 20 && plugin.canvas3d?.debugRegistry?.scenes.length !== 5; i++) await new Promise(resolve => setTimeout(resolve, 10));
        assert(plugin.canvas3d?.debugRegistry instanceof WebGPUDebugRegistry && plugin.canvas3d.debugRegistry.scenes.length === 5, 'The debug extension must register all five helpers without WebGL.');
        assert(plugin.canvas3d?.webgpu && !plugin.canvas3d.webgl, 'The WebGPU plugin must not acquire a WebGL renderer.');
        plugin.canvas3d.webgpu.errors.subscribe(error => errors.push(error.message));
        plugin.events.canvas3d.recreated.subscribe(() => plugin.canvas3d?.webgpu?.errors.subscribe(error => errors.push(error.message)));
        plugin.animationLoop.stop({ noDraw: true });
        const initialCanvas = plugin.canvas3d, initialRenderer = plugin.canvas3dContext!.webgpuRenderer!;
        initialCanvas.resume(); initialCanvas.tick(now());
        const helperTexture = (initialRenderer as unknown as { selectionTexture: GPUTexture }).selectionTexture;
        const helperIds = new Uint32Array((await initialCanvas.webgpu!.readTexture(helperTexture, 0, 0, helperTexture.width, helperTexture.height, 16)).buffer);
        const axisPixel = helperIds.findIndex((v, i) => i % 4 === 2 && v === CameraHelperAxis.X && helperIds[i - 2] > 0) - 2;
        assert(axisPixel >= 0, 'Default native plugins must render pickable camera axes before loading molecular geometry.');
        const axisPick = await initialRenderer.pick(axisPixel / 4 % helperTexture.width, Math.floor(axisPixel / 4 / helperTexture.width));
        assert(axisPick && isCameraAxesLoci(initialCanvas.getLoci(axisPick.id).loci), 'Canvas picking must route camera-axis IDs to existing orientation behavior loci.');
        await verifyHandleCanvas(plugin);
        const pdb = [
            'ATOM      1  N   GLY A   1      -1.200   0.000   0.000  1.00 20.00           N  ',
            'ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  ',
            'ATOM      3  C   GLY A   1       1.200   0.000   0.000  1.00 20.00           C  ',
            'ATOM      4  O   GLY A   1       2.000   0.800   0.000  1.00 20.00           O  ',
            'END',
        ].join('\n');
        const data = await plugin.builders.data.rawData({ data: pdb, label: 'WebGPU verification molecule' });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        const ballAndStick = await plugin.builders.structure.representation.addRepresentation(structure, { type: 'ball-and-stick', color: 'element-symbol' });
        const native = plugin.canvas3d;
        native.resume(); native.tick(now());
        assert(native.camera.state.radius > 1, 'The first asynchronously created molecular geometry must fit the native camera automatically.');
        assert(native.getRenderObjects().length > 0 && native.stats.drawCount > 0, 'PDB ball-and-stick representations must reach the WebGPU renderer.');
        assert(native.getRenderObjects().some(o => o.type === 'spheres'), 'Molecular atoms must select native sphere data when WebGPU is available.');
        for (const key of ['sceneBoundingSpheres', 'visibleSceneBoundingSpheres', 'objectBoundingSpheres', 'instanceBoundingSpheres']) await verifyDebugCanvas(plugin, key);
        const renderer = plugin.canvas3dContext!.webgpuRenderer!;
        const pixels = await renderer.readPixels();
        assert(pixels.array.some((v, i) => i % 4 < 3 && v < 150), 'The PDB structure must produce molecular pixels.');
        let picked;
        for (let y = 0; y < pixels.height && !picked; y += 4) {
            for (let x = 0; x < pixels.width && !picked; x += 4) {
                const o = (y * pixels.width + x) * 4;
                if (Math.min(pixels.array[o], pixels.array[o + 1], pixels.array[o + 2]) < 150) picked = await renderer.pick(x, y);
            }
        }
        assert(picked && !isEmptyLoci(native.getLoci(picked.id).loci), 'WebGPU picking must resolve to molecular loci.');
        assert(native.props.multiSample.mode === 'temporal', 'Default WebGPU canvases must enable temporal multisampling.');
        native.requestDraw(); native.tick(now());
        assert(renderer.multiSampleNeedsFrame, 'A changed scene must start temporal refinement.');
        for (let i = 0; i < 4; i++) native.tick(now());
        assert(!renderer.multiSampleNeedsFrame, 'Default temporal sampling must finish after its four samples.');
        const converged = await renderer.readPixels(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === converged.array[i]), 'An idle canvas must retain its converged frame.');
        native.mark({ repr: ballAndStick.obj!.data.repr, loci: EveryLoci }, MarkerAction.Highlight); native.tick(now());
        assert(!renderer.multiSampleNeedsFrame, 'Marker changes with flicker reduction must complete sampling in one frame.');
        native.setProps({ multiSample: { reduceFlicker: false } }); native.tick(now());
        for (let i = 0; i < 4; i++) native.tick(now());
        native.mark({ repr: ballAndStick.obj!.data.repr, loci: EveryLoci }, MarkerAction.RemoveHighlight); native.tick(now());
        assert(renderer.multiSampleNeedsFrame, 'Disabling flicker reduction must leave marker changes to temporal refinement.');
        // These subsequent material/geometry comparisons isolate single-frame rendering.
        native.setProps({ multiSample: { mode: 'off' } }); native.tick(now());
        const beforeIllumination = await renderer.readPixels();
        native.setProps({ illumination: { enabled: true, maxIterations: 2, denoise: false, steps: 8, bounces: 1 } }); native.tick(now());
        assert(renderer.illuminationNeedsFrame && Number(renderer.illuminationProgress) === 1, 'Canvas settings must activate native progressive illumination.');
        for (let i = 0; i < 3; i++) native.tick(now());
        assert(!renderer.illuminationNeedsFrame && Number(renderer.illuminationProgress) === 4, 'Canvas ticks must finish the configured illumination iterations.');
        const illuminatedPass = native.getImagePass({});
        assert(illuminatedPass instanceof WebGPUImagePass, 'Illuminated exports must use the native image pass.');
        assert(illuminatedPass.props.illumination.enabled, 'Screenshot passes must inherit the current canvas illumination settings.');
        const illuminatedExport = await illuminatedPass.getImageData(RuntimeContext.Synchronous, 160, 120);
        assert(illuminatedExport.data.some((v, i) => i % 4 < 3 && v < 150), 'Native illumination screenshots must retain visible molecular geometry.');
        await illuminatedPass.dispose();
        native.setProps({ illumination: { enabled: false } }); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === beforeIllumination.array[i]), 'Toggling canvas illumination off must restore the previous molecular frame.');
        await verifyStereoCanvas(plugin);

        const originalFog = native.props.cameraFog; native.setProps({ cameraFog: { name: 'off', params: {} } });
        assert(ballAndStick?.obj, 'The PDB ball-and-stick representation must expose its update API.');
        await ballAndStick.obj.data.repr.createOrUpdate({ visuals: ['intra-bond', 'inter-bond'], ignoreLight: true }).run(); native.requestDraw(); native.tick(now());
        const defaultBonds = await renderer.readPixels();
        const defaultBondThemes = native.getRenderObjects().filter(o => o.type === 'cylinders').map(o => createWebGPUThemes(o, createWebGPUGeometry(o)));
        assert(defaultBonds.array.some((v, i) => i % 4 < 3 && v < 150), `Native bonds must render alone. ${JSON.stringify(native.camera.state)} ${JSON.stringify(native.stats)}`);
        await ballAndStick.obj.data.repr.createOrUpdate({ colorMode: 'interpolate' }).run(); native.requestDraw(); native.tick(now());
        const interpolatedMolecule = await renderer.readPixels();
        assert(native.getRenderObjects().some(o => o.type === 'cylinders' && 'dDualColor' in o.values && o.values.dDualColor.ref.value), 'Molecular bond interpolation must enable endpoint themes on native cylinders.');
        assert(native.getRenderObjects().filter(o => o.type === 'cylinders').some((o, i) => createWebGPUThemes(o, createWebGPUGeometry(o)).some((v, j) => v !== defaultBondThemes[i][j])), 'Molecular color updates must rebuild the native endpoint themes.');
        assert(interpolatedMolecule.array.some((v, i) => v !== defaultBonds.array[i]), 'Changing molecular bond color mode must update the rendered PDB structure.');
        const interpolatedExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(interpolatedExport.data.some((v, i) => i % 4 === 0 && interpolatedExport.data[i + 3] > 0 && Math.abs(v - interpolatedExport.data[i + 2]) > 20), 'Molecular screenshot exports must retain endpoint bond colors.');
        showImage('PDB interpolated bond export', interpolatedExport);
        await ballAndStick.obj.data.repr.createOrUpdate({ colorMode: 'default', visuals: ['element-sphere', 'intra-bond', 'inter-bond'], ignoreLight: false }).run(); native.requestDraw(); native.tick(now());
        native.setProps({ cameraFog: originalFog }); native.tick(now());
        const restoredMolecule = await renderer.readPixels();
        assert(restoredMolecule.array.every((v, i) => v === pixels.array[i]), 'Restoring default bond colors must reproduce the original molecular frame.');
        const image = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(image.width === 192 && image.height === 128 && image.data[3] === 0, 'Native screenshot export must support custom size and transparency.');
        showImage('PDB ball-and-stick export', image);
        const outlineExport = await native.getImagePass({ transparentBackground: true, postprocessing: {
            ...PD.getDefaultValues(PostprocessingParams), outline: { name: 'on', params: { scale: 2, threshold: 0.33, color: Color(0xff00ff), includeTransparent: true } },
            antialiasing: { name: 'fxaa', params: { edgeThresholdMin: 0.0312, edgeThresholdMax: 0.063, iterations: 12, subpixelQuality: 0.3 } },
        } }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(outlineExport.data.some((v, i) => i % 4 === 0 && v > 150 && outlineExport.data[i + 2] > 150), 'Screenshot export must apply native outline and FXAA settings.');
        showImage('Molecule export with outline and FXAA', outlineExport);
        const moleculeAlpha = ballAndStick.obj.data.repr.props.alpha, exportBackground = native.props.renderer.backgroundColor;
        native.setProps({ renderer: { backgroundColor: Color(0x204060) }, cameraFog: { name: 'off', params: {} } });
        await ballAndStick.obj.data.repr.createOrUpdate({ alpha: 0.15 }).run(); native.requestDraw(); native.tick(now());
        const outlinePost = { ...PD.getDefaultValues(PostprocessingParams), outline: { name: 'on' as const, params: { scale: 2, threshold: 0.33, color: Color(0xff00ff), includeTransparent: true } }, antialiasing: { name: 'off' as const, params: {} } };
        const transparentOutlineExport = await native.getImagePass({ transparentBackground: true, postprocessing: outlinePost }).getImageData(RuntimeContext.Synchronous, 192, 128);
        const opaqueOutlineExport = await native.getImagePass({ transparentBackground: false, postprocessing: outlinePost }).getImageData(RuntimeContext.Synchronous, 192, 128);
        const exportRgb = Color.toRgb(Color(0x204060));
        for (let i = 0; i < opaqueOutlineExport.data.length; i += 4) {
            const alpha = transparentOutlineExport.data[i + 3] / 255;
            for (let c = 0; c < 3; c++) assert(Math.abs(opaqueOutlineExport.data[i + c] - (transparentOutlineExport.data[i + c] * alpha + exportRgb[c] * (1 - alpha))) <= 2, 'Transparent molecular outline exports must compose consistently over an opaque background.');
            assert(opaqueOutlineExport.data[i + 3] === 255, 'Opaque outlined exports must retain full background alpha.');
        }
        await ballAndStick.obj.data.repr.createOrUpdate({ alpha: moleculeAlpha }).run();
        native.setProps({ renderer: { backgroundColor: exportBackground }, cameraFog: originalFog }); native.requestDraw(); native.tick(now());
        const focusedExport = await native.getImagePass({ transparentBackground: true, postprocessing: {
            ...PD.getDefaultValues(PostprocessingParams), sharpening: { name: 'on', params: { sharpness: 0.8, denoise: true } },
            dof: { name: 'on', params: { blurSize: 5, blurSpread: 1, inFocus: 3, PPM: 1, center: 'scene-center', mode: 'sphere' } },
        } }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(focusedExport.data[3] === 0 && focusedExport.data.some((v, i) => i % 4 === 3 && v > 100), 'Export must support DOF and sharpening while preserving alpha.');
        showImage('Molecule export with DOF and sharpening', focusedExport);
        const animationAmplitudes = native.getRenderObjects().flatMap(o => 'uWiggleAmplitude' in o.values ? [o.values.uWiggleAmplitude] : []);
        assert(animationAmplitudes.length > 0, 'Molecular geometry must expose animation controls.');
        for (const amplitude of animationAmplitudes) ValueCell.update(amplitude, 1);
        const animationTime = now(); native.resetTime(animationTime);
        let frames = 0; const frameSub = native.didDraw.subscribe(() => frames++);
        native.tick(animationTime); const firstFrame = await renderer.readPixels();
        native.tick((animationTime + 1000) as now.Timestamp); const nextFrame = await renderer.readPixels();
        assert(frames >= 2 && firstFrame.array.some((v, i) => v !== nextFrame.array[i]), 'Canvas animation clock must redraw moving molecules without representation updates.');
        frameSub.unsubscribe();
        for (const amplitude of animationAmplitudes) ValueCell.update(amplitude, 0);
        let resetCalled = false;
        native.requestCameraReset({ durationMs: 0, snapshot: (scene, camera) => {
            resetCalled = true;
            assert(scene.boundingSphereVisible.radius > 0, 'Camera reset callbacks must receive current visible molecular bounds.');
            return { ...camera.getFocus(scene.boundingSphereVisible.center, scene.boundingSphereVisible.radius), up: Vec3.create(1, 0, 0) };
        } }); native.tick(now());
        assert(resetCalled && native.camera.state.up[0] === 1, 'Native camera reset callbacks must apply custom snapshot orientation.');
        native.camera.setState({ target: Vec3.create(99, 99, 99), position: Vec3.create(99, 99, 120), radius: 1 }, 0);
        native.requestCameraReset({ durationMs: 0, snapshot: { up: Vec3.create(0, 1, 0) } }); native.tick(now());
        assert(Vec3.distance(native.camera.state.target, native.boundingSphereVisible.center) < 0.001 && native.camera.state.radius > 0, 'Partial snapshots must retain automatic scene fitting.');
        native.requestDraw(); native.tick(now());
        await native.webgpu!.device.queue.onSubmittedWorkDone();
        assert(errors.length === 0, `Plugin WebGPU errors: ${errors.join('; ')}`);
        await verifyBackgroundAssets(plugin);
        const molecularResults = await verifyMolecularRepresentations(plugin);
        const volumeResults = await verifyVolume(plugin);
        await native.webgpu!.device.queue.onSubmittedWorkDone();
        assert(errors.length === 0, `Volume WebGPU errors: ${errors.join('; ')}`);
        return ['native debug extension, all four sphere categories, mesh normals, clip shapes, slice/volume edges, screenshot ownership, clearing, disabling, picking preservation and illumination exclusion', 'plugin background assets, load-triggered redraw, screenshot readiness and shared asset reference cleanup', 'plugin initialization without WebGL', 'PDB parsing and ball-and-stick representation', 'molecular stereo canvas settings, asynchronous per-eye world picking, mono screenshot compatibility and disable restoration', 'molecular bond interpolation updates and transparent screenshot exports', 'molecular loci picking', 'native screenshot export with outline and FXAA', 'canvas animation clock and continuous drawing', 'camera reset callbacks, custom orientation and partial snapshot fitting', ...molecularResults, ...volumeResults];
    } finally {
        plugin.dispose();
        container.remove();
    }
}


const DebugOff = { sceneBoundingSpheres: false, visibleSceneBoundingSpheres: false, objectBoundingSpheres: false, instanceBoundingSpheres: false, meshNormals: false, clipObjects: false, imageEdges: false, directVolumeEdges: false };

async function verifyDebugCanvas(plugin: PluginContext, key: string) {
    const canvas = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!, registry = canvas.debugRegistry;
    assert(registry instanceof WebGPUDebugRegistry, 'Native canvases must provide a CPU-owned debug registry.');
    const original = canvas.props;
    const draw = () => { canvas.requestDraw(); canvas.tick(now()); };
    const readIds = async () => {
        const texture = (renderer as unknown as { selectionTexture: GPUTexture }).selectionTexture;
        return new Uint32Array((await canvas.webgpu!.readTexture(texture, 0, 0, texture.width, texture.height, 16)).buffer);
    };
    try {
        canvas.setProps({ transparentBackground: false, multiSample: { mode: 'off' }, postprocessing: { enabled: false } });
        registry.setProps(DebugOff); draw();
        const baseline = await renderer.readPixels(), ids = await readIds();
        registry.setProps({ [key]: true }); draw();
        const debug = await renderer.readPixels();
        assert(debug.array.some((v, i) => v !== baseline.array[i]), `${key} must draw visible native debug geometry.`);
        assert((await readIds()).every((v, i) => v === ids[i]), `${key} must preserve molecular/volume picking.`);
        const image = canvas.getImagePass({ transparentBackground: false, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } });
        assert(image instanceof WebGPUImagePass, 'Debug exports must use native WebGPU.');
        try {
            const exported = await image.getImageData(RuntimeContext.Synchronous, debug.width, debug.height);
            assert(exported.data.every((v, i) => Math.abs(v - debug.array[i]) <= 1), `${key} must be included in screenshot exports.`);
        } finally { await image.dispose(); }
        const objects = registry.getRenderObjects();
        assert(objects.length > 0, `${key} must survive export disposal.`);
        registry.clear(); draw();
        assert((await renderer.readPixels()).array.every((v, i) => v === debug.array[i]), `${key} must rebuild after clearing its scene.`);
        const traced = renderer.renderTracingInput(canvas.getRenderObjects(), canvas.camera, original.renderer);
        const normal = await canvas.webgpu!.readTexture(traced.normal, 0, 0, debug.width, debug.height, 8);
        registry.setProps(DebugOff);
        const clean = renderer.renderTracingInput(canvas.getRenderObjects(), canvas.camera, original.renderer);
        assert((await canvas.webgpu!.readTexture(clean.normal, 0, 0, debug.width, debug.height, 8)).every((v, i) => v === normal[i]), `${key} must stay out of illumination inputs.`);
        draw();
        assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), `${key} must disappear when disabled.`);
    } finally {
        registry.setProps(DebugOff); canvas.setProps({ transparentBackground: original.transparentBackground, multiSample: original.multiSample, postprocessing: original.postprocessing }); draw();
    }
}

async function verifyHandleCanvas(plugin: PluginContext) {
    const canvas = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const snapshot = canvas.camera.getSnapshot(), original = canvas.props;
    try {
        canvas.camera.setState({ position: Vec3.create(0, 0, 20), target: Vec3(), up: Vec3.unitY, radius: 10, radiusMax: 10 }, 0);
        canvas.setProps({ handle: { handle: { name: 'on', params: HandleHelperParams.handle.map('on').defaultValue } }, multiSample: { mode: 'off' }, postprocessing: { enabled: false } });
        canvas.requestDraw(); canvas.tick(now());
        const texture = (renderer as unknown as { selectionTexture: GPUTexture }).selectionTexture;
        const ids = new Uint32Array((await canvas.webgpu!.readTexture(texture, 0, 0, texture.width, texture.height, 16)).buffer);
        let pick;
        for (let i = 0; i < ids.length && !pick; i += 4) {
            if (ids[i + 2] !== HandleGroup.TranslateObjectX || !ids[i]) continue;
            const candidate = await renderer.pick(i / 4 % texture.width, Math.floor(i / 4 / texture.width));
            if (candidate && isHandleLoci(canvas.getLoci(candidate.id).loci)) pick = candidate;
        }
        assert(pick, 'Enabling native canvas handles must register pickable helper loci.');
        const baseline = await renderer.readPixels();
        const markAll = { loci: EveryLoci };
        canvas.mark(markAll, MarkerAction.Highlight); canvas.tick(now());
        const highlighted = await renderer.readPixels();
        assert(ids.some((v, i) => i % 4 === 0 && v === pick.id.objectId + 1 && highlighted.array.slice(i, i + 3).some((c, k) => c !== baseline.array[i + k])), 'Global highlighting must reach handles even when camera axes are enabled.');
        canvas.mark(markAll, MarkerAction.RemoveHighlight); canvas.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Canvas helper highlights must restore cleanly.');
        const image = canvas.getImagePass({ multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } });
        assert(image instanceof WebGPUImagePass, 'Native canvases must export through WebGPU.');
        const exported = await image.getImageData(RuntimeContext.Synchronous, baseline.width, baseline.height);
        assert(exported.data.every((v, i) => v === baseline.array[i]), 'Canvas screenshot exports must include native handles and axes.');
        await image.dispose();
        canvas.setProps({ handle: { handle: { name: 'off', params: {} } } }); canvas.tick(now());
        const cleared = new Uint32Array((await canvas.webgpu!.readTexture(texture, 0, 0, texture.width, texture.height, 16)).buffer);
        assert(cleared.every((v, i) => i % 4 !== 0 || v !== pick.id.objectId + 1), 'Disabling handles must remove their picking IDs.');
    } finally {
        canvas.setProps({ handle: original.handle, multiSample: original.multiSample, postprocessing: original.postprocessing });
        canvas.camera.setState(snapshot, 0); canvas.requestDraw(); canvas.tick(now());
    }
}

async function verifyStereoCanvas(plugin: PluginContext) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const original = native.props.camera, baseline = await renderer.readPixels();
    const radius = native.camera.state.radius;
    native.setProps({ cameraClipping: { radius: 100 } }); native.tick(now());
    assert(native.camera.state.radius === radius, 'A zero-radius clipping request must preserve the fitted camera, including full viewport settings updates.');
    assert(native.props.cameraClipping.radius === 100 - Math.round(native.camera.transition.target.radius / (native.boundingSphere.radius * native.props.sceneRadiusFactor) * 100), 'Native viewport controls must report current clipping rather than the last input value.');
    const stereo = { ...DefaultStereoCameraProps, eyeSeparation: 0.1, focus: 1 };
    const originalIllumination = native.props.illumination;
    try {
        native.setProps({ camera: { helper: { axes: { name: 'off', params: {} } }, stereo: { name: 'on', params: stereo } } }); native.tick(now());
        const expected = new WebGPUStereoCamera(); expected.update(native.camera, stereo);
        const pixels = await renderer.readPixels();
        assert(pixels.array.some((v, i) => v !== baseline.array[i]), 'Canvas settings must enable two molecular stereo views.');
        for (const eye of [expected.left, expected.right]) {
            let target: { x: number, y: number, depth: number } | undefined;
            for (let y = pixels.height / 4; y < 3 * pixels.height / 4 && !target; y += 8) for (let x = eye.viewport.x; x < eye.viewport.x + eye.viewport.width && !target; x += 8) {
                const picked = await renderer.pick(x, y);
                if (picked && native.getLoci(picked.id).loci.kind === 'element-loci') target = { x, y, depth: picked.depth };
            }
            assert(target, 'Each molecular eye must expose pickable atom or bond geometry.');
            const pending = native.asyncIdentify(Vec2.create(target.x / native.input.pixelRatio, target.y / native.input.pixelRatio));
            assert(pending, 'The native canvas must provide asynchronous stereo identification.');
            let result;
            for (let i = 0; i < 60; i++) {
                const data = pending.tryGet();
                if (data !== 'pending') { result = data; break; }
                await new Promise(resolve => requestAnimationFrame(resolve));
            }
            assert(result && result.position, 'Stereo asynchronous identification must finish with a world-space position.');
            const world = Vec3.scale(Vec3(), eye.unproject(Vec3(), Vec3.create(target.x, pixels.height - target.y, target.depth)), 1 / eye.scale);
            assert(Vec3.distance(result.position, world) < 1e-5 && native.getLoci(result.id).loci.kind === 'element-loci', 'Canvas stereo picks must reconstruct positions with the corresponding eye projection.');
            const wrong = Vec3.scale(Vec3(), native.camera.unproject(Vec3(), Vec3.create(target.x, pixels.height - target.y, target.depth)), 1 / native.camera.scale);
            assert(Vec3.distance(wrong, world) > 0.01, 'The stereo picking fixture must distinguish eye reconstruction from the parent camera.');
        }
        native.setProps({ illumination: { ...originalIllumination, enabled: true, maxIterations: 0 } }); native.tick(now());
        assert(renderer.getPickingCamera(80, 120, native.camera) === native.camera && !renderer.illuminationNeedsFrame, 'Illumination must retain its established single-camera behavior while stereo is configured.');
        native.setProps({ illumination: originalIllumination }); native.tick(now());
        assert(renderer.getPickingCamera(80, 120, native.camera) !== native.camera, 'Disabling illumination must resume configured stereo eye rendering.');
        const exportPass = native.getImagePass({ transparentBackground: true });
        assert(exportPass instanceof WebGPUImagePass, 'Stereo canvas exports must use the native image pass.');
        const stereoExport = await exportPass.getImageData(RuntimeContext.Synchronous, 160, 120); await exportPass.dispose();
        native.setProps({ camera: { stereo: { name: 'off', params: {} } } }); native.tick(now());
        const monoPass = native.getImagePass({ transparentBackground: true });
        assert(monoPass instanceof WebGPUImagePass, 'Mono canvas exports must use the native image pass.');
        const monoExport = await monoPass.getImageData(RuntimeContext.Synchronous, 160, 120); await monoPass.dispose();
        assert(stereoExport.data.every((v, i) => v === monoExport.data[i]), 'Canvas stereo must retain the established mono screenshot API behavior.');
    } finally { native.setProps({ camera: original, illumination: originalIllumination }); native.tick(now()); }
    assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Disabling canvas stereo must restore its original molecular frame.');
}

async function verifyMolecularRepresentations(plugin: PluginContext) {
    const cif = await (await fetch('/fixtures/1crn.cif')).text();
    const data = await plugin.builders.data.rawData({ data: cif, label: 'Crambin representation coverage' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    let native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.setProps({ transparentBackground: true, cameraFog: { name: 'off', params: {} } });
    const results: string[] = [];
    for (const type of ['cartoon', 'backbone', 'spacefill', 'gaussian-surface', 'molecular-surface', 'line', 'point'] as const) {
        native.clear();
        const computeBefore = native.webgpu!.stats.computeDispatches, extractionBefore = native.webgpu!.stats.marchingCubesDispatches;
        const representation = await plugin.builders.structure.representation.addRepresentation(structure, { type, color: 'uniform', colorParams: { value: Color(0xff8800) } });
        native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        assert(native.getRenderObjects().some(o => o.values.drawCount.ref.value > 0), `${type} must create native render geometry.`);
        const pixels = await renderer.readPixels();
        assert(pixels.array.some((v, i) => i % 4 === 3 && v > 0), `${type} must render molecular pixels.`);
        let pick;
        for (let y = 0; y < pixels.height && !pick; y += 5) for (let x = 0; x < pixels.width && !pick; x += 5) {
            if (pixels.array[(y * pixels.width + x) * 4 + 3]) pick = await renderer.pick(x, y);
        }
        assert(pick && !isEmptyLoci(native.getLoci(pick.id).loci), `${type} must retain molecular loci picking.`);
        showPixels(`Crambin ${type}`, pixels);
        results.push(`protein ${type} rendering and loci picking`);
        if (type === 'gaussian-surface') {
            assert(representation?.obj && native.webgpu!.stats.computeDispatches > computeBefore, 'Protein Gaussian surface generation must submit native compute work.');
            assert(native.webgpu!.stats.marchingCubesDispatches > extractionBefore, 'Protein Gaussian surfaces must extract triangles on WebGPU as well as computing density.');
            for (const object of native.getRenderObjects().filter(o => o.values.drawCount.ref.value > 0)) {
                assert(object.type === 'texture-mesh', 'GPU Gaussian surfaces must use native texture meshes.');
                const geometry = value<Record<string, unknown>>(object.values, 'meta', {}).webgpuGeometry;
                assert(geometry instanceof WebGPUTextureMeshGeometry, 'Gaussian surface textures and storage vertices must originate in native compute.');
                const item = (renderer as unknown as { items: Map<number, { buffers: GPUBuffer[] }> }).items.get(object.id);
                assert(item?.buffers[0] === geometry.buffer, 'Rendering must bind the compute output directly without uploading the CPU mirror.');
            }
            const beforeRecovery = await renderer.readPixels(), oldDevice = native.webgpu!.device, source = structure.obj?.data;
            const representationCountBeforeRecovery = native.getRepresentations().length;
            const refs = [...plugin.state.data.cells.keys()].sort(), settings = native.props;
            oldDevice.destroy(); await oldDevice.lost; await plugin.webgpuRecovery;
            native = plugin.canvas3d!; renderer = plugin.canvas3dContext!.webgpuRenderer!;
            assert(native.webgpu!.device !== oldDevice && !native.webgl && structure.obj?.data === source, 'Browser recovery must obtain a fresh native device while retaining molecular data.');
            assert(JSON.stringify([...plugin.state.data.cells.keys()].sort()) === JSON.stringify(refs), 'Browser recovery must preserve state references.');
            native.requestDraw(); native.tick(now());
            const recoveredFrame = await renderer.readPixels();
            assert(native.getRepresentations().length === representationCountBeforeRecovery, 'Device recovery must not redisplay representations removed from the canvas.');
            assert(recoveredFrame.array.every((v, i) => v === beforeRecovery.array[i]), 'Recovered browser GPU surfaces must reproduce their original frame.');
            assert(native.props.transparentBackground === settings.transparentBackground && native.props.camera.mode === settings.camera.mode, 'Recovery must preserve viewer settings.');
            assert((native.debugRegistry as WebGPUDebugRegistry).scenes.length === 5, 'Recovery must restore registered native debug helpers.');
            let recoveredPick;
            for (let y = 0; y < beforeRecovery.height && !recoveredPick; y += 5) for (let x = 0; x < beforeRecovery.width && !recoveredPick; x += 5) recoveredPick = await renderer.pick(x, y);
            assert(recoveredPick && !isEmptyLoci(native.getLoci(recoveredPick.id).loci), 'Recovered GPU picking must retain molecular loci.');
            results.push('forced browser device loss, automatic GPU surface recreation, preserved molecule/state refs/settings, identical frame, native picking and debug helpers');
            const repr = representation.obj.data.repr;
            assert(structure.obj, 'The Gaussian surface must have its source structure.');
            for (const traceOnly of [false, true]) {
                const input = getStructureConformationAndRadius(structure.obj.data, repr.theme.size, { ...PD.getDefaultValues(CommonSurfaceParams), traceOnly });
                assert(OrderedSet.size(input.position.indices) === input.position.id!.length, 'Whole-structure density indices must exclude the one-past-end coordinate.');
                OrderedSet.forEach(input.position.indices, index => {
                    assert(Number.isFinite(input.radius(index)) && Number.isFinite(input.position.x[index]), 'Filtered and unfiltered density inputs must contain valid atom data.');
                });
            }
            const beforeCpu = native.webgpu!.stats.computeDispatches;
            await repr.createOrUpdate({ tryUseGpu: false }).run(); native.requestDraw(); native.tick(now());
            assert(native.webgpu!.stats.computeDispatches === beforeCpu, 'Disabling Gaussian GPU computation must honor the CPU setting.');
            const cpuSurface = await renderer.readPixels();
            assert(cpuSurface.array.some((v, i) => i % 4 === 3 && v > 0), 'CPU Gaussian grids must continue to generate visible native surfaces.');
            await repr.createOrUpdate({ tryUseGpu: true }).run(); native.requestDraw(); native.tick(now());
            assert(native.webgpu!.stats.computeDispatches > beforeCpu, 'Re-enabling Gaussian GPU computation must dispatch a fresh density grid.');
            const restored = await renderer.readPixels();
            assert(restored.array.every((v, i) => v === pixels.array[i]), 'Restoring native Gaussian computation must reproduce the original molecular surface.');
            const beforeParent = native.webgpu!.stats.marchingCubesDispatches;
            await representation.obj.data.repr.createOrUpdate({ includeParent: true }).run(); native.tick(now());
            assert(native.webgpu!.stats.marchingCubesDispatches > beforeParent, 'Parent-inclusive Gaussian surfaces must use native marching-cubes extraction before edge smoothing.');
            assert((await renderer.readPixels()).array.some((v, i) => i % 4 === 3 && v > 0), 'Parent-inclusive Gaussian surfaces must remain visible.');
            await representation.obj.data.repr.createOrUpdate({ includeParent: false }).run(); native.tick(now());
            assert(native.webgpu!.stats.marchingCubesDispatches > beforeParent, 'Disabling parent inclusion must restore native extraction.');
            assert((await renderer.readPixels()).array.every((v, i) => v === restored.array[i]), 'Restoring native extraction after parent smoothing must reproduce the original frame.');


            for (const visual of ['structure-gaussian-surface-mesh', 'gaussian-surface-wireframe', 'structure-gaussian-surface-wireframe']) {
                const beforeVisual = native.webgpu!.stats.computeDispatches;
                await repr.createOrUpdate({ visuals: [visual] }).run(); native.requestDraw(); native.tick(now());
                assert(native.webgpu!.stats.computeDispatches > beforeVisual, `${visual} must generate its density grid on WebGPU.`);
                const visualPixels = await renderer.readPixels();
                assert(visualPixels.array.some((v, i) => i % 4 === 3 && v > 0), `${visual} must render its computed molecular geometry.`);
                let visualPick;
                for (let y = 0; y < visualPixels.height && !visualPick; y += 5) for (let x = 0; x < visualPixels.width && !visualPick; x += 5) {
                    if (visualPixels.array[(y * visualPixels.width + x) * 4 + 3]) visualPick = await renderer.pick(x, y);
                }
                assert(visualPick && !isEmptyLoci(native.getLoci(visualPick.id).loci), `${visual} must retain computed atom IDs for molecular picking.`);
            }
            await repr.createOrUpdate({ visuals: ['gaussian-surface-mesh'] }).run(); native.requestDraw(); native.tick(now());
            const smoothColors = repr.props.smoothColors;
            for (const visual of ['gaussian-surface-mesh', 'structure-gaussian-surface-mesh']) {
                await repr.createOrUpdate({ visuals: [visual], smoothColors: { name: 'on', params: { resolutionFactor: 1, sampleStride: 1 } } }).run();
                native.requestDraw(); native.tick(now()); const beforeMaterial = await renderer.readPixels();
                const substanceDispatches = native.webgpu!.stats.computeDispatches;
                repr.setState({ substance: Substance('every-loci', [{ loci: EveryLoci, material: { metalness: 1, roughness: 0.1, bumpiness: 0 }, clear: false }]) });
                assert(native.webgpu!.stats.computeDispatches > substanceDispatches, 'Synchronous substance updates must submit native smoothing before returning.');
                assert(repr.renderObjects.every(o => (value<WebGPUTextureData>(o.values, 'tSubstanceGrid', undefined as any)).native), 'Smoothed substance data must remain on the GPU.');
                assert(repr.renderObjects.every(o => value<string>(o.values, 'dSubstanceType', '') === 'volumeInstance' && value<unknown>(o.values, 'tSubstanceGrid', undefined) instanceof WebGPUTextureData), 'Unit and whole-structure substance smoothing must produce native material grids without WebGL.');
                native.requestDraw(); native.tick(now()); const spatialSurface = await renderer.readPixels();
                assert(spatialSurface.array.some((v, i) => v !== beforeMaterial.array[i]), 'Smoothed spatial materials must change the rendered protein surface.');
                const materialExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
                assert(materialExport.data.some((v, i) => i % 4 === 3 && v > 0), 'Smoothed spatial materials must render in screenshot exports.');
                repr.setState({ substance: Substance.Empty }); native.requestDraw(); native.tick(now());
                assert((await renderer.readPixels()).array.every((v, i) => v === beforeMaterial.array[i]), 'Clearing a smoothed material layer must restore the original protein frame.');
                const overlayDispatches = native.webgpu!.stats.computeDispatches;
                repr.setState({
                    overpaint: Overpaint('every-loci', [{ loci: EveryLoci, color: Color(0x00cc88), clear: false }]),
                    transparency: Transparency('every-loci', [{ loci: EveryLoci, value: 0.4 }]),
                    emissive: Emissive('every-loci', [{ loci: EveryLoci, value: 0.5 }])
                });
                assert(native.webgpu!.stats.computeDispatches >= overlayDispatches + 3, 'Combined synchronous overlay updates must compute all three native grids.');
                for (const name of ['Overpaint', 'Transparency', 'Emissive']) assert(repr.renderObjects.every(o => value<WebGPUTextureData>(o.values, `t${name}Grid`, undefined as any).native), 'Smoothed overlays must be consumed as GPU textures without readback.');
                for (const name of ['Overpaint', 'Transparency', 'Emissive']) assert(repr.renderObjects.every(o => value<string>(o.values, `d${name}Type`, '') === 'volumeInstance' && value<unknown>(o.values, `t${name}Grid`, undefined) instanceof WebGPUTextureData), 'Smoothed surface overlays must use native spatial grids.');
                native.requestDraw(); native.tick(now()); const spatialOverlays = await renderer.readPixels();
                assert(spatialOverlays.array.some((v, i) => v !== beforeMaterial.array[i]), 'Combined smoothed overlays must update the native protein frame.');
                const overlayExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
                assert(overlayExport.data.some((v, i) => i % 4 === 3 && v > 0 && v < 255), 'Spatial transparency must be retained in smoothed surface exports.');
                repr.setState({ overpaint: Overpaint.Empty, transparency: Transparency.Empty, emissive: Emissive.Empty }); native.requestDraw(); native.tick(now());
                assert((await renderer.readPixels()).array.every((v, i) => v === beforeMaterial.array[i]), 'Clearing smoothed overpaint/transparency/emission layers must restore the original frame.');
                const rapidDispatches = native.webgpu!.stats.computeDispatches;
                for (const color of [0xff3300, 0x3300ff, 0x00cc88]) repr.setState({ overpaint: Overpaint('every-loci', [{ loci: EveryLoci, color: Color(color), clear: false }]) });
                assert(native.webgpu!.stats.computeDispatches >= rapidDispatches + 3, 'Rapid overlay changes must each submit native smoothing synchronously.');
                native.requestDraw(); native.tick(now());
                assert((await renderer.readPixels()).array.some((v, i) => v !== beforeMaterial.array[i]), 'The final queued overlay must remain visible after replacing earlier grid allocations.');
                repr.setState({ overpaint: Overpaint.Empty }); native.requestDraw(); native.tick(now());
                assert((await renderer.readPixels()).array.every((v, i) => v === beforeMaterial.array[i]), 'Clearing rapid overlay updates must restore the protein frame exactly.');

            }
            await repr.createOrUpdate({ smoothColors, visuals: ['gaussian-surface-mesh'] }).run();
            results.push('unit/whole-structure spatial substance/overpaint/transparency/emission smoothing, native grid generation, protein rendering, screenshot export and layer clearing');
            results.push('protein Gaussian density and surface extraction dispatch, CPU/GPU setting updates and deterministic restoration');
            results.push('unit/whole-structure Gaussian meshes and wireframes from native compute, with molecular loci');
        }
        if (type === 'molecular-surface') {
            assert(representation?.obj, 'Molecular surfaces must expose their representation.');
            const repr = representation.obj.data.repr;
            const material = { ...repr.props.material };
            const post = native.props.postprocessing;
            // Color-dependent antialiasing can change edge alpha; compare raw material coverage.
            native.setProps({ postprocessing: { ...post, antialiasing: { name: 'off', params: {} } } }); native.tick(now());
            const materialBaseline = await renderer.readPixels();
            await repr.createOrUpdate({ material: { metalness: 1, roughness: 0.2, bumpiness: 0 } }).run(); native.requestDraw(); native.tick(now());
            const metallic = await renderer.readPixels();
            assert(metallic.array.some((v, i) => v !== materialBaseline.array[i]), 'Changing molecular surface material must update native GGX shading.');
            assert(metallic.array.every((v, i) => i % 4 !== 3 || v === materialBaseline.array[i]), 'Material updates must preserve molecular surface alpha coverage.');
            showPixels('Metallic protein surface', metallic);
            const metalExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            await repr.createOrUpdate({ material }).run(); native.requestDraw(); native.tick(now());
            const matteExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(metalExport.data.some((v, i) => v !== matteExport.data[i]), 'Native screenshot exports must honor the current molecular material.');
            const originalXray = repr.props.xrayShaded;
            await repr.createOrUpdate({ xrayShaded: true }).run(); native.requestDraw(); native.tick(now());
            const xraySurface = await renderer.readPixels();
            assert(xraySurface.array.some((v, i) => i % 4 === 3 && materialBaseline.array[i] === 255 && v > 0 && v < 240), 'Molecular x-ray shading must produce angle-dependent surface opacity.');
            showPixels('Protein x-ray surface', xraySurface);
            const xrayExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(xrayExport.data.some((v, i) => i % 4 === 3 && v > 0 && v < 240), 'Native screenshot exports must retain x-ray surface opacity.');
            await repr.createOrUpdate({ xrayShaded: 'inverted' }).run(); native.requestDraw(); native.tick(now());
            const invertedSurface = await renderer.readPixels();
            assert(invertedSurface.array.some((v, i) => v !== xraySurface.array[i]), 'Inverted molecular x-ray shading must change the surface opacity pattern.');
            await repr.createOrUpdate({ xrayShaded: originalXray }).run(); native.requestDraw(); native.tick(now());
            native.setProps({ postprocessing: post }); native.tick(now());
            results.push('molecular surface material updates, x-ray/inverted opacity and transparent screenshot exports');
        }
        if (type === 'spacefill') {
            assert(representation?.obj, 'Spacefill must expose its representation update API.');
            const sphereRepr = representation.obj.data.repr;
            const post = native.props.postprocessing;
            const barePost = { ...post, occlusion: { name: 'off' as const, params: {} }, antialiasing: { name: 'off' as const, params: {} } };
            native.setProps({ postprocessing: barePost });
            await sphereRepr.createOrUpdate({ alpha: 0.8, alphaThickness: 0 }).run(); native.requestDraw(); native.tick(now());
            const translucentAtoms = await renderer.readPixels();
            await sphereRepr.createOrUpdate({ alphaThickness: 4 }).run(); native.requestDraw(); native.tick(now());
            const thinAtoms = await renderer.readPixels();
            const alphaSum = (array: ArrayLike<number>) => {
                let sum = 0; for (let i = 3; i < array.length; i += 4) sum += array[i]; return sum;
            };
            assert(alphaSum(thinAtoms.array) < alphaSum(translucentAtoms.array), 'Spacefill alpha thickness must reduce native atom coverage according to their physical radii.');
            assert(thinAtoms.array.every((v, i) => i % 4 !== 3 || v <= translucentAtoms.array[i]), 'Atom thickness must never increase alpha, including overlapping atoms.');
            assert(thinAtoms.array.some((v, i) => i % 4 === 3 && v + 10 < translucentAtoms.array[i]), 'Atom thickness must visibly reduce partially covered protein pixels.');
            let thinPick;
            for (let y = 0; y < thinAtoms.height && !thinPick; y += 5) for (let x = 0; x < thinAtoms.width && !thinPick; x += 5) {
                if (thinAtoms.array[(y * thinAtoms.width + x) * 4 + 3]) thinPick = await renderer.pick(x, y);
            }
            assert(thinPick && !isEmptyLoci(native.getLoci(thinPick.id).loci), 'Radius-dependent atom opacity must retain molecular loci selection.');
            const thinExport = await native.getImagePass({ transparentBackground: true, postprocessing: barePost }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(thinExport.data.some((v, i) => i % 4 === 3 && v > 0 && v < 200), 'Spacefill exports must preserve thickness-dependent atom alpha.');
            showImage('Protein atom thickness export', thinExport);
            await sphereRepr.createOrUpdate({ alpha: 1, alphaThickness: 0, tryUseImpostor: false }).run(); native.requestDraw(); native.tick(now());
            assert(native.getRenderObjects().every(o => o.type === 'mesh'), 'Disabling sphere impostors must select the explicit mesh representation.');
            const meshAtoms = await renderer.readPixels();
            assert(meshAtoms.array.some((v, i) => i % 4 === 3 && v > 0), 'Explicit atom meshes must remain visible on WebGPU.');
            await verifyDebugCanvas(plugin, 'meshNormals');
            await sphereRepr.createOrUpdate({ tryUseImpostor: true, visuals: ['structure-element-sphere'] }).run(); native.requestDraw(); native.tick(now());
            assert(native.getRenderObjects().some(o => o.type === 'spheres'), 'Whole-structure atom visuals must also select native sphere data.');
            const structureAtoms = await renderer.readPixels();
            let structurePick;
            for (let y = 0; y < structureAtoms.height && !structurePick; y += 5) for (let x = 0; x < structureAtoms.width && !structurePick; x += 5) {
                if (structureAtoms.array[(y * structureAtoms.width + x) * 4 + 3]) structurePick = await renderer.pick(x, y);
            }
            assert(structurePick && !isEmptyLoci(native.getLoci(structurePick.id).loci), 'Whole-structure native spheres must preserve serial atom loci.');
            await sphereRepr.createOrUpdate({ visuals: ['element-sphere'] }).run(); native.requestDraw(); native.tick(now());
            native.setProps({ postprocessing: post }); native.tick(now());
            const restoredAtoms = await renderer.readPixels();
            assert(restoredAtoms.array.every((v, i) => v === pixels.array[i]), 'Returning to the original native atom representation must reproduce the protein frame.');
            pick = undefined;
            for (let y = 0; y < restoredAtoms.height && !pick; y += 5) for (let x = 0; x < restoredAtoms.width && !pick; x += 5) {
                if (restoredAtoms.array[(y * restoredAtoms.width + x) * 4 + 3]) pick = await renderer.pick(x, y);
            }
            assert(pick, 'Restored atoms must retain picking after representation recreation.');
            results.push('native molecular spheres, thickness-dependent alpha/export, mesh/native switching and whole-structure atom loci');
            native.setProps({ postprocessing: { ...post, occlusion: { name: 'off', params: {} } } }); native.tick(now());
            const baseline = await renderer.readPixels();
            native.setProps({ postprocessing: { ...post, occlusion: { name: 'on', params: { ...PD.getDefaultValues(SsaoParams), samples: 16, radius: 2, bias: 1, blurKernelSize: 7 } } } }); native.tick(now());
            const occluded = await renderer.readPixels();
            let darkened = 0;
            for (let i = 0; i < occluded.array.length; i += 4) if (baseline.array[i + 3] > 200 && occluded.array[i] + 3 < baseline.array[i]) darkened++;
            assert(darkened > 50, 'Native SSAO must darken occluded protein pixels after depth-aware blur.');
            assert(occluded.array[3] === 0, 'SSAO composition must preserve transparent background alpha.');
            showPixels('Protein native SSAO', occluded);
            const multi = { ...PD.getDefaultValues(SsaoParams), samples: 16, bias: 1, multiScale: { name: 'on' as const, params: { levels: [{ radius: 1, bias: 1 }, { radius: 3, bias: 0.5 }], nearThreshold: 0, farThreshold: 10000 } } };
            native.setProps({ postprocessing: { ...post, occlusion: { name: 'on', params: multi } } }); native.tick(now());
            const multiscale = await renderer.readPixels();
            assert(multiscale.array.some((v, i) => v !== baseline.array[i]), 'Native multi-scale SSAO must apply its configured radii.');
            await verifyTransparentSsao(plugin, alpha => sphereRepr.createOrUpdate({ alpha }).run());
            await verifyEnvironmentSsao(plugin, alpha => sphereRepr.createOrUpdate({ alpha }).run());
            native.setProps({ postprocessing: { ...post, occlusion: { name: 'off', params: {} }, shadow: { name: 'on', params: { steps: 16, maxDistance: 3, tolerance: 0.1 } } } }); native.tick(now());
            const shadowed = await renderer.readPixels();
            let shadowPixels = 0;
            for (let i = 0; i < shadowed.array.length; i += 4) if (baseline.array[i + 3] > 200 && shadowed.array[i] + 3 < baseline.array[i]) shadowPixels++;
            assert(shadowPixels > 50, 'Screen-space directional shadows must darken occluded molecular pixels.');
            assert(shadowed.array[3] === 0 && !isEmptyLoci(native.getLoci(pick.id).loci), 'Shadows must preserve transparent backgrounds and molecular picking.');
            showPixels('Protein native shadows', shadowed);
            native.setProps({ postprocessing: { ...post, occlusion: { name: 'off', params: {} }, shadow: { name: 'on', params: { steps: 16, maxDistance: 0, tolerance: 0.1 } } } }); native.tick(now());
            const zeroShadow = await renderer.readPixels();
            assert(zeroShadow.array.every((v, i) => v === baseline.array[i]), 'Zero shadow distance must leave the image unchanged.');
            const lighting = native.props.renderer;
            native.setProps({ renderer: { light: lighting.light.map(l => ({ ...l, azimuth: (l.azimuth + 180) % 360 })) }, postprocessing: { ...post, occlusion: { name: 'off', params: {} }, shadow: { name: 'on', params: { steps: 16, maxDistance: 3, tolerance: 0.1 } } } }); native.tick(now());
            const oppositeLight = await renderer.readPixels();
            assert(oppositeLight.array.some((v, i) => v !== shadowed.array[i]), 'Rotating directional lights must move molecular shadows.');
            native.setProps({ renderer: { light: [...lighting.light, ...lighting.light.map(l => ({ ...l, azimuth: (l.azimuth + 180) % 360, color: Color(0x0088ff) }))] } }); native.tick(now());
            const multipleLights = await renderer.readPixels();
            assert(multipleLights.array.some((v, i) => v !== oppositeLight.array[i]), 'All configured light colors and directions must contribute to shadow visibility.');
            native.setProps({ renderer: { light: [] }, postprocessing: { ...post, occlusion: { name: 'off', params: {} }, shadow: { name: 'off', params: {} } } }); native.tick(now());
            const ambientOnly = await renderer.readPixels();
            native.setProps({ postprocessing: { ...post, occlusion: { name: 'off', params: {} }, shadow: { name: 'on', params: { steps: 16, maxDistance: 3, tolerance: 0.1 } } } }); native.tick(now());
            const ambientShadow = await renderer.readPixels();
            assert(ambientShadow.array.every((v, i) => v === ambientOnly.array[i]), 'Ambient-only lighting must not cast directional shadows.');
            native.setProps({ renderer: lighting });
            const shadowExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            const plainExport = await native.getImagePass({ transparentBackground: true, postprocessing: { ...native.props.postprocessing, shadow: { name: 'off', params: {} } } }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(shadowExport.data.some((v, i) => i % 4 === 0 && v + 3 < plainExport.data[i]), 'Native screenshot exports must include molecular shadows.');
            showImage('Protein shadow export', shadowExport);
            native.setProps({ camera: { mode: 'orthographic' }, postprocessing: { ...post, occlusion: { name: 'off', params: {} }, shadow: { name: 'off', params: {} } } }); native.tick(now());
            const orthoBase = await renderer.readPixels();
            native.setProps({ postprocessing: { ...native.props.postprocessing, shadow: { name: 'on', params: { steps: 16, maxDistance: 3, tolerance: 0.1 } } } }); native.tick(now());
            const orthoShadow = await renderer.readPixels();
            assert(orthoShadow.array.some((v, i) => i % 4 === 0 && v + 3 < orthoBase.array[i]), 'Screen-space shadows must also support orthographic camera depth.');
            native.setProps({ camera: { mode: 'perspective' } });
            native.setProps({ postprocessing: post }); native.tick(now());
            results.push('native SSAO, depth-aware blur, multi-scale radii and alpha preservation');
            results.push('transparent SSAO, opacity thresholds, reduced resolutions, depth pyramids, independent blur, multi-scale controls, analytical fog, picking and screenshot exports');
            results.push('native multi-light shadows, orthographic depth, zero distance, ambient lighting, alpha, picking and screenshot export');
        }
        if (type === 'line') {
            const post = native.props.postprocessing;
            native.setProps({ postprocessing: { ...post, antialiasing: { name: 'off', params: {} } } }); native.tick(now());
            const baseline = await renderer.readPixels();
            native.setProps({ postprocessing: { ...post, antialiasing: { name: 'smaa', params: { edgeThreshold: 0.1, maxSearchSteps: 16 } } } }); native.tick(now());
            const smoothed = await renderer.readPixels();
            assert(smoothed.array.some((v, i) => v !== baseline.array[i]), 'Native SMAA must smooth molecular line geometry.');
            assert(smoothed.array[3] === 0 && !isEmptyLoci(native.getLoci(pick.id).loci), 'SMAA must preserve transparent protein backgrounds and molecular picking.');
            showPixels('Protein line SMAA', smoothed);
            const exportPass = native.getImagePass({ transparentBackground: true });
            const exportSmaa = await exportPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            const exportPlain = await native.getImagePass({ transparentBackground: true, postprocessing: { ...native.props.postprocessing, antialiasing: { name: 'off', params: {} } } }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(exportSmaa.data.some((v, i) => v !== exportPlain.data[i]), 'Native image exports must execute all SMAA passes at the requested dimensions.');
            assert(exportSmaa.data[3] === 0 && exportSmaa.data.some((v, i) => i % 4 === 3 && v > 0 && v < 240), 'SMAA exports must preserve partially covered molecular lines.');
            showImage('Protein line SMAA export', exportSmaa);
            assert('dispose' in exportPass, 'WebGPU must return a native image pass with explicit resource disposal.');
            await exportPass.dispose();
            native.requestDraw(); native.tick(now());
            const afterExportDisposal = await renderer.readPixels();
            assert(afterExportDisposal.array.every((v, i) => v === smoothed.array[i]), 'Disposing an image pass must not destroy the active renderer\'s shared SMAA lookup textures.');
            native.setProps({ postprocessing: post }); native.tick(now());
            results.push('native SMAA for molecular lines and custom-size transparent image exports');
        }
    }
    for (const type of ['gaussian-surface', 'molecular-surface'] as const) {
        native.clear();
        const colored = await plugin.builders.structure.representation.addRepresentation(structure, {
            type, color: 'element-symbol', typeParams: { smoothColors: { name: 'on', params: { resolutionFactor: 1, sampleStride: 1 } } }
        });
        assert(colored?.obj, 'Colored surfaces must expose their representation.');
        const repr = colored.obj.data.repr;
        native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        const prefix = type === 'gaussian-surface' ? 'gaussian-surface' : 'molecular-surface';
        for (const visual of [`${prefix}-mesh`, `structure-${prefix}-mesh`]) {
            await repr.createOrUpdate({ visuals: [visual] }).run(); native.requestDraw(); native.tick(now());
            assert(repr.renderObjects.every(o => ['volume', 'volumeInstance'].includes(value<string>(o.values, 'dColorType', '')) && value<unknown>(o.values, 'tColorGrid', undefined) instanceof WebGPUTextureData), 'Unit and whole-structure smoothed element colors must use native spatial grids.');
            const gridFrame = await renderer.readPixels();
            assert(gridFrame.array.some((v, i) => i % 4 === 3 && v > 0), 'Native spatial color grids must render protein surfaces.');
            assert(repr.renderObjects.every(o => value<WebGPUTextureData>(o.values, 'tColorGrid', undefined as any).native), 'Surface color grids must remain on the GPU.');
            const gridPass = native.getImagePass({ transparentBackground: true });
            const exportGrid = await gridPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(gridPass instanceof WebGPUImagePass, 'Native color grid exports must expose a native image pass.');
            await gridPass.dispose(); native.requestDraw(); native.tick(now());
            assert((await renderer.readPixels()).array.every((v, i) => v === gridFrame.array[i]), "Screenshot renderer disposal must preserve the active scene's borrowed GPU color grids.");
            for (const object of repr.renderObjects) {
                const grid = value<WebGPUTextureData>(object.values, 'tColorGrid', undefined as any);
                const previous = grid.native!.texture;
                const replacement = native.webgpu!.device.createTexture({ size: [previous.width, previous.height], format: 'rgba8unorm', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC });
                const encoder = native.webgpu!.device.createCommandEncoder(); encoder.copyTextureToTexture({ texture: previous }, { texture: replacement }, [previous.width, previous.height]);
                native.webgpu!.device.queue.submit([encoder.finish()]);
                const cellVersion = object.values.tColorGrid.ref.version, textureVersion = grid.version;
                grid.loadGPU(native.webgpu!.device, replacement);
                assert(object.values.tColorGrid.ref.version === cellVersion && grid.version > textureVersion, 'Replacing a native grid must track allocation changes independently of ValueCell versions.');
            }
            native.requestDraw(); native.tick(now());
            assert((await renderer.readPixels()).array.every((v, i) => v === gridFrame.array[i]), 'Renderer caches must rebind replacement GPU grids without stale or destroyed resources.');

            assert(exportGrid.data.some((v, i) => i % 4 === 3 && v > 0), 'Spatial molecular colors must render in screenshot exports.');
            await repr.createOrUpdate({ smoothColors: { name: 'off', params: {} } }).run(); native.requestDraw(); native.tick(now());
            assert(repr.renderObjects.every(o => !value<string>(o.values, 'dColorType', '').startsWith('volume')), 'Disabling smoothing must restore group color data.');
            const groupFrame = await renderer.readPixels();
            assert(groupFrame.array.some((v, i) => v !== gridFrame.array[i]), 'Spatial interpolation must change the protein element-color transition.');
            const smoothingDispatches = native.webgpu!.stats.computeDispatches;
            await repr.createOrUpdate({ smoothColors: { name: 'on', params: { resolutionFactor: 1, sampleStride: 1 } } }).run(); native.requestDraw(); native.tick(now());
            assert(native.webgpu!.stats.computeDispatches > smoothingDispatches, 'Re-enabling molecular color smoothing must submit native compute work before the representation update finishes.');
            assert((await renderer.readPixels()).array.every((v, i) => v === gridFrame.array[i]), 'Restoring native color smoothing must reproduce the original frame.');
            if (visual === `${prefix}-mesh`) showPixels(`Spatial element colors: ${type}`, gridFrame);
        }
        results.push(`${type} unit/whole-structure native spatial color smoothing, element colors, screenshots and smoothing controls`);
    }
    return results;
}

async function verifyVolume(plugin: PluginContext) {
    const space = Tensor.Space([32, 32, 32], [2, 1, 0], Float32Array);
    const data = space.create();
    for (let z = 0; z < 32; z++) for (let y = 0; y < 32; y++) for (let x = 0; x < 32; x++) {
        data[space.dataOffset(x, y, z)] = Math.exp(-((x - 16) ** 2 + (y - 16) ** 2 + (z - 16) ** 2) / 50);
    }
    const transform = Mat4.fromScaling(Mat4(), Vec3.create(0.25, 0.25, 0.25));
    transform[12] = -4; transform[13] = -4; transform[14] = -4;
    const volume: Volume = {
        grid: { cells: Tensor.create(space, data), transform: { kind: 'matrix', matrix: transform }, stats: { min: 0, max: 1, mean: 0.1, sigma: 0.2 } },
        instances: [{ transform: Mat4.identity() }], sourceData: { kind: 'synthetic', name: 'Gaussian test density', data: undefined },
        customProperties: new CustomProperties(), _propertyData: {}, _localPropertyData: {},
    };
    const provider = plugin.representation.volume.registry.get('direct-volume');
    const params = createVolumeRepresentationParams(plugin, volume, {
        type: 'direct-volume', typeParams: { ignoreLight: true, controlPoints: [Vec2.create(0.1, 0), Vec2.create(0.3, 0.4), Vec2.create(1, 0.5)] },
        color: 'uniform', colorParams: { value: Color(0x0088ff) },
    });
    const repr = provider.factory(plugin.representation.volume.themes, provider.getParams);
    repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, params));
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.clear();
    native.setProps({ transparentBackground: true, cameraFog: { name: 'off', params: {} } });
    try {
        for (const dataType of ['byte', 'float', 'halfFloat'] as const) {
            await repr.createOrUpdate({ ...params.type.params, dataType }, volume).run();
            native.add(repr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
            const pixels = await renderer.readPixels();
            const o = (120 * pixels.width + 160) * 4;
            assert(pixels.array[o + 2] > 80 && pixels.array[o + 3] > 80, `Native ${dataType} density must produce blue translucent volume pixels.`);
            const picked = await renderer.pick(160, 120);
            assert(picked && Volume.Cell.isLoci(native.getLoci(picked.id).loci), `Native ${dataType} density picking must resolve to volume-cell loci.`);
            const coordinates = [0, 0, 0];
            space.getCoords(picked.id.groupId, coordinates);
            assert(Math.abs(coordinates[0] - 16) < 2 && Math.abs(coordinates[1] - 16) < 2, 'Volume picking must preserve original grid axis order and cell indexes.');
            if (dataType === 'float') showPixels('Native density volume', pixels);
        }
        const volumeValues = repr.renderObjects[0].values as DirectVolumeValues;
        const gridTexture = volumeValues.tGridTex.ref.value as WebGPUTextureData;
        const originalGrid = gridTexture.data;
        const resourceReference = await renderer.readPixels();
        const uploadsBeforeReload = renderer.volumeResourceStats;
        const emptyGrid = originalGrid.array.slice();
        for (let i = 3; i < emptyGrid.length; i += 4) emptyGrid[i] = 0;
        gridTexture.load({ ...originalGrid, array: emptyGrid });
        native.requestDraw(); native.tick(now());
        assert(!(await renderer.pick(160, 120)), 'Reloading density data without updating its ValueCell must invalidate the native grid and selection.');
        assert(renderer.volumeResourceStats.gridUploads === uploadsBeforeReload.gridUploads + 1, 'A source revision must upload the density grid exactly once across rendering passes.');
        assert(renderer.volumeResourceStats.transferUploads === uploadsBeforeReload.transferUploads && renderer.volumeResourceStats.paletteUploads === uploadsBeforeReload.paletteUploads, 'Density reloads must retain transfer and palette textures.');
        gridTexture.load(originalGrid); native.requestDraw(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === resourceReference.array[i]), 'Restoring density source data must reproduce the original rendered frame.');
        const transferCell = volumeValues.tTransferTex;
        const originalTransfer = transferCell.ref.value;
        const uploadsBeforeTransfer = renderer.volumeResourceStats;
        ValueCell.update(transferCell, { ...originalTransfer, array: new Uint8Array(originalTransfer.array.length) });
        native.requestDraw(); native.tick(now());
        assert(!(await renderer.pick(160, 120)), 'Updating the transfer function must remove fully transparent density from selection.');
        assert(renderer.volumeResourceStats.gridUploads === uploadsBeforeTransfer.gridUploads && renderer.volumeResourceStats.transferUploads === uploadsBeforeTransfer.transferUploads + 1, 'Transfer updates must upload only the transfer texture.');
        ValueCell.update(transferCell, originalTransfer); native.requestDraw(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === resourceReference.array[i]), 'Restoring the transfer function must reproduce the original rendered frame.');
        const savedTransferBox = transferCell.ref;
        const replacementTransfer = ValueCell.create({ ...originalTransfer, array: new Uint8Array(originalTransfer.array.length) });
        while (replacementTransfer.ref.version < savedTransferBox.version) ValueCell.update(replacementTransfer, replacementTransfer.ref.value);
        const uploadsBeforeReplacement = renderer.volumeResourceStats;
        ValueCell.set(transferCell, replacementTransfer.ref); native.requestDraw(); native.tick(now());
        assert(!(await renderer.pick(160, 120)), 'Replacing a transfer value box with an equal revision must invalidate its texture by identity.');
        assert(renderer.volumeResourceStats.transferUploads === uploadsBeforeReplacement.transferUploads + 1 && renderer.volumeResourceStats.gridUploads === uploadsBeforeReplacement.gridUploads, 'Value box replacements must upload only the changed volume texture.');
        ValueCell.set(transferCell, savedTransferBox); native.requestDraw(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === resourceReference.array[i]), 'Restoring a replaced transfer value box must reproduce the original frame.');
        const originalTransforms = volumeValues.aTransform.ref.value;
        const uploadsBeforeFailure = renderer.volumeResourceStats;
        ValueCell.update(transferCell, { ...originalTransfer, array: originalTransfer.array.slice() });
        ValueCell.update(volumeValues.aTransform, new Float32Array(originalTransforms.length));
        let rejectedSingularTransform = false;
        try {
            renderer.render(repr.renderObjects, native.camera, native.props.renderer, true);
        } catch (error) {
            rejectedSingularTransform = error instanceof Error && error.message.includes('singular');
        } finally {
            ValueCell.update(volumeValues.aTransform, originalTransforms);
            ValueCell.update(transferCell, originalTransfer);
        }
        assert(rejectedSingularTransform, 'Invalid volume transforms must reject the resource rebuild.');
        native.requestDraw(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === resourceReference.array[i]), 'A failed rebuild must preserve borrowed density textures for the corrected frame.');
        assert(renderer.volumeResourceStats.gridUploads === uploadsBeforeFailure.gridUploads, 'Failed resource updates must retain the previous uploaded density grid.');
        await verifyDebugCanvas(plugin, 'directVolumeEdges');
        const readVolumePick = async () => {
            const pending = native.asyncIdentify(Vec2.create(160, 120));
            assert(pending, 'Native volume picking must expose asynchronous readback.');
            const start = performance.now();
            for (;;) {
                const data = pending.tryGet();
                if (data !== 'pending') return data;
                assert(performance.now() - start < 5000, 'Asynchronous volume picking must complete.');
                await new Promise(resolve => setTimeout(resolve, 0));
            }
        };
        const referencePosition = (await readVolumePick())?.position;
        assert(referencePosition, 'Volume center must expose a picked molecular position.');
        await verifyNativeRayPicking(plugin);
        const volumeScaleReference = await renderer.readPixels();
        for (const scale of [0.25, 2]) {
            native.camera.scale = scale; native.requestDraw(); native.tick(now());
            const scaled = await renderer.readPixels();
            assert(scaled.array.every((v, i) => Math.abs(v - volumeScaleReference.array[i]) <= 1), 'Volume ray origins and samples must preserve molecular coordinates under camera scaling.');
            const picked = await renderer.pick(160, 120);
            assert(picked && Volume.Cell.isLoci(native.getLoci(picked.id).loci), 'Scaled volumes must retain volume-cell picking.');
            const scaledPosition = (await readVolumePick())?.position;
            assert(scaledPosition && Vec3.distance(referencePosition, scaledPosition) < 0.0001, 'Asynchronous picking must return molecular coordinates independently of camera scale.');
            const pass = native.getImagePass({ transparentBackground: true, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false } });
            const exported = await pass.getImageData(RuntimeContext.Synchronous, 320, 240);
            assert(exported.data.some((v, i) => i % 4 === 3 && v > 80), 'Scaled volume exports must retain visible density.');
            assert(pass instanceof WebGPUImagePass, 'Native volume exports must use the WebGPU image pass.');
            await pass.dispose();
        }
        native.camera.scale = 1; native.requestDraw(); native.tick(now());
        const volumeBeforeEmission = await renderer.readPixels();
        native.setProps({ illumination: { enabled: true, maxIterations: 0, denoise: false } }); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - volumeBeforeEmission.array[i]) <= 1), 'Illumination must preserve the independent transparent direct-volume ray-marching path.');
        native.setProps({ illumination: { enabled: false } }); native.tick(now());

        await verifyVolumeOutlineLayers(plugin);
        await verifyVolumeBackground(plugin);
        const pickingThreshold = native.props.renderer.pickingAlphaThreshold;
        await repr.createOrUpdate({ alpha: 0.01 }, volume).run(); native.tick(now());
        const faintDensity = await renderer.readPixels();
        assert(faintDensity.array.some((v, i) => i % 4 === 3 && v > 0), 'Below-threshold density must remain visible in the color pass.');
        assert(!(await renderer.pick(160, 120)), 'Below-threshold accumulated density opacity must be excluded from selection.');
        native.setProps({ renderer: { pickingAlphaThreshold: 0 } }); native.tick(now());
        assert((await renderer.pick(160, 120))?.id.objectId === repr.renderObjects[0].id, 'Lowering the selection threshold must make faint density selectable.');
        const selectedFaintDensity = await renderer.readPixels();
        assert(selectedFaintDensity.array.every((v, i) => v === faintDensity.array[i]), 'Volume selection threshold must preserve its rendered colors.');
        await repr.createOrUpdate({ alpha: 1 }, volume).run();
        native.setProps({ renderer: { pickingAlphaThreshold: pickingThreshold } }); native.tick(now());
        const uploadsBeforeMaterial = renderer.volumeResourceStats;
        await repr.createOrUpdate({ ignoreLight: false, material: { metalness: 0, roughness: 1, bumpiness: 0 } }, volume).run(); native.tick(now());
        const litDensity = await renderer.readPixels();
        assert(litDensity.array.some((v, i) => v !== volumeBeforeEmission.array[i]), 'Density gradients must use native physical material lighting.');
        await repr.createOrUpdate({ material: { metalness: 1, roughness: 0.2, bumpiness: 0 } }, volume).run(); native.tick(now());
        const metallicDensity = await renderer.readPixels();
        assert(metallicDensity.array.some((v, i) => v !== litDensity.array[i]), 'Volume metalness and roughness must affect native ray-marched lighting.');
        assert(metallicDensity.array.every((v, i) => i % 4 !== 3 || v === volumeBeforeEmission.array[i]), 'Volume material changes must retain density alpha.');
        assert(renderer.volumeResourceStats.gridUploads === uploadsBeforeMaterial.gridUploads, 'Volume material updates must retain the uploaded density grid.');
        showPixels('Native metallic density', metallicDensity);
        await repr.createOrUpdate({ ignoreLight: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } }, volume).run(); native.tick(now());
        const volumeMode = plugin.canvas3dContext!.props.transparency, beforePeeling = await renderer.readPixels();
        try {
            plugin.canvas3dContext!.setProps({ transparency: 'dpoit' }); native.tick(now());
            const peeled = await renderer.readPixels();
            assert(peeled.array.every((v, i) => Math.abs(v - beforePeeling.array[i]) <= 2), 'A single ray-marched density layer must retain its color and alpha under native depth peeling.');
            assert((await renderer.pick(160, 120))?.id.objectId === repr.renderObjects[0].id, 'Depth-peeled densities must retain volume-cell picking.');
            await repr.createOrUpdate({ emissive: 0.5 }, volume).run(); native.tick(now());
            const glow = await renderer.readPixels();
            assert(glow.array.some((v, i) => i % 4 === 2 && v > 2 && peeled.array[i + 1] === 0), 'Depth-peeled volumes must preserve emissive bloom outside their original coverage.');
            const exported = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(exported.data.some((v, i) => i % 4 === 2 && v > 100 && exported.data[i + 1] > 2), 'Depth-peeled density screenshots must preserve emissive color and alpha.');
        } finally {
            await repr.createOrUpdate({ emissive: 0 }, volume).run(); plugin.canvas3dContext!.setProps({ transparency: volumeMode }); native.tick(now());
        }
        assert((await renderer.readPixels()).array.every((v, i) => v === beforePeeling.array[i]), 'Restoring weighted volume rendering must reproduce its original image.');
        const exportBeforeEmission = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
        await repr.createOrUpdate({ emissive: 0.5 }, volume).run(); native.tick(now());
        const volumeGlow = await renderer.readPixels();
        assert(volumeGlow.array.some((v, i) => i % 4 === 2 && v > 2 && volumeBeforeEmission.array[i + 1] === 0), 'Emissive density volumes must produce native bloom outside their original alpha coverage.');
        showPixels('Emissive density volume', volumeGlow);
        const exportGlow = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(exportGlow.data.some((v, i) => i % 4 === 2 && v > 100 && exportGlow.data[i + 1] > 2 && exportBeforeEmission.data[i + 1] === 0), 'Native image exports must include emissive density bloom.');
        showImage('Emissive volume export', exportGlow);
        await repr.createOrUpdate({ emissive: 0 }, volume).run(); native.tick(now());
        const clipObject = { ...Clip.Params.objects.ctor(), type: 'cube' as const, scale: Vec3.create(100, 100, 100) };
        await repr.createOrUpdate({ clip: { variant: 'pixel', objects: [clipObject] } }, volume).run(); native.tick(now());
        await verifyDebugCanvas(plugin, 'clipObjects');
        assert(!(await renderer.pick(160, 120)), 'Clipped density ray samples must not produce colors or picking IDs.');
        await repr.createOrUpdate({ clip: { variant: 'pixel', objects: [{ ...clipObject, invert: true }] } }, volume).run(); native.tick(now());
        assert(await renderer.pick(160, 120), 'Inverted volume clipping must retain density inside the clip shape.');
        await repr.createOrUpdate({ clip: { variant: 'instance', objects: [clipObject] } }, volume).run(); native.tick(now());
        assert(!(await renderer.pick(160, 120)), 'Density instance clipping must reject the full volume.');
        await repr.createOrUpdate({ clip: { variant: 'pixel', objects: [] } }, volume).run(); native.tick(now());
        native.setProps({ viewport: { name: 'static-frame', params: { x: 32, y: 24, width: 240, height: 160 } } });
        native.tick(now());
        assert(native.camera.viewport.width === 240 && native.camera.viewport.x === 32, 'Static viewport parameters must reach the camera.');
        assert((await renderer.pick(152, 136))?.id.objectId === repr.renderObjects[0].id, 'Offset viewport volume rays must stay centered in the active viewport.');
        native.setProps({ viewport: { name: 'canvas', params: {} } }); native.tick(now());
        const uploadsBeforeMarker = renderer.volumeResourceStats;
        native.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
        native.tick(now());
        const pixels = await renderer.readPixels();
        assert(pixels.array[(120 * pixels.width + 160) * 4] > 20, 'Volume highlighting must change native fragment colors.');
        showPixels('Volume highlighting', pixels);
        native.mark({ repr, loci: EveryLoci }, MarkerAction.RemoveHighlight);
        native.tick(now());
        assert(renderer.volumeResourceStats.gridUploads === uploadsBeforeMarker.gridUploads && renderer.volumeResourceStats.transferUploads === uploadsBeforeMarker.transferUploads && renderer.volumeResourceStats.paletteUploads === uploadsBeforeMarker.paletteUploads, 'Highlight updates must retain all uploaded volume textures.');
        const blockerMesh = Mesh.create(new Float32Array([-8, -8, 5, 8, -8, 5, 0, 8, 5]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
        const blockerProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
        const blocker = createRenderObject('mesh', Mesh.Utils.createValuesSimple(blockerMesh, blockerProps, Color(0xff0000), 1), Mesh.Utils.createRenderableState(blockerProps), -1);
        renderer.render([...repr.renderObjects, blocker], native.camera, native.props.renderer, true);
        const occluded = await renderer.readPixels();
        const center = (120 * occluded.width + 160) * 4;
        assert(occluded.array[center] > 240 && occluded.array[center + 2] < 10, 'Opaque geometry in front of a density volume must occlude its ray march.');
        assert((await renderer.pick(160, 120))?.id.objectId === blocker.id, 'Occluded volume fragments must not replace foreground picking IDs.');

        const duplicated: Volume = { ...volume, instances: [
            { transform: Mat4.fromTranslation(Mat4(), Vec3.create(-4, 0, 0)) },
            { transform: Mat4.fromTranslation(Mat4(), Vec3.create(4, 0, 0)) },
        ] };
        await repr.createOrUpdate({}, duplicated).run();
        native.update(repr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        for (const [instance, x] of [-4, 4].entries()) {
            const screen = native.camera.project(Vec4(), Vec3.create(x, 0, 0));
            const picked = await renderer.pick(screen[0], 240 - screen[1]);
            assert(picked?.id.instanceId === instance, 'Native volume instancing must retain logical instance IDs.');
        }
        showPixels('Instanced density volumes', await renderer.readPixels());
        const sliceResults = await verifySlices(plugin, volume);
        const surfaceResults = await verifyIsosurface(plugin, volume);
        const segmentResults = await verifySegments(plugin);
        const paletteResults = await verifyVolumePalettes(plugin, volume);
        const vertexThemeResults = await verifyVolumeVertexThemes(plugin, volume);
        const orbitalResults = await verifyOrbitalVisualization(plugin);
        return ['native ray-based volume picking, camera scaling, cached identify, concurrent requests and unchanged displayed frame', 'volume environment composition, fogged ray opacity and independent cell picking', ...sliceResults, ...surfaceResults, ...segmentResults, ...paletteResults, ...vertexThemeResults, ...orbitalResults, 'direct volume ray marching in byte, float and half-float formats', 'volume GGX materials and gradient lighting', 'volume opacity-threshold selection with independent color rendering', 'volume-cell picking and axis order', 'volume highlighting, emissive bloom/export and native clipping', 'volume depth occlusion and offset viewport', 'instanced density volumes'];
    } finally {
        native.remove(repr); repr.destroy(); volume.customProperties.dispose();
    }
}


async function verifyNativeRayPicking(plugin: PluginContext) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const originalCamera = native.camera.getSnapshot();
    const read = async (ray: Ray3D) => {
        const pending = native.asyncIdentify(ray);
        assert(pending, 'Native ray picking must expose asynchronous readback.');
        const deadline = performance.now() + 10000;
        for (;;) {
            const result = pending.tryGet();
            if (result !== 'pending') return result;
            assert(performance.now() < deadline, 'Native ray picking must finish.');
            await new Promise(resolve => setTimeout(resolve, 0));
        }
    };
    try {
        for (const scale of [1, 0.25, 2]) {
            native.camera.scale = scale; native.requestDraw(); native.tick(now());
            const before = await renderer.readPixels();
            const camera = native.camera.getSnapshot();
            const ray = Ray3D.targetTo(Ray3D(), Ray3D.create(Vec3.scale(Vec3(), native.camera.position, scale), Vec3()), Vec3.scale(Vec3(), native.camera.target, scale));
            const hit = await read(ray);
            assert(hit && Volume.Cell.isLoci(native.getLoci(hit.id).loci), 'Native rays must resolve to the displayed volume cells at every camera scale.');
            const scaledHit = Vec3.scale(Vec3(), hit.position, scale);
            const offset = Vec3.sub(Vec3(), scaledHit, ray.origin);
            const along = Vec3.dot(offset, ray.direction);
            const perpendicular = Vec3.sub(Vec3(), offset, Vec3.scale(Vec3(), ray.direction, along));
            assert(along > 0 && Vec3.magnitude(perpendicular) < 0.0001, 'Native ray picks must return a molecular position on the forward ray.');
            assert(native.identify(ray)?.id.groupId === hit.id.groupId, 'Synchronous ray identify must return its completed cached query.');
            const miss = Ray3D.create(Vec3.add(Vec3(), ray.origin, Vec3.create(100 * scale, 0, 0)), Vec3.clone(ray.direction));
            assert(await read(miss) === undefined, 'Rays missing native geometry must return no hit.');
            assert(await read(Ray3D.create(Vec3(), Vec3())) === undefined, 'Zero-direction rays must return no hit without invalid GPU commands.');
            const after = await renderer.readPixels();
            assert(after.array.every((v, i) => v === before.array[i]), 'Ray queries must leave the displayed color frame intact.');
            assert(JSON.stringify(native.camera.getSnapshot()) === JSON.stringify(camera), 'Ray queries must preserve the main camera.');
            const results = await Promise.all([read(ray), read(miss), read(ray)]);
            assert(results[0] && !results[1] && results[2] && results[0].id.groupId === results[2].id.groupId, 'Concurrent native ray queries must retain independent results.');
        }
    } finally { native.camera.scale = 1; native.camera.setState(originalCamera, 0); native.requestDraw(); native.tick(now()); }
}

async function verifyVolumeVertexThemes(plugin: PluginContext, source: Volume) {
    const volume: Volume = { ...source, grid: { ...source.grid, cells: Tensor.create(source.grid.cells.space, Tensor.Data1(new Float32Array(source.grid.cells.data.length).fill(0.5))), stats: { min: 0, max: 1, mean: 0.5, sigma: 0.1 } }, _propertyData: {}, _localPropertyData: {} };
    const provider = plugin.representation.volume.registry.get('direct-volume');
    const params = createVolumeRepresentationParams(plugin, volume, { type: 'direct-volume', typeParams: { ignoreLight: true, dataType: 'float' }, color: 'uniform', colorParams: { value: Color(0xffffff) } });
    const repr = provider.factory(plugin.representation.volume.themes, provider.getParams);
    repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, params));
    const renderer = await WebGPURenderer.create(plugin.canvas3d!.webgpu!);
    const camera = new Camera({ mode: 'orthographic', position: Vec3.create(0, 0, 20), target: Vec3(), radius: 8, radiusMax: 8, fog: 0 }, { x: 0, y: 0, width: 320, height: 240 }); camera.update();
    const props = PD.getDefaultValues(RendererParams);
    const draw = async () => {
        renderer.render(repr.renderObjects, camera, props, true, 1, { width: 320, height: 240, present: false });
        return renderer.readPixels();
    };
    try {
        await repr.createOrUpdate(params.type.params, volume).run();
        const values = repr.renderObjects[0].values;
        const dimensions = value<number[]>(values, 'uGridDim', []), [nx, ny, nz] = dimensions;
        const count = nx * ny * nz;
        const colors = new Uint8Array(count * 3), overlays = new Uint8Array(count * 4);
        for (let x = 0; x < nx; x++) for (let y = 0; y < ny; y++) for (let z = 0; z < nz; z++) {
            const index = z + y * nz + x * ny * nz;
            const color = [Math.round(x / (nx - 1) * 255), Math.round(y / (ny - 1) * 255), 64];
            colors.set(color, index * 3); overlays.set([...color, 255], index * 4);
        }
        const cartnToUnit = Mat4.invert(Mat4(), value<Mat4>(values, 'uUnitToCartn', Mat4.identity()));
        const referenceByte = (position: number, dimension: number) => {
            const clamped = Math.max(0, Math.min(dimension - 1, position));
            const lo = Math.floor(clamped), hi = Math.min(lo + 1, dimension - 1), fraction = clamped - lo;
            return Math.round(lo / (dimension - 1) * 255) * (1 - fraction) + Math.round(hi / (dimension - 1) * 255) * fraction;
        };
        for (const mode of ['blended', 'wboit', 'dpoit'] as const) {
            renderer.setTransparency(mode);
            ValueCell.update(values.dColorType, 'uniform'); ValueCell.update(values.dOverpaint, false);
            const baseline = await draw();
            ValueCell.update(values.tColor, { array: colors, width: count, height: 1 });
            ValueCell.update(values.dColorType, 'vertex');
            const vertex = await draw();
            let checked = 0;
            for (let y = 35; y < 205; y += 19) for (let x = 45; x < 275; x += 17) {
                const i = (y * 320 + x) * 4, alpha = vertex.array[i + 3] / 255;
                if (alpha < 0.5) continue;
                const world = camera.unproject(Vec3(), Vec3.create(x + 0.5, 240 - y - 0.5, 0.5));
                const unit = Vec3.transformMat4(Vec3(), world, cartnToUnit);
                const expected = [referenceByte(unit[0] * nx, nx), referenceByte(unit[1] * ny, ny), 64];
                for (let c = 0; c < 3; c++) assert(Math.abs(vertex.array[i + c] - expected[c] * alpha) <= 2, `${mode} density vertex colors must match the continuous analytical ramp, channel=${c}, actual=${vertex.array[i + c]}, expected=${expected[c] * alpha}.`);
                checked++;
            }
            assert(checked > 10, 'Vertex-color verification must cover multiple visible rays.');
            assert(vertex.array.every((v, i) => i % 4 !== 3 || v === baseline.array[i]), 'Density vertex colors must retain original opacity.');
            ValueCell.update(values.dColorType, 'vertexInstance');
            const vertexInstance = await draw();
            assert(vertexInstance.array.every((v, i) => v === vertex.array[i]), 'Single-instance vertex colors must agree across both granularities.');
            ValueCell.update(values.dColorType, 'uniform'); ValueCell.update(values.dOverpaint, true);
            ValueCell.update(values.dOverpaintType, 'vertexInstance'); ValueCell.update(values.tOverpaint, { array: overlays, width: count, height: 1 });
            const overpaint = await draw();
            assert(overpaint.array.every((v, i) => Math.abs(v - vertex.array[i]) <= 1), 'Vertex overpaint must interpolate continuously and agree with the independent vertex-color path.');
            const pass = new WebGPUImagePass(plugin.canvas3d!.webgpu!, camera, () => [...repr.renderObjects], { renderer: props, transparentBackground: true, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' }, postprocessing: { ...PD.getDefaultValues(PostprocessingParams), enabled: false }, cameraHelper: { axes: { name: 'off', params: {} } } }, undefined, { transparency: () => mode });
            try {
                const exported = await pass.getImageData(RuntimeContext.Synchronous, 320, 240);
                assert(exported.data.every((v, i) => {
                    const alpha = overpaint.array[i - i % 4 + 3];
                    const expected = i % 4 === 3 ? alpha : alpha ? Math.min(255, Math.round(overpaint.array[i] * 255 / alpha)) : 0;
                    return Math.abs(v - expected) <= 2;
                }), 'Density vertex overpaint exports must agree with canvas colors and transparent alpha.');
            } finally { await pass.dispose(); }
        }
        ValueCell.update(values.dOverpaint, false); ValueCell.update(values.dColorType, 'vertexInstance');
        const instancedColors = new Uint8Array(count * 6);
        for (let i = 0; i < count; i++) { instancedColors.set([255, 0, 0], i * 3); instancedColors.set([0, 255, 0], (count + i) * 3); }
        const duplicate: Volume = { ...volume, instances: [{ transform: Mat4.fromTranslation(Mat4(), Vec3.create(-5, 0, 0)) }, { transform: Mat4.fromTranslation(Mat4(), Vec3.create(5, 0, 0)) }] };
        await repr.createOrUpdate({}, duplicate).run();
        const instanceValues = repr.renderObjects[0].values;
        ValueCell.update(instanceValues.dColorType, 'vertexInstance'); ValueCell.update(instanceValues.tColor, { array: instancedColors, width: count * 2, height: 1 });
        ValueCell.update(instanceValues.dOverpaint, false);
        const frame = await draw();
        for (const [instance, x] of [-5, 5].entries()) {
            const screen = camera.project(Vec4(), Vec3.create(x, 0, 0));
            const pick = await renderer.pick(screen[0], 240 - screen[1]);
            assert(pick?.id.instanceId === instance, 'Vertex-colored densities must preserve per-instance picking.');
            const i = (Math.floor(240 - screen[1]) * 320 + Math.floor(screen[0])) * 4;
            assert(frame.array[i + instance] > 128 && frame.array[i + 1 - instance] === 0, 'Density vertex-instance colors must use the matching instance stride.');
        }
        const instanceOverlays = new Uint8Array(count * 8);
        for (let i = 0; i < count; i++) { instanceOverlays.set([255, 0, 0, 255], i * 4); instanceOverlays.set([0, 255, 0, 255], (count + i) * 4); }
        ValueCell.update(instanceValues.dColorType, 'uniform'); ValueCell.update(instanceValues.uColor, Vec3.create(1, 1, 1));
        ValueCell.update(instanceValues.dOverpaint, true); ValueCell.update(instanceValues.dOverpaintType, 'vertexInstance');
        ValueCell.update(instanceValues.tOverpaint, { array: instanceOverlays, width: count * 2, height: 1 });
        const instanceOverpaint = await draw();
        assert(instanceOverpaint.array.every((v, i) => v === frame.array[i]), 'Density vertex overpaint must use the matching per-instance stride.');
        ValueCell.update(instanceValues.aInstance, new Float32Array([1, 0]));
        const reordered = await draw();
        for (const [physical, x] of [-5, 5].entries()) {
            const logical = 1 - physical, screen = camera.project(Vec4(), Vec3.create(x, 0, 0));
            assert((await renderer.pick(screen[0], 240 - screen[1]))?.id.instanceId === logical, 'Reordered density instances must retain logical picking IDs.');
            const i = (Math.floor(240 - screen[1]) * 320 + Math.floor(screen[0])) * 4;
            assert(reordered.array[i + logical] > 128 && reordered.array[i + 1 - logical] === 0, 'Density vertex overpaint must follow logical instance IDs after reordering.');
        }
        return ['density vertex and vertex-instance colors, analytical trilinear ramps, continuous overpaint, all transparency modes, screenshot parity and reordered instance picking'];
    } finally { renderer.dispose(); repr.destroy(); }
}

async function verifyVolumePalettes(plugin: PluginContext, source: Volume) {
    const volume: Volume = { ...source, grid: { ...source.grid, cells: Tensor.create(source.grid.cells.space, Tensor.Data1(new Float32Array(source.grid.cells.data.length).fill(0.5))), stats: { min: 0, max: 1, mean: 0.5, sigma: 0.1 } }, _propertyData: {}, _localPropertyData: {} };
    const provider = plugin.representation.volume.registry.get('direct-volume');
    const params = createVolumeRepresentationParams(plugin, volume, { type: 'direct-volume', typeParams: { ignoreLight: true, dataType: 'float' }, color: 'volume-value' });
    const repr = provider.factory(plugin.representation.volume.themes, provider.getParams);
    repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, params));
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.clear();
    try {
        await repr.createOrUpdate(params.type.params, volume).run(); native.add(repr);
        native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        const values = repr.renderObjects[0].values;
        ValueCell.update(values.uPaletteDomain, Vec2.create(0, 1));
        for (const filter of ['nearest', 'linear'] as const) {
            const uploadsBeforePalette = renderer.volumeResourceStats;
            ValueCell.update(values.dColorType, 'direct'); ValueCell.update(values.dUsePalette, true);
            ValueCell.update(values.tPalette, { array: new Uint8Array([255, 0, 0, 0, 255, 0, 0, 0, 255, 255, 255, 0]), width: 4, height: 1, filter });
            native.requestDraw(); native.tick(now()); const paletteFrame = await renderer.readPixels();
            assert(renderer.volumeResourceStats.gridUploads === uploadsBeforePalette.gridUploads && renderer.volumeResourceStats.transferUploads === uploadsBeforePalette.transferUploads && renderer.volumeResourceStats.paletteUploads === uploadsBeforePalette.paletteUploads + 1, 'Palette updates must retain density and transfer textures while uploading the changed palette once.');
            const picked = await renderer.pick(160, 120);
            assert(picked && Volume.Cell.isLoci(native.getLoci(picked.id).loci), 'Direct volume palette changes must preserve cell picking.');
            const paletteExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            ValueCell.update(values.dColorType, 'uniform'); ValueCell.update(values.dUsePalette, false);
            ValueCell.update(values.uColor, filter === 'nearest' ? Vec3.create(0, 0, 1) : Vec3.create(0, 0.5, 0.5));
            native.requestDraw(); native.tick(now()); const uniformFrame = await renderer.readPixels();
            assert(paletteFrame.array.some((v, i) => i % 4 === 3 && v > 0), 'Constant-density palette volumes must render visible pixels.');
            assert(paletteFrame.array.every((v, i) => Math.abs(v - uniformFrame.array[i]) <= 1), 'Constant-density volume palettes must match the analytical nearest/linear uniform-color reference.');
            const uniformExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(paletteExport.data.every((v, i) => Math.abs(v - uniformExport.data[i]) <= 1), 'Volume palette filtering must also match the uniform reference in screenshot exports.');
        }
        return ['direct volume nearest/linear palettes, constant-density analytical references, cell picking and screenshot parity'];
    } finally { native.remove(repr); repr.destroy(); }
}

async function verifySlices(plugin: PluginContext, volume: Volume) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.clear();
    const provider = plugin.representation.volume.registry.get('slice');
    const params = createVolumeRepresentationParams(plugin, volume, {
        type: 'slice', typeParams: { dimension: { name: 'z', params: 16 }, isoValue: Volume.IsoValue.absolute(0.1) },
        color: 'uniform', colorParams: { value: Color(0x0088ff) },
    });
    const repr = provider.factory(plugin.representation.volume.themes, provider.getParams);
    repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, params));
    try {
        for (const interpolation of ['nearest', 'catmulrom', 'mitchell', 'bspline'] as const) {
            await repr.createOrUpdate({ ...params.type.params, interpolation }, volume).run();
            native.add(repr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
            const pixels = await renderer.readPixels(), o = (120 * pixels.width + 160) * 4;
            assert(pixels.array[o + 2] > 200 && pixels.array[o] < 10, `Native ${interpolation} slice must display precomputed density colors.`);
            const pick = await renderer.pick(160, 120);
            assert(pick && Volume.Cell.isLoci(native.getLoci(pick.id).loci), 'Slice picking must resolve to volume-cell loci.');
            // Slice groups are mapped by the representation to original tensor cells.
            if (interpolation === 'bspline') showPixels('Native interpolated density slice', pixels);
        }
        await verifyDebugCanvas(plugin, 'imageEdges');
        native.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight); native.tick(now());
        const highlighted = await renderer.readPixels();
        assert(highlighted.array[(120 * highlighted.width + 160) * 4] > 20, 'Slice markers must be sampled per cell.');
        native.mark({ repr, loci: EveryLoci }, MarkerAction.RemoveHighlight);
        await repr.createOrUpdate({ isoValue: Volume.IsoValue.absolute(1.1) }, volume).run();
        native.tick(now());
        assert(!(await renderer.pick(180, 120)), 'Iso-value masking must remove low-density slice pixels and picking.');
        await repr.createOrUpdate({ mode: 'plane', plane: { point: Vec3(), normal: Vec3.create(0.2, 0.1, 1) }, isoValue: Volume.IsoValue.absolute(0.1) }, volume).run();
        native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        assert(await renderer.pick(160, 120), 'Arbitrary plane slices must retain visible cells after trimming in grid coordinates.');
        showPixels('Arbitrary density plane', await renderer.readPixels());
        const paletteParams = createVolumeRepresentationParams(plugin, volume, { type: 'slice', color: 'volume-value', colorParams: { colorList: { kind: 'interpolate', colors: [Color(0xff0000), Color(0x00ff00)] } } });
        repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, paletteParams));
        await repr.createOrUpdate({ mode: 'grid', isoValue: Volume.IsoValue.absolute(0.1) }, volume).run();
        native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        const palettePixels = await renderer.readPixels(), center = (120 * palettePixels.width + 160) * 4;
        assert(palettePixels.array[center + 1] > 180 && palettePixels.array[center] < 80, 'Slice packed colors must decode through the palette in the fragment shader.');
        showPixels('Palette-colored density slice', palettePixels);
        const image = repr.renderObjects[0], post = native.props.postprocessing;
        const originalPalette = value(image.values, 'tPalette', { array: new Uint8Array(3), width: 1, height: 1 });
        const originalDefault = Vec3.clone(value(image.values, 'uPaletteDefault', Vec3()));
        native.setProps({ postprocessing: { ...post, antialiasing: { name: 'off', params: {} } } });
        ValueCell.update(image.values.uPaletteDefault, Vec3.create(1, 0, 1));
        const grayscalePixels = (array: Uint8Array) => {
            const values: number[] = [];
            for (let i = 0; i < array.length; i += 4) if (array[i + 3] === 255 && array[i] === array[i + 1] && array[i] === array[i + 2]) values.push(array[i]);
            return values;
        };
        ValueCell.update(image.values.tPalette, { array: new Uint8Array([0, 0, 0, 255, 255, 255]), width: 2, height: 1, filter: 'nearest' });
        native.requestDraw(); native.tick(now()); const nearestPalette = await renderer.readPixels();
        const nearestGray = grayscalePixels(nearestPalette.array);
        assert(nearestGray.length > 10 && nearestGray.every(v => v === 0 || v === 255), 'Nearest image palettes must produce discrete entries without blending at opaque interior pixels.');
        ValueCell.update(image.values.tPalette, { array: new Uint8Array([0, 0, 0, 255, 255, 255]), width: 2, height: 1, filter: 'linear' });
        native.requestDraw(); native.tick(now()); const linearPalette = await renderer.readPixels();
        assert(grayscalePixels(linearPalette.array).some(v => v > 10 && v < 245), 'Linear image palettes must interpolate between adjacent palette entries.');
        assert(linearPalette.array.some((v, i) => v !== nearestPalette.array[i]), 'Changing palette filtering must update the rendered slice.');
        const palettePick = await renderer.pick(160, 120);
        assert(palettePick && Volume.Cell.isLoci(native.getLoci(palettePick.id).loci), 'Palette filter changes must preserve slice-cell picking.');
        const paletteExport = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(paletteExport.data.some((v, i) => i % 4 === 3 && v > 0), 'Filtered image palettes must render in screenshot exports.');
        ValueCell.update(image.values.tPalette, originalPalette); ValueCell.update(image.values.uPaletteDefault, originalDefault);
        native.setProps({ postprocessing: post }); native.requestDraw(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === palettePixels.array[i]), 'Restoring image palette data and filtering must restore the original slice frame.');
        return ['slice palette decoding and nearest/linear filtering, screenshot export and restoration', 'arbitrary plane slices and grid-space trimming', 'volume slices with four interpolation modes', 'slice-cell picking and highlighting', 'slice iso-value masking'];
    } finally { native.remove(repr); repr.destroy(); }
}


async function verifySegments(plugin: PluginContext) {
    const space = Tensor.Space([16, 16, 16], [1, 2, 0], Float32Array), data = space.create();
    for (let x = 4; x < 12; x++) for (let y = 4; y < 12; y++) for (let z = 4; z < 12; z++) space.set(data, x, y, z, 1);
    const matrix = Mat4.fromScaling(Mat4(), Vec3.create(0.5, 0.5, 0.5)); Mat4.setTranslation(matrix, Vec3.create(-4, -4, -4));
    const volume: Volume = {
        grid: { cells: Tensor.create(space, data), transform: { kind: 'matrix', matrix }, stats: { min: 0, max: 1, mean: 0.125, sigma: 0.3 } },
        instances: [{ transform: Mat4.identity() }], sourceData: { kind: 'synthetic', name: 'Segment verification labels', data: undefined },
        customProperties: new CustomProperties(), _propertyData: {}, _localPropertyData: {},
    };
    const segment = 3 as Volume.SegmentIndex;
    Volume.Segmentation.set(volume, { segments: new Map([[segment, new Set([1])]]), sets: new Map([[1, new Set([segment])]]), bounds: { [segment]: Box3D.create(Vec3.create(4, 4, 4), Vec3.create(12, 12, 12)) }, labels: { [segment]: 'Cube' } });
    Volume.PickingGranularity.set(volume, 'object');
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.clear();
    const provider = plugin.representation.volume.registry.get('segment');
    const params = createVolumeRepresentationParams(plugin, volume, { type: 'segment', typeParams: { segments: [3], ignoreLight: true }, color: 'uniform', colorParams: { value: Color(0x00cc88) } });
    const repr = provider.factory({ ...plugin.representation.volume.themes, webgpu: native.webgpu }, provider.getParams);
    repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, params));
    try {
        const before = native.webgpu!.stats.marchingCubesDispatches;
        await repr.createOrUpdate(params.type.params, volume).run(); native.add(repr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        assert(native.webgpu!.stats.marchingCubesDispatches > before, 'Segment surfaces must use native GPU extraction.');
        const pixels = await renderer.readPixels(), picked = await renderer.pick(160, 120);
        assert(pixels.array[(120 * pixels.width + 160) * 4 + 1] > 150, 'Cropped segment surfaces must render with segment colors.');
        assert(picked && Volume.Segment.isLoci(native.getLoci(picked.id).loci), 'Segment picking must resolve the source segment.');
        Volume.PickingGranularity.set(volume, 'voxel');
        assert(picked && Volume.Cell.isLoci(native.getLoci(picked.id).loci), 'Segment voxel picking must resolve original grid cells.');
        const coordinate = [0, 0, 0]; space.getCoords(picked.id.groupId, coordinate);
        assert(coordinate.every(v => v >= 3 && v <= 12), 'Segment voxel IDs must use global tensor coordinates rather than cropped-grid offsets.');
        const dispatches = native.webgpu!.stats.marchingCubesDispatches, count = repr.renderObjects[0].values.drawCount.ref.value;
        await repr.createOrUpdate({ tryUseGpu: false }, volume).run(); native.tick(now());
        assert(native.webgpu!.stats.marchingCubesDispatches === dispatches && repr.renderObjects[0].values.drawCount.ref.value === count, 'Segment CPU fallback must preserve topology and honor the GPU setting.');
        const cpuPick = await renderer.pick(160, 120);
        assert(cpuPick && Volume.Cell.isLoci(native.getLoci(cpuPick.id).loci), 'CPU segment fallback must preserve source-grid voxel loci.');
        space.getCoords(cpuPick.id.groupId, coordinate); assert(coordinate.every(v => v >= 3 && v <= 12), 'CPU segment voxel IDs must also retain global coordinates.');
        await repr.createOrUpdate({ tryUseGpu: true }, volume).run(); native.tick(now());
        assert(native.webgpu!.stats.marchingCubesDispatches > dispatches && (await renderer.readPixels()).array.every((v, i) => v === pixels.array[i]), 'Restoring GPU segment extraction must reproduce the original frame.');
        const beforeSelection = native.webgpu!.stats.marchingCubesDispatches;
        await repr.createOrUpdate({ segments: [] }, volume).run(); native.tick(now());
        assert(repr.renderObjects.length === 0 && !(await renderer.pick(160, 120)), 'Removing the selected segments must clear geometry and picking.');
        await repr.createOrUpdate({ segments: [3] }, volume).run(); native.tick(now());
        assert(native.webgpu!.stats.marchingCubesDispatches > beforeSelection && (await renderer.readPixels()).array.every((v, i) => v === pixels.array[i]), 'Restoring segment selection must rebuild GPU geometry and preserve the view.');
        const edgeData = space.create(); edgeData.fill(1);
        const edgeVolume: Volume = { ...volume, grid: { ...volume.grid, cells: Tensor.create(space, edgeData) }, _propertyData: {}, _localPropertyData: {} };
        Volume.Segmentation.set(edgeVolume, { segments: new Map([[segment, new Set([1])]]), sets: new Map([[1, new Set([segment])]]), bounds: { [segment]: Box3D.create(Vec3(), Vec3.create(16, 16, 16)) }, labels: { [segment]: 'Boundary cube' } });
        Volume.PickingGranularity.set(edgeVolume, 'voxel');
        await repr.createOrUpdate({}, edgeVolume).run(); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        const edgePixels = await renderer.readPixels(), edgePick = await renderer.pick(160, 120);
        assert(edgePick && Volume.Cell.isLoci(native.getLoci(edgePick.id).loci), 'Boundary segments must remain pickable through the original volume.');
        assert(value(repr.renderObjects[0].values, 'aGroup', new Float32Array()).every(v => Number.isInteger(v) && v >= 0 && v < edgeData.length), 'Padding at every source-grid face must preserve valid original-grid voxel IDs.');
        const vertices = value(repr.renderObjects[0].values, 'aPosition', new Float32Array());
        for (let a = 0; a < 3; a++) {
            let min = Infinity, max = -Infinity;
            for (let i = a; i < vertices.length; i += 3) { min = Math.min(min, vertices[i]); max = Math.max(max, vertices[i]); }
            assert(min > -4.5 && min < -4 && max > 3.5 && max < 4, 'Boundary segment geometry must close within half a source voxel on every face.');
        }
        const edgeCount = repr.renderObjects[0].values.drawCount.ref.value;
        await repr.createOrUpdate({ tryUseGpu: false }, edgeVolume).run(); native.tick(now());
        assert(repr.renderObjects[0].values.drawCount.ref.value === edgeCount, 'Boundary CPU/GPU segment surfaces must preserve topology.');
        await repr.createOrUpdate({ tryUseGpu: true }, edgeVolume).run(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === edgePixels.array[i]), 'Boundary segment GPU restoration must reproduce its original frame.');
        showPixels('Native segment at all volume boundaries', edgePixels);
        showPixels('Native volume segment', pixels);
        return ['native cropped segment extraction, segment colors, segment/voxel loci, selection controls, all grid boundaries and CPU/GPU switching'];
    } finally { native.remove(repr); repr.destroy(); volume.customProperties.dispose(); }
}

async function verifyIsosurface(plugin: PluginContext, volume: Volume) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.clear();
    const provider = plugin.representation.volume.registry.get('isosurface');
    const params = createVolumeRepresentationParams(plugin, volume, {
        type: 'isosurface', typeParams: { isoValue: Volume.IsoValue.absolute(0.3), ignoreLight: true },
        color: 'uniform', colorParams: { value: Color(0xff8800) },
    });
    const repr = provider.factory({ ...plugin.representation.volume.themes, webgpu: native.webgpu }, provider.getParams);
    repr.setTheme(Theme.create(plugin.representation.volume.themes, { volume, locationKinds: provider.locationKinds }, params));
    try {
        const computeBefore = native.webgpu!.stats.computeDispatches;
        await repr.createOrUpdate(params.type.params, volume).run();
        assert(repr.renderObjects.some(o => o.type === 'mesh'), 'GPU-extracted isosurfaces must expose native mesh geometry.');
        native.add(repr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        const pixels = await renderer.readPixels(), o = (120 * pixels.width + 160) * 4;
        assert(native.webgpu!.stats.computeDispatches > computeBefore, 'Native isosurface representations must extract triangles on the GPU.');
        assert(pixels.array[o] > 240 && pixels.array[o + 1] > 100, 'GPU-generated isosurfaces must render on WebGPU.');
        const pick = await renderer.pick(160, 120);
        assert(pick && !isEmptyLoci(native.getLoci(pick.id).loci), 'Isosurface picking must resolve to volume loci.');
        showPixels('Native density isosurface', pixels);
        const count = repr.renderObjects[0].values.drawCount.ref.value;
        await repr.createOrUpdate({ isoValue: Volume.IsoValue.absolute(0.7) }, volume).run();
        native.tick(now());
        assert(repr.renderObjects[0].values.drawCount.ref.value !== count, 'Changing iso-values must regenerate surface geometry.');
        const gpuSurface = await renderer.readPixels(), dispatches = native.webgpu!.stats.computeDispatches;
        await repr.createOrUpdate({ tryUseGpu: false }, volume).run(); native.tick(now());
        assert(native.webgpu!.stats.computeDispatches === dispatches, 'Disabling isosurface GPU extraction must honor the CPU setting.');
        assert((await renderer.readPixels()).array.some((v, i) => i % 4 === 3 && v > 0), 'CPU extraction must remain available for native rendering.');
        await repr.createOrUpdate({ tryUseGpu: true }, volume).run(); native.tick(now());
        assert(native.webgpu!.stats.computeDispatches > dispatches, 'Restoring GPU extraction must regenerate the volume surface.');
        assert((await renderer.readPixels()).array.every((v, i) => v === gpuSurface.array[i]), 'Restored GPU extraction must reproduce the original isosurface frame.');
        const beforeWrap = native.webgpu!.stats.marchingCubesDispatches;
        await repr.createOrUpdate({ wrap: 'on' }, volume).run(); native.tick(now());
        assert(native.webgpu!.stats.marchingCubesDispatches > beforeWrap, 'Wrapped volume surfaces must use native GPU extraction.');
        const wrappedObjects = native.getRenderObjects().filter(object => object.type === 'mesh');
        assert(wrappedObjects.every(object => value(object.values, 'aGroup', new Float32Array()).every(v => v >= 0 && v < volume.grid.cells.data.length)), 'Wrapped surface groups must stay within the original volume cell range.');
        assert((await renderer.pick(160, 120)) && (await renderer.readPixels()).array.some((v, i) => i % 4 === 3 && v > 0), 'Wrapped volume surfaces must remain visible and pickable.');
        await repr.createOrUpdate({ wrap: 'off' }, volume).run(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === gpuSurface.array[i]), 'Disabling periodic wrapping must restore the original native volume surface.');
        for (const floodfill of ['inside', 'outside'] as const) {
            const beforeFill = native.webgpu!.stats.marchingCubesDispatches;
            await repr.createOrUpdate({ floodfill }, volume).run(); native.tick(now());
            assert(native.webgpu!.stats.marchingCubesDispatches > beforeFill, 'Floodfilled volumes must retain native GPU surface extraction.');
            const filledCount = repr.renderObjects[0].values.drawCount.ref.value;
            await repr.createOrUpdate({ tryUseGpu: false }, volume).run(); native.tick(now());
            assert(repr.renderObjects[0].values.drawCount.ref.value === filledCount, 'Floodfilled CPU/GPU surfaces must preserve topology.');
            assert(value(repr.renderObjects[0].values, 'aGroup', new Float32Array()).every(v => Number.isInteger(v) && v >= 0 && v < volume.grid.cells.data.length), 'Floodfill must not alter original volume cell IDs in the CPU fallback.');
            await repr.createOrUpdate({ tryUseGpu: true }, volume).run(); native.tick(now());
        }
        await repr.createOrUpdate({ floodfill: 'off' }, volume).run(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === gpuSurface.array[i]), 'Disabling floodfill must restore the original native volume frame.');
        return ['native GPU marching cubes with volume isosurface rendering', 'isosurface picking, iso-value updates, periodic wrapping and floodfill controls'];
    } finally { native.remove(repr); repr.destroy(); }
}

if (!new URLSearchParams(window.location.search).has('interaction-only')) {
    (window as unknown as { webgpuVerification: ReturnType<typeof verify> }).webgpuVerification = verify();
}

interface WebGPUInteractionVerification {
    snapshot(): Promise<{ position: number[], target: number[], up: number[], distance: number, frameHash: number, clicks: number, hover: boolean, selectedAtoms: number, errors: string[] }>
    target(): Promise<{ x: number, y: number }>
    emptyTarget(): Promise<{ x: number, y: number }>
    setSelectionMode(enabled: boolean): void
    dispose(): void
}

async function prepareInteractionVerification() {
    const container = document.createElement('div');
    container.style.cssText = 'position:fixed;left:12px;top:12px;width:320px;height:240px;z-index:1000';
    const canvas = document.createElement('canvas'); canvas.width = 320; canvas.height = 240;
    canvas.id = 'webgpu-interaction-canvas'; canvas.style.cssText = 'display:block;width:320px;height:240px';
    container.appendChild(canvas); document.body.appendChild(container);
    const plugin = new PluginContext(DefaultPluginSpec());
    assert(plugin.config.get(PluginConfig.General.RenderingBackend) === 'webgpu', 'New plugins must select WebGPU by default.');
    try {
        await plugin.init();
        assert(await plugin.initViewerAsync(canvas, container), 'Interactive WebGPU viewer must initialize.');
        plugin.animationLoop.stop({ noDraw: true });
        const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
        const errors: string[] = []; native.webgpu!.errors.subscribe(error => errors.push(error.message));
        const pdb = [
            'ATOM      1  N   GLY A   1      -3.000   0.000   0.000  1.00 20.00           N  ',
            'ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  ',
            'ATOM      3  C   GLY A   1       3.000   0.000   0.000  1.00 20.00           C  ',
            'ATOM      4  O   GLY A   1       3.000   3.000   1.000  1.00 20.00           O  ',
            'END',
        ].join('\n');
        const data = await plugin.builders.data.rawData({ data: pdb, label: 'Interactive WebGPU molecule' });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill', color: 'element-symbol' });
        plugin.managers.interactivity.setProps({ granularity: 'element' });
        native.setProps({ cameraFog: { name: 'off', params: {} }, trackball: { staticMoving: true } });
        native.resume(); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        let clicks = 0, hover = false;
        const clickSub = native.interaction.click.subscribe(e => { if (!isEmptyLoci(e.current.loci)) clicks++; });
        const hoverSub = native.interaction.hover.subscribe(e => { hover = !isEmptyLoci(e.current.loci); });
        native.animate();
        const api: WebGPUInteractionVerification = {
            async snapshot() {
                // Input and marking can request a redraw during the current frame's interaction tick.
                await new Promise<void>(resolve => window.requestAnimationFrame(() => window.requestAnimationFrame(() => resolve())));
                for (let i = 0; i < 8 && renderer.multiSampleNeedsFrame; i++) await new Promise<void>(resolve => window.requestAnimationFrame(() => resolve()));
                const pixels = await renderer.readPixels();
                let hash = 2166136261;
                for (const byte of pixels.array) hash = Math.imul(hash ^ byte, 16777619);
                return { position: [...native.camera.state.position], target: [...native.camera.state.target], up: [...native.camera.state.up], distance: Vec3.distance(native.camera.state.position, native.camera.state.target), frameHash: hash >>> 0, clicks, hover, selectedAtoms: plugin.managers.structure.selection.stats.elementCount, errors: [...errors] };
            },
            async target() {
                const projected = native.camera.project(Vec4(), Vec3.create(-3, 0, 0));
                const x = projected[0], y = canvas.height - projected[1];
                const pick = await renderer.pick(x, y);
                assert(pick && !isEmptyLoci(native.getLoci(pick.id).loci), 'The interactive target must resolve to a visible atom.');
                const bounds = canvas.getBoundingClientRect();
                return { x: bounds.left + x / native.input.pixelRatio, y: bounds.top + y / native.input.pixelRatio };
            },
            async emptyTarget() {
                const bounds = canvas.getBoundingClientRect();
                for (const [x, y] of [[4, 4], [canvas.width - 5, 4], [4, canvas.height - 5], [canvas.width - 5, canvas.height - 5]]) {
                    if (!(await renderer.pick(x, y))) return { x: bounds.left + x / native.input.pixelRatio, y: bounds.top + y / native.input.pixelRatio };
                }
                throw new Error('The interaction fixture must provide an empty background click target.');
            },
            setSelectionMode(enabled) { plugin.selectionMode = enabled; },
            dispose() { clickSub.unsubscribe(); hoverSub.unsubscribe(); plugin.dispose(); container.remove(); },
        };
        (window as unknown as { webgpuInteraction: WebGPUInteractionVerification }).webgpuInteraction = api;
    } catch (error) { plugin.dispose(); container.remove(); throw error; }
}

(window as unknown as { webgpuPrepareInteraction: typeof prepareInteractionVerification }).webgpuPrepareInteraction = prepareInteractionVerification;

async function verifyGaussianCompute(context: WebGPUContext) {
    const position = {
        indices: OrderedSet.ofSortedArray(new Int32Array([0, 2, 4])),
        x: new Float32Array([-0.43, 100, 0.72, 100, 1.91]),
        y: new Float32Array([0.18, 100, -0.69, 100, 0.87]),
        z: new Float32Array([-0.24, 100, 0.31, 100, -0.58]),
        id: new Float32Array([7, 11, 19]),
    };
    const box = Box3D.create(Vec3.create(-0.43, -0.69, -0.58), Vec3.create(1.91, 0.87, 0.31));
    const radius = (index: number) => [1.03, 10, 0.67, 10, 1.31][index];
    for (const props of [{ resolution: 0.55, radiusOffset: 0.2, smoothness: 1.5 }, { resolution: 0.9, radiusOffset: 0, smoothness: 2.1 }]) {
        const computed = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, position, box, radius, props);
        const batched = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, position, box, radius, props, { maxVoxelsPerBatch: 113 });
        const fullScan = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, position, box, radius, props, { spatialBins: false, maxVoxelsPerBatch: 127 });
        assert(computed.field.data.every((v, i) => v === fullScan.field.data[i]) && computed.idField.data.every((v, i) => v === fullScan.idField.data[i]), 'Spatial bins must preserve all density values and atom IDs exactly against a full atom scan.');
        const cpu = await GaussianDensityCPU(RuntimeContext.Synchronous, position, box, radius, props);
        assert(computed.field.space.dimensions.every((v, i) => v === cpu.field.space.dimensions[i]), 'Native Gaussian dimensions must preserve the CPU grid padding and tensor order.');
        assert(computed.transform.every((v, i) => v === cpu.transform[i]) && computed.maxRadius === cpu.maxRadius && computed.radiusFactor === 1, 'Native Gaussian grids must preserve physical grid transforms and atom radii.');
        assert(computed.field.data.every((v, i) => v === batched.field.data[i]) && computed.idField.data.every((v, i) => v === batched.idField.data[i]), 'Batch boundaries must preserve every Gaussian voxel and atom ID.');
        const dimensions = computed.field.space.dimensions;
        for (let i = 0; i < computed.field.data.length; i++) {
            const coordinate = [Math.floor(i / (dimensions[1] * dimensions[2])), Math.floor(i / dimensions[2]) % dimensions[1], i % dimensions[2]];
            const point = coordinate.map((v, c) => Math.fround(Math.fround(computed.transform[12 + c]) + Math.fround(v * Math.fround(props.resolution))));
            let total = 0, strongest = 0, id = -1;
            for (let a = 0; a < 3; a++) {
                const index = [0, 2, 4][a], r = Math.fround(radius(index) + props.radiusOffset);
                const distanceSquared = (point[0] - position.x[index]) ** 2 + (point[1] - position.y[index]) ** 2 + (point[2] - position.z[index]) ** 2;
                if (distanceSquared > 4 * r * r) continue;
                const contribution = Math.exp(-props.smoothness * distanceSquared / (r * r));
                total = Math.fround(total + contribution);
                if (contribution > strongest) { strongest = contribution; id = position.id[a]; }
            }
            assert(Math.abs(computed.field.data[i] - total) <= 0.00002, 'Native Gaussian values must agree with the analytical truncated exponential kernel.');
            assert(computed.idField.data[i] === id, 'Native Gaussian IDs must identify the strongest contributing atom, with -1 outside the cutoff.');
            // CPU uses a fast exponential approximation; native WGSL follows the existing GPU exp kernel.
            assert(Math.abs(computed.field.data[i] - cpu.field.data[i]) <= 0.065 * Math.max(cpu.field.data[i], 0.05), 'Native Gaussian density must remain within the CPU exponential approximation error.');
        }
    }
    for (const shift of [0, 100000]) {
        const dispersed = { indices: OrderedSet.ofBounds(0, 4), x: [shift, shift + 9, shift + 19, shift + 28], y: [0, 3, -2, 1], z: [0, -1, 2, 0], id: [4, 3, 2, 1] };
        const dispersedBox = Box3D.create(Vec3.create(shift, -2, -1), Vec3.create(shift + 28, 3, 2));
        const props = { resolution: 0.75, radiusOffset: 0, smoothness: 1.5 };
        const indexed = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, dispersed, dispersedBox, () => 1, props);
        const fullScan = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, dispersed, dispersedBox, () => 1, props, { spatialBins: false });
        assert(indexed.field.data.every((v, i) => v === fullScan.field.data[i]) && indexed.idField.data.every((v, i) => v === fullScan.idField.data[i]), 'Spatial binning must preserve density and IDs for separated atoms and large coordinate offsets.');
    }
    const empty = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, { ...position, indices: OrderedSet.ofBounds(0, 0) }, box, radius, { resolution: 1, radiusOffset: 0, smoothness: 1.5 });
    assert(empty.field.data.every(v => v === 0) && empty.idField.data.every(v => v === -1), 'Empty Gaussian selections must produce zero density and unassigned atom IDs.');
    const ties = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, { indices: OrderedSet.ofBounds(0, 2), x: [0, 0], y: [0, 0], z: [0, 0], id: [11, 7] }, Box3D.create(Vec3(), Vec3()), () => 1, { resolution: 0.5, radiusOffset: 0, smoothness: 1.5 }, { maxVoxelsPerBatch: 37 });
    assert(ties.idField.data.every((v, i) => v === (ties.field.data[i] > 0 ? 11 : -1)), 'Equal-density atom ties must retain the first input ID in every batch.');
    const zeroRadius = await GaussianDensityWebGPU(RuntimeContext.Synchronous, context, { indices: OrderedSet.ofBounds(0, 1), x: [0], y: [0], z: [0] }, Box3D.create(Vec3(), Vec3()), () => 0, { resolution: 1, radiusOffset: 0, smoothness: 1.5 });
    assert(zeroRadius.field.data.every(v => v === 0) && zeroRadius.idField.data.every(v => v === -1), 'Zero-radius atoms must contribute no density or invalid floating-point values.');
}

async function verifyMarchingCubesCompute(context: WebGPUContext) {
    for (let cube = 0; cube < 256; cube++) {
        const space = Tensor.Space([2, 2, 2], [0, 1, 2], Float32Array), data = space.create();
        for (let i = 0; i < 8; i++) { const p = CubeVertices[i]; space.set(data, p.i, p.j, p.k, cube & (1 << i) ? 0 : 1); }
        for (const isoLevel of [0.5, -0.5]) {
            const scalarField = Tensor.create(space, isoLevel > 0 ? data : Tensor.Data1(Float32Array.from(data, v => v - 1)));
            const atomSpace = Tensor.Space([2, 2, 2], [2, 1, 0], Float32Array), atomData = atomSpace.create();
            for (let i = 0; i < 8; i++) { const p = CubeVertices[i]; atomSpace.set(atomData, p.i, p.j, p.k, 17 + i); }
            const params = { scalarField, isoLevel, idField: Tensor.create(atomSpace, atomData) };
            const gpu = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, params);
            const cpu = await computeMarchingCubesMesh(params).run();
            assert(gpu.vertexCount === cpu.triangleCount * 3, `Cube case ${cube} must preserve triangle counts.`);
            assert(gpu.readbackBytes === 4 + gpu.vertexCount * 48, 'Compaction must read back only the count and generated vertices, including empty cube cases.');
            if ([0, 1, 23, 254, 255].includes(cube) && isoLevel > 0) {
                const transform = Mat4.fromScaling(Mat4(), Vec3.create(2, 3, 4)); Mat4.setTranslation(transform, Vec3.create(5, -2, 1));
                const surface = await computeMarchingCubesTextureMeshWebGPU(RuntimeContext.Synchronous, context, params, transform, 25, Sphere3D.create(Vec3.create(6, 0, 3), 5));
                const native = surface.meta.webgpuGeometry;
                assert(native instanceof WebGPUTextureMeshGeometry && surface.vertexCount === gpu.vertexCount, 'Native compute output must preserve populated and empty surface topology.');
                assert(surface.meta.webgpuReadbackBytes === 4 + Math.max(80, gpu.vertexCount * 80), 'Native texture extraction must read one final compatibility mirror instead of reading and reuploading compacted triangles.');
                try {
                    const positionData = await native.position.readData(), normalData = await native.normal.readData(), groupData = await native.group.readData();
                    for (const [actual, texture] of [[positionData, native.position], [normalData, native.normal], [groupData, native.group]] as const) assert(actual.array.every((v, i) => v === texture.data.array[i]), 'GPU-generated textures and CPU theme/location mirrors must agree including padding.');
                    const normalTransform = Mat3.directionTransform(Mat3(), transform), normal = Vec3();
                    for (let i = 0; i < gpu.vertexCount; i++) {
                        const position = Vec3.transformMat4(Vec3(), Vec3.fromArray(Vec3(), gpu.positions, i * 3), transform);
                        Vec3.normalize(normal, Vec3.transformMat3(normal, Vec3.fromArray(normal, gpu.normals, i * 3), normalTransform));
                        for (let a = 0; a < 3; a++) {
                            assert(Math.abs(native.geometry.vertices[i * 20 + a] - position[a]) < 0.00001, 'GPU texture positions must preserve affine world transforms.');
                            assert(Math.abs(native.geometry.vertices[i * 20 + a + 4] - normal[a]) < 0.00001, 'GPU normals must use inverse-transpose transforms and normalization.');
                        }
                        assert(native.geometry.vertices[i * 20 + 16] === gpu.groups[i], 'Direct GPU storage vertices must preserve molecular group IDs.');
                    }
                } finally { native.destroy(); surface.doubleBuffer.destroy(); }
            }
            for (let i = 0; i < gpu.vertexCount; i++) for (let a = 0; a < 3; a++) {
                const index = cpu.indexBuffer.ref.value[i];
                assert(gpu.groups[i] === cpu.groupBuffer.ref.value[index], 'Every cube case must preserve atom IDs under the CPU edge conventions.');
                assert(Math.abs(gpu.positions[i * 3 + a] - cpu.vertexBuffer.ref.value[index * 3 + a]) < 0.00001, `Cube case ${cube} must preserve interpolation and winding.`);
                assert(Math.abs(gpu.normals[i * 3 + a] - cpu.normalBuffer.ref.value[index * 3 + a]) < 0.00001, `Cube case ${cube} must preserve interpolated normals.`);
            }
        }
    }
    const sparseSpace = Tensor.Space([18, 6, 4], [2, 0, 1], Float32Array), sparseData = sparseSpace.create();
    for (let x = 0; x < 18; x++) for (let y = 0; y < 6; y++) for (let z = 0; z < 4; z++) sparseSpace.set(sparseData, x, y, z, Math.exp(-((x - 9) ** 2 + (y - 3) ** 2 + (z - 2) ** 2) / 8));
    const sparseParams = { scalarField: Tensor.create(sparseSpace, sparseData), isoLevel: 0.5 };
    const sparse = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, sparseParams, { groupMode: 'cell' });
    assert(sparse.readbackBytes < 17 * 5 * 3 * 736 / 2, 'Sparse surface compaction must eliminate more than half of the fixed-cell readback data.');
    for (const maxCellsPerBatch of [1, 63, 64, 65, 127, 129]) {
        const batched = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, sparseParams, { groupMode: 'cell', maxCellsPerBatch });
        assert(sparse.positions.every((v, i) => v === batched.positions[i]) && sparse.normals.every((v, i) => v === batched.normals[i]) && sparse.groups.every((v, i) => v === batched.groups[i]) && sparse.vertexCount === batched.vertexCount, 'Prefix scans must preserve vertex order across workgroup boundaries, partial groups and empty cells.');
        assert(batched.readbackBytes === Math.ceil(255 / maxCellsPerBatch) * 4 + batched.vertexCount * 48, 'Compaction readback must exclude every unused cell slot in every batch.');
    }
    for (const order of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]) {
        const space = Tensor.Space([3, 4, 5], order, Float32Array), data = space.create();
        for (let x = 0; x < 3; x++) for (let y = 0; y < 4; y++) for (let z = 0; z < 5; z++) space.set(data, x, y, z, Math.sin(x * 2 * Math.PI / 3) + 0.4 * Math.cos(y * 2 * Math.PI / 4) + 0.2 * Math.sin(z * 2 * Math.PI / 5));
        const params = { scalarField: createWrappedTensor(Tensor.create(space, data)), isoLevel: 0.123 };
        const gpu = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, params, { groupMode: 'cell', periodicDimensions: space.dimensions, maxCellsPerBatch: 7 });
        const cpu = await computeMarchingCubesMesh(params).run();
        assert(gpu.vertexCount === cpu.triangleCount * 3, 'Wrapped extraction must preserve CPU surface topology in all tensor axis orders.');
        assert(gpu.positions.some((v, i) => v > space.dimensions[i % 3] - 1), 'Wrapped extraction must include triangles crossing the final periodic boundary.');
        const coordinate = [0, 0, 0];
        for (let i = 0; i < gpu.vertexCount; i++) {
            assert(gpu.groups[i] >= 0 && gpu.groups[i] < data.length && Number.isInteger(gpu.groups[i]), 'Wrapped surface IDs must refer to cells in the original volume.');
            space.getCoords(gpu.groups[i], coordinate);
            for (let a = 0; a < 3; a++) {
                const value = gpu.positions[i * 3 + a], cpuIndex = cpu.indexBuffer.ref.value[i] * 3 + a;
                assert(Math.abs(value - cpu.vertexBuffer.ref.value[cpuIndex]) < 0.00001 && Math.abs(gpu.normals[i * 3 + a] - cpu.normalBuffer.ref.value[cpuIndex]) < 0.00001, 'Periodic GPU interpolation and boundary normals must agree with CPU extraction.');
                assert(value >= coordinate[a] - 0.00001 && value <= coordinate[a] + 1.00001, 'Wrapped picking IDs must identify the original periodic cell containing each triangle.');
            }
        }
    }
    const cavitySpace = Tensor.Space([5, 5, 5], [1, 2, 0], Float32Array), cavityData = cavitySpace.create();
    for (let x = 1; x < 4; x++) for (let y = 1; y < 4; y++) for (let z = 1; z < 4; z++) cavitySpace.set(cavityData, x, y, z, 1);
    cavitySpace.set(cavityData, 2, 2, 2, 0);
    for (const mode of ['inside', 'outside'] as const) {
        const params = { scalarField: Tensor.createFloodfilled(Tensor.create(cavitySpace, cavityData), 0.4, mode), isoLevel: 0.4 };
        const gpu = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, params);
        const cpu = await computeMarchingCubesMesh(params).run();
        assert(gpu.vertexCount === cpu.triangleCount * 3, 'GPU uploads must honor floodfilled tensor accessors without modifying the backing array.');
        for (let i = 0; i < gpu.vertexCount; i++) for (let a = 0; a < 3; a++) {
            const index = cpu.indexBuffer.ref.value[i] * 3 + a;
            assert(Math.abs(gpu.positions[i * 3 + a] - cpu.vertexBuffer.ref.value[index]) < 0.00001 && Math.abs(gpu.normals[i * 3 + a] - cpu.normalBuffer.ref.value[index]) < 0.00001, 'Filled cavity and exterior surfaces must preserve CPU geometry and normals.');
        }
        assert(cavitySpace.get(cavityData, 2, 2, 2) === 0, 'GPU extraction must preserve original scalar samples under tensor views.');
    }
    const cropSpace = Tensor.Space([5, 6, 7], [1, 0, 2], Float32Array), cropData = cropSpace.create();
    for (let x = 0; x < 5; x++) for (let y = 0; y < 6; y++) for (let z = 0; z < 7; z++) cropSpace.set(cropData, x, y, z, x + y * 0.7 + z * 0.3);
    for (const [bottomLeft, topRight] of [[[0, 0, 0], [2, 3, 4]], [[1, 1, 1], [4, 5, 6]], [[3, 2, 3], [5, 6, 7]], [[0, 0, 0], [5, 6, 7]]]) {
        const isoLevel = bottomLeft.reduce((sum, v, a) => sum + (v + topRight[a] - 1) * 0.5 * [1, 0.7, 0.3][a], 0) + 0.125;
        const params = { scalarField: Tensor.create(cropSpace, cropData), isoLevel, bottomLeft, topRight };
        const gpu = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, params, { groupMode: 'cell', maxCellsPerBatch: 7 });
        const cpu = await computeMarchingCubesMesh(params).run();
        assert(gpu.vertexCount > 0 && gpu.vertexCount === cpu.triangleCount * 3, 'Cropped native surfaces must preserve nonempty CPU triangle counts.');
        const coordinate = [0, 0, 0];
        for (let i = 0; i < gpu.vertexCount; i++) {
            cropSpace.getCoords(gpu.groups[i], coordinate);
            for (let a = 0; a < 3; a++) {
                const index = cpu.indexBuffer.ref.value[i] * 3 + a;
                assert(Math.abs(gpu.positions[i * 3 + a] - cpu.vertexBuffer.ref.value[index]) < 0.00001 && Math.abs(gpu.normals[i * 3 + a] - cpu.normalBuffer.ref.value[index]) < 0.00001, 'Cropped extraction must retain global positions and full-grid boundary gradients.');
                assert(coordinate[a] >= bottomLeft[a] && coordinate[a] < topRight[a] - 1, 'Cropped cell IDs must retain original tensor offsets.');
            }
        }
    }
    const beforeEmpty = context.stats.marchingCubesDispatches;
    const emptyRegion = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, { scalarField: Tensor.create(cropSpace, cropData), isoLevel: 4.125, bottomLeft: [1, 1, 1], topRight: [1, 3, 3] });
    assert(emptyRegion.vertexCount === 0 && emptyRegion.readbackBytes === 0 && context.stats.marchingCubesDispatches === beforeEmpty, 'Empty regions must produce no surface without dispatching GPU work.');
    let invalidRegion = false;
    try { await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, { scalarField: Tensor.create(cropSpace, cropData), isoLevel: 4.125, bottomLeft: [-1, 0, 0] }); } catch (error) { invalidRegion = /region/.test(String(error)); }
    assert(invalidRegion, 'Out-of-grid surface regions must fail before GPU access.');
    const idSpace = Tensor.Space([2, 2, 2], [2, 1, 0], Float32Array), idScalar = idSpace.create(), atomIds = idSpace.create();
    idScalar.fill(1); idSpace.set(idScalar, 0, 0, 0, 0); atomIds.fill(17);
    const idAxisSpace = Tensor.Space([2, 2, 2], [0, 1, 2], Float32Array);
    const idParams = { scalarField: Tensor.create(idSpace, idScalar), idField: Tensor.create(idAxisSpace, atomIds), isoLevel: 0.5 };
    const identified = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, idParams);
    assert(identified.vertexCount === 3 && identified.groups.every(v => v === 17), 'GPU extraction must preserve molecular atom IDs from tensors with different axis orders.');
    atomIds.fill(-2);
    const ignored = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, idParams);
    assert(ignored.vertexCount === 0 && ignored.readbackBytes === 4, 'Ignored molecular IDs must remove triangles before GPU compaction.');
    atomIds.fill(17); idAxisSpace.set(atomIds, 0, 0, 0, -1);
    const fallback = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, idParams);
    assert(fallback.groups.every(v => v === 17), 'Missing molecular IDs must fall back to the other edge endpoint.');
    for (const order of [[0, 1, 2], [2, 1, 0], [1, 0, 2]]) {
        const space = Tensor.Space([5, 4, 3], order, Float32Array), data = space.create();
        for (let x = 0; x < 5; x++) for (let y = 0; y < 4; y++) for (let z = 0; z < 3; z++) space.set(data, x, y, z, x + y * 0.7 + z * 0.3);
        const params = { scalarField: Tensor.create(space, data), isoLevel: 2.1 };
        const gpu = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, params, { groupMode: 'cell' });
        const batched = await computeMarchingCubesWebGPU(RuntimeContext.Synchronous, context, params, { groupMode: 'cell', maxCellsPerBatch: 7 });
        assert(gpu.positions.every((v, i) => v === batched.positions[i]) && gpu.normals.every((v, i) => v === batched.normals[i]) && gpu.groups.every((v, i) => v === batched.groups[i]), 'Marching-cubes batches must preserve geometry, normals and original-axis cell IDs.');
        const cpu = await computeMarchingCubesMesh(params).run();
        assert(gpu.vertexCount === cpu.triangleCount * 3, 'Marching cubes must preserve multi-cell topology across tensor axis orders.');
        for (let i = 0; i < gpu.vertexCount; i++) for (let a = 0; a < 3; a++) assert(Math.abs(gpu.positions[i * 3 + a] - cpu.vertexBuffer.ref.value[cpu.indexBuffer.ref.value[i] * 3 + a]) < 0.00001, 'Marching cubes must preserve multi-cell vertex order and positions.');
    }
}

async function verifyOrbitalCompute(context: WebGPUContext) {
    const near = (actual: ArrayLike<number>, expected: ArrayLike<number>, label: string) => {
        assert(actual.length === expected.length, `${label}: grid dimensions must agree.`);
        for (let i = 0; i < actual.length; i++) assert(Math.abs(actual[i] - expected[i]) <= 2e-5 * Math.max(1, Math.abs(expected[i])), `${label}: voxel ${i}: ${actual[i]} versus ${expected[i]}.`);
    };
    const single = initCubeGrid({ basis: { atoms: [{ center: [0, 0, 0], shells: [{ angularMomentum: [0], exponents: [1], coefficients: [[1]] }] }] }, sphericalOrder: 'gaussian', cutoffThreshold: 0, gridSpacing: 0.4, boxExpand: 1 });
    const sOrbital: AlphaOrbital = { alpha: [1], occupancy: 1, energy: 0 };
    near(await computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, single, [sOrbital]), await sphericalCollocation(single, sOrbital, RuntimeContext.Synchronous), 'Single primitive and single alpha coefficient');
    for (const run of [
        () => computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, single, [{ ...sOrbital, alpha: [] }]),
        () => computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, single, [sOrbital], false, { maxVoxelsPerBatch: 0 }),
        () => computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, { ...single, delta: Vec3.create(NaN, 1, 1) }, [sOrbital]),
    ]) {
        const before = context.stats.computeDispatches;
        let rejected = false;
        try { await run(); } catch { rejected = true; }
        assert(rejected && context.stats.computeDispatches === before, 'Invalid orbital inputs must fail before any GPU dispatch.');
    }
    for (const sphericalOrder of ['gaussian', 'cca', 'cca-reverse'] as const) {
        for (const cutoffThreshold of [0, 0.15]) {
            const params: CubeGridComputationParams = {
                basis: { atoms: [
                    { center: [-0.3, 0.2, -0.1], shells: [{ angularMomentum: [0, 1, 2, 3, 4], exponents: [0.5, 1.7], coefficients: [[0.7, -0.2], [0.6, 0.3], [-0.4, 0.8], [0.2, 0.1], [0.3, -0.1]] }] },
                    { center: [0.8, -0.2, 0.4], shells: [{ angularMomentum: [0, 1], exponents: [0.8], coefficients: [[0.3], [-0.5]] }] },
                ] }, sphericalOrder, cutoffThreshold, gridSpacing: 0.6, boxExpand: 2,
            };
            const grid = initCubeGrid(params);
            const alpha = Array.from({ length: 29 }, (_, i) => Math.sin(i * 1.3) * 0.7);
            const orbital: AlphaOrbital = { alpha, occupancy: 2, energy: -0.5 };
            // Isolate each angular momentum as well as testing their sum, so an
            // incorrect polynomial cannot hide behind cancellation.
            for (const l of [-1, 0, 1, 2, 3, 4]) {
                const selected = { ...orbital, alpha: alpha.map((a, i) => l < 0 || (i >= l * l && i < (l + 1) ** 2) ? a : 0) };
                const cpu = await sphericalCollocation(grid, selected, RuntimeContext.Synchronous);
                const gpu = await computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, grid, [selected]);
                near(gpu, cpu, `Orbital ${sphericalOrder} L=${l} cutoff=${cutoffThreshold}`);
                const batched = await computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, grid, [selected], false, { maxVoxelsPerBatch: 65 });
                assert(gpu.every((v, i) => v === batched[i]), 'Partial workgroups and orbital batch boundaries must preserve all values exactly.');
            }
            const orbitals = [orbital, { ...orbital, alpha: alpha.map(v => -v * 0.3), occupancy: 0.5 }, { ...orbital, occupancy: 0 }];
            const expected = new Float32Array(grid.npoints);
            for (const o of orbitals) {
                const cpu = await sphericalCollocation(grid, o, RuntimeContext.Synchronous);
                for (let i = 0; i < cpu.length; i++) expected[i] += o.occupancy * cpu[i] * cpu[i];
            }
            const density = await computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, grid, orbitals, true, { maxVoxelsPerBatch: 127 });
            near(density, expected, 'Electron density with fractional/zero occupancies');
            assert((await computeOrbitalGridWebGPU(RuntimeContext.Synchronous, context, grid, [], true)).every(v => v === 0), 'An empty density must contain zero values.');
            const task = await createSphericalCollocationGrid(params, orbital, undefined, context).run();
            near(task.grid.cells.data, await sphericalCollocation(grid, orbital, RuntimeContext.Synchronous), 'Public orbital task');
            assert(task.isovalues?.negative !== undefined && task.isovalues?.positive !== undefined && Number.isFinite(task.grid.stats.sigma), 'Orbital tasks must preserve signed isovalues and statistics.');
            const densityTask = await createSphericalCollocationDensityGrid(params, orbitals, undefined, context).run();
            near(densityTask.grid.cells.data, expected, 'Public density task');
            assert(densityTask.isovalues?.positive !== undefined && densityTask.grid.stats.min >= 0, 'Density tasks must preserve nonnegative statistics and isovalues.');
            const cpuTask = await createSphericalCollocationGrid(params, orbital).run();
            assert(JSON.stringify(task.grid.transform) === JSON.stringify(cpuTask.grid.transform), 'Orbital task transforms must preserve Bohr-to-Angstrom conversion.');
        }
    }
}

async function verifyOrbitalVisualization(plugin: PluginContext) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    native.clear();
    const before = native.webgpu!.stats.computeDispatches;
    const basis = await plugin.build().toRoot().apply(StaticBasisAndOrbitals, {
        basis: { atoms: [{ center: [0, 0, 0], shells: [{ angularMomentum: [1], exponents: [0.5, 1.5], coefficients: [[0.8, 0.2]] }] }] },
        order: 'gaussian', orbitals: [{ alpha: [0, 1, 0], energy: -0.5, occupancy: 2 }, { alpha: [0, 0, 1], energy: -0.3, occupancy: 0.5 }],
    }).commit();
    try {
        const update = plugin.build();
        const volume = update.to(basis).apply(CreateOrbitalVolume, { index: 0, boxExpand: 4, gridSpacing: [{ atomCount: 0, spacing: 0.25 }] });
        const positive = volume.apply(CreateOrbitalRepresentation3D, { kind: 'positive', color: Color(0x0000ff), alpha: 1 }).selector;
        const negative = volume.apply(CreateOrbitalRepresentation3D, { kind: 'negative', color: Color(0xff0000), alpha: 1 }).selector;
        await update.commit();
        assert(volume.selector.isOk && positive.isOk && negative.isOk, 'Orbital extension state transforms must initialize without WebGL.');
        const positiveRepr = positive.obj!.data.repr, negativeRepr = negative.obj!.data.repr;
        native.add(positiveRepr); native.add(negativeRepr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        assert(native.webgpu!.stats.computeDispatches > before && positiveRepr.renderObjects[0].values.drawCount.ref.value > 0 && negativeRepr.renderObjects[0].values.drawCount.ref.value > 0, 'Both signed orbital surfaces must use native computation and contain triangles.');
        const frame = await renderer.readPixels();
        const coloredPixel = (pixels: typeof frame, channel: number) => {
            const candidates: number[] = [];
            let x = 0, y = 0;
            for (let p = 0; p < pixels.array.length / 4; p++) {
                if (pixels.array[p * 4 + channel] <= 50 || pixels.array[p * 4 + 3] <= 0 || pixels.array[p * 4 + (channel + 2) % 3] >= 10) continue;
                candidates.push(p); x += p % pixels.width; y += Math.floor(p / pixels.width);
            }
            x /= candidates.length; y /= candidates.length;
            candidates.sort((a, b) => (a % pixels.width - x) ** 2 + (Math.floor(a / pixels.width) - y) ** 2 - (b % pixels.width - x) ** 2 - (Math.floor(b / pixels.width) - y) ** 2);
            return candidates[0] ?? -1;
        };
        const blue = coloredPixel(frame, 2), red = coloredPixel(frame, 0);
        assert(blue >= 0 && red >= 0, 'The orbital frame must show both positive and negative lobes.');
        for (const pixel of [blue, red]) {
            const picked = await renderer.pick(pixel % frame.width, Math.floor(pixel / frame.width));
            assert(picked && !isEmptyLoci(native.getLoci(picked.id).loci), 'Both signed orbital lobes must resolve volume loci.');
        }
        const exported = await native.getImagePass({ transparentBackground: true }).getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(exported.data.some((v, i) => i % 4 === 0 && v > 50 && exported.data[i + 2] < 10) && exported.data.some((v, i) => i % 4 === 2 && v > 50 && exported.data[i - 2] < 10), 'Orbital screenshot export must preserve both signed colors.');
        showPixels('Native positive and negative orbital lobes', frame);
        await plugin.build().to(positive).update({ pickable: false }).commit(); native.update(positiveRepr); native.tick(now());
        const pixel = blue;
        assert(!(await renderer.pick(pixel % frame.width, Math.floor(pixel / frame.width))), 'Orbital representation pickable control must exclude its surface.');
        await plugin.build().to(positive).update({ pickable: true }).commit(); native.update(positiveRepr); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === frame.array[i]), 'Pickability changes must preserve orbital colors.');
        const dispatches = native.webgpu!.stats.computeDispatches;
        await plugin.build().to(volume.selector).update({ index: 1 }).commit(); native.update(positiveRepr); native.update(negativeRepr); native.tick(now());
        assert(native.webgpu!.stats.computeDispatches > dispatches && (await renderer.readPixels()).array.some((v, i) => v !== frame.array[i]), 'Orbital index changes must recompute and update both signed surfaces.');
        native.remove(positiveRepr); native.remove(negativeRepr);
        const densityUpdate = plugin.build();
        const density = densityUpdate.to(basis).apply(CreateOrbitalDensityVolume, { boxExpand: 4, gridSpacing: [{ atomCount: 0, spacing: 0.25 }] });
        const densitySurface = density.apply(CreateOrbitalRepresentation3D, { kind: 'positive', color: Color(0x00ff00), alpha: 1 }).selector;
        await densityUpdate.commit();
        assert(densitySurface.isOk && density.selector.isOk && density.selector.obj!.data.grid.stats.min >= 0, 'Density extension state transforms must compute nonnegative grids without WebGL.');
        const repr = densitySurface.obj!.data.repr;
        native.add(repr); native.requestCameraReset({ durationMs: 0 }); native.tick(now());
        const densityFrame = await renderer.readPixels();
        const green = coloredPixel(densityFrame, 1);
        assert(green >= 0 && repr.renderObjects[0].values.drawCount.ref.value > 0, 'Computed electron-density surfaces must render.');
        const densityPixel = green, pick = await renderer.pick(densityPixel % densityFrame.width, Math.floor(densityPixel / densityFrame.width));
        assert(pick && !isEmptyLoci(native.getLoci(pick.id).loci), 'Computed electron-density surfaces must retain volume loci.');
        showPixels('Native orbital electron density', densityFrame);
        native.remove(repr);
        return ['orbital extension state transforms, signed lobes, index updates, pickability, volume loci, screenshots and electron-density surfaces without WebGL'];
    } finally { await plugin.build().delete(basis).commit(); }
}

async function verifyColorSmoothingCompute(context: WebGPUContext) {
    for (const itemSize of [1, 3, 4] as const) for (const colorType of ['group', 'groupInstance'] as const) for (const stride of [1, 2]) {
        const input: ColorSmoothingInput = {
            vertexCount: 5, instanceCount: 2, groupCount: 3, itemSize, colorType,
            positionBuffer: new Float32Array([-0.31, 0.42, 0.19, 0.74, -0.51, -0.28, 0.36, 0.25, 0.87, -0.6, -0.32, 0.7, 0, 0, 0]),
            groupBuffer: new Float32Array([0, 1, 2, 1, 0]), instanceBuffer: new Float32Array([1, 0]),
            transformBuffer: new Float32Array([...Mat4.identity(), ...Mat4.fromTranslation(Mat4(), Vec3.create(2.1, -0.3, 0.2))]),
            colorData: { width: 6, height: 1, array: Uint8Array.from({ length: 6 * itemSize }, (_, i) => (i * 67 + 13) % 256) },
            invariantBoundingSphere: Sphere3D.create(Vec3(), 1.3), boundingSphere: Sphere3D.create(Vec3.create(1, 0, 0), 2.5),
        };
        const options = { resolution: 0.55, stride };
        const reference = calcMeshColorSmoothing(input, options, undefined, new WebGPUTextureData());
        assert(reference.kind === 'volume' && reference.texture instanceof WebGPUTextureData, 'CPU reference must expose its packed smoothing grid.');
        const gpu = await calcMeshColorSmoothingWebGPU(context, input, options);
        assert(gpu.gridDim.every((v, i) => v === reference.gridDim[i]) && gpu.gridTransform.every((v, i) => v === reference.gridTransform[i]) && gpu.type === reference.type, 'Native smoothing must preserve grid layout, transform and theme granularity.');
        const expected = reference.texture.data.array, actual = gpu.texture.data.array;
        assert(actual.length === expected.length && actual.every((v, i) => Math.abs(v - expected[i]) <= 1), `Native ${itemSize}-channel ${colorType} smoothing must match weighted CPU accumulation within one byte.`);
        const direct = calcMeshColorSmoothingTextureWebGPU(context, input, options);
        assert(direct.texture.native?.device === context.device && direct.texture.getByteCount() === actual.length, 'Smoothing must expose an owned native GPU grid with accurate memory size.');
        const directPixels = await context.readTexture(direct.texture.native!.texture, 0, 0, direct.gridTexDim[0], direct.gridTexDim[1], 4);
        assert(directPixels.every((v, i) => v === actual[i]), 'Direct texture compute must preserve all packed bytes, including unused slice padding.');
        let readbackRequired = false;
        try { void direct.texture.data; } catch { readbackRequired = true; }
        assert(readbackRequired, 'GPU-owned texture grids must require explicit readback instead of exposing stale CPU data.');
        const originalTexture = direct.texture.native!.texture;
        const reused = calcMeshColorSmoothingTextureWebGPU(context, input, options, direct.texture);
        assert(reused.texture === direct.texture && reused.texture.native!.texture !== originalTexture, 'Updating a native grid must replace its GPU allocation while preserving the theme texture wrapper.');
        direct.texture.destroy();
        assert(!direct.texture.native, 'Destroying a GPU grid must release the owned texture.');
        const batched = await calcMeshColorSmoothingWebGPU(context, input, options, new WebGPUTextureData(), 65);
        assert(batched.texture.data.array.every((v, i) => v === actual[i]), 'Smoothing batch boundaries and partial workgroups must preserve the grid exactly.');
        const empty = await calcMeshColorSmoothingWebGPU(context, { ...input, vertexCount: 0 }, options);
        assert(empty.texture.data.array.every((v, i) => v === (itemSize === 3 && i % 4 === 3 ? 255 : 0)), 'Empty smoothing grids must contain zero data with the RGB alpha convention.');
        const before = context.stats.computeDispatches;
        let rejected = false;
        try { await calcMeshColorSmoothingWebGPU(context, input, { resolution: 0, stride }); } catch { rejected = true; }
        assert(rejected && context.stats.computeDispatches === before, 'Invalid smoothing parameters must fail before dispatch.');
    }
}

async function verifyOutlineLayers(renderer: WebGPURenderer, camera: Camera, props: RendererProps) {
    const meshProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
    const triangle = (scale: number, z: number, alpha: number) => {
        const mesh = Mesh.create(new Float32Array([-4 * scale, -4 * scale, z, 4 * scale, -4 * scale, z, 0, 4 * scale, z]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
        const values = Mesh.Utils.createValuesSimple(mesh, { ...meshProps, alpha }, Color(0xff0000), 1);
        return createRenderObject('mesh', values, Mesh.Utils.createRenderableState({ ...meshProps, alpha }), -1);
    };
    const opaque = triangle(1, 0, 1), front = triangle(2, 2, 0.25);
    const post = PD.getDefaultValues(PostprocessingParams);
    post.outline = { name: 'on', params: { ...PD.getDefaultValues(OutlineParams), color: Color(0x0000ff) } }; post.antialiasing = { name: 'off', params: {} };
    renderer.render([opaque], camera, props, true, 1, undefined, post);
    const reference = await renderer.readPixels();
    const edges: number[] = [];
    for (let i = 0; i < reference.array.length; i += 4) if (reference.array[i + 2] > 250 && reference.array[i] < 5) edges.push(i);
    assert(edges.length > 0, 'Opaque reference triangles must expose an outline silhouette.');
    front.state.pickable = false;
    for (const includeTransparent of [true, false]) {
        post.outline.params.includeTransparent = includeTransparent;
        renderer.render([opaque, front], camera, props, true, 1, undefined, post);
        const overlapped = await renderer.readPixels();
        for (const offset of edges) {
            assert(Math.abs(overlapped.array[offset] - 64) <= 1 && Math.abs(overlapped.array[offset + 2] - 191) <= 1 && overlapped.array[offset + 3] === 255, 'Opaque outlines underneath a transparent plane must preserve independent depth and blend its color/opacity.');
        }
        const pixel = edges[0] / 4;
        assert(!(await renderer.pick(pixel % reference.width, Math.floor(pixel / reference.width))), 'An unpickable front plane and an outline outside the opaque silhouette must not produce selection IDs.');
    }
    post.outline.params.includeTransparent = true;
    const nearer = triangle(0.9, 2, 0.1), farther = triangle(1, 0, 0.4);
    const threshold = props.pickingAlphaThreshold; props.pickingAlphaThreshold = 0.2;
    try {
        renderer.render([farther, nearer], camera, props, true, 1, undefined, post);
        const layered = await renderer.readPixels();
        let checked = 0;
        for (const offset of edges) if (layered.array[offset] < 2 && layered.array[offset + 2] > 0) {
            assert(Math.abs(layered.array[offset + 3] - 51) <= 1, 'Transparent outlines must double the closest source opacity rather than composited opacity from all overlapping surfaces.'); checked++;
        }
        assert(checked > 10, 'Overlapping transparent surfaces must expose a measurable outline silhouette.');
        assert((await renderer.pick(128, 128))?.id.objectId === farther.id, 'Faint nearest outline depth must remain independent of opacity-threshold selection.');
        nearer.state.colorOnly = true;
        renderer.render([farther, nearer], camera, props, true, 1, undefined, post);
        assert((await renderer.readPixels()).array.every((v, i) => v === layered.array[i]), 'Color-only flags must not remove transparent visual outline depth.');
        assert((await renderer.pick(128, 128))?.id.objectId === farther.id, 'Color-only nearest surfaces must retain selection of geometry behind them.');
    } finally { props.pickingAlphaThreshold = threshold; }
}

async function verifyVolumeOutlineLayers(plugin: PluginContext) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const objects = native.getRenderObjects(), camera = native.camera, props = native.props.renderer;
    const original = await renderer.readPixels();
    const post = PD.getDefaultValues(PostprocessingParams);
    post.antialiasing = { name: 'off', params: {} };
    renderer.render(objects, camera, props, true, 1, undefined, post);
    const baseline = await renderer.readPixels(), picked = await renderer.pick(160, 120);
    post.outline = { name: 'on', params: { ...PD.getDefaultValues(OutlineParams), includeTransparent: false, color: Color(0x0000ff) } };
    renderer.render(objects, camera, props, true, 1, undefined, post);
    assert((await renderer.readPixels()).array.every((v, i) => v === baseline.array[i]), 'Excluding transparent outlines must preserve ray-marched volumes without opaque surfaces.');
    post.outline.params.includeTransparent = true;
    renderer.render(objects, camera, props, true, 1, undefined, post);
    const included = await renderer.readPixels();
    assert(included.array.some((v, i) => i % 4 === 2 && v > baseline.array[i] + 20), 'Ray-marched volumes must populate the independent transparent outline depth and alpha pass.');
    const outlinedPick = await renderer.pick(160, 120);
    assert(picked && outlinedPick?.id.objectId === picked.id.objectId && outlinedPick.id.groupId === picked.id.groupId, 'Volume outlines must preserve selected objects and source-grid cell IDs.');
    const excludedPass = native.getImagePass({ transparentBackground: true, postprocessing: { ...post, outline: { ...post.outline, params: { ...post.outline.params, includeTransparent: false } } } });
    const includedPass = native.getImagePass({ transparentBackground: true, postprocessing: post });
    try {
        const excludedExport = await excludedPass.getImageData(RuntimeContext.Synchronous, 192, 128);
        const includedExport = await includedPass.getImageData(RuntimeContext.Synchronous, 192, 128);
        assert(includedExport.data.some((v, i) => i % 4 === 2 && v > excludedExport.data[i] + 20), 'Volume screenshot exports must include the independent transparent outline pass.');
    } finally {
        assert(excludedPass instanceof WebGPUImagePass && includedPass instanceof WebGPUImagePass, 'Native volume exports must expose disposable native passes.');
        await excludedPass.dispose(); await includedPass.dispose();
        native.requestDraw(); native.tick(now());
    }
    assert((await renderer.readPixels()).array.every((v, i) => v === original.array[i]), 'Volume outline verification and export disposal must restore the original frame.');
}

async function verifyTransparentSsao(plugin: PluginContext, updateAlpha: (alpha: number) => Promise<unknown>) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const originalPost = native.props.postprocessing, originalFog = native.camera.state.fog;
    const params = { ...PD.getDefaultValues(SsaoParams), samples: 16, radius: 2, bias: 1, blurKernelSize: 7, transparentThreshold: 0 };
    const config = { ...PD.getDefaultValues(PostprocessingParams), antialiasing: { name: 'off' as const, params: {} }, occlusion: { name: 'on' as const, params } };
    try {
        await updateAlpha(0.8);
        native.setProps({ postprocessing: config }); native.requestDraw(); native.tick(now());
        const excluded = await renderer.readPixels();
        let point: { x: number, y: number } | undefined, picked;
        for (let y = 0; y < excluded.height && !picked; y += 5) for (let x = 0; x < excluded.width && !picked; x += 5) {
            if (excluded.array[(y * excluded.width + x) * 4 + 3] > 150) { picked = await renderer.pick(x, y); if (picked) point = { x, y }; }
        }
        assert(picked && point, 'Transparent SSAO reference proteins must retain a selectable atom.');
        native.setProps({ postprocessing: { ...config, occlusion: { name: 'off', params: {} } } }); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - excluded.array[i]) <= 1), 'Threshold zero must exclude transparent occlusion without affecting the transparent color pass.');
        native.setProps({ postprocessing: { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.19 } } } }); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === excluded.array[i]), 'The scene transparency minimum must gate transparent SSAO below its threshold.');
        native.setProps({ postprocessing: { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.21 } } } }); native.tick(now());
        const included = await renderer.readPixels();
        let darkened = 0;
        for (let i = 0; i < included.array.length; i += 4) if (included.array[i] + 3 < excluded.array[i] && excluded.array[i + 3] > 150) darkened++;
        assert(darkened > 50, 'Transparent depth normals and alpha-weighted occluders must darken protein cavities.');
        assert(included.array.every((v, i) => i % 4 !== 3 || v === excluded.array[i]), 'Transparent SSAO must preserve every alpha value.');
        const selected = await renderer.pick(point.x, point.y);
        assert(selected?.id.objectId === picked.id.objectId && selected.id.groupId === picked.id.groupId, 'Transparent SSAO must preserve atom and object selection IDs.');
        for (const resolutionScale of [0.5, 0.1]) {
            native.setProps({ postprocessing: { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.21, resolutionScale } } } }); native.tick(now());
            const reduced = await renderer.readPixels();
            assert(reduced.array.some((v, i) => i % 4 === 0 && v + 3 < excluded.array[i]), 'Reduced-resolution SSAO must retain visible transparent occlusion.');
            assert(reduced.array.some((v, i) => v !== included.array[i]), 'Changing SSAO resolution must change its sampling and blur.');
            assert(reduced.array.every((v, i) => i % 4 !== 3 || v === included.array[i]), 'Upsampled transparent SSAO must preserve every full-resolution alpha.');
            if (resolutionScale === 0.5) {
                renderer.render(native.getRenderObjects(), native.camera, native.props.renderer, native.props.transparentBackground, 2, undefined, { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.21 } } });
                assert((await renderer.readPixels()).array.every((v, i) => v === reduced.array[i]), 'A 2x display pixel ratio must calculate SSAO at half resolution without changing geometry coverage.');
            }

        }
        native.setProps({ postprocessing: { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.21 } } } }); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === included.array[i]), 'Restoring full SSAO resolution must recreate its targets and reproduce the original frame.');

        native.setProps({ postprocessing: { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.21, blurKernelSize: 1 } } } }); native.tick(now());
        assert((await renderer.readPixels()).array.some((v, i) => v !== included.array[i]), 'Transparent occlusion must use its own depth-aware blur.');
        const multiScaleParams = { ...params, transparentThreshold: 0.21, multiScale: { name: 'on' as const, params: { levels: [{ radius: 1, bias: 1 }, { radius: 3, bias: 0.5 }], nearThreshold: 0, farThreshold: 10000 } } };
        native.setProps({ postprocessing: { ...config, occlusion: { name: 'on', params: multiScaleParams } } }); native.tick(now());
        assert((await renderer.readPixels()).array.some((v, i) => v !== included.array[i]), 'Transparent SSAO must honor multi-scale radii and biases.');
        const includedPass = native.getImagePass({ transparentBackground: true, postprocessing: { ...config, occlusion: { name: 'on', params: { ...params, transparentThreshold: 0.21 } } } });
        const excludedPass = native.getImagePass({ transparentBackground: true, postprocessing: config });
        try {
            const withAo = await includedPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            const withoutAo = await excludedPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(withAo.data.some((v, i) => i % 4 === 0 && v + 3 < withoutAo.data[i]), 'Screenshot exports must apply transparent occlusion at the export dimensions.');
            assert(withAo.data.every((v, i) => i % 4 !== 3 || v === withoutAo.data[i]), 'Transparent SSAO exports must preserve alpha independently of shading.');
        } finally {
            assert(includedPass instanceof WebGPUImagePass && excludedPass instanceof WebGPUImagePass, 'Transparent SSAO exports must use disposable native passes.');
            await includedPass.dispose(); await excludedPass.dispose();
        }
        const tinted = { ...config, occlusion: { name: 'on' as const, params: { ...params, transparentThreshold: 0.21, color: Color(0xffffff) } } };
        const untinted = { ...tinted, occlusion: { ...tinted.occlusion, params: { ...tinted.occlusion.params, color: Color(0x000000) } } };
        native.camera.setState({ fog: 0 }, 0); native.requestDraw();
        native.setProps({ postprocessing: untinted }); native.tick(now());
        const black = await renderer.readPixels();
        native.setProps({ postprocessing: tinted }); native.tick(now());
        const colored = await renderer.readPixels();
        let strongest = -1, difference = 0;
        for (let i = 0; i < colored.array.length; i += 4) {
            const delta = colored.array[i] - black.array[i];
            if (delta > difference && colored.array[i + 3] > 150) { strongest = i; difference = delta; }
        }
        assert(strongest >= 0 && difference > 10, 'Colored transparent SSAO must tint occluded protein pixels.');
        const pixel = strongest / 4, fogPick = await renderer.pick(pixel % colored.width, Math.floor(pixel / colored.width));
        assert(fogPick, 'Colored occlusion reference pixels must retain a source atom.');
        const view = Vec4.transformMat4(Vec4(), Vec4.create(0, 0, fogPick.depth * 2 - 1, 1), Mat4.invert(Mat4(), native.camera.projection));
        const distance = Math.abs(view[2] / view[3]);
        native.camera.setState({ fog: 100 }, 0); native.requestDraw(); native.tick(now());
        const foggedColored = await renderer.readPixels();
        native.setProps({ postprocessing: untinted }); native.tick(now());
        const foggedBaseline = await renderer.readPixels();
        const t = Math.max(0, Math.min(1, (distance - native.camera.fogNear) / (native.camera.fogFar - native.camera.fogNear)));
        assert(t > 0 && t < 1, 'The colored SSAO fog reference must lie within the fog transition.');
        const expected = difference * (1 - t * t * (3 - 2 * t));
        assert(Math.abs(foggedColored.array[strongest] - foggedBaseline.array[strongest] - expected) <= 3, `Colored transparent occlusion must fade by the analytical source-depth fog smoothstep (actual ${foggedColored.array[strongest] - foggedBaseline.array[strongest]}, expected ${expected}, distance ${distance}, fog ${native.camera.fogNear}–${native.camera.fogFar}).`);
        assert(foggedColored.array.every((v, i) => i % 4 !== 3 || v === excluded.array[i]), 'Fogged transparent occlusion must preserve alpha.');
        native.camera.setState({ fog: originalFog }, 0); native.setProps({ postprocessing: config }); native.requestDraw(); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => v === excluded.array[i]), 'Disabling transparent SSAO and fog must restore the original frame exactly.');
    } finally {
        native.camera.setState({ fog: originalFog }, 0);
        await updateAlpha(1); native.setProps({ postprocessing: originalPost }); native.requestDraw(); native.tick(now());
    }
}


async function verifyDepthPyramid(context: WebGPUContext) {
    const { device } = context, width = 13, height = 9;
    const depth = device.createTexture({ size: { width, height }, format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
    const transparent = device.createTexture({ size: { width, height }, format: 'rgba32float', usage: GPUTextureUsage.COPY_DST | GPUTextureUsage.TEXTURE_BINDING });
    const pyramid = new WebGPUDepthPyramid(context);
    const input = new Float32Array(width * height * 4);
    for (let i = 0; i < width * height; i++) input.set([(i + 1) / 128, (i % 7) / 8, 0, 0], i * 4);
    device.queue.writeTexture({ texture: transparent }, input, { bytesPerRow: width * 16 }, { width, height });
    const module = device.createShaderModule({ code: `
        @vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
            let vertices = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
            return vec4f(vertices[i], 0.0, 1.0);
        }
        @fragment fn fs(@builtin(position) p: vec4f) -> @builtin(frag_depth) f32 { return (floor(p.x) + 2.0 * floor(p.y) + 1.0) / 64.0; }
    ` });
    const pipeline = device.createRenderPipeline({ layout: 'auto', vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'fs', targets: [] }, depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'always' } });
    try {
        const encoder = device.createCommandEncoder();
        const pass = encoder.beginRenderPass({ colorAttachments: [], depthStencilAttachment: { view: depth.createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' } });
        pass.setPipeline(pipeline); pass.draw(3); pass.end(); device.queue.submit([encoder.finish()]);
        for (const [w, h] of [[13, 9], [7, 5], [1, 1], [2, 1], [13, 9]]) {
            const encoder = device.createCommandEncoder();
            const texture = pyramid.render(encoder, depth, transparent, w, h);
            device.queue.submit([encoder.finish()]);
            let previous: Float32Array | undefined, previousWidth = width, previousHeight = height;
            for (let level = 0; level < texture.mipLevelCount; level++) {
                const mw = Math.max(1, w >> level), mh = Math.max(1, h >> level);
                const bytes = await context.readTexture(texture, 0, 0, mw, mh, 16, level);
                const actual = new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4);
                const expected = new Float32Array(mw * mh * 4);
                for (let y = 0; y < mh; y++) for (let x = 0; x < mw; x++) {
                    const sx = Math.min(previousWidth - 1, Math.floor((x + 0.5) * previousWidth / mw));
                    const sy = Math.min(previousHeight - 1, Math.floor((y + 0.5) * previousHeight / mh));
                    const source = (sy * previousWidth + sx) * 4, offset = (y * mw + x) * 4;
                    if (previous) expected.set(previous.subarray(source, source + 4), offset);
                    else expected.set([(sx + 2 * sy + 1) / 64, input[source], input[source + 1], 1], offset);
                }
                assert(actual.every((v, i) => Math.abs(v - expected[i]) < 1e-7), `Native SSAO depth mip ${level} at ${w}×${h} must match analytical opaque depth and nearest transparent depth/alpha.`);
                previous = expected; previousWidth = mw; previousHeight = mh;
            }
        }
    } finally { pyramid.dispose(); depth.destroy(); transparent.destroy(); }
}


async function verifyMarking(renderer: WebGPURenderer, camera: Camera, originalProps: RendererProps) {
    const props = { ...originalProps, highlightStrength: 0, selectStrength: 0 };
    const meshProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
    const triangle = (scale: number, z: number, color: Color) => {
        const mesh = Mesh.create(new Float32Array([-4 * scale, -4 * scale, z, 4 * scale, -4 * scale, z, 0, 4 * scale, z]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
        return createRenderObject('mesh', Mesh.Utils.createValuesSimple(mesh, meshProps, color, 1), Mesh.Utils.createRenderableState(meshProps), -1);
    };
    const selected = triangle(1, 0, Color(0xff0000)), occluder = triangle(2, 2, Color(0x00ff00));
    ValueCell.update(selected.values.uMarker, 2); ValueCell.update(selected.values.markerAverage, 1);
    const marking = { ...PD.getDefaultValues(MarkingParams), selectEdgeColor: Color(0x0000ff), highlightEdgeColor: Color(0xffffff), innerEdgeFactor: 0.5 };
    const post = { ...PD.getDefaultValues(PostprocessingParams), enabled: false, antialiasing: { name: 'off' as const, params: {} } };
    const draw = async (objects = [selected], options = marking, pixelRatio = 1) => {
        renderer.render(objects, camera, props, true, pixelRatio, undefined, post, options);
        return renderer.readPixels();
    };
    const originalFog = camera.state.fog, originalViewport = { ...camera.viewport };
    try {
        camera.setState({ fog: 0 }, 0); camera.update();
        const baseline = await draw([selected], { ...marking, enabled: false });
        const marked = await draw();
        const edges: number[] = [];
        for (let i = 0; i < marked.array.length; i += 4) if (baseline.array[i + 3] === 0 && marked.array[i + 2] > 250 && marked.array[i + 3] === 255) edges.push(i);
        assert(edges.length > 20, 'Selection masks must create blue outer edges on marked primitives.');
        assert(marked.array.some((v, i) => i % 4 === 2 && baseline.array[i + 1] === 255 && Math.abs(v - 128) <= 1), 'Selection inner edges must apply the configured contrast factor.');
        assert((await draw([selected], { ...marking, selectEdgeStrength: 0 })).array.every((v, i) => v === baseline.array[i]), 'Zero selection edge strength must reproduce the unoutlined primitive exactly.');
        const thicker = await draw([selected], { ...marking, edgeScale: 3 });
        assert(thicker.array.filter((v, i) => i % 4 === 3 && baseline.array[i] === 0 && v > 0).length > edges.length * 2, 'Increasing edge scale must widen the selected silhouette.');
        const dpi = await draw([selected], marking, 3);
        assert(dpi.array.every((v, i) => v === thicker.array[i]), 'Marking edge thickness must scale with display pixel ratio.');
        const edgePixel = edges[Math.floor(edges.length / 2)] / 4;
        assert(!(await renderer.pick(edgePixel % marked.width, Math.floor(edgePixel / marked.width))), 'Marking outside geometry coverage must not invent picking IDs.');
        assert((await renderer.pick(128, 128))?.id.objectId === selected.id, 'Marking masks must preserve selected object IDs.');
        const hiddenBaseline = await draw([selected, occluder], { ...marking, enabled: false });
        assert((await draw([selected, occluder], { ...marking, ghostEdgeStrength: 0 })).array.every((v, i) => v === hiddenBaseline.array[i]), 'Disabling ghost edges must hide selections behind unmarked geometry.');
        const hidden = await draw([selected, occluder]);
        let ghostPixels = 0;
        for (let i = 0; i < hidden.array.length; i += 4) if (hiddenBaseline.array[i + 1] === 255 && Math.abs(hidden.array[i + 1] - 179) <= 1 && Math.abs(hidden.array[i + 2] - 77) <= 1) ghostPixels++;
        assert(ghostPixels > 20, 'Hidden outer edges must composite the analytical ghost opacity over their occluder.');
        const solid = await draw([selected, occluder], { ...marking, ghostEdgeStrength: 1 });
        assert(solid.array.some((v, i) => i % 4 === 2 && v === 255 && hiddenBaseline.array[i - 1] === 255), 'Full-strength ghost edges must survive without an unmarked-depth pass.');
        ValueCell.update(selected.values.uMarker, 3);
        const highlighted = await draw();
        assert(edges.every(i => highlighted.array[i] === 255 && highlighted.array[i + 1] === 255 && highlighted.array[i + 2] === 255), 'Highlight must take priority when selection and highlight bits are combined.');
        assert((await draw([selected], { ...marking, highlightEdgeStrength: 0 })).array.every((v, i) => v === baseline.array[i]), 'Highlight edge strength must be independent of selection strength.');
        camera.setState({ fog: 100 }, 0); camera.update();
        const fogged = await draw();
        const distance = 20 * camera.scale;
        const t = Math.max(0, Math.min(1, (distance - camera.fogNear) / (camera.fogFar - camera.fogNear)));
        const expected = 255 * (1 - t * t * (3 - 2 * t));
        assert(t > 0 && t < 1 && edges.every(i => Math.abs(fogged.array[i + 3] - expected) <= 2), 'Marking outer-edge alpha must fade by source-depth fog.');
        camera.setState({ fog: 0 }, 0); camera.update();
        Object.assign(camera.viewport, { x: 32, y: 24, width: 160, height: 128 }); camera.update();
        const offset = await draw();
        assert(offset.array.every((v, i) => { const p = Math.floor(i / 4), x = p % offset.width, y = Math.floor(p / offset.width); return x >= 32 && x < 192 && y >= offset.height - 152 && y < offset.height - 24 || v === 0; }), 'Marking edges must stay within offset viewports.');
        Object.assign(camera.viewport, originalViewport); camera.update();
        const exportPass = new WebGPUImagePass(renderer.context, camera, () => [selected], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: post, marking });
        try {
            const image = await exportPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            exportPass.setProps({ marking: { ...marking, enabled: false } });
            const plain = await exportPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(image.data.some((v, i) => i % 4 === 3 && v > 0 && plain.data[i] === 0), 'Custom-size screenshots must include marking outside geometry coverage.');
            exportPass.setProps({ marking });
            const resized = await exportPass.getImageData(RuntimeContext.Synchronous, 17, 11);
            assert(resized.width === 17 && resized.height === 11 && resized.data.some((v, i) => i % 4 === 3 && v > 0), 'Marking masks and edge targets must resize for small screenshot exports.');
        } finally { await exportPass.dispose(); }
        ValueCell.update(selected.values.uMarker, 0); ValueCell.update(selected.values.markerAverage, 0);
        assert((await draw()).array.every((v, i) => v === baseline.array[i]), 'Clearing all markers must reproduce the original frame and skip stale edge textures.');
    } finally { camera.setState({ fog: originalFog }, 0); Object.assign(camera.viewport, originalViewport); camera.update(); }
}


async function verifyBackground(renderer: WebGPURenderer, camera: Camera, originalProps: RendererProps) {
    const props = { ...originalProps, backgroundColor: Color(0x204060) };
    const post = { ...PD.getDefaultValues(PostprocessingParams), occlusion: { name: 'off' as const, params: {} }, antialiasing: { name: 'off' as const, params: {} } };
    const viewport = { ...camera.viewport }, snapshot = camera.getSnapshot();
    const gradient = { variant: { name: 'horizontalGradient' as const, params: { topColor: Color(0xff0000), bottomColor: Color(0x0000ff), ratio: 0.5, coverage: 'viewport' as const } } };
    const draw = async (background = gradient as typeof post.background, transparent = false, objects: GraphicsRenderObject[] = []) => {
        renderer.render(objects, camera, props, transparent, 1, undefined, { ...post, background });
        return renderer.readPixels();
    };
    const check = (pixels: { width: number, array: Uint8Array }, x: number, y: number, expected: number[], message: string, tolerance = 2) => {
        const actual = Array.from(pixels.array.subarray((y * pixels.width + x) * 4, (y * pixels.width + x) * 4 + 4));
        assert(actual.every((v, i) => Math.abs(v - expected[i]) <= tolerance), `${message}: actual ${actual}, expected ${expected}.`);
    };
    const file = async (name: string, width: number, height: number, rgba: (x: number, y: number) => number[]) => {
        const canvas = document.createElement('canvas'); canvas.width = width; canvas.height = height;
        const data = new Uint8ClampedArray(width * height * 4);
        for (let y = 0; y < height; y++) for (let x = 0; x < width; x++) data.set(rgba(x, y), (y * width + x) * 4);
        canvas.getContext('2d')!.putImageData(new ImageData(data, width, height), 0, 0);
        const blob = await new Promise<Blob>(resolve => canvas.toBlob(blob => resolve(blob!), 'image/png'));
        return { asset: Asset.File(new File([blob], name)), url: canvas.toDataURL('image/png') };
    };
    try {
        camera.setState({ fog: 0 }, 0); camera.update();
        const horizontal = await draw();
        check(horizontal, 16, 16, [255 * (1 - 16.5 / horizontal.width), 0, 255 * 16.5 / horizontal.width, 255], 'Horizontal gradients must follow the analytical top/bottom interpolation');
        const radial = { variant: { name: 'radialGradient' as const, params: { centerColor: Color(0x00ff00), edgeColor: Color(0x0000ff), ratio: 0.8, coverage: 'viewport' as const } } };
        const radialPixels = await draw(radial);
        const d = Math.hypot(128.5 / radialPixels.width - 0.5, 128.5 / radialPixels.width - 0.5) + 0.3;
        check(radialPixels, 128, 128, [0, 255 * (1 - d), 255 * d, 255], 'Radial gradients must honor center, edge and ratio');
        Object.assign(camera.viewport, { x: 32, y: 24, width: 160, height: 128 }); camera.update();
        const adjusted = await draw();
        const y = adjusted.width - 24 - 128 + 16;
        check(adjusted, 64, y, [255 * (1 - 16.5 / 128), 0, 255 * 16.5 / 128, 255], 'Viewport gradients must restart at the viewport boundary');
        check(adjusted, 0, 0, [32, 64, 96, 255], 'Background composition must preserve solid color outside offset viewports');
        const canvasGradient = await draw({ ...gradient, variant: { ...gradient.variant, params: { ...gradient.variant.params, coverage: 'canvas' } } });
        check(canvasGradient, 64, y, [255 * (1 - (y + 0.5) / adjusted.width), 0, 255 * (y + 0.5) / adjusted.width, 255], 'Canvas gradients must use full-canvas coordinates');
        check(await draw(gradient, true), 0, 0, [0, 0, 0, 0], 'Transparent exports must keep pixels outside the viewport transparent');
        Object.assign(camera.viewport, viewport); camera.update();
        const imageDefaults = BackgroundParams.variant.map('image').defaultValue as Extract<typeof post.background.variant, { name: 'image' }>['params'];
        const quadrants = await file('quadrants.png', 4, 4, (x, y) => y < 2 ? x < 2 ? [255, 0, 0, 255] : [0, 255, 0, 255] : x < 2 ? [0, 0, 255, 255] : [255, 255, 255, 255]);
        const image = { variant: { name: 'image' as const, params: { ...imageDefaults, source: { name: 'file' as const, params: quadrants.asset } } } };
        await renderer.updateBackground({ ...post, background: image });
        const imagePixels = await draw(image);
        check(imagePixels, 64, 64, [255, 0, 0, 255], 'Uploaded image backgrounds must retain top-left orientation');
        check(imagePixels, 192, 192, [255, 255, 255, 255], 'Uploaded image backgrounds must retain bottom-right orientation');
        Object.assign(camera.viewport, { x: 32, y: 16, width: 160, height: 64 }); camera.update();
        check(await draw(image), 72, 192, [255, 0, 0, 255], 'Viewport image coverage must use the viewport image origin');
        check(await draw({ variant: { ...image.variant, params: { ...image.variant.params, coverage: 'canvas' } } }), 72, 192, [0, 0, 255, 255], 'Canvas image coverage must retain full-canvas image coordinates');
        Object.assign(camera.viewport, viewport); camera.update();
        const cropColors = [[255, 0, 0, 255], [0, 255, 0, 255], [0, 0, 255, 255], [255, 255, 255, 255]];
        for (const [width, height] of [[8, 4], [4, 8]]) {
            const wide = width > height;
            const stripes = await file('cropping.png', width, height, (x, y) => cropColors[Math.floor((wide ? x : y) / 2)]);
            const cropped = { variant: { ...image.variant, params: { ...image.variant.params, source: { name: 'file' as const, params: stripes.asset } } } };
            await renderer.updateBackground({ ...post, background: cropped });
            const pixels = await draw(cropped);
            check(pixels, wide ? 32 : 64, wide ? 64 : 32, [0, 255, 0, 255], 'Aspect-cover image scaling must crop the leading outer stripe');
            check(pixels, wide ? 224 : 64, wide ? 64 : 224, [0, 0, 255, 255], 'Aspect-cover image scaling must crop the trailing outer stripe');
        }

        const red = await file('red.png', 3, 5, () => [255, 0, 0, 255]);
        const tinted = { variant: { name: 'image' as const, params: { ...imageDefaults, source: { name: 'url' as const, params: red.url }, opacity: 0.25 } } };
        await renderer.updateBackground({ ...post, background: tinted });
        check(await draw(tinted, true), 64, 64, [64, 0, 0, 64], 'Image opacity must produce premultiplied transparent background color');
        check(await draw(tinted), 64, 64, [88, 48, 72, 255], 'Image opacity must composite over the configured solid background');
        const gray = { variant: { ...tinted.variant, params: { ...tinted.variant.params, opacity: 1, saturation: -1, lightness: 0.1 } } };
        check(await draw(gray), 64, 64, [80, 80, 80, 255], 'Background saturation and lightness must match the analytical luminance adjustment');
        const stripes = await file('stripes.png', 8, 8, x => x % 2 ? [255, 255, 255, 255] : [0, 0, 0, 255]);
        const blurred = { variant: { name: 'image' as const, params: { ...imageDefaults, source: { name: 'file' as const, params: stripes.asset }, blur: 1 } } };
        await renderer.updateBackground({ ...post, background: blurred });
        check(await draw(blurred), 16, 32, [128, 128, 128, 255], 'Image blur must sample the native GPU mip pyramid');
        const sharp = await draw({ variant: { ...blurred.variant, params: { ...blurred.variant.params, blur: 0 } } });
        assert(sharp.array[(32 * sharp.width + 16) * 4] < 20, 'Disabling mip blur must restore the image detail.');
        const faces = { px: (await file('px.png', 4, 4, () => [255, 0, 0, 255])).asset, nx: (await file('nx.png', 4, 4, () => [0, 255, 255, 255])).asset,
            py: (await file('py.png', 4, 4, () => [0, 255, 0, 255])).asset, ny: (await file('ny.png', 4, 4, () => [255, 255, 0, 255])).asset,
            pz: (await file('pz.png', 4, 4, () => [255, 255, 255, 255])).asset, nz: (await file('nz.png', 4, 4, () => [0, 0, 255, 255])).asset };
        const sky = { variant: { name: 'skybox' as const, params: { ...(BackgroundParams.variant.map('skybox').defaultValue as Extract<typeof post.background.variant, { name: 'skybox' }>['params']), faces: { name: 'files' as const, params: faces } } } };
        await renderer.updateBackground({ ...post, background: sky });
        check(await draw(sky), 128, 128, [0, 0, 255, 255], 'Cubemap backgrounds must sample the forward negative-Z face');
        check(await draw({ variant: { ...sky.variant, params: { ...sky.variant.params, rotation: { x: 0, y: 90, z: 0 } } } }), 128, 128, [0, 255, 255, 255], 'Cubemap rotation must select the analytical negative-X face');
        for (const [position, up, color] of [
            [[20, 0, 0], [0, 1, 0], [0, 255, 255, 255]], [[-20, 0, 0], [0, 1, 0], [255, 0, 0, 255]],
            [[0, 20, 0], [0, 0, 1], [255, 255, 0, 255]], [[0, -20, 0], [0, 0, 1], [0, 255, 0, 255]],
            [[0, 0, -20], [0, 1, 0], [255, 255, 255, 255]],
        ]) {
            camera.setState({ position: Vec3.create(position[0], position[1], position[2]), up: Vec3.create(up[0], up[1], up[2]) }, 0); camera.update();
            check(await draw(sky), 128, 128, color, 'Camera direction must select the corresponding native cube face');
        }
        camera.setState(snapshot, 0); camera.setState({ fog: 0 }, 0); camera.update();

        camera.setState({ mode: 'orthographic' }, 0); camera.update();
        check(await draw(sky), 128, 128, [0, 0, 255, 255], 'Orthographic cameras must retain perspective skybox directions');
        const urls = { px: quadrants.url, nx: quadrants.url, py: quadrants.url, ny: quadrants.url, pz: quadrants.url, nz: quadrants.url };
        const urlSky = { variant: { ...sky.variant, params: { ...sky.variant.params, faces: { name: 'urls' as const, params: urls }, blur: 1 } } };
        await renderer.updateBackground({ ...post, background: urlSky });
        check(await draw(urlSky), 128, 128, [128, 128, 128, 255], 'URL cubemap faces must build independent native mip chains for blur');
        let invalidCube = false;
        try { await renderer.updateBackground({ ...post, background: { variant: { ...sky.variant, params: { ...sky.variant.params, faces: { name: 'files', params: { ...faces, px: red.asset } } } } } }); } catch { invalidCube = true; }
        assert(invalidCube, 'Unequal or non-square cube faces must fail before GPU texture creation.');

        camera.setState(snapshot, 0); camera.setState({ fog: 0 }, 0); camera.update();
        const meshProps = { ...PD.getDefaultValues(Mesh.Params), ignoreLight: true, doubleSided: true };
        const mesh = Mesh.create(new Float32Array([-4, -4, 0, 4, -4, 0, 0, 4, 0]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array(3), 3, 1);
        const object = createRenderObject('mesh', Mesh.Utils.createValuesSimple(mesh, meshProps, Color(0xff0000), 1), Mesh.Utils.createRenderableState(meshProps), -1);
        const blue = { variant: { ...gradient.variant, params: { ...gradient.variant.params, topColor: Color(0x0000ff), bottomColor: Color(0x0000ff) } } };
        check(await draw(blue, false, [object]), 128, 128, [255, 0, 0, 255], 'Opaque geometry must cover its environment');
        const picked = await renderer.pick(128, 128);
        assert(picked?.id.objectId === object.id, 'Environment composition must preserve geometry picking IDs.');
        camera.setState({ fog: 100 }, 0); camera.update();
        check(await draw(blue, false, [object]), 128, 128, [128, 0, 128, 255], 'Fogged geometry must fade into its environment instead of the solid clear color');
        assert((await renderer.pick(128, 128))?.id.objectId === object.id, 'Environment fog must preserve source selection even while color coverage fades.');
        camera.setState({ fog: 0 }, 0); camera.update();
        ValueCell.update(object.values.alpha, 0.25);
        check(await draw(blue, false, [object]), 128, 128, [64, 0, 191, 255], 'Transparent geometry must composite over environment color without double alpha multiplication');
        const exportPass = new WebGPUImagePass(renderer.context, camera, () => [], { cameraHelper: { axes: { name: 'off', params: {} } }, renderer: props, transparentBackground: true, postprocessing: { ...post, background: tinted } });
        try {
            const image = await exportPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(image.data.every((v, i) => Math.abs(v - [255, 0, 0, 64][i % 4]) <= 1), 'Image screenshots must wait for assets and export straight alpha at custom dimensions.');
            exportPass.setProps({ postprocessing: { ...post, background: sky } });
            const cube = await exportPass.getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(cube.data.some((v, i) => i % 4 === 2 && v > 250), 'Screenshot exports must load all cubemap faces without WebGL.');
        } finally { await exportPass.dispose(); }
        const pending = renderer.updateBackground({ ...post, background: image });
        const restored = await draw(gradient); await pending;
        assert((await draw(gradient)).array.every((v, i) => v === restored.array[i]), 'A stale background asset load must not replace a newer gradient.');
        let failed = false;
        const invalid = { variant: { ...image.variant, params: { ...image.variant.params, source: { name: 'file' as const, params: Asset.File(new File(['invalid image'], 'invalid.png')) } } } };
        try { await renderer.updateBackground({ ...post, background: invalid }); } catch { failed = true; }
        assert(failed, 'Invalid image assets must reject background preparation without a GPU validation error.');
        await draw({ variant: { name: 'off', params: {} } });
    } finally { camera.setState(snapshot, 0); Object.assign(camera.viewport, viewport); camera.update(); }
}


async function verifyBackgroundAssets(plugin: PluginContext) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const post = native.props.postprocessing;
    const canvas = document.createElement('canvas'); canvas.width = 4; canvas.height = 4;
    const ctx = canvas.getContext('2d')!; ctx.fillStyle = '#00ff00'; ctx.fillRect(0, 0, 4, 4);
    const blob = await new Promise<Blob>(resolve => canvas.toBlob(blob => resolve(blob!), 'image/png'));
    const file = new File([blob], 'managed-background.png');
    const stored = Asset.File(file), asset: Asset.File = { kind: 'file', id: stored.id, name: stored.name };
    plugin.managers.asset.set(asset, file, { isStatic: true });
    const background = { variant: { name: 'image' as const, params: { ...(BackgroundParams.variant.map('image').defaultValue as Extract<typeof post.background.variant, { name: 'image' }>['params']), source: { name: 'file' as const, params: asset } } } };
    try {
        native.setProps({ postprocessing: { ...post, enabled: true, background } }); native.tick(now());
        await renderer.updateBackground({ ...post, enabled: true, background });
        // Only the background load notification dirties this frame.
        native.tick(now());
        const pixels = await renderer.readPixels();
        assert(pixels.array[0] === 0 && pixels.array[1] === 255 && pixels.array[2] === 0 && pixels.array[3] === 255, 'Asset readiness must trigger a native canvas redraw without an additional requestDraw.');
        assert(plugin.managers.asset.get(asset)?.refCount === 1, 'Canvas backgrounds must hold their shared asset reference while visible.');
        const pass = native.getImagePass({ transparentBackground: true });
        try {
            await pass.updateBackground();
            const image = await pass.getImageData(RuntimeContext.Synchronous, 192, 128);
            assert(image.data[0] === 0 && image.data[1] === 255 && image.data[2] === 0 && image.data[3] === 255, 'Native screenshots must resolve managed assets without an embedded File object.');
            assert(plugin.managers.asset.get(asset)?.refCount === 2, 'Canvas and screenshot backgrounds must share asset-manager references.');
        } finally { assert(pass instanceof WebGPUImagePass, 'Native background screenshots must expose disposal.'); await pass.dispose(); }
        assert(plugin.managers.asset.get(asset)?.refCount === 1, 'Screenshot disposal must release only its own background asset reference.');
        native.setProps({ postprocessing: post }); native.tick(now());
        assert(plugin.managers.asset.get(asset)?.refCount === 0, 'Clearing the canvas environment must release its retained image asset.');
    } finally { native.setProps({ postprocessing: post }); native.tick(now()); plugin.managers.asset.delete(asset); }
}


async function verifyVolumeBackground(plugin: PluginContext) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const originalPost = native.props.postprocessing, originalFog = native.props.cameraFog;
    const bare = { ...PD.getDefaultValues(PostprocessingParams), enabled: false, antialiasing: { name: 'off' as const, params: {} }, occlusion: { name: 'off' as const, params: {} } };
    const blue = { variant: { name: 'horizontalGradient' as const, params: { topColor: Color(0x0000ff), bottomColor: Color(0x0000ff), ratio: 0.5, coverage: 'viewport' as const } } };
    try {
        native.setProps({ postprocessing: bare, cameraFog: { name: 'off', params: {} } }); native.tick(now());
        const baseline = await renderer.readPixels(), pick = await renderer.pick(160, 120);
        assert(pick, 'Density environment references must retain a selected volume cell.');
        native.setProps({ postprocessing: { ...bare, enabled: true, background: blue } }); native.tick(now());
        const composited = await renderer.readPixels();
        for (let i = 0; i < baseline.array.length; i += 4) {
            assert(Math.abs(composited.array[i] - baseline.array[i]) <= 2 && Math.abs(composited.array[i + 1] - baseline.array[i + 1]) <= 2 && Math.abs(composited.array[i + 2] - Math.min(255, baseline.array[i + 2] + 255 - baseline.array[i + 3])) <= 2 && composited.array[i + 3] === 255, 'Every volume pixel must composite over its blue environment with premultiplied alpha.');
        }
        const selected = await renderer.pick(160, 120);
        assert(selected?.id.objectId === pick.id.objectId && selected.id.groupId === pick.id.groupId, 'Volume environment passes must preserve cell and object selection.');
        const canvas = document.createElement('canvas'); canvas.width = 1; canvas.height = 1;
        const blob = await new Promise<Blob>(resolve => canvas.toBlob(blob => resolve(blob!), 'image/png'));
        const background = { variant: { name: 'image' as const, params: { ...(BackgroundParams.variant.map('image').defaultValue as Extract<typeof bare.background.variant, { name: 'image' }>['params']), opacity: 0, source: { name: 'file' as const, params: Asset.File(new File([blob], 'empty-environment.png')) } } } };
        const post = { ...bare, enabled: true, background };
        await renderer.updateBackground(post); native.setProps({ postprocessing: post }); native.tick(now());
        assert((await renderer.readPixels()).array.every((v, i) => Math.abs(v - baseline.array[i]) <= 1), 'A zero-opacity environment must preserve density color and coverage before fog.');
        native.setProps({ cameraFog: { name: 'on', params: { intensity: 100 } } }); native.tick(now());
        const fogged = await renderer.readPixels();
        assert(fogged.array.some((v, i) => i % 4 === 3 && v + 10 < baseline.array[i]), 'Environment fog must reduce the opacity of ray-marched density samples.');
        const fogPick = await renderer.pick(160, 120);
        assert(fogPick?.id.objectId === pick.id.objectId && fogPick.id.groupId === pick.id.groupId, 'Color fog must not change volume-cell IDs or selection opacity thresholds.');
    } finally { native.setProps({ postprocessing: originalPost, cameraFog: originalFog }); native.tick(now()); }
}


async function verifyEnvironmentSsao(plugin: PluginContext, updateAlpha: (alpha: number) => Promise<unknown>) {
    const native = plugin.canvas3d!, renderer = plugin.canvas3dContext!.webgpuRenderer!;
    const originalPost = native.props.postprocessing, originalFog = native.camera.state.fog;
    const canvas = document.createElement('canvas'); canvas.width = 1; canvas.height = 1;
    const blob = await new Promise<Blob>(resolve => canvas.toBlob(blob => resolve(blob!), 'image/png'));
    const background = { variant: { name: 'image' as const, params: { ...(BackgroundParams.variant.map('image').defaultValue as Extract<typeof originalPost.background.variant, { name: 'image' }>['params']), opacity: 0, source: { name: 'file' as const, params: Asset.File(new File([blob], 'ssao-environment.png')) } } } };
    const post = { ...PD.getDefaultValues(PostprocessingParams), background, antialiasing: { name: 'off' as const, params: {} }, occlusion: { name: 'on' as const, params: { ...PD.getDefaultValues(SsaoParams), samples: 16, radius: 2, blurKernelSize: 7, transparentThreshold: 1 } } };
    const draw = async (color: Color) => { native.setProps({ postprocessing: { ...post, occlusion: { ...post.occlusion, params: { ...post.occlusion.params, color } } } }); native.requestDraw(); native.tick(now()); return renderer.readPixels(); };
    try {
        await renderer.updateBackground(post);
        for (const alpha of [1, 0.8]) {
            await updateAlpha(alpha);
            native.camera.setState({ fog: 0 }, 0);
            const black = await draw(Color(0x000000)), white = await draw(Color(0xffffff));
            let strongest = -1, difference = 0;
            for (let i = 0; i < white.array.length; i += 4) {
                const d = white.array[i] - black.array[i];
                if (d > difference && black.array[i + 3] > 150) { strongest = i; difference = d; }
            }
            assert(strongest >= 0 && difference > 10, 'Environment SSAO references must retain a visibly occluded protein pixel.');
            native.camera.setState({ fog: 100 }, 0);
            const foggedBlack = await draw(Color(0x000000)), foggedWhite = await draw(Color(0xffffff));
            const ratio = foggedBlack.array[strongest + 3] / black.array[strongest + 3];
            assert(ratio > 0 && ratio < 1, 'Environment SSAO references must lie within the geometry fog transition.');
            const expected = difference * ratio, actual = foggedWhite.array[strongest] - foggedBlack.array[strongest];
            assert(Math.abs(actual - expected) <= 3, `Colored environment SSAO must fade once with source coverage for alpha ${alpha}: actual ${actual}, expected ${expected}.`);
            assert(foggedWhite.array.every((v, i) => i % 4 !== 3 || v === foggedBlack.array[i]), 'Environment SSAO tint must preserve fogged geometry coverage.');
        }
    } finally { await updateAlpha(1); native.camera.setState({ fog: originalFog }, 0); native.setProps({ postprocessing: originalPost }); native.requestDraw(); native.tick(now()); }
}
