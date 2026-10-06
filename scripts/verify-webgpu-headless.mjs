/** DOM-free native WebGPU verification. Install dependencies with Bun, run in Node. */
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import { execFile } from 'node:child_process';
import { promisify } from 'node:util';
import { pathToFileURL } from 'node:url';
import { create, globals } from 'webgpu';
import * as pngjs from 'pngjs';
import jpeg from 'jpeg-js';
import { HeadlessPluginContext } from '../lib/commonjs/mol-plugin/headless-plugin-context.js';
import { DefaultPluginSpec } from '../lib/commonjs/mol-plugin/spec.js';
import { RuntimeContext } from '../lib/commonjs/mol-task/index.js';
import { BackgroundParams } from '../lib/commonjs/mol-canvas3d/passes/background.js';
import { WebGPUBackground } from '../lib/commonjs/mol-gl/webgpu/background.js';
import { Camera } from '../lib/commonjs/mol-canvas3d/camera.js';
import { Asset, AssetManager } from '../lib/commonjs/mol-util/assets.js';
import { File_ } from '../lib/commonjs/mol-util/nodejs-shims.js';
import { setCanvasModule } from '../lib/commonjs/mol-geo/geometry/text/font-atlas.js';
import { MultiSampleParams } from '../lib/commonjs/mol-canvas3d/passes/multi-sample.js';
import { ParamDefinition as PD } from '../lib/commonjs/mol-util/param-definition.js';
import { GlbExporter } from '../lib/commonjs/extensions/geo-export/glb-exporter.js';
import { Box3D, Sphere3D } from '../lib/commonjs/mol-math/geometry.js';
import { Vec3, Mat4, Tensor } from '../lib/commonjs/mol-math/linear-algebra.js';
import { Mesh } from '../lib/commonjs/mol-geo/geometry/mesh/mesh.js';
import { Spheres } from '../lib/commonjs/mol-geo/geometry/spheres/spheres.js';
import { SpheresBuilder } from '../lib/commonjs/mol-geo/geometry/spheres/spheres-builder.js';
import { Cylinders } from '../lib/commonjs/mol-geo/geometry/cylinders/cylinders.js';
import { CylindersBuilder } from '../lib/commonjs/mol-geo/geometry/cylinders/cylinders-builder.js';
import { Points } from '../lib/commonjs/mol-geo/geometry/points/points.js';
import { Lines } from '../lib/commonjs/mol-geo/geometry/lines/lines.js';
import { LinesBuilder } from '../lib/commonjs/mol-geo/geometry/lines/lines-builder.js';
import { TextureMesh } from '../lib/commonjs/mol-geo/geometry/texture-mesh/texture-mesh.js';
import { createTransform } from '../lib/commonjs/mol-geo/geometry/transform-data.js';
import { createRenderObject } from '../lib/commonjs/mol-gl/render-object.js';
import { computeMarchingCubesTextureMeshWebGPU, WebGPUTextureMeshGeometry } from '../lib/commonjs/mol-gl/webgpu/texture-mesh.js';
import { WebGPUTextureData } from '../lib/commonjs/mol-gl/webgpu/texture-data.js';
import { GPUTextureUsage } from '../lib/commonjs/mol-gl/webgpu/compat.js';
import { ValueCell } from '../lib/commonjs/mol-util/value-cell.js';
import { StructureElement } from '../lib/commonjs/mol-model/structure.js';
import { OrderedSet } from '../lib/commonjs/mol-data/int.js';
import { Overpaint } from '../lib/commonjs/mol-theme/overpaint.js';
import { Transparency } from '../lib/commonjs/mol-theme/transparency.js';
import { EveryLoci } from '../lib/commonjs/mol-model/loci.js';
import { Color } from '../lib/commonjs/mol-util/color/color.js';
import { packIntToRGBArray } from '../lib/commonjs/mol-util/number-packing.js';
import { ObjExporter } from '../lib/commonjs/extensions/geo-export/obj-exporter.js';
import { UsdzExporter } from '../lib/commonjs/extensions/geo-export/usdz-exporter.js';
import { StlExporter } from '../lib/commonjs/extensions/geo-export/stl-exporter.js';
import { unzip } from '../lib/commonjs/mol-util/zip/zip.js';

async function verifyDeviceRecovery(plugin, structure) {
    const originalCamera = plugin.canvas3d.camera.getSnapshot(), originalSelection = plugin.managers.structure.selection.getSnapshot();
    const surface = await plugin.builders.structure.representation.addRepresentation(structure, { type: 'gaussian-surface', typeParams: { resolution: 1, tryUseGpu: true, smoothColors: { name: 'on', params: { resolutionFactor: 1, sampleStride: 1 } } }, color: 'element-symbol' });
    const repr = surface.obj.data.repr;
    repr.setState({ overpaint: Overpaint('every-loci', [{ loci: EveryLoci, color: Color(0x00cc88), clear: false }]), transparency: Transparency('every-loci', [{ loci: EveryLoci, value: 0.2 }]) });
    const molecule = structure.obj.data;
    const loci = StructureElement.Loci(molecule, [{ unit: molecule.units[0], indices: OrderedSet.ofSingleton(0) }]);
    plugin.managers.interactivity.lociSelects.select({ loci });
    plugin.state.data.setCurrent(surface.ref);
    const before = await plugin.getImageRaw({ width: 103, height: 77 });
    const refs = [...plugin.state.data.cells.keys()], selection = plugin.managers.structure.selection.getSnapshot(), camera = plugin.canvas3d.camera.getSnapshot();
    const oldCanvas = plugin.canvas3d, oldDevice = oldCanvas.webgpu.device;
    oldDevice.destroy(); await oldDevice.lost; await plugin.webgpuRecovery;
    assert.notEqual(plugin.canvas3d, oldCanvas); assert.notEqual(plugin.canvas3d.webgpu.device, oldDevice);
    assert.equal(structure.obj.data, molecule, 'Device recovery must retain molecular data without reparsing/reloading.');
    assert.deepEqual([...plugin.state.data.cells.keys()].sort(), refs.sort(), 'GPU representation recovery must retain state refs.');
    assert.equal(plugin.state.data.behaviors.currentObject.value.ref, surface.ref);
    assert.deepEqual(plugin.managers.structure.selection.getSnapshot(), selection, 'Recovery must retain molecular selection.');
    assert.deepEqual(plugin.canvas3d.camera.getSnapshot(), camera, 'Recovery must retain camera state.');
    assert(surface.obj.data.repr.renderObjects.some(o => o.type === 'texture-mesh' && o.values.meta.ref.value.webgpuGeometry.device === plugin.canvas3d.webgpu.device), 'Recovered surfaces must regenerate native resources on the new device.');
    const errors = []; plugin.canvas3d.webgpu.errors.subscribe(error => errors.push(error.message));
    const after = await plugin.getImageRaw({ width: 103, height: 77 });
    assert.deepEqual(plugin.canvas3d.camera.getSnapshot(), camera, 'Screenshot must preserve recovered camera state.');
    assert(Buffer.from(before.data).equals(Buffer.from(after.data)), 'Recovered molecular surfaces, spatial overlays and selection must render the original image.');
    const provider = plugin.renderer.externalModules.webgpu;
    const transparencyMode = plugin.renderer.context.props.transparency;
    const peelingIterations = plugin.canvas3d.props.dpoitIterations, imagePeelingIterations = plugin.renderer.imagePass.props.dpoitIterations;
    plugin.renderer.context.setProps({ transparency: 'dpoit' });
    plugin.canvas3d.setProps({ dpoitIterations: 3 }); plugin.renderer.imagePass.setProps({ dpoitIterations: 3 });
    plugin.animationLoop.start({ time: 1234 });
    const elapsed = plugin.animationLoop.time;
    const runningCamera = plugin.canvas3d.camera.getSnapshot();
    const runningBefore = await plugin.getImageRaw({ width: 103, height: 77 });
    try {
        plugin.renderer.externalModules.webgpu = { requestAdapter: async () => null };
        plugin.canvas3d.webgpu.device.destroy();
        await assert.rejects(plugin.recoverWebGPU(), /No headless WebGPU adapter/);
        assert.equal(structure.obj.data, molecule, 'Unavailable-device recovery must retain molecular data for retry.');
        assert.deepEqual([...plugin.state.data.cells.keys()].sort(), refs);
    } finally { plugin.renderer.externalModules.webgpu = provider; }
    await plugin.recoverWebGPU();
    assert.equal(plugin.renderer.context.props.transparency, 'dpoit');
    assert.equal(plugin.renderer.context.webgpuRenderer.transparencyMode, 'dpoit', 'Retry must retain the selected native transparency pass.');
    assert.equal(plugin.canvas3d.props.dpoitIterations, 3); assert.equal(plugin.renderer.imagePass.props.dpoitIterations, 3);
    assert(plugin.animationLoop.isAnimating && plugin.animationLoop.time >= elapsed, 'Retry must resume the prior animation clock after an unavailable adapter.');
    plugin.animationLoop.stop({ noDraw: true });
    assert.deepEqual(plugin.canvas3d.camera.getSnapshot(), runningCamera, 'Retry must preserve camera state.');
    const retried = await plugin.getImageRaw({ width: 103, height: 77 });
    assert(Buffer.from(runningBefore.data).equals(Buffer.from(retried.data)), 'Explicit retry must restore the original molecular image and spatial layers.');
    plugin.renderer.context.setProps({ transparency: transparencyMode });
    plugin.canvas3d.setProps({ dpoitIterations: peelingIterations }); plugin.renderer.imagePass.setProps({ dpoitIterations: imagePeelingIterations });
    await plugin.build().delete(surface.ref).commit();
    plugin.managers.structure.selection.setSnapshot(originalSelection);
    plugin.canvas3d.camera.setState(originalCamera, 0); plugin.canvas3d.requestCameraReset({ snapshot: originalCamera, durationMs: 0 });
    await plugin.getImageRaw();
    assert.deepEqual(errors, []);
}

async function verifyHeadlessBackgrounds(context) {
    setCanvasModule(await import('@napi-rs/canvas'));
    assert.equal(typeof createImageBitmap, 'undefined');
    const assets = new AssetManager(), background = new WebGPUBackground(context, assets), camera = new Camera({}, { x: 0, y: 0, width: 17, height: 13 });
    const source = context.device.createTexture({ size: [17, 13], format: context.format, usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT });
    const register = (name, color) => {
        const pixels = Buffer.alloc(4 * 4 * 4);
        for (let i = 0; i < pixels.length; i += 4) pixels.set([...color, 255], i);
        const encoded = pngjs.PNG.sync.write({ width: 4, height: 4, data: pixels });
        const asset = Asset.File(new File_([encoded], name)); assets.set(asset, asset.file); return asset;
    };
    const frame = async props => {
        await background.ready(props);
        const encoder = context.device.createCommandEncoder();
        const output = background.render(encoder, source, camera, props, true, Color(0));
        context.device.queue.submit([encoder.finish()]);
        const rgba = await context.readTexture(output, 0, 0, 17, 13);
        if (context.format === 'bgra8unorm') for (let i = 0; i < rgba.length; i += 4) [rgba[i], rgba[i + 2]] = [rgba[i + 2], rgba[i]];
        return rgba;
    };
    try {
        camera.update();
        const image = PD.getDefaultValues(BackgroundParams); image.variant = { name: 'image', params: { ...BackgroundParams.variant.map('image').defaultValue, source: { name: 'file', params: register('background.png', [32, 96, 160]) }, opacity: 0.5 } };
        const pixels = await frame(image);
        for (let i = 0; i < pixels.length; i += 4) assert.deepEqual(Array.from(pixels.subarray(i, i + 4)), [16, 48, 80, 128], 'Headless image backgrounds must decode and preserve premultiplied opacity.');
        image.variant.params.blur = 1;
        assert.deepEqual(await frame(image), pixels, 'Image blur must sample native mipmaps in headless exports.');
        const colors = [[255, 0, 0], [0, 255, 0], [0, 0, 255], [255, 255, 0], [255, 0, 255], [0, 255, 255]], keys = ['px', 'nx', 'py', 'ny', 'pz', 'nz'];
        const sky = PD.getDefaultValues(BackgroundParams); sky.variant = { name: 'skybox', params: { ...BackgroundParams.variant.map('skybox').defaultValue, faces: { name: 'files', params: Object.fromEntries(keys.map((key, i) => [key, register(`${key}.png`, colors[i])])) } } };
        const center = (6 * 17 + 8) * 4;
        for (let i = 0; i < keys.length; i++) {
            const direction = Vec3.create(...[i < 2 ? (i ? 1 : -1) : 0, i >= 2 && i < 4 ? (i === 3 ? 1 : -1) : 0, i >= 4 ? (i === 5 ? 1 : -1) : 0]);
            camera.setState({ position: direction, target: Vec3(), up: i >= 2 && i < 4 ? Vec3.unitZ : Vec3.unitY }, 0); camera.update();
            const rgba = await frame(sky); assert.deepEqual(Array.from(rgba.subarray(center, center + 4)), [...colors[i], 255], `Headless cubemaps must preserve ${keys[i]} face order.`);
        }
        sky.variant.params.blur = 1;
        assert.deepEqual(Array.from((await frame(sky)).subarray(center, center + 4)), [...colors[5], 255]);
        return image;
    } finally { background.dispose(); source.destroy(); assets.dispose(); }
}

async function verifyGeometryExports(context) {
    const runtime = RuntimeContext.Synchronous;
    const bounds = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(4, 2, 2));
    const owned = [];
    const cpuTexture = array => {
        const texture = new WebGPUTextureData(); owned.push(texture);
        texture.load({ array, width: 2, height: 2 }); return texture;
    };
    const grid = array => {
        const texture = context.device.createTexture({ size: [8, 2], format: 'rgba8unorm', usage: GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC | GPUTextureUsage.TEXTURE_BINDING });
        context.device.queue.writeTexture({ texture }, array, { bytesPerRow: 32 }, [8, 2]);
        const data = new WebGPUTextureData(); owned.push(data); data.loadGPU(context.device, texture); return data;
    };
    const readGlb = async (object, triangulate = false) => {
        const exporter = new GlbExporter(bounds);
        if (triangulate) exporter.setOptions({ linesAsTriangles: true, pointsAsTriangles: true });
        await exporter.add(object, undefined, runtime);
        const bytes = new Uint8Array(await (await exporter.getBlob(runtime)).arrayBuffer());
        const view = new DataView(bytes.buffer);
        const jsonLength = view.getUint32(12, true);
        const json = JSON.parse(Buffer.from(bytes.subarray(20, 20 + jsonLength)).toString().trim());
        const binary = bytes.subarray(28 + jsonLength);
        const accessorBytes = index => {
            const accessor = json.accessors[index], bufferView = json.bufferViews[accessor.bufferView];
            const offset = (bufferView.byteOffset ?? 0) + (accessor.byteOffset ?? 0);
            return binary.subarray(offset, offset + bufferView.byteLength);
        };
        return { json, accessorBytes };
    };
    try {
        const space = Tensor.Space([2, 2, 2], [0, 1, 2], Float32Array), field = space.create();
        for (let x = 0; x < 2; x++) for (let y = 0; y < 2; y++) for (let z = 0; z < 2; z++) space.set(field, x, y, z, x);
        const generated = await computeMarchingCubesTextureMeshWebGPU(runtime, context, { scalarField: Tensor.create(space, field), isoLevel: 0.5 }, Mat4.identity(), 1, Sphere3D.create(Vec3(), 2));
        try {
            const native = generated.meta.webgpuGeometry;
            assert(native instanceof WebGPUTextureMeshGeometry);
            const props = PD.getDefaultValues(TextureMesh.Params);
            const object = createRenderObject('texture-mesh', TextureMesh.Utils.createValuesSimple(generated, props, Color(0xff0000), 1), TextureMesh.Utils.createRenderableState(props), -1);
            const exported = await readGlb(object), primitive = exported.json.meshes[0].primitives[0];
            const positions = exported.accessorBytes(primitive.attributes.POSITION);
            const actual = new Float32Array(positions.buffer, positions.byteOffset, positions.byteLength / 4);
            assert.equal(actual.length, generated.vertexCount * 3);
            for (let i = 0; i < generated.vertexCount; i++) for (let a = 0; a < 3; a++) assert.equal(actual[i * 3 + a], native.geometry.vertices[i * 20 + a], 'Native GPU-generated surface exports must read the actual float textures.');
        } finally {
            generated.meta.webgpuGeometry.destroy(); generated.doubleBuffer.destroy();
        }
        const vertices = new Float32Array([0, 0, 0, 0.75, 0.25, 0.25, 0.25, 0.75, 0.25]);
        const normals = new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]);
        const transforms = new Float32Array(32);
        transforms.set(Mat4.identity()); transforms.set(Mat4.fromTranslation(Mat4(), Vec3.create(2, 0, 0)), 16);
        const transform = createTransform(transforms, 2, undefined, 0, 0);
        const mesh = Mesh.create(vertices, new Uint32Array([0, 1, 2]), normals, new Float32Array([0, 1, 2]), 3, 1);
        const meshProps = PD.getDefaultValues(Mesh.Params);
        const objects = [createRenderObject('mesh', Mesh.Utils.createValuesSimple(mesh, meshProps, Color(0xff0000), 1, transform), Mesh.Utils.createRenderableState(meshProps), -1)];
        const rgbaPositions = new Float32Array(16), rgbaNormals = new Float32Array(16), rgbaGroups = new Uint8Array(16);
        for (let i = 0; i < 3; i++) {
            rgbaPositions.set(vertices.subarray(i * 3, i * 3 + 3), i * 4);
            rgbaNormals.set(normals.subarray(i * 3, i * 3 + 3), i * 4);
            packIntToRGBArray(i, rgbaGroups, i * 4);
        }
        // Padding is deliberately outside the geometry and must never be exported.
        rgbaPositions.set([99, 99, 99, 1], 12);
        for (const groups of [rgbaGroups, Float32Array.from(rgbaGroups, v => v / 255)]) {
            const geometry = TextureMesh.create(3, 3, cpuTexture(rgbaPositions), cpuTexture(groups), cpuTexture(rgbaNormals), Sphere3D.create(Vec3(), 2));
            const props = PD.getDefaultValues(TextureMesh.Params);
            objects.push(createRenderObject('texture-mesh', TextureMesh.Utils.createValuesSimple(geometry, props, Color(0xff0000), 1, transform), TextureMesh.Utils.createRenderableState(props), -1));
        }
        const colorArray = new Uint8Array(64), overpaintArray = new Uint8Array(64), transparencyArray = new Uint8Array(64);
        for (let z = 0; z < 2; z++) for (let y = 0; y < 2; y++) for (let x = 0; x < 4; x++) {
            const offset = (y * 8 + z * 4 + x) * 4;
            colorArray.set(x < 2 ? [255, 0, 0, 255] : [0, 0, 255, 255], offset);
            overpaintArray.set([0, 255, 0, 255], offset); transparencyArray.set([0, 0, 0, 128], offset);
        }
        const colorGrid = grid(colorArray), overpaintGrid = grid(overpaintArray), transparencyGrid = grid(transparencyArray);
        const rampColors = new Uint8Array(64), rampOverpaint = new Uint8Array(64), rampTransparency = new Uint8Array(64);
        for (let z = 0; z < 2; z++) for (let y = 0; y < 2; y++) for (let x = 0; x < 4; x++) {
            const offset = (y * 8 + z * 4 + x) * 4;
            rampColors.set([32 + x * 48, 64 + y * 64, 96 + z * 32, 255], offset);
            rampOverpaint.set([0, 128 + x * 32, 0, 255], offset);
            rampTransparency.set([0, 0, 0, 32 + x * 48 + y * 32 + z * 16], offset);
        }
        const rampColorGrid = grid(rampColors), rampOverpaintGrid = grid(rampOverpaint), rampTransparencyGrid = grid(rampTransparency);
        const linearByte = byte => {
            const srgb = Math.floor(byte) / 255;
            return Math.floor((srgb <= 0.04045 ? srgb / 12.92 : Math.pow((srgb + 0.055) / 1.055, 2.4)) * 255);
        };
        const verifyNativePrimitive = async (object, layered) => {
            if (object.type !== 'points' && object.type !== 'lines') return;
            const exported = await readGlb(object);
            assert.equal(exported.json.meshes.length, 2);
            const expectedPositions = object.type === 'points' ? [0.5, 0.5, 0.5] : [0.25, 0.5, 0.5, 0.75, 0.5, 0.5];
            for (const [instance, mesh] of exported.json.meshes.entries()) {
                const primitive = mesh.primitives[0];
                assert.equal(primitive.mode, object.type === 'points' ? 0 : 1, 'GLB must preserve compact point/line modes when conversion is disabled.');
                const bytes = exported.accessorBytes(primitive.attributes.POSITION), colors = exported.accessorBytes(primitive.attributes.COLOR_0);
                const positions = new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4);
                assert.deepEqual(Array.from(positions), expectedPositions, 'Compact GLB primitives must retain exactly their original positions.');
                for (let i = 0; i < positions.length / 3; i++) {
                    const x = positions[i * 3] + instance * 2, y = positions[i * 3 + 1], z = positions[i * 3 + 2];
                    const reference = layered ? [0, linearByte(128 + x * 32), 0, 255 - Math.floor(32 + x * 48 + y * 32 + z * 16)] : [linearByte(32 + x * 48), linearByte(64 + y * 64), linearByte(96 + z * 32), 255];
                    for (let c = 0; c < 4; c++) assert(Math.abs(colors[i * 4 + c] - reference[c]) <= 2, 'Compact point/line GLB colors and opacity must sample original vertices with the correct instance stride.');
                }
            }
        };
        assert.deepEqual((await colorGrid.readData()).array, colorArray, 'Native grid readback must discard 256-byte row padding.');
        const sphereBuilder = SpheresBuilder.create(1, 1); sphereBuilder.add(0.5, 0.5, 0.5, 0);
        const cylinderBuilder = CylindersBuilder.create(1, 1); cylinderBuilder.add(0.25, 0.5, 0.5, 0.75, 0.5, 0.5, 1, true, true, 2, 0);
        const sphereProps = PD.getDefaultValues(Spheres.Params), cylinderProps = PD.getDefaultValues(Cylinders.Params);
        const pointProps = PD.getDefaultValues(Points.Params), lineProps = PD.getDefaultValues(Lines.Params);
        const lineBuilder = LinesBuilder.create(1, 1); lineBuilder.add(0.25, 0.5, 0.5, 0.75, 0.5, 0.5, 0);
        const primitives = [
            createRenderObject('spheres', Spheres.Utils.createValuesSimple(sphereBuilder.getSpheres(), sphereProps, Color(0xff0000), 0.25, transform), Spheres.Utils.createRenderableState(sphereProps), -1),
            createRenderObject('cylinders', Cylinders.Utils.createValuesSimple(cylinderBuilder.getCylinders(), cylinderProps, Color(0xff0000), 0.1, transform), Cylinders.Utils.createRenderableState(cylinderProps), -1),
            createRenderObject('points', Points.Utils.createValuesSimple(Points.create(new Float32Array([0.5, 0.5, 0.5]), new Float32Array([0]), 1), pointProps, Color(0xff0000), 6, transform), Points.Utils.createRenderableState(pointProps), -1),
            createRenderObject('lines', Lines.Utils.createValuesSimple(lineBuilder.getLines(), lineProps, Color(0xff0000), 6, transform), Lines.Utils.createRenderableState(lineProps), -1),
        ];
        const readPrimitiveGlb = object => readGlb(object, true);
        for (const object of primitives) {
            const values = object.values;
            for (const [name, texture] of [['Color', colorGrid], ['Overpaint', overpaintGrid], ['Transparency', transparencyGrid]]) {
                ValueCell.update(values[`t${name}Grid`], texture);
                ValueCell.update(values[`u${name}GridDim`], [4, 2, 2]); ValueCell.update(values[`u${name}TexDim`], [8, 2]);
                ValueCell.update(values[`u${name}GridTransform`], [0, 0, 0, 1]);
            }
            ValueCell.update(values.dColorType, 'volumeInstance');
            const colored = await readPrimitiveGlb(object);
            assert.equal(colored.json.meshes.length, 2, `${object.type} exports must retain both instances.`);
            let triangleCount = 0, vertexCount = 0;
            for (const [i, exportedMesh] of colored.json.meshes.entries()) {
                const primitive = exportedMesh.primitives[0];
                const count = colored.json.accessors[primitive.attributes.POSITION].count;
                vertexCount += count;
                triangleCount += colored.json.accessors[primitive.indices].count / 3;
                const bytes = colored.accessorBytes(primitive.attributes.POSITION);
                const coordinates = new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4);
                assert(coordinates.every(Number.isFinite), `Instanced ${object.type} must export finite generated geometry.`);
                assert.equal(count, colored.json.accessors[colored.json.meshes[0].primitives[0].attributes.POSITION].count, 'Each primitive instance must contain one copy of its geometry.');
                assert.deepEqual(Array.from(colored.accessorBytes(primitive.attributes.COLOR_0)), Array(count).fill(i ? [0, 0, 255, 255] : [255, 0, 0, 255]).flat(), `${object.type} spatial colors must use generated vertices and the correct instance transform.`);
            }
            ValueCell.update(values.tColorGrid, rampColorGrid);
            const gradient = await readPrimitiveGlb(object);
            for (const [instance, exportedMesh] of gradient.json.meshes.entries()) {
                const primitive = exportedMesh.primitives[0];
                const bytes = gradient.accessorBytes(primitive.attributes.POSITION), colors = gradient.accessorBytes(primitive.attributes.COLOR_0);
                const positions = new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4);
                for (let i = 0; i < positions.length / 3; i++) {
                    const x = positions[i * 3] + instance * 2, y = positions[i * 3 + 1], z = positions[i * 3 + 2];
                    const reference = [32 + x * 48, 64 + y * 64, 96 + z * 32];
                    for (let c = 0; c < 3; c++) assert(Math.abs(colors[i * 4 + c] - linearByte(reference[c])) <= 2, 'Generated primitive color ramps must sample each generated vertex, independently of source-primitive mappings, and export linear glTF colors.');
                }
            }
            await verifyNativePrimitive(object, false);
            ValueCell.update(values.dOverpaint, true); ValueCell.update(values.dOverpaintType, 'volumeInstance');
            ValueCell.update(values.dTransparency, true); ValueCell.update(values.dTransparencyType, 'volumeInstance');
            ValueCell.update(values.tOverpaintGrid, rampOverpaintGrid); ValueCell.update(values.tTransparencyGrid, rampTransparencyGrid);
            const gradientLayers = await readPrimitiveGlb(object);
            for (const [instance, exportedMesh] of gradientLayers.json.meshes.entries()) {
                const primitive = exportedMesh.primitives[0];
                const bytes = gradientLayers.accessorBytes(primitive.attributes.POSITION), colors = gradientLayers.accessorBytes(primitive.attributes.COLOR_0);
                const positions = new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4);
                for (let i = 0; i < positions.length / 3; i++) {
                    const x = positions[i * 3] + instance * 2, y = positions[i * 3 + 1], z = positions[i * 3 + 2];
                    assert.equal(colors[i * 4], 0); assert.equal(colors[i * 4 + 2], 0);
                    assert(Math.abs(colors[i * 4 + 1] - linearByte(128 + x * 32)) <= 2, 'Generated primitive overpaint must sample each generated vertex and instance transform.');
                    assert(Math.abs(colors[i * 4 + 3] - (255 - Math.floor(32 + x * 48 + y * 32 + z * 16))) <= 1, 'Generated primitive transparency must sample each generated vertex and instance transform.');
                }
            }
            await verifyNativePrimitive(object, true);
            ValueCell.update(values.tColorGrid, colorGrid); ValueCell.update(values.tOverpaintGrid, overpaintGrid); ValueCell.update(values.tTransparencyGrid, transparencyGrid);
            const layered = await readPrimitiveGlb(object);
            for (const exportedMesh of layered.json.meshes) {
                const primitive = exportedMesh.primitives[0], count = layered.json.accessors[primitive.attributes.POSITION].count;
                assert.deepEqual(Array.from(layered.accessorBytes(primitive.attributes.COLOR_0)), Array(count).fill([0, 255, 0, 127]).flat(), 'Generated primitive spatial overlays must preserve every exported vertex.');
                assert.equal(layered.json.materials[primitive.material].alphaMode, 'BLEND');
            }
            const obj = new ObjExporter('primitive', bounds); await obj.add(object, undefined, runtime);
            const objData = await obj.getData();
            assert.equal(objData.obj.split('\n').filter(line => line.startsWith('v ')).length, vertexCount);
            assert.equal(objData.obj.split('\n').filter(line => line.startsWith('f ')).length, triangleCount);
            assert(objData.mtl.includes('Kd 0 1 0') && objData.mtl.includes('d 0.5'), 'Generated OBJ primitives must preserve spatial color and opacity.');
            assert(!/NaN|Infinity/.test(objData.obj), `Generated ${object.type} OBJ geometry must be finite: ${objData.obj.split('\n').filter(line => /NaN|Infinity/.test(line)).slice(0, 3).join('; ')}`);
            const usd = new UsdzExporter(bounds, 2); await usd.add(object, undefined, runtime);
            const archive = await unzip(runtime, (await usd.getData(runtime)).usdz);
            const model = Buffer.from(archive['model.usda']).toString();
            assert(model.includes('inputs:diffuseColor = (0,1,0)') && model.includes('inputs:opacity = 0.5'), 'Generated USDZ primitives must preserve spatial color and opacity.');
            assert(!/NaN|Infinity/.test(model), 'Generated USDZ geometry must be finite.');
            const faceCounts = Array.from(model.matchAll(/int\[\] faceVertexCounts = \[([^\]]*)\]/g), match => match[1].split(',').filter(v => v.trim()).length);
            assert.equal(faceCounts.reduce((sum, count) => sum + count, 0), triangleCount, 'USDZ must retain every generated primitive triangle.');
            const stl = new StlExporter(bounds); await stl.add(object, undefined, runtime);
            const bytes = (await stl.getData()).stl, view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
            assert.equal(view.getUint32(80, true), triangleCount); assert.equal(bytes.length, 84 + triangleCount * 50);
            for (let i = 0; i < triangleCount; i++) for (let c = 0; c < 12; c++) assert(Number.isFinite(view.getFloat32(84 + i * 50 + c * 4, true)), 'Generated STL normals and vertices must be finite.');
        }
        for (const object of objects) {
            const values = object.values;
            ValueCell.update(values.dColorType, 'group');
            ValueCell.update(values.tColor, { array: new Uint8Array([255, 0, 0, 0, 255, 0, 0, 0, 255]), width: 3, height: 1 });
            const grouped = await readGlb(object);
            const groupedPrimitive = grouped.json.meshes[0].primitives[0];
            assert.deepEqual(Array.from(grouped.accessorBytes(groupedPrimitive.attributes.COLOR_0)), [255, 0, 0, 255, 0, 255, 0, 255, 0, 0, 255, 255], 'Packed byte and normalized float texture groups must select their original colors.');
            const exportedPositions = grouped.accessorBytes(groupedPrimitive.attributes.POSITION);
            assert.deepEqual(Array.from(new Float32Array(exportedPositions.buffer, exportedPositions.byteOffset, exportedPositions.byteLength / 4)), Array.from(vertices), 'Geometry exports must preserve vertex coordinates and discard texture padding.');
            for (const [name, texture] of [['Color', colorGrid], ['Overpaint', overpaintGrid], ['Transparency', transparencyGrid]]) {
                ValueCell.update(values[`t${name}Grid`], texture);
                ValueCell.update(values[`u${name}GridDim`], [4, 2, 2]);
                ValueCell.update(values[`u${name}TexDim`], [8, 2]);
                ValueCell.update(values[`u${name}GridTransform`], [0, 0, 0, 1]);
            }
            ValueCell.update(values.dColorType, 'volumeInstance');
            let exported = await readGlb(object);
            assert.equal(exported.json.meshes.length, 2);
            for (let i = 0; i < 2; i++) {
                const primitive = exported.json.meshes[i].primitives[0];
                assert.equal(exported.json.accessors[primitive.attributes.POSITION].count, 3);
                assert.deepEqual(Array.from(exported.accessorBytes(primitive.attributes.COLOR_0)), Array(3).fill(i ? [0, 0, 255, 255] : [255, 0, 0, 255]).flat());
            }
            ValueCell.update(values.dOverpaint, true); ValueCell.update(values.dOverpaintType, 'volumeInstance');
            ValueCell.update(values.dTransparency, true); ValueCell.update(values.dTransparencyType, 'volumeInstance');
            exported = await readGlb(object);
            for (const exportedMesh of exported.json.meshes) {
                const primitive = exportedMesh.primitives[0];
                assert.deepEqual(Array.from(exported.accessorBytes(primitive.attributes.COLOR_0)), Array(3).fill([0, 255, 0, 127]).flat());
                assert.equal(exported.json.materials[primitive.material].alphaMode, 'BLEND');
            }
            const obj = new ObjExporter('native', bounds); await obj.add(object, undefined, runtime);
            const objData = await obj.getData();
            assert.equal(objData.obj.split('\n').filter(line => line.startsWith('v ')).length, 6);
            assert(objData.mtl.includes('Kd 0 1 0') && objData.mtl.includes('d 0.5'), 'OBJ must preserve spatial color/opacity.');
            const usd = new UsdzExporter(bounds, 2); await usd.add(object, undefined, runtime);
            const archive = await unzip(runtime, (await usd.getData(runtime)).usdz);
            const model = Buffer.from(archive['model.usda']).toString();
            assert(model.includes('inputs:diffuseColor = (0,1,0)') && model.includes('inputs:opacity = 0.5'));
            const stl = new StlExporter(bounds); await stl.add(object, undefined, runtime);
            const stlData = (await stl.getData()).stl;
            assert.equal(new DataView(stlData.buffer, stlData.byteOffset, stlData.byteLength).getUint32(80, true), 2);
            assert.equal(stlData.byteLength, 184, 'STL must contain exactly two 50-byte triangle records.');
        }
    } finally {
        for (const texture of owned) texture.destroy();
    }
}


Object.assign(globalThis, globals);
assert.equal(typeof document, 'undefined');
assert.equal(typeof ImageData, 'undefined');
let gpu = create(process.platform === 'darwin' ? ['backend=metal'] : []);
let plugin;
try {
    plugin = await HeadlessPluginContext.create({ webgpu: gpu, pngjs, 'jpeg-js': jpeg }, DefaultPluginSpec(), { width: 161, height: 129 });
    await plugin.init();
    await plugin.canvas3dInitialized;
    assert(plugin.canvas3d.webgpu && !plugin.canvas3d.webgl);
    assert.equal(plugin.canvas3dContext.canvas, undefined);
    const errors = [];
    plugin.canvas3d.webgpu.errors.subscribe(error => errors.push(error.message));
    const pdb = ['ATOM      1  N   GLY A   1      -1.200   0.000   0.000  1.00 20.00           N  ', 'ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  ', 'ATOM      3  C   GLY A   1       1.200   0.000   0.000  1.00 20.00           C  ', 'ATOM      4  O   GLY A   1       2.000   0.800   0.000  1.00 20.00           O  ', 'END'].join('\n');
    const data = await plugin.builders.data.rawData({ data: pdb });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'ball-and-stick', color: 'element-symbol' });
    const raw = await plugin.getImageRaw();
    assert.equal(raw.data.length, 161 * 129 * 4);
    assert(raw.data.some((v, i) => i % 4 < 3 && v < 150), 'Headless exports must contain visible molecular pixels.');
    assert(raw.data.every((v, i) => i % 4 !== 3 || v === 255));
    await verifyDeviceRecovery(plugin, structure);
    await verifyGeometryExports(plugin.canvas3d.webgpu);
    const background = await verifyHeadlessBackgrounds(plugin.canvas3d.webgpu);
    const sampling = plugin.renderer.imagePass.props.multiSample;
    plugin.renderer.imagePass.setProps({ multiSample: { ...sampling, mode: 'off' } });
    const backdrop = await plugin.getImageRaw({ width: 103, height: 77 }, { background, occlusion: { name: 'off', params: {} }, shadow: { name: 'off', params: {} }, outline: { name: 'off', params: {} }, antialiasing: { name: 'off', params: {} }, sharpening: { name: 'off', params: {} }, dof: { name: 'off', params: {} }, bloom: { name: 'off', params: {} } });
    for (let i = 0; i < 4; i++) assert(Math.abs(backdrop.data[i] - [144, 176, 208, 255][i]) <= 1, `The actual headless molecular screenshot must composite its image background: ${Array.from(backdrop.data.subarray(0, 4))}.`);
    plugin.renderer.imagePass.setProps({ multiSample: sampling });
    const objects = plugin.canvas3d.getRenderObjects();
    assert(objects.some(o => o.type === 'spheres') && objects.some(o => o.type === 'cylinders'));
    const exporter = new GlbExporter(Box3D.fromSphere3D(Box3D(), plugin.canvas3d.boundingSphereVisible));
    for (const object of objects) await exporter.add(object, undefined, RuntimeContext.Synchronous);
    const glb = new Uint8Array(await (await exporter.getBlob(RuntimeContext.Synchronous)).arrayBuffer());
    const header = new DataView(glb.buffer);
    assert.equal(header.getUint32(0, true), 0x46546c67); assert.equal(header.getUint32(4, true), 2); assert.equal(header.getUint32(8, true), glb.length);
    const geometry = JSON.parse(Buffer.from(glb.subarray(20, 20 + header.getUint32(12, true))).toString().trim());
    assert(geometry.meshes.length > 0 && geometry.meshes.every(mesh => mesh.primitives.every(primitive => geometry.accessors[primitive.attributes.POSITION].count > 0)), 'Native spheres and cylinders must export populated GLB meshes without WebGL.');
    const png = await plugin.getImagePng();
    assert(Buffer.from(raw.data).equals(png.data), 'PNG bytes must match native raw readback.');
    const encodedJpeg = await plugin.getImageJpeg();
    const decodedJpeg = jpeg.decode(encodedJpeg.data);
    assert.equal(decodedJpeg.width, raw.width); assert.equal(decodedJpeg.height, raw.height);
    await fs.mkdir('tmp/webgpu', { recursive: true });
    await plugin.saveImage('tmp/webgpu/headless-native.png');
    assert(Buffer.from(raw.data).equals(pngjs.PNG.sync.read(await fs.readFile('tmp/webgpu/headless-native.png')).data));
    plugin.canvas3d.setProps({ transparentBackground: true, multiSample: { mode: 'off' } });
    plugin.renderer.imagePass.setProps({ transparentBackground: true, multiSample: { ...PD.getDefaultValues(MultiSampleParams), mode: 'off' } });
    const transparent = await plugin.getImageRaw({ width: 103, height: 77 });
    assert.equal(transparent.data.length, 103 * 77 * 4);
    assert(transparent.data.some((v, i) => i % 4 === 3 && v === 0));
    assert(transparent.data.some((v, i) => i % 4 === 3 && v > 0));
    const crop = await plugin.renderer.imagePass.getImageRaw(RuntimeContext.Synchronous, 103, 77, { x: 20, y: 10, width: 41, height: 33 });
    for (let y = 0; y < 33; y++) assert(Buffer.from(crop.data.subarray(y * 41 * 4, (y + 1) * 41 * 4)).equals(Buffer.from(transparent.data.subarray(((y + 10) * 103 + 20) * 4, ((y + 10) * 103 + 61) * 4))));
    for (const mode of ['smaa', 'fxaa']) {
        const frame = await plugin.getImageRaw({ width: 103, height: 77 }, { antialiasing: { name: mode, params: mode === 'smaa' ? { edgeThreshold: 0.1, maxSearchSteps: 16 } : { edgeThresholdMin: 0.0312, edgeThresholdMax: 0.125, iterations: 12, subpixelQuality: 0.3 } } });
        assert(frame.data.some((v, i) => i % 4 === 3 && v > 0 && v < 255), `${mode} must smooth transparent molecular silhouettes without DOM APIs.`);
    }
    await plugin.canvas3d.webgpu.device.queue.onSubmittedWorkDone();
    assert.deepEqual(errors, []);
    const component = { kind: 'component', params: { selector: 'polymer' }, children: [{ kind: 'representation', params: { type: 'cartoon' }, children: [{ kind: 'color', params: { color: '#ff8800' } }] }, { kind: 'label', params: { text: 'Native WebGPU' } }] };
    const mvs = { metadata: { version: '1', timestamp: '2026-10-04T00:00:00Z' }, root: { kind: 'root', children: [{ kind: 'download', params: { url: pathToFileURL(`${process.cwd()}/examples/1crn.cif`).href }, children: [{ kind: 'parse', params: { format: 'mmcif' }, children: [{ kind: 'structure', params: { type: 'model' }, children: [component] }] }] }] } };
    await fs.writeFile('tmp/webgpu/cli-label.mvsj', JSON.stringify(mvs));
    component.children.pop();
    await fs.writeFile('tmp/webgpu/cli-plain.mvsj', JSON.stringify(mvs));
    const cli = promisify(execFile);
    await cli(process.execPath, ['lib/commonjs/cli/mvs/mvs-render.js', '-i', 'tmp/webgpu/cli-label.mvsj', 'tmp/webgpu/cli-label.mvsj', 'tmp/webgpu/cli-plain.mvsj', '-o', 'tmp/webgpu/cli-label.png', 'tmp/webgpu/cli-label.jpg', 'tmp/webgpu/cli-plain.png', '--size', '161x129', '--molj'], { timeout: 60000 });
    const labeled = pngjs.PNG.sync.read(await fs.readFile('tmp/webgpu/cli-label.png'));
    const plain = pngjs.PNG.sync.read(await fs.readFile('tmp/webgpu/cli-plain.png'));
    assert.equal(labeled.width, 161); assert.equal(labeled.height, 129);
    assert(labeled.data.some((v, i) => i % 4 === 0 && v > 100 && labeled.data[i + 1] > 40 && labeled.data[i + 2] < 40), 'Default CLI rendering must show the orange protein.');
    assert(!labeled.data.equals(plain.data), 'Font-backed MVS labels must change the exported frame.');
    assert.equal(jpeg.decode(await fs.readFile('tmp/webgpu/cli-label.jpg')).width, 161);
    const state = JSON.parse(await fs.readFile('tmp/webgpu/cli-label.molj', 'utf8'));
    assert(state.entries?.[0]?.snapshot?.data?.tree?.transforms?.length > 0, 'The CLI must save a usable Mol* state alongside its image.');
    const second = structuredClone(mvs.root);
    const recolor = node => { if (node.kind === 'color') node.params.color = '#0088ff'; for (const child of node.children ?? []) recolor(child); };
    recolor(second);
    await fs.writeFile('tmp/webgpu/cli-animation.mvsj', JSON.stringify({ kind: 'multiple', metadata: mvs.metadata, snapshots: [{ root: mvs.root, metadata: { key: 'orange', duration_ms: 100 } }, { root: second, metadata: { key: 'blue', duration_ms: 100 } }] }));
    await cli(process.execPath, ['lib/commonjs/cli/mvs/mvs-render.js', '-i', 'tmp/webgpu/cli-animation.mvsj', '-o', 'tmp/webgpu/cli-animation.mp4', '--size', '161x129'], { timeout: 60000 });
    const movie = await fs.readFile('tmp/webgpu/cli-animation.mp4');
    assert(movie.length > 1000 && movie.toString('ascii', 4, 8) === 'ftyp', 'CLI snapshot animation must produce an MP4 container.');
    const probe = JSON.parse((await cli('ffprobe', ['-v', 'error', '-select_streams', 'v:0', '-show_entries', 'stream=width,height,nb_frames', '-of', 'json', 'tmp/webgpu/cli-animation.mp4'])).stdout).streams[0];
    assert.equal(probe.width, 160); assert.equal(probe.height, 128); assert.equal(Number(probe.nb_frames), 7);
    const decoded = (await cli('ffmpeg', ['-v', 'error', '-i', 'tmp/webgpu/cli-animation.mp4', '-f', 'rawvideo', '-pix_fmt', 'rgba', 'pipe:1'], { encoding: 'buffer', maxBuffer: 16 * 1024 * 1024 })).stdout;
    const frameBytes = 160 * 128 * 4;
    assert.equal(decoded.length, frameBytes * 7);
    const colors = frame => {
        let orange = 0, blue = 0;
        for (let i = frame * frameBytes; i < (frame + 1) * frameBytes; i += 4) {
            if (decoded[i] > 100 && decoded[i + 1] > 40 && decoded[i + 2] < 50) orange++;
            if (decoded[i + 2] > 100 && decoded[i + 1] > 40 && decoded[i] < 50) blue++;
        }
        return { orange, blue };
    };
    assert(colors(0).orange > 100 && colors(6).blue > 100, 'Decoded animation frames must show both native molecular snapshots.');
    console.log(JSON.stringify({ results: ['DOM-free WebGPU plugin initialization, PDB representations, raw/PNG/JPEG exports, odd/custom sizes, transparent alpha, cropping, full sampling and native SMAA/FXAA', 'native GPU-generated texture-mesh GLB export and native mesh/texture-mesh GLB/OBJ/USDZ/STL exports, exact grid samples, packed groups, per-instance spatial colors/overpaint/transparency, aligned GPU readback, instanced sphere/cylinder/point/line GLB/OBJ/USDZ/STL exports and compact GLB point/line modes, analytical generated-vertex spatial gradients, finite geometry and matching triangle counts without WebGL', 'forced native device loss and unavailable-adapter retry, animation-clock continuation, GPU surface reconstruction, preserved molecule/state refs/camera/selection/spatial layers and byte-identical restored exports', 'DOM-free image and six-face cubemap backgrounds, opacity, face orientation and mipmap blur', 'default WebGPU MVS CLI, protein and font-backed labels, PNG/JPEG output, state snapshots and independently decoded odd-size MP4 animation'], errors }, null, 2));
} catch (error) {
    console.error(error.stack); process.exitCode = 1;
} finally {
    plugin?.dispose();
    await plugin?.renderer.dispose();
    plugin = undefined; gpu = undefined;
}
