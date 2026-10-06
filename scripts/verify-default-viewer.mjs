/** Production Viewer defaults and explicit legacy backend verification. Run with Bun. */
import { chromium } from '@playwright/test';
import { createServer } from 'node:http';
import { readFile, mkdir } from 'node:fs/promises';
import { resolve, extname } from 'node:path';
import { PNG } from 'pngjs';
import jpeg from 'jpeg-js';
import { inflateRawSync } from 'node:zlib';

function readZipFiles(bytes) {
    let end = bytes.length - 22;
    while (end >= 0 && bytes.readUInt32LE(end) !== 0x06054b50) end--;
    if (end < 0) throw new Error('Export archive must have a ZIP central directory.');
    const files = new Map();
    let offset = bytes.readUInt32LE(end + 16);
    for (let i = 0; i < bytes.readUInt16LE(end + 10); i++) {
        if (bytes.readUInt32LE(offset) !== 0x02014b50) throw new Error('Invalid ZIP entry.');
        const method = bytes.readUInt16LE(offset + 10), size = bytes.readUInt32LE(offset + 20), length = bytes.readUInt16LE(offset + 28);
        const name = bytes.subarray(offset + 46, offset + 46 + length).toString();
        const local = bytes.readUInt32LE(offset + 42);
        if (bytes.readUInt32LE(local) !== 0x04034b50) throw new Error('Invalid ZIP local header.');
        const start = local + 30 + bytes.readUInt16LE(local + 26) + bytes.readUInt16LE(local + 28);
        const compressed = bytes.subarray(start, start + size);
        if (method !== 0 && method !== 8) throw new Error('Unexpected export compression.');
        const data = method === 8 ? inflateRawSync(compressed) : compressed;
        if (data.length !== bytes.readUInt32LE(offset + 24)) throw new Error('Export archive entry size must match decoded bytes.');
        files.set(name, data);
        offset += 46 + length + bytes.readUInt16LE(offset + 30) + bytes.readUInt16LE(offset + 32);
    }
    return files;
}

async function verifyGeometryDownloads(page, results) {
    const directory = resolve('tmp/webgpu'); await mkdir(directory, { recursive: true });
    await page.getByRole('button', { name: 'Export Geometry', exact: true }).click();
    const controls = page.locator('.msp-transform-wrapper').filter({ has: page.locator('#exportgeometry') });
    let current = 'glTF 2.0 Binary (.glb)';
    const save = async (label, extension) => {
        if (label !== current) {
            await controls.getByRole('button', { name: current, exact: true }).click();
            await page.getByRole('button', { name: label, exact: true }).click(); current = label;
        }
        const pending = page.waitForEvent('download');
        await controls.getByRole('button', { name: 'Save', exact: true }).click();
        const file = await pending;
        if (!file.suggestedFilename().endsWith(`.${extension}`)) throw new Error('Geometry download must use the selected format.');
        const path = resolve(directory, `viewer-geometry.${extension}`);
        await file.saveAs(path); return readFile(path);
    };
    const glb = await save(current, 'glb');
    if (glb.readUInt32LE(0) !== 0x46546c67 || glb.readUInt32LE(4) !== 2 || glb.readUInt32LE(8) !== glb.length) throw new Error('GLB download must be a complete glTF 2.0 binary.');
    const jsonLength = glb.readUInt32LE(12), model = JSON.parse(glb.subarray(20, 20 + jsonLength).toString().trim());
    const binary = glb.subarray(28 + jsonLength);
    let triangles = 0;
    for (const mesh of model.meshes) for (const primitive of mesh.primitives) {
        const positions = model.accessors[primitive.attributes.POSITION], normals = model.accessors[primitive.attributes.NORMAL];
        if (positions.count <= 0 || positions.type !== 'VEC3' || positions.componentType !== 5126 || normals.count !== positions.count) throw new Error('GLB meshes must have populated position/normal accessors.');
        const view = model.bufferViews[positions.bufferView], start = (view.byteOffset || 0) + (positions.byteOffset || 0), stride = view.byteStride || 12;
        for (let i = 0; i < positions.count; i++) for (let c = 0; c < 3; c++) if (!Number.isFinite(binary.readFloatLE(start + i * stride + c * 4))) throw new Error('Exported molecular positions must be finite.');
        const indices = model.accessors[primitive.indices];
        if (!indices || indices.count % 3 !== 0) throw new Error('GLB meshes must contain complete indexed triangles.');
        triangles += indices.count / 3;
    }
    if (triangles < 100) throw new Error('Molecular geometry download must contain populated surfaces and primitives.');
    const stl = await save('Stl (.stl)', 'stl');
    if (stl.readUInt32LE(80) !== triangles || stl.length !== 84 + triangles * 50) throw new Error('STL triangle counts and packed records must agree with GLB geometry.');
    for (let i = 0; i < triangles; i++) for (let c = 0; c < 12; c++) if (!Number.isFinite(stl.readFloatLE(84 + i * 50 + c * 4))) throw new Error('STL normals and vertices must be finite.');
    const objFiles = readZipFiles(await save('Wavefront (.obj)', 'zip'));
    const obj = [...objFiles].find(([name]) => name.endsWith('.obj')), mtl = [...objFiles].find(([name]) => name.endsWith('.mtl'));
    if (!obj || !mtl || (obj[1].toString().match(/^f /gm) || []).length !== triangles || !/^newmtl /m.test(mtl[1].toString())) throw new Error('OBJ download must contain matching triangle geometry and material definitions.');
    const usdFiles = readZipFiles(await save('Universal Scene Description (.usdz)', 'usdz'));
    const usda = [...usdFiles].find(([name]) => name.endsWith('.usda'));
    if (!usda || !usda[1].toString().startsWith('#usda 1.0') || !usda[1].toString().includes('faceVertexIndices') || !usda[1].toString().includes('def Material')) throw new Error('USDZ download must contain a complete geometry/material scene.');
    await page.evaluate(() => {
        if (window.recoveryGpuErrors.length) throw new Error(window.recoveryGpuErrors.join('\n'));
        if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('Geometry downloads must not acquire a WebGL context.');
    });
    results.push('webgpu: trusted GLB/STL/OBJ/USDZ geometry downloads, independently decoded archives, finite vertices and matching triangle counts');
}

async function verifyVolumeControls(page, results) {
    await page.evaluate(async () => {
        const viewer = window.verifiedViewer, plugin = viewer.plugin;
        await plugin.clear();
        plugin.canvas3dContext.setProps({ transparency: 'wboit' });
        // The default transfer ramp is deliberately faint; include its cells in picking.
        plugin.canvas3d.setProps({ multiSample: { mode: 'off' }, renderer: { pickingAlphaThreshold: 0.01 } });
        await viewer.loadVolumeFromUrl({ url: `${location.origin}/fixtures/density.cube`, format: 'cube', isBinary: false }, [{ type: 'absolute', value: 0.25, color: 0x3377aa }]);
        window.verifiedVolumeRepresentation = () => {
            const volume = plugin.managers.volume.hierarchy.current.volumes[0];
            return volume && volume.representations[0] && volume.representations[0].cell;
        };
    });
    const wait = async iso => {
        await page.waitForFunction(value => {
            const plugin = window.verifiedViewer.plugin, cell = window.verifiedVolumeRepresentation();
            return cell && cell.status === 'ok' && cell.transform.params.type.params.isoValue.absoluteValue === value && !plugin.behaviors.state.isBusy.value && !plugin.canvas3d.camera.transition.inTransition;
        }, iso);
        return page.evaluate(async () => {
            const plugin = window.verifiedViewer.plugin, canvas = plugin.canvas3d, renderer = plugin.canvas3dContext.webgpuRenderer;
            await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
            const texture = renderer.selectionTexture, volume = plugin.managers.volume.hierarchy.current.volumes[0].cell.obj.data;
            const ids = new Uint32Array((await canvas.webgpu.readTexture(texture, 0, 0, texture.width, texture.height, 16)).buffer);
            const objects = new Set(window.verifiedVolumeRepresentation().obj.data.repr.renderObjects.map(o => o.id + 1));
            let count = 0;
            for (let i = 0; i < ids.length; i += 4) if (objects.has(ids[i])) {
                if (count === 0) {
                    const pick = await renderer.pick(i / 4 % texture.width, Math.floor(i / 4 / texture.width));
                    const loci = pick && canvas.getLoci(pick.id).loci;
                    if (!loci || !['volume-loci', 'isosurface-loci', 'cell-loci'].includes(loci.kind) || loci.volume !== volume) throw new Error('Density controls must retain picking into the displayed volume.');
                }
                count++;
            }
            return { count, pixels: Array.from((await renderer.readPixels()).array), camera: Array.from(canvas.camera.state.position), target: Array.from(canvas.camera.state.target) };
        });
    };
    const initial = await wait(0.25);
    if (initial.count < 100) throw new Error('Loaded density must produce a visible pickable native surface.');
    await page.getByRole('button', { name: 'Actions', exact: true }).click();
    const isoRow = page.locator('.msp-control-row').filter({ has: page.getByText('Iso Value', { exact: true }) }).first();
    const isoInput = isoRow.getByRole('textbox').first();
    await isoInput.fill('0.7'); await isoInput.press('Enter');
    const smaller = await wait(0.7);
    if (smaller.count < 20 || smaller.count >= initial.count * 0.8) throw new Error('Increasing the Gaussian density isovalue must shrink the displayed surface.');
    if (JSON.stringify(smaller.camera) !== JSON.stringify(initial.camera) || JSON.stringify(smaller.target) !== JSON.stringify(initial.target)) throw new Error('Isovalue controls must preserve camera position and target.');
    await isoInput.fill('0.25'); await isoInput.press('Enter');
    const restored = await wait(0.25);
    if (restored.pixels.some((v, i) => v !== initial.pixels[i])) throw new Error('Restoring the density isovalue must restore the original frame.');
    await page.getByRole('button', { name: 'Hide component', exact: true }).click();
    await page.waitForFunction(() => window.verifiedVolumeRepresentation().state.isHidden);
    const hidden = await wait(0.25);
    if (hidden.count !== 0) throw new Error('Hiding the density must remove its pixels and picking IDs.');
    await page.getByRole('button', { name: 'Show component', exact: true }).click();
    await page.waitForFunction(() => !window.verifiedVolumeRepresentation().state.isHidden);
    const shown = await wait(0.25);
    if (shown.pixels.some((v, i) => v !== initial.pixels[i])) throw new Error('Showing the density must restore its original frame.');
    await page.evaluate(() => {
        if (window.recoveryGpuErrors.length) throw new Error(window.recoveryGpuErrors.join('\n'));
        if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('Density controls must not acquire a WebGL context.');
    });
    results.push('webgpu: production CUBE loading, trusted density isovalue/visibility controls, shrinking surfaces, restored pixels, stable camera and volume picking');
    await verifyDirectVolumeControls(page, results);
}

async function verifyDirectVolumeControls(page, results) {
    await page.getByRole('button', { name: 'Isosurface', exact: true }).click();
    await page.getByRole('button', { name: 'Direct Volume', exact: true }).click();
    await page.getByRole('button', { name: 'Update', exact: true }).click();
    const snapshot = async (key, value) => {
        await page.waitForFunction(({ key, value }) => {
            const plugin = window.verifiedViewer.plugin, cell = window.verifiedVolumeRepresentation();
            return cell && cell.status === 'ok' && cell.transform.params.type.name === 'direct-volume' && (!key || cell.transform.params.type.params[key] === value) && !plugin.behaviors.state.isBusy.value && !plugin.canvas3d.camera.transition.inTransition;
        }, { key, value });
        return page.evaluate(async () => {
            const plugin = window.verifiedViewer.plugin, canvas = plugin.canvas3d, renderer = plugin.canvas3dContext.webgpuRenderer;
            await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
            const pixels = await renderer.readPixels(), picking = await renderer.readPickingPixels(), ids = picking.array;
            const cell = window.verifiedVolumeRepresentation(), objects = new Set(cell.obj.data.repr.renderObjects.map(o => o.id + 1));
            let count = 0;
            for (let i = 0; i < ids.length; i += 4) if (objects.has(ids[i])) {
                if (count === 0) {
                    const picked = await renderer.pick(i / 4 % picking.width, Math.floor(i / 4 / picking.width));
                    const loci = picked && canvas.getLoci(picked.id).loci;
                    if (!loci || loci.kind !== 'cell-loci' || loci.volume !== plugin.managers.volume.hierarchy.current.volumes[0].cell.obj.data) throw new Error('Direct-volume controls must preserve picking into the displayed density cells.');
                }
                count++;
            }
            if (window.recoveryGpuErrors.length) throw new Error(window.recoveryGpuErrors.join('\n'));
            return { count, pixels: Array.from(pixels.array), camera: Array.from(canvas.camera.state.position), target: Array.from(canvas.camera.state.target) };
        });
    };
    const initial = await snapshot();
    if (initial.count < 100) throw new Error('Trusted representation switching must create visible ray-marched density.');
    await page.locator('.msp-mapped-parameter-group').filter({ has: page.getByRole('button', { name: 'Direct Volume', exact: true }) }).last().getByRole('button', { name: 'Type Properties', exact: true }).click();
    const opacityRow = page.locator('.msp-control-row').filter({ has: page.getByText('Opacity', { exact: true }) }).first();
    const opacity = opacityRow.getByRole('textbox').first();
    await opacity.fill('0'); await opacity.press('Enter');
    const transparent = await snapshot('alpha', 0);
    if (transparent.count !== 0 || !transparent.pixels.some((v, i) => v !== initial.pixels[i])) throw new Error('Zero density opacity must remove volume colors and picking.');
    await opacity.fill('1'); await opacity.press('Enter');
    const restored = await snapshot('alpha', 1);
    if (restored.pixels.some((v, i) => v !== initial.pixels[i])) throw new Error('Restoring direct-volume opacity must recover the original frame.');
    const dataRow = page.locator('.msp-control-row').filter({ has: page.getByText('Data Type', { exact: true }) }).first();
    await page.getByRole('button', { name: 'Advanced Options', exact: true }).first().click();
    for (const [name, type] of [['Float', 'float'], ['Half Float', 'halfFloat'], ['Byte', 'byte']]) {
        await dataRow.getByRole('button').first().click(); await page.getByRole('button', { name, exact: true }).click();
        const frame = await snapshot('dataType', type);
        if (frame.count < 100) throw new Error(`The ${type} density control must preserve visible pickable ray marching.`);
        if (JSON.stringify(frame.camera) !== JSON.stringify(initial.camera) || JSON.stringify(frame.target) !== JSON.stringify(initial.target)) throw new Error('Density format controls must preserve the camera.');
    }
    const byte = await snapshot('dataType', 'byte');
    if (byte.pixels.some((v, i) => v !== initial.pixels[i])) throw new Error('Restoring byte density must restore its original frame.');
    await page.evaluate(() => {
        if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('Direct-volume controls must not acquire WebGL.');
    });
    results.push('webgpu: trusted isosurface/direct-volume switching, opacity zero/restoration, byte/float/half-float formats, stable camera, volume-cell picking and zero WebGL requests');
}

async function verifyScreenshotDownloads(page, results) {
    const directory = resolve('tmp/webgpu'); await mkdir(directory, { recursive: true });
    await page.evaluate(() => {
        const plugin = window.verifiedViewer.plugin, helper = plugin.helpers.viewportScreenshot;
        // SMAA uses color contrast to filter alpha; isolate the exposure update.
        plugin.canvas3d.setProps({ postprocessing: { antialiasing: { name: 'off', params: {} } } });
        helper.behaviors.values.next({ ...helper.values, resolution: { name: 'custom', params: { width: 192, height: 128 } }, transparent: false, format: { name: 'png', params: {} } });
    });
    const open = () => page.getByRole('button', { name: 'Screenshot / State Snapshot', exact: true }).click();
    const panel = page.locator('.msp-viewport-controls-panel');
    const download = async filename => {
        const pending = page.waitForEvent('download');
        await panel.getByRole('button', { name: 'Download', exact: true }).click();
        const file = await pending;
        if (!file.suggestedFilename().endsWith(extname(filename))) throw new Error('Screenshot download must use the selected file format.');
        await file.saveAs(resolve(directory, filename));
        return readFile(resolve(directory, filename));
    };
    await open();
    await panel.getByRole('button', { name: /^Auto-crop/ }).click();
    await page.waitForFunction(() => !window.verifiedViewer.plugin.helpers.viewportScreenshot.cropParams.auto);
    const transparency = panel.locator('.msp-control-row').filter({ has: page.getByText('Transparent', { exact: true }) });
    await transparency.getByRole('button', { name: 'Off', exact: true }).click();
    await page.waitForFunction(() => window.verifiedViewer.plugin.helpers.viewportScreenshot.values.transparent);
    const png = PNG.sync.read(await download('viewer-export.png'));
    if (png.width !== 192 || png.height !== 128 || !png.data.some((v, i) => i % 4 === 3 && v === 0) || !png.data.some((v, i) => i % 4 !== 3 && v > 30 && png.data[i - i % 4 + 3] > 128)) throw new Error('Transparent PNG download must retain custom dimensions, empty alpha and visible molecular colors.');
    await page.evaluate(() => {
        const plugin = window.verifiedViewer.plugin;
        window.screenshotExposure = plugin.canvas3d.props.renderer.exposure;
        plugin.canvas3d.setProps({ renderer: { exposure: 0 }, dpoitIterations: 3 });
    });
    await open();
    const dark = PNG.sync.read(await download('viewer-export-updated.png'));
    if (dark.data.some((v, i) => i % 4 !== 3 && v !== 0) || dark.data.some((v, i) => i % 4 === 3 && Math.abs(v - png.data[i]) > 1)) throw new Error('Reused WebGPU screenshot passes must honor updated exposure while retaining transparency.');
    await page.evaluate(() => {
        const plugin = window.verifiedViewer.plugin, helper = plugin.helpers.viewportScreenshot;
        if (helper.imagePass.props.dpoitIterations !== 3) throw new Error('Reused screenshot passes must inherit updated depth-peeling iteration counts.');
        plugin.canvas3d.setProps({ renderer: { exposure: window.screenshotExposure } });
        plugin.canvas3dContext.setProps({ transparency: 'dpoit' });
    });
    await open();
    await transparency.getByRole('button', { name: 'On', exact: true }).click();
    const format = panel.locator('.msp-control-row').filter({ has: page.getByText('Format', { exact: true }) });
    await format.getByRole('button', { name: 'PNG', exact: true }).click();
    await page.getByRole('button', { name: 'JPEG', exact: true }).click();
    await page.waitForFunction(() => !window.verifiedViewer.plugin.helpers.viewportScreenshot.values.transparent && window.verifiedViewer.plugin.helpers.viewportScreenshot.values.format.name === 'jpeg');
    const decoded = jpeg.decode(await download('viewer-export.jpg'), { useTArray: true });
    if (decoded.width !== 192 || decoded.height !== 128 || !decoded.data.some((v, i) => i % 4 !== 3 && v < 180)) throw new Error('JPEG download must retain the configured size and visible molecular geometry.');
    await page.evaluate(() => {
        if (window.recoveryGpuErrors.length) throw new Error(window.recoveryGpuErrors.join('\n'));
        if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('Screenshot controls must not acquire a WebGL context.');
    });
    results.push('webgpu: trusted transparent PNG and opaque JPEG downloads, decoded pixels, custom sizes and live renderer/depth-peeling settings');
}

async function verifyTrajectory(page, results) {
    await page.evaluate(async () => {
        const viewer = window.verifiedViewer, plugin = viewer.plugin;
        await plugin.clear();
        plugin.canvas3d.setProps({ multiSample: { mode: 'off' } });
        await viewer.loadStructureFromUrl(`${location.origin}/fixtures/trajectory.pdb`, 'pdb', false);
        window.trajectoryModel = () => [...plugin.state.data.cells.values()].find(cell => cell.transform.transformer.definition.name === 'model-from-trajectory');
        if (!window.trajectoryModel() || window.trajectoryModel().transform.params.modelIndex !== 0) throw new Error('Trajectory must start at its first model.');
    });
    const wait = async index => {
        await page.waitForFunction(i => window.trajectoryModel().transform.params.modelIndex === i && !window.verifiedViewer.plugin.behaviors.state.isBusy.value && !window.verifiedViewer.plugin.canvas3d.camera.transition.inTransition, index);
        return page.evaluate(async () => {
            const plugin = window.verifiedViewer.plugin, canvas = plugin.canvas3d;
            await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
            const pixels = await plugin.canvas3dContext.webgpuRenderer.readPixels();
            let picked;
            for (let y = 0; y < pixels.height && !picked; y += 12) for (let x = 0; x < pixels.width && !picked; x += 12) {
                const id = await plugin.canvas3dContext.webgpuRenderer.pick(x, y);
                if (!id) continue;
                const loci = canvas.getLoci(id.id).loci;
                if (loci.kind === 'element-loci') picked = loci.elements[0].unit.model.modelNum;
            }
            const index = window.trajectoryModel().transform.params.modelIndex;
            if (picked !== index + 1) throw new Error(`Trajectory picking must resolve to the rendered model: index=${index}, model=${picked}.`);
            return { pixels: Array.from(pixels.array), camera: JSON.stringify(canvas.camera.getSnapshot()), ref: window.trajectoryModel().transform.ref };
        });
    };
    const first = await wait(0);
    await page.getByRole('button', { name: 'Next Model', exact: true }).click();
    const next = await wait(1);
    const firstCamera = JSON.parse(first.camera), nextCamera = JSON.parse(next.camera);
    // Translating float32 coordinates can change computed bounds by a few ulps.
    for (const key of ['radius', 'radiusMax']) {
        if (Math.abs(firstCamera[key] - nextCamera[key]) > 1e-6) throw new Error('Trajectory advancement must retain camera clipping for rigidly translated frames.');
        delete firstCamera[key]; delete nextCamera[key];
    }
    if (next.ref !== first.ref || JSON.stringify(nextCamera) !== JSON.stringify(firstCamera)) throw new Error('Trajectory advancement must retain model state references and camera orientation/position.');
    if (next.pixels.filter((v, i) => v !== first.pixels[i]).length < 100) throw new Error('Next Model must update visible molecular geometry.');
    await page.getByRole('button', { name: 'Previous Model', exact: true }).click();
    const previous = await wait(0);
    if (previous.pixels.some((v, i) => v !== first.pixels[i])) throw new Error('Previous Model must restore the original rendered frame.');
    await page.getByRole('button', { name: 'Previous Model', exact: true }).click();
    await wait(2);
    await page.getByRole('button', { name: 'Next Model', exact: true }).click();
    await wait(0);
    await page.getByRole('button', { name: 'Next Model', exact: true }).click();
    await wait(1);
    await page.getByRole('button', { name: 'First Model', exact: true }).click();
    await wait(0);
    results.push('webgpu: trusted trajectory next/previous/first controls, wraparound, changing geometry, stable camera and per-model picking');
    await page.evaluate(() => {
        const plugin = window.verifiedViewer.plugin, animation = plugin.managers.animation;
        animation.updateParams({ current: 'built-in.animate-model-index' });
        animation.updateCurrentParams({ mode: { name: 'loop', params: { direction: 'forward' } }, duration: { name: 'sequential', params: { maxFps: 5 } } });
        window.trajectoryAnimatedFrames = new Set();
        window.trajectoryAnimationSubscription = plugin.state.data.events.changed.subscribe(() => window.trajectoryAnimatedFrames.add(window.trajectoryModel().transform.params.modelIndex));
    });
    await page.getByRole('button', { name: 'Select Animation', exact: true }).click();
    await page.locator('.msp-animation-viewport-controls-select').getByRole('button', { name: 'Start', exact: true }).click();
    await page.waitForFunction(() => window.trajectoryAnimatedFrames.size === 3 && window.verifiedViewer.plugin.managers.animation.state.animationState === 'playing');
    await page.evaluate(async () => {
        const plugin = window.verifiedViewer.plugin, device = plugin.canvas3d.webgpu.device;
        const ref = window.trajectoryModel().transform.ref, time = plugin.animationLoop.time;
        device.destroy(); await device.lost; await plugin.webgpuRecovery;
        if (plugin.canvas3d.webgpu.device === device || plugin.canvas3d.webgl) throw new Error('Trajectory recovery must create a fresh WebGPU device.');
        if (window.trajectoryModel().transform.ref !== ref || !plugin.animationLoop.isAnimating || plugin.animationLoop.time < time || plugin.managers.animation.state.animationState !== 'playing') throw new Error('Trajectory recovery must retain model references and active playback.');
        window.recoveryGpuErrors = []; plugin.canvas3d.webgpu.errors.subscribe(error => window.recoveryGpuErrors.push(error.message));
        window.trajectoryAnimatedFrames.clear();
    });
    await page.waitForFunction(() => window.trajectoryAnimatedFrames.size === 3 && window.verifiedViewer.plugin.managers.animation.state.animationState === 'playing');
    await page.locator('.msp-animation-viewport-controls').getByRole('button', { name: 'Stop', exact: true }).click();
    await page.waitForFunction(() => window.verifiedViewer.plugin.managers.animation.state.animationState === 'stopped' && !window.verifiedViewer.plugin.behaviors.state.isBusy.value);
    await page.evaluate(() => window.trajectoryAnimationSubscription.unsubscribe());
    const stopped = await page.evaluate(() => window.trajectoryModel().transform.params.modelIndex);
    await wait(stopped);
    // Observe completed redraws to verify that stopping leaves the model stationary.
    await page.evaluate(async () => { for (let i = 0; i < 20; i++) await new Promise(resolve => requestAnimationFrame(resolve)); });
    if (await page.evaluate(() => window.trajectoryModel().transform.params.modelIndex) !== stopped) throw new Error('Stop must prevent further trajectory advancement.');
    await page.getByRole('button', { name: 'First Model', exact: true }).click();
    await wait(0);
    await page.evaluate(() => {
        if (window.recoveryGpuErrors.length) throw new Error(window.recoveryGpuErrors.join('\n'));
        if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('Trajectory playback must not acquire a WebGL context.');
    });
    results.push('webgpu: trusted trajectory animation start/stop, all frames, device-loss playback recovery, stationary stopped model and zero WebGL requests');
}

async function main() {
const root = resolve('build/viewer');
const html = (await readFile(`${root}/index.html`, 'utf8')).replace('}).then(viewer => {', '}).then(viewer => { window.verifiedViewer = viewer;');
const protein = await readFile('examples/1crn.cif');
const atoms = (await readFile('examples/trajectory/protein.pdb', 'utf8')).split('\n').filter(line => line.startsWith('ATOM'));
const trajectory = [0, 2, -2].map((offset, index) => `MODEL     ${String(index + 1).padStart(4)}\n${atoms.map(line => line.slice(0, 30) + (Number(line.slice(30, 38)) + offset).toFixed(3).padStart(8) + line.slice(38)).join('\n')}\nENDMDL`).join('\n') + '\nEND\n';
const densityValues = [];
for (let x = 0; x < 32; x++) for (let y = 0; y < 32; y++) for (let z = 0; z < 32; z++) densityValues.push(Math.exp(-((x - 16) ** 2 + (y - 16) ** 2 + (z - 16) ** 2) / 50).toExponential(6));
const density = ['Gaussian density', 'WebGPU production control fixture', '1 -8 -8 -8', '32 0.5 0 0', '32 0 0.5 0', '32 0 0 0.5', '6 0 0 0 0', densityValues.join(' ')].join('\n') + '\n';
const server = createServer(async (request, response) => {
    try {
        const path = new URL(request.url, 'http://localhost').pathname;
        if (path === '/fixtures/1crn.cif') { response.end(protein); return; }
        if (path === '/fixtures/trajectory.pdb') { response.end(trajectory); return; }
        if (path === '/fixtures/density.cube') { response.end(density); return; }
        if (path === '/viewer/' || path === '/viewer/index.html') { response.setHeader('Content-Type', 'text/html'); response.end(html); return; }
        const file = resolve(root, path.replace(/^\/viewer\//, ''));
        if (!path.startsWith('/viewer/') || !file.startsWith(`${root}/`)) { response.writeHead(404); response.end(); return; }
        response.setHeader('Content-Type', { '.js': 'text/javascript', '.css': 'text/css', '.ico': 'image/x-icon' }[extname(file)] || 'application/octet-stream');
        response.end(await readFile(file));
    } catch (error) { response.writeHead(404); response.end(); }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
let browser;
try {
    browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--enable-unsafe-webgpu'] });
    const results = [];
    for (const backend of process.env.MOLSTAR_VERIFY_VOLUME_ONLY ? ['webgpu'] : ['webgpu', 'webgl']) {
        const page = await browser.newPage();
        const errors = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.addInitScript(() => {
            window.requestedContexts = [];
            const getContext = HTMLCanvasElement.prototype.getContext;
            HTMLCanvasElement.prototype.getContext = function (kind, ...args) { window.requestedContexts.push(kind); return getContext.call(this, kind, ...args); };
        });
        try {
            await page.goto(`http://127.0.0.1:${server.address().port}/viewer/${backend === 'webgl' ? '?renderer=webgl' : ''}`, { waitUntil: 'domcontentloaded' });
            await page.waitForFunction(() => !!(window.verifiedViewer && window.verifiedViewer.plugin.canvas3d));
            const result = await page.evaluate(async expected => {
                const viewer = window.verifiedViewer, plugin = viewer.plugin;
                const actual = plugin.canvas3d.webgpu ? 'webgpu' : plugin.canvas3d.webgl ? 'webgl' : 'missing';
                if (actual !== expected) throw new Error(`Production viewer selected ${actual}, expected ${expected}.`);
                const gpuErrors = []; window.verifiedGpuErrors = gpuErrors;
                if (plugin.canvas3d.webgpu) plugin.canvas3d.webgpu.errors.subscribe(error => gpuErrors.push(error.message));
                await viewer.loadStructureFromUrl(`${location.origin}/fixtures/1crn.cif`, 'mmcif', false);
                plugin.canvas3d.commit(true);
                plugin.canvas3d.requestDraw();
                const objects = plugin.canvas3d.getRenderObjects();
                if (!objects.length) throw new Error('The production viewer must create protein representations.');
                if (expected === 'webgpu') {
                    const deadline = performance.now() + 5000;
                    let visible = false;
                    while (performance.now() < deadline && !visible) {
                        await new Promise(resolve => requestAnimationFrame(resolve));
                        const pixels = await plugin.canvas3dContext.webgpuRenderer.readPixels();
                        visible = pixels.array.some((v, i) => i % 4 !== 3 && v < 200 && pixels.array[i - i % 4 + 3] > 0);
                    }
                    plugin.animationLoop.stop({ noDraw: true });
                    if (!visible) throw new Error(`The default production viewer must render visible protein pixels. GPU errors: ${gpuErrors.join('; ')}`);
                    if (plugin.canvas3d.camera.state.radius <= 1) throw new Error('The production viewer must fit its camera to the loaded protein.');
                    await plugin.canvas3d.webgpu.device.queue.onSubmittedWorkDone();
                    if (gpuErrors.length) throw new Error(gpuErrors.join('\n'));
                    if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('The default viewer must not acquire a WebGL context.');
                }
                return `${actual}: production UI initialization and protein loading${expected === 'webgpu' ? ', native GPU readback and zero WebGL context requests' : ', explicit URL override'}`;
            }, backend);
            if (backend === 'webgpu' && process.env.MOLSTAR_VERIFY_VOLUME_ONLY) {
                await page.evaluate(() => {
                    window.recoveryGpuErrors = window.verifiedGpuErrors;
                    window.verifiedViewer.plugin.animationLoop.start();
                });
                await verifyVolumeControls(page, results);
                results.push(result);
                continue;
            }
            if (backend === 'webgpu') {
                await page.evaluate(() => window.verifiedViewer.plugin.animationLoop.start());
            await page.evaluate(async () => {
                const plugin = window.verifiedViewer.plugin, device = plugin.canvas3d.webgpu.device;
                const refs = [...plugin.state.data.cells.keys()].sort(), time = plugin.animationLoop.time;
                device.destroy(); await device.lost; await plugin.webgpuRecovery;
                if (plugin.canvas3d.webgpu.device === device || plugin.canvas3d.webgl) throw new Error('Production recovery must create a fresh native device.');
                if (JSON.stringify([...plugin.state.data.cells.keys()].sort()) !== JSON.stringify(refs)) throw new Error('Production recovery must preserve the scene state.');
                if (!plugin.animationLoop.isAnimating || plugin.animationLoop.time < time) throw new Error('Recovery must resume the animation clock.');
                window.recoveryGpuErrors = []; plugin.canvas3d.webgpu.errors.subscribe(error => window.recoveryGpuErrors.push(error.message));
            });
            results.push('webgpu: forced production device loss, automatic scene recovery and animation clock continuation');
                await page.waitForFunction(() => !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.multiSampleNeedsFrame);
                const axisTarget = await page.evaluate(async () => {
                    const plugin = window.verifiedViewer.plugin, canvas = plugin.canvas3d, renderer = plugin.canvas3dContext.webgpuRenderer;
                    const texture = renderer.selectionTexture;
                    const ids = new Uint32Array((await canvas.webgpu.readTexture(texture, 0, 0, texture.width, texture.height, 16)).buffer);
                    for (let i = 0; i < ids.length; i += 4) if (ids[i] && ids[i + 2] === 1) {
                        const loci = canvas.getLoci({ objectId: ids[i] - 1, instanceId: ids[i + 1], groupId: ids[i + 2] }).loci;
                        if (loci.kind !== 'data-loci' || loci.tag !== 'camera-axes') continue;
                        const rect = plugin.canvas3dContext.canvas.getBoundingClientRect();
                        window.verifiedAxisCamera = JSON.stringify(canvas.camera.state.up);
                        return { x: rect.left + (i / 4 % texture.width + 0.5) / canvas.input.pixelRatio, y: rect.top + (Math.floor(i / 4 / texture.width) + 0.5) / canvas.input.pixelRatio };
                    }
                    throw new Error('The default production viewer must expose a pickable X camera axis.');
                });
                await page.mouse.click(axisTarget.x, axisTarget.y);
                await page.waitForFunction(() => JSON.stringify(window.verifiedViewer.plugin.canvas3d.camera.state.up) !== window.verifiedAxisCamera);
                await page.waitForFunction(() => !window.verifiedViewer.plugin.canvas3d.camera.transition.inTransition && !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.multiSampleNeedsFrame);
                results.push('webgpu: default camera axes and trusted axis click changes camera orientation');
                await page.getByRole('button', { name: 'Settings / Controls Info', exact: true }).click();
                const settingsPanel = page.locator('.msp-viewport-controls-panel');
                const cameraGroup = settingsPanel.locator('.msp-mapped-parameter-group').filter({ has: page.getByText('Camera', { exact: true }) });
                await cameraGroup.getByRole('button', { name: 'More Options', exact: true }).click();
                const stereoRow = settingsPanel.locator('.msp-control-row').filter({ has: page.getByText('Stereo', { exact: true }) });
                await stereoRow.getByRole('button', { name: 'Off', exact: true }).click();
                await page.waitForFunction(() => window.verifiedViewer.plugin.canvas3d.props.camera.stereo.name === 'on' && window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.stereoActive && !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.multiSampleNeedsFrame);
                await page.evaluate(async () => {
                    const plugin = window.verifiedViewer.plugin, canvas = plugin.canvas3d, renderer = plugin.canvas3dContext.webgpuRenderer;
                    const pixels = await renderer.readPixels(), split = Math.floor(pixels.width / 2), found = [false, false];
                    for (let y = 0; y < pixels.height; y += 12) for (let x = 0; x < pixels.width; x += 12) {
                        const eye = x < split ? 0 : 1;
                        if (found[eye]) continue;
                        const pick = await renderer.pick(x, y);
                        if (pick && canvas.getLoci(pick.id).loci.kind === 'element-loci') found[eye] = true;
                    }
                    if (!found.every(Boolean)) throw new Error(`The production stereo setting must create two pickable molecular views: found=${found}, dimensions=${pixels.width}x${pixels.height}, camera=${JSON.stringify(canvas.camera.state)}, GPU errors=${window.verifiedGpuErrors.join('; ')}.`);
                    if (window.verifiedGpuErrors.length) throw new Error(window.verifiedGpuErrors.join('\n'));
                });
                await stereoRow.getByRole('button', { name: 'On', exact: true }).click();
                await page.waitForFunction(() => window.verifiedViewer.plugin.canvas3d.props.camera.stereo.name === 'off' && !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.stereoActive && !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.multiSampleNeedsFrame);
                await settingsPanel.getByRole('button', { name: 'Settings / Controls Info', exact: true }).click();
                results.push('webgpu: trusted production stereo settings toggle, two molecular views and per-eye picking');
                await page.getByRole('button', { name: 'Illumination', exact: true }).click();
                await page.waitForFunction(() => {
                    const plugin = window.verifiedViewer.plugin, renderer = plugin.canvas3dContext.webgpuRenderer;
                    return plugin.canvas3d.props.illumination.enabled && renderer.illuminationProgress === Math.pow(2, plugin.canvas3d.props.illumination.maxIterations) && !renderer.illuminationNeedsFrame;
                }, undefined, { timeout: 60000 });
                await page.evaluate(async () => {
                    const plugin = window.verifiedViewer.plugin;
                    await plugin.canvas3d.webgpu.device.queue.onSubmittedWorkDone();
                    const pixels = await plugin.canvas3dContext.webgpuRenderer.readPixels();
                    if (!pixels.array.some((v, i) => i % 4 !== 3 && v < 200 && pixels.array[i - i % 4 + 3] > 0)) throw new Error('Illumination UI must retain visible protein pixels.');
                    if (window.verifiedGpuErrors.length) throw new Error(window.verifiedGpuErrors.join('\n'));
                    if (window.requestedContexts.some(kind => kind === 'webgl' || kind === 'webgl2' || kind === 'experimental-webgl')) throw new Error('Illumination must not request a WebGL context.');
                });
                await page.getByRole('button', { name: 'Illumination', exact: true }).click();
                await page.waitForFunction(() => !window.verifiedViewer.plugin.canvas3d.props.illumination.enabled && !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.illuminationNeedsFrame);
                results.push('webgpu: trusted production illumination toggle, default 32-iteration convergence, visible protein and zero GPU errors');
                await page.evaluate(() => window.verifiedViewer.plugin.canvas3d.setProps({ multiSample: { mode: 'on', sampleLevel: 2 } }));
                await page.getByRole('button', { name: 'Illumination', exact: true }).click();
                await page.waitForFunction(() => {
                    const plugin = window.verifiedViewer.plugin, renderer = plugin.canvas3dContext.webgpuRenderer;
                    return plugin.canvas3d.props.illumination.enabled && renderer.illuminationProgress === 32 && !renderer.illuminationNeedsFrame;
                }, undefined, { timeout: 60000 });
                await page.evaluate(async () => {
                    const plugin = window.verifiedViewer.plugin;
                    await plugin.canvas3d.webgpu.device.queue.onSubmittedWorkDone();
                    const pixels = await plugin.canvas3dContext.webgpuRenderer.readPixels();
                    const color = plugin.canvas3d.props.renderer.backgroundColor, background = [color >> 16 & 255, color >> 8 & 255, color & 255];
                    const visible = i => pixels.array[i + 3] > 0 && Math.max(...background.map((v, c) => Math.abs(v - pixels.array[i + c]))) >= 16;
                    let visiblePixels = 0;
                    for (let i = 0; i < pixels.array.length; i += 4) if (visible(i)) visiblePixels++;
                    if (visiblePixels < 100) {
                        const renderer = plugin.canvas3dContext.webgpuRenderer, gpu = plugin.canvas3d.webgpu;
                        const ids = new Uint32Array((await gpu.readTexture(renderer.selectionTexture, 0, 0, pixels.width, pixels.height, 16)).buffer);
                        let canonical = 0; for (let i = 0; i < ids.length; i += 4) if (ids[i]) canonical++;
                        const depths = new Float32Array((await gpu.readTexture(renderer.tracingInput.textures.depth, 0, 0, pixels.width, pixels.height)).buffer);
                        const hold = await gpu.readTexture(renderer.multiSample.hold, 0, 0, pixels.width, pixels.height);
                        let held = 0; for (let i = 0; i < hold.length; i += 4) if (Math.min(hold[i], hold[i + 1], hold[i + 2]) < 239) held++;
                        throw new Error(`Supersampled illumination contrast=${visiblePixels}, canonical=${canonical}, depth=${depths.filter(v => v < 1).length}, held=${held}, camera=${JSON.stringify(plugin.canvas3d.camera.state)}, quality=${JSON.stringify(renderer.illuminationQuality)}.`);
                    }
                    let picked;
                    for (let y = 0; y < pixels.height && !picked; y += 16) for (let x = 0; x < pixels.width && !picked; x += 16) {
                        if (visible((y * pixels.width + x) * 4)) {
                            const candidate = await plugin.canvas3dContext.webgpuRenderer.pick(x, y);
                            if (candidate && plugin.canvas3d.getLoci(candidate.id).loci.kind === 'element-loci') picked = candidate;
                        }
                    }
                    if (!picked || plugin.canvas3d.getLoci(picked.id).loci.kind === 'empty-loci') throw new Error('Supersampled production pixels must resolve to canonical protein geometry.');
                    if (window.verifiedGpuErrors.length) throw new Error(window.verifiedGpuErrors.join('\n'));
                });
                await page.getByRole('button', { name: 'Illumination', exact: true }).click();
                await page.waitForFunction(() => !window.verifiedViewer.plugin.canvas3d.props.illumination.enabled && !window.verifiedViewer.plugin.canvas3dContext.webgpuRenderer.illuminationNeedsFrame);
                results.push('webgpu: production supersampled illumination convergence and disable');
                await page.evaluate(() => { if (window.recoveryGpuErrors.length) throw new Error(window.recoveryGpuErrors.join('\n')); });
                await verifyTrajectory(page, results);
                await verifyGeometryDownloads(page, results);
                await verifyScreenshotDownloads(page, results);
                await verifyVolumeControls(page, results);
            }
            if (errors.length) throw new Error(errors.join('\n'));
            results.push(result);
        } finally { await page.close(); }
    }
    console.log(JSON.stringify({ results }, null, 2));
} finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
}

}
main().catch(error => { console.error(error); process.exitCode = 1; });
