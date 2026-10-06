/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { copyWebGPUCameraState } from '../mol-gl/webgpu/camera';
import { WebGPURayPick } from '../mol-gl/webgpu/ray-pick';
import { BehaviorSubject, Subject, Subscription } from 'rxjs';
import { PickingId } from '../mol-geo/geometry/picking';
import { RendererStats } from '../mol-gl/renderer';
import { value } from '../mol-gl/webgpu/geometry';
import { BoundaryHelper } from '../mol-math/geometry/boundary-helper';
import { Sphere3D } from '../mol-math/geometry';
import { Ray3D } from '../mol-math/geometry/primitives/ray3d';
import { Mat4, Vec2, Vec3 } from '../mol-math/linear-algebra';
import { degToRad, radToDeg } from '../mol-math/misc';
import { EmptyLoci, isEmptyLoci } from '../mol-model/loci';
import { Representation } from '../mol-repr/representation';
import { MarkerAction } from '../mol-util/marker-action';
import { now } from '../mol-util/now';
import { deepClone } from '../mol-util/object';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { produce } from '../mol-util/produce';
import { Camera } from './camera';
import { Canvas3D, Canvas3DContext, Canvas3DParams, Canvas3DProps, Canvas3DAttribs, Canvas3DCameraResetOptions, DefaultCanvas3DAttribs } from './canvas3d';
import { TrackballControls } from './controls/trackball';
import { Canvas3dInteractionHelper } from './helper/interaction-events';
import { createWebGPUCameraHelper } from './helper/camera-helper';
import { createWebGPUHandleHelper } from './helper/handle-helper';
import { createWebGPUPointerHelper } from './helper/pointer-helper';
import { AsyncPickData, PickData } from './passes/pick';
import { WebGPUImagePass } from './passes/webgpu-image';
import { WebGPUDebugRegistry } from './helper/debug-registry';

/** Shared Mol* input, camera, representations and loci, backed by native WebGPU. */
export function createWebGPUCanvas3D(ctx: Canvas3DContext, props: Partial<Canvas3DProps>, attribs: Partial<Canvas3DAttribs>): Canvas3D {
    const { webgpu, webgpuRenderer: renderer, input } = ctx;
    if (!webgpu || !renderer) throw new Error('WebGPU Canvas3D requires an initialized WebGPU context.');
    const canvas = ctx.canvas ?? webgpu.canvas;
    const p = PD.merge(Canvas3DParams, PD.getDefaultValues(Canvas3DParams), props);
    const a = { ...deepClone(DefaultCanvas3DAttribs), ...deepClone(attribs) };
    const representations = new Map<Representation.Any, Subscription>();
    const imagePasses: WebGPUImagePass[] = [];
    const boundingSphere = Sphere3D(), boundingSphereVisible = Sphere3D();
    const scene = { boundingSphere, boundingSphereVisible };
    const camera = new Camera({ position: Vec3.create(0, 0, 100), mode: p.camera.mode, fov: degToRad(p.camera.fov), fog: p.cameraFog.name === 'on' ? p.cameraFog.params.intensity : 0 }, { x: 0, y: 0, width: canvas.width || 128, height: canvas.height || 128 });
    const controls = TrackballControls.create(input, camera, scene, p.trackball, a.trackball);
    const cameraHelper = createWebGPUCameraHelper(() => input.pixelRatio, p.camera.helper);
    renderer.setCameraHelper(cameraHelper);
    const handleHelper = createWebGPUHandleHelper(() => input.pixelRatio, p.handle);
    const pointerHelper = createWebGPUPointerHelper(p.pointer);
    renderer.setWorldHelpers(handleHelper, pointerHelper);
    const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
    const commited = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
    const commitQueueSize = new BehaviorSubject(0), reprCount = new BehaviorSubject(0), resized = new BehaviorSubject(0);
    const stats: RendererStats = { programCount: 1, shaderCount: 1, attributeCount: 0, elementsCount: 0, framebufferCount: 0, renderbufferCount: 0, textureCount: 0, vertexArrayCount: 0, drawCount: 0, instanceCount: 0, instancedDrawCount: 0 };
    let dirty = true, boundsDirty = false, paused = false, disposed = false, lost = false, reset: Canvas3DCameraResetOptions | undefined;
    let hasFittedScene = false;
    let frameHandle: number | undefined;
    let lastPick: { x: number, y: number, data: PickData | undefined } | undefined;
    let lastRayPick: { ray: Ray3D, data: PickData | undefined } | undefined;
    const rayPick = new WebGPURayPick(webgpu, { handle: handleHelper, pointer: pointerHelper });
    const objects = () => [...representations.keys()].flatMap(repr => [...repr.renderObjects]);
    const debugRegistry = new WebGPUDebugRegistry(webgpu, {
        boundingSphere, boundingSphereVisible,
        has: object => objects().includes(object),
        forEach: callback => { for (const object of objects()) callback({ values: object.values }, object); },
    });
    renderer.setDebugRegistry(debugRegistry);

    function updateBounds() {
        const helper = new BoundaryHelper('98');
        const list = objects().filter(object => object.values.drawCount.ref.value > 0);
        for (const visible of [false, true]) {
            helper.reset();
            const spheres = list.filter(object => !visible || object.state.visible).map(object => object.values.boundingSphere.ref.value).filter(sphere => sphere.radius > 0);
            for (const sphere of spheres) helper.includeSphere(sphere);
            helper.finishedIncludeStep();
            for (const sphere of spheres) helper.radiusSphere(sphere);
            const target = visible ? boundingSphereVisible : boundingSphere;
            if (spheres.length) helper.getSphere(target);
            else { Vec3.set(target.center, 0, 0, 0); target.radius = 0; }
        }
        boundsDirty = false;
        // Empty geometry can be a temporary result of iso-value/floodfill edits.
        // Keep the camera clipping radius so restoring geometry restores its view.
        if (boundingSphere.radius > 0) camera.setState({ radiusMax: boundingSphere.radius * p.sceneRadiusFactor }, 0);
    }

    function requestDraw() { dirty = true; lastPick = undefined; lastRayPick = undefined; }
    function commit() {
        if (!boundsDirty && !reset) return;
        if (boundsDirty) {
            updateBounds();
            // A representation can be attached before its asynchronous geometry exists.
            // Defer the first fit until that geometry arrives, while preserving later empty edits.
            if (!p.camera.manualReset && !hasFittedScene && boundingSphereVisible.radius > 0 && !reset) reset = { durationMs: 0 };
        }
        reprCount.next(representations.size);
        if (reset) {
            const radius = boundingSphereVisible.radius;
            if (radius > 0) {
                hasFittedScene = true;
                const adjust = controls.props.autoAdjustMinMaxDistance;
                if (adjust.name === 'on') {
                    const minDistance = adjust.params.minDistanceFactor * radius + adjust.params.minDistancePadding;
                    const maxDistance = Math.max(adjust.params.maxDistanceFactor * radius, adjust.params.maxDistanceMin);
                    controls.setProps({ minDistance, maxDistance });
                }
                const focus = camera.getFocus(boundingSphereVisible.center, radius);
                const next = typeof reset.snapshot === 'function' ? reset.snapshot(scene, camera) : reset.snapshot;
                camera.setState({ ...focus, ...next, radiusMax: boundingSphere.radius * p.sceneRadiusFactor }, reset.durationMs ?? p.cameraResetDurationMs,
                    { keyframes: reset.keyframes, easing: reset.easing ?? p.cameraResetEasing, trajectory: reset.trajectory ?? p.cameraResetTrajectory });
            }
            reset = undefined;
            requestDraw();
        }
        commited.next(now());
    }

    function getLoci(id: PickingId | undefined): Representation.Loci {
        if (!id) return { loci: EmptyLoci };
        const axes = cameraHelper.getLoci(id);
        if (!isEmptyLoci(axes)) return { loci: axes, repr: Representation.Empty };
        const handle = handleHelper.getLoci(id);
        if (!isEmptyLoci(handle)) return { loci: handle, repr: Representation.Empty };
        for (const repr of representations.keys()) {
            const loci = repr.getLoci(id);
            if (!isEmptyLoci(loci)) return { loci, repr };
        }
        return { loci: EmptyLoci };
    }

    async function pick(target: Vec2 | Ray3D): Promise<PickData | undefined> {
        if (disposed || lost) return undefined;
        if (!Array.isArray(target)) {
            const ray = Ray3D.clone(target);
            const data = await rayPick.pick(ray, camera, objects(), p.renderer, p.pickPadding, renderer!.animationTime);
            if (disposed || lost) return undefined;
            lastRayPick = { ray, data };
            return data;
        }
        const x = target[0] * input.pixelRatio, y = target[1] * input.pixelRatio;
        // Snapshot matrices before the asynchronous GPU read completes.
        const source = renderer!.getPickingCamera(x, y, camera);
        const snapshot = new Camera(source.getSnapshot(), { ...source.viewport });
        copyWebGPUCameraState(snapshot, source);
        Object.assign(snapshot.viewOffset, source.viewOffset);
        Mat4.copy(snapshot.view, source.view); Mat4.copy(snapshot.projection, source.projection);
        Mat4.copy(snapshot.projectionView, source.projectionView); Mat4.copy(snapshot.inverseProjectionView, source.inverseProjectionView);
        const result = await renderer!.pick(x, y);
        if (disposed || lost) return undefined;
        const data = result ? { id: result.id, position: Vec3.scale(Vec3(), snapshot.unproject(Vec3(), Vec3.create(x, canvas!.height - y, result.depth)), 1 / snapshot.scale) } : undefined;
        if (!disposed) lastPick = { x: target[0], y: target[1], data };
        return data;
    }
    function identify(target: Vec2 | Ray3D) {
        if (disposed || lost) return undefined;
        if (!Array.isArray(target)) return lastRayPick && Vec3.equals(lastRayPick.ray.origin, target.origin) && Vec3.equals(lastRayPick.ray.direction, target.direction) ? lastRayPick.data : undefined;
        return lastPick && lastPick.x === target[0] && lastPick.y === target[1] ? lastPick.data : undefined;
    }
    function asyncIdentify(target: Vec2 | Ray3D): AsyncPickData {
        let data: 'pending' | PickData | undefined = 'pending';
        pick(target).then(result => { data = result; }).catch(error => { data = undefined; if (!disposed && !lost) console.error('WebGPU picking failed', error); });
        return { tryGet: () => data };
    }
    const interactionHelper = new Canvas3dInteractionHelper(identify, asyncIdentify, getLoci, input, camera, controls, p.interaction, pick);

    function handleResize() {
        const v = p.viewport;
        if (v.name === 'canvas') Object.assign(camera.viewport, { x: 0, y: 0, width: Math.max(1, canvas!.width), height: Math.max(1, canvas!.height) });
        else if (v.name === 'static-frame') Object.assign(camera.viewport, v.params);
        else Object.assign(camera.viewport, { x: Math.round(v.params.x * canvas!.width), y: Math.round(v.params.y * canvas!.height), width: Math.max(1, Math.round(v.params.width * canvas!.width)), height: Math.max(1, Math.round(v.params.height * canvas!.height)) });
        Object.assign(controls.viewport, camera.viewport);
        requestDraw(); resized.next(Date.now());
    }
    let markingChanged = false;
    let startTime: number | undefined;
    function tick(t: now.Timestamp, options?: { manualDraw?: boolean, updateControls?: boolean }) {
        if (disposed || lost) return;
        commit();
        renderer!.setDpoitIterations(p.dpoitIterations);
        if (startTime === undefined) startTime = t;
        renderer!.setTime((t - startTime) / 1000);
        if (options?.updateControls !== false) controls.update(t);
        camera.transition.tick(t);
        const changed = camera.update();
        const sceneChanged = dirty || changed || controls.isAnimating || (p.renderer.enableAnimation && objects().some(o => value(o.values, 'uWiggleAmplitude', 0) > 0 || value(o.values, 'uTumbleAmplitude', 0) > 0 || value(o.values, 'dWiggle', false)));
        if (!paused && !options?.manualDraw && (sceneChanged || renderer!.multiSampleNeedsFrame || renderer!.illuminationNeedsFrame)) {
            const forceOn = p.multiSample.reduceFlicker && markingChanged && !changed && !controls.isAnimating;
            if (p.camera.stereo.name === 'on' && !p.illumination.enabled) {
                renderer!.renderStereo(objects(), camera, p.camera.stereo.params, p.renderer, p.transparentBackground, input.pixelRatio, undefined, p.postprocessing, p.marking, p.multiSample, sceneChanged, forceOn);
            } else {
                renderer!.render(objects(), camera, p.renderer, p.transparentBackground, input.pixelRatio, undefined, p.postprocessing, p.marking, p.multiSample, sceneChanged, forceOn, p.illumination);
            }
            Object.assign(stats, { drawCount: renderer!.stats.drawCount, instanceCount: renderer!.stats.instanceCount, instancedDrawCount: renderer!.stats.triangleCount * 3 });
            dirty = false; markingChanged = false;
            if (api.notifyDidDraw) didDraw.next(t);
        }
        interactionHelper.tick(t);
    }
    function animate() {
        if (frameHandle !== undefined || disposed) return;
        paused = false;
        const frame = (t: number) => { frameHandle = undefined; tick(t as now.Timestamp); if (!disposed && !paused) frameHandle = requestFrame(frame); };
        frameHandle = requestFrame(frame);
    }
    const lossSub = ctx.contextLost?.subscribe(() => { lost = true; });
    const backgroundSub = renderer.backgroundChanged.subscribe(requestDraw);
    const changedSub = ctx.changed?.subscribe(() => { handleResize(); });
    const xr = { isSupported: new BehaviorSubject(false), isPresenting: new BehaviorSubject(false), requestFailed: new Subject<string>(), request: async () => { xr.requestFailed.next('WebXR is not supported by the WebGPU backend.'); }, end: async () => {} };
    const api: Canvas3D = {
        webgpu,
        add: repr => {
            const isNew = !representations.has(repr);
            if (isNew) representations.set(repr, repr.updated.subscribe(() => { if (!repr.state.syncManually) { boundsDirty = true; requestDraw(); } }));
            boundsDirty = true;
            if (isNew && representations.size === 1 && !p.camera.manualReset) reset = { durationMs: 0 };
            requestDraw();
        },
        remove: repr => { representations.get(repr)?.unsubscribe(); representations.delete(repr); if (!representations.size) hasFittedScene = false; boundsDirty = true; requestDraw(); },
        commit,
        tick,
        update: () => { boundsDirty = true; requestDraw(); },
        clear: () => { for (const sub of representations.values()) sub.unsubscribe(); representations.clear(); hasFittedScene = false; boundsDirty = true; commit(); requestDraw(); },
        syncVisibility: () => { boundsDirty = true; requestDraw(); },
        requestDraw,
        resetTime: t => { startTime = t; controls.start(t); interactionHelper.resetTime(t); },
        animate,
        pause: noDraw => { if (frameHandle !== undefined) cancelFrame(frameHandle); frameHandle = undefined; paused = !!noDraw; },
        resume: () => { paused = false; requestDraw(); },
        requestAnimationFrame: callback => requestFrame(callback as FrameRequestCallback),
        cancelAnimationFrame: handle => cancelFrame(handle),
        identify,
        asyncIdentify,
        mark: (loci, action: MarkerAction) => { const cameraChanged = cameraHelper.mark(loci.loci, action); markingChanged = handleHelper.mark(loci.loci, action) || cameraChanged || markingChanged; for (const repr of representations.keys()) if (!loci.repr || loci.repr === repr) markingChanged = repr.mark(loci.loci, action) || markingChanged; requestDraw(); },
        getLoci,
        getRepresentations: () => Array.from(representations.keys()),
        notifyDidDraw: true,
        didDraw, commited, commitQueueSize, reprCount, resized,
        handleResize,
        requestResize: handleResize,
        requestCameraReset: options => { reset = options ?? {}; requestDraw(); },
        camera, boundingSphere, boundingSphereVisible,
        setProps: (properties, noDraw) => {
            const updated = typeof properties === 'function' ? produce(deepClone(p), properties as (draft: Canvas3DProps) => void) : properties;
            if (!updated) return;
            Object.assign(p, PD.merge(Canvas3DParams, p, updated));
            if (updated.camera) camera.setState({ mode: p.camera.mode, fov: degToRad(p.camera.fov) }, 0);
            if (updated.camera?.helper) cameraHelper.setProps(p.camera.helper);
            if (updated.handle) handleHelper.setProps(p.handle);
            if (updated.pointer) pointerHelper.setProps(p.pointer);
            if (updated.cameraFog) camera.setState({ fog: p.cameraFog.name === 'on' ? p.cameraFog.params.intensity : 0 }, 0);
            if (updated.cameraClipping) camera.setState({ clipFar: p.cameraClipping.far, minNear: p.cameraClipping.minNear }, 0);
            if (updated.cameraClipping?.radius !== undefined) {
                const radius = boundingSphere.radius * p.sceneRadiusFactor * (100 - p.cameraClipping.radius) / 100;
                // The established setter leaves zero-radius clipping requests unchanged.
                if (radius > 0) camera.setState({ radius: Math.max(0.01, radius) }, 0);
            }
            if (updated.trackball) controls.setProps(p.trackball);
            interactionHelper.setProps(p.interaction);
            if (updated.viewport) handleResize();
            if (!noDraw) requestDraw();
        },
        setAttribs: attribs => { Object.assign(a, attribs); controls.setAttribs(a.trackball); },
        getImagePass: props => { const pass = new WebGPUImagePass(webgpu, camera, objects, { dpoitIterations: p.dpoitIterations, renderer: p.renderer, postprocessing: p.postprocessing, marking: p.marking, illumination: p.illumination, cameraHelper: p.camera.helper, ...props }, ctx.assetManager, { handle: handleHelper, pointer: pointerHelper, debug: debugRegistry, transparency: () => ctx.props.transparency }); imagePasses.push(pass); return pass; },
        debugRegistry,
        getRenderObjects: objects,
        get props() {
            const current = deepClone(p);
            current.trackball = deepClone(controls.props);
            current.camera.mode = camera.state.mode;
            current.camera.fov = Math.round(radToDeg(camera.state.fov));
            current.camera.helper = deepClone(cameraHelper.props);
            current.cameraFog = camera.state.fog > 0 ? { name: 'on', params: { intensity: camera.state.fog } } : { name: 'off', params: {} };
            current.cameraClipping = { far: camera.state.clipFar, minNear: camera.state.minNear,
                radius: boundingSphere.radius > 0 ? 100 - Math.round(camera.transition.target.radius / (boundingSphere.radius * p.sceneRadiusFactor) * 100) : 0 };
            return current;
        },
        get attribs() { return deepClone(a); },
        input, stats, interaction: interactionHelper.events,
        xr,
        dispose: () => {
            if (disposed) return;
            disposed = true;
            if (frameHandle !== undefined) cancelFrame(frameHandle);
            for (const sub of representations.values()) sub.unsubscribe();
            representations.clear();
            lossSub?.unsubscribe(); changedSub?.unsubscribe(); backgroundSub.unsubscribe();
            controls.dispose(); interactionHelper.dispose(); renderer.dispose();
            rayPick.dispose();
            cameraHelper.scene.clear(); handleHelper.scene.clear(); pointerHelper.scene.clear();
            debugRegistry.clear();
            for (const pass of imagePasses) void pass.dispose();
            for (const subject of [didDraw, commited, commitQueueSize, reprCount, resized, xr.isSupported, xr.isPresenting, xr.requestFailed]) subject.complete();
        },
    };
    handleResize();
    return api;
}

const requestFrame = (callback: FrameRequestCallback): number => typeof window === 'undefined' ? setImmediate(() => callback(now())) as unknown as number : window.requestAnimationFrame(callback);
const cancelFrame = (handle: number) => typeof window === 'undefined' ? clearImmediate(handle as unknown as NodeJS.Immediate) : window.cancelAnimationFrame(handle);
