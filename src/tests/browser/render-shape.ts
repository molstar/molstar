/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Representation } from '../../mol-repr/representation';
import { Canvas3D, Canvas3DContext } from '../../mol-canvas3d/canvas3d';
import { lociLabel } from '../../mol-theme/label';
import { MarkerAction } from '../../mol-util/marker-action';
import { EveryLoci } from '../../mol-model/loci';
import { RuntimeContext, Progress } from '../../mol-task';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mat4, Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { Sphere } from '../../mol-geo/primitive/sphere';
import { ColorNames } from '../../mol-util/color/names';
import { Shape } from '../../mol-model/shape';
import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { AssetManager } from '../../mol-util/assets';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { ParamDefinition } from '../../mol-util/param-definition';

const parent = document.getElementById('app')!;
parent.style.width = '100%';
parent.style.height = '100%';

const canvas = document.createElement('canvas');
parent.appendChild(canvas);
resizeCanvas(canvas, parent);

const assetManager = new AssetManager();

const canvas3dContext = Canvas3DContext.fromCanvas(canvas, assetManager);
const canvas3d = Canvas3D.create(canvas3dContext);
resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
canvas3dContext.syncPixelScale();
canvas3d.requestResize();
canvas3d.animate();

const info = document.createElement('div');
info.style.position = 'absolute';
info.style.fontFamily = 'sans-serif';
info.style.fontSize = '24pt';
info.style.bottom = '20px';
info.style.right = '20px';
info.style.color = 'white';
parent.appendChild(info);

let prevReprLoci = Representation.Loci.Empty;
canvas3d.input.move.subscribe(({ x, y }) => {
    const pickingId = canvas3d.identify(Vec2.create(x, y))?.id;
    let label = '';
    if (pickingId) {
        const reprLoci = canvas3d.getLoci(pickingId);
        label = lociLabel(reprLoci.loci);
        if (!Representation.Loci.areEqual(prevReprLoci, reprLoci)) {
            canvas3d.mark(prevReprLoci, MarkerAction.RemoveHighlight);
            canvas3d.mark(reprLoci, MarkerAction.Highlight);
            prevReprLoci = reprLoci;
        }
    } else {
        canvas3d.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        prevReprLoci = Representation.Loci.Empty;
    }
    info.innerText = label;
});

canvas3d.input.resize.subscribe(() => {
    resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
    canvas3dContext.syncPixelScale();
    canvas3d.requestResize();
});

const commonData = {
    // version, increment to trigger update
    version: 0,
    // centers of spheres
    centers: [
        0, 0, 0,
        0, 3, 0,
        1, 0, 4
    ],
    // color per group
    colors: [ColorNames.tomato, ColorNames.springgreen, ColorNames.springgreen],
    // size per group
    sizes: [1, 0.5, 0.2],
    // labels per group and instance
    labels: [
        'Sphere 0, Instance A',
        'Sphere 1, Instance A',
        'Sphere 2, Instance A',
        'Sphere 0, Instance B',
        'Sphere 1, Instance B',
        'Sphere 2, Instance B',
    ],
    // transforms
    transforms: [
        Mat4.identity(),
        Mat4.fromTranslation(Mat4(), Vec3.create(3, 0, 0))
    ],
};

const meshData = {
    ...commonData,
    props: {
        sphereDetail: 3,
    },
};
type MeshData = typeof meshData
type MeshProps = ParamDefinition.Values<Mesh.Params>;

const spheresData = {
    ...commonData,
};
type SpheresData = typeof spheresData
type SpheresProps = ParamDefinition.Values<Spheres.Params>;

/**
 * Create a mesh of spheres at given centers
 * - asynchronous (using async/await)
 * - progress tracking (via `ctx.update`)
 * - re-use storage from an existing mesh if given
 */
async function getSphereMesh(ctx: RuntimeContext, data: MeshData, props: MeshProps, mesh?: Mesh) {
    console.log('getSphereMesh');
    const { centers, sizes } = data;
    const builderState = MeshBuilder.createState(centers.length / 3, centers.length / 3, mesh);
    const t = Mat4.identity();
    const v = Vec3();
    const sphere = Sphere(data.props.sphereDetail);
    builderState.currentGroup = 0;
    for (let i = 0, il = centers.length / 3; i < il; ++i) {
        // for production, calls to update should be guarded by `if (ctx.shouldUpdate)`
        await ctx.update({ current: i, max: il, message: `adding mesh sphere ${i}` });
        builderState.currentGroup = i;
        Mat4.setIdentity(t);
        Mat4.scaleUniformly(t, t, sizes[i]);
        Mat4.setTranslation(t, Vec3.fromArray(v, centers, i * 3));
        MeshBuilder.addPrimitive(builderState, t, sphere);
    }
    return MeshBuilder.getMesh(builderState);
}

async function getSpheres(ctx: RuntimeContext, data: SpheresData, props: SpheresProps, spheres?: Spheres) {
    console.log('getSpheres');
    const { centers } = data;
    const builder = SpheresBuilder.create(centers.length / 3, centers.length / 3, spheres);
    for (let i = 0, il = centers.length / 3; i < il; ++i) {
        // for production, calls to update should be guarded by `if (ctx.shouldUpdate)`
        await ctx.update({ current: i, max: il, message: `adding sphere ${i}` });
        builder.add(centers[i * 3], centers[i * 3 + 1], centers[i * 3 + 2], i);
    }
    return builder.getSpheres();
}

/**
 * Get mesh shape from `MyData` object
 */
async function getMeshShape(ctx: RuntimeContext, data: MeshData, props: MeshProps, shape?: Shape<Mesh>) {
    const currentData = shape ? shape.sourceData as MeshData : undefined;
    if (shape && currentData?.version === data.version) {
        return shape;
    }

    await ctx.update('async creation of mesh shape from myData');
    const { centers, colors, sizes, labels, transforms } = data;
    const mesh = await getSphereMesh(ctx, data, props, shape && shape.geometry);
    const groupCount = centers.length / 3;
    return Shape.create(
        'test', { ...data }, mesh,
        (groupId: number) => colors[groupId], // color: per group, same for instances
        (groupId: number) => sizes[groupId], // size: per group, same for instances
        (groupId: number, instanceId: number) => labels[instanceId * groupCount + groupId], // label: per group and instance
        transforms
    );
}

/**
 * Get spheres shape from `MyData` object
 */
async function getSpheresShape(ctx: RuntimeContext, data: SpheresData, props: SpheresProps, shape?: Shape<Spheres>) {
    const currentData = shape ? shape.sourceData as SpheresData : undefined;
    if (shape && currentData?.version === data.version) {
        return shape;
    }

    await ctx.update('async creation of spheres shape from myData');
    const { centers, colors, sizes, labels, transforms } = data;
    const mesh = await getSpheres(ctx, data, props, shape && shape.geometry);
    const groupCount = centers.length / 3;
    return Shape.create(
        'test', { ...data }, mesh,
        (groupId: number) => colors[groupId], // color: per group, same for instances
        (groupId: number) => sizes[groupId] / 3, // size: per group, same for instances
        (groupId: number, instanceId: number) => labels[instanceId * groupCount + groupId], // label: per group and instance
        transforms
    );
}

// Init ShapeRepresentation containers
const meshRepr = ShapeRepresentation(getMeshShape, Mesh.Utils);
const spheresRepr = ShapeRepresentation(getSpheresShape, Spheres.Utils);

async function changeMeshColor() {
    meshData.colors[0] = ColorNames.darkmagenta;
    // Calling `createOrUpdate` with `data` will trigger color and transform update
    await meshRepr.createOrUpdate({}, meshData).run();
}

async function changeMeshDetail() {
    meshData.version += 1; // to trigger geometry update
    meshData.props.sphereDetail = (meshData.props.sphereDetail + 1) % 4;
    await meshRepr.createOrUpdate({}, meshData).run();
    setTimeout(changeMeshDetail, 500);
}

async function toggleSpheresVisibility() {
    spheresRepr.setState({ visible: !spheresRepr.state.visible });
    setTimeout(toggleSpheresVisibility, 1000);
}

export async function init() {
    // Create shape from meshData and add to canvas3d
    await meshRepr.createOrUpdate({ alpha: 0.5 }, meshData).run((p: Progress) => console.log(Progress.format(p)));
    console.log('mesh shape', meshRepr);
    canvas3d.add(meshRepr);

    // Create shape from spheresData and add to canvas3d
    await spheresRepr.createOrUpdate({}, spheresData).run((p: Progress) => console.log(Progress.format(p)));
    console.log('spheres shape', spheresRepr);
    canvas3d.add(spheresRepr);

    canvas3d.requestCameraReset();

    // Change color after 1s
    setTimeout(changeMeshColor, 1000);

    // Start changing mesh sphereDetail after 1s
    setTimeout(changeMeshDetail, 1000);

    // Start toggling spheres visibility after 1s
    setTimeout(toggleSpheresVisibility, 1000);
}
init();
