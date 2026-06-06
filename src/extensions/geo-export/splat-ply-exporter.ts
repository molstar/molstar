/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Exports `spheres`/`points` render objects (e.g. a Spacefill representation) as a standard
 * 3D Gaussian-Splatting PLY: one isotropic gaussian per atom/point, encoded with the INRIA
 * conventions used by SuperSplat / PlayCanvas / gsplat (log-scale, logit-opacity, float f_dc
 * spherical-harmonic DC color term). Other geometry kinds are skipped.
 */

import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { SpheresValues } from '../../mol-gl/renderable/spheres';
import { PointsValues } from '../../mol-gl/renderable/points';
import { BaseValues, SizeValues } from '../../mol-gl/renderable/schema';
import { TextureImage } from '../../mol-gl/renderable/util';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { RuntimeContext } from '../../mol-task';
import { Box3D } from '../../mol-math/geometry';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Color } from '../../mol-util/color/color';
import { unpackRGBToInt } from '../../mol-util/number-packing';
import { sizeDataFactor } from '../../mol-geo/geometry/size-data';
import { MeshExporter, MeshGeoData } from './mesh-exporter';

// SH DC band constant Y_0^0 = sqrt(1 / (4*pi)); 3DGS stores color as f_dc = (c - 0.5) / C0.
const C0 = 0.28209479177387814;

export type SplatPlyData = {
    ply: Uint8Array
}

const v3fromArray = Vec3.fromArray;

export class SplatPlyExporter extends MeshExporter<SplatPlyData> {
    readonly fileExtension = 'ply';

    private centerTransform: Mat4;
    private readonly px: number[] = [];
    private readonly py: number[] = [];
    private readonly pz: number[] = [];
    private readonly radii: number[] = [];
    private readonly colors: Color[] = [];
    private readonly alphas: number[] = [];

    // not used: splats are collected in `add` directly from sphere/point centers
    protected async addMeshWithColors(): Promise<void> {}

    private static decodeSize(values: BaseValues & SizeValues, instanceIndex: number, group: number, vertexIndex: number): number {
        const tSize = values.tSize.ref.value as TextureImage<Uint8Array>;
        const fromTex = (i: number) => unpackRGBToInt(tSize.array[i * 3], tSize.array[i * 3 + 1], tSize.array[i * 3 + 2]) / sizeDataFactor;
        let size = 0;
        switch (values.dSizeType.ref.value) {
            case 'uniform': size = values.uSize.ref.value; break;
            case 'instance': size = fromTex(instanceIndex); break;
            case 'group': size = fromTex(group); break;
            case 'groupInstance': size = fromTex(instanceIndex * values.uGroupCount.ref.value + group); break;
            case 'vertex': size = fromTex(vertexIndex); break;
            case 'vertexInstance': size = fromTex(instanceIndex * values.uVertexCount.ref.value + vertexIndex); break;
        }
        return size * values.uSizeFactor.ref.value;
    }

    private pushSplat(world: Vec3, radius: number, color: Color, alpha: number) {
        this.px.push(world[0]);
        this.py.push(world[1]);
        this.pz.push(world[2]);
        this.radii.push(radius);
        this.colors.push(color);
        this.alphas.push(alpha);
    }

    private addCenters(values: SpheresValues | PointsValues, centers: Float32Array, count: number, groupOf: (i: number) => number) {
        const isVolumeColor = values.dColorType.ref.value.startsWith('volume');
        const t = Mat4();
        const c = Vec3();
        const w = Vec3();
        const aTransform = values.aTransform.ref.value;
        const instanceCount = values.instanceCount.ref.value;
        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            Mat4.fromArray(t, aTransform, instanceIndex * 16);
            Mat4.mul(t, this.centerTransform, t);
            for (let i = 0; i < count; ++i) {
                v3fromArray(c, centers, i * 3);
                Vec3.transformMat4(w, c, t);
                const group = groupOf(i);
                const radius = SplatPlyExporter.decodeSize(values, instanceIndex, group, i);
                const geoData: MeshGeoData = { values, groups: undefined, vertexCount: count, instanceIndex, isGeoTexture: false, mode: 'points' };
                // volume color types need GPU-interpolated colors we don't read here; fall back to white
                const color = isVolumeColor ? Color(0xffffff) : MeshExporter.getColor(i, geoData);
                const transparency = isVolumeColor ? 0 : MeshExporter.getTransparency(i, geoData);
                this.pushSplat(w, radius, color, 1 - transparency);
            }
        }
    }

    add(renderObject: GraphicsRenderObject, webgl: WebGLContext, ctx: RuntimeContext): Promise<void> | undefined {
        if (!renderObject.state.visible && !this.options.includeHidden) return;
        if (renderObject.values.drawCount.ref.value === 0) return;
        if (renderObject.values.instanceCount.ref.value === 0) return;

        if (renderObject.type === 'spheres') {
            const values = renderObject.values as SpheresValues;
            const centers = values.centerBuffer.ref.value;
            const groups = values.groupBuffer.ref.value;
            const sphereCount = values.uVertexCount.ref.value / 6;
            this.addCenters(values, centers, sphereCount, i => groups[i]);
        } else if (renderObject.type === 'points') {
            const values = renderObject.values as PointsValues;
            const centers = values.aPosition.ref.value;
            const groups = values.aGroup.ref.value;
            const pointCount = values.uVertexCount.ref.value;
            this.addCenters(values, centers, pointCount, i => groups[i]);
        }
        // mesh / lines / cylinders / texture-mesh are not gaussian splats -> skip
        return;
    }

    async getData() {
        const n = this.px.length;
        const header =
            'ply\n' +
            'format binary_little_endian 1.0\n' +
            `element vertex ${n}\n` +
            'property float x\nproperty float y\nproperty float z\n' +
            'property float scale_0\nproperty float scale_1\nproperty float scale_2\n' +
            'property float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\n' +
            'property float opacity\n' +
            'property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n' +
            'end_header\n';
        const headerBytes = new TextEncoder().encode(header);

        const stride = 14 * 4; // 14 float32 per splat
        const body = new ArrayBuffer(n * stride);
        const dv = new DataView(body);

        for (let i = 0; i < n; ++i) {
            let o = i * stride;
            dv.setFloat32(o, this.px[i], true); o += 4;
            dv.setFloat32(o, this.py[i], true); o += 4;
            dv.setFloat32(o, this.pz[i], true); o += 4;

            // isotropic gaussian: scale stored as log(radius) (viewers apply exp)
            const logS = Math.log(Math.max(this.radii[i], 1e-6));
            dv.setFloat32(o, logS, true); o += 4;
            dv.setFloat32(o, logS, true); o += 4;
            dv.setFloat32(o, logS, true); o += 4;

            // identity quaternion, INRIA layout rot_0 = w
            dv.setFloat32(o, 1, true); o += 4;
            dv.setFloat32(o, 0, true); o += 4;
            dv.setFloat32(o, 0, true); o += 4;
            dv.setFloat32(o, 0, true); o += 4;

            // opacity stored as logit (viewers apply sigmoid)
            const alpha = Math.min(Math.max(this.alphas[i], 1e-3), 1 - 1e-3);
            dv.setFloat32(o, Math.log(alpha / (1 - alpha)), true); o += 4;

            // color as SH DC term: f_dc = (color_normalized - 0.5) / C0
            const [r, g, b] = Color.toRgb(this.colors[i]);
            dv.setFloat32(o, (r / 255 - 0.5) / C0, true); o += 4;
            dv.setFloat32(o, (g / 255 - 0.5) / C0, true); o += 4;
            dv.setFloat32(o, (b / 255 - 0.5) / C0, true); o += 4;
        }

        const ply = new Uint8Array(headerBytes.length + body.byteLength);
        ply.set(headerBytes, 0);
        ply.set(new Uint8Array(body), headerBytes.length);
        return { ply };
    }

    async getBlob(ctx: RuntimeContext): Promise<Blob> {
        return new Blob([(await this.getData()).ply], { type: 'application/octet-stream' });
    }

    constructor(boundingBox: Box3D) {
        super();
        const tmpV = Vec3();
        Vec3.add(tmpV, boundingBox.min, boundingBox.max);
        Vec3.scale(tmpV, tmpV, -0.5);
        this.centerTransform = Mat4.fromTranslation(Mat4(), tmpV);
    }
}
