/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadPositions, } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { AttributeSpec, DefineSpec, TextureSpec, UniformSpec, Values, ValueSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { background_frag } from '../../mol-gl/shader/background.frag';
import { background_vert } from '../../mol-gl/shader/background.vert';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createCubeTexture, createNullTexture, createTexture, CubeFaces, ImageTexture, Texture } from '../../mol-gl/webgl/texture';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { ValueCell } from '../../mol-util/value-cell';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { isTimingMode } from '../../mol-util/debug';
import { Camera, ICamera } from '../camera';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { Color } from '../../mol-util/color';

const SkyboxParams = {
    size: PD.Select(512, [[256, '256x256'], [512, '512x512'], [1024, '1024x1024'], [2048, '2048x2048'], [4096, '4096x4096']] as const), // TODO: remove
    faces: PD.MappedStatic('urls', {
        urls: PD.Group({
            nx: PD.Text('', { label: 'Negative X' }),
            ny: PD.Text('', { label: 'Negative Y' }),
            nz: PD.Text('', { label: 'Negative Z' }),
            px: PD.Text('', { label: 'Positive X' }),
            py: PD.Text('', { label: 'Positive Y' }),
            pz: PD.Text('', { label: 'Positive Z' }),
        }, { isExpanded: true, label: 'URLs' }),
        // TODO: files
    })
};
type SkyboxProps = PD.Values<typeof SkyboxParams>

const ImageParams = {
    source: PD.MappedStatic('url', {
        url: PD.Text(''),
        // TODO: file
    })
};
type ImageProps = PD.Values<typeof ImageParams>

export const BackgroundParams = {
    variant: PD.MappedStatic('off', {
        off: PD.EmptyGroup(),
        skybox: PD.Group(SkyboxParams),
        image: PD.Group(ImageParams),
        horizontalGradient: PD.Group({
            topColor: PD.Color(Color(0xDDDDDD)),
            bottomColor: PD.Color(Color(0xEEEEEE)),
            ratio: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }),
        }),
        radialGradient: PD.Group({
            centerColor: PD.Color(Color(0xDDDDDD)),
            edgeColor: PD.Color(Color(0xEEEEEE)),
            ratio: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }),
        }),
    }),
    opacity: PD.Numeric(1, { min: 0.0, max: 1.0, step: 0.01 }, { hideIf: p => p?.variant === 'off' }),
};
export type BackgroundProps = PD.Values<typeof BackgroundParams>

export class BackgroundPass {
    private renderable: BackgroundRenderable;

    private skybox: ImageTexture | undefined;
    private skyboxProps: SkyboxProps | undefined;

    private image: ImageTexture | undefined;
    private imageProps: ImageProps | undefined;

    private readonly camera = new Camera();
    private readonly target = Vec3();
    private readonly position = Vec3();
    private readonly dir = Vec3();

    readonly texture: Texture;

    constructor(private webgl: WebGLContext, width: number, height: number) {
        this.renderable = getBackgroundRenderable(webgl, width, height);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;

        if (width !== w || height !== h) {
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
        }
    }

    private updateSkybox(camera: ICamera, props: SkyboxProps) {
        const tf = this.skyboxProps?.faces;
        const f = props.faces.params;
        if (!f.nx || !f.ny || !f.nz || !f.px || !f.py || !f.pz) {
            this.skybox = undefined;
            this.skyboxProps = undefined;
            return;
        }
        if (!this.skyboxProps || !tf || areSkyboxTexturePropsEqual(this.skyboxProps.faces.params, this.skyboxProps.size, props.faces.params, props.size)) {
            this.skybox = getSkyboxTexture(this.webgl, props.faces.params, props.size);
            ValueCell.update(this.renderable.values.tSkybox, this.skybox);
            this.renderable.update();
            this.skyboxProps = { ...props };
        }
        if (!this.skybox) return;

        let cam = camera;
        if (camera.state.mode === 'orthographic') {
            this.camera.setState({ ...camera.state, mode: 'perspective' });
            this.camera.update();
            cam = this.camera;
        }

        const m = this.renderable.values.uViewDirectionProjectionInverse.ref.value;
        Vec3.sub(this.dir, cam.state.position, cam.state.target);
        Vec3.setMagnitude(this.dir, this.dir, 0.1);
        Vec3.copy(this.position, this.dir);
        Mat4.lookAt(m, this.position, this.target, cam.state.up);
        Mat4.mul(m, cam.projection, m);
        Mat4.invert(m, m);
        ValueCell.update(this.renderable.values.uViewDirectionProjectionInverse, m);

        ValueCell.updateIfChanged(this.renderable.values.dVariant, 'skybox');
        this.renderable.update();
    }

    updateImage(props: ImageProps) {
        if (!props.source.params) {
            this.image = undefined;
            this.imageProps = undefined;
            return;
        }
        if (!this.imageProps || !this.imageProps.source.params || !props.source.params !== !this.imageProps.source.params) {
            this.image = getImageTexture(this.webgl, props.source.params);
            ValueCell.update(this.renderable.values.tImage, this.image);
            this.renderable.update();
            this.imageProps = { ...props };
        }
        if (!this.image) return;

        ValueCell.updateIfChanged(this.renderable.values.dVariant, 'image');
        this.renderable.update();
    }

    updateImageScaling() {
        const v = this.renderable.values;
        const [w, h] = v.uTexSize.ref.value;
        const iw = this.image?.getWidth() || 0;
        const ih = this.image?.getHeight() || 0;
        const r = w / h;
        const ir = iw / ih;
        // responsive scaling with offset
        if (r < ir) {
            ValueCell.update(v.uImageScale, Vec2.set(v.uImageScale.ref.value, iw * h / ih, h));
        } else {
            ValueCell.update(v.uImageScale, Vec2.set(v.uImageScale.ref.value, w, ih * w / iw));
        }
        const [rw, rh] = v.uImageScale.ref.value;
        const sr = rw / rh;
        if (sr > r) {
            ValueCell.update(v.uImageOffset, Vec2.set(v.uImageOffset.ref.value, (1 - r / sr) / 2, 0));
        } else {
            ValueCell.update(v.uImageOffset, Vec2.set(v.uImageOffset.ref.value, 0, (1 - sr / r) / 2));
        }
    }

    updateGradient(colorA: Color, colorB: Color, ratio: number, variant: 'horizontalGradient' | 'radialGradient') {
        ValueCell.update(this.renderable.values.uGradientColorA, Color.toVec3Normalized(this.renderable.values.uGradientColorA.ref.value, colorA));
        ValueCell.update(this.renderable.values.uGradientColorB, Color.toVec3Normalized(this.renderable.values.uGradientColorB.ref.value, colorB));
        ValueCell.updateIfChanged(this.renderable.values.uGradientRatio, ratio);
        ValueCell.updateIfChanged(this.renderable.values.dVariant, variant);
        this.renderable.update();
    }

    update(camera: ICamera, props: BackgroundProps) {
        if (props.variant.name === 'off') {
            this.skyboxProps = undefined;
            return;
        } else if (props.variant.name === 'skybox') {
            this.imageProps = undefined;
            this.updateSkybox(camera, props.variant.params);
        } else if (props.variant.name === 'image') {
            this.skyboxProps = undefined;
            this.updateImage(props.variant.params);
        } else if (props.variant.name === 'horizontalGradient') {
            this.imageProps = undefined;
            this.skyboxProps = undefined;
            this.updateGradient(props.variant.params.topColor, props.variant.params.bottomColor, props.variant.params.ratio, props.variant.name);
        } else if (props.variant.name === 'radialGradient') {
            this.imageProps = undefined;
            this.skyboxProps = undefined;
            this.updateGradient(props.variant.params.centerColor, props.variant.params.edgeColor, props.variant.params.ratio, props.variant.name);
        }
        ValueCell.updateIfChanged(this.renderable.values.uOpacity, props.opacity);
    }

    isEnabled(props: BackgroundProps) {
        return !!(
            (this.skyboxProps && this.skybox?.isLoaded) ||
            (this.imageProps && this.image?.isLoaded) ||
            props.variant.name === 'horizontalGradient' ||
            props.variant.name === 'radialGradient'
        );
    }

    private isReady() {
        return !!(
            (this.skyboxProps && this.skybox?.isLoaded) ||
            (this.imageProps && this.image?.isLoaded) ||
            this.renderable.values.dVariant.ref.value === 'horizontalGradient' ||
            this.renderable.values.dVariant.ref.value === 'radialGradient'
        );
    }

    render() {
        if (!this.isReady()) return;

        if (this.renderable.values.dVariant.ref.value === 'image') {
            this.updateImageScaling();
        }

        if (isTimingMode) this.webgl.timer.mark('BackgroundPass.render');
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('BackgroundPass.render');
    }

    //

    static areTexturePropsEqual(propsNew: BackgroundProps, propsOld: BackgroundProps) {
        if (propsNew.variant.name === 'skybox') {
            if (propsOld.variant.name !== 'skybox') return false;
            return areSkyboxTexturePropsEqual(propsNew.variant.params.faces.params, propsNew.variant.params.size, propsOld.variant.params.faces.params, propsOld.variant.params.size);
        } else if (propsNew.variant.name === 'image') {
            if (propsOld.variant.name !== 'image') return false;
            return areImageTexturePropsEqual(propsNew.variant.params.source.params, propsOld.variant.params.source.params);
        } else {
            return true;
        }
    }

    static loadTexture(ctx: WebGLContext, props: BackgroundProps, onload?: () => void) {
        if (props.variant.name === 'skybox') {
            getSkyboxTexture(ctx, props.variant.params.faces.params, props.variant.params.size, onload);
        } else if (props.variant.name === 'image') {
            getImageTexture(ctx, props.variant.params.source.params, onload);
        }
    }
}

//

const SkyboxName = 'background-skybox';

function getSkyboxHash(faces: CubeFaces, size: number) {
    return `${SkyboxName}_${faces.nx}|${faces.ny}|${faces.nz}|${faces.px}|${faces.py}|${faces.pz}|${size}`;
}

function areSkyboxTexturePropsEqual(facesA: CubeFaces, sizeA: number, facesB: CubeFaces, sizeB: number) {
    return sizeA === sizeB && facesA.nx === facesB.nx && facesA.ny === facesB.ny && facesA.nz === facesB.nz && facesA.px === facesB.px && facesA.py === facesB.py && facesA.pz === facesB.pz;
}

function getSkyboxTexture(ctx: WebGLContext, faces: CubeFaces, size: number, onload?: () => void): ImageTexture {
    const hash = getSkyboxHash(faces, size);
    if (!ctx.namedTextures[hash]) {
        ctx.namedTextures[hash] = createCubeTexture(ctx.gl, faces, size, onload);
    } else if (onload) {
        onload();
    }
    return ctx.namedTextures[hash] as ImageTexture;
}

//

const ImageName = 'background-image';

function getImageHash(source: string) {
    return `${ImageName}_${source}`;
}

function areImageTexturePropsEqual(sourceA: string, sourceB: string) {
    return sourceA === sourceB;
}

function getImageTexture(ctx: WebGLContext, source: string, onload?: () => void): ImageTexture {
    const hash = getImageHash(source);
    if (!ctx.namedTextures[hash]) {
        const texture = {
            ...createTexture(ctx.gl, ctx.extensions, 'image-uint8', 'rgba', 'ubyte', 'linear'),
            isLoaded: false,
        };
        const img = new Image();
        img.onload = () => {
            texture.load(img);
            texture.isLoaded = true;
            onload?.();
        };
        img.src = source;
        ctx.namedTextures[hash] = texture;
    } else if (onload) {
        onload();
    }
    return ctx.namedTextures[hash] as ImageTexture;
}

//

const BackgroundSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),
    aPosition: AttributeSpec('float32', 2, 0),
    tSkybox: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tImage: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uImageScale: UniformSpec('v2'),
    uImageOffset: UniformSpec('v2'),
    uTexSize: UniformSpec('v2'),
    uViewDirectionProjectionInverse: UniformSpec('m4'),
    uGradientColorA: UniformSpec('v3'),
    uGradientColorB: UniformSpec('v3'),
    uGradientRatio: UniformSpec('f'),
    uOpacity: UniformSpec('f'),
    dVariant: DefineSpec('string', ['skybox', 'image', 'verticalGradient', 'horizontalGradient', 'radialGradient']),
};
const SkyboxShaderCode = ShaderCode('background', background_vert, background_frag);
type BackgroundRenderable = ComputeRenderable<Values<typeof BackgroundSchema>>

function getBackgroundRenderable(ctx: WebGLContext, width: number, height: number): BackgroundRenderable {
    const values: Values<typeof BackgroundSchema> = {
        drawCount: ValueCell.create(6),
        instanceCount: ValueCell.create(1),
        aPosition: ValueCell.create(QuadPositions),
        tSkybox: ValueCell.create(createNullTexture()),
        tImage: ValueCell.create(createNullTexture()),
        uImageScale: ValueCell.create(Vec2()),
        uImageOffset: ValueCell.create(Vec2()),
        uTexSize: ValueCell.create(Vec2.create(width, height)),
        uViewDirectionProjectionInverse: ValueCell.create(Mat4()),
        uGradientColorA: ValueCell.create(Vec3()),
        uGradientColorB: ValueCell.create(Vec3()),
        uGradientRatio: ValueCell.create(0.5),
        uOpacity: ValueCell.create(1),
        dVariant: ValueCell.create('skybox'),
    };

    const schema = { ...BackgroundSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', SkyboxShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}
