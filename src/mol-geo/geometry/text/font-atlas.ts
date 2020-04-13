/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { edt } from '../../../mol-math/geometry/distance-transform';
import { createTextureImage, TextureImage } from '../../../mol-gl/renderable/util';

const TextAtlasCache: { [k: string]: FontAtlas } = {};

export function getFontAtlas (props: Partial<FontAtlasProps>) {
    const hash = JSON.stringify(props);
    if (TextAtlasCache[hash] === undefined) {
        TextAtlasCache[hash] = new FontAtlas(props);
    }
    return TextAtlasCache[hash];
}

export type FontFamily = 'sans-serif' | 'monospace' | 'serif' | 'cursive'
export type FontStyle = 'normal' | 'italic' | 'oblique'
export type FontVariant = 'normal' | 'small-caps'
export type FontWeight = 'normal' | 'bold'

export const FontAtlasParams = {
    fontFamily: PD.Select('sans-serif', [['sans-serif', 'Sans Serif'], ['monospace', 'Monospace'], ['serif', 'Serif'], ['cursive', 'Cursive']] as [FontFamily, string][]),
    fontQuality: PD.Select(3, [[0, 'lower'], [1, 'low'], [2, 'medium'], [3, 'high'], [4, 'higher']]),
    fontStyle: PD.Select('normal', [['normal', 'Normal'], ['italic', 'Italic'], ['oblique', 'Oblique']] as [FontStyle, string][]),
    fontVariant: PD.Select('normal', [['normal', 'Normal'], ['small-caps', 'Small Caps']] as [FontVariant, string][]),
    fontWeight: PD.Select('normal', [['normal', 'Normal'], ['bold', 'Bold']] as [FontWeight, string][]),
};
export type FontAtlasParams = typeof FontAtlasParams
export type FontAtlasProps = PD.Values<FontAtlasParams>

export type FontAtlasMap = {
    x: number, y: number, w: number, h: number,
    nw: number, nh: number // normalized to lineheight
}

export class FontAtlas {
    readonly props: Readonly<FontAtlasProps>
    readonly mapped: { [k: string]: FontAtlasMap } = {}
    readonly placeholder: FontAtlasMap
    readonly texture: TextureImage<Uint8Array>

    private scratchW = 0
    private scratchH = 0
    private currentX = 0
    private currentY = 0
    private readonly scratchData: Uint8Array

    private readonly cutoff = 0.5
    readonly buffer: number
    private readonly radius: number

    private gridOuter: Float64Array
    private gridInner: Float64Array
    private f: Float64Array
    private d: Float64Array
    private z: Float64Array
    private v: Int16Array

    private scratchCanvas: HTMLCanvasElement
    private scratchContext: CanvasRenderingContext2D

    readonly lineHeight: number

    private readonly maxWidth: number
    private readonly middle: number

    constructor (props: Partial<FontAtlasProps> = {}) {
        const p = { ...PD.getDefaultValues(FontAtlasParams), ...props };
        this.props = p;

        // create measurements
        const fontSize = 32 * (p.fontQuality + 1);
        this.buffer = fontSize / 8;
        this.radius = fontSize / 3;
        this.lineHeight = Math.round(fontSize + 2 * this.buffer + this.radius);
        this.maxWidth = Math.round(this.lineHeight * 0.75);

        // create texture (for ~350 characters)
        this.texture = createTextureImage(350 * this.lineHeight * this.maxWidth, 1, Uint8Array);

        // prepare scratch canvas
        this.scratchCanvas = document.createElement('canvas');
        this.scratchCanvas.width = this.maxWidth;
        this.scratchCanvas.height = this.lineHeight;

        this.scratchContext = this.scratchCanvas.getContext('2d')!;
        this.scratchContext.font = `${p.fontStyle} ${p.fontVariant} ${p.fontWeight} ${fontSize}px ${p.fontFamily}`;
        this.scratchContext.fillStyle = 'black';
        this.scratchContext.textBaseline = 'middle';

        // SDF scratch values
        this.scratchData = new Uint8Array(this.lineHeight * this.maxWidth);

        // temporary arrays for the distance transform
        this.gridOuter = new Float64Array(this.lineHeight * this.maxWidth);
        this.gridInner = new Float64Array(this.lineHeight * this.maxWidth);
        this.f = new Float64Array(Math.max(this.lineHeight, this.maxWidth));
        this.d = new Float64Array(Math.max(this.lineHeight, this.maxWidth));
        this.z = new Float64Array(Math.max(this.lineHeight, this.maxWidth) + 1);
        this.v = new Int16Array(Math.max(this.lineHeight, this.maxWidth));

        this.middle = Math.ceil(this.lineHeight / 2);

        // replacement Character
        this.placeholder = this.get(String.fromCharCode(0xFFFD));
    }

    get (char: string) {
        if (this.mapped[char] === undefined) {
            this.draw(char);

            const { array, width, height } = this.texture;
            const data = this.scratchData;

            if (this.currentX + this.scratchW > width) {
                this.currentX = 0;
                this.currentY += this.scratchH;
            }
            if (this.currentY + this.scratchH > height) {
                console.warn('canvas to small');
                return this.placeholder;
            }

            this.mapped[char] = {
                x: this.currentX, y: this.currentY,
                w: this.scratchW, h: this.scratchH,
                nw: this.scratchW / this.lineHeight, nh: this.scratchH / this.lineHeight
            };

            for (let y = 0; y < this.scratchH; ++y) {
                for (let x = 0; x < this.scratchW; ++x) {
                    array[width * (this.currentY + y) + this.currentX + x] = data[y * this.scratchW + x];
                }
            }

            this.currentX += this.scratchW;
        }

        return this.mapped[char];
    }

    draw (char: string) {
        const h = this.lineHeight;
        const ctx = this.scratchContext;
        const data = this.scratchData;

        // measure text
        const m = ctx.measureText(char);
        const w = Math.min(this.maxWidth, Math.ceil(m.width + 2 * this.buffer));
        const n = w * h;

        ctx.clearRect(0, 0, w, h); // clear scratch area
        ctx.fillText(char, this.buffer, this.middle); // draw text
        const imageData = ctx.getImageData(0, 0, w, h);

        for (let i = 0; i < n; i++) {
            const a = imageData.data[i * 4 + 3] / 255; // alpha value
            this.gridOuter[i] = a === 1 ? 0 : a === 0 ? Number.MAX_SAFE_INTEGER : Math.pow(Math.max(0, 0.5 - a), 2);
            this.gridInner[i] = a === 1 ? Number.MAX_SAFE_INTEGER : a === 0 ? 0 : Math.pow(Math.max(0, a - 0.5), 2);
        }

        edt(this.gridOuter, w, h, this.f, this.d, this.v, this.z);
        edt(this.gridInner, w, h, this.f, this.d, this.v, this.z);

        for (let i = 0; i < n; i++) {
            const d = this.gridOuter[i] - this.gridInner[i];
            data[i] = Math.max(0, Math.min(255, Math.round(255 - 255 * (d / this.radius + this.cutoff))));
        }

        this.scratchW = w;
        this.scratchH = h;
    }
}