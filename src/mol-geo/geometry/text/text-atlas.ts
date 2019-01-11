/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { edt } from 'mol-math/geometry/distance-transform';

const TextAtlasCache: { [k: string]: TextAtlas } = {}

export function getTextAtlas (props: Partial<TextAtlasProps>) {
  const hash = JSON.stringify(props)
  if (TextAtlasCache[ hash ] === undefined) {
    TextAtlasCache[ hash ] = new TextAtlas(props)
  }
  return TextAtlasCache[ hash ]
}

export type TextFonts = 'sans-serif' | 'monospace' | 'serif' | 'cursive'
export type TextStyles = 'normal' | 'italic' | 'oblique'
export type TextVariants = 'normal' | 'small-caps'
export type TextWeights = 'normal' | 'bold'

export const TextAtlasParams = {
  fontFamily: PD.Select('sans-serif', [['sans-serif', 'Sans Serif'], ['monospace', 'Monospace'], ['serif', 'Serif'], ['cursive', 'Cursive']] as [TextFonts, string][]),
  fontSize: PD.Numeric(36, { min: 4, max: 96, step: 1 }),
  fontStyle: PD.Select('normal', [['normal', 'Normal'], ['italic', 'Italic'], ['oblique', 'Oblique']] as [TextStyles, string][]),
  fontVariant: PD.Select('normal', [['normal', 'Normal'], ['small-caps', 'Small Caps']] as [TextVariants, string][]),
  fontWeight: PD.Select('normal', [['normal', 'Normal'], ['bold', 'Bold']] as [TextWeights, string][]),

  width: PD.Numeric(1024),
  height: PD.Numeric(1024)
}
export type TextAtlasParams = typeof TextAtlasParams
export type TextAtlasProps = PD.Values<TextAtlasParams>

export type TextAtlasMap = { x: number, y: number, w: number, h: number }

export class TextAtlas {
    readonly props: Readonly<TextAtlasProps>
    readonly mapped: { [k: string]: TextAtlasMap } = {}
    readonly placeholder: TextAtlasMap
    readonly context: CanvasRenderingContext2D

    private canvas: HTMLCanvasElement

    private scratchW = 0
    private scratchH = 0
    private currentX = 0
    private currentY = 0

    private readonly cutoff = 0.25
    private padding: number
    private radius: number

    private gridOuter: Float64Array
    private gridInner: Float64Array
    private f: Float64Array
    private d: Float64Array
    private z: Float64Array
    private v: Int16Array

    private scratchCanvas: HTMLCanvasElement
    private scratchContext: CanvasRenderingContext2D

    private lineHeight: number
    private maxWidth: number
    private middle: number

    constructor (props: Partial<TextAtlasProps> = {}) {
        const p = { ...PD.getDefaultValues(TextAtlasParams), ...props }
        this.props = p

        this.padding = p.fontSize / 8
        this.radius = p.fontSize / 3
        this.lineHeight = Math.round(p.fontSize + 6 * this.padding)
        this.maxWidth = this.lineHeight * 1.5

        // Prepare scratch canvas
        this.scratchCanvas = document.createElement('canvas')
        this.scratchCanvas.width = this.maxWidth
        this.scratchCanvas.height = this.lineHeight

        this.scratchContext = this.scratchCanvas.getContext('2d')!
        this.scratchContext.font = `${p.fontStyle} ${p.fontVariant} ${p.fontWeight} ${p.fontSize}px ${p.fontFamily}`
        this.scratchContext.fillStyle = 'black'
        this.scratchContext.textBaseline = 'middle'

        // temporary arrays for the distance transform
        this.gridOuter = new Float64Array(this.lineHeight * this.maxWidth)
        this.gridInner = new Float64Array(this.lineHeight * this.maxWidth)
        this.f = new Float64Array(Math.max(this.lineHeight, this.maxWidth))
        this.d = new Float64Array(Math.max(this.lineHeight, this.maxWidth))
        this.z = new Float64Array(Math.max(this.lineHeight, this.maxWidth) + 1)
        this.v = new Int16Array(Math.max(this.lineHeight, this.maxWidth))

        this.middle = Math.ceil(this.lineHeight / 2)

        //

        this.canvas = document.createElement('canvas')
        this.canvas.width = p.width
        this.canvas.height = p.height
        this.context = this.canvas.getContext('2d')!

        // Replacement Character
        this.placeholder = this.map(String.fromCharCode(0xFFFD))

        // Basic Latin (subset)
        for (let i = 0x0020; i <= 0x007E; ++i) this.map(String.fromCharCode(i))

        // TODO: to slow to always prepare them
        // // Latin-1 Supplement (subset)
        // for (let i = 0x00A1; i <= 0x00FF; ++i) this.map(String.fromCharCode(i))

        // Degree sign
        this.map(String.fromCharCode(0x00B0))

        // // Greek and Coptic (subset)
        // for (let i = 0x0391; i <= 0x03C9; ++i) this.map(String.fromCharCode(i))

        // // Cyrillic (subset)
        // for (let i = 0x0400; i <= 0x044F; ++i) this.map(String.fromCharCode(i))

        // Angstrom Sign
        this.map(String.fromCharCode(0x212B))
    }

    map (text: string) {
        if (this.mapped[text] === undefined) {
            this.draw(text)

            if (this.currentX + this.scratchW > this.props.width) {
                this.currentX = 0
                this.currentY += this.scratchH
            }
            if (this.currentY + this.scratchH > this.props.height) {
                console.warn('canvas to small')
            }

            this.mapped[text] = {
                x: this.currentX, y: this.currentY,
                w: this.scratchW, h: this.scratchH
            }

            this.context.drawImage(
                this.scratchCanvas,
                0, 0, this.scratchW, this.scratchH,
                this.currentX, this.currentY, this.scratchW, this.scratchH
            )

            this.currentX += this.scratchW
        }

        return this.mapped[text]
    }

    get (text: string) {
        return this.mapped[text] || this.placeholder
    }

    draw (text: string) {
        const h = this.lineHeight
        const ctx = this.scratchContext

        // Measure text
        const m = ctx.measureText(text)
        const w = Math.min(this.maxWidth, Math.ceil(m.width + 2 * this.padding + this.radius / 2))
        const n = w * h

        ctx.clearRect(0, 0, w, h) // clear scratch area
        ctx.fillText(text, this.padding + this.radius / 4, this.middle) // draw text

        const imageData = ctx.getImageData(0, 0, w, h)
        const data = imageData.data

        for (let i = 0; i < n; i++) {
            const a = imageData.data[i * 4 + 3] / 255 // alpha value
            this.gridOuter[i] = a === 1 ? 0 : a === 0 ? Number.MAX_SAFE_INTEGER : Math.pow(Math.max(0, 0.5 - a), 2)
            this.gridInner[i] = a === 1 ? Number.MAX_SAFE_INTEGER : a === 0 ? 0 : Math.pow(Math.max(0, a - 0.5), 2)
        }

        edt(this.gridOuter, w, h, this.f, this.d, this.v, this.z)
        edt(this.gridInner, w, h, this.f, this.d, this.v, this.z)

        for (let i = 0; i < n; i++) {
            const d = this.gridOuter[i] - this.gridInner[i];
            data[i * 4 + 3] = Math.max(0, Math.min(255, Math.round(255 - 255 * (d / this.radius + this.cutoff))));
        }

        ctx.putImageData(imageData, 0, 0)
        this.scratchW = w
        this.scratchH = h
    }
}