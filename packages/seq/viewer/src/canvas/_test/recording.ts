/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A fake 2D context that logs what it is asked to do, and a canvas that holds one: the Canvas2D backend is tested
// against these, with no DOM.

import type { CanvasLike, Context2D } from '../backend.js';

const Recorded = [
  'fillStyle',
  'strokeStyle',
  'lineWidth',
  'font',
  'textAlign',
  'textBaseline',
  'imageSmoothingEnabled',
];

const show = (a: unknown) => (typeof a === 'number' ? String(Math.round(a * 100) / 100) : String(a));

export class RecordingContext implements Context2D {
  readonly log: string[] = [];
  private readonly _props: Record<string, unknown> = {};

  constructor() {
    for (const p of Recorded)
      Object.defineProperty(this, p, {
        get: () => this._props[p],
        set: (v) => {
          this._props[p] = v;
          this.log.push(`${p}=${typeof v === 'string' ? v : typeof v === 'number' ? v : '[object]'}`);
        },
      });
  }

  declare fillStyle: unknown;
  declare strokeStyle: unknown;
  declare lineWidth: number;
  declare font: string;
  declare textAlign: string;
  declare textBaseline: string;
  declare imageSmoothingEnabled: boolean;

  private call(name: string, ...args: unknown[]) {
    this.log.push(`${name}(${args.map(show).join(', ')})`);
  }

  save() {
    this.call('save');
  }
  restore() {
    this.call('restore');
  }
  setTransform(a: number, b: number, c: number, d: number, e: number, f: number) {
    this.call('setTransform', a, b, c, d, e, f);
  }
  clearRect(x: number, y: number, w: number, h: number) {
    this.call('clearRect', x, y, w, h);
  }
  fillRect(x: number, y: number, w: number, h: number) {
    this.call('fillRect', x, y, w, h);
  }
  strokeRect(x: number, y: number, w: number, h: number) {
    this.call('strokeRect', x, y, w, h);
  }
  beginPath() {
    this.call('beginPath');
  }
  closePath() {
    this.call('closePath');
  }
  moveTo(x: number, y: number) {
    this.call('moveTo', x, y);
  }
  lineTo(x: number, y: number) {
    this.call('lineTo', x, y);
  }
  bezierCurveTo(a: number, b: number, c: number, d: number, e: number, f: number) {
    this.call('bezierCurveTo', a, b, c, d, e, f);
  }
  rect(x: number, y: number, w: number, h: number) {
    this.call('rect', x, y, w, h);
  }
  fill() {
    this.call('fill');
  }
  stroke() {
    this.call('stroke');
  }
  clip() {
    this.call('clip');
  }
  setLineDash(s: number[]) {
    this.call('setLineDash', s.join('/'));
  }
  fillText(t: string, x: number, y: number) {
    this.call('fillText', JSON.stringify(t), x, y);
  }
  createLinearGradient(x0: number, y0: number, x1: number, y1: number) {
    this.call('createLinearGradient', x0, y0, x1, y1);
    return { addColorStop: (o: number, c: string) => this.call('addColorStop', o, c) };
  }
  createImageData(w: number, h: number) {
    this.call('createImageData', w, h);
    return { data: new Uint8ClampedArray(w * h * 4) };
  }
  putImageData(_image: unknown, x: number, y: number) {
    this.call('putImageData', x, y);
  }
  drawImage(_source: unknown, x: number, y: number, w: number, h: number) {
    this.call('drawImage', x, y, w, h);
  }

  /** The calls made, with the log cleared. */
  take(): string[] {
    return this.log.splice(0);
  }

  /** The calls of one name. */
  calls(name: string): string[] {
    return this.log.filter((l) => l.startsWith(`${name}(`));
  }
}

export class FakeCanvas implements CanvasLike {
  width = 300;
  height = 150;
  style = { width: '', height: '' };
  readonly context = new RecordingContext();
  getContext() {
    return this.context;
  }
}
