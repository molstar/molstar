/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { inflateSync } from 'node:zlib';
import { printJson } from '@molstar/core/repr/codecs/json';
import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// The render CLI's test inputs: p53 as tracks (its letters, domains and conservation, all inline), a spec of its own
// folder with a file source, and a PNG's pixels decoded without a dependency.

/** Human p53 (UniProt P04637), the first 120 residues. */
export const P53 =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMF';

/** p53's tracks: letters, domains and a conservation value a residue. */
export function p53Tracks(): Node {
  const b = createSeqSpec({ title: 'p53 tracks' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: {
      start: new Int32Array([0, 63, 93]),
      end: new Int32Array([61, 94, 120]),
      id: ['TAD', 'PRR', 'DBD'],
      label: ['Transactivation', 'Proline-rich', 'DNA-binding'],
    },
  });
  const values = new Float32Array(P53.length);
  for (let i = 0; i < values.length; i++) values[i] = (i % 10) / 10;
  const conservation = b.data.values('conservation', { target: p53, values });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    s.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: conservation });
  });
  return b.build();
}

/** p53's letters alone: no columns, so its S-expression is whole (a printed column is elided). */
export function p53Letters(): Node {
  const b = createSeqSpec({ title: 'p53 letters' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  b.views.view('main', (v) => {
    v.section('protein', { of: p53 }).track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** p53 whose letters are read from a FASTA file at `uri` (relative to the spec's folder). */
export function p53FromFile(uri: string): Node {
  const b = createSeqSpec({ title: 'p53 from a file' });
  const src = b.data.source('fa', { format: 'fasta', uri });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', source: src });
  b.views.view('main', (v) => {
    v.section('protein', { of: p53 }).track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

export const JsonText = (n: Node) => printJson(n);
export const SexprText = (n: Node) => printSexpr(n);

/** A scratch folder, removed by `dispose`. */
export function scratch() {
  const dir = mkdtempSync(join(tmpdir(), 'seqx-render-'));
  return {
    dir,
    file: (name: string, content: string | Uint8Array) => {
      const p = join(dir, name);
      writeFileSync(p, content);
      return p;
    },
    dispose: () => rmSync(dir, { recursive: true, force: true }),
  };
}

export interface Decoded {
  readonly width: number;
  readonly height: number;
  /** RGBA, row-major. */
  readonly rgba: Uint8Array;
  /** The pixel at (x, y) as [r, g, b, a]. */
  at(x: number, y: number): [number, number, number, number];
}

/** An 8-bit RGBA PNG's pixels (any filter, one or more IDAT). */
export function decodePng(png: Uint8Array): Decoded {
  const v = new DataView(png.buffer, png.byteOffset, png.byteLength);
  if (v.getUint32(0) !== 0x89504e47) throw new Error('not a PNG');
  let at = 8;
  let width = 0;
  let height = 0;
  const data: Uint8Array[] = [];
  while (at < png.length) {
    const n = v.getUint32(at);
    const type = String.fromCharCode(...png.subarray(at + 4, at + 8));
    const body = png.subarray(at + 8, at + 8 + n);
    if (type === 'IHDR') {
      width = v.getUint32(at + 8);
      height = v.getUint32(at + 12);
      if (png[at + 16] !== 8 || png[at + 17] !== 6) throw new Error('only 8-bit RGBA');
    } else if (type === 'IDAT') data.push(body);
    at += 12 + n;
  }
  const raw = inflateSync(Buffer.concat(data));
  const row = width * 4;
  const out = new Uint8Array(row * height);
  for (let y = 0; y < height; y++) {
    const f = raw[y * (row + 1)];
    for (let x = 0; x < row; x++) {
      const a = x >= 4 ? out[y * row + x - 4] : 0;
      const b = y > 0 ? out[(y - 1) * row + x] : 0;
      const c = x >= 4 && y > 0 ? out[(y - 1) * row + x - 4] : 0;
      const cur = raw[y * (row + 1) + 1 + x];
      let p = 0;
      if (f === 1) p = a;
      else if (f === 2) p = b;
      else if (f === 3) p = (a + b) >> 1;
      else if (f === 4) {
        const pa = Math.abs(b - c);
        const pb = Math.abs(a - c);
        const pc = Math.abs(a + b - 2 * c);
        p = pa <= pb && pa <= pc ? a : pb <= pc ? b : c;
      }
      out[y * row + x] = (cur + p) & 0xff;
    }
  }
  return {
    width,
    height,
    rgba: out,
    at: (x, y) => {
      const i = (y * width + x) * 4;
      return [out[i], out[i + 1], out[i + 2], out[i + 3]];
    },
  };
}
