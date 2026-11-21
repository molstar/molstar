/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome drawing (_spec/README.md "Render", Blocks and Rings layout): runs of small blocks (their extent, band and
// count), the groups of a section of several assemblies (a bracket under a row, an arc and a name round a ring), the
// closing run past `maxAssemblyBlocks`, the ring's title over assemblies. A chunk of its own, landed with the genome
// chunk (`withGenome`, registered by `layout/blocks.ts`): a reader that never opens a genome pays nothing for it.

import { piecesTotal, type ResolvedSection } from '@molstar/seqx-core/view';
import type { RenderConfig } from '../config.js';
import { fontSize } from '../config.js';
import { Anchor, grown, NoColor, type PassWriter } from '../draw-lists.js';
import { elide } from '../text.js';
import type { SeqTokens } from '../tokens.js';
import type { SectionLayout } from '../tracks.js';
import { type BlockPlacement, blockLength, isChip, type ResolvedBlock } from './blocks.js';
import { angleOf, type RingPlacement } from './rings.js';
import type { PlaceContext } from './types.js';

// the anchors as numbers read once: a draw reads no property of an imported object
const AnchorStart: number = Anchor.start;
const AnchorMiddle: number = Anchor.middle;
const AnchorEnd: number = Anchor.end;

/** A count as a run's label says it ("2,841 contigs"). */
export function runLabel(count: number): string {
  return `${count.toLocaleString('en-US')} ${count === 1 ? 'contig' : 'contigs'}`;
}

/** The closing run's label past `maxAssemblyBlocks` ("+96,312 sequences not shown"). */
export function moreLabel(count: number): string {
  return `+${count.toLocaleString('en-US')} ${count === 1 ? 'sequence' : 'sequences'} not shown`;
}

/** The runs in view as laid (a row's): each its blocks, extent and label ('' when it has no room). */
export function layRuns(p: BlockPlacement, ctx: PlaceContext, blocks: readonly ResolvedBlock[]): void {
  const row = p.row;
  p.runN = 0;
  if (p.stacked) return;
  for (let b = 0; b < blocks.length; b++) {
    if (row.run[b] !== b) continue;
    let e = b;
    while (e + 1 < blocks.length && row.run[e + 1] === b) e++;
    // its extent in view: its first block in view to its last
    let x0 = Infinity;
    let x1 = -Infinity;
    for (let k = b; k <= e; k++) {
      if (row.ve[k] <= row.vs[k]) continue;
      x0 = Math.min(x0, row.x[k]);
      x1 = Math.max(x1, row.x[k] + row.w[k]);
    }
    if (!(x1 > x0)) continue;
    const r = p.runN++;
    if (r >= p.runFirst.length) {
      p.runFirst = grown(p.runFirst, r + 1);
      p.runLast = grown(p.runLast, r + 1);
      p.runX = grown(p.runX, r + 1);
      p.runW = grown(p.runW, r + 1);
    }
    p.runFirst[r] = b;
    p.runLast[r] = e;
    p.runX[r] = x0;
    p.runW[r] = x1 - x0;
    const text = runLabel(e - b + 1);
    p.runText[r] = ctx.text.measure(ctx.config.labelFont, text) + 8 <= x1 - x0 ? text : '';
  }
}

/**
 * The row under a row of blocks (_spec/README.md "Render", Blocks): a bracket a group in view with its assembly's name
 * (a section of several assemblies, EC-210), and the closing run's label past `maxAssemblyBlocks` at its end. Its
 * height (0 when it has neither).
 */
export function layFoot(
  p: BlockPlacement,
  section: ResolvedSection,
  ctx: PlaceContext,
  blocks: readonly ResolvedBlock[],
): number {
  const row = p.row;
  const groups = section.groups;
  p.groupN = 0;
  p.groupText.length = 0;
  p.moreText = section.more ? moreLabel(section.more) : '';
  if (groups && groups.length > 1)
    for (const g of groups) {
      let x0 = Infinity;
      let x1 = -Infinity;
      for (let k = g.first; k < g.first + g.count && k < blocks.length; k++) {
        if (row.ve[k] <= row.vs[k]) continue;
        x0 = Math.min(x0, row.x[k]);
        x1 = Math.max(x1, row.x[k] + row.w[k]);
      }
      if (!(x1 > x0)) continue;
      const i = p.groupN++;
      if (i >= p.groupX.length) {
        p.groupX = grown(p.groupX, i + 1);
        p.groupW = grown(p.groupW, i + 1);
      }
      p.groupX[i] = x0;
      p.groupW[i] = x1 - x0;
      p.groupText[i] = elide(ctx.text, ctx.config.labelFont, g.label, Math.max(0, x1 - x0 - 8));
    }
  return p.groupN || p.moreText ? fontSize(ctx.config.labelFont) + 10 : 0;
}

/**
 * The runs' bands under the tracks' data (the data pass, once a section): each a band in `seq-wash` over its extent with a
 * hairline at each end (_spec/README.md "Render", Blocks); through `args`, nothing boxed.
 */
export function drawRuns(s: SectionLayout, p: BlockPlacement, out: PassWriter): void {
  const g = out.args;
  const h = p.tracksBottom - p.tracksTop;
  if (!(h > 0)) return;
  for (let r = 0; r < p.runN; r++) {
    g[0] = p.runX[r];
    g[1] = p.tracksTop;
    g[2] = p.runW[r];
    g[3] = h;
    g[4] = s.tokens.colors['seq-wash'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    out.rectArgs();
    for (let e = 0; e < 2; e++) {
      g[0] = e ? p.runX[r] + p.runW[r] - 1 : p.runX[r];
      g[1] = p.tracksTop;
      g[2] = 1;
      g[3] = h;
      g[4] = s.tokens.colors['msui-text-muted'];
      g[7] = 0;
      out.rectArgs();
    }
  }
}

/**
 * The runs' labels in the label row, the groups' brackets with their assemblies' names and the closing run's label in
 * the row under the tracks (the labels pass, once a section).
 */
export function drawRunLabels(s: SectionLayout, p: BlockPlacement, out: PassWriter, size: number): void {
  const g = out.args;
  const muted = s.tokens.colors['msui-text-muted'];
  const font = s.config.labelFont;
  for (let r = 0; r < p.runN; r++) {
    const text = p.runText[r];
    if (!text) continue;
    g[0] = p.runX[r] + 4;
    g[1] = p.tracksTop - p.labelRow + size;
    g[2] = 0;
    g[3] = Math.max(0, p.runW[r] - 8);
    g[4] = muted;
    out.textArgs(text, AnchorStart, font);
  }
  if (!p.footRow) return;
  const y = p.tracksBottom + 4;
  for (let i = 0; i < p.groupN; i++) {
    // a bracket: a hairline across the group, a tick up at each end, its name in the middle below it
    g[0] = p.groupX[i];
    g[1] = y;
    g[2] = p.groupW[i];
    g[3] = 1;
    g[4] = muted;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    out.rectArgs();
    for (let e = 0; e < 2; e++) {
      g[0] = e ? p.groupX[i] + p.groupW[i] - 1 : p.groupX[i];
      g[1] = y - 4;
      g[2] = 1;
      g[3] = 5;
      g[4] = muted;
      g[7] = 0;
      out.rectArgs();
    }
    if (!p.groupText[i]) continue;
    g[0] = p.groupX[i] + p.groupW[i] / 2;
    g[1] = y + size + 3;
    g[2] = 0;
    g[3] = p.groupW[i];
    g[4] = s.tokens.colors['msui-text'];
    out.textArgs(p.groupText[i], AnchorMiddle, font);
  }
  if (p.moreText) {
    const a = p.row;
    g[0] = a.x0 + a.width;
    g[1] = y + size + 3;
    g[2] = 0;
    g[3] = a.width;
    g[4] = muted;
    out.textArgs(p.moreText, AnchorEnd, font);
  }
}

/** A length in bases as a genome's title says it: `3.1 Gb`, `48 Mb`, `16.6 kb`, `950 bp`. */
export function basesText(n: number): string {
  const units: [number, string][] = [
    [1e9, 'Gb'],
    [1e6, 'Mb'],
    [1e3, 'kb'],
  ];
  for (const [k, u] of units) if (n >= k) return `${+(n / k).toPrecision(n / k >= 100 ? 3 : 2)} ${u}`;
  return `${Math.round(n)} bp`;
}

/**
 * A ring's runs at the outer radius's zoom (px a unit): two or more consecutive blocks of one group, no chips, each
 * under `blocks.runPx` (_spec/README.md "Render", Rings layout); each block's run start or −1, undefined with none.
 */
export function ringRuns(
  p: RingPlacement,
  blocks: readonly ResolvedBlock[] | undefined,
  config: RenderConfig,
  zoom: number,
): Int32Array | undefined {
  const small = config.blocks.runPx;
  if (!blocks || !(small > 0)) return undefined;
  const n = blocks.length;
  p.run = grown(p.run, n);
  const tiny = (b: number) => !isChip(blocks[b]) && blockLength(blocks[b]) * zoom < small;
  let any = false;
  let b = 0;
  while (b < n) {
    if (!tiny(b)) {
      p.run[b++] = -1;
      continue;
    }
    let e = b + 1;
    while (e < n && tiny(e) && (blocks[e].group ?? -1) === (blocks[b].group ?? -1)) e++;
    const v = e - b > 1 ? b : -1;
    if (v >= 0) any = true;
    for (let k = b; k < e; k++) p.run[k] = v;
    b = e;
  }
  return any ? p.run : undefined;
}

/** A block of a run laid at frame x: its run's extent grown (a run's first block opens it). */
export function addRun(p: RingPlacement, start: number, b: number, x: number, w: number): void {
  if (start === b) {
    const r = p.runN++;
    if (r >= p.runFirst.length) {
      p.runFirst = grown(p.runFirst, r + 1);
      p.runLast = grown(p.runLast, r + 1);
      p.runX = grown(p.runX, r + 1);
      p.runW = grown(p.runW, r + 1);
    }
    p.runFirst[r] = b;
    p.runX[r] = x;
    p.runW[r] = 0;
  }
  const r = p.runN - 1;
  p.runLast[r] = b;
  p.runW[r] = x + w - p.runX[r];
}

/**
 * A ring's runs' labels where they have room (with `outside`) and its groups' arcs (a section of several assemblies,
 * EC-210): each group's blocks' frame extent, its name among the outside labels.
 */
export function ringRunsAndGroups(
  p: RingPlacement,
  section: ResolvedSection,
  blocks: readonly ResolvedBlock[] | undefined,
  outside: boolean,
): void {
  if (outside)
    for (let k = 0; k < p.runN; k++)
      p.labels.add(
        runLabel(p.runLast[k] - p.runFirst[k] + 1),
        angleOf(p, p.runX[k] + p.runW[k] / 2),
        Infinity,
        -1,
        undefined,
      );
  const groups = section.groups;
  if (!groups || groups.length < 2 || !blocks) return;
  const s = p.strips;
  for (const g of groups) {
    let x0 = Infinity;
    let x1 = -Infinity;
    for (let i = 0; i < s.n; i++) {
      const b = s.block[i];
      if (b < g.first || b >= g.first + g.count) continue;
      x0 = Math.min(x0, s.axes[i].x0);
      x1 = Math.max(x1, s.axes[i].x0 + s.axes[i].width);
    }
    for (let r = 0; r < p.runN; r++)
      if (p.runFirst[r] >= g.first && p.runFirst[r] < g.first + g.count) {
        x0 = Math.min(x0, p.runX[r]);
        x1 = Math.max(x1, p.runX[r] + p.runW[r]);
      }
    if (!(x1 > x0)) continue;
    const k = p.groupN++;
    if (k >= p.groupX.length) {
      p.groupX = grown(p.groupX, k + 1);
      p.groupW = grown(p.groupW, k + 1);
    }
    p.groupX[k] = x0;
    p.groupW[k] = x1 - x0;
    if (outside) p.labels.add(g.label, angleOf(p, (x0 + x1) / 2), Infinity, -2, undefined);
  }
}

/** A ring's title over assemblies: their names, and the blocks' length in bases; the closing run's count after. */
export function groupTitle(section: ResolvedSection): [string, string] {
  let bases = 0;
  for (const b of section.blocks ?? []) bases += piecesTotal(b.axis.ranges ?? []) * b.members.length;
  const more = section.more ? ` (+${section.more.toLocaleString('en-US')} not shown)` : '';
  return [(section.groups ?? []).map((g) => g.label).join(' · '), `${basesText(bases)}${more}`];
}

/**
 * A ring's runs (a band in `seq-wash` over each run's arc through the rings, a hairline at each end) and its groups'
 * arcs (a hairline just outside the outer ring), through the warp in the data pass (_spec/README.md "Render", Rings
 * layout); through `args`, nothing boxed.
 */
export function drawRingRuns(p: RingPlacement, page: PassWriter, tokens: SeqTokens): void {
  const out = p.warp.writer(page);
  const g = out.args;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  for (let r = 0; r < p.runN; r++) {
    g[0] = p.runX[r];
    g[1] = 0;
    g[2] = p.runW[r];
    g[3] = p.depth;
    g[4] = tokens.colors['seq-wash'];
    out.rectArgs();
    for (let e = 0; e < 2; e++) {
      g[0] = e ? p.runX[r] + p.runW[r] - 0.75 : p.runX[r];
      g[1] = 0;
      g[2] = 0.75;
      g[3] = p.depth;
      g[4] = tokens.colors['msui-text-muted'];
      out.rectArgs();
    }
  }
  for (let k = 0; k < p.groupN; k++) {
    g[0] = p.groupX[k];
    g[1] = -4;
    g[2] = p.groupW[k];
    g[3] = 1.5;
    g[4] = tokens.colors['msui-text-muted'];
    out.rectArgs();
  }
  p.warp.close(page);
}
