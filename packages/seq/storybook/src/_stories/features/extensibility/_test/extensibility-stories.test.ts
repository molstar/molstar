/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  DemoExtension,
  DomainsLayer,
  Letters,
  SiteStarts,
  extensibilityComposition,
  extensibilitySpec,
  issueLines,
} from '../compose.js';
import { FocusSites, LollipopsType, WarmScheme } from '../demo-extension.js';
import { defineViewerExtension } from '@molstar/seqx-viewer/extension';

// The Extensibility story, headless on a fake clock, composed as the story composes it (compose.ts): the viewer
// without the example extension draws the sites as their fallback and says so; `extend` draws them as lollipops on the
// canvas and in the SVG export (one draw list), the warm scheme colors the domains, the binding runs the command; the
// undo brings the fallback back. And the extension's params are checked where it is composed.

const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  vi.restoreAllMocks();
  vi.useRealTimers();
});

async function open() {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  vi.spyOn(DomHost, 'requestAnimationFrame').mockImplementation((fn) => {
    setTimeout(() => fn(performance.now()), 16);
    return undefined as never;
  });
  const v = createSeqViewer(extensibilityComposition(false));
  viewers.push(v);
  v.resize(900);
  await v.load(extensibilitySpec());
  await settle(v);
  return v;
}

async function settle(v: SeqViewer) {
  await vi.advanceTimersByTimeAsync(100);
  await v.settle();
  v.draw();
}

const types = (v: SeqViewer) =>
  v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers.map((l) => l.type.name)));
const codes = (v: SeqViewer) => issueLines(v).map((l) => l.split(' ')[1]);

describe('the Extensibility story', () => {
  it('has each site (S6, S15, S33, S46) on a serine of p53, its label its residue number', async () => {
    expect(SiteStarts.map((i) => Letters[i])).toEqual(['S', 'S', 'S', 'S']);
    const v = await open();
    await v.extend(DemoExtension);
    await settle(v);
    const layer = v
      .frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers))
      .find((l) => l.type.name === 'demo:lollipops')!;
    const items = layer.type.items!(layer.lay);
    // each site's label names its residue (1-based) and its extent is that residue
    expect(items.map((i) => [i.label, i.start])).toEqual(SiteStarts.map((x) => [`S${x + 1}`, x]));
    expect(items.map((i) => Letters[i.start])).toEqual(['S', 'S', 'S', 'S']);
  });

  it('draws the fallback without the extension, the lollipops with it, the fallback again after its undo', async () => {
    const v = await open();
    expect(types(v)).toEqual(['letters', 'ranges', 'ranges']);
    expect(codes(v)).toEqual(expect.arrayContaining(['seq.view.layer.extension', 'seq.render.theme.scheme']));
    const without = await v.renderSvg();
    const undo = await v.extend(DemoExtension);
    await settle(v);
    expect(types(v)).toEqual(['letters', 'demo:lollipops', 'ranges']);
    expect(issueLines(v)).toEqual(['No issues: every layer is drawn as its own type.']);
    // the SVG export draws the same registry: the lollipops' heads (rounded) in it, not in the export without it
    const withIt = await v.renderSvg();
    expect(withIt.issues.filter((i) => i.severity !== 'info')).toEqual([]);
    expect(withIt.svg).not.toEqual(without.svg);
    // each head keeps the label the fallback had (F20), in the export too; the Domains legend lists what the data has
    for (const site of ['S6', 'S15', 'S33', 'S46']) expect(withIt.svg).toContain(site);
    expect(recordDrawLists(v.frame!.lists, ['labels'])).toContain('"S6"');
    const domains = v.legends.legends.value.find((l) => l.title === 'Domains')!;
    expect(domains.entries.kind).toBe('swatches');
    expect(domains.entries.kind === 'swatches' ? domains.entries.items.map((i) => i.label) : []).toEqual([
      'transactivation',
      'proline-rich',
    ]);
    // the binding runs the extension's command
    expect(v.bindings.resolve('double-click', { alt: true }, 'cell')).toBe(FocusSites);
    expect(v.runCommand(FocusSites, {})).toBe('Zoomed to the sites (S6 to S46)');
    undo();
    await settle(v);
    expect(types(v)).toEqual(['letters', 'ranges', 'ranges']);
    expect(codes(v)).toContain('seq.view.layer.extension');
    expect((await v.renderSvg()).issues.map((i) => i.code)).toContain('seq.view.layer.extension');
    expect(v.commands.some((c) => c.id === FocusSites)).toBe(false);
  });

  /** The point on S6's head and the hover pass's rects after the pointer rests there. */
  async function hoverS6(ext: typeof DemoExtension) {
    const v = await open();
    await v.extend(ext);
    await settle(v);
    const layer = v
      .frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers))
      .find((l) => l.type.name === 'demo:lollipops')!;
    const s6 = layer.type.items!(layer.lay).find((i) => i.key === 'S6')!;
    // the head is at the extent's middle, at the stems' top: a point on it, whichever layout the story gives
    const at = [...Array(60).keys()]
      .map((k) => ({ x: (s6.x0 + s6.x1) / 2, y: k * 2 }))
      .find((p) => v.hit(p.x, p.y)?.element?.key === 'S6')!;
    expect(at).toBeDefined();
    expect(recordDrawLists(v.frame!.lists, ['hover'])).not.toContain('rect');
    const mark = { x: 0, y: 0, width: 0, height: 0 };
    expect(layer.type.mark!(layer.lay, v.hit(at.x, at.y)!.element!, mark)).toBe(true);
    v.input.move(at.x, at.y);
    await settle(v);
    const ink = (v.frame!.layout.tokens.colors['seq-hover'] >>> 8).toString(16).padStart(6, '0');
    const rects = recordDrawLists(v.frame!.lists, ['hover'])
      .split('\n')
      .filter((l) => l.startsWith('rect') && l.includes(`stroke #${ink}`));
    return { rects, mark };
  }

  it('outlines a hovered lollipop as a built-in item is: the hover ink, outside its head and round', async () => {
    const { rects, mark } = await hoverS6(DemoExtension);
    // the section's own outline (the head's) is one of them: it contains the head with the offset to spare
    const boxes = rects.map((r) => r.split(' ').slice(1, 5).map(Number));
    const around = boxes.find(
      ([x, y, w, h]) => x < mark.x && y < mark.y && x + w > mark.x + mark.width && y + h > mark.y + mark.height,
    );
    expect(around, rects.join('\n')).toBeDefined();
    // concentric with the head, not widened to the band: a head's width plus its offset and line, not the extent's
    expect(around![2]).toBeLessThan(mark.width + 6);
    expect(rects.some((r) => r.includes('r'))).toBe(true);
  });

  it('draws no outline for a type that draws its own hover (`ownHover`)', async () => {
    const own = defineViewerExtension({
      ...DemoExtension,
      registry: [
        { kind: 'layer-type', type: { ...LollipopsType, ownHover: true } },
        { kind: 'scheme', scheme: WarmScheme },
      ],
    });
    const { rects } = await hoverS6(DemoExtension);
    const { rects: none } = await hoverS6(own);
    expect(rects.length - none.length).toBe(1);
  });

  it('colors the domains by the warm scheme only where it is composed', async () => {
    const spec = extensibilitySpec();
    const color = (ctx: ReturnType<typeof createSeqContext>) =>
      ctx.renderSvg(spec, { size: { width: 600, height: 'auto' } }).then((r) => r.issues.map((i) => i.code));
    expect(await color(createSeqContext())).toContain('seq.render.theme.scheme');
    expect(await color(createSeqContext({ extensions: [DemoExtension] }))).not.toContain('seq.render.theme.scheme');
    expect(DomainsLayer).toMatch(/domains$/);
    // by category, never by the order met: TAD1 (transactivation) and TAD2 (proline-rich) their own colors
    const { svg } = await createSeqContext({ extensions: [DemoExtension] }).renderSvg(spec, {
      size: { width: 600, height: 'auto' },
    });
    const fills = [...svg.matchAll(/<path d="M [^"]*" fill="(#[0-9a-f]{6})"/g)].map((m) => m[1]);
    expect(fills.filter((f) => f === '#d9480f' || f === '#f59f00')).toEqual(['#d9480f', '#f59f00']);
  });

  it('checks its params where it is composed, and carries them where not', () => {
    const spec = extensibilitySpec();
    const set = (root: Node, at: string, key: string, value: unknown) => {
      const r = applyPatch(root, { ops: [{ op: 'set', at: at.split('/').filter(Boolean), key, value }] } as never);
      if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
      return r.root;
    };
    const bad = set(
      set(spec, '/views/main/protein/sites/lollipops', 'head', 40),
      '/views/main/protein',
      'demo:stems',
      'low',
    );
    expect(createSeqContext({ extensions: [DemoExtension] }).check(spec)).toEqual([]);
    expect(
      createSeqContext({ extensions: [DemoExtension] })
        .check(bad)
        .map((i) => i.code),
    ).toEqual(['repr.schema.param.enum', 'repr.schema.param.range']);
    expect(createSeqContext().check(bad)).toEqual([]);
  });
});
