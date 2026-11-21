/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { DemoExtension, DomainsLayer, extensibilityComposition, extensibilitySpec, issueLines } from '../compose.js';
import { FocusSites } from '../demo-extension.js';

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
    // the binding runs the extension's command
    expect(v.bindings.resolve('double-click', { alt: true }, 'cell')).toBe(FocusSites);
    expect(v.runCommand(FocusSites, {})).toBe('The sites are shown');
    undo();
    await settle(v);
    expect(types(v)).toEqual(['letters', 'ranges', 'ranges']);
    expect(codes(v)).toContain('seq.view.layer.extension');
    expect((await v.renderSvg()).issues.map((i) => i.code)).toContain('seq.view.layer.extension');
    expect(v.commands.some((c) => c.id === FocusSites)).toBe(false);
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
