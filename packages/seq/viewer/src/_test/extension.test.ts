/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A viewer's extensions (_spec/README.md "Composition", Extensions): composed, an extension's layer type is drawn on
// the canvas and its behaviors, commands and bindings are the viewer's; registered late (`extend`), a layer drawn as
// its fallback is drawn as its type, and as the fallback again after the undo, its parts gone with it.

import type { Node } from '@molstar/core/repr/tree';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { type LayerType, RangesType, TracksLayout } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { type HeadersBehavior, HeadersProvider } from '../behaviors/headers.js';
import { HoverProvider } from '../behaviors/providers.js';
import { defineViewerExtension } from '../extension.js';
import type { SeqViewerExtension, SeqViewerSpec } from '../viewer-spec.js';
import { createSeqViewer, type SeqViewer, type ViewerCommand } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

/** Ranges drawn as one bar each (a ranges layer under another name): what the frame shows is its name. */
const Bars: LayerType = { ...(RangesType as LayerType), name: 'demo:bars', shape: 'ranges', fallback: 'ranges' };

function spec(): Node {
  const b = createSeqSpec({ title: 'ext', extensions: [{ id: 'demo', version: '0.1.0' }] });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const a = b.data.ranges('a', { target: p, rows: { start: new Int32Array([2]), end: new Int32Array([9]) } });
  b.views.view('main', (v) => {
    v.section('s', { of: p })
      .track('t')
      .layer('l', { type: 'demo:bars', fallback: 'ranges', data: a } as never);
  });
  return b.build();
}

function demo() {
  const disposed = vi.fn();
  const ran = vi.fn();
  const command = {
    id: 'demo:hello',
    title: 'Hello',
    category: 'Demo',
    run: () => (ran(), { message: 'hello' }),
  } as unknown as ViewerCommand;
  const base: SeqViewerExtension = {
    id: 'demo',
    version: '0.1.0',
    requires: '0.11',
    registry: [{ kind: 'layer-type', type: Bars }],
    behaviors: [{ provider: { name: 'demo:behavior', create: () => ({ dispose: disposed }) } }],
    commands: [command],
    bindings: [{ gesture: 'double-click', modifiers: ['alt'], target: 'any', action: 'demo:hello' }],
  };
  const ext = defineViewerExtension(base);
  /** Another descriptor of it, made again with some fields changed. */
  const with_ = (over: Partial<SeqViewerExtension>) => defineViewerExtension({ ...base, ...over });
  return { ext, with_, disposed, ran };
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(extra: Partial<SeqViewerSpec> = {}) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames ??= installFrames();
  const v = createSeqViewer({
    registry: [
      { kind: 'format-runtime', runtime: SeqRuntime },
      { kind: 'layer-type', type: RangesType },
      { kind: 'layout', layout: TracksLayout },
    ],
    behaviors: [{ provider: HoverProvider }],
    ...extra,
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const drawn = (v: SeqViewer) =>
  v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers.map((l) => l.type.name)));
const alt = { alt: true };

describe('a viewer’s extensions', () => {
  it('composed: its type drawn, its behaviors, commands and bindings the viewer’s', async () => {
    const { ext, ran } = demo();
    const v = await make({ extensions: [ext] });
    expect(drawn(v)).toEqual(['demo:bars']);
    expect(v.view.resolved.value!.issues).toEqual([]);
    expect(v.behaviors.has('demo:behavior')).toBe(true);
    expect(v.bindings.resolve('double-click', alt, 'cell')).toBe('demo:hello');
    v.runCommand('demo:hello', {});
    expect(ran).toHaveBeenCalledOnce();
  });

  it('registered late: drawn as its type, then as its fallback again after the undo, its parts gone', async () => {
    const { ext, disposed } = demo();
    const v = await make();
    expect(drawn(v)).toEqual(['ranges']);
    expect(v.view.resolved.value!.issues.map((i) => i.code)).toEqual(['seq.view.layer.extension']);
    expect(v.bindings.resolve('double-click', alt, 'cell')).not.toBe('demo:hello');
    const before = v.context.version;
    const undo = await v.extend(ext);
    const again = await v.extend(ext);
    await v.settle();
    frames!.flush();
    expect(drawn(v)).toEqual(['demo:bars']);
    expect(v.context.version).toBe(before + 1);
    expect([v.behaviors.has('demo:behavior'), v.commands.some((c) => c.id === 'demo:hello')]).toEqual([true, true]);
    expect(v.bindings.resolve('double-click', alt, 'cell')).toBe('demo:hello');
    // counted: the first undo leaves it, the last takes it away
    undo();
    undo();
    expect(drawn(v)).toEqual(['demo:bars']);
    again();
    await v.settle();
    frames!.flush();
    expect(drawn(v)).toEqual(['ranges']);
    expect(disposed).toHaveBeenCalledOnce();
    expect([v.behaviors.has('demo:behavior'), v.commands.some((c) => c.id === 'demo:hello')]).toEqual([false, false]);
    expect(v.bindings.resolve('double-click', alt, 'cell')).not.toBe('demo:hello');
    expect(v.context.extensions).toEqual([]);
  });

  it('its settings compose on a late extend: the composition’s pairs under its namespace read, gone with its undo', async () => {
    const size = { key: 'demo.size' };
    const v = await make({ config: [[size, 7]] as never });
    expect(v.context.renderConfig.settings?.['demo.size']).toBeUndefined();
    const { with_ } = demo();
    const undo = await v.extend(with_({ config: [size] }));
    expect(v.context.renderConfig.settings?.['demo.size']).toBe(7);
    undo();
    expect(v.context.renderConfig.settings?.['demo.size']).toBeUndefined();
  });

  it('refuses what it cannot register late, and one that does not check, leaving the viewer as it was', async () => {
    const v = await make();
    const { with_ } = demo();
    expect(() => with_({ requires: '0.4' })).toThrow(/seq.context.extension.version/);
    const late = [
      with_({ registry: [{ kind: 'format-runtime', runtime: { ...SeqRuntime, name: 'demo:x' } }] }),
      // a part over the canvas would never be attached on a viewer already mounted
      with_({ behaviors: [{ provider: { name: 'demo:part', create: () => ({ attach: () => () => {} }) } }] }),
    ];
    for (const ext of late) await expect(v.extend(ext)).rejects.toThrow(/composed, not registered late/);
    expect([v.context.extensions, v.behaviors.has('demo:part'), drawn(v)]).toEqual([[], false, ['ranges']]);
    // an extension without its viewer maker loses no part silently
    await expect(v.extend({ ...with_({}), composeViewer: undefined })).rejects.toThrow(/defineViewerExtension/);
  });

  it('its header entries registered late go on the viewer’s headers, and with its undo; one naming no command refused', async () => {
    const v = await make({ behaviors: [{ provider: HoverProvider }, { provider: HeadersProvider }] });
    const headers = v.behavior<HeadersBehavior>('headers')!;
    const { with_ } = demo();
    const act = { name: 'demo:act', on: 'track' as const, label: 'Hello', command: 'demo:hello' };
    const undo = await v.extend(with_({ headerActions: [act] }));
    expect(headers.entries.includes(act)).toBe(true);
    undo();
    expect(headers.entries.includes(act)).toBe(false);
    // an entry naming a command nobody has: refused, the viewer as it was
    const bad = with_({ headerActions: [{ ...act, name: 'demo:bad', command: 'demo:none' }] });
    await expect(v.extend(bad)).rejects.toThrow(/no command demo:none/);
    expect([v.context.extensions, v.commands.some((c) => c.id === 'demo:hello')]).toEqual([[], false]);
  });

  it('a failed registration takes away only what it made, never the viewer’s own behaviors', async () => {
    const v = await make();
    const { with_ } = demo();
    const hover = v.behaviors.get('hover');
    const clash = with_({
      behaviors: [
        { provider: { name: 'demo:first', create: () => ({}) } },
        { provider: { name: 'hover', create: () => ({}) } },
      ],
    });
    await expect(v.extend(clash)).rejects.toThrow(/listed twice/);
    expect([v.behaviors.get('hover') === hover, v.behaviors.has('demo:first'), v.context.extensions]).toEqual([
      true,
      false,
      [],
    ]);
  });
});
