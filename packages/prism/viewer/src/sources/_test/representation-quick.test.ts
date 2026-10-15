/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, headlessViewer, until } from '../../_test/harness.js';
import type { SceneNode } from '../../model/scene.js';
import type { PrismViewer } from '../../viewer.js';

// A structure representation's quick settings (model/representation-quick.ts) over a headless viewer: the common types
// and themes that apply, a theme or one color picked, another type keeping the color and opacity, the size as a multiple
// of the type's default, a drag one undo step.

async function cartoon(file = cbs()) {
  const v = await headlessViewer();
  await v.api.load.files({ files: [file] });
  await until(() => v.views.scene.roots.value.length > 0);
  const find = (nodes: readonly SceneNode[]): SceneNode | undefined => {
    for (const n of nodes) {
      if (n.kind === 'representation' && n.label.value === 'Cartoon') return n;
      const inner = find(n.children.value);
      if (inner) return inner;
    }
    return undefined;
  };
  await until(() => !!find(v.views.scene.roots.value));
  const row = find(v.views.scene.roots.value)!;
  return { v, row, quick: v.views.scene.inspect(row).quick! };
}

const params = (v: PrismViewer, ref: string) => v.plugin.state.data.cells.get(ref)!.transform.params as any;
const idle = (v: PrismViewer) => until(() => v.api.running.value.length === 0);

describe('the quick settings of a structure representation', () => {
  it('offers the common types and themes that apply, the current ones marked', async () => {
    const { quick } = await cartoon();
    expect(quick.types.value.map((t) => t.id)).toContain('cartoon');
    expect(quick.types.value.map((t) => t.id)).toContain('molecular-surface');
    expect(quick.themes.value.map((t) => t.id)).toContain('uniform');
    // An X-ray structure: no pLDDT.
    expect(quick.themes.value.map((t) => t.id)).not.toContain('plddt-confidence');
    expect(quick.type.value).toBe('cartoon');
    expect(quick.opacity.value).toBe(1);
    expect(quick.size.value).toBe(1);
  });

  it('offers no pLDDT for an experimental structure whose file does not say how it was determined', async () => {
    // Crambin's file predates the determination methodology Mol*'s pLDDT theme looks for; it has an X-ray method.
    const { quick } = await cartoon(crambin());
    expect(quick.themes.value.map((t) => t.id)).not.toContain('plddt-confidence');
  });

  it('sets a theme, one color, another type keeping the color and opacity, and the size', async () => {
    const { v, row, quick } = await cartoon();
    quick.change({ theme: 'element-symbol', color: undefined });
    quick.settle();
    await until(() => params(v, row.id).colorTheme.name === 'element-symbol');

    quick.change({ theme: 'uniform', color: 0x3362b2 });
    quick.settle();
    await until(() => params(v, row.id).colorTheme.name === 'uniform');
    expect(params(v, row.id).colorTheme.params.value).toBe(0x3362b2);
    expect(quick.color.value).toBe(0x3362b2);

    quick.change({ opacity: 0.5 });
    quick.settle();
    await until(() => params(v, row.id).type.params.alpha === 0.5);

    quick.change({ type: 'molecular-surface' });
    quick.settle();
    await until(() => params(v, row.id).type.name === 'molecular-surface');
    expect(params(v, row.id).type.params.alpha).toBe(0.5);
    expect(params(v, row.id).colorTheme).toMatchObject({ name: 'uniform', params: { value: 0x3362b2 } });

    quick.change({ type: 'cartoon' });
    quick.settle();
    await until(() => params(v, row.id).type.name === 'cartoon');
    const base = params(v, row.id).type.params.sizeFactor;
    quick.change({ size: 2 });
    quick.settle();
    await until(() => params(v, row.id).type.params.sizeFactor === Math.round(base * 2 * 1000) / 1000);
    await until(() => quick.size.value === 2);
  });

  it('makes a slider drag one undo step', async () => {
    const { v, row, quick } = await cartoon();
    await idle(v);
    const steps = v.api.history.steps.value.length;
    for (const a of [0.9, 0.8, 0.7, 0.6]) quick.change({ opacity: a });
    expect(quick.opacity.value).toBe(0.6);
    quick.settle();
    await until(() => params(v, row.id).type.params.alpha === 0.6);
    await until(() => v.api.history.steps.value.length === steps + 1);
  });
});
