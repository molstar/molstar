/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs, headlessViewer, until } from '../../../_test/harness.js';
import type { AnnotationProvider } from '../../../api/annotations.js';
import { type AnnotationOffer, createAnnotationOffer, type HomeModel } from '../../../model/home.js';
import type { ViewerExtension } from '../../../spec.js';
import type { PrismViewer } from '../../../viewer.js';
import { HomeTab } from '../home-tab.js';

// Home's result line after a load (parity.md §6a): an extension's annotation a new structure carries is offered once;
// the built-ins (B-factors, secondary structure) never are. A stand-in extension plays model-archive's pLDDT.

const Score: AnnotationProvider = {
  id: 'test-score',
  title: 'test scores',
  description: 'A stand-in for pLDDT',
  isApplicable: (s) => s.polymerResidueCount > 0,
  apply: (ctx, structure) => ctx.api.representation.color({ target: 'polymer', theme: 'uncertainty', structure }),
};
const Scores: ViewerExtension = { id: 'test-scores', title: 'Scores', annotations: [Score] };

afterEach(() => {
  document.body.innerHTML = '';
});

const structures = (v: PrismViewer) => v.sources.structures;

/** Follows the offer as a view would; returns the latest value and the unsubscribe. */
function watch(model: ReturnType<typeof createAnnotationOffer>) {
  const seen: { value: AnnotationOffer | undefined } = { value: undefined };
  const sub = model.offer.subscribe((o) => {
    seen.value = o;
  });
  return { seen, stop: () => sub.unsubscribe() };
}

describe("Home's annotation offer", () => {
  it("offers an extension's annotation of a structure that loads, and Show applies it as one step", async () => {
    const v = await headlessViewer({ extensions: [Scores] });
    const model = createAnnotationOffer(v.api.as({ kind: 'person' }), structures(v));
    const { seen, stop } = watch(model);
    expect(seen.value).toBeUndefined();
    expect((await v.api.load.files({ files: [cbs()] })).ok).toBe(true);
    await until(() => !!seen.value);
    expect(seen.value).toMatchObject({ annotation: 'test-score', title: 'test scores' });
    expect(seen.value?.structure).toBe(structures(v).value[0].ref);

    const steps = v.managers.history.steps.value.length;
    model.accept();
    expect(seen.value).toBeUndefined();
    await until(() => v.managers.history.steps.value.length === steps + 1);
    stop();
  });

  it('never offers the built-ins, and a dismissed offer does not come back', async () => {
    const plain = await headlessViewer({ files: [cbs()] });
    const none = createAnnotationOffer(plain.api.as({ kind: 'person' }), structures(plain));
    const a = watch(none);
    // B-factors and secondary structure are what 1CBS carries: listed, not offered.
    const listed = await plain.api.annotations.list({});
    expect(listed.ok && JSON.stringify(listed.value)).toMatch(/secondary-structure/);
    await new Promise((r) => setTimeout(r, 50));
    expect(a.seen.value).toBeUndefined();
    a.stop();

    const v = await headlessViewer({ files: [cbs()], extensions: [Scores] });
    const model = createAnnotationOffer(v.api.as({ kind: 'person' }), structures(v));
    const b = watch(model);
    await until(() => !!b.seen.value);
    model.dismiss();
    expect(b.seen.value).toBeUndefined();
    // Watched again (Home shown again): the structure was asked already.
    b.stop();
    const c = watch(model);
    await new Promise((r) => setTimeout(r, 50));
    expect(c.seen.value).toBeUndefined();
    c.stop();
  });

  it("still offers the first structure's annotation when another arrives while it is being listed", async () => {
    const list = observableState<readonly { ref: string; label: string }[]>([{ ref: 'a', label: 'A' }]);
    let finishA!: () => void;
    const lists: Record<string, () => Promise<unknown>> = {
      a: () =>
        new Promise((resolve) => {
          finishA = () => resolve({ ok: true, value: [{ annotations: [{ id: 'plddt', title: 'pLDDT' }] }] });
        }),
      b: async () => ({ ok: true, value: [{ annotations: [{ id: 'secondary-structure', title: 'Secondary' }] }] }),
    };
    const api = {
      area: () => true,
      annotations: { list: ({ structure }: { structure: string }) => lists[structure](), providers: { value: [] } },
    };
    const model = createAnnotationOffer(api as never, list.asReadonly);
    const { seen, stop } = watch(model);
    await until(() => !!finishA);
    // B arrives while A is listed and carries nothing; A's result must not be dropped.
    list.set([
      { ref: 'a', label: 'A' },
      { ref: 'b', label: 'B' },
    ]);
    await new Promise((r) => setTimeout(r, 10));
    finishA();
    await until(() => !!seen.value);
    expect(seen.value).toMatchObject({ structure: 'a', annotation: 'plddt' });
    stop();
  });

  it('shows on Home as a line with Show and Dismiss', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [Scores] });
    const model = createAnnotationOffer(v.api.as({ kind: 'person' }), structures(v));
    const home: HomeModel = {
      actions: observableState([]).asReadonly,
      formats: [],
      openFiles: () => {},
      annotationOffer: model,
    };
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <HomeTab home={home} />, host);
    await until(() => host.textContent?.includes('has test scores') ?? false);
    const dismiss = host.querySelector<HTMLButtonElement>('button[data-tip="Dismiss"]')!;
    dismiss.click();
    await until(() => !host.textContent?.includes('has test scores'));
  });
});
