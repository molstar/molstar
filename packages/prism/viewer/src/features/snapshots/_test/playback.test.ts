// @vitest-environment jsdom

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, until } from '../../../_test/harness.js';
import type { RunRecord } from '../../../api/run.js';
import type { PrismViewer } from '../../../viewer.js';
import { SnapshotsFeature } from '../index.js';
import { createMarkdownHost } from '../markdown.js';
import { createSnapshotStepper, openPlaybackOnLoad } from '../model.js';

// The playback card's descriptions (Mol*'s markdown commands) and its opening on a load (parity.md §4), headless.

const viewer = () => headlessViewer({ extensions: [SnapshotsFeature], files: [cbs()] });
const snap = (v: PrismViewer) => v.api.snapshot!;
const current = (v: PrismViewer) => v.plugin.managers.snapshot.state.current;

describe('snapshot.link', () => {
  it("runs a description's command through the plugin's markdown extensions, with no history step", async () => {
    const v = await viewer();
    const a = await snap(v).add({ name: 'A', key: 'a' });
    expectOk(a);
    expectOk(await snap(v).add({ name: 'B', key: 'b' }));
    const steps = v.managers.history.steps.value.length;
    const r = await snap(v).link({ href: '!apply-snapshot=a' });
    expectOk(r);
    expect(r.did).toBe('Ran apply-snapshot.');
    await until(() => current(v) === (a.value as { id: string }).id);
    expect(v.managers.history.steps.value.length).toBe(steps);
    // Hover events go to the extensions too (highlight-refs lights its refs up while the pointer is on the link).
    expectOk(await snap(v).link({ href: '!highlight-refs=nothing', event: 'enter' }));
    expectOk(await snap(v).link({ href: '!highlight-refs=nothing', event: 'leave' }));
  });

  it('refuses a link that is not a command', async () => {
    const v = await viewer();
    expect(await snap(v).link({ href: 'https://molstar.org' })).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('the descriptions’ markdown host', () => {
  it('claims command links only, and runs them through snapshot.link', async () => {
    const v = await viewer();
    const runs = recordRuns(v);
    const host = createMarkdownHost(v.api, snap(v).markdown);
    expect(host.link?.('#key')).toBeUndefined();
    expect(host.link?.('https://molstar.org')).toBeUndefined();
    host.link?.('!center-camera')?.('click');
    await until(() => runs.some((r) => r.operation === 'snapshot.link'));
    expect(runs.find((r) => r.operation === 'snapshot.link')?.values).toMatchObject({ href: '!center-camera' });
  });

  it("renders v5's inline commands and leaves other images to the plugin's resolvers", async () => {
    const v = await viewer();
    const host = createMarkdownHost(v.api, snap(v).markdown);
    const swatch = host.image?.('!color-swatch=red', '') as HTMLElement;
    expect(swatch.style.backgroundColor).toBe('red');
    const palette = host.image?.('!color-palette-name=viridis', '') as HTMLElement;
    expect(palette.style.background).toContain('gradient');
    expect(host.image?.('picture.png', 'P')).toBeUndefined();
    // The audio player is the card's own bar, not media in the text.
    expect(host.media).toBeUndefined();
  });

  it('gives the stepper the host and whether the description is plain text', async () => {
    const v = await viewer();
    expectOk(await snap(v).add({ name: 'A', description: 'line one\nline two' }));
    const host = createMarkdownHost(v.api, snap(v).markdown);
    const stepper = createSnapshotStepper(v.api, snap(v), host);
    expect(stepper.markdown).toBe(host);
    expect(stepper.plainText?.value).toBe(false);
    const entry = v.plugin.managers.snapshot.state.entries.get(0)!;
    v.plugin.managers.snapshot.update(entry, { descriptionFormat: 'plaintext' });
    expect(stepper.plainText?.value).toBe(true);
  });
});

describe('the playback card on a load', () => {
  const fake = (v: PrismViewer) => {
    const started = new Subject<RunRecord>();
    const finished = new Subject<RunRecord>();
    const open = vi.fn(async () => ({ ok: true }));
    const api = { events: { started, finished }, view: { open } } as unknown as Parameters<
      typeof openPlaybackOnLoad
    >[0];
    openPlaybackOnLoad(api, snap(v));
    let id = 0;
    const run = async (operation: string, work: () => Promise<unknown>, ok = true) => {
      const record = { id: ++id, operation, title: operation, values: {}, origin: { kind: 'person' }, started: 0 };
      started.next(record as unknown as RunRecord);
      await work();
      finished.next({ ...record, result: { ok } } as unknown as RunRecord);
    };
    return { open, run };
  };

  it('opens once a load brings two or more snapshots, and not for a load that leaves the list', async () => {
    const v = await viewer();
    const { open, run } = fake(v);
    await run('load.mvs', async () => {
      await snap(v).add({ name: 'A' });
      await snap(v).add({ name: 'B' });
    });
    expect(open).toHaveBeenCalledTimes(1);
    expect(open).toHaveBeenCalledWith({ card: 'snapshots.overlay.playback' });
    await run('load.entry', async () => {});
    // A snapshot the person adds is not a load.
    await run('snapshot.add', async () => void (await snap(v).add({ name: 'C' })));
    expect(open).toHaveBeenCalledTimes(1);
  });

  it('opens when the views come after a load that brought a story, unless a load is still running', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    await snap(v).add({ name: 'B' });
    const open = vi.fn(async () => ({ ok: true }));
    const events = { started: new Subject<RunRecord>(), finished: new Subject<RunRecord>() };
    const running = (operations: string[]) => ({ value: operations.map((operation) => ({ operation })) });
    openPlaybackOnLoad({ events, view: { open }, running: running(['load.mvs']) } as never, snap(v));
    expect(open).not.toHaveBeenCalled();
    // That load ends (its snapshots were in before the views came): the card opens.
    const record = { id: 7, operation: 'load.mvs', title: '', values: {}, origin: { kind: 'person' }, started: 0 };
    events.finished.next({ ...record, result: { ok: true } } as never);
    expect(open).toHaveBeenCalledTimes(1);
    open.mockClear();
    openPlaybackOnLoad({ events, view: { open }, running: running([]) } as never, snap(v));
    expect(open).toHaveBeenCalledWith({ card: 'snapshots.overlay.playback' });
  });

  it('opens for one snapshot with a description; not for one without, nor for a failed load', async () => {
    const v = await viewer();
    const { open, run } = fake(v);
    await run('session.open', async () => void (await snap(v).add({ name: 'A' })));
    expect(open).not.toHaveBeenCalled();
    await run(
      'load.mvs',
      async () => {
        await snap(v).clear();
        await snap(v).add({ name: 'B', description: 'The **story**.' });
      },
      false,
    );
    expect(open).not.toHaveBeenCalled();
    await run('load.mvs', async () => {
      await snap(v).clear();
      await snap(v).add({ name: 'C', description: 'The **story**.' });
    });
    expect(open).toHaveBeenCalledTimes(1);
  });
});
