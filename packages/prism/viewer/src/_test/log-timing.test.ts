/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from './harness.js';

// How long things took, in the Log only (journal.md §3): a run's message carries its duration as the log's detail, and
// the plugin's transforms are logged with theirs, the latest 100.

describe('timings in the log', () => {
  it("gives a person's run its duration in the log, not in the status line's text", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const person = v.api.as({ kind: 'person' });
    expectOk(await person.styles.apply({ style: 'spacefill' }));
    await until(() => !!v.managers.notifications.current.value);
    const shown = v.managers.notifications.current.value!;
    expect(shown.text).not.toMatch(/\d+ ms|\d s$/);
    expect(shown.detail).toMatch(/^\d+ ms$|^\d+(\.\d)? s$/);
  });

  it("logs the plugin's transforms with their time, the latest 100", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const transforms = () => v.managers.notifications.history.value.filter((m) => m.source === 'Transforms');
    await until(() => transforms().length > 0);
    expect(transforms()[0].text).toMatch(/^(Created|Updated) /);
    expect(transforms()[0].detail).toMatch(/^[\d.]+ m?s$/);
    expect(v.managers.notifications.current.value?.source).toBeUndefined();
    for (let i = 0; i < 6; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
    await Promise.resolve();
    expect(transforms()).toHaveLength(100);
  });
});
