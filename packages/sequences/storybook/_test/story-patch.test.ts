/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' `storyPatch` (lib/seqstar/storybook/src/_stories/support/seqstar-patch.ts): a patch applied through it keeps the issues
// that came back (warning or error) for the banner and logs them; a patch that applies cleanly leaves nothing.

import type { Patch } from '@molstar/seq-core/format';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { pf00042 } from '../../core/_test/model/alignment-fixtures.js';
import { keyedFamilyDoc, MEMBERS_LAYER, MODES, modePatch } from '../../viewer/_test/keyed-family.js';
import { opened } from '../../viewer/_test/layout-fixtures.js';
import { patchNotices, storyPatch, storyPropose } from '../src/_stories/support/seqstar-patch.js';

afterEach(() => vi.restoreAllMocks());
const doc = () => keyedFamilyDoc(pf00042());
const reps = MODES.find((m) => m.key === 'reps90') ?? MODES[0];
/** What the story's first version sent: the track's ref, which no layer has. */
const wrong: Patch = {
  format: 'seq_view_spec_patch',
  ops: [{ op: 'set', ref: 'members', params: { order: reps.order, filter: reps.filter } }],
};

describe('storyPatch', () => {
  it('keeps and logs the issues of a patch that names the wrong ref, and still returns the result', async () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const t = await opened(doc());
    const r = await storyPatch(t.ctx, wrong, 'Representatives at 90%');
    const notices = patchNotices(t.ctx).list;
    expect(notices).toHaveLength(1);
    expect(notices[0].label).toBe('Representatives at 90%');
    expect(notices[0].issues.length).toBeGreaterThan(0);
    expect(notices[0].issues.every((i) => i.severity === 'warning' || i.severity === 'error')).toBe(true);
    expect(notices[0].issues.map((i) => i.message).join(' ')).toContain('order');
    expect(r.issues.length).toBeGreaterThanOrEqual(notices[0].issues.length);
    expect(warn).toHaveBeenCalledTimes(1);
    expect(String(warn.mock.calls[0][0])).toContain('Representatives at 90%');
  });

  it('a clean patch leaves no notice and logs nothing', async () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const t = await opened(doc());
    const r = await storyPatch(t.ctx, modePatch(reps) as Patch, reps.label);
    expect(r.ok).toBe(true);
    expect(patchNotices(t.ctx).list).toEqual([]);
    expect(warn).not.toHaveBeenCalled();
    expect(MEMBERS_LAYER).toBe('members.rows');
  });

  it('notices are listed, can be dismissed one by one or cleared, and tell their subscribers', async () => {
    vi.spyOn(console, 'warn').mockImplementation(() => {});
    const t = await opened(doc());
    const store = patchNotices(t.ctx);
    let told = 0;
    const off = store.subscribe(() => told++);
    await storyPatch(t.ctx, wrong, 'one');
    await storyPatch(t.ctx, wrong, 'two');
    expect(store.list.map((n) => n.label)).toEqual(['one', 'two']);
    store.dismiss(store.list[0].id);
    expect(store.list.map((n) => n.label)).toEqual(['two']);
    store.clear();
    expect(store.list).toEqual([]);
    expect(told).toBe(4);
    off();
    await storyPatch(t.ctx, wrong, 'three');
    expect(told).toBe(4);
    // one list per context
    const other = await opened(doc());
    expect(patchNotices(other.ctx).list).toEqual([]);
  });

  it('storyPropose reports the issues of a proposed patch the same way', async () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const t = await opened(doc());
    const r = await storyPropose(t.ctx, wrong, 'Style', 'host:story');
    expect(r.ok).toBe(false);
    expect(patchNotices(t.ctx).list).toHaveLength(1);
    expect(warn).toHaveBeenCalledTimes(1);
  });
});
