/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { firstValueFrom } from 'rxjs';
import { mockInspector } from '../../testing/mock-inspector.js';
import { findLimit } from '../model.js';

describe('InspectorModel', () => {
  it('selects the latest row until one is chosen; move steps through the rows and stops at the ends', async () => {
    const model = mockInspector();
    const rows = await firstValueFrom(model.rows);
    expect(rows[0].key).toBe('turn:turn-1');
    expect(rows[1].key).toBe('req:turn-1');
    expect(await firstValueFrom(model.selected)).toBe(rows.at(-1)!.key);
    model.move(-1);
    expect(await firstValueFrom(model.selected)).toBe(rows.at(-2)!.key);
    model.move(5);
    expect(await firstValueFrom(model.selected)).toBe(rows.at(-1)!.key);
    model.select('turn:turn-1');
    model.move(-1);
    expect(await firstValueFrom(model.selected)).toBe('turn:turn-1');
  });

  it('marks nested steps (a tool’s helper) and keeps rows across re-reads of the same trace', async () => {
    const model = mockInspector();
    const rows = await firstValueFrom(model.rows);
    const nested = rows.filter((r) => r.kind === 'step' && r.nested);
    expect(nested).toHaveLength(1);
    expect(await firstValueFrom(model.rows)).toBe(rows);
  });

  it('applies a page focus: a view, a turn (the latest when unknown), a step, a tool, a group', async () => {
    const model = mockInspector();
    expect(model.focus('tools')).toBe(true);
    expect(model.view.value).toBe('tools');
    expect(model.focus('turn:turn-2')).toBe(true);
    expect(model.view.value).toBe('trace');
    expect(await firstValueFrom(model.selected)).toBe('turn:turn-2');
    expect(model.focus('turn:not-in-the-trace')).toBe(true);
    expect(await firstValueFrom(model.selected)).toBe('turn:turn-3');
    expect(model.focus('step:s3')).toBe(true);
    expect(await firstValueFrom(model.selected)).toBe('step:s3');
    expect(model.focus('tool:bounds.get')).toBe(true);
    expect(model.view.value).toBe('tools');
    expect(await firstValueFrom(model.selectedTool)).toBe('tool:bounds.get');
    expect(model.focus('group:volumes')).toBe(true);
    expect(await firstValueFrom(model.selectedTool)).toBe('group:volumes');
    expect(model.focus('limits')).toBe(true);
    expect(model.view.value).toBe('limits');
    expect(model.focus('nonsense')).toBe(false);
    expect(model.focus('tool:')).toBe(false);
    expect(model.view.value).toBe('limits');
  });

  it('opens a tool from a tool call: the Tools view, the tool selected', async () => {
    const model = mockInspector();
    model.openTool('selection.set');
    expect(model.view.value).toBe('tools');
    const row = await firstValueFrom(model.toolRow);
    expect(row?.kind === 'tool' && row.tool.name).toBe('selection.set');
  });

  it('sets a limit, refuses what the row cannot take, and resets to the default', () => {
    const model = mockInspector();
    const value = (key: string) => findLimit(model.limits.value, key)?.value;
    expect(model.setLimit('modelCalls', 40)).toBe(true);
    expect(value('modelCalls')).toBe(40);
    expect(model.setLimit('modelCalls', 1000)).toBe(false);
    expect(model.setLimit('modelCalls', 0)).toBe(false);
    expect(model.setLimit('modelCalls', Number.NaN)).toBe(false);
    expect(value('modelCalls')).toBe(40);
    expect(model.setLimit('modelCalls', undefined)).toBe(true);
    expect(value('modelCalls')).toBe(24);
    // Rounded to the row's step.
    expect(model.setLimit('output', 40_400)).toBe(true);
    expect(value('output')).toBe(40_000);
    // A choice: only its values.
    expect(model.setLimit('thinking', 5)).toBe(false);
    expect(model.setLimit('thinking', 3)).toBe(true);
    expect(value('thinking')).toBe(3);
    // Information, not a limit.
    expect(model.setLimit('credits.turn', 2)).toBe(false);
    expect(model.setLimit('no-such-limit', 2)).toBe(false);
  });

  it('has no rows before the agent has run', async () => {
    const model = mockInspector({ empty: true });
    expect(await firstValueFrom(model.rows)).toEqual([]);
    expect(await firstValueFrom(model.selected)).toBeUndefined();
    model.move(1);
    expect(model.focus('turn:turn-1')).toBe(true);
    expect(await firstValueFrom(model.selected)).toBeUndefined();
  });
});
