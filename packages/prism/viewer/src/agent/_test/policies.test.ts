/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '@molstar/prism-ui/core';
import { describe, expect, it } from 'vitest';
import type { OperationDefinition } from '../../api/area.js';
import { allowedAddress, createAgentPolicy, createDomainRules } from '../policies.js';

// What goes without asking (agent.md §4).

const op = (o: Partial<OperationDefinition<any, any>>) =>
  ({ id: 'x.y', title: 'X', params: {}, run: async () => ({}) as never, ...o }) as OperationDefinition<any, any>;

describe('the agent policy', () => {
  it('lets changes go while the turn can be undone, and asks for removals and uploads', () => {
    const p = createAgentPolicy({ checkpointing: observableState(true) });
    expect(p.approves(op({}), {})).toBe(true);
    expect(p.approves(op({ approval: 'destructive' }), {})).toBe(false);
    expect(p.approves(op({ approval: 'upload' }), {})).toBe(false);
  });

  it('asks once for the turn when undo is off', () => {
    const checkpointing = observableState(false);
    const p = createAgentPolicy({ checkpointing });
    expect(p.approves(op({}), {})).toBe(false);
    expect(p.approves(op({ history: false }), {})).toBe(true);
    p.allowedForTurn();
    expect(p.approves(op({}), {})).toBe(true);
    p.beginTurn();
    expect(p.approves(op({}), {})).toBe(false);
  });

  it('reaches the known services and allowed domains by https; others ask', () => {
    const rules = createDomainRules();
    const p = createAgentPolicy({ checkpointing: observableState(true), rules });
    const net = op({ approval: 'network' });
    expect(p.approves(net, { url: 'https://files.rcsb.org/download/1cbs.cif' })).toBe(true);
    expect(p.approves(net, { url: 'https://example.org/x.cif' })).toBe(false);
    expect(p.approves(net, { url: 'http://files.rcsb.org/x.cif' })).toBe(false);
    rules.allow('example.org', 'session');
    expect(p.approves(net, { url: 'https://data.example.org/x.cif' })).toBe(true);
    expect(allowedAddress('https://evil-rcsb.org/x', rules)).toBe(false);
  });

  it('keeps the browser rules through its host', () => {
    let saved: readonly string[] = [];
    const rules = createDomainRules(['kept.org'], (b) => {
      saved = b;
    });
    rules.allow('session.org', 'session');
    rules.allow('browser.org', 'browser');
    expect(saved).toEqual(['kept.org', 'browser.org']);
    rules.remove('kept.org');
    expect(saved).toEqual(['browser.org']);
  });
});
