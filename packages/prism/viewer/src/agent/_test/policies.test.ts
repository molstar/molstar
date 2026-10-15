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
import { addressesInRequest, allowedAddress, createAgentPolicy, createDomainRules, mayAsk } from '../policies.js';

// What goes without asking (agent.md §4).

const op = (o: Partial<OperationDefinition<any, any>>) =>
  ({ id: 'x.y', title: 'X', params: {}, run: async () => ({}) as never, ...o }) as OperationDefinition<any, any>;

describe('the agent policy', () => {
  it('lets changes go while the turn can be undone, and asks for removals and uploads', () => {
    const p = createAgentPolicy({ checkpointing: observableState(true).asReadonly });
    expect(p.approves(op({}), {})).toBe(true);
    expect(p.approves(op({ approval: 'destructive' }), {})).toBe(false);
    expect(p.approves(op({ approval: 'upload' }), {})).toBe(false);
  });

  it('asks once for the turn when undo is off', () => {
    const checkpointing = observableState(false).asReadonly;
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
    const p = createAgentPolicy({ checkpointing: observableState(true).asReadonly, rules });
    const net = op({ approval: 'network', addresses: ({ url }) => [url] });
    expect(p.approves(net, { url: 'https://files.rcsb.org/download/1cbs.cif' })).toBe(true);
    expect(p.approves(net, { url: 'https://example.org/x.cif' })).toBe(false);
    expect(p.approves(net, { url: 'http://files.rcsb.org/x.cif' })).toBe(false);
    rules.allow('example.org', 'session');
    expect(p.approves(net, { url: 'https://data.example.org/x.cif' })).toBe(true);
    expect(allowedAddress('https://evil-rcsb.org/x', rules)).toBe(false);
  });

  it('asks for an address it cannot read whole, and for a network operation that does not say what it reaches', () => {
    const p = createAgentPolicy({ checkpointing: observableState(true).asReadonly });
    const net = op({ approval: 'network', addresses: ({ url }) => [url] });
    // The forms that slipped past a check of top-level `scheme://` strings (G6).
    for (const url of [' https://evil.example/x.cif', '//evil.example/x.cif', 'evil.example/x.cif', 'rcsb.org/x.cif']) {
      expect(p.approves(net, { url }), url).toBe(false);
    }
    expect(p.approves(op({ approval: 'network' }), { url: 'https://files.rcsb.org/x.cif' })).toBe(false);
    // The known services only (the operation says so).
    expect(p.approves(op({ approval: 'network', addresses: () => [] }), { id: '1cbs' })).toBe(true);
    // Several addresses: each one counts.
    const two = op({ approval: 'network', addresses: ({ a, b }) => [a, b] });
    expect(p.approves(two, { a: 'https://files.rcsb.org/x', b: 'https://evil.example/y' })).toBe(false);
  });

  it('lets the hosts the request names go for that turn only', () => {
    const p = createAgentPolicy({ checkpointing: observableState(true).asReadonly });
    const net = op({ approval: 'network', addresses: ({ url }) => [url] });
    p.beginTurn(addressesInRequest('load https://lab.example.edu/files/model.cif please, not http://plain.example'));
    expect(p.approves(net, { url: 'https://lab.example.edu/files/other.cif' })).toBe(true);
    expect(p.approves(net, { url: 'https://plain.example/x' })).toBe(false);
    p.beginTurn();
    expect(p.approves(net, { url: 'https://lab.example.edu/files/other.cif' })).toBe(false);
  });

  it('asks for a change no undo restores, and every call that may ask is known up front', () => {
    const p = createAgentPolicy({ checkpointing: observableState(true).asReadonly });
    expect(p.approves(op({ undoable: false }), {})).toBe(false);
    expect(mayAsk(op({ read: true }))).toBe(false);
    expect(mayAsk(op({ history: false }))).toBe(false);
    expect(mayAsk(op({ read: true, approval: 'network' }))).toBe(true);
    expect(mayAsk(op({ undoable: false }))).toBe(true);
    // A file made for the person (an export: no history step, nothing to undo) asks too.
    expect(mayAsk(op({ history: false, undoable: false }))).toBe(true);
    expect(p.approves(op({ history: false, undoable: false }), {})).toBe(false);
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
