/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from '../../mol-task';
import { Asset, AssetManager } from '../assets';
import { ajaxGet } from '../data-source';
import { GraphQLClient } from '../graphql-client';

describe('graphql transport', () => {
    it('adds JSON headers to GraphQL requests', async () => {
        const assetManager = new AssetManager();
        const dispose = jest.fn();

        const resolveSpy = jest.spyOn(assetManager, 'resolve').mockImplementation((asset, type) => {
            expect(type).toBe('json');
            expect(Asset.isUrl(asset)).toBe(true);
            if (!Asset.isUrl(asset)) throw new Error('expected URL asset');

            expect(asset.url).toBe('https://example.org/graphql');
            expect(asset.headers).toEqual({
                'Content-Type': 'application/json; charset=utf-8',
                'Accept': 'application/json'
            });
            expect(asset.body).toContain('"query"');
            expect(asset.body).toContain('"variables"');

            return {
                id: 0,
                name: 'mock',
                run: async () => ({ data: { data: { ok: true } }, dispose }),
                runAsChild: async () => ({ data: { data: { ok: true } }, dispose }),
                runInContext: async () => ({ data: { data: { ok: true } }, dispose }),
            } as any;
        });

        const client = new GraphQLClient('https://example.org/graphql', assetManager);
        const result = await client.request(RuntimeContext.Synchronous, 'query Test { test }', { id: '1' });

        expect(resolveSpy).toHaveBeenCalledTimes(1);
        expect(result.data).toEqual({ ok: true });
        result.dispose();
        expect(dispose).toHaveBeenCalledTimes(1);
    });

    it('preserves POST body and headers in Node.js HTTP requests', async () => {
        const fetchSpy = jest.spyOn(globalThis, 'fetch').mockResolvedValue({
            status: 200,
            bytes: async () => new TextEncoder().encode(JSON.stringify({ ok: true }))
        } as any);

        const result = await ajaxGet({
            url: 'https://example.org/graphql',
            type: 'json',
            body: '{"query":"{ test }"}',
            headers: {
                'Content-Type': 'application/json',
                'Accept': 'application/json'
            }
        }).run();

        expect(fetchSpy).toHaveBeenCalledWith('https://example.org/graphql', {
            signal: expect.any(AbortSignal),
            method: 'POST',
            body: '{"query":"{ test }"}',
            headers: {
                'Content-Type': 'application/json',
                'Accept': 'application/json'
            }
        });
        expect(result).toEqual({ ok: true });

        fetchSpy.mockRestore();
    });
});
