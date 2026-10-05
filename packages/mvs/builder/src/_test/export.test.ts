import { unzipSync } from 'fflate';
import { createMVSX } from '../export.js';
import { MVSData } from '../mvs-data.js';

function dataWithAsset(url: string) {
  const builder = MVSData.createBuilder();
  builder.download({ url });
  return builder.getState();
}

afterEach(() => jest.restoreAllMocks());

test('toMVSX defaults to platform fetch when options are omitted', async () => {
  const fetchAsset = jest.spyOn(globalThis, 'fetch').mockResolvedValue(new Response('asset content'));
  const files = unzipSync(await MVSData.toMVSX(dataWithAsset('https://example.org/asset.txt')));
  expect(fetchAsset).toHaveBeenCalledWith('https://example.org/asset.txt');
  expect(new TextDecoder().decode(files['./assets/0-asset.txt'])).toBe('asset content');
});

test('toMVSX uses custom fetch with resolved URIs and reuses cached assets', async () => {
  const platformFetch = jest.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('unexpected platform fetch'));
  const fetchAsset = jest.fn(async () => new Response('local asset'));
  const cache = {};
  const data = dataWithAsset('asset.txt');
  const options = { baseUri: 'file:///data/scene.mvsj', fetch: fetchAsset, cache };
  const files = unzipSync(await MVSData.toMVSX(data, options));
  fetchAsset.mockRejectedValue(new Error('source is now unavailable'));
  expect(unzipSync(await MVSData.toMVSX(data, options))).toEqual(files);
  expect(fetchAsset).toHaveBeenCalledTimes(1);
  expect(fetchAsset).toHaveBeenCalledWith('file:///data/asset.txt');
  expect(platformFetch).not.toHaveBeenCalled();
  expect(new TextDecoder().decode(files['./assets/0-asset.txt'])).toBe('local asset');
  const exported = JSON.parse(new TextDecoder().decode(files['index.mvsj']));
  expect(exported.root.children[0].params.url).toBe('./assets/0-asset.txt');
  expect(data.root.children![0].params.url).toBe('asset.txt');
});

test('toMVSX treats empty cached content as a cache hit', async () => {
  const fetchAsset = jest.fn(async () => {
    throw new Error('unexpected fetch');
  });
  const url = 'https://example.org/empty.txt';
  const files = unzipSync(await MVSData.toMVSX(dataWithAsset(url), { cache: { [url]: '' }, fetch: fetchAsset }));
  expect(files['./assets/0-empty.txt']).toHaveLength(0);
  expect(fetchAsset).not.toHaveBeenCalled();
});

test('toMVSX bypasses fetching for explicit assets and skipped external URIs', async () => {
  const fetchAsset = jest.fn(async () => {
    throw new Error('unexpected fetch');
  });
  const url = 'https://example.org/asset.txt';
  const data = dataWithAsset(url);
  const explicit = unzipSync(await MVSData.toMVSX(data, { assets: { [url]: 'explicit' }, fetch: fetchAsset }));
  expect(new TextDecoder().decode(explicit['./assets/0-asset.txt'])).toBe('explicit');
  const skipped = unzipSync(await MVSData.toMVSX(data, { skipExternal: true, fetch: fetchAsset }));
  expect(Object.keys(skipped)).toEqual(['index.mvsj']);
  expect(fetchAsset).not.toHaveBeenCalled();
});

test('toMVSX reports failed responses without caching them', async () => {
  const fetchAsset = jest.fn(async () => new Response(null, { status: 404, statusText: 'Not Found' }));
  const cache = {};
  const url = 'https://example.org/missing.txt';
  await expect(MVSData.toMVSX(dataWithAsset(url), { fetch: fetchAsset, cache })).rejects.toThrow(
    `Failed to fetch MVS asset ${url}: 404 Not Found`,
  );
  expect(cache).toEqual({});
});

test('MVSX output is deterministic across different clock times with explicit mtime', async () => {
  jest.useFakeTimers();
  try {
    const data = MVSData.createBuilder().getState();
    const assets = [{ name: 'asset.txt', content: 'archive asset' }];
    // ZIP stores local date fields; use matching fields for reproducibility across time zones.
    const options = { zip: { mtime: new Date(1980, 0, 1), level: 0 as const } };
    jest.setSystemTime(new Date('2025-01-01T12:00:00Z'));
    const first = await createMVSX(data, assets, options);
    jest.setSystemTime(new Date('2026-10-04T18:00:00Z'));
    const second = await createMVSX(data, assets, options);
    expect(second).toEqual(first);
    const files = unzipSync(first);
    expect(new TextDecoder().decode(files['asset.txt'])).toBe('archive asset');
    expect(JSON.parse(new TextDecoder().decode(files['index.mvsj']))).toEqual(data);
  } finally {
    jest.useRealTimers();
  }
});

test('MVSX output uses the current time when ZIP options are omitted', async () => {
  jest.useFakeTimers();
  try {
    const data = MVSData.createBuilder().getState();
    jest.setSystemTime(new Date('2025-01-01T12:00:00Z'));
    const first = await createMVSX(data, []);
    jest.setSystemTime(new Date('2026-10-04T18:00:00Z'));
    const second = await createMVSX(data, []);
    expect(second).not.toEqual(first);
    expect(unzipSync(second)).toEqual(unzipSync(first));
  } finally {
    jest.useRealTimers();
  }
});
