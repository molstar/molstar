/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** version from package.json, to be filled in during deployment */
const VERSION = '__MOLSTAR_VERSION__';

const CACHE_NAME = `molstar-viewer-${VERSION}`;

// The static resources that the app needs to function.
const APP_STATIC_RESOURCES = [
    'favicon.ico',
    'index.html',
    'molstar.css',
    'molstar.js',
    'manifest.webmanifest',
    'logo-144.png',
    'pwa.js'
];

async function cacheStaticResources() {
    const cache = await caches.open(CACHE_NAME);
    await cache.addAll(APP_STATIC_RESOURCES);
    await self.skipWaiting(); // Ensures the new service worker takes control immediately.
}

async function deleteOldCaches() {
    const keys = await caches.keys();
    await Promise.all(
        keys.map((key) => {
            if (key !== CACHE_NAME) {
                return caches.delete(key);
            }
        }),
    );
    await self.clients.claim(); // Ensures the new service worker takes control immediately.
}

async function respondWithCacheFirst(request) {
    // Try to match the request with the cache
    const cachedResponse = await caches.match(request);
    return cachedResponse || fetch(request);
}

self.addEventListener('install', (event) => {
    // console.log(`Service Worker version ${VERSION} installed.`);
    event.waitUntil(cacheStaticResources());
});

self.addEventListener('activate', (event) => {
    // console.log(`Service Worker version ${VERSION} activated.`);
    event.waitUntil(deleteOldCaches());
});

self.addEventListener('fetch', (event) => {
    event.respondWith(respondWithCacheFirst(event.request));
});
