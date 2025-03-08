/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 *
 * @description Service worker for the viewer app.
 * The service worker:
 * - caches the static resources that the app needs to function
 * - intercepts server requests and responds with cached responses instead of going to the network
 * - deletes old caches on activation
 * - responds with cached resources on fetch
 * - responds with a network error if fetching fails
 * - responds with a cache error if the cache match fails
 *
 * To not complicate the build process, this file is not transpiled but written in JavaScript.
 * It is coppied to the build directory as is and no imports must be used.
 */

/** version from package.json, to be filled in at build time */
const VERSION = '';

const CACHE_NAME = `molstar-viewer-${VERSION}`;

// The static resources that the app needs to function.
const APP_STATIC_RESOURCES = [
    'favicon.ico',
    'index.html',
    'molstar.css',
    'molstar.js',
    'manifest.webmanifest',
    'logo-144.png'
];

// On install, cache the static resources.
self.addEventListener('install', (event) => {
    console.log(`Service Worker version ${VERSION} installed.`);
    event.waitUntil(
        (async () => {
            const cache = await caches.open(CACHE_NAME);
            await cache.addAll(APP_STATIC_RESOURCES);
            await self.skipWaiting(); // Ensures the new service worker takes control immediately.
        })(),
    );
});

// On activate, delete old caches.
self.addEventListener('activate', (event) => {
    console.log(`Service Worker version ${VERSION} activated.`);
    event.waitUntil(
        (async () => {
            const keys = await caches.keys();
            await Promise.all(
                keys.map((key) => {
                    if (key !== CACHE_NAME) {
                        return caches.delete(key);
                    }
                }),
            );
            await self.clients.claim(); // Ensures the new service worker takes control immediately.
        })(),
    );
});

// On fetch, respond with cached resources.
self.addEventListener('fetch', (event) => {
    event.respondWith(
        (async () => {
            try {
                // Try to match the request with the cache
                const cachedResponse = await caches.match(event.request);
                if (cachedResponse) {
                    return cachedResponse;
                }
                // Fetch from network if not found in cache
                const networkResponse = await fetch(event.request);
                // Check if the network response is valid
                if (!networkResponse || networkResponse.status !== 200 || networkResponse.type !== 'basic') {
                    return networkResponse;
                }
                // Clone the network response
                const responseToCache = networkResponse.clone();
                // Open the cache and put the network response in it
                const cache = await caches.open(CACHE_NAME);
                await cache.put(event.request, responseToCache);
                return networkResponse;
            } catch (error) {
                console.error('Fetching failed:', error);
                return new Response('Network error occurred', {
                    status: 408,
                    statusText: 'Network error occurred'
                });
            }
        })()
    );
});
