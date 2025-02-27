/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 * @description Service worker for the viewer app.
 * The service worker:
 * - caches the static resources that the app needs to function
 * - intercepts server requests and responds with cached responses instead of going to the network
 * - deletes old caches on activation
 * - responds with cached resources on fetch
 * - responds with a network error if fetching fails
 * - responds with a cache error if the cache match fails
 */
/// <reference lib="webworker" />

// Added to make the file into a module in order to avoid the following error:
// 'Cannot use import statement outside a module'
export { };

// Top-level async function to import package.json and set CACHE_NAME
(async () => {

    // Use dynamic import to load the JSON file
    const packageJson = await import('../../../package.json');
    const version = packageJson.version;

    // Previous attempts to load the JSON file and import the version.
    //// Use require to load the JSON file
    //const packageJson = require('../../../package.json');
    //const version = packageJson.version;
    //import { version } from '../../../package.json';

    const CACHE_NAME = `molstar-viewer-` + version;

    // The static resources that the app needs to function.
    const APP_STATIC_RESOURCES = [
        "favicon.ico",
        "circle.ico",
        "circle.svg",
        "wheel.svg",
        "tire.svg",
        "index.html",
        "molstar.css",
        "molstar.js",
        "manifest.webmanifest"
    ];

    // Cast self to ServiceWorkerGlobalScope
    const swSelf = self as unknown as ServiceWorkerGlobalScope;

    // On install, cache the static resources.
    swSelf.addEventListener("install", (event: ExtendableEvent) => {
        console.log(`Service Worker version ${version} installed.`);
        event.waitUntil(
            (async () => {
                const cache = await caches.open(CACHE_NAME);
                await cache.addAll(APP_STATIC_RESOURCES);
                await swSelf.skipWaiting();
            })(),
        );
    });

    // On activate, delete old caches.
    swSelf.addEventListener("activate", (event: ExtendableEvent) => {
        console.log(`Service Worker version ${version} activated.`);
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
                await swSelf.clients.claim();
            })(),
        );
    });


    // On fetch, respond with cached resources.
    swSelf.addEventListener("fetch", (event: FetchEvent) => {
        event.respondWith(
            caches.match(event.request).then((response) => {
                // Return the cached response if found, otherwise fetch from network
                return response || fetch(event.request).then((networkResponse) => {
                    // Check if the network response is valid
                    if (!networkResponse || networkResponse.status !== 200 || networkResponse.type !== 'basic') {
                        return networkResponse;
                    }

                    // Clone the network response
                    const responseToCache = networkResponse.clone();

                    // Open the cache and put the network response in it
                    caches.open(CACHE_NAME).then((cache) => {
                        cache.put(event.request, responseToCache);
                    });

                    return networkResponse;
                }).catch((error) => {
                    console.error('Fetching failed:', error);
                    return new Response('Network error occurred', {
                        status: 408,
                        statusText: 'Network error occurred'
                    });
                });
            }).catch((error) => {
                console.error('Cache match failed:', error);
                return new Response('Cache error occurred', {
                    status: 408,
                    statusText: 'Cache error occurred'
                });
            })
        );
    });
})();