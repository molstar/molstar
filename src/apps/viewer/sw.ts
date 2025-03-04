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

// Hardcoded version number as all ways of injecting the version number failed :(
//const version = '1.0.0'; // Replace with your actual version number

// Declare the VERSION constant injected by Webpack
declare const VERSION: string;
//// Use the injected version number
const version = VERSION;
//const version = process.env.VERSION;
console.log(`Service Worker version: ${version}`);
const CACHE_NAME = `molstar-viewer-${version}`;

//const version = process.env.VERSION;
//const CACHE_NAME = `molstar-viewer-${version}`;

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
            await swSelf.skipWaiting(); // Ensures the new service worker takes control immediately.
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
            await swSelf.clients.claim(); // Ensures the new service worker takes control immediately.
        })(),
    );
});

// On fetch, respond with cached resources.
swSelf.addEventListener("fetch", (event: FetchEvent) => {
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