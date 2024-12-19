// The docking-viewer service worker:
// - caches the static resources that the app needs to function
// - intercepts server requests and responds with cached responses instead of going to the network
// - deletes old caches on activation
import { version } from './package.json';

const CACHE_NAME = `molstar-docking-viewer-${version}`;

// The static resources that the app needs to function.
const APP_STATIC_RESOURCES = [
    "./",
    "src/apps/viewer/favicon.ico",
    "src/icons/circle.ico",
    "src/icons/circle.svg",
    "src/icons/tire.svg",
    "src/icons/wheel.svg",
    "src/apps/docking-viewer/index.html"
];

// On install, cache the static resources
self.addEventListener("install", (event) => {
    event.waitUntil(
        (async () => {
            const cache = await caches.open(CACHE_NAME);
            cache.addAll(APP_STATIC_RESOURCES);
        })(),
    );
});

// delete old caches on activate
self.addEventListener("activate", (event) => {
    event.waitUntil(
        (async () => {
            const names = await caches.keys();
            await Promise.all(
                names.map((name) => {
                    if (name !== CACHE_NAME) {
                        return caches.delete(name);
                    }
                }),
            );
            await clients.claim();
        })(),
    );
});

// On fetch, intercept server requests
// and respond with cached responses instead of going to network
self.addEventListener("fetch", (event) => {
    // As a single page app, direct app to always go to cached home page.
    if (event.request.mode === "navigate") {
        event.respondWith(caches.match("/"));
        return;
    }

    // For all other requests, go to the cache first, and then the network.
    event.respondWith(
        (async () => {
            const cache = await caches.open(CACHE_NAME);
            const cachedResponse = await cache.match(event.request.url);
            if (cachedResponse) {
                // Return the cached response if it's available.
                return cachedResponse;
            }
            // If resource isn't in the cache, return a 404.
            return new Response(null, { status: 404 });
        })(),
    );
});