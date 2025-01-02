// Generated by running pwa.js, only modify this directly in testing as it is automatically (re)generated in the build process.
// The viewer service worker:
// - caches the static resources that the app needs to function
// - intercepts server requests and responds with cached responses instead of going to the network
// - deletes old caches on activation

const CACHE_NAME = `molstar-viewer-4.10.0`;

// The static resources that the app needs to function.
const APP_STATIC_RESOURCES = [
    "./",
    "src/apps/viewer/favicon.ico",
    "src/icons/circle.ico",
    "src/icons/circle.svg",
    "src/icons/tire.svg",
    "src/icons/wheel.svg",
    "src/apps/viewer/index-pwa.html"
];
    
// On install, cache the static resources.
self.addEventListener("install", (event) => {
    event.waitUntil(
        (async () => {
            const cache = await caches.open("molstar-viewer-4.10.0");
            await cache.addAll(["/"]);
            await self.skipWaiting();
        })(),
    );
});

// On activate, delete old caches.
self.addEventListener("activate", (event) => {
    event.waitUntil(
        (async () => {
            const keys = await caches.keys();
            await Promise.all(
                keys.map((key) => {
                    if (key !== "molstar-viewer-4.10.0") {
                        return caches.delete(key);
                    }
                }),
            );
            await clients.claim();
        })(),
    );
});

// On fetch, intercept server requests and respond with cached responses instead of going to network.
self.addEventListener("fetch", (event) => {
    // As a single page app, direct app to always go to cached home page.
    if (event.request.mode === "navigate") {
        event.respondWith(caches.match("/"));
        return;
    }

    // For all other requests, go to the cache first, and then the network.
    event.respondWith(
        (async () => {
            const cache = await caches.open("molstar-viewer-4.10.0");
            const cachedResponse = await cache.match(event.request.url);
            if (cachedResponse) {
                // If available, return the cached response.
                return cachedResponse;
            }
            return new Response(null, { status: 404 });
        })(),
    );
});