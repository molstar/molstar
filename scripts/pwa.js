#!/usr/bin/env node
/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
const fs = require('fs');
const path = require('path');
const { version } = require('../package.json');

// Generate the service worker file.
function generateServiceWorkerFile(name) {
    const cacheName = `molstar-${name}-${version}`;
    const entryPoint = `index-pwa.html`;
    const filePath = path.join(__dirname, '..', 'build', `${name}`, `sw.js`);
    // Set faviconpath.
    // Initialise content.
    const content = `// Generated by running pwa.js, only modify this directly in testing as it is automatically (re)generated in the build process.
// The service worker:
// - caches the static resources that the app needs to function
// - intercepts server requests and responds with cached responses instead of going to the network
// - deletes old caches on activation

const CACHE_NAME = \`${cacheName}\`;

// The static resources that the app needs to function.
const APP_STATIC_RESOURCES = [
    "favicon.ico",
    "circle.ico",
    "circle.svg",
    "wheel.svg",
    "tire.svg",
    "${entryPoint}",
    "molstar.css",
    "molstar.js",
    "manifest.webmanifest"
];
    
// On install, cache the static resources.
self.addEventListener("install", (event) => {
    event.waitUntil(
        (async () => {
            const cache = await caches.open(CACHE_NAME);
            await cache.addAll(APP_STATIC_RESOURCES);
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
                    if (key !== CACHE_NAME) {
                        return caches.delete(key);
                    }
                }),
            );
            await clients.claim();
        })(),
    );
});

// On fetch, respond with cached resources.
self.addEventListener("fetch", (event) => {
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
});`;
    // Write content to file.
    fs.writeFile(filePath, content, 'utf8', (err) => {
        if (err) {
            console.error('Error writing the service worker file:', err);
            return;
        }
        console.log(`${filePath} has been created successfully.`);
    });
}

// Generate the manifest file
function generateManifestFile(name) {
    const filePath = path.join(__dirname, '..', 'build', `${name}`, 'manifest.webmanifest');
    const content=`{
        "id": "/molstar/build/viewer/index-pwa.html",
        "name": "Mol* ${name}",
        "short_name": "Mol*",
        "description": "Mol* ${name}",
        "start_url": "./index-pwa.html",
        "theme_color": "#eeffee",
        "background_color": "#eeffee",
        "display": "standalone",
        "icons": [
            {
                "src": "favicon.ico",
                "sizes": "48x48"
            },
            {
                "src": "circle.ico",
                "sizes": "72x72"
            },
            {
                "src": "circle.svg",
                "sizes": "72x72 96x96",
                "purpose": "maskable"
            },
            {
                "src": "wheel.svg",
                "sizes": "144x144"
            },
            {
                "src": "tire.svg",
                "sizes": "192x192"
            }
        ]
    }`;
    // Write content to file.
    fs.writeFile(filePath, content, 'utf8', (err) => {
        if (err) {
            console.error('Error writing the service worker file:', err);
            return;
        }
        console.log(`${filePath} has been created successfully.`);
    });
}

// Read the existing index.html file and create an index-pwa.html file injected with what a PWA needs.
function processHtmlFile(name) {
    
    // Define the source and destination file paths
    let inFilePath = path.join(__dirname, '..', 'src', 'apps', `${name}`, 'index.html');
    let outFilePath = path.join(__dirname, '..', 'build', `${name}`, 'index-pwa.html');

    // Read the content of the source file
    fs.readFile(inFilePath, 'utf8', (err, data) => {
        if (err) {
            console.error('Error reading the source file:', err);
            return;
        }

        const manifest = `    <link rel="manifest" href="./manifest.webmanifest" type="application/manifest+json">`;

        // Insert the manifest before the </head> tag
        let modifiedContent = data.replace('</head>', `${manifest}\n    </head>`);

        // Define the service worker registration script
        const sWScript = `    <!-- Register the app's service worker. -->
        <script type="module">
            if ('serviceWorker' in navigator) {
                window.addEventListener('load', function () {
                    navigator.serviceWorker.register(new URL('./sw.js', import.meta.url))
                    .then(function (registration) {
                        // Registration was successful
                        console.log('ServiceWorker registration successful with scope: ', registration.scope);
                    }, function (err) {
                        // registration failed :(
                        console.log('ServiceWorker registration failed: ', err);
                    });
                });
            }
        </script>
        `;

        // Insert the service worker script before the closing </body> tag
        modifiedContent = modifiedContent.replace('</body>', `${sWScript}</body>`);

        // Write the modified content to the destination file
        fs.writeFile(outFilePath, modifiedContent, 'utf8', (err) => {
            if (err) {
                console.error('Error writing to the destination file:', err);
                return;
            }
            console.log(`${outFilePath} has been created successfully.`);
        });
    });
}

function createPWA() {
    let name = 'viewer';
    // (re)generate the service worker file
    generateServiceWorkerFile(name);
    // (re)generate the manifest files
    generateManifestFile(name);
    // (re)generate the html files
    processHtmlFile(name);
}

// Create the Progressive Web App
createPWA();