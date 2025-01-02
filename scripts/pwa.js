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
function generateServiceWorkerFile(name, dirname) {
    const cacheName = `molstar-${name}-${version}`;
    const filePath = path.join(__dirname, '..', `sw-${name}.js`);
    // Set faviconpath.
    let faviconPath;
    if (dirname == 'example') {
        faviconPath = "src/apps/viewer/favicon.ico";
    } else {
        faviconPath = `src/${dirname}/viewer/favicon.ico`;
    }
    // Initialise content.
    const content = `// Generated by running pwa.js, only modify this directly in testing as it is automatically (re)generated in the build process.
// The ${name} service worker:
// - caches the static resources that the app needs to function
// - intercepts server requests and responds with cached responses instead of going to the network
// - deletes old caches on activation

const CACHE_NAME = \`${cacheName}\`;

// The static resources that the app needs to function.
const APP_STATIC_RESOURCES = [
    "./",
    "${faviconPath}",
    "src/icons/circle.ico",
    "src/icons/circle.svg",
    "src/icons/tire.svg",
    "src/icons/wheel.svg",
    "src/${dirname}/${name}/index-pwa.html"
];
    
// On install, cache the static resources.
self.addEventListener("install", (event) => {
    event.waitUntil(
        (async () => {
            const cache = await caches.open("${cacheName}");
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
                    if (key !== "${cacheName}") {
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
            const cache = await caches.open("${cacheName}");
            const cachedResponse = await cache.match(event.request.url);
            if (cachedResponse) {
                // If available, return the cached response.
                return cachedResponse;
            }
            return new Response(null, { status: 404 });
        })(),
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
function generateManifestFile(name, dirname) {
    const filePath = path.join(__dirname, '..', `manifest-${name}.webmanifest`);
    const content=`{
    "name": "Mol* ${name}",
    "short_name": "Mol*",
    "description": "Mol* ${name}",
    "start_url": "molstar/build/${name}/index-pwa.html",
    "theme_color": "#eeffee",
    "background_color": "#eeffee",
    "display": "standalone",
    "icons": [
        {
            "src": "src/${dirname}/${name}/favicon.ico",
            "sizes": "48x48"
        },
        {
            "src": "src/icons/circle.svg",
            "sizes": "72x72 96x96",
            "purpose": "maskable"
        },
        {
            "src": "src/icons/tire.svg",
            "sizes": "128x128 256x256"
        },
        {
            "src": "src/icons/wheel.svg",
            "sizes": "512x512"
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
function processHtmlFile(name, dirname) {
    
    // Define the source and destination file paths
    let inFilePath = path.join(__dirname, '..', `src/${dirname}/`, `${name}/index.html`);
    let outFilePath = path.join(__dirname, '..', `src/${dirname}/`, `${name}/index-pwa.html`);

    // Read the content of the source file
    fs.readFile(inFilePath, 'utf8', (err, data) => {
        if (err) {
            console.error('Error reading the source file:', err);
            return;
        }

        const manifest = `    <link rel="manifest" href="/molstar/manifest-${name}.webmanifest" type="application/manifest+json">`;

        // Insert the manifest before the </head> tag
        let modifiedContent = data.replace('</head>', `${manifest}\n    </head>`);

        // Define the service worker registration script
        const sWScript = `    <!-- Register the app's service worker. -->
            <script type="module">
                if ('serviceWorker' in navigator) {
                    window.addEventListener('load', function () {
                        //navigator.serviceWorker.register('/sw.js')
                        navigator.serviceWorker.register(new URL('/sw-${name}.js', import.meta.url))
                        //navigator.serviceWorker.register(new URL('/sw.js', import.meta.url), { scope: '/' })
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


function createPWA(name, dir) {
    // (re)generate the service worker file
    generateServiceWorkerFile(name, dir);
    // (re)generate the manifest files
    generateManifestFile(name);
    // (re)generate the html files
    processHtmlFile(name, dir);
}

// Apps
createPWA('viewer', 'apps');
createPWA('docking-viewer', 'apps');
createPWA('mesoscale-explorer', 'apps');

// Examples
createPWA('alpha-orbitals', 'examples');
createPWA('alphafolddb-pae', 'examples');
createPWA('basic-wrapper', 'examples');
//createPWA('domain-annotation-server', 'examples');
//createPWA('gbl-export', 'examples');
//createPWA('image-renderer', 'examples');
createPWA('lighting', 'examples');
createPWA('proteopedia-wrapper', 'examples');