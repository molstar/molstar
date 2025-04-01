/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

window.addEventListener('molstarViewerCreated', e => {
    const viewer = e.detail.viewer;

    // Handle incoming files
    if ('launchQueue' in window) {
        launchQueue.setConsumer((launchParams) => {
            if (!launchParams.files.length) return;

            const files = [];
            for (const fileHandle of launchParams.files) {
                files.push(fileHandle.getFile());
            }

            Promise.all(files).then((files) => {
                viewer.loadFiles(files);
            });
        });
    }
});

// Register Progressive Web App service worker.
if ('serviceWorker' in navigator) {
    window.addEventListener('load', function () {
        navigator.serviceWorker.register('./sw.js')
            .then(function (registration) {
                // Registration was successful
                if (molstar.isDebugMode) {
                    console.log('ServiceWorker registration successful with scope: ', registration.scope);
                }
            }, function (err) {
                // registration failed :(
                if (molstar.isDebugMode) {
                    console.error('ServiceWorker registration failed: ', err);
                }
            });
    });
}
