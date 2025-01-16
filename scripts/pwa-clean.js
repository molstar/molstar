#!/usr/bin/env node
/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
const fs = require('fs');
const path = require('path');

// General function to delete files based on prefix and suffix
function deleteFiles(dir, prefix, suffix) {
    fs.readdir(dir, (err, files) => {
        if (err) {
            console.error('Error reading directory:', err);
            return;
        }
        files.forEach(file => {
            const filePath = path.join(dir, file);
            if (file.startsWith(prefix) && file.endsWith(suffix)) {
                fs.unlink(filePath, err => {
                    if (err) {
                        console.error('Error deleting file:', err);
                        return;
                    }
                    console.log(`Deleted: ${filePath}`);
                });
            }
        });
    });
}

// Function to delete manifest files
function deleteManifestFiles(dir) {
    deleteFiles(dir, 'manifest-', '.webmanifest');
}

// Function to delete service worker files
function deleteServiceWorkerFiles(dir) {
    deleteFiles(dir, 'sw-', '.js');
}

deleteManifestFiles(path.join(__dirname, '..'));
deleteServiceWorkerFiles(__dirname);