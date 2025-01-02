#!/usr/bin/env node
/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
const fs = require('fs');
const path = require('path');

// Function to delete manifest files
function deleteManifestFiles(dir) {
    fs.readdir(dir, (err, files) => {
        if (err) {
            console.error('Error reading directory:', err);
            return;
        }

        files.forEach(file => {
            const filePath = path.join(dir, file);
            if (file.startsWith('manifest-') && file.endsWith('.webmanifest') ||
                file.startsWith('sw-') && file.endsWith('.js')) {
                // Delete manifest files
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

deleteManifestFiles(path.join(__dirname, '..'));