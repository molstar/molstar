
The files in this directory are used for deploying the Mol* Viewer as a PWA (Progressive Web App) at https://molstar.org/viewer/. They may serve as an example for creating your own PWA but wont work as-is. See `/script/deploy.js` for where these files are copied and how they are transformed during deployment.


## PWA features

- The Service Worker will cache static resources so the Viewer can be used without internet access. This works without installing, i.e., also in Firefox.
- Once installed, file types listed in the Manifest can be opened from, e.g., the Windows File Explorer.


## Notes for development

In Chrome you can see a list of installed PWAs at chrome://apps/. A right-click opens a menu with an option uninstall.

The Chrome Dev Tools have a section 'Application' to inspect and manage PWA aspects like the Manifest and Service Workers.
