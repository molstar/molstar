## Downloads by `fetch`

- `ajaxGet` (and so `PluginContext.fetch` and URL assets) downloads by `fetch` instead of `XMLHttpRequest`, one path
  for browsers and Node.js; `file://` URLs in Node.js are still read through `setFSModule`. Progress comes from the
  response stream against `Content-Length` (bytes only when the length is unknown or the body is compressed); an
  abort, seen between chunks as with XHR, cancels the request
- `AjaxGetParams.signal`: an `AbortSignal` that cancels the download also while it waits for the server or on a stalled
  body; the task then rejects as aborted
- Every download is read as bytes and decoded in JS, as opened files are: a `string` download past the 512 MB string
  limit comes back as a `ChunkedBigString` instead of failing. Text is always UTF-8 (XHR followed a `charset` in
  `Content-Type`); JSON and XML parse from it
- In a browser a network failure rejects with an `Error` carrying the same message, not a bare string (Node.js keeps
  its own error); a `json` body that does not parse gives `null`, now in Node.js too; a `zip` download is unzipped in
  the browser too, as its type says (it gave the raw bytes); malformed XML gives the parser's error document (XHR gave
  `null`)
- A `file:` URL in a browser (a page opened from disk, an Electron host) no longer loads: browsers' `fetch` has no
  `file:` scheme, where XHR read it when the browser allowed file access
