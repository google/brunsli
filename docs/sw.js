let brunsliModule = null;

self.addEventListener('install', event => {self.skipWaiting();});
self.addEventListener('activate', event => {event.waitUntil(clients.claim());});
self.addEventListener('fetch', onFetch);

function shouldUpgradeRequest(request) {
  if (request.method != 'GET') return false;
  if (request.referrer.includes('no-jxl')) return false;
  if (request.url.includes('no-jxl')) return false;
  if (request.headers.get('Accept').indexOf('image/') == -1) return false;
  if (request.destination != "image") return false;
  return true;
}
async function waitForBrunsli() {
  if (brunsliModule && brunsliModule.ready) return;
  while (true) {
    if (brunsliModule && brunsliModule.then && !brunsliModule.warm) {
      brunsliModule.warm = true;
      brunsliModule.then(wasmModule => brunsliModule.ready = true);
    }
    if (brunsliModule && brunsliModule.ready) break;
    await new Promise((resolve, reject) => setTimeout(resolve, 1));
  }
}
async function wrapRequest(request) {
  let modifiedRequestHeaders = new Headers(request.headers);
  modifiedRequestHeaders.append('Accept', 'image/x-j');
  let modifiedRequest = new Request(request, {headers: modifiedRequestHeaders});
  let originalResponse = await fetch(modifiedRequest);
  let contentType = originalResponse.headers.get('Content-Type');
  if (contentType != 'image/x-j') return originalResponse;
  let data = await originalResponse.arrayBuffer();
  await waitForBrunsli();
  let jpeg = decodeBrunsli(data);
  if (!jpeg) return new Response(data, { headers: originalResponse.headers });
  let modifiedResponseHeaders = new Headers(originalResponse.headers);
  modifiedResponseHeaders.set('Content-Type', 'image/jpeg');
  return new Response(jpeg, { headers: modifiedResponseHeaders });
}
async function onFetch(event) {
  if (!shouldUpgradeRequest(event.request)) return;
  event.respondWith(wrapRequest(event.request));
}
function decodeBrunsli(arrayBuffer) {
  let size = arrayBuffer.byteLength;
  let bytesView = new Uint8Array(arrayBuffer);
  let data = brunsliModule._malloc(size);
  try {
    brunsliModule.HEAPU8.set(bytesView, data);
    let jpeg = brunsliModule._BrunsliToJpeg(data, size);
    if (!jpeg) return null;
    let jpegData = brunsliModule._GetJpegData(jpeg);
    let jpegSize = brunsliModule._GetJpegLength(jpeg);
    let jpegDataView = new Uint8Array(
        brunsliModule.HEAPU8.subarray(jpegData, jpegData + jpegSize));
    let output = new Uint8Array(jpegSize);
    output.set(jpegDataView);
    brunsliModule._FreeJpeg(jpeg);
    return output;
  } finally {
    brunsliModule._free(data);
  }
  return data;
}

importScripts('brunslidec-wasm.js');
brunsliModule = BrunsliDecModule();
brunsliModule.warm = false;
brunsliModule.ready = false;
