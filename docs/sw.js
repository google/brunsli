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
  await waitForBrunsli();
  let jxlStream = await originalResponse.body;
  let jpegStream = jxlStream.pipeThrough(new BrunsliTransformer());
  let modifiedResponseHeaders = new Headers(originalResponse.headers);
  modifiedResponseHeaders.set('Content-Type', 'image/jpeg');
  return new Response(jpegStream, { headers: modifiedResponseHeaders });
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
class BrunsliTransformer {
  constructor() {
    let instance = brunsliModule._BrunsliDecoderInit();
    let increment = brunsliModule.HEAPU32[(instance + 4) >> 2];
    let input = brunsliModule.HEAPU32[(instance + 8) >> 2];
    let out = brunsliModule.HEAPU32[(instance + 16) >> 2];
    let readController = null;
    let writeController = null;

    let cleanup = () => {
      brunsliModule._BrunsliDecoderCleanup(instance);
      instance = 0;
    }

    this.readable = new ReadableStream({
      start(controller) {
        readController = controller;
      },
      cancel(reason) {
        cleanup();
      }
    });

    this.writable = new WritableStream({
      start(controller) {
        writeController = controller;
      },
      write(chunk, controller) {
        if (!instance) return;
        let pos = 0;
        while (true) {
          let end = Math.min(pos + increment, chunk.length);
          let delta = end - pos;
          if (delta > 0) brunsliModule.HEAPU8.set(chunk.subarray(pos, end), input);
          brunsliModule.HEAPU32[(instance + 12) >> 2] = delta;
          pos += delta;
          let result = brunsliModule._BrunsliDecoderProcess(instance);
          let outLen = brunsliModule.HEAPU32[(instance + 20) >> 2];
          if (outLen > 0) {
            let outChunk = new Uint8Array(outLen);
            outChunk.set(brunsliModule.HEAPU8.subarray(out, out + outLen));
            readController.enqueue(outChunk);
          }
          if (result == 0) {
            if ((delta == 0) && (outLen == 0)) break;
            continue;
          } else if (result == 1) {
            cleanup();
            readController.close();
            break;
          } else {
            cleanup();
            readController.error("corrupted Brunsli file");
            break;
          }
        }
      },
      close(controller) {
        cleanup();
        // TODO: what if input is truncated?
        readController.close();
      },
      abort(reason) {
        cleanup();
        readController.error(reason);
      }
    });
  }
}

importScripts('brunslidec-wasm.js');
brunsliModule = BrunsliDecModule();
brunsliModule.warm = false;
brunsliModule.ready = false;
