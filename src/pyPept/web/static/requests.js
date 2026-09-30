// Only use this for calculations that do not write library definitions.
// A busy single-worker server may reject overlapping previews and renders.
async function fetchCalculation(url, options = {}) {
  for (let attempt = 0; ; attempt++) {
    if (options.signal?.aborted) throw options.signal.reason || new Error('Request cancelled');
    const response = await fetch(url, options);
    const retryAfter = Number(response.headers?.get('Retry-After'));
    if (response.status !== 503 || attempt >= 2 || !Number.isFinite(retryAfter) ||
        retryAfter <= 0 || retryAfter > 2) return response;
    await response.body?.cancel();
    await new Promise((resolve, reject) => {
      const signal = options.signal;
      let timer;
      const aborted = () => {
        clearTimeout(timer);
        reject(signal.reason || new Error('Request cancelled'));
      };
      if (signal?.aborted) { aborted(); return; }
      signal?.addEventListener('abort', aborted, { once: true });
      timer = setTimeout(() => {
        signal?.removeEventListener('abort', aborted);
        resolve();
      }, Math.min(retryAfter * (attempt + 1), 2) * 1000);
    });
  }
}

