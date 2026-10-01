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

// A response may already be queued when abort() runs. Only the current request
// in each group may change the UI, even if a cancelled fetch still resolves.
const requests = (() => {
  const activeRequests = new Map();
  function cancelRequests(...keys) {
    for (const key of keys) {
      const pending = activeRequests.get(key);
      if (!pending) continue;
      activeRequests.delete(key);
      pending.controller.abort();
      pending.onEnd();
    }
  }
  function startRequest(key, onEnd = () => {}) {
    cancelRequests(key);
    let settle;
    const done = new Promise(resolve => { settle = resolve; });
    const pending = { controller: new AbortController(), done,
      onEnd() { try { onEnd(); } finally { settle(); } },
    };
    activeRequests.set(key, pending);
    return {
      signal: pending.controller.signal,
      done,
      current: () => activeRequests.get(key) === pending,
      finish() {
        if (activeRequests.get(key) !== pending) return;
        activeRequests.delete(key);
        pending.onEnd();
      },
    };
  }

  return {
    start: startRequest, cancel: cancelRequests,
    has: key => activeRequests.has(key),
    pending: (...keys) => keys.map(key => activeRequests.get(key)?.done).find(Boolean),
  };
})();

async function readResponse(response) {
  const data = await response.json();
  if (!response.ok && !data.error) {
    const detail = data.detail;
    data.error = Array.isArray(detail)
      ? detail.map(item => item.msg).join('; ')
      : detail || 'The request could not be completed.';
  }
  return data;
}

