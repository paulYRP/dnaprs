(function () {
  "use strict";
  // Shared bounded cache. Chunk scripts also work in downloaded local reports.
  const cache = new Map();
  const pending = new Map();
  window.dnaprsTableChunk = function (key, rows) {
    cache.delete(key);
    cache.set(key, rows);
    while (cache.size > 4) cache.delete(cache.keys().next().value);
  };
  function load(meta, number) {
    const key = meta.key + ":" + (number + 1);
    if (cache.has(key)) {
      const rows = cache.get(key);
      window.dnaprsTableChunk(key, rows);
      return Promise.resolve(rows);
    }
    if (pending.has(key)) return pending.get(key);
    const promise = new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.src = meta.base + meta.chunks[number].file;
      const timeout = setTimeout(
        () => finish(new Error("Table loading timed out. Use the complete download or try again.")),
        30000,
      );
      function finish(error) {
        clearTimeout(timeout);
        script.remove();
        if (error) reject(error);
        else if (cache.has(key)) resolve(cache.get(key));
        else reject(new Error("The table chunk is unavailable. Use the complete download."));
      }
      script.onload = () => finish();
      script.onerror = () => finish(new Error("Unable to load table rows. Use the complete download."));
      document.head.appendChild(script);
    }).finally(() => pending.delete(key));
    pending.set(key, promise);
    return promise;
  }
  function matches(value, operator, query) {
    const text = value == null ? "" : String(value);
    if (operator === "contains") return text.toLowerCase().includes(query.toLowerCase());
    if (operator === "equals") return text === query;
    if (!text.trim() || !query.trim()) return false;
    const left = Number(text),
      right = Number(query);
    if (!Number.isFinite(left) || !Number.isFinite(right)) return false;
    return { lt: left < right, lte: left <= right, gt: left > right, gte: left >= right }[operator] || false;
  }
  const pause = () => new Promise((resolve) => setTimeout(resolve, 0));
  function initialise(root) {
    const meta = JSON.parse(root.querySelector("[data-table-metadata]").textContent);
    const q = (role) => root.querySelector("[data-" + role + "]");
    const status = q("status"),
      body = root.querySelector("tbody");
    let page = 1,
      indices = null,
      renderToken = 0,
      filterToken = 0,
      busy = false;
    const rowsControl = q("rows");
    rowsControl.value = root.dataset.pageSize;
    const size = () => Number(rowsControl.value);
    const total = () => (indices ? indices.length : meta.total);
    const pages = () => Math.max(1, Math.ceil(total() / size()));
    const heading = document.createElement("tr");
    meta.columns.forEach((column, i) => {
      const th = document.createElement("th");
      th.textContent = column;
      th.scope = "col";
      heading.appendChild(th);
      const option = document.createElement("option");
      option.textContent = column;
      option.value = i;
      q("column").appendChild(option);
    });
    root.querySelector("thead").appendChild(heading);
    function chunkFor(index) {
      let low = 0,
        high = meta.chunks.length - 1;
      while (low < high) {
        const mid = Math.ceil((low + high) / 2);
        if (meta.chunks[mid].start <= index) low = mid;
        else high = mid - 1;
      }
      return low;
    }
    function controls() {
      q("page").value = page;
      q("page").max = pages();
      q("first").disabled = q("previous").disabled = page <= 1;
      q("last").disabled = q("next").disabled = page >= pages();
      q("apply").disabled = busy;
      q("find").disabled = busy || !meta.identifier;
      q("cancel").disabled = !busy;
    }
    async function render() {
      if (!q("open").open) return;
      const token = ++renderToken;
      page = Math.max(1, Math.min(pages(), page));
      controls();
      const start = (page - 1) * size(),
        end = Math.min(total(), start + size());
      status.textContent = "Loading rows...";
      try {
        const fragment = document.createDocumentFragment();
        for (let i = start; i < end; i++) {
          const index = indices ? indices[i] : i;
          const number = chunkFor(index);
          const rows = await load(meta, number);
          if (token !== renderToken) return;
          const tr = document.createElement("tr");
          rows[index - meta.chunks[number].start].forEach((value) => {
            const td = document.createElement("td");
            td.textContent = value == null ? "" : String(value);
            tr.appendChild(td);
          });
          fragment.appendChild(tr);
        }
        body.replaceChildren(fragment);
        status.textContent = total()
          ? "Showing " +
            (start + 1) +
            "â€“" +
            end +
            " of " +
            total() +
            (indices ? " matching records (" + meta.total + " total)." : " records.")
          : "No matching records.";
      } catch (error) {
        if (token === renderToken) {
          body.replaceChildren();
          status.textContent = error.message;
        }
      }
    }
    async function filter() {
      const query = q("query").value,
        column = Number(q("column").value),
        operator = q("operator").value;
      if (!query.trim()) {
        status.textContent = "Enter a filter value.";
        return;
      }
      const token = ++filterToken;
      ++renderToken;
      busy = true;
      controls();
      // One bit per source row while scanning, then compact positions for matches.
      const mask = new Uint8Array(Math.ceil(meta.total / 8));
      let count = 0;
      try {
        for (let number = 0; number < meta.chunks.length; number++) {
          if (token !== filterToken) return;
          const descriptor = meta.chunks[number];
          if (
            operator === "equals" &&
            meta.columns[column] === meta.identifier &&
            descriptor.idMin != null &&
            descriptor.idMax != null &&
            /^[ -~]*$/.test(query) &&
            (query < descriptor.idMin || query > descriptor.idMax)
          )
            continue;
          status.textContent = "Filtering chunk " + (number + 1) + " of " + meta.chunks.length + "...";
          const rows = await load(meta, number);
          for (let i = 0; i < rows.length; i++) {
            if (i % 250 === 0) {
              await pause();
              if (token !== filterToken) return;
            }
            if (matches(rows[i][column], operator, query)) {
              const index = meta.chunks[number].start + i;
              mask[index >> 3] |= 1 << (index & 7);
              count++;
            }
          }
        }
        const selected = new Uint32Array(count);
        let offset = 0;
        for (let i = 0; i < meta.total; i++) {
          if (i % 50000 === 0) {
            await pause();
            if (token !== filterToken) return;
          }
          if (mask[i >> 3] & (1 << (i & 7))) selected[offset++] = i;
        }
        indices = selected;
        page = 1;
        q("open").open = true;
        await render();
      } catch (error) {
        if (token === filterToken) status.textContent = error.message;
      } finally {
        if (token === filterToken) {
          busy = false;
          controls();
        }
      }
    }
    q("apply").onclick = filter;
    q("find").onclick = () => {
      q("column").value = String(meta.columns.indexOf(meta.identifier));
      q("operator").value = "equals";
      return filter();
    };
    q("cancel").onclick = () => {
      ++filterToken;
      busy = false;
      controls();
      status.textContent = "Filtering cancelled. Previous selection retained.";
    };
    q("clear").onclick = () => {
      ++filterToken;
      busy = false;
      indices = null;
      page = 1;
      q("query").value = "";
      controls();
      render();
    };
    for (const [name, next] of Object.entries({
      first: () => 1,
      previous: () => page - 1,
      next: () => page + 1,
      last: pages,
    })) {
      q(name).onclick = () => {
        page = next();
        q("open").open = true;
        render();
      };
    }
    q("page").onchange = () => {
      page = Math.floor(Number(q("page").value)) || 1;
      render();
    };
    rowsControl.onchange = () => {
      page = 1;
      render();
    };
    q("open").ontoggle = () => {
      if (q("open").open) render();
      else {
        ++renderToken;
        ++filterToken;
        busy = false;
        controls();
        body.replaceChildren();
        for (const key of cache.keys()) if (key.startsWith(meta.key + ":")) cache.delete(key);
      }
    };
    controls();
  }
  function start() {
    document.querySelectorAll("[data-paged-table]").forEach(initialise);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start, { once: true });
  else start();
})();
