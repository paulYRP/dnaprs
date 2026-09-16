(function () {
  "use strict";
  var pending = false;

  function loadChunk(path) {
    return new Promise(function (resolve, reject) {
      var script = document.createElement("script");
      var timer = window.setTimeout(function () {
        finish(new Error("Download timed out."));
      }, 60000);
      function finish(error, value) {
        window.clearTimeout(timer);
        script.onload = script.onerror = null;
        script.downloadBytes = null;
        script.remove();
        if (error) reject(error);
        else resolve(value);
      }
      script.onload = function () {
        try {
          if (typeof script.downloadBytes !== "string") throw new Error("Missing download data.");
          var binary = window.atob(script.downloadBytes);
          var bytes = new Uint8Array(binary.length);
          for (var i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
          finish(null, bytes);
        } catch (error) {
          finish(error);
        }
      };
      script.onerror = function () {
        finish(new Error("Download data could not be loaded."));
      };
      script.src = path;
      document.head.appendChild(script);
    });
  }

  async function prepareFile(entry, progress) {
    var parts = [];
    var total = 0;
    for (var n = 1; n <= entry.chunks; n++) {
      var bytes = await loadChunk(entry.base + String(n).padStart(6, "0") + ".js");
      var expected = Math.min(entry.chunkBytes, entry.bytes - total);
      if (bytes.length !== expected) throw new Error("Incomplete download data.");
      total += bytes.length;
      parts.push(new Blob([bytes]));
      progress(Math.round((100 * total) / entry.bytes));
    }
    if (total !== entry.bytes) throw new Error("Incomplete download data.");
    return new Blob(parts, { type: "application/octet-stream" });
  }

  function statusFor(link) {
    var status = link.parentNode.querySelector("[data-download-status]");
    if (!status) {
      status = document.createElement("span");
      status.dataset.downloadStatus = "true";
      status.setAttribute("role", "status");
      status.setAttribute("aria-live", "polite");
      link.parentNode.appendChild(status);
    }
    return status;
  }

  document.addEventListener("click", async function (event) {
    var link = event.target.closest("a[download]");
    if (!link || link.dataset.preparedDownload === "true") return;
    var href = link.getAttribute("href");
    if (!href) return;
    var index = window.dnaprsDownloadIndex || {};
    var key = decodeURIComponent(href);
    // Quarto can make local links absolute after loading the page.
    if (!index[key]) {
      var absolute = new URL(href, document.baseURI).href;
      key = Object.keys(index).find(function (candidate) {
        return new URL(candidate.split("/").map(encodeURIComponent).join("/"), document.baseURI).href === absolute;
      });
    }
    var entry = index[key];
    if (!entry) return;
    event.preventDefault();
    var status = statusFor(link);
    if (pending) {
      status.textContent = " Wait for the current download to finish preparing.";
      return;
    }
    pending = true;
    // Snapshot the selection; changing figures during preparation cannot rename this file.
    var filename = link.getAttribute("download") || key.split("/").pop();
    link.setAttribute("aria-busy", "true");
    status.textContent = " Preparing download...";
    try {
      var blob = await prepareFile(entry, function (percent) {
        status.textContent = " Preparing download: " + percent + "%";
      });
      var url = URL.createObjectURL(blob);
      var save = document.createElement("a");
      save.href = url;
      save.download = filename;
      save.dataset.preparedDownload = "true";
      save.hidden = true;
      document.body.appendChild(save);
      try {
        save.click();
      } finally {
        save.remove();
        window.setTimeout(function () {
          URL.revokeObjectURL(url);
        }, 60000);
      }
      status.textContent = "";
    } catch (error) {
      status.textContent = " Download could not be prepared. Keep the complete report folder and try again.";
    } finally {
      link.removeAttribute("aria-busy");
      pending = false;
    }
  });
})();
