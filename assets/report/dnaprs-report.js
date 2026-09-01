(function () {
  "use strict";

  var root = document.documentElement;
  var minimumZoom = 0.25;
  var maximumZoom = 4;
  var zoomStep = 1.15;
  var sharedZoom = readStoredZoom();

  function clamp(value, minimum, maximum) {
    return Math.min(maximum, Math.max(minimum, value));
  }

  function readStoredZoom() {
    var value = null;
    try {
      value = Number(window.sessionStorage.getItem("dnaprs-figure-zoom"));
    } catch (error) {
      value = null;
    }
    return Number.isFinite(value) && value > 0 ? clamp(value, minimumZoom, maximumZoom) : 1;
  }

  function storeZoom(value) {
    sharedZoom = clamp(value, minimumZoom, maximumZoom);
    try {
      window.sessionStorage.setItem("dnaprs-figure-zoom", String(sharedZoom));
    } catch (error) {
      // Session storage can be unavailable for a local report.
    }
    document.dispatchEvent(new CustomEvent("dnaprs:figure-zoom", { detail: { zoom: sharedZoom } }));
  }

  function preferredTheme() {
    var stored = null;
    try {
      stored = window.localStorage.getItem("dnaprs-theme");
    } catch (error) {
      stored = null;
    }
    if (stored === "light" || stored === "dark") {
      return stored;
    }
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }

  function themeIcon(theme) {
    return theme === "dark"
      ? "<path d='M12 3a9 9 0 1 0 9 9 7 7 0 0 1-9-9Z'/>"
      : "<path d='M12 8a4 4 0 1 0 0 8 4 4 0 0 0 0-8Zm0-6h1v4h-2V2h1Zm0 16h1v4h-2v-4h1ZM2 11h4v2H2v-2Zm16 0h4v2h-4v-2ZM4.2 5.6l1.4-1.4 2.8 2.8L7 8.4 4.2 5.6Zm11.4 11.4 1.4-1.4 2.8 2.8-1.4 1.4-2.8-2.8Zm0-10 2.8-2.8 1.4 1.4L17 8.4 15.6 7ZM4.2 18.4 7 15.6 8.4 17l-2.8 2.8-1.4-1.4Z'/>";
  }

  function setTheme(theme, button) {
    var next = theme === "dark" ? "light" : "dark";
    root.dataset.theme = theme;
    button.setAttribute("aria-label", "Use " + next + " appearance");
    button.setAttribute("title", "Use " + next + " appearance");
    button.querySelector("svg").innerHTML = themeIcon(theme);
    document.querySelectorAll("[data-dnaprs-logo]").forEach(function (logo) {
      logo.src = theme === "dark" ? "assets/nf-core-dnaprs_logo_dark.svg" : "assets/nf-core-dnaprs_logo_light.svg";
    });
  }

  function initializeTheme() {
    document.querySelectorAll("#quarto-header .navbar-logo, #quarto-header .navbar-logo-image, #quarto-header img").forEach(function (logo) {
      logo.setAttribute("data-dnaprs-logo", "");
    });
    var button = document.createElement("button");
    button.type = "button";
    button.className = "dnaprs-theme-button";
    button.setAttribute("data-dnaprs-theme", "");
    button.innerHTML = "<svg viewBox='0 0 24 24' aria-hidden='true'></svg>";
    document.body.appendChild(button);
    setTheme(preferredTheme(), button);
    button.addEventListener("click", function () {
      var theme = root.dataset.theme === "dark" ? "light" : "dark";
      setTheme(theme, button);
      try {
        window.localStorage.setItem("dnaprs-theme", theme);
      } catch (error) {
        // The selected theme remains active for the current page.
      }
    });
  }

  function updateBodyLock() {
    var expanded = document.querySelector(".dnaprs-table-viewer.is-expanded, .dnaprs-figure-card.is-expanded");
    document.body.classList.toggle("dnaprs-lock-scroll", Boolean(expanded));
  }

  function setExpandState(container, button, expanded, label) {
    container.classList.toggle("is-expanded", expanded);
    button.setAttribute("aria-label", expanded ? "Close expanded " + label : "Expand " + label);
    button.setAttribute("title", expanded ? "Close expanded " + label : "Expand " + label);
    button.textContent = expanded ? "Close" : "⛶";
    updateBodyLock();
  }

  function initializeTable(viewer) {
    if (viewer.dataset.tableInitialized === "true") return;
    viewer.dataset.tableInitialized = "true";
    var rows = Array.from(viewer.querySelectorAll("tbody tr"));
    viewer.querySelectorAll("tbody td").forEach(function (cell) {
      var status = cell.textContent.trim().toUpperCase().replace(/_/g, " ");
      if (["PASS", "REVIEW", "FAIL", "NOT RUN"].indexOf(status) >= 0) {
        cell.classList.add("dnaprs-status", "dnaprs-status-" + status.toLowerCase().replace(" ", "-"));
      }
    });
    var search = viewer.querySelector("[data-table-search]");
    var size = viewer.querySelector("[data-table-size]");
    var previous = viewer.querySelector("[data-table-previous]");
    var next = viewer.querySelector("[data-table-next]");
    var status = viewer.querySelector("[data-table-status]");
    var expand = viewer.querySelector("[data-table-expand]");
    var page = 1;
    var filtered = rows.slice();

    function render() {
      var pageSize = Number(size.value) || 25;
      var pages = Math.max(1, Math.ceil(filtered.length / pageSize));
      page = Math.min(Math.max(1, page), pages);
      rows.forEach(function (row) {
        row.hidden = true;
      });
      filtered.slice((page - 1) * pageSize, page * pageSize).forEach(function (row) {
        row.hidden = false;
      });
      status.textContent = filtered.length
        ? "Page " + page + " of " + pages + " | " + filtered.length + " records"
        : "No matching records";
      previous.disabled = page <= 1;
      next.disabled = page >= pages;
    }

    function applySearch() {
      var query = search.value.trim().toLocaleLowerCase();
      filtered = rows.filter(function (row) {
        return !query || row.textContent.toLocaleLowerCase().indexOf(query) >= 0;
      });
      page = 1;
      render();
    }

    search.addEventListener("input", applySearch);
    size.addEventListener("change", function () {
      page = 1;
      render();
    });
    previous.addEventListener("click", function () {
      page -= 1;
      render();
    });
    next.addEventListener("click", function () {
      page += 1;
      render();
    });
    expand.addEventListener("click", function () {
      setExpandState(viewer, expand, !viewer.classList.contains("is-expanded"), "table");
    });
    render();
  }

  function initializeLogText(target) {
    var panel = target.closest("details");
    if (!panel) return;

    function loadText() {
      if (!panel.open || target.dataset.logLoaded === "true" || target.dataset.logLoading === "true") return;
      target.dataset.logLoading = "true";
      target.textContent = "Loading completed execution record...";
      fetch(target.dataset.logSource, { cache: "no-store" })
        .then(function (response) {
          if (!response.ok) throw new Error("HTTP " + response.status);
          return response.text();
        })
        .then(function (value) {
          target.textContent = value || "The completed execution record is empty.";
          target.dataset.logLoaded = "true";
        })
        .catch(function () {
          target.textContent = "The browser could not load this local text record. Use the download link below, or serve the report directory with a local web server.";
        })
        .finally(function () {
          target.dataset.logLoading = "false";
        });
    }

    panel.addEventListener("toggle", loadText);
    loadText();
  }

  function initializeGallery(gallery) {
    if (gallery.dataset.galleryInitialized === "true") return;
    var data = gallery.querySelector("[data-figure-data]");
    var figures = [];
    try {
      figures = JSON.parse(data.textContent);
    } catch (error) {
      return;
    }
    if (!figures.length) return;
    gallery.dataset.galleryInitialized = "true";

    var select = gallery.querySelector("[data-figure-select]");
    var previous = gallery.querySelector("[data-figure-previous]");
    var next = gallery.querySelector("[data-figure-next]");
    var count = gallery.querySelector("[data-figure-count]");
    var title = gallery.querySelector("[data-figure-title]");
    var description = gallery.querySelector("[data-figure-description]");
    var inspection = gallery.querySelector("[data-figure-inspection]");
    var image = gallery.querySelector("[data-figure-image]");
    var canvas = gallery.querySelector(".dnaprs-figure-canvas");
    var stage = gallery.querySelector("[data-figure-stage]");
    var svg = gallery.querySelector("[data-figure-svg]");
    var tiff = gallery.querySelector("[data-figure-tiff]");
    var png = gallery.querySelector("[data-figure-png]");
    var jpeg = gallery.querySelector("[data-figure-jpeg]");
    var source = gallery.querySelector("[data-figure-source]");
    var card = gallery.querySelector(".dnaprs-figure-card");
    var expand = gallery.querySelector("[data-figure-expand]");
    var zoomOut = gallery.querySelector("[data-figure-zoom-out]");
    var zoomReset = gallery.querySelector("[data-figure-zoom-reset]");
    var zoomIn = gallery.querySelector("[data-figure-zoom-in]");
    var zoomStatus = gallery.querySelector("[data-figure-zoom-status]");
    var zoom = sharedZoom;
    var geometryFrame = 0;

    function updateZoomStatus() {
      zoomStatus.textContent = Math.abs(zoom - 1) < 0.001 ? "Fit" : Math.round(zoom * 100) + "%";
    }

    function renderGeometry() {
      if (!image.naturalWidth || !image.naturalHeight) return;
      var availableWidth = Math.max(1, canvas.clientWidth - 16);
      var availableHeight = Math.max(1, canvas.clientHeight - 16);
      var fitScale = Math.min(availableWidth / image.naturalWidth, availableHeight / image.naturalHeight);
      var imageWidth = Math.max(1, image.naturalWidth * fitScale * zoom);
      var imageHeight = Math.max(1, image.naturalHeight * fitScale * zoom);
      image.style.width = imageWidth + "px";
      image.style.height = imageHeight + "px";
      stage.style.width = Math.max(availableWidth, imageWidth) + "px";
      stage.style.height = Math.max(availableHeight, imageHeight) + "px";
      updateZoomStatus();
    }

    function scheduleGeometry() {
      if (geometryFrame) window.cancelAnimationFrame(geometryFrame);
      geometryFrame = window.requestAnimationFrame(function () {
        geometryFrame = 0;
        renderGeometry();
      });
    }

    function changeZoom(nextZoom, offsetX, offsetY) {
      var oldWidth = Math.max(1, canvas.scrollWidth);
      var oldHeight = Math.max(1, canvas.scrollHeight);
      var anchorX = (canvas.scrollLeft + offsetX) / oldWidth;
      var anchorY = (canvas.scrollTop + offsetY) / oldHeight;
      zoom = clamp(nextZoom, minimumZoom, maximumZoom);
      storeZoom(zoom);
      renderGeometry();
      window.requestAnimationFrame(function () {
        canvas.scrollLeft = anchorX * canvas.scrollWidth - offsetX;
        canvas.scrollTop = anchorY * canvas.scrollHeight - offsetY;
      });
    }

    function resetZoom() {
      zoom = 1;
      storeZoom(zoom);
      renderGeometry();
      canvas.scrollTo(0, 0);
    }

    function render() {
      var index = clamp(Number(select.value) - 1, 0, figures.length - 1);
      var figure = figures[index];
      title.textContent = figure.title;
      description.textContent = figure.description;
      inspection.textContent = figure.inspection;
      image.src = figure.svg || figure.png;
      image.alt = figure.title + ". " + figure.description;
      svg.href = figure.svg;
      tiff.href = figure.tiff;
      png.href = figure.png;
      jpeg.href = figure.jpeg;
      svg.download = figure.svg.split("/").pop();
      tiff.download = figure.tiff.split("/").pop();
      png.download = figure.png.split("/").pop();
      jpeg.download = figure.jpeg.split("/").pop();
      if (figure.source_table) {
        source.hidden = false;
        source.href = figure.source_table;
        source.download = figure.source_table.split("/").pop();
      } else {
        source.hidden = true;
        source.removeAttribute("href");
        source.removeAttribute("download");
      }
      count.textContent = index + 1 + " of " + figures.length + " figures";
      previous.disabled = index === 0;
      next.disabled = index === figures.length - 1;
      zoom = sharedZoom;
      canvas.scrollTo(0, 0);
      if (image.complete && image.naturalWidth) scheduleGeometry();
    }

    function change(offset) {
      select.value = String(clamp(Number(select.value) + offset, 1, figures.length));
      render();
    }

    previous.addEventListener("click", function () {
      change(-1);
    });
    next.addEventListener("click", function () {
      change(1);
    });
    select.addEventListener("change", render);
    image.addEventListener("load", scheduleGeometry);
    zoomOut.addEventListener("click", function () {
      changeZoom(zoom / zoomStep, canvas.clientWidth / 2, canvas.clientHeight / 2);
    });
    zoomIn.addEventListener("click", function () {
      changeZoom(zoom * zoomStep, canvas.clientWidth / 2, canvas.clientHeight / 2);
    });
    zoomReset.addEventListener("click", resetZoom);
    canvas.addEventListener(
      "wheel",
      function (event) {
        if (event.ctrlKey || event.metaKey) {
          event.preventDefault();
          var bounds = canvas.getBoundingClientRect();
          changeZoom(
            zoom * (event.deltaY < 0 ? zoomStep : 1 / zoomStep),
            event.clientX - bounds.left,
            event.clientY - bounds.top,
          );
        } else if (event.shiftKey && !event.deltaX) {
          event.preventDefault();
          canvas.scrollLeft += event.deltaY;
        }
      },
      { passive: false },
    );
    canvas.addEventListener("keydown", function (event) {
      if (event.key === "+" || event.key === "=") {
        event.preventDefault();
        changeZoom(zoom * zoomStep, canvas.clientWidth / 2, canvas.clientHeight / 2);
      } else if (event.key === "-") {
        event.preventDefault();
        changeZoom(zoom / zoomStep, canvas.clientWidth / 2, canvas.clientHeight / 2);
      } else if (event.key === "0") {
        event.preventDefault();
        resetZoom();
      }
    });
    canvas.addEventListener("dblclick", resetZoom);
    gallery.addEventListener("keydown", function (event) {
      if (event.target.matches("input,select,button,a")) return;
      if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
        event.preventDefault();
        change(event.key === "ArrowRight" ? 1 : -1);
      }
    });
    expand.addEventListener("click", function () {
      setExpandState(card, expand, !card.classList.contains("is-expanded"), "figure");
      scheduleGeometry();
    });
    document.addEventListener("dnaprs:figure-zoom", function (event) {
      zoom = clamp(Number(event.detail.zoom) || 1, minimumZoom, maximumZoom);
      scheduleGeometry();
    });
    window.addEventListener("resize", scheduleGeometry);
    if (window.ResizeObserver) {
      new ResizeObserver(scheduleGeometry).observe(canvas);
    }
    render();
  }

  function closeExpandedPanels(event) {
    if (event.key !== "Escape") return;
    document.querySelectorAll(".dnaprs-table-viewer.is-expanded").forEach(function (viewer) {
      setExpandState(viewer, viewer.querySelector("[data-table-expand]"), false, "table");
    });
    document.querySelectorAll(".dnaprs-figure-card.is-expanded").forEach(function (card) {
      setExpandState(card, card.querySelector("[data-figure-expand]"), false, "figure");
    });
  }

  function initializeAll() {
    initializeTheme();
    document.querySelectorAll("[data-table-viewer]").forEach(initializeTable);
    document.querySelectorAll("[data-figure-gallery]").forEach(initializeGallery);
    document.querySelectorAll("[data-log-source]").forEach(initializeLogText);
    document.addEventListener("keydown", closeExpandedPanels);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initializeAll, { once: true });
  } else {
    initializeAll();
  }
})();
