const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

// A DOM-contract test, not a visual browser test.
class Element {
    constructor() {
        this.dataset = {};
        this.style = {};
        this.listeners = {};
        this.value = "1";
        this.clientWidth = 1000;
        this.clientHeight = 600;
        this.scrollWidth = 1000;
        this.scrollHeight = 600;
        this.naturalWidth = 1600;
        this.naturalHeight = 900;
    }
    addEventListener(type, handler) {
        (this.listeners[type] ||= []).push(handler);
    }
    emit(type, event = {}) {
        for (const handler of this.listeners[type] || []) handler(event);
    }
    cloneNode() {
        return new Element();
    }
    removeAttribute(name) {
        delete this[name];
    }
    setAttribute(name, value) {
        this[name] = value;
    }
    replaceWith(value) {
        elements["[data-figure-image]"] = value;
    }
    scrollTo() {}
}
const elements = {};
const gallery = new Element();
gallery.querySelector = (selector) => (elements[selector] ||= new Element());
const figures = ["MDD", "GAD"].map((trait) => ({
    title: `${trait} - QQ plot`,
    description: "Full-data ranks",
    inspection: "Inspect the diagonal",
    preview: `${trait}-preview.png`,
    svg: `${trait}.svg`,
    png: `${trait}.png`,
    tiff: `${trait}.tiff`,
    jpeg: `${trait}.jpeg`,
    source_table: `${trait}.tsv`,
}));
gallery.querySelector("[data-figure-data]").textContent = JSON.stringify(figures);
const script = fs.readFileSync("assets/report/dnaprs-report.js", "utf8");
// Evaluate the actual gallery function, without unrelated theme and table DOM.
const start = script.indexOf("  function initializeGallery(");
const end = script.indexOf("  function closeExpandedPanels(", start);
const context = {
    gallery,
    sharedZoom: 1,
    minimumZoom: 0.25,
    maximumZoom: 4,
    zoomStep: 1.15,
    clamp: (n, min, max) => Math.min(max, Math.max(min, n)),
    storeZoom() {},
    setExpandState() {},
    document: { addEventListener() {} },
    window: {
        requestAnimationFrame(fn) {
            fn();
            return 1;
        },
        cancelAnimationFrame() {},
        addEventListener() {},
    },
};
vm.runInNewContext(`${script.slice(start, end)}\ninitializeGallery(gallery);`, context);
const get = (name) => gallery.querySelector(`[data-figure-${name}]`);
const first = get("image");
assert.equal(first.src, figures[0].preview);
assert.equal(first.loading, "eager");
assert.match(get("status").textContent, /Loading/);
get("select").value = "2";
get("select").emit("change");
const second = get("image");
first.emit("load");
assert.equal(get("image"), second);
assert.match(get("status").textContent, /Loading/);
second.emit("error");
assert.match(get("status").textContent, /could not be loaded/);
for (const type of ["svg", "png", "tiff", "jpeg"]) assert.equal(get(type).href, figures[1][type]);
assert.equal(get("source").href, figures[1].source_table);
get("previous").emit("click");
const retry = get("image");
retry.emit("load");
assert.equal(retry.hidden, false);
assert.match(get("status").textContent, /Web preview/);
get("zoom-in").emit("click");
assert.match(get("zoom-status").textContent, /115%/);
get("zoom-reset").emit("click");
assert.equal(get("zoom-status").textContent, "Fit");
console.log("Report gallery passed: selected preview, stale loads, errors, downloads and zoom.");
