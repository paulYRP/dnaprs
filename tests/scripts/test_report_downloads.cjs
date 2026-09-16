const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const { setTimeout: delay } = require("node:timers/promises");

// Exercise the real offline handler without fetching local files or starting a server.
async function testDownloads() {
    let clickHandler;
    let requests = 0;
    let missing = false;
    const saves = [];
    const urls = new Map();
    const status = {};
    const original = Buffer.from([0, 255, 128, 65, 13, 10, 0, 1, 2]);
    const chunks = [original.subarray(0, 4), original.subarray(4, 8), original.subarray(8)];
    const context = {
        Blob,
        Uint8Array,
        URL: class extends URL {
            static createObjectURL(blob) {
                const key = `blob:${urls.size}`;
                urls.set(key, blob);
                return key;
            }
            static revokeObjectURL(key) {
                urls.delete(key);
            }
        },
        window: {
            atob: (text) => Buffer.from(text, "base64").toString("binary"),
            setTimeout(fn, ms) {
                const timer = setTimeout(fn, ms);
                timer.unref();
                return timer;
            },
            clearTimeout,
            dnaprsDownloadIndex: {
                "downloads/a & b.pgen": { base: "assets/downloads/test/", chunks: 3, bytes: 9, chunkBytes: 4 },
            },
        },
        document: {
            baseURI: "file:///report/index.html",
            addEventListener(type, fn) {
                assert.equal(type, "click");
                clickHandler = fn;
            },
            createElement(tag) {
                return {
                    dataset: {},
                    remove() {},
                    setAttribute() {},
                    click() {
                        assert.equal(tag, "a");
                        saves.push({ name: this.download, blob: urls.get(this.href) });
                    },
                };
            },
            head: {
                appendChild(script) {
                    requests++;
                    const number = Number(script.src.match(/(\d+)\.js$/)[1]);
                    setTimeout(() => {
                        if (missing) script.onerror();
                        else {
                            script.downloadBytes = chunks[number - 1].toString("base64");
                            script.onload();
                        }
                    }, 1);
                },
            },
            body: { appendChild() {} },
        },
    };
    vm.runInNewContext(fs.readFileSync("assets/report/dnaprs-downloads.js", "utf8"), context);
    assert.equal(requests, 0, "Opening the report must not preload download bytes");
    const link = {
        dataset: {},
        name: "a & b.pgen",
        busy: false,
        getAttribute(key) {
            return key === "href" ? "file:///report/downloads/a%20%26%20b.pgen" : this.name;
        },
        setAttribute() {
            this.busy = true;
        },
        removeAttribute() {
            this.busy = false;
        },
        parentNode: {
            querySelector() {
                return status;
            },
        },
    };
    const event = { target: { closest: () => link }, preventDefault() {} };
    const first = clickHandler(event);
    link.name = "new-selection.pgen";
    await clickHandler(event);
    await first;
    assert.equal(saves.length, 1);
    assert.equal(saves[0].name, "a & b.pgen");
    assert.deepEqual(Buffer.from(await saves[0].blob.arrayBuffer()), original);
    assert.equal(requests, 3);
    assert.equal(link.busy, false);
    assert.equal(status.textContent, "", "A successful handoff must clear preparation text");
    missing = true;
    await clickHandler(event);
    assert.match(status.textContent, /could not be prepared/);
    assert.equal(saves.length, 1);
    assert.equal(link.busy, false);
    missing = false;
    await clickHandler(event);
    assert.equal(saves.length, 2, "A failed download must be retryable");
    assert.equal(status.textContent, "", "A successful retry must clear the error text");
    await delay(5);
    console.log("Offline downloads passed: lazy loading, exact binary bytes, names, duplicate clicks and retry.");
}
testDownloads().catch((error) => {
    console.error(error);
    process.exitCode = 1;
});
