const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const crypto = require("node:crypto");
const vm = require("node:vm");

async function checkDownloads() {
    const site = path.resolve(process.argv[2] || "_site");
    const context = { window: {} };
    vm.runInNewContext(fs.readFileSync(path.join(site, "assets/download-index.js"), "utf8"), context);
    const index = context.window.dnaprsDownloadIndex;
    const entries = Object.entries(index);
    assert.ok(entries.length);
    let total = 0;
    for (const [href, metadata] of entries) {
        const expected = crypto.createHash("sha256");
        for await (const bytes of fs.createReadStream(path.join(site, href))) expected.update(bytes);
        const actual = crypto.createHash("sha256");
        let size = 0;
        for (let number = 1; number <= metadata.chunks; number++) {
            const script = fs.readFileSync(
                path.join(site, metadata.base, String(number).padStart(6, "0") + ".js"),
                "utf8",
            );
            const match = script.match(/^document\.currentScript\.downloadBytes="([A-Za-z0-9+/=]*)";\r?\n?$/);
            assert.ok(match, `Invalid download asset for ${href}`);
            const bytes = Buffer.from(match[1], "base64");
            assert.equal(bytes.length, Math.min(metadata.chunkBytes, metadata.bytes - size));
            size += bytes.length;
            actual.update(bytes);
        }
        assert.equal(size, metadata.bytes, href);
        assert.equal(actual.digest("hex"), expected.digest("hex"), href);
        total += size;
    }
    console.log(
        `Download integrity passed: ${entries.length} files, ${total} bytes; every SHA-256 matches its original.`,
    );
}
checkDownloads().catch((error) => {
    console.error(error);
    process.exitCode = 1;
});
