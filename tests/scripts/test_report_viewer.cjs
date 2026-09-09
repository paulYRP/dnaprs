// Isolated DOM contract tests; no browser, network or participant data.
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const assert = require("node:assert/strict");
const project = path.resolve(process.argv[2] || "report-table-test");
const html = fs.readFileSync(path.join(project, "index.qmd"), "utf8");
const metadata = JSON.parse(html.match(/data-table-metadata>(.*?)<\/script>/s)[1]);
class Element {
    constructor() {
        this.children = [];
        this.value = "";
        this.open = false;
        this.disabled = false;
        this.textContent = "";
    }
    appendChild(child) {
        if (child.fragment) this.children.push(...child.children);
        else this.children.push(child);
        return child;
    }
    replaceChildren(...children) {
        this.children = [];
        children.forEach((child) => this.appendChild(child));
    }
    remove() {}
}
const roles = new Map();
for (const name of [
    "table-metadata",
    "status",
    "rows",
    "page",
    "first",
    "previous",
    "next",
    "last",
    "apply",
    "find",
    "cancel",
    "clear",
    "open",
    "query",
    "column",
    "operator",
])
    roles.set(name, new Element());
roles.get("table-metadata").textContent = JSON.stringify(metadata);
roles.get("column").value = "0";
roles.get("operator").value = "equals";
const body = new Element(),
    head = new Element();
const root = {
    dataset: { pageSize: "25" },
    querySelector(selector) {
        if (selector === "tbody") return body;
        if (selector === "thead") return head;
        return roles.get(selector.slice(6, -1));
    },
};
let loads = 0,
    maximum = 0,
    failLoad = false;
const document = {
    readyState: "complete",
    querySelectorAll: () => [root],
    createElement: () => new Element(),
    createDocumentFragment: () => Object.assign(new Element(), { fragment: true }),
    head: {
        appendChild(script) {
            loads++;
            setImmediate(() => {
                if (failLoad) {
                    script.onerror();
                    return;
                }
                vm.runInContext(fs.readFileSync(path.join(project, script.src), "utf8"), context);
                script.onload();
            });
        },
    },
};
const maps = [];
class ObservedMap extends Map {
    constructor(...args) {
        super(...args);
        maps.push(this);
    }
}
const context = vm.createContext({
    document,
    window: {},
    setTimeout,
    clearTimeout,
    console,
    Uint8Array,
    Uint32Array,
    Map: ObservedMap,
});
vm.runInContext(fs.readFileSync("assets/report/dnaprs-tables.js", "utf8"), context);
const original = context.window.dnaprsTableChunk;
context.window.dnaprsTableChunk = (key, rows) => {
    maximum = Math.max(maximum, rows.length);
    original(key, rows);
    assert.ok(maps[0].size <= 4, "shared chunk cache exceeds four entries");
};
const delay = () => new Promise((resolve) => setTimeout(resolve, 10));
async function settled() {
    for (let i = 0; i < 300; i++) {
        await delay();
        if (!roles.get("status").textContent.startsWith("Loading")) return;
    }
    throw new Error("Viewer did not settle");
}
(async () => {
    assert.equal(loads, 0, "closed tables must not load chunks");
    roles.get("open").open = true;
    roles.get("open").ontoggle();
    await settled();
    assert.equal(body.children.length, 25);
    assert.equal(body.children[0].children[0].textContent, "rs10");
    assert.equal(body.children[2].children[2].textContent, "</script><b>literal</b>");
    roles.get("last").onclick();
    await settled();
    assert.equal(body.children.length, 12);
    assert.match(roles.get("status").textContent, /76â€“87 of 87/);
    roles.get("first").onclick();
    await settled();
    assert.equal(body.children.length, 25);
    roles.get("query").value = "rs10";
    await roles.get("find").onclick();
    assert.match(roles.get("status").textContent, /58 matching/);
    roles.get("column").value = "1";
    roles.get("operator").value = "gte";
    roles.get("query").value = "80";
    await roles.get("apply").onclick();
    assert.equal(body.children.length, 7);
    roles.get("query").value = "999";
    await roles.get("apply").onclick();
    assert.equal(body.children.length, 0);
    roles.get("clear").onclick();
    await settled();
    assert.equal(body.children.length, 25);
    roles.get("query").value = "1";
    const filtering = roles.get("apply").onclick();
    roles.get("cancel").onclick();
    await filtering;
    assert.match(roles.get("status").textContent, /cancelled/);
    roles.get("open").open = false;
    roles.get("open").ontoggle();
    assert.equal(body.children.length, 0);
    failLoad = true;
    roles.get("open").open = true;
    roles.get("open").ontoggle();
    await settled();
    assert.match(roles.get("status").textContent, /complete download/);
    assert.ok(maximum <= 17);
    assert.ok(loads < 100);
    console.log(
        "Viewer contract tests passed: lazy loading, paging, repeated IDs, numeric filters, cancellation, text escaping and missing chunks.",
    );
})().catch((error) => {
    console.error(error);
    process.exitCode = 1;
});
