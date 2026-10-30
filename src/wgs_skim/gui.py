"""`wgs_skim gui`: tick the functions and options in a web page, save them as a preset (~/.wgs_skim/presets/NAME.yaml)
and get the command that runs it. Standard library only; listens on 127.0.0.1 (on a server:
`ssh -L 8765:localhost:8765 server`, then open http://localhost:8765 on your computer).

    wgs_skim gui [--port 8765]
"""
import json
import re
import secrets
import subprocess
import sys
import tempfile
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import yaml

from .paths import HOME, TAXA, presets, workflow_root

PRESET_DIR = HOME / "presets"
TOKEN = secrets.token_urlsafe(32)    # per session: other web pages in the browser cannot call this server


def default_keys():
    """the workflow options a preset may set (config/config.yaml)"""
    return set(yaml.safe_load((workflow_root() / "config" / "config.yaml").read_text()))


def check_preset(text):
    """(options, problems) of preset YAML text"""
    try:
        opts = yaml.safe_load(text) or {}
    except yaml.YAMLError as e:
        return {}, [f"not valid YAML: {e}"]
    if not isinstance(opts, dict):
        return {}, ["a preset is a set of key: value lines"]
    unknown = sorted(set(opts) - default_keys())
    return opts, [f"unknown option {k} (see `wgs_skim config`)" for k in unknown]


def dry_run(text, reads):
    """the jobs a run with this preset would start (snakemake -n), as text"""
    if not reads or not Path(reads).is_dir():
        return "give a reads directory that exists to preview the jobs"
    with tempfile.TemporaryDirectory(prefix="wgs_skim_gui_") as tmp:
        Path(tmp, "preset.yaml").write_text(text)
        r = subprocess.run([sys.executable, "-m", "wgs_skim", "run", "-r", reads, "-o", f"{tmp}/out", "-n",
                            "--preset", f"{tmp}/preset.yaml"], capture_output=True, text=True, timeout=300)
    out = r.stdout + r.stderr
    m = re.search(r"Job stats:\n(.*?\ntotal\s+\d+)", out, re.S)
    if m:
        return m.group(1)
    err = [x for x in out.splitlines() if "rror" in x or x.startswith("[wgs_skim]")]
    return "\n".join(err[-8:]) or out[-2000:]


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>wgs_skim presets</title>
<style>
:root{--bg:#f7f8fa;--card:#fff;--fg:#1d2433;--muted:#5b6475;--line:#e1e5ec;--accent:#2563eb;--ok:#15803d;--bad:#b91c1c}
@media (prefers-color-scheme:dark){:root{--bg:#111418;--card:#1a1e24;--fg:#e6e9ef;--muted:#9aa3b2;--line:#2b313b;
--accent:#60a5fa;--ok:#4ade80;--bad:#f87171}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.5 system-ui,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px}h1{font-size:20px;margin:0 0 4px}
p.sub{color:var(--muted);margin:0 0 20px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card h2{font-size:15px;margin:0 0 10px}
label.chk{display:flex;gap:8px;align-items:flex-start;padding:4px 0;cursor:pointer}
label.chk small{display:block;color:var(--muted)}
.opt{display:grid;grid-template-columns:150px 1fr;gap:6px 10px;align-items:center;margin-top:8px}
.opt label{color:var(--muted)}
input[type=text],input[type=number],select,textarea{width:100%;padding:6px 8px;border:1px solid var(--line);
border-radius:6px;background:var(--bg);color:var(--fg);font:inherit}
textarea{font-family:ui-monospace,monospace;font-size:12.5px;min-height:220px}
.bar{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:16px 0}
button{padding:7px 14px;border-radius:6px;border:1px solid var(--accent);background:var(--accent);color:#fff;
font:inherit;cursor:pointer}button.sec{background:transparent;color:var(--accent)}
pre{background:var(--bg);border:1px solid var(--line);border-radius:6px;padding:10px;overflow:auto;
font-size:12.5px;white-space:pre-wrap;margin:6px 0 0}
.msg.ok{color:var(--ok)}.msg.bad{color:var(--bad)}
</style></head><body><main>
<h1>wgs_skim presets</h1>
<p class="sub">Tick what to run and save it as a preset; then <code>wgs_skim run -r rawdata/ -o results/ --preset NAME</code>.
Presets are saved in <code>__PRESET_DIR__</code>. Options left empty keep the defaults (and your ~/.wgs_skim/config.yaml).</p>

<div class="bar"><label>Load preset <select id="load"><option value="">(new)</option>__PRESETS__</select></label></div>

<div class="grid">
<section class="card"><h2>Reads</h2>
<label class="chk"><input type="checkbox" id="qc" checked><span>Quality control (fastp), once<small>untick when the
reads are already clean: tools read them where they are</small></span></label>
<div class="opt"><label>Taxon preset</label><select id="taxon"><option value="">(default: annelid)</option>__TAXA__</select></div>
</section>

<section class="card"><h2>Mitogenome</h2>
<label class="chk"><input type="checkbox" data-step="mito" checked><span>MEANGS + GetOrganelle<small>reseeded with the
previous result while not circular</small></span></label>
<div class="opt"><label>MEANGS seed</label><select id="meangs"><option value="">yes</option><option value="false">no (GetOrganelle's own seeds)</option></select>
<label>Reseeded rounds</label><input type="number" id="mt_iter_rounds" min="0" placeholder="3 (0 = one run)"></div>
<label class="chk"><input type="checkbox" data-step="annotate" checked><span>MitoM annotation</span></label>
<label class="chk"><input type="checkbox" data-step="topology" checked><span>Read-based circular / linear check</span></label>
<label class="chk"><input type="checkbox" data-step="mitofinder"><span>MitoFinder<small>on the SPAdes assembly (or the
transcripts), with gene tables of all samples and the BLAST contamination check</small></span></label>
<div class="opt"><label>MitoFinder reference</label><input type="text" id="mitofinder_ref" placeholder="GenBank file or directory (empty: each sample's MitoM annotation)"></div>
</section>

<section class="card"><h2>Nuclear markers</h2>
<label class="chk"><input type="checkbox" data-step="nr" checked><span>rDNA (18S-ITS-28S) and histones (H1-H4)
<small>GetOrganelle, MitoM annotation (.gbf)</small></span></label>
<label class="chk"><input type="checkbox" data-step="nuclear" checked><span>Copies in the SPAdes assembly, copy numbers,
comparison with GetOrganelle</span></label>
</section>

<section class="card"><h2>Whole genome</h2>
<label class="chk"><input type="checkbox" data-step="spades" checked><span>SPAdes assembly</span></label>
<label class="chk"><input type="checkbox" data-step="busco" checked><span>BUSCO</span></label>
<div class="opt"><label>BUSCO lineage</label><input type="text" id="busco_lineage" placeholder="name or directory (empty: from the taxon)"></div>
<label class="chk"><input type="checkbox" data-step="genome_size" checked><span>Genome size and heterozygosity (GenomeScope)</span></label>
</section>

<section class="card"><h2>Transcriptome (RNA-seq)</h2>
<label class="chk"><input type="checkbox" data-step="rnaspades"><span>rnaSPAdes</span></label>
<label class="chk"><input type="checkbox" data-step="trinity"><span>Trinity</span></label>
<p class="sub" style="margin:6px 0 0">BUSCO and MitoFinder above apply to the transcripts too.</p>
</section>

<section class="card"><h2>Species check and report</h2>
<label class="chk"><input type="checkbox" data-step="blast" checked><span>BLAST (top hits with GenBank taxonomy)</span></label>
<div class="opt"><label>Mitogenome queries</label><input type="text" id="blast_mt" placeholder="cox1 (full, cox1, 16s, cytb ... or none)">
<label>rDNA queries</label><input type="text" id="blast_nr" placeholder="18S,ITS,28S (or none)">
<label>Local database</label><input type="text" id="blast_db" placeholder="empty: NCBI online / your config">
<label>E-mail (NCBI)</label><input type="text" id="email" placeholder="you@example.org"></div>
<label class="chk"><input type="checkbox" data-step="summary" checked><span>Summary table and HTML report</span></label>
</section>
</div>

<div class="bar">
<label>Preset name <input type="text" id="name" placeholder="e.g. my_skim" style="width:180px"></label>
<label>Description <input type="text" id="desc" placeholder="one line" style="width:320px"></label>
</div>
<section class="card"><h2>Preset (you can edit it before saving)</h2>
<textarea id="yaml" spellcheck="false"></textarea>
<div class="bar"><button id="save">Save preset</button>
<label>Reads directory <input type="text" id="reads" placeholder="/path/to/rawdata" style="width:300px"></label>
<button class="sec" id="preview">Preview jobs (dry run)</button></div>
<div id="msg" class="msg"></div><pre id="out" hidden></pre>
</section>
</main>
<script>
const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
const OPTS = ["taxon", "meangs", "mt_iter_rounds", "mitofinder_ref", "busco_lineage", "blast_mt", "blast_nr", "blast_db", "email"];
function build() {
  const steps = $$("[data-step]").filter(x => x.checked).map(x => x.dataset.step);
  const lines = [`# ${$("#desc").value || ($("#name").value || "my preset")}`,
                 `steps: "${steps.join(",")}"`];
  if (!$("#qc").checked) lines.push("qc: false");
  for (const k of OPTS) {
    const v = $("#" + k).value.trim();
    if (v === "") continue;
    lines.push(k === "mt_iter_rounds" || k === "meangs" ? `${k}: ${v}` : `${k}: "${v.replace(/"/g, '\\"')}"`);
  }
  $("#yaml").value = lines.join("\n") + "\n";
}
document.addEventListener("input", e => { if (e.target.id !== "yaml" && e.target.id !== "reads") build(); });
const TOKEN = "__TOKEN__";
async function post(url, body) {
  const r = await fetch(url, {method: "POST", headers: {"Content-Type": "application/json", "X-Token": TOKEN},
                              body: JSON.stringify(body)});
  return r.json();
}
function say(text, ok) { $("#msg").textContent = text; $("#msg").className = "msg " + (ok ? "ok" : "bad"); }
$("#save").onclick = async () => {
  const res = await post("/save", {name: $("#name").value, yaml: $("#yaml").value});
  say(res.message, res.ok);
  if (res.ok) { $("#out").hidden = false; $("#out").textContent = res.command; }
};
$("#preview").onclick = async () => {
  say("dry run ...", true);
  const res = await post("/preview", {yaml: $("#yaml").value, reads: $("#reads").value});
  say(res.message, res.ok); $("#out").hidden = false; $("#out").textContent = res.jobs;
};
$("#load").onchange = async () => {
  const name = $("#load").value; if (!name) return;
  const res = await (await fetch("/preset?name=" + encodeURIComponent(name), {headers: {"X-Token": TOKEN}})).json();
  const o = res.options || {}, steps = String(o.steps || "").split(",").map(s => s.trim());
  $$("[data-step]").forEach(x => x.checked = steps.includes(x.dataset.step) || (steps[0] === "" && x.dataset.step !== "rnaspades" && x.dataset.step !== "trinity" && (x.dataset.step !== "mitofinder" || o.mitofinder)));
  $("#qc").checked = o.qc !== false;
  for (const k of OPTS) $("#" + k).value = o[k] === undefined ? "" : String(o[k]);
  $("#name").value = name; $("#desc").value = res.description || "";
  $("#yaml").value = res.text;
};
build();
</script></body></html>"""


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        data = body.encode() if isinstance(body, str) else json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", f"{ctype}; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, fmt, *args):            # quiet
        pass

    def _allowed(self, need_token=True):
        """same machine, this server's own page: Host is 127.0.0.1 / localhost (DNS rebinding), the session token
        is sent (other sites' requests), POST bodies are JSON"""
        port = self.server.server_address[1]
        if self.headers.get("Host") not in (f"127.0.0.1:{port}", f"localhost:{port}"):
            return False
        if need_token and not secrets.compare_digest(self.headers.get("X-Token", ""), TOKEN):
            return False
        if self.command == "POST" and self.headers.get("Content-Type", "").split(";")[0] != "application/json":
            return False
        return True

    def do_GET(self):
        if self.path.split("?")[0] == "/" and self._allowed(need_token=False):
            if self.path != f"/?t={TOKEN}":
                return self._send(403, "open the address printed by `wgs_skim gui` (it carries a session token)",
                                  "text/plain")
            opts = "".join(f'<option value="{n}">{n}</option>' for n in presets())
            taxa = "".join(f'<option value="{t}">{t}</option>' for t in TAXA)
            page = PAGE.replace("__PRESETS__", opts).replace("__TAXA__", taxa).replace("__PRESET_DIR__", str(PRESET_DIR))
            page = page.replace("__TOKEN__", TOKEN)
            return self._send(200, page, "text/html")
        if not self._allowed():
            return self._send(403, {"error": "forbidden"})
        if self.path.startswith("/preset?name="):
            from urllib.parse import unquote
            name = unquote(self.path.split("=", 1)[1])
            f = presets().get(name)
            if not f:
                return self._send(404, {"error": "no such preset"})
            text = f.read_text()
            desc = next((x.lstrip("# ").strip() for x in text.splitlines() if x.startswith("#")), "")
            return self._send(200, {"text": text, "options": yaml.safe_load(text) or {}, "description": desc})
        self._send(404, {"error": "not found"})

    def do_POST(self):
        if not self._allowed():
            return self._send(403, {"ok": False, "message": "forbidden"})
        body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        if self.path == "/save":
            name = str(body.get("name", "")).strip()
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
                return self._send(200, {"ok": False, "message": "preset name: letters, digits, _ . - only"})
            if (workflow_root() / "config" / "presets" / f"{name}.yaml").exists():
                return self._send(200, {"ok": False, "message": f"'{name}' is a built-in preset: choose another name"})
            _, problems = check_preset(body.get("yaml", ""))
            if problems:
                return self._send(200, {"ok": False, "message": "; ".join(problems)})
            PRESET_DIR.mkdir(parents=True, exist_ok=True)
            f = PRESET_DIR / f"{name}.yaml"
            f.write_text(body["yaml"])
            return self._send(200, {"ok": True, "message": f"saved {f}",
                                    "command": f"wgs_skim run -r /path/to/rawdata -o results/ --preset {name}\n"
                                               f"wgs_skim run -r /path/to/rawdata -o results/ --preset {name} -n"
                                               "    # dry run: what would run"})
        if self.path == "/preview":
            _, problems = check_preset(body.get("yaml", ""))
            if problems:
                return self._send(200, {"ok": False, "message": "; ".join(problems), "jobs": ""})
            jobs = dry_run(body.get("yaml", ""), str(body.get("reads", "")).strip())
            return self._send(200, {"ok": True, "message": "jobs of a run with this preset:", "jobs": jobs})
        self._send(404, {"error": "not found"})


def serve(port=8765, open_browser=True):
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{port}/?t={TOKEN}"
    print(f"[wgs_skim] preset editor at {url}  (Ctrl-C to stop; on a server: ssh -L {port}:localhost:{port} <server>, "
          f"then open the address with localhost)", flush=True)
    if open_browser:
        try:
            webbrowser.open(url)
        except Exception:
            pass
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0
